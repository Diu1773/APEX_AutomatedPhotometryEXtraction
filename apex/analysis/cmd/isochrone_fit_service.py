"""High-level isochrone-fit service: dataframe in → posterior + figures out.

This is the single entry point the GUI (and any caller) uses to run the automatic
MCMC isochrone fit. It composes the Qt-free pieces — :mod:`isochrone_data`
(membership + array extraction + isochrone loading), :func:`fit_isochrone_mcmc`
(the sampler), and :mod:`isochrone_plots` (the paper-quality figures) — so the GUI
worker only has to call one function off the main thread and then drop the
returned Matplotlib figures into a canvas.

Qt-free by construction (analysis layer); the GUI imports it, never vice-versa.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from apex.analysis.cmd import isochrone_data as D
from apex.analysis.cmd.isochrone_fitter_v2 import FitBounds, IsochroneFitterV2

__all__ = ["IsochroneFitConfig", "IsochroneFitOutput", "fit_cluster_isochrone"]


@dataclass
class IsochroneFitConfig:
    """User-facing knobs for the automatic isochrone fit."""

    colors: List[Tuple[str, str]]          # e.g. [("g", "r")] or [("g", "r"), ("r", "i")]
    mag_band: str = "g"
    iso_file: Optional[Path] = None        # default bundled PARSEC grid if None
    age_bounds: Tuple[float, float] = (8.5, 10.2)
    mh_bounds: Tuple[float, float] = (-1.0, 0.5)
    dm_bounds: Tuple[float, float] = (5.0, 18.0)
    ecolor_bounds: Tuple[float, float] = (0.0, 0.5)
    # Gaussian priors (mean, sigma); a tiny sigma ≈ a hard "fix".
    mh_prior: Optional[Tuple[float, float]] = None
    ecolor_prior: Optional[Tuple[float, float]] = None
    use_membership: bool = True
    data_snr_min: float = 20.0
    fit_snr_min: float = 5.0
    max_stars: int = 500
    n_walkers: int = 32
    n_burn: int = 600
    n_steps: int = 2000
    f_bin: float = 0.3
    f_field: float = 0.1
    err_floor: float = 0.02
    seed: int = 2024


@dataclass
class IsochroneFitOutput:
    """Everything the GUI needs to render and persist a fit."""

    result: Any                             # McmcResult
    summary: Dict[str, Any]
    member_meta: Dict[str, Any]
    n_stars: int
    cmd_figure: Any = None                  # Matplotlib Figure
    corner_figure: Any = None               # Matplotlib Figure
    warnings: List[str] = field(default_factory=list)


def _summary_dict(result, e_color_to_ebv: float) -> Dict[str, Any]:
    def trip(p):
        return [float(p[0]), float(p[1]), float(p[2])] if p is not None else None

    out = {
        "log_age": trip(result.log_age_p),
        "metallicity": trip(result.mh_p),
        "distance_mod": trip(result.dm_p),
        "e_color": trip(result.e_color_p),
        "age_gyr": trip(result.age_gyr_p),
        "distance_pc": trip(result.distance_pc_p),
        "acceptance_fraction": float(result.acceptance_fraction),
        "convergence_ok": bool(result.convergence_ok),
    }
    e_bv = trip(getattr(result, "e_bv_p", None))
    if e_bv is None and out["e_color"] is not None:
        e_bv = [v * e_color_to_ebv for v in out["e_color"]]
    out["e_bv"] = e_bv
    return out


def fit_cluster_isochrone(
    df: pd.DataFrame,
    config: IsochroneFitConfig,
    *,
    make_figures: bool = True,
    progress_cb=None,
) -> IsochroneFitOutput:
    """Run the automatic MCMC isochrone fit on a Step-10 wide-CMD dataframe.

    Parameters
    ----------
    df :
        Step-10 wide CMD dataframe (``mag_std_*`` columns; Gaia astrometry
        optional, used only for membership).
    config :
        :class:`IsochroneFitConfig` knobs.
    make_figures :
        Build the CMD + corner figures (set False for a headless number-only run).
    progress_cb :
        Optional ``callable(fraction: float, message: str)`` for GUI feedback.
    """
    from apex.analysis.cmd.isochrone_mcmc import MultiBandIsochrone, fit_isochrone_mcmc

    def _say(frac, msg):
        if progress_cb is not None:
            progress_cb(float(frac), str(msg))

    warnings: List[str] = []
    colors = [tuple(c) for c in config.colors]
    mag_band = config.mag_band
    bands = sorted({b for pair in colors for b in pair} | {mag_band})

    _say(0.05, "Selecting members" if config.use_membership else "Selecting stars")
    member_meta: Dict[str, Any] = {"applied": False, "reason": "disabled"}
    work = df
    if config.use_membership:
        mmask, member_meta = D.membership_mask(df)
        if member_meta.get("applied"):
            work = df.iloc[np.flatnonzero(mmask)].reset_index(drop=True)
        else:
            warnings.append(f"Membership not applied: {member_meta.get('reason')}")

    _say(0.15, "Extracting photometry")
    obs_mc, err_mc, labels, data_meta = D.extract_multicolor_arrays(
        work, colors, mag_band,
        data_snr_min=config.data_snr_min, max_stars=config.max_stars, seed=config.seed,
    )
    n_stars = int(len(obs_mc))
    if n_stars < 20:
        raise ValueError(f"Too few stars after selection (n={n_stars}); relax SNR/membership.")

    iso_file = Path(config.iso_file) if config.iso_file else D.default_iso_file(colors[0])
    if not Path(iso_file).exists():
        raise FileNotFoundError(f"Isochrone file not found: {iso_file}")

    _say(0.25, "Loading isochrone grid")
    iso_data, bcol, _ = D.load_multiband_isochrone(
        iso_file, bands, config.age_bounds, config.mh_bounds,
    )
    b1, b2 = colors[0]
    R_color = D.EXTINCTION_R[b1] - D.EXTINCTION_R[b2]

    mb = MultiBandIsochrone(
        iso_data, col_age=2, col_mh=1, col_mass=iso_data.shape[1] - 1,
        band_cols={b: bcol[b] for b in bands},
        R_band={b: D.EXTINCTION_R[b] for b in bands},
    )
    # IsochroneFitterV2 (built on colour 0 / mag) supplies the IMF + single-colour
    # interpolation used to seed the walkers; the multi-band likelihood drives the fit.
    fitter = IsochroneFitterV2(
        Path(iso_file), col_mh=1, col_age=2, col_g=bcol[b1], col_r=bcol[b2],
        col_mag=bcol[mag_band], col_mass=iso_data.shape[1] - 1,
        R_band1=D.EXTINCTION_R[b1], R_band2=D.EXTINCTION_R[b2],
        R_mag=D.EXTINCTION_R[mag_band], iso_data=iso_data, use_bilinear=False,
    )

    bounds = FitBounds(
        log_age=config.age_bounds, metallicity=config.mh_bounds,
        distance_mod=config.dm_bounds, extinction_gr=config.ecolor_bounds,
    )

    obs_c, obs_m = obs_mc[:, 0], obs_mc[:, -1]
    color_err, mag_err = err_mc[:, 0], err_mc[:, -1]

    _say(0.35, "Grid-scan seed")
    try:
        grid = fitter.fit_grid_scan(obs_c, obs_m, color_err, mag_err,
                                    bounds=bounds, snr_min=config.fit_snr_min)
        gb = grid.best_fit
        init = [gb.log_age, gb.metallicity, gb.distance_mod, gb.extinction_gr]
    except Exception as exc:  # pragma: no cover - seed is best-effort
        init = None
        warnings.append(f"Grid-scan seed failed ({exc}); using bounds midpoint.")

    priors: Dict[str, Tuple[float, float]] = {}
    if config.mh_prior is not None:
        priors["mh"] = (float(config.mh_prior[0]), float(config.mh_prior[1]))
    if config.ecolor_prior is not None:
        priors["e_color"] = (float(config.ecolor_prior[0]), float(config.ecolor_prior[1]))

    _say(0.45, f"Running MCMC ({config.n_walkers}×{config.n_steps})")
    result = fit_isochrone_mcmc(
        obs_c, obs_m, color_err, mag_err,
        bounds=bounds,
        interp_fn=fitter._interpolate_isochrone,
        imf_fn=fitter._imf_weight,
        n_walkers=config.n_walkers, n_burn=config.n_burn, n_steps=config.n_steps,
        init_from_gridscan=init, priors=priors or None,
        f_bin=config.f_bin, f_field=config.f_field,
        R_color=R_color, seed=config.seed, progress=False,
        colors=colors, mag=mag_band, multiband_iso=mb,
        obs_multicolor=obs_mc, err_multicolor=err_mc,
        err_floor=config.err_floor,
    )

    summary = _summary_dict(result, 1.0 / R_color if R_color else 1.0)
    if not result.convergence_ok:
        warnings.append(
            f"MCMC not converged (acceptance {result.acceptance_fraction:.2f}); "
            "treat the posterior with caution — see methodology §10."
        )

    out = IsochroneFitOutput(
        result=result, summary=summary, member_meta=member_meta, n_stars=n_stars,
        warnings=warnings,
    )

    if make_figures:
        _say(0.9, "Rendering figures")
        from apex.analysis.cmd import isochrone_plots as P

        bands_app, mass = mb.apparent(*result.median_theta[:2],
                                      float(result.median_theta[2]),
                                      float(result.median_theta[3]) / R_color
                                      if R_color else float(result.median_theta[3]))
        iso_color = iso_mag = None
        if bands_app is not None:
            iso_color = bands_app[b1] - bands_app[b2]
            iso_mag = bands_app[mag_band]
            order = np.argsort(iso_mag)
            iso_color, iso_mag = iso_color[order], iso_mag[order]
        ann = {
            "age_gyr": summary["age_gyr"], "mh": summary["metallicity"],
            "dm": summary["distance_mod"],
        }
        if summary.get("e_bv"):
            ann["e_bv"] = summary["e_bv"]
        out.cmd_figure = P.cmd_figure(
            obs_c, obs_m,
            iso_color if iso_color is not None else np.array([]),
            iso_mag if iso_mag is not None else np.array([]),
            color_label=labels[0], mag_label=mag_band,
            title=f"Auto isochrone fit ({'+'.join(labels[:-1])})",
            annotations=ann,
        )
        if result.flat_chain is not None and result.labels is not None:
            out.corner_figure = P.corner_figure(result.flat_chain, list(result.labels))

    _say(1.0, "Done")
    return out
