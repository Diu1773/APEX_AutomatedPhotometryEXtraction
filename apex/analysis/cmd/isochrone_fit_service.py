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
    dm_prior: Optional[Tuple[float, float]] = None        # (m-M)0 Gaussian prior
    parallax_distance_prior: bool = False                 # auto dm prior from Gaia parallax members
    parallax_dm_sigma: float = 0.05                       # sigma (mag) for the auto parallax dm prior
    parallax_dm_window: float = 0.06                      # hard (m-M)0 half-window around the parallax dm
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
    """Everything the GUI needs to render and persist a fit.

    The plot *data* (obs/iso arrays, chain) is always populated so a caller can
    build the figures on the main thread (Matplotlib/pyplot is not thread-safe).
    ``cmd_figure``/``corner_figure`` are only filled when ``make_figures=True``.
    """

    result: Any                             # McmcResult
    summary: Dict[str, Any]
    member_meta: Dict[str, Any]
    n_stars: int
    # Plot data (always populated) so figures can be built on the main thread.
    obs_color: Any = None                   # np.ndarray (colour 0)
    obs_mag: Any = None                     # np.ndarray (magnitude axis)
    iso_color: Any = None                   # np.ndarray (best-fit track, colour 0)
    iso_mag: Any = None                     # np.ndarray (best-fit track, magnitude)
    color_label: str = ""
    mag_label: str = ""
    annotations: Dict[str, Any] = field(default_factory=dict)
    cmd_figure: Any = None                  # Matplotlib Figure (if make_figures)
    corner_figure: Any = None               # Matplotlib Figure (if make_figures)
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
        "convergence_detail": str(getattr(result, "convergence_detail", "")),
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

    # Distance constraint from the Gaia parallax of the membership clump. The gri
    # degeneracy makes the likelihood prefer a spurious SHORT distance by a large
    # margin (~200+ logL on M67), so a soft Gaussian prior alone is NOT enough — we
    # must also TIGHTEN the dm bounds to a window around the parallax distance. The
    # cluster distance is known to ~0.03 mag from the member-averaged parallax, so a
    # ±0.1 mag hard window is well justified and is what actually breaks the rail.
    dm_lo, dm_hi = config.dm_bounds
    dm_prior_pair: Optional[Tuple[float, float]] = None
    if config.dm_prior is not None:
        dm_prior_pair = (float(config.dm_prior[0]), float(config.dm_prior[1]))
        dm_c, w = dm_prior_pair[0], max(5.0 * dm_prior_pair[1], 0.15)
        dm_lo, dm_hi = max(dm_lo, dm_c - w), min(dm_hi, dm_c + w)
    elif config.parallax_distance_prior and member_meta.get("applied"):
        plx = float(member_meta.get("parallax_center") or 0.0)
        if plx > 0.0:
            dm_plx = 5.0 * np.log10(1000.0 / plx) - 5.0
            w = float(config.parallax_dm_window)
            dm_lo, dm_hi = max(dm_lo, dm_plx - w), min(dm_hi, dm_plx + w)
            dm_prior_pair = (dm_plx, float(config.parallax_dm_sigma))
            member_meta["dm_from_parallax"] = round(dm_plx, 3)
            member_meta["dm_window"] = [round(dm_lo, 3), round(dm_hi, 3)]

    # Reddening constraint: when a reddening (E(colour)) prior is supplied (e.g. a
    # dust-map value), TIGHTEN the E bounds around it for the same reason — the
    # degeneracy overrides a soft Gaussian, so dm AND E both need hard windows for
    # the fit to recover the correct age (verified: pinning both gives M67 age 3.99
    # Gyr; pinning dm only leaves age/E trading off, age drifts young).
    ec_lo, ec_hi = config.ecolor_bounds
    if config.ecolor_prior is not None:
        ec_c, ec_w = float(config.ecolor_prior[0]), max(2.0 * float(config.ecolor_prior[1]), 0.008)
        ec_lo, ec_hi = max(ec_lo, ec_c - ec_w), min(ec_hi, ec_c + ec_w)

    # [M/H] constraint: a soft Gaussian (±few logL) is overridden by a confident
    # likelihood (same as dm/E), so when an [M/H] prior is supplied — e.g. a
    # high-resolution spectroscopic [Fe/H] — also tighten the [M/H] bounds to
    # ±2.5σ. The window scales with the prior's σ, so a precise spectroscopic
    # value pins [M/H] hard while a loose 'assumed solar' stays broad.
    mh_lo, mh_hi = config.mh_bounds
    if config.mh_prior is not None:
        mh_c, mh_w = float(config.mh_prior[0]), max(2.0 * float(config.mh_prior[1]), 0.04)
        mh_lo, mh_hi = max(mh_lo, mh_c - mh_w), min(mh_hi, mh_c + mh_w)

    bounds = FitBounds(
        log_age=config.age_bounds, metallicity=(mh_lo, mh_hi),
        distance_mod=(dm_lo, dm_hi), extinction_gr=(ec_lo, ec_hi),
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
    # Distance prior (computed above with the matching dm-bounds tightening).
    if dm_prior_pair is not None:
        priors["dm"] = dm_prior_pair

    _say(0.45, f"Running MCMC ({config.n_walkers}×{config.n_steps})")

    def _mcmc_progress(frac, label):
        # Map the sampler's 0..1 onto the 0.45..0.88 slice of the overall bar.
        _say(0.45 + 0.43 * float(frac), f"{label} {int(frac * 100)}%")

    result = fit_isochrone_mcmc(
        obs_c, obs_m, color_err, mag_err,
        bounds=bounds,
        interp_fn=fitter._interpolate_isochrone,
        imf_fn=fitter._imf_weight,
        n_walkers=config.n_walkers, n_burn=config.n_burn, n_steps=config.n_steps,
        init_from_gridscan=init, priors=priors or None,
        f_bin=config.f_bin, f_field=config.f_field,
        R_color=R_color, seed=config.seed, progress=False, progress_cb=_mcmc_progress,
        colors=colors, mag=mag_band, multiband_iso=mb,
        obs_multicolor=obs_mc, err_multicolor=err_mc,
        err_floor=config.err_floor,
    )

    summary = _summary_dict(result, 1.0 / R_color if R_color else 1.0)
    if not result.convergence_ok:
        warnings.append(
            f"MCMC not converged — "
            f"{getattr(result, 'convergence_detail', '') or 'reason unavailable'}"
            f" (acceptance {result.acceptance_fraction:.2f}); "
            "treat the posterior with caution — see methodology §10."
        )

    out = IsochroneFitOutput(
        result=result, summary=summary, member_meta=member_meta, n_stars=n_stars,
        warnings=warnings,
    )

    # Always compute the best-fit track + annotations as plain arrays/dicts so the
    # caller can render figures on the main thread (pyplot is not thread-safe).
    median = result.median_theta if result.median_theta is not None else [
        result.log_age_med, result.mh_med, result.dm_med, result.e_color_med]
    e_bv_med = float(median[3]) / R_color if R_color else float(median[3])
    # Overlay = the SINGLE nearest-grid isochrone (no 4-corner mass-blend, which
    # fans into a mess of near-horizontal lines), reddened + distance-shifted —
    # exactly what the CMD viewer draws. mass-ordered -> one clean MS→giant curve.
    try:
        ages_g = np.unique(np.round(iso_data[:, 2], 4))
        mhs_g = np.unique(np.round(iso_data[:, 1], 4))
        a0 = float(ages_g[np.argmin(np.abs(ages_g - float(median[0])))])
        m0 = float(mhs_g[np.argmin(np.abs(mhs_g - float(median[1])))])
        sel = (np.round(iso_data[:, 2], 4) == a0) & (np.round(iso_data[:, 1], 4) == m0)
        sub = iso_data[sel]
        # KEEP FILE (EEP/evolutionary) ORDER — do NOT sort. PARSEC isochrones are
        # EEP-parametrised, so the same initial mass recurs at different phases;
        # sorting by mass interleaves the phases and fans the line into a mess.
        dm = float(median[2])

        def _app(band):
            return sub[:, bcol[band]] + dm + D.EXTINCTION_R[band] * e_bv_med
        out.iso_color = _app(b1) - _app(b2)
        out.iso_mag = _app(mag_band)
    except Exception:
        bands_app, _mass = mb.apparent(float(median[0]), float(median[1]),
                                       float(median[2]), e_bv_med)
        if bands_app is not None:
            out.iso_color = bands_app[b1] - bands_app[b2]
            out.iso_mag = bands_app[mag_band]
    # Display ALL members (not just the SNR/max-stars fit subsample) so the full
    # CMD — bright giants and the faint MS — is visible; the fit itself still used
    # the subsample (obs_mc).
    try:
        dv1, _e1, _s1 = D.magnitude_column(work, b1)
        dv2, _e2, _s2 = D.magnitude_column(work, b2)
        dvm, _em, _sm = D.magnitude_column(work, mag_band)
        dcol = dv1 - dv2
        dfin = np.isfinite(dcol) & np.isfinite(dvm)
        out.obs_color, out.obs_mag = dcol[dfin], dvm[dfin]
    except Exception:
        out.obs_color, out.obs_mag = obs_c, obs_m
    # Truncate the displayed isochrone just above the brightest observed member.
    # The post-turn-off EEP blocks (thermal-pulse AGB, post-AGB → WD) oscillate
    # wildly in colour and, for young isochrones, dip back into the plotted
    # magnitude range as a tangle, while representing a negligible, essentially
    # unobserved population. Keep the EEP sequence from the faint MS up to the first
    # point that climbs ~0.5 mag brighter than the brightest member, then stop — so
    # the overlay is one clean MS→turn-off→(sub)giant curve bounded by the data.
    try:
        if out.iso_mag is not None and out.obs_mag is not None and len(out.iso_mag) > 5:
            ic = np.asarray(out.iso_color, float)
            im = np.asarray(out.iso_mag, float)
            bright_limit = float(np.nanmin(out.obs_mag)) - 0.5
            faint_limit = float(np.nanmax(out.obs_mag)) + 0.5
            # EEP order runs faint MS -> turn-off -> RGB -> loops/AGB/post-AGB/WD.
            # Locate the turn-off (bluest point WITHIN the plotted magnitude range, so
            # the very blue/faint white-dwarf tail can't masquerade as it), then cut at
            # the first post-turn-off point that climbs brighter than the brightest
            # member: that drops the RGB tip / thermal-pulse-AGB / post-AGB tangle while
            # protecting the MS+turn-off. For giant-less (young) clusters the crossing
            # falls on the MS at/before the turn-off, so nothing is truncated.
            inrange = np.where((im >= bright_limit) & (im <= faint_limit))[0]
            crossed = np.where(im < bright_limit)[0]
            if len(inrange) and len(crossed):
                to_idx = int(inrange[int(np.nanargmin(ic[inrange]))])
                cut = int(crossed[0])
                if cut > to_idx and cut >= 5:
                    out.iso_color, out.iso_mag = ic[:cut], im[:cut]
    except Exception:
        pass
    out.color_label, out.mag_label = labels[0], mag_band
    ann = {"age_gyr": summary["age_gyr"], "mh": summary["metallicity"],
           "dm": summary["distance_mod"]}
    if summary.get("e_bv"):
        ann["e_bv"] = summary["e_bv"]
    out.annotations = ann

    if make_figures:
        _say(0.9, "Rendering figures")
        from apex.analysis.cmd import isochrone_plots as P

        out.cmd_figure = P.cmd_figure(
            out.obs_color, out.obs_mag,
            out.iso_color if out.iso_color is not None else np.array([]),
            out.iso_mag if out.iso_mag is not None else np.array([]),
            color_label=out.color_label, mag_label=out.mag_label,
            title=f"Auto isochrone fit ({'+'.join(labels[:-1])})",
            annotations=out.annotations,
        )
        if result.flat_chain is not None and result.labels is not None:
            out.corner_figure = P.corner_figure(result.flat_chain, list(result.labels))

    _say(1.0, "Done")
    return out
