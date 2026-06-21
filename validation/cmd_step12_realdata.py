"""Headless real-data validation runner for CMD Step12 isochrone fitting.

The runner reads an APEX result directory, extracts Step10 wide CMD photometry,
loads a compact PARSEC isochrone subset, runs the Step12 grid-scan fitter, and
writes a JSON report plus an optional CMD overlay plot.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from apex.analysis.cmd.isochrone_fitter_v2 import FitBounds, IsochroneFitterV2
from apex.utils.gaia_transforms import build_color_pairs, filter_bands_from_columns


PARSEC_ISO_COL = {
    "u": 28, "g": 29, "r": 30, "i": 31, "z": 32,
    "B": 30, "V": 31, "R": 32, "I": 33,
}

EXTINCTION_R = {
    "U": 4.902, "B": 4.035, "V": 3.116, "R": 2.634, "I": 1.903,
    "g": 3.303, "r": 2.285, "i": 1.698, "z": 1.263,
}

M5_SANITY_RANGES = {
    "log_age": (9.95, 10.15),
    "metallicity": (-1.75, -0.95),
    "distance_mod": (14.05, 14.75),
    "e_bv": (0.00, 0.08),
}


def parse_pair(text: str) -> tuple[str, str]:
    parts = [p.strip() for p in str(text).split("-")]
    if len(parts) != 2 or not all(parts):
        raise argparse.ArgumentTypeError("color must look like g-r or B-V")
    return parts[0], parts[1]


def parse_bounds(values: list[float]) -> tuple[float, float]:
    if len(values) != 2:
        raise argparse.ArgumentTypeError("bounds require two floats")
    lo, hi = float(values[0]), float(values[1])
    if hi <= lo:
        raise argparse.ArgumentTypeError("bounds must be increasing")
    return lo, hi


def resolve_wide_csv(result_dir: Path) -> Path:
    candidates = [
        result_dir / "cmd_zeropoint" / "median_by_ID_filter_wide_cmd.csv",
        result_dir / "cmd_zeropoint" / "median_by_ID_filter_wide.csv",
        result_dir / "median_by_ID_filter_wide_cmd.csv",
        result_dir / "median_by_ID_filter_wide.csv",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError("Step10 wide CMD CSV not found under result_dir")


def infer_color_mag(df: pd.DataFrame) -> tuple[tuple[str, str], str]:
    bands = filter_bands_from_columns(df.columns, "mag_std_")
    pairs = build_color_pairs(bands, adjacent_only=True)
    for pair, mag in [("g-r", "g"), ("B-V", "V"), ("V-I", "V"), ("r-i", "r")]:
        b1, b2 = pair.split("-")
        if (b1, b2) in pairs and mag in bands:
            return (b1, b2), mag
    if not pairs:
        raise ValueError("No mag_std_* bands found in Step10 CSV")
    b1, b2 = pairs[0]
    return (b1, b2), b1


def default_iso_file(color: tuple[str, str]) -> Path:
    if any(b in {"B", "V", "R", "I"} for b in color):
        return ROOT / "isochrone" / "PARSEC" / "johnson" / "iso_data_cmd39_johnson.dat"
    return ROOT / "isochrone" / "PARSEC" / "sdss" / "iso_data_cmd39.dat"


def magnitude_column(df: pd.DataFrame, band: str) -> tuple[np.ndarray, np.ndarray, str]:
    std = f"mag_std_{band}"
    inst = f"mag_inst_{band}"
    if std in df.columns:
        values = pd.to_numeric(df[std], errors="coerce").to_numpy(float)
        source = std
    elif inst in df.columns:
        values = pd.to_numeric(df[inst], errors="coerce").to_numpy(float)
        source = inst
    else:
        raise KeyError(f"Missing magnitude column for band {band}")

    err_candidates = [f"mag_std_err_{band}", f"mag_inst_err_{band}", f"mag_inst_med_err_{band}"]
    err = None
    for col in err_candidates:
        if col in df.columns:
            err = pd.to_numeric(df[col], errors="coerce").to_numpy(float)
            break
    if err is None:
        err = np.full(len(values), 0.03, dtype=float)
    err = np.where(np.isfinite(err) & (err > 0), err, 0.03)
    return values, err, source


def membership_mask(
    df: pd.DataFrame, *, nsig: float = 3.0, plx_nsig: float = 3.0,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Boolean mask of probable cluster members from Gaia proper motion + parallax.

    Auto-detects the cluster clump in (pmra, pmdec) as the densest 2-D histogram
    peak, refines the centre by iterative robust (MAD) clipping, then keeps stars
    within ``nsig`` robust sigma of the PM centroid AND ``plx_nsig`` of the
    members' median parallax. Stars lacking astrometry are treated as non-members
    (conservative: unconfirmable membership is excluded). Returns an all-True mask
    as a graceful no-op when the astrometry columns are absent.
    """
    need = {"pmra", "pmdec", "parallax"}
    if not need.issubset(df.columns):
        return np.ones(len(df), dtype=bool), {"applied": False, "reason": "no_astrometry"}
    pmra = pd.to_numeric(df["pmra"], errors="coerce").to_numpy(float)
    pmdec = pd.to_numeric(df["pmdec"], errors="coerce").to_numpy(float)
    plx = pd.to_numeric(df["parallax"], errors="coerce").to_numpy(float)
    fin = np.isfinite(pmra) & np.isfinite(pmdec) & np.isfinite(plx)
    if int(fin.sum()) < 20:
        return np.ones(len(df), dtype=bool), {"applied": False, "reason": "too_few_astrometric"}

    x, y, p = pmra[fin], pmdec[fin], plx[fin]
    lo_x, hi_x = np.percentile(x, [2, 98])
    lo_y, hi_y = np.percentile(y, [2, 98])
    bw = 0.5
    bx = np.arange(lo_x, hi_x + bw, bw)
    by = np.arange(lo_y, hi_y + bw, bw)
    if len(bx) < 2 or len(by) < 2:
        return np.ones(len(df), dtype=bool), {"applied": False, "reason": "degenerate_pm_range"}
    hist, xe, ye = np.histogram2d(x, y, bins=[bx, by])
    ip, jp = np.unravel_index(int(np.argmax(hist)), hist.shape)
    cx = 0.5 * (xe[ip] + xe[ip + 1])
    cy = 0.5 * (ye[jp] + ye[jp + 1])

    sel = np.hypot(x - cx, y - cy) < 1.5
    sx = sy = 0.3
    for _ in range(5):
        if int(sel.sum()) < 10:
            break
        cx = float(np.median(x[sel]))
        cy = float(np.median(y[sel]))
        sx = max(1.4826 * float(np.median(np.abs(x[sel] - cx))), 0.2)
        sy = max(1.4826 * float(np.median(np.abs(y[sel] - cy))), 0.2)
        sel = (np.abs(x - cx) < nsig * sx) & (np.abs(y - cy) < nsig * sy)

    pc = float(np.median(p[sel]))
    ps = max(1.4826 * float(np.median(np.abs(p[sel] - pc))), 0.05)
    member_fin = sel & (np.abs(p - pc) < plx_nsig * ps)

    mask = np.zeros(len(df), dtype=bool)
    mask[np.flatnonzero(fin)[member_fin]] = True
    meta = {
        "applied": True,
        "pm_center": [round(cx, 3), round(cy, 3)],
        "pm_sigma": [round(sx, 3), round(sy, 3)],
        "parallax_center": round(pc, 3),
        "parallax_sigma": round(ps, 3),
        "n_astrometric": int(fin.sum()),
        "n_members": int(mask.sum()),
        "nsig": nsig,
        "plx_nsig": plx_nsig,
    }
    return mask, meta


def extract_multicolor_arrays(
    df: pd.DataFrame,
    colors: list[tuple[str, str]],
    mag_band: str,
    data_snr_min: float,
    max_stars: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, list[str], dict[str, Any]]:
    """Extract a (N, D) observed array + errors for the multi-colour fit.

    Axis order: one column per colour (b1 - b2) then the magnitude band, matching
    :func:`fit_isochrone_mcmc`'s multi-colour convention.  Errors are propagated
    in quadrature for colours.  Prefers calibrated ``mag_std_*`` columns over
    instrumental ``mag_inst_*`` (handled by :func:`magnitude_column`) because the
    multi-colour fit is sensitive to per-band zeropoint offsets that a single
    colour absorbs into E(colour).
    """
    bands_needed = sorted({b for pair in colors for b in pair} | {mag_band})
    vals: dict[str, np.ndarray] = {}
    errs: dict[str, np.ndarray] = {}
    sources: dict[str, str] = {}
    for b in bands_needed:
        v, e, src = magnitude_column(df, b)
        vals[b], errs[b], sources[b] = v, e, src

    axes = []
    axis_errs = []
    labels = []
    for b1, b2 in colors:
        axes.append(vals[b1] - vals[b2])
        axis_errs.append(np.sqrt(errs[b1] ** 2 + errs[b2] ** 2))
        labels.append(f"{b1}-{b2}")
    axes.append(vals[mag_band])
    axis_errs.append(errs[mag_band])
    labels.append(mag_band)

    obs = np.column_stack(axes)
    err = np.column_stack(axis_errs)

    mask = np.isfinite(obs).all(axis=1) & np.isfinite(err).all(axis=1)
    for b in bands_needed:
        snr_col = f"snr_{b}"
        if snr_col in df.columns and data_snr_min > 0:
            snr = pd.to_numeric(df[snr_col], errors="coerce").to_numpy(float)
            mask &= np.isfinite(snr) & (snr >= data_snr_min)

    idx = np.flatnonzero(mask)
    if max_stars > 0 and len(idx) > max_stars:
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(idx, size=max_stars, replace=False))

    meta = {
        "n_total_rows": int(len(df)),
        "n_after_filters": int(len(idx)),
        "axes": labels,
        "columns": sources,
    }
    return obs[idx], err[idx], labels, meta


def load_multiband_isochrone(
    iso_file: Path,
    bands: list[str],
    age_bounds: tuple[float, float],
    mh_bounds: tuple[float, float],
    cache_dir: Path | None = None,
) -> tuple[np.ndarray, dict[int, int], dict[str, Any]]:
    """Load a compact PARSEC subset carrying every requested band.

    Returns ``(iso_data, band_col_in_subset, meta)`` where ``iso_data`` columns
    are ``[pad, mh, log_age, band_0, band_1, ..., mass]`` and
    ``band_col_in_subset`` maps band-name -> column index in ``iso_data``.
    """
    missing = [b for b in bands if b not in PARSEC_ISO_COL]
    if missing:
        raise KeyError(f"No PARSEC column mapping for: {', '.join(missing)}")

    band_src_cols = {b: PARSEC_ISO_COL[b] for b in bands}
    usecols = sorted({1, 2, 5, *band_src_cols.values()})

    cache_path = None
    if cache_dir is not None:
        stat = iso_file.stat()
        key_text = json.dumps(
            {
                "iso_file": str(iso_file.resolve()),
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "bands": bands,
                "age_bounds": age_bounds,
                "mh_bounds": mh_bounds,
                "kind": "multiband",
            },
            sort_keys=True,
        )
        digest = hashlib.sha1(key_text.encode("utf-8")).hexdigest()[:16]
        cache_path = cache_dir / f"parsec_multiband_{digest}.npy"
        if cache_path.exists():
            data = np.load(cache_path)
            band_col = {b: 3 + k for k, b in enumerate(bands)}
            return data, band_col, {"cache_hit": True, "cache_path": str(cache_path)}

    chunks: list[np.ndarray] = []
    for chunk in pd.read_csv(
        iso_file, sep=r"\s+", comment="#", header=None,
        usecols=usecols, chunksize=200_000, engine="c",
    ):
        age = pd.to_numeric(chunk[2], errors="coerce").to_numpy(float)
        mh = pd.to_numeric(chunk[1], errors="coerce").to_numpy(float)
        mask = (
            np.isfinite(age) & np.isfinite(mh)
            & (age >= age_bounds[0]) & (age <= age_bounds[1])
            & (mh >= mh_bounds[0]) & (mh <= mh_bounds[1])
        )
        if not mask.any():
            continue
        sub = chunk.loc[mask]
        cols = [np.zeros(len(sub)), sub[1].to_numpy(float), sub[2].to_numpy(float)]
        cols += [sub[band_src_cols[b]].to_numpy(float) for b in bands]
        cols.append(sub[5].to_numpy(float))
        chunks.append(np.column_stack(cols))

    if not chunks:
        raise ValueError("No isochrone rows survived the selected age/[M/H] bounds")
    iso_data = np.vstack(chunks)
    band_col = {b: 3 + k for k, b in enumerate(bands)}
    meta = {"cache_hit": False, "cache_path": str(cache_path) if cache_path else None}
    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.save(cache_path, iso_data)
    return iso_data, band_col, meta


def extract_cmd_arrays(
    df: pd.DataFrame,
    color: tuple[str, str],
    mag_band: str,
    data_snr_min: float,
    max_stars: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, pd.DataFrame, dict[str, Any]]:
    b1, b2 = color
    v1, e1, c1_source = magnitude_column(df, b1)
    v2, e2, c2_source = magnitude_column(df, b2)
    vm, em, mag_source = magnitude_column(df, mag_band)

    obs_color = v1 - v2
    obs_mag = vm
    obs_color_err = np.sqrt(e1**2 + e2**2)
    obs_mag_err = em

    mask = (
        np.isfinite(obs_color)
        & np.isfinite(obs_mag)
        & np.isfinite(obs_color_err)
        & np.isfinite(obs_mag_err)
    )
    for band in sorted(set([b1, b2, mag_band])):
        snr_col = f"snr_{band}"
        if snr_col in df.columns and data_snr_min > 0:
            snr = pd.to_numeric(df[snr_col], errors="coerce").to_numpy(float)
            mask &= np.isfinite(snr) & (snr >= data_snr_min)

    idx = np.flatnonzero(mask)
    if max_stars > 0 and len(idx) > max_stars:
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(idx, size=max_stars, replace=False))

    filtered = df.iloc[idx].copy()
    meta = {
        "n_total_rows": int(len(df)),
        "n_after_filters": int(len(idx)),
        "columns": {
            "band1": c1_source,
            "band2": c2_source,
            "mag": mag_source,
        },
    }
    return (
        obs_color[idx],
        obs_mag[idx],
        obs_color_err[idx],
        obs_mag_err[idx],
        filtered,
        meta,
    )


def load_compact_isochrone(
    iso_file: Path,
    color: tuple[str, str],
    mag_band: str,
    age_bounds: tuple[float, float],
    mh_bounds: tuple[float, float],
    cache_dir: Path | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    b1, b2 = color
    missing = [b for b in (b1, b2, mag_band) if b not in PARSEC_ISO_COL]
    if missing:
        raise KeyError(f"No PARSEC column mapping for: {', '.join(missing)}")

    col_b1 = PARSEC_ISO_COL[b1]
    col_b2 = PARSEC_ISO_COL[b2]
    col_mag = PARSEC_ISO_COL[mag_band]
    usecols = sorted({1, 2, 5, col_b1, col_b2, col_mag})

    cache_path = None
    if cache_dir is not None:
        stat = iso_file.stat()
        key_text = json.dumps(
            {
                "iso_file": str(iso_file.resolve()),
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "color": color,
                "mag_band": mag_band,
                "age_bounds": age_bounds,
                "mh_bounds": mh_bounds,
                "cols": [col_b1, col_b2, col_mag],
            },
            sort_keys=True,
        )
        digest = hashlib.sha1(key_text.encode("utf-8")).hexdigest()[:16]
        cache_path = cache_dir / f"parsec_subset_{digest}.npy"
        if cache_path.exists():
            return np.load(cache_path), {"cache_hit": True, "cache_path": str(cache_path)}

    chunks: list[np.ndarray] = []
    for chunk in pd.read_csv(
        iso_file,
        sep=r"\s+",
        comment="#",
        header=None,
        usecols=usecols,
        chunksize=200_000,
        engine="c",
    ):
        age = pd.to_numeric(chunk[2], errors="coerce").to_numpy(float)
        mh = pd.to_numeric(chunk[1], errors="coerce").to_numpy(float)
        mask = (
            np.isfinite(age)
            & np.isfinite(mh)
            & (age >= age_bounds[0])
            & (age <= age_bounds[1])
            & (mh >= mh_bounds[0])
            & (mh <= mh_bounds[1])
        )
        if not mask.any():
            continue
        sub = chunk.loc[mask]
        compact = np.column_stack(
            [
                np.zeros(len(sub), dtype=float),
                sub[1].to_numpy(float),
                sub[2].to_numpy(float),
                sub[col_b1].to_numpy(float),
                sub[col_b2].to_numpy(float),
                sub[col_mag].to_numpy(float),
                sub[5].to_numpy(float),
            ]
        )
        chunks.append(compact)

    if not chunks:
        raise ValueError("No isochrone rows survived the selected age/[M/H] bounds")
    iso_data = np.vstack(chunks)
    meta = {"cache_hit": False, "cache_path": str(cache_path) if cache_path else None}
    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.save(cache_path, iso_data)
    return iso_data, meta


def build_fitter(
    iso_file: Path,
    iso_data: np.ndarray,
    color: tuple[str, str],
    mag_band: str,
    fit_fraction: float,
) -> IsochroneFitterV2:
    b1, b2 = color
    return IsochroneFitterV2(
        iso_file,
        col_mh=1,
        col_age=2,
        col_g=3,
        col_r=4,
        col_mag=5,
        col_mass=6,
        fit_fraction=fit_fraction,
        R_band1=EXTINCTION_R[b1],
        R_band2=EXTINCTION_R[b2],
        R_mag=EXTINCTION_R[mag_band],
        iso_data=iso_data,
        use_bilinear=False,
    )


def range_check(value: float, bounds: tuple[float, float]) -> dict[str, Any]:
    lo, hi = bounds
    return {
        "value": float(value),
        "range": [float(lo), float(hi)],
        "passed": bool(lo <= value <= hi),
    }


def maybe_m5_checks(result, color: tuple[str, str]) -> dict[str, Any]:
    b1, b2 = color
    e_bv = float(result.extinction_gr) / (EXTINCTION_R[b1] - EXTINCTION_R[b2])
    return {
        "note": "Broad M5 sanity ranges; use literature-cited ranges for publication validation.",
        "log_age": range_check(float(result.log_age), M5_SANITY_RANGES["log_age"]),
        "metallicity": range_check(float(result.metallicity), M5_SANITY_RANGES["metallicity"]),
        "distance_mod": range_check(float(result.distance_mod), M5_SANITY_RANGES["distance_mod"]),
        "e_bv": range_check(e_bv, M5_SANITY_RANGES["e_bv"]),
    }


def save_plot(
    output: Path,
    obs_color: np.ndarray,
    obs_mag: np.ndarray,
    iso_color: np.ndarray,
    iso_mag: np.ndarray,
    color_label: str,
    mag_label: str,
) -> None:
    # Delegates to the shared, Qt-free figure module so the runner, CLI and GUI
    # all produce the same paper-quality CMD figure (single source of truth).
    from apex.analysis.cmd.isochrone_plots import save_cmd_plot

    save_cmd_plot(
        output, obs_color, obs_mag, iso_color, iso_mag,
        color_label=color_label, mag_label=mag_label,
        title="APEX Step12 real-data validation",
    )


def save_corner_plot(
    output: Path,
    flat_chain: np.ndarray,
    labels: list[str],
    truths: list[float] | None = None,
) -> bool:
    """Write a corner plot of the posterior. Returns False if ``corner`` absent."""
    from apex.analysis.cmd.isochrone_plots import corner_figure, save_figure

    try:
        import corner  # type: ignore  # noqa: F401
        used = True
    except Exception:
        used = False
    fig = corner_figure(flat_chain, labels, truths=truths)
    save_figure(fig, output, dpi=130)
    return used


def run_mcmc_fit(
    fitter: IsochroneFitterV2,
    obs_color: np.ndarray,
    obs_mag: np.ndarray,
    color_err: np.ndarray,
    mag_err: np.ndarray,
    bounds: FitBounds,
    color: tuple[str, str],
    args: argparse.Namespace,
) -> tuple[Any, Any, np.ndarray, np.ndarray]:
    """Run grid-scan init then MCMC; return (mcmc_result, grid_best, iso_c, iso_m)."""
    from apex.analysis.cmd.isochrone_mcmc import fit_isochrone_mcmc

    b1, b2 = color
    R_color = EXTINCTION_R[b1] - EXTINCTION_R[b2]

    # 1) Grid scan to initialise walkers.
    grid = fitter.fit_grid_scan(
        obs_color, obs_mag, color_err, mag_err,
        bounds=bounds, snr_min=args.fit_snr_min,
        initial_dm=args.initial_dm, initial_egr=args.initial_ecolor,
        local_maxiter=args.local_maxiter,
    )
    gb = grid.best_fit
    init = [gb.log_age, gb.metallicity, gb.distance_mod, gb.extinction_gr]

    priors: dict[str, tuple[float, float]] = {}
    if args.mh_prior is not None:
        priors["mh"] = (float(args.mh_prior[0]), float(args.mh_prior[1]))
    if args.ecolor_prior is not None:
        priors["e_color"] = (float(args.ecolor_prior[0]), float(args.ecolor_prior[1]))

    # Apply the same SNR mask the fitter uses, so MCMC sees the same stars.
    obs_color_err = np.clip(color_err, 0.01, None)
    obs_mag_err = np.clip(mag_err, 0.01, None)
    snr = np.minimum(1.0 / obs_color_err, 1.0 / obs_mag_err)
    mask = (snr > args.fit_snr_min) & np.isfinite(obs_color) & np.isfinite(obs_mag)

    mcmc = fit_isochrone_mcmc(
        obs_color[mask], obs_mag[mask], obs_color_err[mask], obs_mag_err[mask],
        bounds=bounds,
        interp_fn=fitter._interpolate_isochrone,
        imf_fn=fitter._imf_weight,
        n_walkers=args.n_walkers,
        n_burn=args.n_burn,
        n_steps=args.n_steps,
        init_from_gridscan=init,
        priors=priors or None,
        f_bin=args.f_bin,
        f_field=args.f_field,
        sample_f_bin=args.sample_f_bin,
        R_color=R_color,
        seed=args.seed,
        progress=False,
    )

    iso_c, iso_m, _ = fitter._interpolate_isochrone(*mcmc.median_theta)
    return mcmc, gb, iso_c, iso_m


def run_mcmc_fit_multicolor(
    fitter: IsochroneFitterV2,
    mb,  # MultiBandIsochrone
    obs_mc: np.ndarray,
    err_mc: np.ndarray,
    colors: list[tuple[str, str]],
    mag_band: str,
    bounds: FitBounds,
    args: argparse.Namespace,
) -> tuple[Any, Any]:
    """Grid-scan init (on the first colour) then multi-colour MCMC."""
    from apex.analysis.cmd.isochrone_mcmc import fit_isochrone_mcmc

    b1, b2 = colors[0]
    R_color = EXTINCTION_R[b1] - EXTINCTION_R[b2]

    # Grid scan on the first colour to seed walkers (reuses the single-colour
    # fitter built on colour 0).  Falls back to bounds midpoint on failure.
    obs_c = obs_mc[:, 0]
    obs_m = obs_mc[:, -1]
    color_err = err_mc[:, 0]
    mag_err = err_mc[:, -1]
    try:
        grid = fitter.fit_grid_scan(
            obs_c, obs_m, color_err, mag_err,
            bounds=bounds, snr_min=args.fit_snr_min,
            initial_dm=args.initial_dm, initial_egr=args.initial_ecolor,
            local_maxiter=args.local_maxiter,
        )
        gb = grid.best_fit
        init = [gb.log_age, gb.metallicity, gb.distance_mod, gb.extinction_gr]
    except Exception:
        gb = None
        init = None

    priors: dict[str, tuple[float, float]] = {}
    if args.mh_prior is not None:
        priors["mh"] = (float(args.mh_prior[0]), float(args.mh_prior[1]))
    if args.ecolor_prior is not None:
        priors["e_color"] = (float(args.ecolor_prior[0]), float(args.ecolor_prior[1]))

    mcmc = fit_isochrone_mcmc(
        obs_c, obs_m, color_err, mag_err,   # ignored in multi-colour mode
        bounds=bounds,
        interp_fn=fitter._interpolate_isochrone,
        imf_fn=fitter._imf_weight,
        n_walkers=args.n_walkers,
        n_burn=args.n_burn,
        n_steps=args.n_steps,
        init_from_gridscan=init,
        priors=priors or None,
        f_bin=args.f_bin,
        f_field=args.f_field,
        sample_f_bin=args.sample_f_bin,
        R_color=R_color,
        seed=args.seed,
        progress=False,
        colors=colors,
        mag=mag_band,
        multiband_iso=mb,
        obs_multicolor=obs_mc,
        err_multicolor=err_mc,
        err_floor=getattr(args, "err_floor", 0.02),
    )
    return mcmc, gb


def parse_colors(text: str) -> list[tuple[str, str]]:
    """Parse '--colors g-r,r-i' into [('g','r'), ('r','i')]."""
    out = []
    for tok in str(text).split(","):
        tok = tok.strip()
        if not tok:
            continue
        out.append(parse_pair(tok))
    if not out:
        raise argparse.ArgumentTypeError("--colors needs at least one pair")
    return out


def run_validation(args: argparse.Namespace) -> dict[str, Any]:
    result_dir = Path(args.result_dir)
    wide_csv = resolve_wide_csv(result_dir)
    df = pd.read_csv(wide_csv)

    member_meta: dict[str, Any] = {"applied": False, "reason": "not_requested"}
    if getattr(args, "membership", False):
        mmask, member_meta = membership_mask(
            df, nsig=args.membership_nsig, plx_nsig=args.membership_plx_nsig)
        if member_meta.get("applied"):
            df = df.iloc[np.flatnonzero(mmask)].reset_index(drop=True)

    color = parse_pair(args.color) if args.color else None
    mag_band = args.mag
    if color is None or not mag_band:
        inferred_color, inferred_mag = infer_color_mag(df)
        color = color or inferred_color
        mag_band = mag_band or inferred_mag

    iso_file = Path(args.iso_file) if args.iso_file else default_iso_file(color)
    if not iso_file.exists():
        raise FileNotFoundError(f"Isochrone file not found: {iso_file}")

    age_bounds = parse_bounds(args.age_bounds)
    mh_bounds = parse_bounds(args.mh_bounds)
    dm_bounds = parse_bounds(args.dm_bounds)
    ecolor_bounds = parse_bounds(args.ecolor_bounds)

    t0 = time.perf_counter()
    obs_color, obs_mag, color_err, mag_err, fit_df, data_meta = extract_cmd_arrays(
        df,
        color,
        mag_band,
        data_snr_min=args.data_snr_min,
        max_stars=args.max_stars,
        seed=args.seed,
    )
    cache_dir = Path(args.cache_dir) if args.cache_dir else None
    iso_data, iso_cache_meta = load_compact_isochrone(
        iso_file,
        color,
        mag_band,
        age_bounds,
        mh_bounds,
        cache_dir=cache_dir,
    )
    fitter = build_fitter(iso_file, iso_data, color, mag_band, args.fit_fraction)
    bounds = FitBounds(
        log_age=age_bounds,
        metallicity=mh_bounds,
        distance_mod=dm_bounds,
        extinction_gr=ecolor_bounds,
    )
    color_label = f"{color[0]}-{color[1]}"
    b1, b2 = color
    R_color = EXTINCTION_R[b1] - EXTINCTION_R[b2]

    report: dict[str, Any] = {
        "scenario": "cmd_step12_realdata",
        "cluster": args.cluster,
        "method": args.method,
        "result_dir": str(result_dir),
        "wide_csv": str(wide_csv),
        "iso_file": str(iso_file),
        "data": data_meta,
        "membership": member_meta,
        "selection": {
            "color_index": color_label,
            "mag_band": mag_band,
            "data_snr_min": float(args.data_snr_min),
            "fit_snr_min": float(args.fit_snr_min),
            "fit_fraction": float(args.fit_fraction),
        },
        "isochrone_subset": {
            "rows": int(len(iso_data)),
            **iso_cache_meta,
            "age_bounds": list(age_bounds),
            "metallicity_bounds": list(mh_bounds),
            "distance_mod_bounds": list(dm_bounds),
            "extinction_color_bounds": list(ecolor_bounds),
        },
    }

    multicolor = bool(getattr(args, "multicolor", False) or getattr(args, "colors", None))
    if args.method == "mcmc" and multicolor:
        from apex.analysis.cmd.isochrone_mcmc import MultiBandIsochrone

        colors = (parse_colors(args.colors) if getattr(args, "colors", None)
                  else [("g", "r"), ("r", "i")])
        bands = sorted({b for pair in colors for b in pair} | {mag_band})
        obs_mc, err_mc, axis_labels, mc_meta = extract_multicolor_arrays(
            df, colors, mag_band,
            data_snr_min=args.data_snr_min, max_stars=args.max_stars, seed=args.seed,
        )
        # Empirical model-color correction: subtract the measured (data - model) r-i
        # residual line A + B*(color0) from the second color axis, putting the data on
        # the isochrone's color system so multi-color fitting is self-consistent.
        ri_res = getattr(args, "ri_residual", None)
        if ri_res is not None and obs_mc.shape[1] >= 2:
            A, B = float(ri_res[0]), float(ri_res[1])
            obs_mc[:, 1] = obs_mc[:, 1] - (A + B * obs_mc[:, 0])
            report["ri_residual_correction"] = {"intercept": A, "slope": B}
        mb_iso_data, band_col, mb_cache_meta = load_multiband_isochrone(
            iso_file, bands, age_bounds, mh_bounds, cache_dir=cache_dir,
        )
        mb = MultiBandIsochrone(
            mb_iso_data, col_age=2, col_mh=1, col_mass=mb_iso_data.shape[1] - 1,
            band_cols=band_col,
            R_band={b: EXTINCTION_R[b] for b in bands},
        )
        mcmc, grid_best = run_mcmc_fit_multicolor(
            fitter, mb, obs_mc, err_mc, colors, mag_band, bounds, args,
        )
        elapsed = time.perf_counter() - t0
        report["elapsed_sec"] = float(elapsed)
        report["multicolor"] = {
            "colors": [f"{a}-{b}" for a, b in colors],
            "axes": axis_labels,
            "n_stars": int(obs_mc.shape[0]),
            "band_columns_used": mc_meta["columns"],
            "isochrone_rows": int(len(mb_iso_data)),
            **mb_cache_meta,
        }

        def _pc(p):
            if p is None:
                return None
            lo, mid, hi = p
            return {"median": float(mid), "p16": float(lo), "p84": float(hi),
                    "err_minus": float(mid - lo), "err_plus": float(hi - mid)}

        report["fit"] = {
            "log_age": _pc(mcmc.log_age_p),
            "age_gyr": _pc(mcmc.age_gyr_p),
            "metallicity": _pc(mcmc.mh_p),
            "distance_mod": _pc(mcmc.dm_p),
            "distance_pc": _pc(mcmc.distance_pc_p),
            "extinction_color": _pc(mcmc.e_color_p),
            "e_bv": _pc(mcmc.e_bv_p),
            "f_bin": _pc(mcmc.f_bin_p),
            "extinction_label": f"E({color_label})",
            "grid_init": (None if grid_best is None else {
                "log_age": float(grid_best.log_age),
                "metallicity": float(grid_best.metallicity),
                "distance_mod": float(grid_best.distance_mod),
                "extinction_color": float(grid_best.extinction_gr),
            }),
        }
        report["mcmc_diagnostics"] = {
            "n_walkers": int(mcmc.n_walkers),
            "n_burn": int(mcmc.n_burn),
            "n_steps": int(mcmc.n_steps),
            "acceptance_fraction": float(mcmc.acceptance_fraction),
            "autocorr_time": ([float(t) for t in mcmc.autocorr_time]
                              if mcmc.autocorr_time is not None else None),
            "convergence_ok": bool(mcmc.convergence_ok),
            "f_field": float(args.f_field),
            "f_bin": float(args.f_bin),
            "sample_f_bin": bool(args.sample_f_bin),
        }
        report["passed"] = bool(mcmc.convergence_ok)
        if args.output_corner is not None and mcmc.flat_chain is not None:
            used = save_corner_plot(Path(args.output_corner), mcmc.flat_chain,
                                    list(mcmc.labels))
            report["corner_plot"] = str(args.output_corner)
            report["corner_used_corner_lib"] = bool(used)
        # Median isochrone overlay (first colour) for the optional plot.
        iso_color, iso_mag, _ = fitter._interpolate_isochrone(*mcmc.median_theta)
    elif args.method == "mcmc":
        mcmc, grid_best, iso_color, iso_mag = run_mcmc_fit(
            fitter, obs_color, obs_mag, color_err, mag_err,
            bounds, color, args,
        )
        elapsed = time.perf_counter() - t0
        report["elapsed_sec"] = float(elapsed)

        def _pc(p):  # (p16, p50, p84) -> {median, ci_low, ci_high}
            if p is None:
                return None
            lo, mid, hi = p
            return {
                "median": float(mid),
                "p16": float(lo),
                "p84": float(hi),
                "err_minus": float(mid - lo),
                "err_plus": float(hi - mid),
            }

        report["fit"] = {
            "log_age": _pc(mcmc.log_age_p),
            "age_gyr": _pc(mcmc.age_gyr_p),
            "metallicity": _pc(mcmc.mh_p),
            "distance_mod": _pc(mcmc.dm_p),
            "distance_pc": _pc(mcmc.distance_pc_p),
            "extinction_color": _pc(mcmc.e_color_p),
            "e_bv": _pc(mcmc.e_bv_p),
            "f_bin": _pc(mcmc.f_bin_p),
            "extinction_label": f"E({color_label})",
            "grid_init": {
                "log_age": float(grid_best.log_age),
                "metallicity": float(grid_best.metallicity),
                "distance_mod": float(grid_best.distance_mod),
                "extinction_color": float(grid_best.extinction_gr),
            },
        }
        report["mcmc_diagnostics"] = {
            "n_walkers": int(mcmc.n_walkers),
            "n_burn": int(mcmc.n_burn),
            "n_steps": int(mcmc.n_steps),
            "acceptance_fraction": float(mcmc.acceptance_fraction),
            "autocorr_time": (
                [float(t) for t in mcmc.autocorr_time]
                if mcmc.autocorr_time is not None else None
            ),
            "convergence_ok": bool(mcmc.convergence_ok),
            "f_field": float(args.f_field),
            "f_bin": float(args.f_bin),
            "sample_f_bin": bool(args.sample_f_bin),
        }
        report["passed"] = bool(mcmc.convergence_ok)

        if args.output_corner is not None and mcmc.flat_chain is not None:
            used = save_corner_plot(
                Path(args.output_corner), mcmc.flat_chain, list(mcmc.labels),
            )
            report["corner_plot"] = str(args.output_corner)
            report["corner_used_corner_lib"] = bool(used)
    else:
        grid = fitter.fit_grid_scan(
            obs_color,
            obs_mag,
            color_err,
            mag_err,
            bounds=bounds,
            snr_min=args.fit_snr_min,
            initial_dm=args.initial_dm,
            initial_egr=args.initial_ecolor,
            local_maxiter=args.local_maxiter,
        )
        result = grid.best_fit
        result.extinction_label = f"E({color_label})"
        result.color_label = color_label
        result.mag_label = mag_band
        iso_color, iso_mag = fitter.get_best_fit_isochrone(result)
        elapsed = time.perf_counter() - t0
        report["elapsed_sec"] = float(elapsed)

        e_bv = float(result.extinction_gr) / R_color
        distance_pc = float(result.distance_pc)
        report["fit"] = {
            "log_age": float(result.log_age),
            "age_gyr": float(result.age_gyr),
            "metallicity": float(result.metallicity),
            "distance_mod": float(result.distance_mod),
            "distance_pc": distance_pc,
            "extinction_label": result.extinction_label,
            "extinction_color": float(result.extinction_gr),
            "e_bv": e_bv,
            "chi2": float(result.chi2),
            "reduced_chi2": float(result.reduced_chi2),
            "n_stars": int(result.n_stars),
            "converged": bool(result.converged),
            "grid_cells_evaluated": int(grid.n_evaluated),
        }
        if str(args.cluster).strip().lower() == "m5":
            checks = maybe_m5_checks(result, color)
            report["sanity_checks"] = checks
            report["passed"] = bool(
                result.converged
                and result.n_stars >= 20
                and all(v.get("passed", True) for k, v in checks.items() if isinstance(v, dict))
            )
        else:
            report["passed"] = bool(result.converged and result.n_stars >= 20)

    if args.output_plot:
        save_plot(
            Path(args.output_plot),
            obs_color,
            obs_mag,
            iso_color,
            iso_mag,
            color_label,
            mag_band,
        )
        report["plot"] = str(args.output_plot)

    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result_dir", help="APEX result directory, e.g. E:/observed_Analysis/M5/light/result")
    parser.add_argument("--cluster", default="M5")
    parser.add_argument("--iso-file")
    parser.add_argument("--color", help="Color index, e.g. g-r or B-V")
    parser.add_argument("--mag", help="Magnitude band, e.g. g or V")
    parser.add_argument("--age-bounds", nargs=2, type=float, default=[9.80, 10.25])
    parser.add_argument("--mh-bounds", nargs=2, type=float, default=[-1.80, -0.80])
    parser.add_argument("--dm-bounds", nargs=2, type=float, default=[13.80, 14.90])
    parser.add_argument("--ecolor-bounds", nargs=2, type=float, default=[0.00, 0.14])
    parser.add_argument("--initial-dm", type=float, default=14.35)
    parser.add_argument("--initial-ecolor", type=float, default=0.03)
    parser.add_argument("--fit-fraction", type=float, default=0.65)
    parser.add_argument("--fit-snr-min", type=float, default=5.0)
    parser.add_argument("--data-snr-min", type=float, default=10.0)
    parser.add_argument("--max-stars", type=int, default=5000)
    parser.add_argument("--local-maxiter", type=int, default=80)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--output-plot", type=Path)
    parser.add_argument("--output-corner", type=Path,
                        help="Path for the MCMC corner plot (method=mcmc only)")
    parser.add_argument("--cache-dir", default=str(ROOT / "validation" / "reports" / "iso_cache"))
    # --- fit method selection (backward compatible: default = grid) ----------
    parser.add_argument("--method", choices=["grid", "mcmc"], default="grid",
                        help="grid (legacy chi2 grid scan) or mcmc (generative posterior)")
    parser.add_argument("--n-walkers", type=int, default=32)
    parser.add_argument("--n-burn", type=int, default=400)
    parser.add_argument("--n-steps", type=int, default=1500)
    parser.add_argument("--f-bin", type=float, default=0.3,
                        help="Unresolved-binary mixture fraction")
    parser.add_argument("--f-field", type=float, default=0.1,
                        help="Field/outlier mixture fraction")
    parser.add_argument("--sample-f-bin", action="store_true",
                        help="Sample the binary fraction as a 5th nuisance parameter")
    parser.add_argument("--mh-prior", nargs=2, type=float, default=None,
                        metavar=("MU", "SIGMA"),
                        help="Gaussian prior on [M/H] (e.g. spectroscopic [Fe/H])")
    parser.add_argument("--ecolor-prior", nargs=2, type=float, default=None,
                        metavar=("MU", "SIGMA"),
                        help="Gaussian prior on E(color) (e.g. a reddening estimate)")
    parser.add_argument("--membership", action="store_true",
                        help="Pre-filter to probable cluster members using Gaia proper motion "
                             "+ parallax (auto-detected clump). Strongly recommended for real "
                             "clusters: field interlopers otherwise frustrate the generative fit.")
    parser.add_argument("--membership-nsig", type=float, default=3.0,
                        help="Robust-sigma radius in proper-motion space for membership (default 3).")
    parser.add_argument("--membership-plx-nsig", type=float, default=3.0,
                        help="Robust-sigma window in parallax for membership (default 3).")
    parser.add_argument("--err-floor", type=float, default=0.02,
                        help="Systematic/model-error floor (mag) added in quadrature to per-star "
                             "errors in the likelihood. Raise (~0.04-0.05) to tolerate a real "
                             "model-vs-data color-system offset on the second color axis.")
    parser.add_argument("--ri-residual", nargs=2, type=float, default=None, metavar=("A", "B"),
                        help="Empirical model color correction: subtract A + B*(first color) "
                             "from the second color axis (data->isochrone color system).")
    parser.add_argument("--multicolor", action="store_true",
                        help="Multi-colour MCMC over g-r AND r-i (breaks the "
                             "age/[M/H]/reddening degeneracy). Implies method=mcmc.")
    parser.add_argument("--colors", default=None,
                        help="Comma-separated colour list for multi-colour mode, "
                             "e.g. 'g-r,r-i'. Implies --multicolor.")
    args = parser.parse_args(argv)
    if args.colors or args.multicolor:
        args.method = "mcmc"

    report = run_validation(args)
    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
