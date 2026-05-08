from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

try:
    from scipy.spatial import cKDTree as KDTree
except Exception:  # pragma: no cover - scipy is a runtime dependency for APEX
    KDTree = None


def _get_float(params: Any, name: str, default: float) -> float:
    try:
        value = getattr(params, name, default)
        value = float(value)
        return value if np.isfinite(value) else float(default)
    except Exception:
        return float(default)


def _numeric_col(df: pd.DataFrame, name: str, default: float = np.nan) -> np.ndarray:
    if name not in df.columns:
        return np.full(len(df), default, dtype=float)
    return pd.to_numeric(df[name], errors="coerce").to_numpy(float)


def _bool_status(df: pd.DataFrame, name: str, value: str) -> np.ndarray:
    if name not in df.columns:
        return np.zeros(len(df), dtype=bool)
    return df[name].astype(str).str.strip().str.lower().eq(value).to_numpy(bool)


def _nearest_neighbor_px(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    out = np.full(len(x), np.nan, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    idx_valid = np.flatnonzero(valid)
    if len(idx_valid) < 2:
        return out
    xy = np.column_stack([x[valid], y[valid]])
    if KDTree is not None:
        dists, _ = KDTree(xy).query(xy, k=2)
        out[idx_valid] = dists[:, 1]
        return out

    # Slow fallback for stripped-down environments.
    delta = xy[:, None, :] - xy[None, :, :]
    dist = np.hypot(delta[:, :, 0], delta[:, :, 1])
    np.fill_diagonal(dist, np.inf)
    out[idx_valid] = np.nanmin(dist, axis=1)
    return out


def _edge_margin_px(x: np.ndarray, y: np.ndarray, image_shape: tuple[int, int] | None) -> np.ndarray:
    if image_shape is None or len(image_shape) < 2:
        return np.full(len(x), np.nan, dtype=float)
    h, w = int(image_shape[0]), int(image_shape[1])
    return np.minimum.reduce([x, y, (w - 1) - x, (h - 1) - y])


def _flux_percentile(flux: np.ndarray) -> np.ndarray:
    out = np.full(len(flux), np.nan, dtype=float)
    valid = np.isfinite(flux) & (flux > 0)
    idx = np.flatnonzero(valid)
    if len(idx) == 0:
        return out
    order = np.argsort(flux[idx], kind="mergesort")
    ranks = np.empty(len(idx), dtype=float)
    if len(idx) == 1:
        ranks[order] = 100.0
    else:
        ranks[order] = np.linspace(0.0, 100.0, len(idx))
    out[idx] = ranks
    return out


def _flag_strings(flag_pairs: list[tuple[str, np.ndarray]], n: int) -> list[str]:
    out: list[str] = []
    for i in range(n):
        labels = [name for name, mask in flag_pairs if bool(mask[i])]
        out.append(";".join(labels))
    return out


def compute_source_quality(
    sources: pd.DataFrame,
    *,
    frame_fwhm_px: float | None = None,
    image_shape: tuple[int, int] | None = None,
    params: Any = None,
) -> pd.DataFrame:
    """Add source-quality and role-candidate columns to a Step4 detection table.

    Step4 owns the frame-local measurements. Later steps can tighten these
    first-pass labels for registration, aperture correction, and EPSF building.
    """
    out = sources.copy()
    n = len(out)
    if n == 0:
        return out

    x = _numeric_col(out, "x")
    y = _numeric_col(out, "y")
    peak = _numeric_col(out, "peak_adu")
    dao_flux = _numeric_col(out, "dao_flux")
    dao_peak = _numeric_col(out, "dao_peak")
    fwhm_src = _numeric_col(out, "fwhm_px")
    elong = _numeric_col(out, "elongation")
    sharp = _numeric_col(out, "sharpness")
    roundness = _numeric_col(out, "roundness")
    if not np.isfinite(roundness).any():
        r1 = np.abs(_numeric_col(out, "roundness1"))
        r2 = np.abs(_numeric_col(out, "roundness2"))
        roundness = np.nanmax(np.vstack([r1, r2]), axis=0)

    if frame_fwhm_px is None or not np.isfinite(frame_fwhm_px) or frame_fwhm_px <= 0:
        ok_fwhm = fwhm_src[np.isfinite(fwhm_src) & (fwhm_src > 0)]
        frame_fwhm = float(np.nanmedian(ok_fwhm)) if len(ok_fwhm) else _get_float(params, "fwhm_seed_px", 6.0)
    else:
        frame_fwhm = float(frame_fwhm_px)
    frame_fwhm = max(frame_fwhm, 1.0)

    sat_adu = _get_float(params, "saturation_adu", 60000.0)
    datamax_adu = _get_float(params, "datamax_adu", sat_adu)
    elong_max = _get_float(params, "ref_cat_max_elong", 1.5)
    round_abs_max = _get_float(params, "ref_cat_max_abs_round", 0.4)
    sharp_min = _get_float(params, "ref_cat_sharp_min", 0.2)
    sharp_max = _get_float(params, "ref_cat_sharp_max", 1.0)
    fwhm_ratio_lo = _get_float(params, "source_quality_fwhm_ratio_lo", 0.6)
    fwhm_ratio_hi = _get_float(params, "source_quality_fwhm_ratio_hi", 1.6)
    anchor_iso_mult = _get_float(params, "source_quality_anchor_neighbor_fwhm_mult", 2.0)
    apcorr_iso_mult = _get_float(params, "apcorr_isolation_factor", 2.5)
    epsf_iso_mult = _get_float(params, "psf_isolation_fwhm_mult", 3.0)
    anchor_min_pct = _get_float(params, "source_quality_anchor_flux_pct", 60.0)
    apcorr_min_pct = _get_float(params, "source_quality_apcorr_flux_pct", 60.0)
    epsf_pct_lo = _get_float(params, "psf_flux_percentile_lo", 75.0)
    epsf_pct_hi = _get_float(params, "psf_flux_percentile_hi", 95.0)
    psf_seed_min_pct = _get_float(params, "source_quality_psf_seed_flux_pct", 30.0)
    edge_min_mult = _get_float(params, "source_quality_edge_fwhm_mult", 1.0)

    flux = np.where(np.isfinite(dao_flux) & (dao_flux > 0), dao_flux, np.nan)
    flux = np.where(np.isfinite(flux), flux, dao_peak)
    flux = np.where(np.isfinite(flux), flux, peak)
    flux_pct = _flux_percentile(flux)
    peak_fraction = peak / sat_adu if sat_adu > 0 else np.full(n, np.nan)
    nn_px = _nearest_neighbor_px(x, y)
    nn_fwhm = nn_px / frame_fwhm
    edge_margin = _edge_margin_px(x, y, image_shape)
    measured_fwhm = np.isfinite(fwhm_src) & (fwhm_src > 0)
    # Most Step4 runs measure FWHM only for QC candidates. For unmeasured
    # sources, use the frame median for role selection instead of rejecting
    # otherwise good isolated stars.
    fwhm_ratio = np.where(measured_fwhm, fwhm_src / frame_fwhm, 1.0)

    invalid_xy = ~(np.isfinite(x) & np.isfinite(y))
    near_edge = _bool_status(out, "fwhm_status", "edge")
    if image_shape is not None:
        near_edge = near_edge | (np.isfinite(edge_margin) & (edge_margin < edge_min_mult * frame_fwhm))
    saturated = _bool_status(out, "fwhm_status", "saturated")
    saturated = saturated | (np.isfinite(peak) & (peak >= sat_adu))
    nonlinear = np.isfinite(peak) & (datamax_adu > 0) & (peak >= datamax_adu)
    fwhm_missing = _bool_status(out, "fwhm_status", "not_measured") | _bool_status(out, "fwhm_status", "failed")
    fwhm_outlier = (
        measured_fwhm
        & np.isfinite(fwhm_ratio)
        & ((fwhm_ratio < fwhm_ratio_lo) | (fwhm_ratio > fwhm_ratio_hi))
    )
    bad_elong = np.isfinite(elong) & (elong_max > 0) & (elong > elong_max)
    bad_round = np.isfinite(roundness) & (round_abs_max > 0) & (np.abs(roundness) > round_abs_max)
    bad_sharp = np.zeros(n, dtype=bool)
    if np.isfinite(sharp_min):
        bad_sharp = bad_sharp | (np.isfinite(sharp) & (sharp < sharp_min))
    if np.isfinite(sharp_max) and sharp_max > 0:
        bad_sharp = bad_sharp | (np.isfinite(sharp) & (sharp > sharp_max))
    faint = np.isfinite(flux_pct) & (flux_pct < psf_seed_min_pct)
    close_neighbor = np.isfinite(nn_fwhm) & (nn_fwhm < anchor_iso_mult)

    base_ok = (
        ~invalid_xy
        & ~near_edge
        & ~saturated
        & ~nonlinear
        & ~fwhm_outlier
        & ~bad_elong
        & ~bad_round
        & ~bad_sharp
    )
    bright_enough = np.isfinite(flux_pct)
    peak_safe = ~np.isfinite(peak_fraction) | (peak_fraction < 0.8)

    anchor_candidate = (
        base_ok
        & bright_enough
        & (flux_pct >= anchor_min_pct)
        & (nn_fwhm >= anchor_iso_mult)
        & peak_safe
    )
    apcorr_candidate = (
        base_ok
        & bright_enough
        & (flux_pct >= apcorr_min_pct)
        & (nn_fwhm >= apcorr_iso_mult)
        & peak_safe
    )

    epsf_sharp_lo = _get_float(params, "psf_epsf_sharp_lo", 0.3)
    epsf_sharp_hi = _get_float(params, "psf_epsf_sharp_hi", 0.8)
    epsf_round_max = _get_float(params, "psf_epsf_round_abs_max", 0.5)
    epsf_elong_max = _get_float(params, "psf_epsf_elong_max", 1.3)
    epsf_shape_ok = np.ones(n, dtype=bool)
    epsf_shape_ok &= ~np.isfinite(sharp) | ((sharp >= epsf_sharp_lo) & (sharp <= epsf_sharp_hi))
    epsf_shape_ok &= ~np.isfinite(roundness) | (np.abs(roundness) <= epsf_round_max)
    epsf_shape_ok &= ~np.isfinite(elong) | (elong <= epsf_elong_max)
    epsf_candidate = (
        base_ok
        & epsf_shape_ok
        & bright_enough
        & (flux_pct >= epsf_pct_lo)
        & (flux_pct <= epsf_pct_hi)
        & (nn_fwhm >= epsf_iso_mult)
        & peak_safe
    )

    psf_seed_candidate = (
        base_ok
        & bright_enough
        & (flux_pct >= psf_seed_min_pct)
        & (~np.isfinite(peak_fraction) | (peak_fraction < 0.95))
    )

    score = np.full(n, 100.0, dtype=float)
    penalties = [
        (invalid_xy, 100.0),
        (saturated | nonlinear, 100.0),
        (near_edge, 35.0),
        (fwhm_missing, 5.0),
        (fwhm_outlier, 25.0),
        (bad_elong, 20.0),
        (bad_round, 20.0),
        (bad_sharp, 20.0),
        (close_neighbor, 20.0),
        (faint, 15.0),
    ]
    for mask, penalty in penalties:
        score[mask] -= penalty
    score = np.clip(score, 0.0, 100.0)

    flag_pairs = [
        ("invalid_xy", invalid_xy),
        ("edge", near_edge),
        ("saturated", saturated),
        ("nonlinear", nonlinear),
        ("fwhm_missing", fwhm_missing),
        ("fwhm_outlier", fwhm_outlier),
        ("elongation", bad_elong),
        ("roundness", bad_round),
        ("sharpness", bad_sharp),
        ("close_neighbor", close_neighbor),
        ("faint", faint),
    ]

    out["nearest_neighbor_px"] = nn_px
    out["nearest_neighbor_fwhm"] = nn_fwhm
    out["edge_margin_px"] = edge_margin
    out["fwhm_ratio_to_frame"] = fwhm_ratio
    out["flux_for_quality"] = flux
    out["flux_percentile"] = flux_pct
    out["peak_fraction_to_sat"] = peak_fraction
    out["quality_score"] = score
    out["quality_flags"] = _flag_strings(flag_pairs, n)
    out["anchor_candidate"] = anchor_candidate
    out["apcorr_candidate"] = apcorr_candidate
    out["epsf_candidate"] = epsf_candidate
    out["psf_seed_candidate"] = psf_seed_candidate
    return out
