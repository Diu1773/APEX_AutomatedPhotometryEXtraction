"""Frame-level quality control for Step 4 source detection.

The source detector owns per-frame measurements.  This module turns those
measurements into a conservative PASS/REVIEW/FAIL decision that can be reused
by the GUI, downstream gates, and paper-quality validation plots.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd


PASS = "PASS"
REVIEW = "REVIEW"
FAIL = "FAIL"

_SEEING_AIRMASS_EXPONENT = 0.6  # Kolmogorov seeing: FWHM proportional to X^(3/5).


@dataclass(frozen=True)
class FrameQCThresholds:
    """Conservative default thresholds for automatic frame QC."""

    fwhm_z_review: float = 3.0
    fwhm_z_fail: float = 4.5
    fwhm_model_ratio_review: float = 1.30
    fwhm_model_ratio_fail: float = 1.60
    sky_z_review: float = 4.0
    sky_z_fail: float = 5.5
    sky_noise_ratio_review: float = 2.5
    sky_noise_ratio_fail: float = 6.0
    nsrc_z_review: float = 3.0
    nsrc_z_fail: float = 4.5
    elong_review: float = 1.22
    elong_fail: float = 1.35
    quality_score_review: float = 70.0
    quality_score_fail: float = 45.0
    sat_rate_review: float = 0.08
    sat_rate_fail: float = 0.30
    min_anchor_review: int = 3
    min_psf_seed_review: int = 10
    min_epsf_review: int = 3


def _get_float(params: Any, name: str, default: float) -> float:
    try:
        value = float(getattr(params, name, default))
        return value if np.isfinite(value) else float(default)
    except Exception:
        return float(default)


def _get_int(params: Any, name: str, default: int) -> int:
    try:
        value = int(getattr(params, name, default))
        return value if value >= 0 else int(default)
    except Exception:
        return int(default)


def _num(df: pd.DataFrame, name: str, default: float = np.nan) -> pd.Series:
    if name not in df.columns:
        return pd.Series(np.full(len(df), default, dtype=float), index=df.index)
    return pd.to_numeric(df[name], errors="coerce")


def _robust_center_scale(values: np.ndarray) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return np.nan, np.nan
    med = float(np.nanmedian(finite))
    mad = float(np.nanmedian(np.abs(finite - med)))
    if np.isfinite(mad) and mad > 0:
        return med, mad / 0.6745
    std = float(np.nanstd(finite))
    if np.isfinite(std) and std > 0:
        return med, std
    return med, np.nan


def _robust_z(values: np.ndarray) -> tuple[np.ndarray, float, float]:
    med, scale = _robust_center_scale(values)
    z = np.full(len(values), np.nan, dtype=float)
    if np.isfinite(med) and np.isfinite(scale) and scale > 0:
        z = (values - med) / scale
    return z, med, scale


def _finite_median(values: np.ndarray, default: float = np.nan) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return float(default)
    return float(np.nanmedian(finite))


def _filter_groups(df: pd.DataFrame):
    if "filter" not in df.columns:
        yield "", df.index
        return
    filters = df["filter"].fillna("").astype(str)
    for _, idx in filters.groupby(filters, sort=False).groups.items():
        yield "", pd.Index(idx)


def _seeing_model_for_group(fwhm: np.ndarray, airmass: np.ndarray) -> np.ndarray:
    model = np.full(len(fwhm), np.nan, dtype=float)
    good = np.isfinite(fwhm) & (fwhm > 0) & np.isfinite(airmass) & (airmass > 0)
    if int(good.sum()) >= 3 and float(np.nanmax(airmass[good]) - np.nanmin(airmass[good])) > 0.02:
        fwhm0 = _finite_median(fwhm[good] / np.power(airmass[good], _SEEING_AIRMASS_EXPONENT))
        if np.isfinite(fwhm0) and fwhm0 > 0:
            model = fwhm0 * np.power(np.clip(airmass, 1e-6, None), _SEEING_AIRMASS_EXPONENT)
    if not np.isfinite(model).any():
        med = _finite_median(fwhm[np.isfinite(fwhm) & (fwhm > 0)])
        if np.isfinite(med) and med > 0:
            model[:] = med
    return model


def _group_envelope(values: np.ndarray, high_z: float | None = None, low_z: float | None = None):
    _, med, scale = _robust_z(values)
    high = np.full(len(values), np.nan, dtype=float)
    low = np.full(len(values), np.nan, dtype=float)
    if np.isfinite(med):
        if high_z is not None and np.isfinite(scale) and scale > 0:
            high[:] = med + high_z * scale
        if low_z is not None and np.isfinite(scale) and scale > 0:
            low[:] = med - low_z * scale
    return med, high, low


def evaluate_frame_qc(
    frame_df: pd.DataFrame,
    params: Any = None,
    thresholds: FrameQCThresholds | None = None,
) -> pd.DataFrame:
    """Return ``frame_df`` with automatic frame QC columns added.

    The decision policy is intentionally conservative.  ``FAIL`` is reserved
    for frames that are very likely to poison later WCS/photometry steps;
    ambiguous frames are left as ``REVIEW`` so the user can keep them.
    """

    if frame_df is None or frame_df.empty:
        return pd.DataFrame() if frame_df is None else frame_df.copy()

    thr = thresholds or FrameQCThresholds(
        elong_fail=max(1.0, _get_float(params, "fwhm_elong_max", FrameQCThresholds.elong_fail)),
    )
    out = frame_df.copy()

    fwhm = _num(out, "fwhm_med").to_numpy(float)
    sky = _num(out, "sky_med").to_numpy(float)
    sky_sigma = _num(out, "sky_sigma").to_numpy(float)
    nsrc = _num(out, "n_sources", 0.0).to_numpy(float)
    elong = _num(out, "elong_med").to_numpy(float)
    quality = _num(out, "quality_score_median").to_numpy(float)
    sat = _num(out, "sat_star_count", 0.0).to_numpy(float)
    n_anchor = _num(out, "n_anchor_candidates").to_numpy(float)
    n_psf_seed = _num(out, "n_psf_seed_candidates").to_numpy(float)
    n_epsf = _num(out, "n_epsf_candidates").to_numpy(float)
    gain = _num(out, "gain_e_per_adu").to_numpy(float)
    rdnoise = _num(out, "rdnoise_e").to_numpy(float)
    airmass = _num(out, "airmass").to_numpy(float)

    n = len(out)
    fwhm_z = np.full(n, np.nan, dtype=float)
    sky_z = np.full(n, np.nan, dtype=float)
    sky_sigma_z = np.full(n, np.nan, dtype=float)
    nsrc_z = np.full(n, np.nan, dtype=float)
    fwhm_model = np.full(n, np.nan, dtype=float)
    fwhm_high_cut = np.full(n, np.nan, dtype=float)
    sky_high_cut = np.full(n, np.nan, dtype=float)
    nsrc_trend = np.full(n, np.nan, dtype=float)
    nsrc_low_cut = np.full(n, np.nan, dtype=float)

    for _, idx in _filter_groups(out):
        loc = out.index.get_indexer(idx)
        loc = loc[loc >= 0]
        if loc.size == 0:
            continue
        fwhm_z[loc], _, _ = _robust_z(fwhm[loc])
        sky_z[loc], _, _ = _robust_z(sky[loc])
        sky_sigma_z[loc], _, _ = _robust_z(sky_sigma[loc])
        nsrc_z[loc], _, _ = _robust_z(nsrc[loc])
        fwhm_model[loc] = _seeing_model_for_group(fwhm[loc], airmass[loc])
        _, f_high, _ = _group_envelope(fwhm[loc], high_z=thr.fwhm_z_review)
        _, sky_high, _ = _group_envelope(sky[loc], high_z=thr.sky_z_review)
        n_med, _, n_low = _group_envelope(nsrc[loc], low_z=thr.nsrc_z_review)
        fwhm_high_cut[loc] = f_high
        sky_high_cut[loc] = sky_high
        nsrc_trend[loc] = n_med
        nsrc_low_cut[loc] = n_low

    with np.errstate(divide="ignore", invalid="ignore"):
        fwhm_model_ratio = fwhm / fwhm_model
        sky_e = sky * gain
        sky_sigma_e = sky_sigma * gain
        sky_sigma_expected_e = np.sqrt(np.clip(sky_e, 0.0, None) + np.square(rdnoise))
        sky_noise_ratio = sky_sigma_e / sky_sigma_expected_e
        sat_rate = sat / np.maximum(nsrc, 1.0)

    fwhm_px_min = _get_float(params, "fwhm_px_min", 0.0)
    fwhm_px_max = _get_float(params, "fwhm_px_max", 0.0)
    min_sources_abs = max(_get_int(params, "gate_nsrc_min", 0), 1)
    critical_sources_abs = max(min_sources_abs, 3)

    statuses: list[str] = []
    reasons: list[str] = []
    scores: list[float] = []

    for i in range(n):
        fail: list[str] = []
        review: list[str] = []
        score = 100.0

        if not np.isfinite(nsrc[i]) or nsrc[i] < critical_sources_abs:
            fail.append("no_or_few_sources")
            score -= 55.0
        if not np.isfinite(fwhm[i]) or fwhm[i] <= 0:
            fail.append("fwhm_missing")
            score -= 55.0

        if np.isfinite(elong[i]):
            if elong[i] > thr.elong_fail:
                fail.append("high_elongation")
                score -= 35.0
            elif elong[i] > thr.elong_review:
                review.append("elongation_warning")
                score -= 15.0

        if np.isfinite(fwhm[i]) and fwhm[i] > 0:
            if fwhm_px_min > 0 and fwhm[i] < fwhm_px_min:
                review.append("fwhm_below_config")
                score -= 10.0
            if fwhm_px_max > 0 and fwhm[i] > fwhm_px_max:
                review.append("fwhm_above_config")
                score -= 15.0

        high_fwhm_z = np.isfinite(fwhm_z[i]) and fwhm_z[i] > thr.fwhm_z_fail
        high_fwhm_model = np.isfinite(fwhm_model_ratio[i]) and fwhm_model_ratio[i] > thr.fwhm_model_ratio_fail
        if high_fwhm_z or high_fwhm_model:
            fail.append("high_fwhm")
            score -= 35.0
        elif (
            (np.isfinite(fwhm_z[i]) and fwhm_z[i] > thr.fwhm_z_review)
            or (np.isfinite(fwhm_model_ratio[i]) and fwhm_model_ratio[i] > thr.fwhm_model_ratio_review)
        ):
            review.append("fwhm_warning")
            score -= 15.0

        low_nsrc_ratio_fail = (
            np.isfinite(nsrc[i])
            and np.isfinite(nsrc_trend[i])
            and nsrc_trend[i] >= 20
            and nsrc[i] < 0.25 * nsrc_trend[i]
        )
        low_nsrc_fail = (np.isfinite(nsrc_z[i]) and nsrc_z[i] < -thr.nsrc_z_fail) or low_nsrc_ratio_fail
        low_nsrc_review = np.isfinite(nsrc_z[i]) and nsrc_z[i] < -thr.nsrc_z_review
        if low_nsrc_fail:
            fail.append("low_n_sources")
            score -= 30.0
        elif low_nsrc_review:
            review.append("n_sources_warning")
            score -= 12.0

        sky_fail = (
            np.isfinite(sky_z[i])
            and sky_z[i] > thr.sky_z_fail
            and np.isfinite(sky_noise_ratio[i])
            and sky_noise_ratio[i] > thr.sky_noise_ratio_fail
        )
        sky_review = (
            (np.isfinite(sky_z[i]) and sky_z[i] > thr.sky_z_review)
            or (np.isfinite(sky_sigma_z[i]) and sky_sigma_z[i] > thr.sky_z_review)
            or (np.isfinite(sky_noise_ratio[i]) and sky_noise_ratio[i] > thr.sky_noise_ratio_review)
        )
        if sky_fail and low_nsrc_review:
            fail.append("high_sky_noise")
            score -= 30.0
        elif sky_review:
            review.append("sky_warning")
            score -= 10.0

        if np.isfinite(quality[i]):
            if quality[i] < thr.quality_score_fail:
                fail.append("low_source_quality")
                score -= 25.0
            elif quality[i] < thr.quality_score_review:
                review.append("source_quality_warning")
                score -= 10.0

        if np.isfinite(sat_rate[i]):
            if sat_rate[i] > thr.sat_rate_fail and nsrc[i] >= 10:
                fail.append("many_saturated_sources")
                score -= 25.0
            elif sat_rate[i] > thr.sat_rate_review:
                review.append("saturation_warning")
                score -= 10.0

        if np.isfinite(n_anchor[i]) and n_anchor[i] < thr.min_anchor_review:
            review.append("few_anchor_candidates")
            score -= 8.0
        if np.isfinite(n_psf_seed[i]) and n_psf_seed[i] < thr.min_psf_seed_review:
            review.append("few_psf_seed_candidates")
            score -= 8.0
        if np.isfinite(n_epsf[i]) and n_epsf[i] < thr.min_epsf_review:
            review.append("few_epsf_candidates")
            score -= 8.0

        if fail:
            statuses.append(FAIL)
            reason = fail + review
        elif review:
            statuses.append(REVIEW)
            reason = review
        else:
            statuses.append(PASS)
            reason = []

        reasons.append(",".join(dict.fromkeys(reason)))
        scores.append(float(np.clip(score, 0.0, 100.0)))

    out["fwhm_z"] = fwhm_z
    out["sky_z"] = sky_z
    out["sky_sigma_z"] = sky_sigma_z
    out["nsrc_z"] = nsrc_z
    out["fwhm_model_px"] = fwhm_model
    out["fwhm_model_ratio"] = fwhm_model_ratio
    out["fwhm_high_cut_px"] = fwhm_high_cut
    out["sky_high_cut_adu"] = sky_high_cut
    out["sky_e"] = sky_e
    out["sky_sigma_e"] = sky_sigma_e
    out["sky_sigma_expected_e"] = sky_sigma_expected_e
    out["sky_noise_ratio"] = sky_noise_ratio
    out["n_sources_trend"] = nsrc_trend
    out["n_sources_low_cut"] = nsrc_low_cut
    out["sat_star_rate"] = sat_rate
    out["qc_status"] = statuses
    out["qc_score"] = scores
    out["qc_reasons"] = reasons
    return out


def summarize_frame_qc(frame_df: pd.DataFrame) -> dict[str, int]:
    """Return PASS/REVIEW/FAIL counts for a QC-evaluated frame table."""

    counts = {PASS: 0, REVIEW: 0, FAIL: 0}
    if frame_df is None or frame_df.empty or "qc_status" not in frame_df.columns:
        return counts
    series = frame_df["qc_status"].fillna("").astype(str).str.upper()
    for key in counts:
        counts[key] = int((series == key).sum())
    return counts
