"""Robustly anchor per-frame PSF fluxes to Step 7 total aperture fluxes."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PSFApertureScale:
    scale: float
    applied: bool
    n_matched: int
    n_candidates: int
    n_used: int
    median_delta_mag_raw: float
    scatter_mag: float
    reason: str


def _bool_series(frame: pd.DataFrame, column: str, default: bool) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(default, index=frame.index, dtype=bool)
    values = frame[column]
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(default).astype(bool)
    text = values.astype(str).str.strip().str.lower()
    mapped = text.map({
        "true": True, "1": True, "yes": True, "on": True,
        "false": False, "0": False, "no": False, "off": False,
    })
    return mapped.fillna(default).astype(bool)


def _numeric(frame: pd.DataFrame, column: str, default: float = np.nan) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(default, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce")


def estimate_psf_aperture_scale(
    psf: pd.DataFrame,
    aperture: pd.DataFrame,
    *,
    min_snr: float = 50.0,
    min_stars: int = 8,
    min_neighbor_fwhm: float = 4.0,
    max_scatter_mag: float = 0.10,
    require_apcorr_candidate: bool = True,
    scale_bounds: tuple[float, float] = (0.5, 2.0),
) -> tuple[PSFApertureScale, pd.DataFrame]:
    """Estimate one multiplicative scale from clean common stars.

    The returned table contains every positive-UID match and records which
    stars were eligible and retained by robust clipping.
    """
    empty = PSFApertureScale(1.0, False, 0, 0, 0, np.nan, np.nan, "no_matches")
    if not isinstance(psf, pd.DataFrame) or not isinstance(aperture, pd.DataFrame):
        return empty, pd.DataFrame()
    if "det_uid" not in psf.columns or "det_uid" not in aperture.columns:
        return empty, pd.DataFrame()

    psf_work = pd.DataFrame({
        "det_uid": _numeric(psf, "det_uid"),
        "flux_psf_e": _numeric(psf, "flux_psf_e"),
        "flags_psf": _numeric(psf, "flags_psf", 0.0),
        "neighbor_dist_fwhm": _numeric(psf, "neighbor_dist_fwhm"),
        "crowding_unreliable_psf": _bool_series(psf, "crowding_unreliable_psf", False),
    })
    aperture_work = pd.DataFrame({
        "det_uid": _numeric(aperture, "det_uid"),
        "flux_aperture_e": _numeric(aperture, "flux_e"),
        "snr_aperture": _numeric(aperture, "snr"),
        "detected_flag": _bool_series(aperture, "detected_flag", True),
        "bad_phot_flag": _bool_series(aperture, "bad_phot_flag", False),
        "off_frame_flag": _bool_series(aperture, "off_frame_flag", False),
        "step4_apcorr_candidate": _bool_series(
            aperture, "step4_apcorr_candidate", True
        ),
    })
    psf_work = psf_work.loc[np.isfinite(psf_work["det_uid"]) & (psf_work["det_uid"] >= 0)].copy()
    aperture_work = aperture_work.loc[
        np.isfinite(aperture_work["det_uid"]) & (aperture_work["det_uid"] >= 0)
    ].copy()
    psf_work["det_uid"] = psf_work["det_uid"].astype(np.int64)
    aperture_work["det_uid"] = aperture_work["det_uid"].astype(np.int64)
    psf_work = psf_work.drop_duplicates("det_uid", keep="first")
    aperture_work = aperture_work.drop_duplicates("det_uid", keep="first")
    matched = psf_work.merge(aperture_work, on="det_uid", how="inner", validate="one_to_one")
    if matched.empty:
        return empty, matched

    eligible = (
        np.isfinite(matched["flux_psf_e"])
        & (matched["flux_psf_e"] > 0)
        & np.isfinite(matched["flux_aperture_e"])
        & (matched["flux_aperture_e"] > 0)
        & np.isfinite(matched["snr_aperture"])
        & (matched["snr_aperture"] >= max(0.0, float(min_snr)))
        & (matched["flags_psf"].fillna(-1) == 0)
        & ~matched["crowding_unreliable_psf"]
        & matched["detected_flag"]
        & ~matched["bad_phot_flag"]
        & ~matched["off_frame_flag"]
    )
    if min_neighbor_fwhm > 0:
        eligible &= (
            np.isfinite(matched["neighbor_dist_fwhm"])
            & (matched["neighbor_dist_fwhm"] >= float(min_neighbor_fwhm))
        )
    if require_apcorr_candidate:
        eligible &= matched["step4_apcorr_candidate"]

    matched["eligible"] = eligible
    matched["delta_mag_raw"] = np.nan
    matched.loc[eligible, "delta_mag_raw"] = -2.5 * np.log10(
        matched.loc[eligible, "flux_psf_e"] / matched.loc[eligible, "flux_aperture_e"]
    )
    candidate_index = np.flatnonzero(eligible.to_numpy(dtype=bool))
    matched["used"] = False
    min_stars = max(3, int(min_stars))
    if candidate_index.size < min_stars:
        result = PSFApertureScale(
            1.0, False, len(matched), int(candidate_index.size), 0,
            np.nan, np.nan, "too_few_candidates",
        )
        return result, matched

    values = matched.loc[eligible, "delta_mag_raw"].to_numpy(dtype=float)
    keep = np.isfinite(values)
    for _ in range(4):
        current = values[keep]
        if current.size < min_stars:
            break
        median = float(np.median(current))
        scatter = float(1.4826 * np.median(np.abs(current - median)))
        clip_radius = max(0.03, 3.0 * scatter)
        next_keep = np.isfinite(values) & (np.abs(values - median) <= clip_radius)
        if np.array_equal(next_keep, keep):
            break
        keep = next_keep

    used_values = values[keep]
    used_index = candidate_index[keep]
    matched.loc[used_index, "used"] = True
    if used_values.size < min_stars:
        result = PSFApertureScale(
            1.0, False, len(matched), int(candidate_index.size), int(used_values.size),
            np.nan, np.nan, "too_few_after_clipping",
        )
        return result, matched

    median_delta = float(np.median(used_values))
    scatter = float(1.4826 * np.median(np.abs(used_values - median_delta)))
    scale = float(10.0 ** (median_delta / 2.5))
    lower, upper = sorted((float(scale_bounds[0]), float(scale_bounds[1])))
    reason = "ok"
    applied = True
    if not np.isfinite(scatter) or scatter > max(0.0, float(max_scatter_mag)):
        applied = False
        reason = "scatter_exceeds_limit"
    elif not np.isfinite(scale) or scale < lower or scale > upper:
        applied = False
        reason = "scale_out_of_bounds"
    if not applied:
        scale = 1.0

    result = PSFApertureScale(
        scale=scale,
        applied=applied,
        n_matched=len(matched),
        n_candidates=int(candidate_index.size),
        n_used=int(used_values.size),
        median_delta_mag_raw=median_delta,
        scatter_mag=scatter,
        reason=reason,
    )
    return result, matched


def apply_psf_aperture_scale(
    psf: pd.DataFrame,
    result: PSFApertureScale,
    *,
    zeropoint: float,
    exptime: float,
) -> pd.DataFrame:
    """Preserve raw columns and apply an accepted scale to PSF flux output."""
    output = psf.copy()
    flux = _numeric(output, "flux_psf_e")
    flux_error = _numeric(output, "flux_psf_err_e")
    output["flux_psf_raw_e"] = flux
    output["flux_psf_err_raw_e"] = flux_error
    output["mag_psf_raw"] = _numeric(output, "mag_psf")
    scale = float(result.scale) if result.applied else 1.0
    output["flux_psf_e"] = flux * scale
    output["flux_psf_err_e"] = flux_error * scale
    corrected_flux = output["flux_psf_e"].to_numpy(dtype=float)
    valid_mag = np.isfinite(corrected_flux) & (corrected_flux > 0) & (float(exptime) > 0)
    corrected_mag = np.full(len(output), np.nan, dtype=float)
    corrected_mag[valid_mag] = float(zeropoint) - 2.5 * np.log10(
        corrected_flux[valid_mag] / float(exptime)
    )
    original_mag = _numeric(output, "mag_psf").to_numpy(dtype=float)
    output["mag_psf"] = np.where(np.isfinite(original_mag), corrected_mag, np.nan)
    output["psf_aperture_scale"] = scale
    output["psf_aperture_scale_applied"] = bool(result.applied)
    output["psf_aperture_scale_n"] = int(result.n_used)
    output["psf_aperture_scale_scatter_mag"] = float(result.scatter_mag)
    output["psf_aperture_scale_reason"] = result.reason
    return output
