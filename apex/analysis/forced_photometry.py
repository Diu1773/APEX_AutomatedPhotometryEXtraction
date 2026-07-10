"""Qt-free forced aperture photometry (headless port of ForcedPhotWorker.run()).

This module is a VERBATIM RELOCATION of the compute body of
``apex.gui.workflow.step7_forced_aperture_phot.ForcedPhotWorker`` (its ``run`` /
``_run_impl`` and all the per-frame measurement helpers it calls) with only the
Qt coupling substituted by optional callbacks. The numerical algorithm, parameter
lookups, rounding, and file writes are unchanged.

LAYER RULE: this module must NOT import apex.gui or PyQt5.
"""

from __future__ import annotations

import hashlib  # noqa: F401  (kept for parity with the GUI module imports)
import json
import math
import time
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from threading import Lock
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.wcs import WCS
from scipy.spatial import cKDTree as KDTree

from apex.utils.step_paths import (
    step2_cropped_dir,
    step4_dir,
    step5_wcs_dir,  # noqa: F401
    step6_refbuild_dir,
    step7_forced_phot_dir,
    crop_is_active,
)
from apex.utils.qc_utils import filter_files_by_qc  # noqa: F401
from apex.utils.photometry_utils import (
    phot_vectorized,
    refine_local_centroid,
)
from apex.utils.cache_utils import astap_wcs_candidates, parse_astap_wcs_file
from apex.utils.constants import (
    get_parallel_workers,
    MAG_ERR_COEFF,
    MAD_TO_SIGMA,
    INSTRUMENTAL_ZMAG,
    EXPTIME_HEADER_KEYS,
)
from apex.utils.noise_params import resolve_effective_noise_params

_GC_N_STEPS = 14   # number of radii in the growth curve


# Qt-free relocation of run_control.format_duration (a GUI module that imports
# PyQt5); used only for the orchestration log strings below. Copied byte-for-byte.
def _fmt_duration(seconds: float) -> str:
    """Format elapsed\\remaining seconds as MM:SS or H:MM:SS."""
    try:
        value = float(seconds)
    except Exception:
        value = 0.0
    if not math.isfinite(value):
        value = 0.0
    total = int(max(0, round(value)))
    hours, rem = divmod(total, 3600)
    minutes, secs = divmod(rem, 60)
    if hours > 0:
        return f"{hours:d}:{minutes:02d}:{secs:02d}"
    return f"{minutes:02d}:{secs:02d}"


# ── Scalar helpers ─────────────────────────────────────────────────────────────

def _to_float(val, default: float) -> float:
    try:
        v = float(val) if val is not None else float(default)
        return v if np.isfinite(v) else float(default)
    except Exception:
        return float(default)


def _to_int(val, default: int) -> int:
    try:
        return int(float(val)) if val is not None else int(default)
    except Exception:
        return int(default)


def _exptime_from_header(header, default: float = 1.0) -> tuple[float, bool]:
    """Exposure time (seconds) from a FITS header, for count-rate magnitudes.

    Returns ``(exptime, found)``. ``found`` is False when no usable EXPTIME-like
    key is present and ``default`` (1.0 s) is used; the caller should warn,
    because a frame silently defaulting to 1.0 s inside a mixed-exposure set
    would land on the wrong magnitude scale.
    """
    if header is None:
        return float(default), False
    for key in EXPTIME_HEADER_KEYS:
        if key in header:
            try:
                val = float(header[key])
            except (TypeError, ValueError):
                continue
            if np.isfinite(val) and val > 0:
                return val, True
    return float(default), False


def _catalog_series(df: pd.DataFrame, col: str, fallback) -> pd.Series:
    if col in df.columns:
        return df[col].reset_index(drop=True)
    return pd.Series(fallback)


def _normalize_filter_value(value) -> str:
    from apex.utils.astro_utils import normalize_filter_name
    import re
    raw = str(value or "").strip()
    if not raw:
        return ""
    # Try direct normalization first (preserves Johnson/SDSS case distinction)
    direct = normalize_filter_name(raw)
    if direct:
        return direct
    # Extract a token from compound strings like "Bessel_R" or "sdss-r_001"
    tokens = [t for t in re.split(r"[^a-zA-Z0-9']+", raw) if t]
    for token in reversed(tokens or [raw]):
        candidate = normalize_filter_name(token)
        if candidate:
            return candidate
    return raw


def _numeric_col(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.nan, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce")


def _bool_col(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(False, index=df.index, dtype=bool)
    s = df[col]
    if s.dtype == bool:
        return s.fillna(False).astype(bool)
    return s.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y", "on"})


def _finite_values(values) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr)]


def _finite_median(values) -> float:
    arr = _finite_values(values)
    return float(np.median(arr)) if arr.size else float("nan")


def _finite_percentile(values, q: float) -> float:
    arr = _finite_values(values)
    return float(np.percentile(arr, q)) if arr.size else float("nan")


def _finite_max(values) -> float:
    arr = _finite_values(values)
    return float(np.max(arr)) if arr.size else float("nan")


def _robust_frame_shift(dx, dy, *, sigma: float = 3.0, max_iter: int = 3) -> dict:
    """Robust shift-only registration from anchor residuals."""
    dx_arr = np.asarray(dx, dtype=float)
    dy_arr = np.asarray(dy, dtype=float)
    keep = np.isfinite(dx_arr) & np.isfinite(dy_arr)
    if int(keep.sum()) == 0:
        return {
            "dx": float("nan"),
            "dy": float("nan"),
            "rms": float("nan"),
            "p95": float("nan"),
            "keep": keep,
        }

    for _ in range(max(1, int(max_iter))):
        if int(keep.sum()) < 3:
            break
        med_dx = float(np.median(dx_arr[keep]))
        med_dy = float(np.median(dy_arr[keep]))
        resid = np.hypot(dx_arr - med_dx, dy_arr - med_dy)
        resid_keep = resid[keep]
        med_resid = float(np.median(resid_keep))
        mad = float(np.median(np.abs(resid_keep - med_resid)))
        scale = MAD_TO_SIGMA * mad
        if not np.isfinite(scale) or scale <= 0:
            scale = float(np.percentile(resid_keep, 68)) if resid_keep.size else float("nan")
        if not np.isfinite(scale) or scale <= 0:
            break
        new_keep = keep & (resid <= max(0.25, sigma * scale))
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep

    med_dx = float(np.median(dx_arr[keep]))
    med_dy = float(np.median(dy_arr[keep]))
    resid = np.hypot(dx_arr - med_dx, dy_arr - med_dy)
    resid_keep = resid[keep]
    rms = float(np.sqrt(np.mean(resid_keep * resid_keep))) if resid_keep.size else float("nan")
    p95 = float(np.percentile(resid_keep, 95)) if resid_keep.size else float("nan")
    return {"dx": med_dx, "dy": med_dy, "rms": rms, "p95": p95, "keep": keep}


def _safe_rate(num: int, den: int) -> float:
    return float(num) / float(den) if den > 0 else float("nan")


def _step4_quality_check(row: pd.Series, P, sat_adu: float, datamax_adu: float) -> Tuple[bool, str, bool]:
    """Reuse Step4 detection-quality columns when they are available."""
    reasons: List[str] = []
    used = False

    if "anchor_candidate" in row.index:
        used = True
        anchor_ok = str(row.get("anchor_candidate")).strip().lower() in {"1", "true", "t", "yes", "y"}
        if not anchor_ok:
            qflags = str(row.get("quality_flags", "") or "").strip()
            reasons.append(f"anchor_candidate_false:{qflags}" if qflags else "anchor_candidate_false")

    def _finite_row_float(name: str) -> float:
        nonlocal used
        if name not in row.index:
            return float("nan")
        try:
            v = float(row.get(name))
        except Exception:
            return float("nan")
        if np.isfinite(v):
            used = True
        return v

    elong = _finite_row_float("elongation")
    elong_max = _to_float(getattr(P, "ref_cat_max_elong", 1.5), 1.5)
    if np.isfinite(elong) and elong_max > 0 and elong > elong_max:
        reasons.append("elongation")

    roundness = _finite_row_float("roundness")
    if not np.isfinite(roundness):
        r1 = _finite_row_float("roundness1")
        r2 = _finite_row_float("roundness2")
        vals = [abs(v) for v in (r1, r2) if np.isfinite(v)]
        roundness = max(vals) if vals else float("nan")
    round_max = _to_float(getattr(P, "ref_cat_max_abs_round", 0.4), 0.4)
    if np.isfinite(roundness) and round_max > 0 and abs(roundness) > round_max:
        reasons.append("roundness")

    sharp = _finite_row_float("sharpness")
    sharp_min = _to_float(getattr(P, "ref_cat_sharp_min", 0.2), 0.2)
    sharp_max = _to_float(getattr(P, "ref_cat_sharp_max", 1.0), 1.0)
    if np.isfinite(sharp):
        if np.isfinite(sharp_min) and sharp < sharp_min:
            reasons.append("sharpness_low")
        if np.isfinite(sharp_max) and sharp_max > 0 and sharp > sharp_max:
            reasons.append("sharpness_high")

    peak = _finite_row_float("peak_adu")
    if np.isfinite(peak):
        if np.isfinite(sat_adu) and sat_adu > 0 and peak >= sat_adu:
            reasons.append("saturated_peak")
        elif np.isfinite(datamax_adu) and datamax_adu > 0 and peak >= datamax_adu:
            reasons.append("nonlinear_peak")

    status = str(row.get("fwhm_status", "") or "").strip().lower()
    if status:
        used = True
    if status in {"saturated", "edge"}:
        reasons.append(status)

    return (len(reasons) == 0), ";".join(reasons), used


def _frame_centering_stats(
    phot_df: pd.DataFrame,
    *,
    fname: str,
    filt: str,
    status: str,
    fwhm_px: float,
    wcs_ok: bool,
    apcorr: float,
    max_shift_px: float,
    outlier_px: float,
    out_path: str,
) -> dict:
    """Summarize forced-photometry centering reliability for one frame."""
    n_master = int(len(phot_df))
    if n_master == 0:
        return {
            "file": fname,
            "filter": filt,
            "status": status,
            "centering_status": "NO_SOURCES",
            "n_master": 0,
            "path": out_path,
        }

    detected = _bool_col(phot_df, "detected_flag")
    forced = _bool_col(phot_df, "forced_flag")
    recentered = _bool_col(phot_df, "recentered_flag")
    outlier = _bool_col(phot_df, "centroid_outlier")
    off_frame = _bool_col(phot_df, "off_frame_flag")
    x_fit = _numeric_col(phot_df, "x_fit")
    y_fit = _numeric_col(phot_df, "y_fit")
    on_frame = (~off_frame) & x_fit.notna() & y_fit.notna()

    match_offset = _numeric_col(phot_df, "match_offset_px")
    reg_resid = _numeric_col(phot_df, "registration_resid_px")
    reg_anchor = _bool_col(phot_df, "registration_anchor")
    frame_dx = _numeric_col(phot_df, "frame_dx_px")
    frame_dy = _numeric_col(phot_df, "frame_dy_px")
    center_error = _numeric_col(phot_df, "center_error_px")
    centroid_shift = _numeric_col(phot_df, "centroid_shift_px")
    mag_err = _numeric_col(phot_df, "mag_err")
    snr = _numeric_col(phot_df, "snr")

    n_on_frame = int(on_frame.sum())
    n_detected = int(detected.sum())
    n_forced = int(forced.sum())
    n_forced_on_frame = int((forced & on_frame).sum())
    n_recentered = int(recentered.sum())
    n_center_tested = int(center_error.notna().sum())
    n_outlier = int(outlier.sum())
    n_registration_anchor = int(reg_anchor.sum())

    det_match_offset = match_offset[detected & match_offset.notna()]
    reg_anchor_resid = reg_resid[reg_anchor & reg_resid.notna()]
    cen_error = center_error[center_error.notna()]
    cen_shift = centroid_shift[recentered & centroid_shift.notna()]

    high_mag_err = mag_err[outlier & mag_err.notna()]
    low_mag_err = mag_err[(~outlier) & detected & mag_err.notna()]
    high_snr = snr[outlier & snr.notna()]
    low_snr = snr[(~outlier) & detected & snr.notna()]
    high_mag_med = _finite_median(high_mag_err)
    low_mag_med = _finite_median(low_mag_err)
    high_shift_delta = (
        float(high_mag_med - low_mag_med)
        if np.isfinite(high_mag_med) and np.isfinite(low_mag_med)
        else float("nan")
    )

    center_p90 = _finite_percentile(cen_error, 90)
    outlier_rate = _safe_rate(n_outlier, n_center_tested)
    detected_rate = _safe_rate(n_detected, n_on_frame)
    forced_rate = _safe_rate(n_forced_on_frame, n_on_frame)
    recentered_rate = _safe_rate(n_recentered, n_detected)

    centering_status = "OK"
    advice = "centering stable"
    if n_on_frame <= 0:
        centering_status = "NO_ON_FRAME"
        advice = "projected master positions are outside this frame"
    elif n_center_tested < max(5, min(20, int(0.05 * n_on_frame))):
        centering_status = "LOW_MATCH"
        advice = "too few detected matches to verify centers"
    elif np.isfinite(outlier_rate) and outlier_rate >= 0.20:
        centering_status = "CHECK"
        advice = "many centers exceed the outlier threshold"
    elif np.isfinite(center_p90) and center_p90 >= max(max_shift_px * 0.85, outlier_px * 1.5):
        centering_status = "CHECK"
        advice = "p90 center error is close to the recenter limit"
    elif np.isfinite(outlier_rate) and outlier_rate >= 0.05:
        centering_status = "REVIEW"
        advice = "some sources have large center shifts"
    elif np.isfinite(forced_rate) and forced_rate >= 0.50:
        centering_status = "REVIEW"
        advice = "many sources are forced-only; inspect downstream scatter"

    return {
        "file": fname,
        "filter": filt,
        "status": status,
        "centering_status": centering_status,
        "advice": advice,
        "wcs_ok": bool(wcs_ok),
        "fwhm_px": float(fwhm_px) if np.isfinite(fwhm_px) else float("nan"),
        "max_recenter_shift_px": float(max_shift_px),
        "centroid_outlier_px": float(outlier_px),
        "apcorr": float(apcorr) if np.isfinite(apcorr) else float("nan"),
        "n_master": n_master,
        "n_on_frame": n_on_frame,
        "n_detected": n_detected,
        "n_forced": n_forced,
        "n_forced_on_frame": n_forced_on_frame,
        "n_recentered": n_recentered,
        "n_center_tested": n_center_tested,
        "n_centroid_outlier": n_outlier,
        "n_registration_anchor": n_registration_anchor,
        "detected_rate": detected_rate,
        "forced_rate": forced_rate,
        "recentered_rate": recentered_rate,
        "centroid_outlier_rate": outlier_rate,
        "frame_dx_px": _finite_median(frame_dx),
        "frame_dy_px": _finite_median(frame_dy),
        "registration_resid_med_px": _finite_median(reg_anchor_resid),
        "registration_resid_p95_px": _finite_percentile(reg_anchor_resid, 95),
        "match_offset_med_px": _finite_median(det_match_offset),
        "match_offset_p90_px": _finite_percentile(det_match_offset, 90),
        "match_offset_max_px": _finite_max(det_match_offset),
        "center_error_med_px": _finite_median(cen_error),
        "center_error_p90_px": center_p90,
        "center_error_max_px": _finite_max(cen_error),
        "centroid_shift_med_px": _finite_median(cen_shift),
        "centroid_shift_p90_px": _finite_percentile(cen_shift, 90),
        "centroid_shift_max_px": _finite_max(cen_shift),
        "mag_err_med": _finite_median(mag_err),
        "mag_err_p90": _finite_percentile(mag_err, 90),
        "mag_err_high_shift_med": high_mag_med,
        "mag_err_low_shift_med": low_mag_med,
        "mag_err_high_shift_delta": high_shift_delta,
        "snr_high_shift_med": _finite_median(high_snr),
        "snr_low_shift_med": _finite_median(low_snr),
        "path": out_path,
    }


def run_forced_photometry(
    file_list,
    params,
    data_dir,
    cache_dir,
    *,
    result_dir=None,
    output_dir=None,
    progress_cb=None,
    log_cb=None,
    worker_status_cb=None,
    apcorr_cb=None,
    center_stats_cb=None,
    error_cb=None,
    should_stop=None,
    logger=None,
) -> dict:
    """Headless forced aperture photometry over ``file_list``.

    Verbatim relocation of ``ForcedPhotWorker.run`` / ``_run_impl``: same
    algorithm, parameter lookups, rounding, and on-disk outputs. Qt signals are
    replaced by optional callbacks; the summary dict is returned instead of
    emitted.
    """
    # Bind call/state shims so the relocated body keeps using the same names.
    data_dir = Path(data_dir)
    cache_dir = Path(cache_dir)
    result_dir = Path(result_dir) if result_dir is not None else Path(getattr(params.P, "result_dir", data_dir))
    output_dir = Path(output_dir) if output_dir is not None else step7_forced_phot_dir(result_dir)
    file_list = list(file_list)
    _wcs_header_cache: Dict[str, fits.Header] = {}
    _wcs_cache_lock = Lock()
    _results_lock = Lock()
    max_workers = get_parallel_workers(params)

    def _stop_requested() -> bool:
        return should_stop() if should_stop else False

    def _log(msg: str):
        (log_cb(msg) if log_cb else None)

    # ── Path resolution ────────────────────────────────────────────────────────

    def _resolve_fits_path(fname: str) -> Optional[Path]:
        if crop_is_active(result_dir):
            p = step2_cropped_dir(result_dir) / fname
            if p.exists():
                return p
        try:
            orig = Path(params.get_file_path(fname))
            if orig.exists():
                return orig
        except Exception:
            pass
        cand = data_dir / fname
        return cand if cand.exists() else None

    def _load_header(fname: str) -> Optional[fits.Header]:
        path = _resolve_fits_path(fname)
        if path is None or not path.exists():
            return None
        key = str(path)
        with _wcs_cache_lock:
            hdr = _wcs_header_cache.get(key)
        if hdr is not None:
            return hdr
        try:
            hdr = fits.getheader(path)
        except Exception:
            return None
        with _wcs_cache_lock:
            _wcs_header_cache[key] = hdr
        return hdr

    def _load_wcs_for_frame(fname: str) -> Optional[WCS]:
        fits_path = _resolve_fits_path(fname)
        candidates: List[Path] = []
        if fits_path is not None and fits_path.exists():
            candidates.append(fits_path)
        try:
            orig = Path(params.get_file_path(fname))
            if orig not in candidates and orig.exists():
                candidates.append(orig)
        except Exception:
            pass

        for path in candidates:
            try:
                key = str(path)
                with _wcs_cache_lock:
                    hdr = _wcs_header_cache.get(key)
                if hdr is None:
                    hdr = fits.getheader(path)
                    with _wcs_cache_lock:
                        _wcs_header_cache[key] = hdr
                w = WCS(hdr, relax=True)
                if w.has_celestial:
                    return w
                for wcs_path in astap_wcs_candidates(path):
                    if not wcs_path.exists():
                        continue
                    wcs_dict = parse_astap_wcs_file(wcs_path)
                    if not wcs_dict:
                        continue
                    hdr2 = fits.Header()
                    for k, v in wcs_dict.items():
                        try:
                            hdr2[k] = v
                        except Exception:
                            pass
                    w2 = WCS(hdr2, relax=True)
                    if w2.has_celestial:
                        return w2
            except Exception:
                continue
        return None

    def _load_image(fname: str) -> Optional[np.ndarray]:
        path = _resolve_fits_path(fname)
        if path is None:
            return None
        try:
            with fits.open(path) as hdul:
                for hdu in hdul:
                    if hdu.data is not None and hdu.data.ndim == 2:
                        # float32, not float64: halves per-frame memory (a 61 MP
                        # frame is 244 MB vs 488 MB). Aperture sums stay exact to
                        # ~1e-7 relative — far below the photometric noise floor.
                        return hdu.data.astype(np.float32)
        except Exception:
            pass
        return None

    def _load_fwhm(fname: str, fallback: float = 6.0) -> float:
        for cand in [
            cache_dir / f"detect_{fname}.json",
            step4_dir(result_dir) / f"detect_{fname}.json",
        ]:
            if not cand.exists():
                continue
            try:
                meta = json.loads(cand.read_text(encoding="utf-8"))
                for key in ("fwhm_med_rad_px", "fwhm_med_px", "fwhm_px", "fwhm_med"):
                    v = meta.get(key)
                    if v is not None:
                        f = float(v)
                        if np.isfinite(f) and f > 0:
                            return f
            except Exception:
                continue
        return fallback

    def _load_sky_stats(fname: str) -> Tuple[float, float]:
        for cand in [
            cache_dir / f"detect_{fname}.json",
            step4_dir(result_dir) / f"detect_{fname}.json",
        ]:
            if not cand.exists():
                continue
            try:
                meta = json.loads(cand.read_text(encoding="utf-8"))
                sky_med = float(meta.get("sky_med") or meta.get("bkg_median") or np.nan)
                sky_sig = float(meta.get("sky_sigma") or meta.get("bkg_rms") or np.nan)
                return sky_med, sky_sig
            except Exception:
                continue
        return np.nan, np.nan

    def _load_detect_positions(fname: str) -> Optional[pd.DataFrame]:
        for cand in [
            cache_dir / f"detect_{fname}.csv",
            step4_dir(result_dir) / f"detect_{fname}.csv",
        ]:
            if not cand.exists() or cand.stat().st_size == 0:
                continue
            try:
                df = pd.read_csv(cand)
                xc = "x" if "x" in df.columns else ("xcenter" if "xcenter" in df.columns else None)
                yc = "y" if "y" in df.columns else ("ycenter" if "ycenter" in df.columns else None)
                if xc is None or yc is None:
                    continue
                out = pd.DataFrame({
                    "det_uid": (
                        pd.to_numeric(df["det_uid"], errors="coerce").fillna(
                            pd.Series(range(len(df)))
                        ).astype(int)
                        if "det_uid" in df.columns
                        else pd.Series(range(len(df)))
                    ),
                    "x": pd.to_numeric(df[xc], errors="coerce"),
                    "y": pd.to_numeric(df[yc], errors="coerce"),
                })
                for col in (
                    "elongation", "roundness", "roundness1", "roundness2",
                    "sharpness", "dao_peak", "dao_flux", "peak_adu",
                    "fwhm_px", "sep_flag", "nearest_neighbor_px",
                    "nearest_neighbor_fwhm", "edge_margin_px",
                    "fwhm_ratio_to_frame", "flux_for_quality",
                    "flux_percentile", "peak_fraction_to_sat",
                    "quality_score",
                ):
                    if col in df.columns:
                        out[col] = pd.to_numeric(df[col], errors="coerce")
                for col in ("fwhm_status", "source_type", "quality_flags"):
                    if col in df.columns:
                        out[col] = df[col].astype(str).reset_index(drop=True)
                for col in ("anchor_candidate", "apcorr_candidate", "epsf_candidate", "psf_seed_candidate"):
                    if col in df.columns:
                        out[col] = (
                            df[col]
                            .astype(str)
                            .str.strip()
                            .str.lower()
                            .isin({"1", "true", "t", "yes", "y"})
                            .reset_index(drop=True)
                        )
                return out.dropna(subset=["x", "y"]).reset_index(drop=True)
            except Exception:
                continue
        return None

    def _load_master_catalog(filt: str) -> Optional[pd.DataFrame]:
        refbuild_dir = step6_refbuild_dir(result_dir)
        candidates = [
            refbuild_dir / f"ref_catalog_{filt}.tsv",
            refbuild_dir / f"ref_catalog_{filt.lower()}.tsv",
            refbuild_dir / f"ref_catalog_{filt.upper()}.tsv",
            refbuild_dir / "ref_catalog.tsv",
        ]
        for path in candidates:
            if not path.exists():
                continue
            try:
                df = pd.read_csv(path, sep="\t")
                has_radec = {"ra_deg", "dec_deg"} <= set(df.columns)
                has_xy = {"x_ref", "y_ref"} <= set(df.columns)
                if not has_radec and not has_xy:
                    continue
                if "master_id" not in df.columns:
                    if "ID" in df.columns:
                        df.insert(0, "master_id", df["ID"])
                    elif "source_id" in df.columns:
                        df.insert(0, "master_id", df["source_id"])
                    else:
                        df.insert(0, "master_id", range(1, len(df) + 1))
                if "source_id" not in df.columns:
                    df["source_id"] = df["master_id"]
                if "ID" not in df.columns:
                    df["ID"] = df["master_id"]
                for col in ("crowding_flag", "neighbor_dist_px"):
                    if col not in df.columns:
                        df[col] = False if col == "crowding_flag" else np.nan
                return df.reset_index(drop=True)
            except Exception as exc:
                _log(f"[FORCED] catalog load failed ({path.name}): {exc}")
        return None

    # ── Per-frame photometry ───────────────────────────────────────────────────

    def _phot_frame(
        fname: str,
        master_df: pd.DataFrame,
        img: np.ndarray,
        header: Optional[fits.Header],
        wcs: Optional[WCS],
        fwhm_px: float,
        P,
        status_cb=None,
    ) -> Tuple[pd.DataFrame, float, dict]:
        """Measure forced aperture photometry + growth-curve apcorr for one frame.

        Returns
        -------
        (phot_df, apcorr, growth_curve_data)
        """
        def _status(label: str, pct: int):
            if callable(status_cb):
                try:
                    status_cb(label, int(pct))
                except Exception:
                    pass

        n = len(master_df)
        h, w = img.shape

        # Photometry parameters
        r_ap_scale   = _to_float(getattr(P, "forced_r_ap_scale",     0.8), 0.8)
        ref_ap_scale = _to_float(getattr(P, "forced_ref_ap_scale",   2.4), 2.4)
        min_r_ap     = _to_float(getattr(P, "min_r_ap_px",           4.0), 4.0)
        ann_scale    = _to_float(getattr(P, "fitsky_annulus_scale",   4.0), 4.0)
        dann_scale   = _to_float(getattr(P, "fitsky_dannulus_scale",  2.0), 2.0)
        ann_gap      = _to_float(getattr(P, "annulus_min_gap_px",     6.0), 6.0)
        max_shift    = _to_float(getattr(P, "max_recenter_shift",     2.0), 2.0)
        outlier_px   = _to_float(getattr(P, "centroid_outlier_px",    1.0), 1.0)
        recenter_enabled = bool(getattr(P, "recenter_aperture", True))
        noise        = resolve_effective_noise_params(P, header)
        gain         = noise.gain_e_per_adu
        rn_e         = noise.rdnoise_e
        exptime, _exptime_found = _exptime_from_header(header)
        if not _exptime_found:
            _log(
                f"[FORCED][{fname}] WARNING: no usable EXPTIME in header "
                f"(tried {', '.join(EXPTIME_HEADER_KEYS)}); using {exptime:.3g}s for "
                f"count-rate magnitudes — in a mixed-exposure set this frame may be "
                f"on the wrong scale."
            )
        noise_info = {
            "gain_e_per_adu": float(noise.gain_e_per_adu),
            "rdnoise_e": float(noise.rdnoise_e),
            "binning_x": int(noise.bin_x),
            "binning_y": int(noise.bin_y),
            "gain_source": noise.gain_source,
            "rdnoise_source": noise.rdnoise_source,
        }
        sat_adu      = _to_float(getattr(P, "saturation_adu",     65000.0), 65000.0)
        datamax_adu  = _to_float(getattr(P, "datamax_adu",        55000.0), 55000.0)
        sigma_clip   = _to_float(getattr(P, "phot_sigma_clip",        3.0), 3.0)
        max_iter     = _to_int(  getattr(P, "phot_max_iter",            5), 5)
        sky_mode     = str(getattr(P, "sky_sigma_mode", "local") or "local").strip().lower()
        sky_incl_rn  = bool(getattr(P, "sky_sigma_includes_rn", True))
        min_n_sky    = _to_int(getattr(P, "sky_sigma_min_n_sky", 50), 50)
        cbox_scale   = _to_float(getattr(P, "center_cbox_scale", 1.5), 1.5)
        apcorr_min_snr = _to_float(getattr(P, "apcorr_min_snr", 40.0), 40.0)
        apcorr_min_n   = _to_int(getattr(P, "apcorr_use_min_n", 12), 12)
        apcorr_max_sources = _to_int(getattr(P, "apcorr_max_sources", 250), 250)

        r_ap  = max(min_r_ap, r_ap_scale * fwhm_px)
        r_ref = max(r_ap + 2.0, ref_ap_scale * fwhm_px)
        r_in  = max(r_ref + ann_gap, ann_scale * fwhm_px)
        r_out = r_in + max(ann_gap, dann_scale * fwhm_px)

        sky_med_frame, sky_sigma_frame = _load_sky_stats(fname)
        sky_e = sky_sigma_frame * gain if np.isfinite(sky_sigma_frame) else np.nan

        # --- Project master positions to pixels via WCS ---
        has_radec = {"ra_deg", "dec_deg"} <= set(master_df.columns)
        has_xy    = {"x_ref", "y_ref"} <= set(master_df.columns)

        x_pred = np.full(n, np.nan)
        y_pred = np.full(n, np.nan)

        _status("Project centers", 28)
        if has_radec and wcs is not None:
            ra  = pd.to_numeric(master_df["ra_deg"],  errors="coerce").to_numpy(float)
            dec = pd.to_numeric(master_df["dec_deg"], errors="coerce").to_numpy(float)
            valid = np.isfinite(ra) & np.isfinite(dec)
            if valid.any():
                try:
                    xy = wcs.all_world2pix(
                        np.column_stack([ra[valid], dec[valid]]), 0
                    )
                    x_pred[valid] = xy[:, 0]
                    y_pred[valid] = xy[:, 1]
                except Exception as exc:
                    _log(f"[FORCED][{fname}] WCS project failed: {exc}")

        if has_xy:
            x_ref_col = pd.to_numeric(master_df["x_ref"], errors="coerce").to_numpy(float)
            y_ref_col = pd.to_numeric(master_df["y_ref"], errors="coerce").to_numpy(float)
            missing = ~np.isfinite(x_pred)
            x_pred[missing] = x_ref_col[missing]
            y_pred[missing] = y_ref_col[missing]

        # --- Match to detections for recentering ---
        _status("Match detections", 34)
        det_df = _load_detect_positions(fname)
        detected_flag = np.zeros(n, dtype=bool)
        det_uid_arr   = np.full(n, -1, dtype=np.int64)
        det_x_arr = np.full(n, np.nan)
        det_y_arr = np.full(n, np.nan)
        match_dx = np.full(n, np.nan)
        match_dy = np.full(n, np.nan)
        match_offset = np.full(n, np.nan)
        step4_quality_ok = np.ones(n, dtype=bool)
        step4_quality_used = np.zeros(n, dtype=bool)
        step4_quality_reason = np.full(n, "", dtype=object)
        step4_quality_score = np.full(n, np.nan)
        step4_quality_flags = np.full(n, "", dtype=object)
        step4_anchor_candidate = np.zeros(n, dtype=bool)
        step4_apcorr_candidate = np.zeros(n, dtype=bool)
        registration_match_radius = _to_float(
            getattr(P, "registration_match_radius_px", max(max_shift * 2.0, 4.0)),
            max(max_shift * 2.0, 4.0),
        )
        registration_match_radius = max(registration_match_radius, max_shift)

        if det_df is not None and len(det_df) > 0:
            det_xy = det_df[["x", "y"]].to_numpy(float)
            master_xy_pred = np.column_stack([x_pred, y_pred])
            valid_pred = np.isfinite(x_pred) & np.isfinite(y_pred)
            if valid_pred.any():
                tree = KDTree(det_xy)
                dists, idxs = tree.query(master_xy_pred[valid_pred], k=1)
                match_mask = dists <= registration_match_radius
                valid_indices = np.where(valid_pred)[0]
                # Check column presence once before the loop (avoid per-iteration overhead).
                _has_uid    = "det_uid"          in det_df.columns
                _has_qscore = "quality_score"    in det_df.columns
                _has_qflags = "quality_flags"    in det_df.columns
                _has_anchor = "anchor_candidate" in det_df.columns
                _has_apcorr = "apcorr_candidate" in det_df.columns
                for vi, matched, idx, dist in zip(valid_indices, match_mask, idxs, dists):
                    if not matched:
                        continue
                    # Single iloc per matched star — reuse the row Series for all accesses.
                    row = det_df.iloc[idx]
                    x_d = float(row["x"])
                    y_d = float(row["y"])
                    detected_flag[vi] = bool(dist <= max_shift)
                    det_uid_arr[vi]   = int(row["det_uid"]) if _has_uid else -1
                    det_x_arr[vi]     = x_d
                    det_y_arr[vi]     = y_d
                    match_dx[vi]      = x_d - float(x_pred[vi])
                    match_dy[vi]      = y_d - float(y_pred[vi])
                    match_offset[vi]  = float(dist)
                    ok, reason, used  = _step4_quality_check(row, P, sat_adu, datamax_adu)
                    step4_quality_ok[vi]     = ok
                    step4_quality_used[vi]   = used
                    step4_quality_reason[vi] = reason
                    if _has_qscore:
                        step4_quality_score[vi]   = float(row["quality_score"])
                    if _has_qflags:
                        step4_quality_flags[vi]   = str(row["quality_flags"])
                    if _has_anchor:
                        step4_anchor_candidate[vi] = bool(row["anchor_candidate"])
                    if _has_apcorr:
                        step4_apcorr_candidate[vi] = bool(row["apcorr_candidate"])

        # --- Frame-level registration from good anchor stars ---
        crowding_flag = (
            master_df["crowding_flag"].to_numpy(bool)
            if "crowding_flag" in master_df.columns
            else np.zeros(n, dtype=bool)
        )
        has_step4_anchor = det_df is not None and "anchor_candidate" in det_df.columns
        has_step4_apcorr = det_df is not None and "apcorr_candidate" in det_df.columns
        matched_for_registration = np.isfinite(match_dx) & np.isfinite(match_dy)
        registration_anchor = (
            matched_for_registration
            & ~crowding_flag
            & step4_quality_ok
        )
        if has_step4_anchor:
            registration_anchor &= step4_anchor_candidate
        n_anchor_initial = int(registration_anchor.sum())
        registration_method = "wcs_only"
        frame_dx = float("nan")
        frame_dy = float("nan")
        registration_rms = float("nan")
        registration_p95 = float("nan")
        registration_resid = np.full(n, np.nan)
        min_reg_anchors = max(3, _to_int(getattr(P, "registration_min_anchors", 6), 6))

        if n_anchor_initial >= min_reg_anchors:
            reg = _robust_frame_shift(match_dx[registration_anchor], match_dy[registration_anchor])
            anchor_indices = np.flatnonzero(registration_anchor)
            keep_local = np.asarray(reg["keep"], dtype=bool)
            kept_indices = anchor_indices[keep_local]
            if len(kept_indices) >= min_reg_anchors:
                registration_anchor[:] = False
                registration_anchor[kept_indices] = True
                frame_dx = float(reg["dx"])
                frame_dy = float(reg["dy"])
                registration_rms = float(reg["rms"])
                registration_p95 = float(reg["p95"])
                registration_method = "anchor_shift"
                _log(
                    f"[FORCED][{fname}] registration shift "
                    f"dx={frame_dx:.3f} dy={frame_dy:.3f} px | "
                    f"anchors={len(kept_indices)}/{n_anchor_initial} | rms={registration_rms:.3f}"
                )
            else:
                registration_anchor[:] = False
                _log(
                    f"[WARN][FORCED][{fname}] registration rejected after clipping: "
                    f"anchors={len(kept_indices)}/{n_anchor_initial} < {min_reg_anchors}"
                )
        elif has_step4_anchor and int((matched_for_registration & ~crowding_flag & step4_quality_ok).sum()) >= min_reg_anchors:
            fallback_anchor = (
                matched_for_registration
                & ~crowding_flag
                & step4_quality_ok
            )
            reg = _robust_frame_shift(match_dx[fallback_anchor], match_dy[fallback_anchor])
            anchor_indices = np.flatnonzero(fallback_anchor)
            keep_local = np.asarray(reg["keep"], dtype=bool)
            kept_indices = anchor_indices[keep_local]
            if len(kept_indices) >= min_reg_anchors:
                registration_anchor[:] = False
                registration_anchor[kept_indices] = True
                frame_dx = float(reg["dx"])
                frame_dy = float(reg["dy"])
                registration_rms = float(reg["rms"])
                registration_p95 = float(reg["p95"])
                registration_method = "quality_shift"
                _log(
                    f"[WARN][FORCED][{fname}] too few Step4 anchors ({n_anchor_initial}); "
                    f"using quality-filtered shift anchors={len(kept_indices)}"
                )
            else:
                registration_anchor[:] = False
                _log(
                    f"[WARN][FORCED][{fname}] quality registration rejected after clipping: "
                    f"anchors={len(kept_indices)} < {min_reg_anchors}"
                )

        x_reg = x_pred.copy()
        y_reg = y_pred.copy()
        if np.isfinite(frame_dx) and np.isfinite(frame_dy):
            x_reg = x_pred + frame_dx
            y_reg = y_pred + frame_dy
            registration_resid = np.hypot(det_x_arr - x_reg, det_y_arr - y_reg)
        registration_match_offset = np.where(np.isfinite(registration_resid), registration_resid, match_offset)
        detected_flag = detected_flag | (
            np.isfinite(registration_match_offset)
            & (registration_match_offset <= max_shift)
        )

        # --- Recenter detected sources ---
        _status(f"Recenter {int(detected_flag.sum())}", 40)
        x_fit = x_reg.copy()
        y_fit = y_reg.copy()
        centroid_shift = np.full(n, np.nan)
        recentered_flag = np.zeros(n, dtype=bool)
        recenter_method = np.full(n, registration_method, dtype=object)

        if not recenter_enabled:
            recenter_method[detected_flag] = "disabled"
        else:
            n_detected_total = int(detected_flag.sum())
            recenter_seen = 0
            for i in range(n):
                if not detected_flag[i]:
                    continue
                recenter_seen += 1
                if recenter_seen == 1 or recenter_seen % 250 == 0:
                    pct = 40 + min(12, int(12 * recenter_seen / max(n_detected_total, 1)))
                    _status(f"Recenter {recenter_seen}/{n_detected_total}", pct)
                xi, yi = x_reg[i], y_reg[i]
                if not (np.isfinite(xi) and np.isfinite(yi)):
                    recenter_method[i] = "invalid_position"
                    continue

                seed_x = det_x_arr[i] if np.isfinite(det_x_arr[i]) else xi
                seed_y = det_y_arr[i] if np.isfinite(det_y_arr[i]) else yi
                accepted = False
                try:
                    xc, yc = refine_local_centroid(img, seed_x, seed_y, fwhm_px, cbox_scale)
                    shift = np.hypot(xc - xi, yc - yi)
                    if np.isfinite(shift) and shift <= max_shift:
                        x_fit[i] = xc
                        y_fit[i] = yc
                        centroid_shift[i] = shift
                        recentered_flag[i] = True
                        recenter_method[i] = "local_centroid"
                        accepted = True
                except Exception:
                    pass

                if (
                    not accepted
                    and np.isfinite(det_x_arr[i])
                    and np.isfinite(det_y_arr[i])
                    and np.isfinite(registration_match_offset[i])
                    and registration_match_offset[i] <= max_shift
                ):
                    x_fit[i] = det_x_arr[i]
                    y_fit[i] = det_y_arr[i]
                    centroid_shift[i] = registration_match_offset[i]
                    recentered_flag[i] = True
                    recenter_method[i] = "detected_xy"
                    accepted = True

                if not accepted:
                    recenter_method[i] = "failed"

        center_error = np.where(np.isfinite(centroid_shift), centroid_shift, registration_match_offset)
        centroid_outlier = (
            detected_flag &
            np.isfinite(center_error) &
            (center_error > outlier_px)
        )

        # --- Vectorized photometry at r_ap ---
        _status("Aperture photometry", 56)
        positions = np.column_stack([x_fit, y_fit])
        phot_kw = dict(
            gain=gain, rn_param_e=rn_e, sky_frame_e=sky_e,
            sky_sigma_mode=sky_mode, sky_sigma_includes_rn=sky_incl_rn,
            min_n_sky_for_local=min_n_sky, sat_adu=sat_adu,
            datamax_adu=datamax_adu, sigma_clip_val=sigma_clip, maxiters=max_iter,
        )
        (
            flux_arr, flux_err_arr, snr_arr, sky_arr,
            is_sat_arr, is_nl_arr, sky_std_arr, n_sky_arr,
        ) = phot_vectorized(
            img, positions, r_ap, r_in, r_out,
            return_sky_stats=True,
            **phot_kw,
        )

        # --- Growth-curve aperture correction (vectorized: N_radii calls, not N_stars×N_radii) ---
        apcorr_mask = (
            detected_flag &
            ~crowding_flag &
            step4_quality_ok &
            ((not has_step4_apcorr) | step4_apcorr_candidate) &
            ~centroid_outlier &
            np.isfinite(flux_arr) & (flux_arr > 0) &
            np.isfinite(snr_arr) & (snr_arr >= apcorr_min_snr) &
            ~is_sat_arr &
            ~is_nl_arr
        )
        n_step4_quality_reject = int((detected_flag & ~step4_quality_ok).sum())
        n_step4_apcorr_reject = int((detected_flag & ~step4_apcorr_candidate).sum()) if has_step4_apcorr else 0
        n_center_outlier_reject = int(centroid_outlier.sum())
        n_apcorr_candidates = int(apcorr_mask.sum())
        if apcorr_max_sources > 0 and n_apcorr_candidates > apcorr_max_sources:
            candidate_idx = np.flatnonzero(apcorr_mask)
            candidate_snr = np.asarray(snr_arr[candidate_idx], dtype=float)
            rank = np.argsort(np.nan_to_num(candidate_snr, nan=-np.inf))[::-1]
            keep_idx = candidate_idx[rank[:apcorr_max_sources]]
            capped_mask = np.zeros(n, dtype=bool)
            capped_mask[keep_idx] = True
            apcorr_mask = capped_mask
            _log(
                f"[FORCED][{fname}] apcorr refs capped "
                f"{n_apcorr_candidates}->{int(apcorr_mask.sum())} "
                f"(apcorr_max_sources={apcorr_max_sources})"
            )
        gc_positions = positions[apcorr_mask]
        n_gc_stars = int(apcorr_mask.sum())
        _status(f"Apcorr refs {n_gc_stars}/{n_apcorr_candidates}", 68)

        apcorr = 1.0
        growth_curve_data: dict = {"noise": dict(noise_info)}

        if n_gc_stars >= apcorr_min_n:
            gc_radii = np.linspace(max(2.0, r_ap * 0.4), r_ref * 1.15, _GC_N_STEPS)
            gc_flux_matrix, gc_ferr_matrix = _growth_curve_fixed_sky(
                img, gc_positions, gc_radii, r_in, r_out, phot_kw,
                sky_median=sky_arr[apcorr_mask],
                sky_std=sky_std_arr[apcorr_mask],
                n_sky=n_sky_arr[apcorr_mask],
                status_cb=lambda ri, total: _status(
                    f"Apcorr radius {ri}/{total}",
                    70 + int(10 * ri / max(total, 1)),
                ),
            )

            # Normalize each star by its flux at the outermost radius
            f_outer = gc_flux_matrix[-1, :]
            valid_ref = np.isfinite(f_outer) & (f_outer > 0)
            if valid_ref.sum() >= apcorr_min_n:
                gc_norm = gc_flux_matrix[:, valid_ref] / f_outer[valid_ref]  # (N_radii, N_valid)
                gc_curve = np.nanmedian(gc_norm, axis=1)                     # (N_radii,)

                # Per-radius median magnitude error from flux/flux_err
                f_pos = gc_flux_matrix[:, valid_ref]
                fe_pos = gc_ferr_matrix[:, valid_ref]
                with np.errstate(divide="ignore", invalid="ignore"):
                    mag_err_matrix = (2.5 / np.log(10.0)) * (fe_pos / f_pos)
                mag_err_matrix = np.where(
                    np.isfinite(mag_err_matrix) & (f_pos > 0), mag_err_matrix, np.nan
                )
                med_mag_err = np.nanmedian(mag_err_matrix, axis=1)           # (N_radii,)

                # apcorr = 1 / enclosed_fraction at r_ap
                r_ap_idx = int(np.argmin(np.abs(gc_radii - r_ap)))
                enc_frac = float(gc_curve[r_ap_idx]) if np.isfinite(gc_curve[r_ap_idx]) else 0.0
                apcorr = float(1.0 / enc_frac) if enc_frac > 0.05 else 1.0

                # Optimal aperture radius = minimum of mag_err U-curve
                err_finite = np.where(np.isfinite(med_mag_err), med_mag_err, np.inf)
                r_opt_idx = int(np.argmin(err_finite)) if np.any(np.isfinite(med_mag_err)) else r_ap_idx
                r_opt_px  = float(gc_radii[r_opt_idx])
                err_opt   = float(med_mag_err[r_opt_idx]) if np.isfinite(med_mag_err[r_opt_idx]) else float("nan")

                growth_curve_data = {
                    "noise":          dict(noise_info),
                    "radii_px":       gc_radii.tolist(),
                    "enclosed_frac":  [float(v) for v in gc_curve],
                    "mag_err":        [float(v) for v in med_mag_err],
                    "n_stars":        int(valid_ref.sum()),
                    "n_apcorr_candidates": int(n_apcorr_candidates),
                    "n_apcorr_used":  int(n_gc_stars),
                    "n_step4_quality_reject": int(n_step4_quality_reject),
                    "n_step4_apcorr_reject": int(n_step4_apcorr_reject),
                    "n_center_outlier_reject": int(n_center_outlier_reject),
                    "apcorr_max_sources": int(apcorr_max_sources),
                    "r_ap_px":        float(r_ap),
                    "r_ref_px":       float(r_ref),
                    "r_opt_px":       r_opt_px,
                    "mag_err_opt":    err_opt,
                    "fwhm_px":        float(fwhm_px),
                    "fname":          fname,
                }
            else:
                # Too few valid reference stars — fall back to simple ratio
                _status("Apcorr fallback", 80)
                apcorr = _simple_apcorr(flux_arr, apcorr_mask, img, positions, r_ref, r_in, r_out, phot_kw)
        elif n_gc_stars > 0:
            _status("Apcorr fallback", 80)
            apcorr = _simple_apcorr(flux_arr, apcorr_mask, img, positions, r_ref, r_in, r_out, phot_kw)

        # --- Apply apcorr ---
        _status("Finalize photometry", 84)
        flux_corr     = flux_arr * apcorr
        flux_err_corr = flux_err_arr * apcorr

        # IRAF-style instrumental magnitude: count rate (e-/s) with a fixed
        # zeropoint baked in (see INSTRUMENTAL_ZMAG). Normalizing by exptime
        # makes magnitudes from frames of different exposure time directly
        # combinable; the +zmag constant makes mag_inst a real (positive)
        # magnitude. The per-frame ZP in Step 10 (median(gaia - mag_inst))
        # re-absorbs both automatically, so no downstream retuning is needed.
        rate = np.where(
            np.isfinite(flux_corr) & (flux_corr > 0),
            flux_corr / exptime,
            np.nan,
        )
        mag_inst = np.where(
            np.isfinite(rate) & (rate > 0),
            INSTRUMENTAL_ZMAG - 2.5 * np.log10(np.where(rate > 0, rate, np.nan)),
            np.nan,
        )
        mag_err = np.where(
            np.isfinite(flux_corr) & (flux_corr > 0) &
            np.isfinite(flux_err_corr) & (flux_err_corr >= 0),
            MAG_ERR_COEFF * flux_err_corr / np.where(flux_corr > 0, flux_corr, np.nan),
            np.nan,
        )
        bad_phot = is_sat_arr | is_nl_arr | ~np.isfinite(flux_corr)
        off_frame = (
            ~np.isfinite(x_fit) | ~np.isfinite(y_fit) |
            (x_fit < 0) | (x_fit >= w) |
            (y_fit < 0) | (y_fit >= h)
        )
        centering_quality = np.full(n, "forced_only", dtype=object)
        centering_quality[off_frame] = "off_frame"
        centering_quality[detected_flag & ~centroid_outlier] = "ok"
        centering_quality[centroid_outlier] = "warn_shift"
        # recenter_enabled is a scalar bool; ``~`` on a Python bool is a deprecated
        # bitwise inversion (NOT logical negation). When recentring is off, every
        # detected star is "matched but not recentred".
        if not recenter_enabled:
            centering_quality[detected_flag] = "matched_no_recenter"

        master_id  = _catalog_series(master_df, "master_id", np.arange(1, n + 1))
        source_id  = _catalog_series(master_df, "source_id", master_id)
        display_id = _catalog_series(master_df, "ID", master_id)

        out = pd.DataFrame({
            "master_id":         master_id,
            "source_id":         source_id,
            "ID":                display_id,
            "det_uid":           det_uid_arr,
            "x":                 x_fit,
            "y":                 y_fit,
            "x_pred":            x_pred,
            "y_pred":            y_pred,
            "x_reg":             x_reg,
            "y_reg":             y_reg,
            "x_fit":             x_fit,
            "y_fit":             y_fit,
            "det_x":             det_x_arr,
            "det_y":             det_y_arr,
            "match_dx_px":       match_dx,
            "match_dy_px":       match_dy,
            "match_offset_px":   match_offset,
            "frame_dx_px":       frame_dx,
            "frame_dy_px":       frame_dy,
            "registration_method": registration_method,
            "registration_anchor": registration_anchor,
            "registration_resid_px": registration_resid,
            "registration_rms_px": registration_rms,
            "registration_p95_px": registration_p95,
            "detected_flag":     detected_flag,
            "forced_flag":       ~detected_flag,
            "centroid_shift_px": centroid_shift,
            "center_error_px":   center_error,
            "recentered_flag":   recentered_flag,
            "recenter_method":   recenter_method,
            "step4_quality_ok":  step4_quality_ok,
            "step4_quality_used": step4_quality_used,
            "step4_quality_reason": step4_quality_reason,
            "step4_quality_score": step4_quality_score,
            "step4_quality_flags": step4_quality_flags,
            "step4_anchor_candidate": step4_anchor_candidate,
            "step4_apcorr_candidate": step4_apcorr_candidate,
            "centroid_outlier":  centroid_outlier,
            "centering_quality": centering_quality,
            "off_frame_flag":    off_frame,
            "flux":              flux_corr,
            "flux_e":            flux_corr,
            "flux_net_adu":      flux_corr / max(gain, 1e-12),
            "flux_err":          flux_err_corr,
            "gain_e_per_adu":    np.full(n, gain),
            "exptime":           np.full(n, float(exptime)),
            "rdnoise_e":         np.full(n, rn_e),
            "binning_x":         np.full(n, noise.bin_x),
            "binning_y":         np.full(n, noise.bin_y),
            "gain_source":       np.full(n, noise.gain_source, dtype=object),
            "rdnoise_source":    np.full(n, noise.rdnoise_source, dtype=object),
            "mag_inst":          mag_inst,
            "mag_err":           mag_err,
            "snr":               snr_arr,
            "sky":               sky_arr,
            "sky_std":           sky_std_arr,
            "bkg_median_adu":    sky_arr,
            "bkg_std_adu":       sky_std_arr,
            "n_sky":             n_sky_arr,
            "apcorr":            apcorr,
            "apcorr_applied":    np.full(n, bool(np.isfinite(apcorr) and abs(apcorr - 1.0) > 1e-6)),
            "apcorr_ref_count":  np.full(n, n_gc_stars),
            "is_saturated":      is_sat_arr,
            "is_nonlinear":      is_nl_arr,
            "bad_phot_flag":     bad_phot,
        })
        for col in ("ra_deg", "dec_deg", "gaia_source_id"):
            if col in master_df.columns and col not in out.columns:
                out[col] = master_df[col].reset_index(drop=True)

        return out, apcorr, growth_curve_data

    # ── Main run ──────────────────────────────────────────────────────────────

    run_t0 = time.time()
    P = params.P
    out_dir = output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    fwhm_guess = _to_float(getattr(P, "fwhm_pix_guess", 6.0), 6.0)

    def _get_filter(fname: str) -> str:
        fits_path = _resolve_fits_path(fname)
        if fits_path is not None and fits_path.exists():
            try:
                hdr = fits.getheader(fits_path)
                for key in ("FILTER", "FILTER1", "FILTER2", "FILTNAM", "FILTERID"):
                    filt = _normalize_filter_value(hdr.get(key))
                    if filt:
                        return filt
            except Exception:
                pass
        import re
        m = re.search(r"[-_]([ugrizbvUGRIZBV])[-_.]", str(fname), re.IGNORECASE)
        return _normalize_filter_value(m.group(1)) if m else "unknown"

    filter_map: Dict[str, List[str]] = {}
    for fname in file_list:
        filt = _get_filter(fname)
        filter_map.setdefault(filt, []).append(fname)

    try:
        (out_dir / "filter_frames.json").write_text(
            json.dumps(filter_map, indent=2, ensure_ascii=False), encoding="utf-8"
        )
    except Exception as exc:
        _log(f"[FORCED] filter_frames write failed: {exc}")

    master_catalogs: Dict[str, Optional[pd.DataFrame]] = {}
    for filt in filter_map:
        cat = _load_master_catalog(filt)
        if cat is None:
            cat = _load_master_catalog("r") or _load_master_catalog("unknown")
        master_catalogs[filt] = cat
        if cat is not None:
            _log(f"[FORCED] Loaded master catalog for filter '{filt}': {len(cat)} sources")
        else:
            _log(f"[FORCED] WARNING: No master catalog found for filter '{filt}'")

    try:
        master_rows = []
        for filt, cat in master_catalogs.items():
            if cat is None or cat.empty:
                continue
            tmp = cat.copy()
            tmp.insert(0, "filter", filt)
            master_rows.append(tmp)
        if master_rows:
            pd.concat(master_rows, ignore_index=True).drop_duplicates(
                subset=[c for c in ("filter", "source_id") if c in master_rows[0].columns]
            ).to_csv(out_dir / "master_sources.csv", index=False)
    except Exception as exc:
        _log(f"[FORCED] master_sources write failed: {exc}")

    total = len(file_list)
    index_rows: List[dict] = []
    apcorr_rows: List[dict] = []
    center_stats_rows: List[dict] = []
    n_done = [0]

    max_workers = max(1, int(max_workers))
    from queue import Queue as _Queue
    slot_queue: _Queue = _Queue()
    for _i in range(max_workers):
        slot_queue.put(_i)

    def _process_one(filt: str, fname: str, master_df: pd.DataFrame):
        if _stop_requested():
            return None

        slot = slot_queue.get()
        t0 = time.time()
        try:
            (worker_status_cb(slot, fname, "Loading", 5) if worker_status_cb else None)
            img = _load_image(fname)
            t_load = time.time()
            if img is None:
                (worker_status_cb(slot, fname, "No image", 100) if worker_status_cb else None)
                return ("no_image", fname, filt, None, None, None, None)

            (worker_status_cb(slot, fname, "WCS / FWHM", 20) if worker_status_cb else None)
            header = _load_header(fname)
            wcs = _load_wcs_for_frame(fname)
            fwhm_px = _load_fwhm(fname, fallback=fwhm_guess)
            t_wcs = time.time()

            (worker_status_cb(slot, fname, "Photometry", 45) if worker_status_cb else None)
            try:
                def _status(label: str, pct: int):
                    (worker_status_cb(slot, fname, label, pct) if worker_status_cb else None)

                phot_df, apcorr_val, gc_data = _phot_frame(
                    fname, master_df, img, header, wcs, fwhm_px, P, status_cb=_status
                )
                t_phot = time.time()
            except Exception as exc:
                tb = traceback.format_exc()
                (worker_status_cb(slot, fname, "Error", 100) if worker_status_cb else None)
                return ("error", fname, filt, None, None, None,
                        {"wcs_ok": wcs is not None, "exc": str(exc), "tb": tb,
                         "n_master": len(master_df),
                         "elapsed_s": time.time() - t0})

            (worker_status_cb(slot, fname, "Saving", 85) if worker_status_cb else None)
            out_path = out_dir / f"photometry_{fname}.tsv"
            center_stats = _frame_centering_stats(
                phot_df,
                fname=fname,
                filt=filt,
                status="ok",
                fwhm_px=fwhm_px,
                wcs_ok=wcs is not None,
                apcorr=apcorr_val,
                max_shift_px=_to_float(getattr(P, "max_recenter_shift", 2.0), 2.0),
                outlier_px=_to_float(getattr(P, "centroid_outlier_px", 1.0), 1.0),
                out_path=str(out_path),
            )
            try:
                phot_df.insert(0, "file", fname)
                phot_df.insert(1, "filter", filt)
                phot_df.insert(2, "FILTER", filt)
                phot_df.to_csv(out_path, sep="\t", index=False, float_format="%.6f")
            except Exception as exc:
                _log(f"[FORCED] [{fname}] write failed: {exc}")
            t_save = time.time()

            (worker_status_cb(slot, fname, "Done", 100) if worker_status_cb else None)
            return ("ok", fname, filt, phot_df, apcorr_val, gc_data,
                    {"wcs_ok": wcs is not None, "fwhm_px": fwhm_px,
                     "n_master": len(master_df), "out_path": str(out_path),
                     "noise": gc_data.get("noise", {}) if isinstance(gc_data, dict) else {},
                     "center_stats": center_stats,
                     "elapsed_s": t_save - t0,
                     "timing": {
                         "load_s": t_load - t0,
                         "wcs_fwhm_s": t_wcs - t_load,
                         "phot_apcorr_s": t_phot - t_wcs,
                         "save_s": t_save - t_phot,
                     }})
        finally:
            slot_queue.put(slot)

    tasks: List[Tuple[str, str, pd.DataFrame]] = []
    for filt, files in filter_map.items():
        master_df = master_catalogs.get(filt)
        if master_df is None or master_df.empty:
            _log(f"[FORCED] Skipping filter '{filt}' — no master catalog")
            with _results_lock:
                n_done[0] += len(files)
            continue
        for fname in files:
            tasks.append((filt, fname, master_df))

    _log(f"[FORCED] Starting parallel forced photometry with {max_workers} workers ({len(tasks)} frames)")

    executor = ThreadPoolExecutor(max_workers=max_workers)
    try:
        future_map = {
            executor.submit(_process_one, filt, fname, mdf): (filt, fname)
            for (filt, fname, mdf) in tasks
        }
        for future in as_completed(future_map):
            filt, fname = future_map[future]
            if _stop_requested():
                break

            with _results_lock:
                n_done[0] += 1
                cur = n_done[0]
            (progress_cb(cur, total, fname) if progress_cb else None)

            try:
                result = future.result()
            except Exception as exc:
                _log(f"[FORCED] [{fname}] task failed: {exc}\n{traceback.format_exc()}")
                with _results_lock:
                    index_rows.append({
                        "file": fname, "filter": filt, "status": "error",
                        "n_master": 0, "n_detected": 0, "n_forced": 0, "apcorr": np.nan,
                    })
                continue

            if result is None:
                continue

            status, fname_r, filt_r, phot_df, apcorr_val, gc_data, info = result

            if status == "no_image":
                _log(f"[FORCED] [{fname_r}] Image not found, skipping")
                with _results_lock:
                    index_rows.append({"file": fname_r, "filter": filt_r, "status": "no_image",
                                       "n_master": 0, "n_detected": 0, "n_forced": 0, "apcorr": np.nan})
                continue

            if status == "error":
                _log(f"[FORCED] [{fname_r}] phot_frame failed: {info['exc']}\n{info['tb']}")
                with _results_lock:
                    index_rows.append({"file": fname_r, "filter": filt_r, "status": "error",
                                       "n_master": info["n_master"], "n_detected": 0, "n_forced": 0,
                                       "wcs_ok": info["wcs_ok"], "apcorr": np.nan})
                continue

            # status == "ok"
            n_det    = int(phot_df["detected_flag"].sum())
            n_forced = int(phot_df["forced_flag"].sum())
            n_valid  = int(np.isfinite(phot_df["flux"]).sum())
            fwhm_px  = info["fwhm_px"]
            center_stats = info.get("center_stats") or {}
            elapsed = max(0.0, time.time() - run_t0)
            eta = ""
            if cur > 0 and total > cur:
                eta = f"  ETA={_fmt_duration((elapsed / float(cur)) * float(total - cur))}"
            timing = info.get("timing") or {}
            frame_dt = _fmt_duration(float(info.get("elapsed_s", 0.0) or 0.0))
            phot_dt = _fmt_duration(float(timing.get("phot_apcorr_s", 0.0) or 0.0))
            save_dt = _fmt_duration(float(timing.get("save_s", 0.0) or 0.0))

            _log(
                f"[FORCED] [{cur}/{total}] {fname_r} fwhm={fwhm_px:.1f}px  "
                f"det={n_det}  forced={n_forced}  "
                f"valid={n_valid}  apcorr={apcorr_val:.4f}  "
                f"dt={frame_dt} phot={phot_dt} save={save_dt}  "
                f"elapsed={_fmt_duration(elapsed)}{eta}"
            )

            if gc_data:
                gc_data["filter"] = filt_r
                gc_data["apcorr"] = apcorr_val
                gc_data["n_detected"] = n_det
                (apcorr_cb(gc_data) if apcorr_cb else None)

            with _results_lock:
                index_row = {
                    "file":         fname_r,
                    "filter":       filt_r,
                    "status":       "ok",
                    "n_master":     info["n_master"],
                    "n_detected":   n_det,
                    "n_forced":     n_forced,
                    "n_valid_phot": n_valid,
                    "fwhm_px":      fwhm_px,
                    "wcs_ok":       info["wcs_ok"],
                    "apcorr":       apcorr_val,
                    "path":         info["out_path"],
                }
                noise_info = info.get("noise") or {}
                if noise_info:
                    index_row.update(noise_info)
                for key, value in center_stats.items():
                    if key not in index_row:
                        index_row[key] = value
                index_rows.append(index_row)
                apcorr_rows.append({
                    "file":         fname_r,
                    "filter":       filt_r,
                    "apcorr":       apcorr_val,
                    "n_apcorr_stars": (
                        int(gc_data.get("n_apcorr_used", gc_data.get("n_stars", n_det)))
                        if isinstance(gc_data, dict) else n_det
                    ),
                    "n_apcorr_candidates": (
                        int(gc_data.get("n_apcorr_candidates", n_det))
                        if isinstance(gc_data, dict) else n_det
                    ),
                    "n_step4_quality_reject": (
                        int(gc_data.get("n_step4_quality_reject", 0))
                        if isinstance(gc_data, dict) else 0
                    ),
                    "n_step4_apcorr_reject": (
                        int(gc_data.get("n_step4_apcorr_reject", 0))
                        if isinstance(gc_data, dict) else 0
                    ),
                    "n_center_outlier_reject": (
                        int(gc_data.get("n_center_outlier_reject", 0))
                        if isinstance(gc_data, dict) else 0
                    ),
                })
                if center_stats:
                    center_stats_rows.append(center_stats)
                    (center_stats_cb(center_stats) if center_stats_cb else None)
    finally:
        executor.shutdown(wait=not _stop_requested(), cancel_futures=True)

    if _stop_requested():
        _log("[FORCED] Stopped by user.")
        return {
            "index_rows": index_rows,
            "apcorr_rows": apcorr_rows,
            "center_stats_rows": center_stats_rows,
        }

    index_write_ok = False
    if index_rows:
        idx_df = pd.DataFrame(index_rows)
        try:
            idx_df.to_csv(out_dir / "photometry_index.csv",  index=False, float_format="%.6f")
            idx_df.to_csv(out_dir / "frame_stats.csv", index=False, float_format="%.6f")
            index_write_ok = True
        except Exception as exc:
            _log(f"[FORCED] photometry_index write failed: {exc}")

    if apcorr_rows:
        try:
            pd.DataFrame(apcorr_rows).to_csv(
                out_dir / "apcorr_summary.csv", index=False, float_format="%.6f"
            )
        except Exception as exc:
            _log(f"[FORCED] apcorr_summary write failed: {exc}")

    if center_stats_rows:
        try:
            pd.DataFrame(center_stats_rows).to_csv(
                out_dir / "centering_stats.csv", index=False, float_format="%.6f"
            )
        except Exception as exc:
            _log(f"[FORCED] centering_stats write failed: {exc}")

    valid_tsv_names = {
        f"photometry_{str(r.get('file', ''))}.tsv"
        for r in index_rows
        if r.get("status") == "ok" and r.get("file")
    }
    stale_removed = 0
    if index_write_ok and index_rows:
        for old_tsv in out_dir.glob("photometry_*.tsv"):
            if old_tsv.name in valid_tsv_names:
                continue
            try:
                old_tsv.unlink()
                stale_removed += 1
            except Exception as exc:
                _log(f"[FORCED] stale output cleanup failed for {old_tsv.name}: {exc}")
        if stale_removed:
            _log(f"[FORCED] Removed {stale_removed} stale per-frame TSV outputs.")

    n_ok    = sum(1 for r in index_rows if r.get("status") == "ok")
    n_det_t = sum(r.get("n_detected", 0) for r in index_rows)
    n_frc_t = sum(r.get("n_forced",   0) for r in index_rows)
    _log(
        f"[FORCED] Done. {n_ok}/{total} frames OK. "
        f"Total detected={n_det_t} forced={n_frc_t} "
        f"elapsed={_fmt_duration(time.time() - run_t0)}"
    )
    return {
        "index_rows": index_rows,
        "apcorr_rows": apcorr_rows,
        "center_stats_rows": center_stats_rows,
    }


# ── Module-level helper (keeps _phot_frame readable) ───────────────────────────

def _growth_curve_fixed_sky(
    img: np.ndarray,
    positions: np.ndarray,
    radii: np.ndarray,
    r_in: float,
    r_out: float,
    phot_kw: dict,
    sky_median: Optional[np.ndarray] = None,
    sky_std: Optional[np.ndarray] = None,
    n_sky: Optional[np.ndarray] = None,
    status_cb=None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Growth-curve aperture fluxes with one annulus-stat pass.

    The annulus does not change while scanning aperture radii, so recomputing
    sigma-clipped sky statistics for every radius is the dominant Step7 cost.
    """
    from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats, aperture_photometry

    positions = np.asarray(positions, dtype=float)
    radii = np.asarray(radii, dtype=float)
    n_r = len(radii)
    n_s = len(positions)
    flux_matrix = np.full((n_r, n_s), np.nan, dtype=float)
    ferr_matrix = np.full((n_r, n_s), np.nan, dtype=float)
    if n_r == 0 or n_s == 0:
        return flux_matrix, ferr_matrix

    gain = _to_float(phot_kw.get("gain"), 1.0)
    rn_e = _to_float(phot_kw.get("rn_param_e"), 7.5)
    sky_frame_e = _to_float(phot_kw.get("sky_frame_e"), np.nan)
    sigma_clip_val = _to_float(phot_kw.get("sigma_clip_val"), 3.0)
    maxiters = _to_int(phot_kw.get("maxiters"), 5)
    sky_mode = str(phot_kw.get("sky_sigma_mode", "local") or "local").strip().lower()
    sky_incl_rn = bool(phot_kw.get("sky_sigma_includes_rn", True))
    min_n_sky = _to_int(phot_kw.get("min_n_sky_for_local"), 50)

    ann_area = float(np.pi * max(r_out * r_out - r_in * r_in, 1.0))
    reuse_sky = (
        sky_median is not None
        and len(np.asarray(sky_median)) == n_s
    )
    if reuse_sky:
        bkg_med = np.asarray(sky_median, dtype=float)
        bkg_std = (
            np.asarray(sky_std, dtype=float)
            if sky_std is not None and len(np.asarray(sky_std)) == n_s
            else np.full(n_s, np.nan, dtype=float)
        )
        n_sky_arr = (
            np.asarray(n_sky, dtype=float)
            if n_sky is not None and len(np.asarray(n_sky)) == n_s
            else np.full(n_s, ann_area, dtype=float)
        )
    else:
        sc = SigmaClip(sigma=sigma_clip_val, maxiters=maxiters)
        ann = CircularAnnulus(positions, r_in=r_in, r_out=r_out)
        ann_stats = ApertureStats(img, ann, sigma_clip=sc)
        bkg_med = np.asarray(ann_stats.median, dtype=float)
        bkg_std = np.asarray(ann_stats.std, dtype=float)
        n_sky_arr = np.asarray(ann_stats.sum_aper_area, dtype=float)

    bkg_med = np.where(np.isfinite(bkg_med), bkg_med, 0.0)
    bkg_std = np.where(np.isfinite(bkg_std), bkg_std, np.nan)
    n_sky_arr = np.where(np.isfinite(n_sky_arr) & (n_sky_arr > 0), n_sky_arr, ann_area)

    sigma_local_e2 = (np.maximum(bkg_std, 0.0) * gain) ** 2
    sigma_frame_e2 = float(sky_frame_e) ** 2 if np.isfinite(sky_frame_e) else np.nan
    use_local = np.isfinite(sigma_local_e2) & (bkg_std > 0) & (n_sky_arr >= min_n_sky)
    if sky_mode == "frame":
        sigma_pix_e2 = np.where(np.isfinite(sigma_frame_e2), sigma_frame_e2, sigma_local_e2)
    elif sky_mode == "max":
        sigma_pix_e2 = np.fmax(sigma_local_e2, sigma_frame_e2)
    else:
        sigma_pix_e2 = np.where(use_local, sigma_local_e2, sigma_frame_e2)
    sigma_pix_e2 = np.where(np.isfinite(sigma_pix_e2), sigma_pix_e2, 0.0)
    if sky_incl_rn:
        sigma_pix_e2 = np.maximum(sigma_pix_e2 - rn_e ** 2, 0.0)

    for ri, radius in enumerate(radii):
        if callable(status_cb):
            try:
                status_cb(ri + 1, n_r)
            except Exception:
                pass
        ap = CircularAperture(positions, r=float(radius))
        phot_tbl = aperture_photometry(img, ap, method="exact")
        ap_sum = np.asarray(phot_tbl["aperture_sum"], dtype=float)
        ap_area = float(ap.area)
        flux_net_adu = ap_sum - bkg_med * ap_area
        flux_e = flux_net_adu * gain

        var_source = np.maximum(flux_e, 0.0)
        var_bkg_ap = ap_area * sigma_pix_e2
        var_bkg_est = (ap_area ** 2 / np.maximum(n_sky_arr, 1)) * sigma_pix_e2
        var_rn = ap_area * rn_e ** 2 if sky_incl_rn else 0.0
        ferr_e = np.sqrt(np.maximum(var_source + var_bkg_ap + var_bkg_est + var_rn, 0.0))
        flux_matrix[ri, :] = flux_e
        ferr_matrix[ri, :] = ferr_e

    return flux_matrix, ferr_matrix


def _simple_apcorr(
    flux_arr: np.ndarray,
    apcorr_mask: np.ndarray,
    img: np.ndarray,
    positions: np.ndarray,
    r_ref: float,
    r_in: float,
    r_out: float,
    phot_kw: dict,
) -> float:
    """Fallback aperture correction: median(flux_ref / flux_ap) with sigma-clip."""
    flux_ref_arr, *_ = phot_vectorized(img, positions, r_ref, r_in, r_out, **phot_kw)
    apcorr_values = flux_ref_arr[apcorr_mask] / flux_arr[apcorr_mask]
    apcorr_values = apcorr_values[np.isfinite(apcorr_values) & (apcorr_values > 0)]
    if len(apcorr_values) == 0:
        return 1.0
    _, med, _ = sigma_clipped_stats(apcorr_values, sigma=3.0)
    return float(med) if np.isfinite(med) else 1.0
