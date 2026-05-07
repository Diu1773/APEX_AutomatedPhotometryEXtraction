"""
Step 7 (UI index 6): Forced Aperture Photometry

Master-catalog-driven forced aperture photometry.
For every frame: project master RA/Dec → pixel coords via WCS,
re-center on detected stars (±max_recenter_shift px), measure flux
at a fixed small aperture, compute aperture correction from a
growth curve measured on bright isolated stars (vectorized over N stars,
loop only over N_radii ≈ 14 steps).

Outputs (step7_forced_phot/):
  photometry_{fname}.tsv   — per-frame, per-source measurements
  photometry_index.csv     — summary row per frame
  apcorr_summary.csv       — per-frame aperture correction values
"""

from __future__ import annotations

import json
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

from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QGroupBox,
    QTextEdit, QProgressBar, QTableWidget, QTableWidgetItem,
    QHeaderView, QAbstractItemView, QWidget, QMessageBox,
    QTabWidget, QSplitter, QFormLayout, QCheckBox, QDoubleSpinBox,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal

try:
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
    from matplotlib.figure import Figure
    _HAVE_MPL = True
except ImportError:
    _HAVE_MPL = False

from .step_window_base import StepWindowBase
from .run_control import RunControlBar
from .log_panel import WorkflowLogWindow, WorkerStatusPanel, append_timestamped_log, show_raised
from apex.utils.step_paths import (
    step2_cropped_dir,
    step4_dir,
    step5_wcs_dir,
    step6_refbuild_dir,
    step7_forced_phot_dir,
    crop_is_active,
)
from apex.utils.photometry_utils import (
    phot_vectorized,
    refine_local_centroid,
)
from apex.utils.cache_utils import astap_wcs_candidates, parse_astap_wcs_file
from apex.utils.constants import get_parallel_workers

MAG_ERR_COEFF = 2.5 / np.log(10)

_GC_N_STEPS = 14   # number of radii in the growth curve


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


def _catalog_series(df: pd.DataFrame, col: str, fallback) -> pd.Series:
    if col in df.columns:
        return df[col].reset_index(drop=True)
    return pd.Series(fallback)


def _normalize_filter_value(value) -> str:
    raw = str(value or "").strip().lower()
    if not raw:
        return ""
    import re
    tokens = [t for t in re.split(r"[^a-z0-9]+", raw) if t]
    for token in reversed(tokens or [raw]):
        if token in {"ha", "halpha", "h-alpha"}:
            return "ha"
        if token in {"u", "g", "r", "i", "z", "b", "v"}:
            return token
    if raw in {"u", "g", "r", "i", "z", "b", "v"}:
        return raw
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


def _safe_rate(num: int, den: int) -> float:
    return float(num) / float(den) if den > 0 else float("nan")


def _step4_quality_check(row: pd.Series, P, sat_adu: float, datamax_adu: float) -> Tuple[bool, str, bool]:
    """Reuse Step4 detection-quality columns when they are available."""
    reasons: List[str] = []
    used = False

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
    x_pred = _numeric_col(phot_df, "x_pred")
    y_pred = _numeric_col(phot_df, "y_pred")
    on_frame = (~off_frame) & x_pred.notna() & y_pred.notna()

    match_offset = _numeric_col(phot_df, "match_offset_px")
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

    det_match_offset = match_offset[detected & match_offset.notna()]
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
        "detected_rate": detected_rate,
        "forced_rate": forced_rate,
        "recentered_rate": recentered_rate,
        "centroid_outlier_rate": outlier_rate,
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


# ── ForcedPhotWorker ───────────────────────────────────────────────────────────

class ForcedPhotWorker(QThread):
    """Per-frame forced aperture photometry worker."""

    progress      = pyqtSignal(int, int, str)   # current, total, fname
    log           = pyqtSignal(str)
    apcorr_update = pyqtSignal(dict)            # emitted per-frame with gc_data
    center_stats_update = pyqtSignal(dict)       # emitted per-frame with centering diagnostics
    worker_status = pyqtSignal(int, str, str, int)  # worker_id, fname, status, progress(0-100)
    finished      = pyqtSignal(dict)
    error         = pyqtSignal(str, str)

    def __init__(
        self,
        params,
        data_dir: Path,
        result_dir: Path,
        cache_dir: Path,
        file_list: List[str],
        output_dir: Optional[Path] = None,
    ):
        super().__init__()
        self.params = params
        self.data_dir = Path(data_dir)
        self.result_dir = Path(result_dir)
        self.cache_dir = Path(cache_dir)
        self.file_list = list(file_list)
        self.output_dir = Path(output_dir) if output_dir is not None else step7_forced_phot_dir(result_dir)
        self._stop_requested = False
        self._wcs_header_cache: Dict[str, fits.Header] = {}
        self._wcs_cache_lock = Lock()
        self._results_lock = Lock()
        self.max_workers = get_parallel_workers(params)

    def stop(self):
        self._stop_requested = True

    def _log(self, msg: str):
        self.log.emit(msg)

    # ── Path resolution ────────────────────────────────────────────────────────

    def _resolve_fits_path(self, fname: str) -> Optional[Path]:
        if crop_is_active(self.result_dir):
            p = step2_cropped_dir(self.result_dir) / fname
            if p.exists():
                return p
        try:
            orig = Path(self.params.get_file_path(fname))
            if orig.exists():
                return orig
        except Exception:
            pass
        cand = self.data_dir / fname
        return cand if cand.exists() else None

    def _load_wcs_for_frame(self, fname: str) -> Optional[WCS]:
        fits_path = self._resolve_fits_path(fname)
        candidates: List[Path] = []
        if fits_path is not None and fits_path.exists():
            candidates.append(fits_path)
        try:
            orig = Path(self.params.get_file_path(fname))
            if orig not in candidates and orig.exists():
                candidates.append(orig)
        except Exception:
            pass

        for path in candidates:
            try:
                key = str(path)
                with self._wcs_cache_lock:
                    hdr = self._wcs_header_cache.get(key)
                if hdr is None:
                    hdr = fits.getheader(path)
                    with self._wcs_cache_lock:
                        self._wcs_header_cache[key] = hdr
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

    def _load_image(self, fname: str) -> Optional[np.ndarray]:
        path = self._resolve_fits_path(fname)
        if path is None:
            return None
        try:
            with fits.open(path) as hdul:
                for hdu in hdul:
                    if hdu.data is not None and hdu.data.ndim == 2:
                        return hdu.data.astype(float)
        except Exception:
            pass
        return None

    def _load_fwhm(self, fname: str, fallback: float = 6.0) -> float:
        for cand in [
            self.cache_dir / f"detect_{fname}.json",
            step4_dir(self.result_dir) / f"detect_{fname}.json",
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

    def _load_sky_stats(self, fname: str) -> Tuple[float, float]:
        for cand in [
            self.cache_dir / f"detect_{fname}.json",
            step4_dir(self.result_dir) / f"detect_{fname}.json",
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

    def _load_detect_positions(self, fname: str) -> Optional[pd.DataFrame]:
        for cand in [
            self.cache_dir / f"detect_{fname}.csv",
            step4_dir(self.result_dir) / f"detect_{fname}.csv",
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
                    "fwhm_px", "sep_flag",
                ):
                    if col in df.columns:
                        out[col] = pd.to_numeric(df[col], errors="coerce")
                for col in ("fwhm_status", "source_type"):
                    if col in df.columns:
                        out[col] = df[col].astype(str).reset_index(drop=True)
                return out.dropna(subset=["x", "y"]).reset_index(drop=True)
            except Exception:
                continue
        return None

    def _load_master_catalog(self, filt: str) -> Optional[pd.DataFrame]:
        refbuild_dir = step6_refbuild_dir(self.result_dir)
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
                self._log(f"[FORCED] catalog load failed ({path.name}): {exc}")
        return None

    # ── Per-frame photometry ───────────────────────────────────────────────────

    def _phot_frame(
        self,
        fname: str,
        master_df: pd.DataFrame,
        img: np.ndarray,
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
        gain         = _to_float(getattr(P, "gain_e_per_adu",         1.0), 1.0)
        rn_e         = _to_float(getattr(P, "rdnoise_e",              7.5), 7.5)
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

        sky_med_frame, sky_sigma_frame = self._load_sky_stats(fname)
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
                    self._log(f"[FORCED][{fname}] WCS project failed: {exc}")

        if has_xy:
            x_ref_col = pd.to_numeric(master_df["x_ref"], errors="coerce").to_numpy(float)
            y_ref_col = pd.to_numeric(master_df["y_ref"], errors="coerce").to_numpy(float)
            missing = ~np.isfinite(x_pred)
            x_pred[missing] = x_ref_col[missing]
            y_pred[missing] = y_ref_col[missing]

        # --- Match to detections for recentering ---
        _status("Match detections", 34)
        det_df = self._load_detect_positions(fname)
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

        if det_df is not None and len(det_df) > 0:
            det_xy = det_df[["x", "y"]].to_numpy(float)
            master_xy_pred = np.column_stack([x_pred, y_pred])
            valid_pred = np.isfinite(x_pred) & np.isfinite(y_pred)
            if valid_pred.any():
                tree = KDTree(det_xy)
                dists, idxs = tree.query(master_xy_pred[valid_pred], k=1)
                match_mask = dists <= max_shift
                valid_indices = np.where(valid_pred)[0]
                for k, (vi, matched, idx, dist) in enumerate(zip(valid_indices, match_mask, idxs, dists)):
                    if matched:
                        detected_flag[vi] = True
                        det_uid_arr[vi] = int(det_df["det_uid"].iloc[idx])
                        dx = float(det_df["x"].iloc[idx]) - float(x_pred[vi])
                        dy = float(det_df["y"].iloc[idx]) - float(y_pred[vi])
                        det_x_arr[vi] = float(det_df["x"].iloc[idx])
                        det_y_arr[vi] = float(det_df["y"].iloc[idx])
                        match_dx[vi] = dx
                        match_dy[vi] = dy
                        match_offset[vi] = float(dist)
                        ok, reason, used = _step4_quality_check(
                            det_df.iloc[idx], P, sat_adu, datamax_adu
                        )
                        step4_quality_ok[vi] = ok
                        step4_quality_used[vi] = used
                        step4_quality_reason[vi] = reason

        # --- Recenter detected sources ---
        _status(f"Recenter {int(detected_flag.sum())}", 40)
        x_fit = x_pred.copy()
        y_fit = y_pred.copy()
        centroid_shift = np.full(n, np.nan)
        recentered_flag = np.zeros(n, dtype=bool)
        recenter_method = np.full(n, "forced_wcs", dtype=object)

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
                xi, yi = x_pred[i], y_pred[i]
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
                    and np.isfinite(match_offset[i])
                    and match_offset[i] <= max_shift
                ):
                    x_fit[i] = det_x_arr[i]
                    y_fit[i] = det_y_arr[i]
                    centroid_shift[i] = match_offset[i]
                    recentered_flag[i] = True
                    recenter_method[i] = "detected_xy"
                    accepted = True

                if not accepted:
                    recenter_method[i] = "failed"

        center_error = np.where(np.isfinite(centroid_shift), centroid_shift, match_offset)
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
        crowding_flag = (
            master_df["crowding_flag"].to_numpy(bool)
            if "crowding_flag" in master_df.columns
            else np.zeros(n, dtype=bool)
        )
        apcorr_mask = (
            detected_flag &
            ~crowding_flag &
            step4_quality_ok &
            ~centroid_outlier &
            np.isfinite(flux_arr) & (flux_arr > 0) &
            np.isfinite(snr_arr) & (snr_arr >= apcorr_min_snr) &
            ~is_sat_arr &
            ~is_nl_arr
        )
        n_step4_quality_reject = int((detected_flag & ~step4_quality_ok).sum())
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
            self._log(
                f"[FORCED][{fname}] apcorr refs capped "
                f"{n_apcorr_candidates}->{int(apcorr_mask.sum())} "
                f"(apcorr_max_sources={apcorr_max_sources})"
            )
        gc_positions = positions[apcorr_mask]
        n_gc_stars = int(apcorr_mask.sum())
        _status(f"Apcorr refs {n_gc_stars}/{n_apcorr_candidates}", 68)

        apcorr = 1.0
        growth_curve_data: dict = {}

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
                    "radii_px":       gc_radii.tolist(),
                    "enclosed_frac":  [float(v) for v in gc_curve],
                    "mag_err":        [float(v) for v in med_mag_err],
                    "n_stars":        int(valid_ref.sum()),
                    "n_apcorr_candidates": int(n_apcorr_candidates),
                    "n_apcorr_used":  int(n_gc_stars),
                    "n_step4_quality_reject": int(n_step4_quality_reject),
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

        mag_inst = np.where(
            np.isfinite(flux_corr) & (flux_corr > 0),
            -2.5 * np.log10(np.where(flux_corr > 0, flux_corr, np.nan)),
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
            ~np.isfinite(x_pred) | ~np.isfinite(y_pred) |
            (x_pred < 0) | (x_pred >= w) |
            (y_pred < 0) | (y_pred >= h)
        )
        centering_quality = np.full(n, "forced_only", dtype=object)
        centering_quality[off_frame] = "off_frame"
        centering_quality[detected_flag & ~centroid_outlier] = "ok"
        centering_quality[centroid_outlier] = "warn_shift"
        centering_quality[detected_flag & ~recenter_enabled] = "matched_no_recenter"

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
            "x_fit":             x_fit,
            "y_fit":             y_fit,
            "det_x":             det_x_arr,
            "det_y":             det_y_arr,
            "match_dx_px":       match_dx,
            "match_dy_px":       match_dy,
            "match_offset_px":   match_offset,
            "detected_flag":     detected_flag,
            "forced_flag":       ~detected_flag,
            "centroid_shift_px": centroid_shift,
            "center_error_px":   center_error,
            "recentered_flag":   recentered_flag,
            "recenter_method":   recenter_method,
            "step4_quality_ok":  step4_quality_ok,
            "step4_quality_used": step4_quality_used,
            "step4_quality_reason": step4_quality_reason,
            "centroid_outlier":  centroid_outlier,
            "centering_quality": centering_quality,
            "off_frame_flag":    off_frame,
            "flux":              flux_corr,
            "flux_e":            flux_corr,
            "flux_net_adu":      flux_corr / max(gain, 1e-12),
            "flux_err":          flux_err_corr,
            "mag_inst":          mag_inst,
            "mag_err":           mag_err,
            "snr":               snr_arr,
            "sky":               sky_arr,
            "apcorr":            apcorr,
            "is_saturated":      is_sat_arr,
            "is_nonlinear":      is_nl_arr,
            "bad_phot_flag":     bad_phot,
        })
        for col in ("ra_deg", "dec_deg", "gaia_source_id"):
            if col in master_df.columns and col not in out.columns:
                out[col] = master_df[col].reset_index(drop=True)

        return out, apcorr, growth_curve_data

    # ── Main run ──────────────────────────────────────────────────────────────

    def run(self):
        try:
            self._run_impl()
        except Exception as exc:
            self._log(f"[FORCED] Unhandled error: {exc}\n{traceback.format_exc()}")
            self.error.emit("ForcedPhot", str(exc))
            self.finished.emit({})

    def _run_impl(self):  # noqa: C901
        P = self.params.P
        out_dir = self.output_dir
        out_dir.mkdir(parents=True, exist_ok=True)

        fwhm_guess = _to_float(getattr(P, "fwhm_pix_guess", 6.0), 6.0)

        def _get_filter(fname: str) -> str:
            fits_path = self._resolve_fits_path(fname)
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
            return m.group(1).lower() if m else "unknown"

        filter_map: Dict[str, List[str]] = {}
        for fname in self.file_list:
            filt = _get_filter(fname)
            filter_map.setdefault(filt, []).append(fname)

        try:
            (out_dir / "filter_frames.json").write_text(
                json.dumps(filter_map, indent=2, ensure_ascii=False), encoding="utf-8"
            )
        except Exception as exc:
            self._log(f"[FORCED] filter_frames write failed: {exc}")

        master_catalogs: Dict[str, Optional[pd.DataFrame]] = {}
        for filt in filter_map:
            cat = self._load_master_catalog(filt)
            if cat is None:
                cat = self._load_master_catalog("r") or self._load_master_catalog("unknown")
            master_catalogs[filt] = cat
            if cat is not None:
                self._log(f"[FORCED] Loaded master catalog for filter '{filt}': {len(cat)} sources")
            else:
                self._log(f"[FORCED] WARNING: No master catalog found for filter '{filt}'")

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
            self._log(f"[FORCED] master_sources write failed: {exc}")

        total = len(self.file_list)
        index_rows: List[dict] = []
        apcorr_rows: List[dict] = []
        center_stats_rows: List[dict] = []
        n_done = [0]

        max_workers = max(1, int(self.max_workers))
        from queue import Queue as _Queue
        slot_queue: _Queue = _Queue()
        for _i in range(max_workers):
            slot_queue.put(_i)

        def _process_one(filt: str, fname: str, master_df: pd.DataFrame):
            if self._stop_requested:
                return None

            slot = slot_queue.get()
            try:
                self.worker_status.emit(slot, fname, "Loading", 5)
                img = self._load_image(fname)
                if img is None:
                    self.worker_status.emit(slot, fname, "No image", 100)
                    return ("no_image", fname, filt, None, None, None, None)

                self.worker_status.emit(slot, fname, "WCS / FWHM", 20)
                wcs = self._load_wcs_for_frame(fname)
                fwhm_px = self._load_fwhm(fname, fallback=fwhm_guess)

                self.worker_status.emit(slot, fname, "Photometry", 45)
                try:
                    def _status(label: str, pct: int):
                        self.worker_status.emit(slot, fname, label, pct)

                    phot_df, apcorr_val, gc_data = self._phot_frame(
                        fname, master_df, img, wcs, fwhm_px, P, status_cb=_status
                    )
                except Exception as exc:
                    tb = traceback.format_exc()
                    self.worker_status.emit(slot, fname, "Error", 100)
                    return ("error", fname, filt, None, None, None,
                            {"wcs_ok": wcs is not None, "exc": str(exc), "tb": tb,
                             "n_master": len(master_df)})

                self.worker_status.emit(slot, fname, "Saving", 85)
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
                    self._log(f"[FORCED] [{fname}] write failed: {exc}")

                self.worker_status.emit(slot, fname, "Done", 100)
                return ("ok", fname, filt, phot_df, apcorr_val, gc_data,
                        {"wcs_ok": wcs is not None, "fwhm_px": fwhm_px,
                         "n_master": len(master_df), "out_path": str(out_path),
                         "center_stats": center_stats})
            finally:
                slot_queue.put(slot)

        tasks: List[Tuple[str, str, pd.DataFrame]] = []
        for filt, files in filter_map.items():
            master_df = master_catalogs.get(filt)
            if master_df is None or master_df.empty:
                self._log(f"[FORCED] Skipping filter '{filt}' — no master catalog")
                with self._results_lock:
                    n_done[0] += len(files)
                continue
            for fname in files:
                tasks.append((filt, fname, master_df))

        self._log(f"[FORCED] Starting parallel forced photometry with {max_workers} workers ({len(tasks)} frames)")

        executor = ThreadPoolExecutor(max_workers=max_workers)
        try:
            future_map = {
                executor.submit(_process_one, filt, fname, mdf): (filt, fname)
                for (filt, fname, mdf) in tasks
            }
            for future in as_completed(future_map):
                filt, fname = future_map[future]
                if self._stop_requested:
                    break

                with self._results_lock:
                    n_done[0] += 1
                    cur = n_done[0]
                self.progress.emit(cur, total, fname)

                try:
                    result = future.result()
                except Exception as exc:
                    self._log(f"[FORCED] [{fname}] task failed: {exc}\n{traceback.format_exc()}")
                    with self._results_lock:
                        index_rows.append({
                            "file": fname, "filter": filt, "status": "error",
                            "n_master": 0, "n_detected": 0, "n_forced": 0, "apcorr": np.nan,
                        })
                    continue

                if result is None:
                    continue

                status, fname_r, filt_r, phot_df, apcorr_val, gc_data, info = result

                if status == "no_image":
                    self._log(f"[FORCED] [{fname_r}] Image not found, skipping")
                    with self._results_lock:
                        index_rows.append({"file": fname_r, "filter": filt_r, "status": "no_image",
                                           "n_master": 0, "n_detected": 0, "n_forced": 0, "apcorr": np.nan})
                    continue

                if status == "error":
                    self._log(f"[FORCED] [{fname_r}] phot_frame failed: {info['exc']}\n{info['tb']}")
                    with self._results_lock:
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

                self._log(
                    f"[FORCED] [{cur}/{total}] {fname_r} fwhm={fwhm_px:.1f}px  "
                    f"det={n_det}  forced={n_forced}  "
                    f"valid={n_valid}  apcorr={apcorr_val:.4f}"
                )

                if gc_data:
                    gc_data["filter"] = filt_r
                    gc_data["apcorr"] = apcorr_val
                    gc_data["n_detected"] = n_det
                    self.apcorr_update.emit(gc_data)

                with self._results_lock:
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
                        "n_center_outlier_reject": (
                            int(gc_data.get("n_center_outlier_reject", 0))
                            if isinstance(gc_data, dict) else 0
                        ),
                    })
                    if center_stats:
                        center_stats_rows.append(center_stats)
                        self.center_stats_update.emit(center_stats)
        finally:
            executor.shutdown(wait=not self._stop_requested, cancel_futures=True)

        if self._stop_requested:
            self._log("[FORCED] Stopped by user.")
            self.finished.emit({
                "index_rows": index_rows,
                "apcorr_rows": apcorr_rows,
                "center_stats_rows": center_stats_rows,
            })
            return

        if index_rows:
            idx_df = pd.DataFrame(index_rows)
            try:
                idx_df.to_csv(out_dir / "photometry_index.csv",  index=False, float_format="%.6f")
                idx_df.to_csv(out_dir / "frame_stats.csv", index=False, float_format="%.6f")
            except Exception as exc:
                self._log(f"[FORCED] photometry_index write failed: {exc}")

        if apcorr_rows:
            try:
                pd.DataFrame(apcorr_rows).to_csv(
                    out_dir / "apcorr_summary.csv", index=False, float_format="%.6f"
                )
            except Exception as exc:
                self._log(f"[FORCED] apcorr_summary write failed: {exc}")

        if center_stats_rows:
            try:
                pd.DataFrame(center_stats_rows).to_csv(
                    out_dir / "centering_stats.csv", index=False, float_format="%.6f"
                )
            except Exception as exc:
                self._log(f"[FORCED] centering_stats write failed: {exc}")

        n_ok    = sum(1 for r in index_rows if r.get("status") == "ok")
        n_det_t = sum(r.get("n_detected", 0) for r in index_rows)
        n_frc_t = sum(r.get("n_forced",   0) for r in index_rows)
        self._log(
            f"[FORCED] Done. {n_ok}/{total} frames OK. "
            f"Total detected={n_det_t} forced={n_frc_t}"
        )
        self.finished.emit({
            "index_rows": index_rows,
            "apcorr_rows": apcorr_rows,
            "center_stats_rows": center_stats_rows,
        })


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


# ── ForcedPhotWindow ───────────────────────────────────────────────────────────

class ForcedPhotWindow(StepWindowBase):
    """Step 7: Forced Aperture Photometry (master-driven)."""

    def __init__(self, params, file_manager, project_state, main_window):
        self.file_manager = file_manager
        self.worker = None
        self.results: dict = {}
        self._gc_accumulator: List[dict] = []   # growth curves collected across frames
        self._gc_per_frame: Dict[str, dict] = {}   # fname → growth curve dict
        self._center_stats_rows: List[dict] = []
        self._stats_seen = False

        super().__init__(
            step_index=6,
            step_name="Forced Aperture Phot",
            params=params,
            project_state=project_state,
            main_window=main_window,
        )

        self.setup_step_ui()
        self.restore_state()

    def setup_step_ui(self):
        info = QLabel(
            "Project master catalog positions onto each frame via WCS and measure forced "
            "aperture photometry.\nAperture correction is derived from a growth curve "
            f"measured on bright isolated stars ({_GC_N_STEPS} radii, vectorized)."
        )
        info.setStyleSheet("QLabel { background-color: #E8F5E9; padding: 10px; border-radius: 5px; }")
        info.setWordWrap(True)
        self.content_layout.addWidget(info)

        # Prerequisites
        status_group = QGroupBox("Prerequisites")
        status_layout = QVBoxLayout(status_group)
        self.status_label = QLabel("Checking...")
        status_layout.addWidget(self.status_label)
        self.content_layout.addWidget(status_group)
        self._check_prerequisites()

        # Controls
        ctrl_layout = QHBoxLayout()
        self.run_bar = RunControlBar(
            "Run Forced Photometry", "Show Log",
            run_cb=self.run_forced_phot,
            stop_cb=self.stop_forced_phot,
            log_cb=self._show_log,
        )
        ctrl_layout.addWidget(self.run_bar)
        self.btn_run  = self.run_bar.btn_run
        self.btn_stop = self.run_bar.btn_stop
        self.content_layout.addLayout(ctrl_layout)

        center_group = QGroupBox("Centering / Recenter")
        center_form = QFormLayout(center_group)
        self.chk_recenter = QCheckBox("Use detected-source recentering")
        self.chk_recenter.setChecked(bool(getattr(self.params.P, "recenter_aperture", True)))
        self.chk_recenter.setToolTip(
            "When a projected master position matches a Step4 detection, seed local centroiding "
            "from that detection before aperture photometry."
        )
        center_form.addRow("Mode", self.chk_recenter)

        self.spin_max_shift = QDoubleSpinBox()
        self.spin_max_shift.setRange(0.10, 20.0)
        self.spin_max_shift.setDecimals(2)
        self.spin_max_shift.setSingleStep(0.10)
        self.spin_max_shift.setSuffix(" px")
        self.spin_max_shift.setValue(float(getattr(self.params.P, "max_recenter_shift", 2.0)))
        self.spin_max_shift.setToolTip(
            "Maximum allowed distance between projected master position and matched detection. "
            "Too large can jump to a neighbor in crowded fields."
        )
        center_form.addRow("Match / recenter limit", self.spin_max_shift)

        self.spin_outlier_px = QDoubleSpinBox()
        self.spin_outlier_px.setRange(0.05, 10.0)
        self.spin_outlier_px.setDecimals(2)
        self.spin_outlier_px.setSingleStep(0.05)
        self.spin_outlier_px.setSuffix(" px")
        self.spin_outlier_px.setValue(float(getattr(self.params.P, "centroid_outlier_px", 1.0)))
        self.spin_outlier_px.setToolTip(
            "Stats warning threshold for center_error_px. Sources above this value are "
            "marked centroid_outlier in the forced photometry TSV."
        )
        center_form.addRow("Outlier threshold", self.spin_outlier_px)

        self.center_param_label = QLabel()
        self.center_param_label.setWordWrap(True)
        center_form.addRow("Current", self.center_param_label)
        self.chk_recenter.stateChanged.connect(self._on_center_params_changed)
        self.spin_max_shift.valueChanged.connect(self._on_center_params_changed)
        self.spin_outlier_px.valueChanged.connect(self._on_center_params_changed)
        self._update_center_param_label()
        self.content_layout.addWidget(center_group)

        # ── Progress bar — placed ABOVE the tabs so it is always visible ──────
        prog_layout = QHBoxLayout()
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        prog_layout.addWidget(self.progress_bar)
        self.progress_label = QLabel("Ready")
        self.progress_label.setMinimumWidth(300)
        prog_layout.addWidget(self.progress_label)
        self.content_layout.addLayout(prog_layout)

        # ── Tab widget ─────────────────────────────────────────────────────────
        self.tabs = QTabWidget()
        self.content_layout.addWidget(self.tabs)

        # ── Tab 0: Centering Stats ─────────────────────────────────────────────
        tab_stats = QWidget()
        tab_stats_layout = QVBoxLayout(tab_stats)

        self.center_summary_label = QLabel("Run forced photometry to populate centering stats.")
        self.center_summary_label.setWordWrap(True)
        self.center_summary_label.setStyleSheet(
            "QLabel { background-color: #F5F5F5; padding: 8px; border-radius: 4px; }"
        )
        tab_stats_layout.addWidget(self.center_summary_label)

        center_plot_group = QGroupBox("Selected Frame Center Error")
        center_plot_layout = QVBoxLayout(center_plot_group)
        if _HAVE_MPL:
            self._center_fig = Figure(figsize=(5, 3.5), tight_layout=True)
            self._center_ax_hist = self._center_fig.add_subplot(121)
            self._center_ax_scatter = self._center_fig.add_subplot(122)
            self._center_canvas = FigureCanvasQTAgg(self._center_fig)
            self._center_canvas.setMinimumHeight(260)
            center_plot_layout.addWidget(self._center_canvas)
            self._init_center_axes()
        else:
            center_plot_layout.addWidget(QLabel("matplotlib not available — install it for center-error plots"))
            self._center_canvas = None
        tab_stats_layout.addWidget(center_plot_group)

        center_table_group = QGroupBox("Per-frame Centering Diagnostics")
        center_table_layout = QVBoxLayout(center_table_group)
        self.center_stats_table = QTableWidget()
        _center_headers = [
            "File", "Filter", "Status", "Centering",
            "Detected %", "Forced %", "Recentered %",
            "Match p50", "Match p90", "Center p50", "Center p90",
            "Outlier %", "Med mag_err", "High-shift Δmag_err", "Advice",
        ]
        self.center_stats_table.setColumnCount(len(_center_headers))
        self.center_stats_table.setHorizontalHeaderLabels(_center_headers)
        _center_tooltips = [
            "FITS 파일명",
            "필터명",
            "처리 상태",
            "Centering health: OK / REVIEW / CHECK / LOW_MATCH",
            "Projected master positions matched to Step4 detections, divided by on-frame sources",
            "On-frame sources without a matched detection; their center remains WCS/master-driven",
            "Matched detections that actually moved to a refined/detected center",
            "Median projected-position to detection offset (px)",
            "90th percentile projected-position to detection offset (px)",
            "Median final center error proxy (px)",
            "90th percentile final center error proxy (px)",
            "Fraction of tested centers above centroid_outlier_px",
            "Median reported photometric magnitude error",
            "Median mag_err(outliers) - median mag_err(non-outlier detections)",
            "Suggested inspection target",
        ]
        for col, tip in enumerate(_center_tooltips):
            item = self.center_stats_table.horizontalHeaderItem(col)
            if item is not None:
                item.setToolTip(tip)
        self.center_stats_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.center_stats_table.horizontalHeader().setStretchLastSection(True)
        self.center_stats_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.center_stats_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.center_stats_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.center_stats_table.itemSelectionChanged.connect(self._on_center_stats_row_selected)
        center_table_layout.addWidget(self.center_stats_table)
        tab_stats_layout.addWidget(center_table_group)

        self._stats_tab_index = self.tabs.addTab(tab_stats, "Stats")

        # ── Tab 1: Apcorr / Growth Curve ──────────────────────────────────────
        tab0 = QWidget()
        tab0_layout = QVBoxLayout(tab0)

        # Growth curve plot
        gc_group = QGroupBox("Growth Curve (selected frame)")
        gc_layout = QVBoxLayout(gc_group)
        if _HAVE_MPL:
            self._gc_fig    = Figure(figsize=(5, 4.5), tight_layout=True)
            self._gc_ax_mag = self._gc_fig.add_subplot(211)
            self._gc_ax_err = self._gc_fig.add_subplot(212)
            self._gc_canvas = FigureCanvasQTAgg(self._gc_fig)
            self._gc_canvas.setMinimumHeight(320)
            gc_layout.addWidget(self._gc_canvas)
            self._init_gc_axes()
        else:
            gc_layout.addWidget(QLabel("matplotlib not available — install it for the growth curve plot"))
            self._gc_canvas = None
        tab0_layout.addWidget(gc_group)

        # Per-frame apcorr table (below the plot)
        ap_group = QGroupBox("Per-frame Aperture Correction")
        ap_layout = QVBoxLayout(ap_group)
        self.apcorr_table = QTableWidget()
        self.apcorr_table.setColumnCount(10)
        _ap_headers = [
            "File", "Filter", "Apcorr", "N stars",
            "Candidates", "Step4 reject", "Center reject",
            "r_opt (px)", "r_opt/FWHM", "min mag_err",
        ]
        self.apcorr_table.setHorizontalHeaderLabels(_ap_headers)
        _ap_tooltips = [
            "FITS 파일명",
            "필터명",
            "Aperture Correction = 1 / enclosed_fraction(at r_ap).\n"
            "측정 flux에 곱하는 무차원 보정 계수 (보통 1.05~1.3).",
            "Apcorr 산정에 사용된 isolated bright star 개수",
            "Apcorr 조건을 통과한 후보 개수. 계산은 apcorr_max_sources까지만 사용.",
            "Step4 detect 품질 컬럼(shape/peak/status)으로 제외된 검출 source 수",
            "centroid_outlier_px보다 중심 오차가 커서 apcorr reference에서 제외된 source 수",
            "Optimal aperture 반지름 (px): U-shape mag_err 곡선의 최저점 = SNR 최대 구경.",
            "최적 구경 배수 = r_opt / FWHM. 현재 r_ap/FWHM과 다르면 forced_r_ap_scale 튜닝 고려.",
            "최저점 median mag_err.",
        ]
        for col, tip in enumerate(_ap_tooltips):
            item = self.apcorr_table.horizontalHeaderItem(col)
            if item is not None:
                item.setToolTip(tip)
        self.apcorr_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.apcorr_table.horizontalHeader().setStretchLastSection(True)
        self.apcorr_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.apcorr_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.apcorr_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.apcorr_table.setMaximumHeight(180)
        self.apcorr_table.itemSelectionChanged.connect(self._on_apcorr_row_selected)
        ap_layout.addWidget(self.apcorr_table)
        tab0_layout.addWidget(ap_group)

        self._apcorr_tab_index = self.tabs.addTab(tab0, "Apcorr")

        # ── Tab 2: Photometry Results ──────────────────────────────────────────
        tab1 = QWidget()
        tab1_layout = QVBoxLayout(tab1)

        self.results_table = QTableWidget()
        _res_headers = [
            "File", "Filter", "Status", "WCS",
            "FWHM (px)", "r_ap (px)", "r_ref (px)",
            "Apcorr", "N master", "N detected", "N forced", "N valid",
        ]
        self.results_table.setColumnCount(len(_res_headers))
        self.results_table.setHorizontalHeaderLabels(_res_headers)
        _res_tooltips = [
            "FITS 파일명",
            "필터명",
            "처리 상태 (ok / no_image / error)",
            "WCS solve 성공 여부",
            "프레임 FWHM (px) — step4 detection에서",
            "측광 aperture 반지름 (px). forced_r_ap_scale × FWHM",
            "Reference aperture 반지름 (px). forced_ref_ap_scale × FWHM",
            "Aperture correction = 1 / enclosed_fraction(at r_ap)",
            "Master catalog source 수",
            "프레임에서 검출된 source 수 (recenter용)",
            "Forced 측광 source 수",
            "유효 flux (finite) source 수",
        ]
        for col, tip in enumerate(_res_tooltips):
            it = self.results_table.horizontalHeaderItem(col)
            if it is not None:
                it.setToolTip(tip)
        self.results_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.results_table.horizontalHeader().setStretchLastSection(True)
        self.results_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        tab1_layout.addWidget(self.results_table)

        self.tabs.addTab(tab1, "Results")

        # ── Floating log window with Workers panel ─────────────────────────────
        _workers_group = QGroupBox("Workers")
        _workers_group.setMinimumWidth(430)
        _wg_layout = QVBoxLayout(_workers_group)
        _wg_layout.setContentsMargins(5, 5, 5, 5)
        self._worker_panel = WorkerStatusPanel(_workers_group)
        _wg_layout.addWidget(self._worker_panel)

        self._log_win = WorkflowLogWindow(
            self, "Forced Phot Log & Workers", width=900, height=500,
            side_widget=_workers_group,
        )
        self.log_text = self._log_win.log_text

    # ── Prerequisites ──────────────────────────────────────────────────────────

    def _check_prerequisites(self):
        refbuild_dir = step6_refbuild_dir(self.params.P.result_dir)
        wcs_dir      = step5_wcs_dir(self.params.P.result_dir)
        has_wcs = wcs_dir.exists() and any(wcs_dir.glob("wcs_solve_summary.csv"))
        has_cat = refbuild_dir.exists() and any(refbuild_dir.glob("ref_catalog*.tsv"))
        parts = []
        if not has_wcs:
            parts.append("WCS (step5_wcs/) not found")
        if not has_cat:
            parts.append("Master catalog (step6_refbuild/) not found")
        if parts:
            self.status_label.setText("Missing: " + "; ".join(parts))
            self.status_label.setStyleSheet("QLabel { color: #f44336; }")
        else:
            self.status_label.setText("Prerequisites OK — ready to run")
            self.status_label.setStyleSheet("QLabel { color: #4CAF50; }")

    # ── Centering controls ─────────────────────────────────────────────────────

    def _sync_centering_params(self):
        if not hasattr(self, "chk_recenter"):
            return
        self.params.P.recenter_aperture = bool(self.chk_recenter.isChecked())
        self.params.P.max_recenter_shift = float(self.spin_max_shift.value())
        self.params.P.centroid_outlier_px = float(self.spin_outlier_px.value())

    def _update_center_param_label(self):
        if not hasattr(self, "center_param_label"):
            return
        mode = "enabled" if self.chk_recenter.isChecked() else "disabled"
        self.center_param_label.setText(
            f"recenter={mode}  |  match_limit={self.spin_max_shift.value():.2f}px  "
            f"|  outlier={self.spin_outlier_px.value():.2f}px"
        )

    def _on_center_params_changed(self, *_):
        self._sync_centering_params()
        self._update_center_param_label()
        self.persist_params()

    # ── Load existing results ──────────────────────────────────────────────────

    def _try_load_existing_results(self):
        idx_path = step7_forced_phot_dir(self.params.P.result_dir) / "photometry_index.csv"
        if not idx_path.exists():
            return
        try:
            df = pd.read_csv(idx_path)
            self._populate_results_table(df.to_dict(orient="records"))
            self.progress_label.setText(f"Loaded {len(df)} frames from previous run")
            self.update_navigation_buttons()
        except Exception:
            pass

        apc_path = step7_forced_phot_dir(self.params.P.result_dir) / "apcorr_summary.csv"
        if apc_path.exists():
            try:
                apc_df = pd.read_csv(apc_path)
                for _, row in apc_df.iterrows():
                    ri = self.apcorr_table.rowCount()
                    self.apcorr_table.insertRow(ri)
                    for ci, val in enumerate([
                        str(row.get("file", "")),
                        str(row.get("filter", "")),
                        f"{row.get('apcorr', float('nan')):.4f}" if pd.notna(row.get("apcorr")) else "—",
                        str(row.get("n_apcorr_stars", "")),
                        str(row.get("n_apcorr_candidates", "")),
                        str(row.get("n_step4_quality_reject", "")),
                        str(row.get("n_center_outlier_reject", "")),
                        "—", "—", "—",
                    ]):
                        item = QTableWidgetItem(val)
                        item.setTextAlignment(Qt.AlignCenter)
                        self.apcorr_table.setItem(ri, ci, item)
            except Exception:
                pass

        stats_path = step7_forced_phot_dir(self.params.P.result_dir) / "centering_stats.csv"
        if stats_path.exists():
            try:
                stats_df = pd.read_csv(stats_path)
                rows = stats_df.to_dict(orient="records")
                self._populate_center_stats_table(rows)
                self._refresh_center_summary()
            except Exception:
                pass

    # ── Table helpers ──────────────────────────────────────────────────────────

    def _populate_results_table(self, rows: list):
        def _fmt_f(v, fmt="{:.2f}"):
            try:
                f = float(v)
                if not np.isfinite(f):
                    return "—"
                return fmt.format(f)
            except (TypeError, ValueError):
                return "—"

        self.results_table.setRowCount(len(rows))
        for ri, row in enumerate(rows):
            wcs_ok = row.get("wcs_ok")
            wcs_str = "—" if wcs_ok is None else ("OK" if wcs_ok else "FAIL")
            # r_ap / r_ref might not be in row dict; pull from gc_per_frame if available
            gc = self._gc_per_frame.get(str(row.get("file", "")), {})
            r_ap_val  = row.get("r_ap_px",  gc.get("r_ap_px"))
            r_ref_val = row.get("r_ref_px", gc.get("r_ref_px"))
            vals = [
                str(row.get("file", "")),
                str(row.get("filter", "")),
                str(row.get("status", "")),
                wcs_str,
                _fmt_f(row.get("fwhm_px"), "{:.2f}"),
                _fmt_f(r_ap_val, "{:.1f}"),
                _fmt_f(r_ref_val, "{:.1f}"),
                _fmt_f(row.get("apcorr"), "{:.4f}"),
                str(row.get("n_master", "")),
                str(row.get("n_detected", "")),
                str(row.get("n_forced", "")),
                str(row.get("n_valid_phot", "")),
            ]
            for ci, val in enumerate(vals):
                item = QTableWidgetItem(val)
                item.setTextAlignment(Qt.AlignCenter)
                self.results_table.setItem(ri, ci, item)

    def _fmt_table_float(self, v, fmt="{:.2f}") -> str:
        try:
            f = float(v)
            if not np.isfinite(f):
                return "—"
            return fmt.format(f)
        except (TypeError, ValueError):
            return "—"

    def _fmt_pct(self, v) -> str:
        try:
            f = float(v)
            if not np.isfinite(f):
                return "—"
            return f"{100.0 * f:.1f}%"
        except (TypeError, ValueError):
            return "—"

    def _populate_center_stats_table(self, rows: list):
        self._center_stats_rows = list(rows or [])
        self.center_stats_table.setRowCount(len(self._center_stats_rows))
        for ri, row in enumerate(self._center_stats_rows):
            self._set_center_stats_row(ri, row)
        if self._center_stats_rows and self.center_stats_table.rowCount() > 0:
            self.center_stats_table.selectRow(0)

    def _append_center_stats_row(self, row: dict):
        self._center_stats_rows.append(row)
        ri = self.center_stats_table.rowCount()
        self.center_stats_table.insertRow(ri)
        self._set_center_stats_row(ri, row)
        self.center_stats_table.scrollToBottom()

    def _set_center_stats_row(self, ri: int, row: dict):
        vals = [
            str(row.get("file", "")),
            str(row.get("filter", "")),
            str(row.get("status", "")),
            str(row.get("centering_status", "")),
            self._fmt_pct(row.get("detected_rate")),
            self._fmt_pct(row.get("forced_rate")),
            self._fmt_pct(row.get("recentered_rate")),
            self._fmt_table_float(row.get("match_offset_med_px"), "{:.3f}"),
            self._fmt_table_float(row.get("match_offset_p90_px"), "{:.3f}"),
            self._fmt_table_float(row.get("center_error_med_px"), "{:.3f}"),
            self._fmt_table_float(row.get("center_error_p90_px"), "{:.3f}"),
            self._fmt_pct(row.get("centroid_outlier_rate")),
            self._fmt_table_float(row.get("mag_err_med"), "{:.4f}"),
            self._fmt_table_float(row.get("mag_err_high_shift_delta"), "{:+.4f}"),
            str(row.get("advice", "")),
        ]
        for ci, val in enumerate(vals):
            item = QTableWidgetItem(val)
            item.setTextAlignment(Qt.AlignCenter if ci < len(vals) - 1 else Qt.AlignLeft | Qt.AlignVCenter)
            if ci == len(vals) - 1:
                item.setToolTip(str(row.get("advice", "")))
            self.center_stats_table.setItem(ri, ci, item)

    def _refresh_center_summary(self):
        rows = self._center_stats_rows
        if not rows:
            self.center_summary_label.setText("Run forced photometry to populate centering stats.")
            return
        statuses: Dict[str, int] = {}
        for row in rows:
            key = str(row.get("centering_status", "") or "UNKNOWN")
            statuses[key] = statuses.get(key, 0) + 1
        center_p90 = _finite_median([r.get("center_error_p90_px", np.nan) for r in rows])
        forced_med = _finite_median([r.get("forced_rate", np.nan) for r in rows])
        outlier_max = _finite_max([r.get("centroid_outlier_rate", np.nan) for r in rows])
        status_text = ", ".join(f"{k}:{v}" for k, v in sorted(statuses.items()))
        self.center_summary_label.setText(
            f"Frames={len(rows)}  |  {status_text}  |  "
            f"median center p90={self._fmt_table_float(center_p90, '{:.3f}')} px  |  "
            f"median forced={self._fmt_pct(forced_med)}  |  "
            f"max outlier={self._fmt_pct(outlier_max)}"
        )

    def _init_center_axes(self):
        ax_hist = self._center_ax_hist
        ax_scatter = self._center_ax_scatter
        ax_hist.set_title("Center error", fontsize=9)
        ax_hist.set_xlabel("center_error_px", fontsize=9)
        ax_hist.set_ylabel("N", fontsize=9)
        ax_hist.grid(True, alpha=0.3)
        ax_scatter.set_title("Error vs mag_err", fontsize=9)
        ax_scatter.set_xlabel("center_error_px", fontsize=9)
        ax_scatter.set_ylabel("mag_err", fontsize=9)
        ax_scatter.grid(True, alpha=0.3)

    def _update_center_plot(self, fname: str, row: Optional[dict] = None):
        if not _HAVE_MPL or self._center_canvas is None:
            return
        ax_hist = self._center_ax_hist
        ax_scatter = self._center_ax_scatter
        ax_hist.cla()
        ax_scatter.cla()
        self._init_center_axes()
        if not fname:
            self._center_canvas.draw_idle()
            return

        path = step7_forced_phot_dir(self.params.P.result_dir) / f"photometry_{fname}.tsv"
        if not path.exists():
            ax_hist.set_title(f"{fname} not found", fontsize=9)
            self._center_canvas.draw_idle()
            return

        try:
            df = pd.read_csv(path, sep="\t")
        except Exception as exc:
            ax_hist.set_title(f"Read failed: {exc}", fontsize=9)
            self._center_canvas.draw_idle()
            return

        center_error = _numeric_col(df, "center_error_px")
        mag_err = _numeric_col(df, "mag_err")
        outlier = _bool_col(df, "centroid_outlier")
        finite_center = center_error[np.isfinite(center_error)]
        if len(finite_center) > 0:
            ax_hist.hist(finite_center.to_numpy(float), bins=30, color="#1565C0", alpha=0.78)
        outlier_default = float(getattr(self.params.P, "centroid_outlier_px", 1.0))
        outlier_px = _to_float(row.get("centroid_outlier_px") if row else None, outlier_default)
        if np.isfinite(outlier_px):
            ax_hist.axvline(outlier_px, color="#E53935", ls="--", lw=1.2, label=f"outlier={outlier_px:.2f}px")
            ax_hist.legend(fontsize=7, frameon=False)

        finite_scatter = center_error.notna() & mag_err.notna()
        if finite_scatter.any():
            idx = np.where(finite_scatter.to_numpy())[0]
            if len(idx) > 5000:
                idx = np.linspace(0, len(idx) - 1, 5000).astype(int)
                idx = np.where(finite_scatter.to_numpy())[0][idx]
            x = center_error.iloc[idx].to_numpy(float)
            y = mag_err.iloc[idx].to_numpy(float)
            colors = np.where(outlier.iloc[idx].to_numpy(bool), "#E53935", "#2E7D32")
            ax_scatter.scatter(x, y, s=8, c=colors, alpha=0.55, linewidths=0)

        title = str(fname)
        if row:
            title += (
                f" | p90={self._fmt_table_float(row.get('center_error_p90_px'), '{:.3f}')}px"
                f" | outliers={self._fmt_pct(row.get('centroid_outlier_rate'))}"
            )
        ax_hist.set_title(title, fontsize=9)
        self._center_canvas.draw_idle()

    def _on_center_stats_row_selected(self):
        rows = self.center_stats_table.selectionModel().selectedRows() if self.center_stats_table.selectionModel() else []
        if not rows:
            return
        ri = rows[0].row()
        if ri < 0 or ri >= len(self._center_stats_rows):
            return
        row = self._center_stats_rows[ri]
        fname = str(row.get("file", ""))
        if fname:
            self._update_center_plot(fname, row)

    def _init_gc_axes(self):
        ax_mag = self._gc_ax_mag
        ax_err = self._gc_ax_err
        ax_mag.set_ylabel("Inst Magnitude", fontsize=9)
        ax_mag.set_title("Growth Curve", fontsize=9)
        ax_mag.grid(True, alpha=0.3)
        ax_err.set_xlabel("Aperture radius (px)", fontsize=9)
        ax_err.set_ylabel("Median mag_err", fontsize=9)
        ax_err.set_title("Error vs Aperture (U-shape)", fontsize=9)
        ax_err.grid(True, alpha=0.3)

    def _update_gc_plot(self, gc: dict):
        """Redraw the growth curve figure for a single frame."""
        if not _HAVE_MPL or self._gc_canvas is None:
            return
        ax_mag = self._gc_ax_mag
        ax_err = self._gc_ax_err
        ax_mag.cla()
        ax_err.cla()
        self._init_gc_axes()

        if not gc:
            self._gc_canvas.draw_idle()
            return

        fwhm  = float(gc.get("fwhm_px", 0.0) or 0.0)
        r_ap  = float(gc.get("r_ap_px",  0.0))
        r_ref = float(gc.get("r_ref_px", 0.0))
        r_opt = float(gc.get("r_opt_px", 0.0) or 0.0)
        fname = str(gc.get("fname", ""))
        apcorr_val = float(gc.get("apcorr", np.nan))

        radii = gc.get("radii_px", [])
        encs  = gc.get("enclosed_frac", [])
        errs  = gc.get("mag_err", [])

        if radii and encs:
            arr = np.asarray(encs, dtype=float)
            arr = np.where((arr > 0) & np.isfinite(arr), arr, np.nan)
            mag = -2.5 * np.log10(arr)
            ax_mag.plot(radii, mag, "-o", color="#1565C0",
                        lw=1.5, markersize=5,
                        markeredgecolor="white", markeredgewidth=0.5)
        if radii and errs:
            ax_err.plot(radii, errs, "-s", color="#E53935",
                        lw=1.8, markersize=5,
                        markeredgecolor="white", markeredgewidth=0.5)

        for ax in (ax_mag, ax_err):
            if r_opt > 0:
                ax.axvline(r_opt, color="#7B1FA2", lw=1.6, ls="-",
                           alpha=0.85, label=f"r_opt={r_opt:.1f}px")
            if r_ap > 0:
                ax.axvline(r_ap, color="#E53935", lw=1.2, ls="--",
                           alpha=0.8, label=f"r_ap={r_ap:.1f}px")
            if r_ref > 0:
                ax.axvline(r_ref, color="#43A047", lw=1.2, ls="--",
                           alpha=0.8, label=f"r_ref={r_ref:.1f}px")
            if fwhm > 0:
                ax.axvline(fwhm, color="#6D4C41", lw=1.0, ls=":",
                           alpha=0.8, label=f"FWHM={fwhm:.2f}px")

        ax_mag.invert_yaxis()
        title = fname
        if np.isfinite(apcorr_val):
            title += f"  |  apcorr={apcorr_val:.4f}"
        ax_mag.set_title(title, fontsize=9)
        ax_mag.legend(fontsize=7, frameon=False, loc="best")
        ax_err.legend(fontsize=7, frameon=False, loc="best")
        self._gc_canvas.draw_idle()

    def _on_apcorr_row_selected(self):
        rows = self.apcorr_table.selectionModel().selectedRows() if self.apcorr_table.selectionModel() else []
        if not rows:
            return
        ri = rows[0].row()
        item = self.apcorr_table.item(ri, 0)
        if item is None:
            return
        fname = item.text()
        gc = self._gc_per_frame.get(fname)
        if gc:
            self._update_gc_plot(gc)

    # ── Log window ─────────────────────────────────────────────────────────────

    def _show_log(self):
        show_raised(self._log_win)

    # ── Worker control ─────────────────────────────────────────────────────────

    def run_forced_phot(self):
        self._sync_centering_params()
        P = self.params.P
        try:
            result_dir = Path(P.result_dir)
            data_dir   = Path(P.data_dir)
            cache_dir  = Path(P.cache_dir)
        except Exception as exc:
            QMessageBox.critical(self, "Config Error", f"Cannot read paths from params: {exc}")
            return

        file_list = []
        try:
            file_list = list(self.file_manager.get_file_list())
        except Exception:
            pass
        if not file_list:
            QMessageBox.warning(self, "No Files", "File list is empty. Run File Selection first.")
            return

        refbuild_dir = step6_refbuild_dir(result_dir)
        if not any(refbuild_dir.glob("ref_catalog*.tsv")):
            QMessageBox.warning(
                self, "No Master Catalog",
                "Master catalog not found in step6_refbuild/. Run Master Catalog Build first."
            )
            return

        self.log_text.clear()
        self._gc_accumulator.clear()
        self._gc_per_frame.clear()
        self._center_stats_rows.clear()
        self._stats_seen = False
        self.apcorr_table.setRowCount(0)
        self.center_stats_table.setRowCount(0)
        self._update_gc_plot({})
        self._update_center_plot("")
        self._refresh_center_summary()
        if hasattr(self, "_worker_panel"):
            self._worker_panel.clear()

        self.worker = ForcedPhotWorker(
            params=self.params,
            data_dir=data_dir,
            result_dir=result_dir,
            cache_dir=cache_dir if cache_dir.is_absolute() else result_dir / cache_dir,
            file_list=file_list,
        )
        self.worker.progress.connect(self._on_progress)
        self.worker.log.connect(self._on_log)
        self.worker.apcorr_update.connect(self._on_apcorr_update)
        self.worker.center_stats_update.connect(self._on_center_stats_update)
        self.worker.worker_status.connect(self._worker_panel.update_worker)
        self.worker.finished.connect(self._on_finished)
        self.worker.error.connect(self._on_error)
        self.worker.start()

        self.run_bar.set_running(True)
        self.progress_bar.setValue(0)
        self.progress_label.setText("Running…")

    def stop_forced_phot(self):
        if self.worker and self.worker.isRunning():
            self.worker.stop()
            self.worker.wait(5000)
        self.run_bar.set_running(False)
        self.progress_label.setText("Stopped")

    # ── Slots ──────────────────────────────────────────────────────────────────

    def _on_progress(self, current: int, total: int, fname: str):
        self.progress_bar.setMaximum(total)
        self.progress_bar.setValue(current)
        self.progress_label.setText(f"[{current}/{total}] {fname}")

    def _on_log(self, msg: str):
        append_timestamped_log(self.log_text, msg)

    def _on_apcorr_update(self, gc_data: dict):
        """Called per frame when growth curve data is available."""
        fname = str(gc_data.get("fname", ""))
        self._gc_accumulator.append(gc_data)
        self._gc_per_frame[fname] = gc_data

        # Append row to apcorr_table
        ri = self.apcorr_table.rowCount()
        self.apcorr_table.insertRow(ri)
        fwhm  = float(gc_data.get("fwhm_px", 0) or 0)
        r_opt = float(gc_data.get("r_opt_px", 0) or 0)
        err_opt = float(gc_data.get("mag_err_opt", float("nan")))
        apcorr_val = float(gc_data.get("apcorr", np.nan))
        r_opt_scale = (r_opt / fwhm) if (fwhm > 0 and r_opt > 0) else float("nan")
        for ci, val in enumerate([
            fname,
            str(gc_data.get("filter", "")),
            f"{apcorr_val:.4f}" if np.isfinite(apcorr_val) else "—",
            str(gc_data.get("n_stars", "")),
            str(gc_data.get("n_apcorr_candidates", "")),
            str(gc_data.get("n_step4_quality_reject", "")),
            str(gc_data.get("n_center_outlier_reject", "")),
            f"{r_opt:.1f}" if r_opt > 0 else "—",
            f"{r_opt_scale:.2f}" if np.isfinite(r_opt_scale) else "—",
            f"{err_opt:.4f}" if np.isfinite(err_opt) else "—",
        ]):
            item = QTableWidgetItem(val)
            item.setTextAlignment(Qt.AlignCenter)
            self.apcorr_table.setItem(ri, ci, item)
        self.apcorr_table.scrollToBottom()

        # If the user hasn't selected a row yet, show the most recent frame
        if not self.apcorr_table.selectionModel().hasSelection():
            self._update_gc_plot(gc_data)

    def _on_center_stats_update(self, row: dict):
        self._append_center_stats_row(row)
        self._refresh_center_summary()
        if not self._stats_seen:
            self._stats_seen = True
            self.tabs.setCurrentIndex(self._stats_tab_index)
            if self.center_stats_table.rowCount() > 0:
                self.center_stats_table.selectRow(0)

    def _on_finished(self, results: dict):
        self.results = results
        rows = results.get("index_rows", [])
        self._populate_results_table(rows)
        center_rows = results.get("center_stats_rows", [])
        if center_rows:
            self._populate_center_stats_table(center_rows)
            self._refresh_center_summary()
        self.run_bar.set_running(False)
        n_ok = sum(1 for r in rows if r.get("status") == "ok")
        self.progress_label.setText(f"Done — {n_ok}/{len(rows)} frames OK")
        self.update_navigation_buttons()

    def _on_error(self, error_type: str, msg: str):
        QMessageBox.critical(self, error_type, msg)
        self.run_bar.set_running(False)
        self.progress_label.setText("Error — see log")

    # ── StepWindowBase overrides ───────────────────────────────────────────────

    def save_state(self):
        self._sync_centering_params()

    def validate_step(self) -> bool:
        idx_path = step7_forced_phot_dir(self.params.P.result_dir) / "photometry_index.csv"
        return idx_path.exists()

    def restore_state(self):
        self._try_load_existing_results()
        self._check_prerequisites()
