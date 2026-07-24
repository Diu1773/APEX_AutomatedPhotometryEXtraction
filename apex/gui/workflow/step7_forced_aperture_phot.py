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

import hashlib
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
    QScrollArea, QFrame,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal

try:
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
    from matplotlib.figure import Figure
    _HAVE_MPL = True
except ImportError:
    _HAVE_MPL = False

from .step_window_base import StepWindowBase
from .run_control import RunControlBar, format_duration as _fmt_duration, progress_status_text
from .log_panel import WorkflowLogWindow, WorkerStatusPanel, append_timestamped_log, show_raised
from .ui_helpers import (
    create_output_reuse_checkbox,
    set_table_row_background,
    status_row_background,
)
from apex.utils.step_paths import (
    step2_cropped_dir,
    step4_dir,
    step5_wcs_dir,
    step6_refbuild_dir,
    step7_forced_phot_dir,
    crop_is_active,
)
from apex.utils.qc_utils import filter_files_by_qc, should_use_frame_quality_qc
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
_FORCED_SIGNATURE_FILE = "forced_phot_signature.json"
# v2: mag_inst redefined as IRAF-style count-rate magnitude
# (INSTRUMENTAL_ZMAG - 2.5*log10(flux_e/exptime)); invalidates v1 caches.
_FORCED_SIGNATURE_VERSION = 2
_FORCED_SIGNATURE_PARAMS = (
    "phot_use_qc_pass_only",
    "fwhm_pix_guess",
    "recenter_aperture",
    "max_recenter_shift",
    "centroid_outlier_px",
    "registration_match_radius_px",
    "registration_min_anchors",
    "forced_r_ap_scale",
    "forced_ref_ap_scale",
    "min_r_ap_px",
    "fitsky_annulus_scale",
    "fitsky_dannulus_scale",
    "annulus_min_gap_px",
    "gain_e_per_adu",
    "rdnoise_e",
    "noise_use_fits_header",
    "noise_reference_binning",
    "noise_scale_by_binning",
    "saturation_adu",
    "datamax_adu",
    "phot_sigma_clip",
    "phot_max_iter",
    "sky_sigma_mode",
    "sky_sigma_includes_rn",
    "sky_sigma_min_n_sky",
    "center_cbox_scale",
    "apcorr_min_snr",
    "apcorr_use_min_n",
    "apcorr_max_sources",
    "ref_cat_max_elong",
    "ref_cat_max_abs_round",
    "ref_cat_sharp_min",
    "ref_cat_sharp_max",
    "ref_cat_min_peak_adu",
)


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

    # ── Main run ──────────────────────────────────────────────────────────────

    def run(self):
        """Run forced photometry on all frames.

        Thin wrapper that delegates the compute to the Qt-free
        ``apex.analysis.forced_photometry.run_forced_photometry``, re-emitting its
        callbacks through the existing Qt signals. Behavior is identical to the
        prior inline implementation.
        """
        from apex.analysis.forced_photometry import run_forced_photometry

        try:
            summary = run_forced_photometry(
                self.file_list,
                self.params,
                self.data_dir,
                self.cache_dir,
                result_dir=self.result_dir,
                output_dir=self.output_dir,
                progress_cb=self.progress.emit,
                log_cb=self.log.emit,
                worker_status_cb=self.worker_status.emit,
                apcorr_cb=self.apcorr_update.emit,
                center_stats_cb=self.center_stats_update.emit,
                error_cb=self.error.emit,
                should_stop=lambda: self._stop_requested,
            )
        except Exception as exc:
            self._log(f"[FORCED] Unhandled error: {exc}\n{traceback.format_exc()}")
            self.error.emit("ForcedPhot", str(exc))
            self.finished.emit({})
            return
        self.finished.emit(summary)


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
        self._run_started_ts: Optional[float] = None
        self._current_forced_signature: dict | None = None
        self._current_forced_files: list[str] = []
        self._reuse_initial_cache_status = True
        self._initial_cache_status: tuple[bool, str] | None = None

        super().__init__(
            step_index=6,
            step_name="Forced Aperture Phot",
            params=params,
            project_state=project_state,
            main_window=main_window,
        )

        self.setup_step_ui()
        self.restore_state()
        self._reuse_initial_cache_status = False
        self._initial_cache_status = None

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
        self.chk_use_existing_output = create_output_reuse_checkbox(
            not bool(getattr(self.params.P, "force_rephot", False)),
            "When enabled, Step 7 loads existing photometry_index.csv and per-frame TSVs if the "
            "selected frame set, inputs, and photometry parameters match the saved signature. "
            "Disable to force forced photometry to rerun.",
        )
        ctrl_layout.addWidget(self.chk_use_existing_output)
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
        self.content_layout.addWidget(self.tabs, stretch=1)

        # ── Tab 0: Centering Stats ─────────────────────────────────────────────
        # The tab content (summary + plot + table) needs ~700px of vertical
        # space.  Cramped windows (small laptops, 720p captures) cannot
        # afford that, so wrap the whole tab in a QScrollArea — both plot
        # and table keep their own minimum height, and the user scrolls
        # when the tab is shorter than the content.
        tab_stats = QWidget()
        tab_stats_outer = QVBoxLayout(tab_stats)
        tab_stats_outer.setContentsMargins(0, 0, 0, 0)

        stats_scroll = QScrollArea()
        stats_scroll.setWidgetResizable(True)
        stats_scroll.setFrameShape(QFrame.NoFrame)
        stats_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        stats_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)

        stats_inner = QWidget()
        tab_stats_layout = QVBoxLayout(stats_inner)
        tab_stats_layout.setContentsMargins(6, 6, 6, 6)

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
            # Hard minimum so X-axis ticks/labels aren't cropped away;
            # maximum so the plot doesn't devour the whole tab on tall
            # windows and push the diagnostics table out of view.
            # 230~280 keeps two histogram axes legible while leaving
            # room for ~5 table rows in a 900px-tall window.
            self._center_canvas.setMinimumHeight(230)
            self._center_canvas.setMaximumHeight(280)
            center_plot_layout.addWidget(self._center_canvas)
            self._init_center_axes()
        else:
            center_plot_layout.addWidget(QLabel("matplotlib not available — install it for center-error plots"))
            self._center_canvas = None
        tab_stats_layout.addWidget(center_plot_group)

        center_table_group = QGroupBox("Per-frame Centering Diagnostics")
        center_table_layout = QVBoxLayout(center_table_group)
        self.center_stats_table = QTableWidget()
        # Enough room for 4–5 rows + header + horizontal scrollbar.
        self.center_stats_table.setMinimumHeight(220)
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

        stats_scroll.setWidget(stats_inner)
        tab_stats_outer.addWidget(stats_scroll)

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
        self.apcorr_table.setColumnCount(11)
        _ap_headers = [
            "File", "Filter", "Apcorr", "N stars",
            "Candidates", "Step4 reject", "Apcorr reject", "Center reject",
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
            "Step4 anchor 품질 컬럼(shape/peak/status)으로 제외된 검출 source 수",
            "Step4 apcorr_candidate=false라서 aperture correction reference에서 제외된 source 수",
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
        out_dir = step7_forced_phot_dir(self.params.P.result_dir)
        idx_path = out_dir / "photometry_index.csv"
        loaded = False
        self._gc_accumulator.clear()
        self._gc_per_frame.clear()
        self._center_stats_rows.clear()
        if hasattr(self, "apcorr_table"):
            self.apcorr_table.setRowCount(0)
        if hasattr(self, "center_stats_table"):
            self.center_stats_table.setRowCount(0)
        if not idx_path.exists():
            return False
        try:
            df = pd.read_csv(idx_path)
            rows = df.to_dict(orient="records")
            self.results = {"index_rows": rows}
            self._populate_results_table(rows)
            self.progress_label.setText(f"Loaded {len(df)} frames from previous run")
            self.update_navigation_buttons()
            loaded = True
        except Exception:
            pass

        apc_path = out_dir / "apcorr_summary.csv"
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
                        str(row.get("n_step4_apcorr_reject", "")),
                        str(row.get("n_center_outlier_reject", "")),
                        "—", "—", "—",
                    ]):
                        item = QTableWidgetItem(val)
                        item.setTextAlignment(Qt.AlignCenter)
                        self.apcorr_table.setItem(ri, ci, item)
            except Exception:
                pass

        stats_path = out_dir / "centering_stats.csv"
        if stats_path.exists():
            try:
                stats_df = pd.read_csv(stats_path)
                rows = stats_df.to_dict(orient="records")
                self._populate_center_stats_table(rows)
                self._refresh_center_summary()
            except Exception:
                pass
        else:
            self._refresh_center_summary()
        return loaded

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
            status = str(row.get("status", "") or "").strip().lower()
            centering = str(row.get("centering_status", "") or "").strip().upper()
            ok = status == "ok" and wcs_ok is not False
            warning = ok and centering in {"CHECK", "LOW_MATCH", "REVIEW"}
            set_table_row_background(self.results_table, ri, status_row_background(ok, warning=warning))

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
        self._sync_cache_params()
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

        qc_info = None
        use_qc = should_use_frame_quality_qc(
            result_dir,
            P,
            "phot_use_qc_pass_only",
            default=False,
        )
        file_list, qc_info = filter_files_by_qc(result_dir, file_list, require_qc=use_qc)
        if use_qc and not file_list:
            QMessageBox.warning(self, "No Frames", "No frames remain after Step 4 QC filtering.")
            return

        self.log_text.clear()
        if use_qc and qc_info is not None:
            if qc_info.get("applied"):
                append_timestamped_log(
                    self.log_text,
                    f"[FORCED][QC] Frame QC filter: {qc_info['kept']}/{qc_info['total']} kept."
                )
            elif qc_info.get("path") is None:
                append_timestamped_log(self.log_text, "[FORCED][QC] frame_quality.csv not found; using all frames.")
            else:
                append_timestamped_log(
                    self.log_text,
                    f"[FORCED][QC] frame_quality.csv ignored ({qc_info['reason']}); using all frames."
                )
        signature = self._build_forced_output_signature(file_list)
        self._current_forced_signature = signature
        self._current_forced_files = list(file_list)
        use_existing = bool(getattr(self, "chk_use_existing_output", None) and self.chk_use_existing_output.isChecked())
        if use_existing:
            cache_complete, cache_reason = self._existing_output_covers(file_list, signature)
            if cache_complete:
                append_timestamped_log(
                    self.log_text,
                    f"[FORCED][CACHE] Existing Step 7 output matches current inputs ({len(file_list)} frames); loading cached result."
                )
                self._try_load_existing_results()
                self.progress_bar.setMaximum(len(file_list))
                self.progress_bar.setValue(len(file_list))
                self.progress_label.setText(f"Cached Step 7 output loaded ({len(file_list)} frames)")
                self.update_navigation_buttons()
                self._current_forced_signature = None
                self._current_forced_files = []
                return
            append_timestamped_log(
                self.log_text,
                f"[FORCED][CACHE] Existing Step 7 output not reusable ({cache_reason}); recomputing all selected frames."
            )
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
        self._run_started_ts = time.monotonic()
        self.progress_bar.setValue(0)
        self.progress_label.setText(
            progress_status_text(0, len(file_list), self._run_started_ts, message="Starting...")
        )

    def stop_forced_phot(self):
        if self.worker and self.worker.isRunning():
            self.worker.stop()
            self.worker.wait(5000)
        self.run_bar.set_running(False)
        elapsed_txt = ""
        if self._run_started_ts is not None:
            elapsed_txt = f" | elapsed {_fmt_duration(time.monotonic() - float(self._run_started_ts))}"
        self.progress_label.setText(f"Stopped{elapsed_txt}")

    # ── Slots ──────────────────────────────────────────────────────────────────

    @staticmethod
    def _signature_value(value):
        if isinstance(value, np.generic):
            value = value.item()
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, float):
            return float(value) if np.isfinite(value) else None
        if isinstance(value, (bool, int, str)) or value is None:
            return value
        if isinstance(value, (list, tuple, set)):
            return [ForcedPhotWindow._signature_value(v) for v in value]
        if isinstance(value, dict):
            return {
                str(k): ForcedPhotWindow._signature_value(v)
                for k, v in sorted(value.items(), key=lambda item: str(item[0]))
            }
        try:
            if pd.isna(value):
                return None
        except Exception:
            pass
        return str(value)

    @staticmethod
    def _file_signature(path: Path | None) -> dict | None:
        if path is None:
            return None
        try:
            p = Path(path)
            if not p.exists():
                return None
            st = p.stat()
            try:
                path_text = str(p.resolve())
            except Exception:
                path_text = str(p)
            return {
                "path": path_text,
                "size": int(st.st_size),
                "mtime_ns": int(st.st_mtime_ns),
            }
        except Exception:
            return None

    @staticmethod
    def _first_existing(candidates: list[Path]) -> Path | None:
        for path in candidates:
            try:
                if path.stat().st_size > 0:
                    return path
            except OSError:
                continue
        return None

    def _resolved_cache_dir(self) -> Path:
        result_dir = Path(self.params.P.result_dir)
        cache_dir = Path(self.params.P.cache_dir)
        return cache_dir if cache_dir.is_absolute() else result_dir / cache_dir

    def _resolve_fits_path_for_signature(self, fname: str) -> Path | None:
        result_dir = Path(self.params.P.result_dir)
        if crop_is_active(result_dir):
            cropped = step2_cropped_dir(result_dir) / fname
            if cropped.exists():
                return cropped
        try:
            original = Path(self.params.get_file_path(fname))
            if original.exists():
                return original
        except Exception:
            pass
        try:
            data_path = Path(self.params.P.data_dir) / fname
            if data_path.exists():
                return data_path
        except Exception:
            pass
        return None

    def _build_forced_output_signature(self, files: list[str]) -> dict:
        result_dir = Path(self.params.P.result_dir)
        cache_dir = self._resolved_cache_dir()
        s4_dir = step4_dir(result_dir)
        s5_dir = step5_wcs_dir(result_dir)
        s6_dir = step6_refbuild_dir(result_dir)

        ref_inputs = [
            s6_dir / "ref_build_signature.json",
            s6_dir / "ref_build_meta.json",
            s6_dir / "ref_frame_stats.csv",
            s6_dir / "master_catalog.tsv",
        ]
        ref_inputs.extend(sorted(s6_dir.glob("ref_catalog*.tsv")))
        ref_signatures = []
        for path in ref_inputs:
            sig = self._file_signature(path)
            if sig is not None:
                ref_signatures.append(sig)

        frame_inputs = []
        for fname in files:
            fits_path = self._resolve_fits_path_for_signature(str(fname))
            wcs_path = self._first_existing(list(astap_wcs_candidates(fits_path)) if fits_path else [])
            detect_json = self._first_existing([
                cache_dir / f"detect_{fname}.json",
                s4_dir / f"detect_{fname}.json",
            ])
            detect_csv = self._first_existing([
                cache_dir / f"detect_{fname}.csv",
                s4_dir / f"detect_{fname}.csv",
            ])
            frame_inputs.append({
                "file": str(fname),
                "fits": self._file_signature(fits_path),
                "wcs_sidecar": self._file_signature(wcs_path),
                "detect_json": self._file_signature(detect_json),
                "detect_csv": self._file_signature(detect_csv),
            })

        payload = {
            "signature_version": _FORCED_SIGNATURE_VERSION,
            "step": "step7_forced_phot",
            "frames": [str(f) for f in files],
            "params": {
                k: self._signature_value(getattr(self.params.P, k, None))
                for k in _FORCED_SIGNATURE_PARAMS
            },
            "inputs": {
                "frame_quality": self._file_signature(s4_dir / "frame_quality.csv"),
                "wcs_summary": self._file_signature(s5_dir / "wcs_solve_summary.csv"),
                "reference_catalogs": ref_signatures,
                "frames": frame_inputs,
            },
        }
        encoded = json.dumps(payload, sort_keys=True, ensure_ascii=False, allow_nan=False)
        payload["signature_hash"] = hashlib.sha1(encoded.encode("utf-8")).hexdigest()
        return payload

    def _stored_forced_signature(self) -> dict | None:
        sig_path = step7_forced_phot_dir(self.params.P.result_dir) / _FORCED_SIGNATURE_FILE
        if not sig_path.exists():
            return None
        try:
            data = json.loads(sig_path.read_text(encoding="utf-8"))
            return data if isinstance(data, dict) else None
        except Exception:
            return None

    def _write_forced_signature(self, signature: dict) -> None:
        out_dir = step7_forced_phot_dir(self.params.P.result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / _FORCED_SIGNATURE_FILE).write_text(
            json.dumps(signature, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False),
            encoding="utf-8",
        )

    def _remove_forced_signature(self) -> None:
        sig_path = step7_forced_phot_dir(self.params.P.result_dir) / _FORCED_SIGNATURE_FILE
        try:
            sig_path.unlink()
        except FileNotFoundError:
            pass
        except Exception:
            pass

    def _forced_signature_matches(self, signature: dict) -> tuple[bool, str]:
        stored = self._stored_forced_signature()
        if not stored:
            return False, "missing signature"
        if stored.get("signature_version") != _FORCED_SIGNATURE_VERSION:
            return False, "signature version mismatch"
        if stored.get("signature_hash") != signature.get("signature_hash"):
            return False, "signature hash mismatch"
        return True, "ok"

    def _sync_cache_params(self) -> None:
        use_existing = bool(
            getattr(self, "chk_use_existing_output", None)
            and self.chk_use_existing_output.isChecked()
        )
        self.params.P.force_rephot = not use_existing
        if hasattr(self, "persist_params"):
            self.persist_params()

    def _existing_output_covers(self, file_list: list[str], signature: dict) -> tuple[bool, str]:
        out_dir = step7_forced_phot_dir(self.params.P.result_dir)
        sig_ok, sig_reason = self._forced_signature_matches(signature)
        if not sig_ok:
            return False, sig_reason
        idx_path = out_dir / "photometry_index.csv"
        if not idx_path.exists():
            return False, "photometry_index.csv missing"
        try:
            df = pd.read_csv(idx_path)
        except Exception as exc:
            return False, f"photometry_index.csv unreadable: {exc}"
        if "file" not in df.columns:
            return False, "photometry_index.csv has no file column"

        status = df.get("status", pd.Series(["ok"] * len(df), index=df.index)).fillna("").astype(str).str.lower()
        ok_files = set(df.loc[status.eq("ok"), "file"].fillna("").astype(str).tolist())
        missing = []
        for fname in file_list:
            fname_s = str(fname)
            if fname_s not in ok_files:
                missing.append(fname_s)
                continue
            if not (out_dir / f"photometry_{fname_s}.tsv").exists():
                missing.append(fname_s)
        if missing:
            preview = ", ".join(missing[:3])
            if len(missing) > 3:
                preview += f", +{len(missing) - 3} more"
            return False, f"missing/incomplete frames {len(missing)}/{len(file_list)}: {preview}"
        return True, "ok"

    def _current_forced_cache_status(self) -> tuple[bool, str]:
        if self._reuse_initial_cache_status and self._initial_cache_status is not None:
            return self._initial_cache_status

        def _finish(valid: bool, reason: str) -> tuple[bool, str]:
            result = (valid, reason)
            if self._reuse_initial_cache_status:
                self._initial_cache_status = result
            return result

        try:
            files = list(self.file_manager.get_file_list())
        except Exception as exc:
            return _finish(False, f"file list unavailable: {exc}")
        if not files:
            return _finish(False, "no current frames")
        use_qc = should_use_frame_quality_qc(
            Path(self.params.P.result_dir),
            self.params.P,
            "phot_use_qc_pass_only",
            default=False,
        )
        files, _ = filter_files_by_qc(
            Path(self.params.P.result_dir),
            files,
            require_qc=use_qc,
        )
        if not files:
            return _finish(False, "no current frames after QC")
        signature = self._build_forced_output_signature(list(files))
        valid, reason = self._existing_output_covers(list(files), signature)
        return _finish(valid, reason)

    def _on_progress(self, current: int, total: int, fname: str):
        self.progress_bar.setMaximum(total)
        self.progress_bar.setValue(current)
        self.progress_label.setText(
            progress_status_text(current, total, self._run_started_ts, message=fname)
        )

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
            str(gc_data.get("n_step4_apcorr_reject", "")),
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
        expected = len(self._current_forced_files)
        if self._current_forced_signature:
            if expected > 0 and n_ok == expected and len(rows) == expected:
                try:
                    self._write_forced_signature(self._current_forced_signature)
                    append_timestamped_log(self.log_text, "[FORCED][CACHE] Output signature saved for future reuse.")
                except Exception as exc:
                    append_timestamped_log(self.log_text, f"[FORCED][CACHE] Signature write failed: {exc}")
            else:
                self._remove_forced_signature()
                append_timestamped_log(self.log_text, "[FORCED][CACHE] Signature not saved: output incomplete.")
        self._current_forced_signature = None
        self._current_forced_files = []
        elapsed_txt = ""
        if self._run_started_ts is not None:
            elapsed_txt = f" | elapsed {_fmt_duration(time.monotonic() - float(self._run_started_ts))}"
        self.progress_label.setText(f"Done — {n_ok}/{len(rows)} frames OK{elapsed_txt}")
        # Second-stage photometric QC: per-frame transparency offsets from the
        # matched bright stars. Catches clouds that image-level Step 4 QC is
        # blind to (fig6 validation). Never blocks the step on failure.
        if n_ok > 0:
            try:
                from apex.analysis.photometric_qc import (
                    run_photometric_qc,
                    summarize_photometric_qc,
                )

                qc_df = run_photometric_qc(self.params.P.result_dir)
                counts = summarize_photometric_qc(qc_df)
                if not qc_df.empty:
                    append_timestamped_log(
                        self.log_text,
                        "[FORCED][QC] Transparency QC: "
                        f"PASS {counts['PASS']} / REVIEW {counts['REVIEW']} / "
                        f"FAIL {counts['FAIL']} / SKIP {counts['SKIP']} -> phot_quality.csv",
                    )
                    flagged = qc_df[qc_df["phot_qc_status"].isin(("REVIEW", "FAIL"))]
                    for _, row in flagged.head(8).iterrows():
                        append_timestamped_log(
                            self.log_text,
                            f"[FORCED][QC] {row['phot_qc_status']}: {row['file']} "
                            f"offset={row['transparency_offset_mag']:+.3f} mag "
                            f"({row['phot_qc_reasons']})",
                        )
                    if len(flagged) > 8:
                        append_timestamped_log(
                            self.log_text,
                            f"[FORCED][QC] ... and {len(flagged) - 8} more flagged frames.",
                        )
            except Exception as exc:
                append_timestamped_log(
                    self.log_text, f"[FORCED][QC] Photometric QC skipped: {exc}"
                )
        self.update_navigation_buttons()

    def _on_error(self, error_type: str, msg: str):
        QMessageBox.critical(self, error_type, msg)
        self.run_bar.set_running(False)
        self._current_forced_signature = None
        self._current_forced_files = []
        elapsed_txt = ""
        if self._run_started_ts is not None:
            elapsed_txt = f" | elapsed {_fmt_duration(time.monotonic() - float(self._run_started_ts))}"
        self.progress_label.setText(f"Error — see log{elapsed_txt}")

    # ── StepWindowBase overrides ───────────────────────────────────────────────

    def save_state(self):
        self._sync_centering_params()
        self._sync_cache_params()

    def validate_step(self) -> bool:
        valid, _ = self._current_forced_cache_status()
        return valid

    def restore_state(self):
        valid, reason = self._current_forced_cache_status()
        if valid:
            self._try_load_existing_results()
            append_timestamped_log(
                self.log_text,
                "[FORCED][CACHE] Restored Step 7 output matching current inputs.",
            )
        elif (step7_forced_phot_dir(self.params.P.result_dir) / "photometry_index.csv").exists():
            append_timestamped_log(
                self.log_text,
                f"[FORCED][CACHE] Previous output not restored ({reason}).",
            )
        self._check_prerequisites()
