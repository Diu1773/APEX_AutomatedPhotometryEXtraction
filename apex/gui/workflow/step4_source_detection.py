"""
Step 4: Source Detection Window
Parallel source detection with segmentation and peak finding
APEX source detection workflow.
"""

from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QGroupBox, QMessageBox, QTextEdit, QComboBox, QDialog,
    QFormLayout, QLineEdit, QDialogButtonBox, QSplitter, QApplication,
    QProgressBar, QCheckBox, QSpinBox, QDoubleSpinBox, QTableWidget,
    QTableWidgetItem, QHeaderView, QAbstractItemView, QSlider, QGridLayout,
    QWidget, QTabWidget, QScrollArea,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt5.QtGui import QColor
from apex.gui.widgets.fits_viewer import FITSViewerWidget, OverlayMarker
from pathlib import Path
from typing import Optional
import copy
import shutil
import json
from collections import OrderedDict
import numpy as np
import pandas as pd
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from astropy.io import fits
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.time import Time
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from scipy.ndimage import gaussian_filter, median_filter
from scipy.spatial import cKDTree as KDTree

from .step_window_base import StepWindowBase
from .run_control import RunControlBar, format_duration, progress_status_text
from .ui_helpers import (
    create_cache_action_button,
    create_collapsible_section,
    create_detection_cache_checkbox,
    create_parameter_button,
    configure_parameter_dialog,
)
from apex.core.cache_manager import StepCacheManager
from .log_panel import WorkflowLogWindow, WorkerStatusPanel, append_timestamped_log, show_raised
from apex.utils.constants import get_parallel_workers
from apex.utils.fast_stats import finite_nanmedian, finite_nanstd, robust_median_mad
from apex.utils.step_paths import (
    step2_cropped_dir,
    crop_is_active,
    step4_dir,
)
from apex.utils.cache_utils import (
    DETECTION_CACHE_SCHEMA_VERSION,
    build_detection_cache_signature,
    detection_cache_signature_matches,
    norm_path_key,
    normalize_detect_engine as _normalize_detect_engine_util,
)
from apex.utils.astro_utils import compute_airmass_from_header
from apex.utils.qc_utils import is_passed_value
from apex.utils.source_quality import compute_source_quality

_DETECT_MODE_PRESETS = {
    "normal": {
        "detect_sigma": 3.2,
        "minarea_pix": 3,
        "deblend_enable": True,
        "deblend_nthresh": 64,
        "deblend_cont": 0.004,
        "deblend_max_labels": 4000,
        "deblend_label_hard_max": 7000,
        "dao_refine_enable": True,
        "peak_pass_enable": False,
        "peak_nsigma": 3.4,
        "peak_max_add": 500,
        "fwhm_qc_max_sources": 20,
        "bkg2d_downsample": 8,
        "deblend_mode": "exponential",   # faster than linear (photutils>=1.0.1)
    },
    "crowded": {
        "detect_sigma": 3.0,
        "minarea_pix": 2,
        "deblend_enable": True,
        "deblend_nthresh": 96,
        "deblend_cont": 0.0015,
        "deblend_max_labels": 6000,
        "deblend_label_hard_max": 10000,
        "dao_refine_enable": True,
        "peak_pass_enable": False,
        "peak_nsigma": 3.8,
        "peak_max_add": 300,
        "fwhm_qc_max_sources": 20,
        "bkg2d_downsample": 8,
        "deblend_mode": "exponential",
    },
    "faint": {
        "detect_sigma": 2.6,
        "minarea_pix": 2,
        "deblend_enable": True,
        "deblend_nthresh": 48,
        "deblend_cont": 0.0060,
        "deblend_max_labels": 5000,
        "deblend_label_hard_max": 9000,
        "dao_refine_enable": False,
        "peak_pass_enable": False,
        "peak_nsigma": 2.8,
        "peak_max_add": 1200,
        "fwhm_qc_max_sources": 20,
        "bkg2d_downsample": 8,
        "deblend_mode": "exponential",
    },
}


def _normalize_detect_mode(value) -> str:
    s = str(value or "normal").strip().lower()
    aliases = {
        "default": "normal",
        "standard": "normal",
        "dense": "crowded",
        "cluster": "crowded",
        "dim": "faint",
        "deep": "faint",
        "custom": "custom",
    }
    s = aliases.get(s, s)
    if s in _DETECT_MODE_PRESETS:
        return s
    if s == "custom":
        return s
    return "normal"


def _get_detect_mode_preset(mode: str) -> dict:
    mode_key = _normalize_detect_mode(mode)
    return dict(_DETECT_MODE_PRESETS.get(mode_key, _DETECT_MODE_PRESETS["normal"]))


def _get_detect_mode_from_params(params_obj) -> str:
    return _normalize_detect_mode(getattr(params_obj, "detect_mode", "normal"))


def _normalize_detect_engine(value) -> str:
    return _normalize_detect_engine_util(value)


class DetectionWorker(QThread):
    """Worker thread for parallel source detection"""
    progress = pyqtSignal(int, int, str, int)  # current, total, message, active_workers
    file_done = pyqtSignal(str, dict)  # filename, result
    finished = pyqtSignal(dict)  # summary
    error = pyqtSignal(str, str)  # filename, error message
    worker_status = pyqtSignal(int, str, str, int)  # worker_id, filename, status, progress(0-100)

    def __init__(self, file_list, params, data_dir, cache_dir, use_cropped=False, filter_sigma_map=None):
        super().__init__()
        self.file_list = file_list
        self.params = params
        self.data_dir = Path(data_dir)
        self.result_dir = Path(getattr(params.P, "result_dir", data_dir))
        self.cache_dir = Path(cache_dir)
        self.use_cropped = use_cropped
        self.filter_sigma_map = filter_sigma_map or {}
        self._stop_requested = False
        self._worker_file_map = {}  # thread_id -> (worker_num, filename)
        self._executor = None

    def stop(self):
        self._stop_requested = True

    def run(self):
        """Run detection on all files"""
        # Suppress tqdm and other console output from photutils
        import warnings
        import os
        os.environ['TQDM_DISABLE'] = '1'
        warnings.filterwarnings('ignore', category=UserWarning)

        try:
            from photutils.segmentation import SourceCatalog, detect_sources, deblend_sources
            from photutils.detection import find_peaks
            from photutils.detection import DAOStarFinder

            print(f"[DetectionWorker] Starting with {len(self.file_list)} files")
            print(f"[DetectionWorker] Data dir: {self.data_dir}, use_cropped: {self.use_cropped}")

            results = {}
            total = len(self.file_list)
            P = self.params.P
            step4_out = step4_dir(self.result_dir)
            step4_out.mkdir(parents=True, exist_ok=True)

            # Get parameters
            detect_sigma_base = float(getattr(P, 'detect_sigma', 3.2))
            minarea_pix = int(getattr(P, 'minarea_pix', 3))
            deblend_enable = getattr(P, 'deblend_enable', True)
            deblend_nthresh = int(getattr(P, 'deblend_nthresh', 64))
            deblend_cont = float(getattr(P, 'deblend_cont', 0.004))
            deblend_max_labels = int(getattr(P, 'deblend_max_labels', 4000))
            deblend_label_hard_max = int(getattr(P, 'deblend_label_hard_max', 7000))
            deblend_mode = str(getattr(P, 'deblend_mode', 'exponential')).strip().lower()
            bkg2d_enable = getattr(P, 'bkg2d_in_detect', True)
            bkg2d_box = int(getattr(P, 'bkg2d_box', 64))
            dao_refine_enable = getattr(P, 'dao_refine_enable', False)
            dao_fwhm_px = float(getattr(P, 'dao_fwhm_px', getattr(P, 'fwhm_seed_px', 6.0)))
            dao_sharp_lo = float(getattr(P, 'dao_sharp_lo', 0.2))
            dao_sharp_hi = float(getattr(P, 'dao_sharp_hi', 1.0))
            dao_round_lo = float(getattr(P, 'dao_round_lo', -0.5))
            dao_round_hi = float(getattr(P, 'dao_round_hi', 0.5))
            dao_match_tol = float(getattr(P, 'dao_match_tol_px', 2.0))
            peak_enable = bool(getattr(P, 'peak_pass_enable', False))
            peak_nsigma = float(getattr(P, 'peak_nsigma', 3.2))
            peak_scales = getattr(P, 'peak_kernel_scales', "0.9,1.3")
            peak_min_sep = float(getattr(P, 'peak_min_sep_px', 4.0))
            peak_max_add = int(getattr(P, 'peak_max_add', 600))
            peak_max_elong = float(getattr(P, 'peak_max_elong', 1.6))
            peak_sharp_lo = float(getattr(P, 'peak_sharp_lo', 0.12))
            peak_skip_if_nsrc_ge = int(getattr(P, 'peak_skip_if_nsrc_ge', 4500))
            if isinstance(peak_scales, str):
                peak_scales = [float(s) for s in peak_scales.split(",") if s.strip()]

            max_workers = get_parallel_workers(self.params)

            print(f"[DetectionWorker] max_workers={max_workers}, sigma_base={detect_sigma_base}")

            # Track worker assignments
            import threading
            worker_counter = [0]
            worker_counter_lock = threading.Lock()
            thread_to_worker = {}

            def get_worker_id():
                """Get or assign worker ID for current thread"""
                tid = threading.get_ident()
                if tid not in thread_to_worker:
                    with worker_counter_lock:
                        thread_to_worker[tid] = worker_counter[0]
                        worker_counter[0] += 1
                return thread_to_worker[tid]

            def detect_single(filename):
                """Detect sources in a single file"""
                if self._stop_requested:
                    return filename, None

                worker_id = get_worker_id()
                short_name = filename[:25] + "..." if len(filename) > 28 else filename

                try:
                    # Stage 1: Loading
                    self.worker_status.emit(worker_id, short_name, "Loading", 10)

                    # Determine path
                    if self.use_cropped:
                        cropped_dir = step2_cropped_dir(self.result_dir)
                        file_path = cropped_dir / filename
                    else:
                        file_path = self.params.get_file_path(filename)
                    source_sig = build_detection_cache_signature(
                        Path(file_path),
                        use_cropped=bool(self.use_cropped),
                        detect_engine=getattr(P, 'detect_engine', 'segm'),
                    )

                    # Load FITS
                    with fits.open(file_path) as hdul:
                        data = hdul[0].data.astype(float)
                        header = hdul[0].header

                    # Get filter - preserve original case from header
                    filt = header.get('FILTER', '').strip()
                    filt_lower = filt.lower()

                    # Get sigma for this filter from flexible mapping
                    nsig = detect_sigma_base
                    if filt in self.filter_sigma_map:
                        nsig = float(self.filter_sigma_map[filt])
                    elif filt_lower in self.filter_sigma_map:
                        nsig = float(self.filter_sigma_map[filt_lower])
                    elif filt.upper() in self.filter_sigma_map:
                        nsig = float(self.filter_sigma_map[filt.upper()])

                    # Stage 2: Background
                    self.worker_status.emit(worker_id, short_name, "Background", 25)
                    if self._stop_requested:
                        return filename, None

                    def refine_centroid(img, x, y, seed_fwhm_px):
                        h, w = img.shape
                        r = max(int(round(3.5 * max(seed_fwhm_px, 2.0))), 8)
                        xi, yi = int(round(x)), int(round(y))
                        x0, x1 = max(0, xi - r), min(w, xi + r + 1)
                        y0, y1 = max(0, yi - r), min(h, yi + r + 1)
                        if (x1 - x0) < 9 or (y1 - y0) < 9:
                            return None
                        cut = img[y0:y1, x0:x1]
                        _, med, _ = sigma_clipped_stats(cut, sigma=3.0, maxiters=5, mask=~np.isfinite(cut))
                        z = cut - med
                        z[~np.isfinite(z)] = 0.0
                        z[z < 0] = 0.0
                        s = np.nansum(z)
                        if s <= 0:
                            return None
                        yy, xx = np.mgrid[y0:y1, x0:x1]
                        xc = float(np.nansum(xx * z) / s)
                        yc = float(np.nansum(yy * z) / s)
                        return xc, yc, float(med), float(0.0)

                    def shape_metrics(img, xc, yc, seed_fwhm_px):
                        h, w = img.shape
                        r = max(int(round(2.5 * max(seed_fwhm_px, 2.0))), 7)
                        xi, yi = int(round(xc)), int(round(yc))
                        x0, x1 = max(0, xi - r), min(w, xi + r + 1)
                        y0, y1 = max(0, yi - r), min(h, yi + r + 1)
                        if (x1 - x0) < 7 or (y1 - y0) < 7:
                            return np.inf, 0.0
                        cut = img[y0:y1, x0:x1].astype(float)
                        _, med, _ = sigma_clipped_stats(cut, sigma=3.0, maxiters=5, mask=~np.isfinite(cut))
                        z = cut - med
                        z[~np.isfinite(z)] = 0.0
                        z[z < 0] = 0.0
                        s = float(np.nansum(z))
                        if s <= 0:
                            return np.inf, 0.0
                        yy, xx = np.mgrid[y0:y1, x0:x1]
                        xbar = float(np.nansum(xx * z) / s)
                        ybar = float(np.nansum(yy * z) / s)
                        dx = xx - xbar
                        dy = yy - ybar
                        ixx = float(np.nansum(z * dx * dx) / s)
                        iyy = float(np.nansum(z * dy * dy) / s)
                        ixy = float(np.nansum(z * dx * dy) / s)
                        wv, _ = np.linalg.eigh(np.array([[ixx, ixy], [ixy, iyy]]))
                        wv = np.sort(np.maximum(wv, 1e-12))
                        elong = float(np.sqrt(wv[1] / wv[0]))
                        cx, cy = int(round(xbar)), int(round(ybar))

                        def box_sum(R):
                            xa0, xa1 = max(0, cx - R), min(w, cx + R + 1)
                            ya0, ya1 = max(0, cy - R), min(h, cy + R + 1)
                            b = img[ya0:ya1, xa0:xa1] - med
                            b = np.where(np.isfinite(b), b, 0.0)
                            return float(np.nansum(b[b > 0]))

                        s33 = box_sum(1)
                        s55 = box_sum(2)
                        sharp = (s33 / max(s55, 1e-9)) if s55 > 0 else 0.0
                        return elong, sharp

                    def radial_fwhm(img, xc, yc, seed_fwhm_px, ann_in_scale, ann_out_scale,
                                    ann_min_gap_px, ann_min_width_px, sigma=3.0, maxiters=5, dr=0.5):
                        r_in = max(ann_in_scale * seed_fwhm_px, ann_min_gap_px, 12.0)
                        r_out = max(r_in + ann_out_scale * seed_fwhm_px, r_in + ann_min_width_px, 20.0)
                        h, w = img.shape
                        rmax = int(max(6.0 * seed_fwhm_px, r_out + 6.0))
                        xi, yi = int(round(xc)), int(round(yc))
                        x0, x1 = max(0, xi - rmax), min(w, xi + rmax + 1)
                        y0, y1 = max(0, yi - rmax), min(h, yi + rmax + 1)
                        if (x1 - x0) < 5 or (y1 - y0) < 5:
                            return np.nan, 0.0, 0.0
                        cut = img[y0:y1, x0:x1].astype(float, copy=False)
                        yy, xx = np.mgrid[y0:y1, x0:x1]
                        rr = np.hypot(xx - xc, yy - yc)

                        ann_mask = (rr <= r_out) & (rr > r_in)
                        vals = cut[ann_mask]
                        vals = vals[np.isfinite(vals)]
                        if vals.size:
                            sc = SigmaClip(sigma=sigma, maxiters=maxiters)
                            vv = sc(vals)
                            vv = vv.compressed() if np.ma.isMaskedArray(vv) else np.asarray(vv)
                            vv = vv[np.isfinite(vv)]
                            sky_med = finite_nanmedian(vv, default=finite_nanmedian(vals, default=0.0))
                            sky_std = finite_nanstd(vv, ddof=1, default=0.0) if vv.size > 1 else finite_nanstd(vals, ddof=1, default=0.0) if vals.size > 1 else 0.0
                        else:
                            sky_med = 0.0
                            sky_std = 0.0

                        val = cut - sky_med
                        edges = np.arange(0.0, rmax + dr, dr)
                        centers = 0.5 * (edges[:-1] + edges[1:])
                        prof = np.full_like(centers, np.nan, dtype=float)
                        _n_bins = len(centers)
                        _bi = np.digitize(rr.ravel(), edges) - 1
                        _vf = val.ravel()
                        _mask = np.isfinite(_vf) & (_bi >= 0) & (_bi < _n_bins)
                        if _mask.any():
                            _cnt = np.bincount(_bi[_mask], minlength=_n_bins)
                            _sm = np.bincount(_bi[_mask], weights=_vf[_mask], minlength=_n_bins)
                            _ok = _cnt > 0
                            prof[_ok] = _sm[_ok] / _cnt[_ok]
                        if not np.isfinite(prof).any():
                            return np.nan, sky_med, sky_std
                        k_peak = int(max(2, np.round(1.5 * seed_fwhm_px / max(dr, 1e-9))))
                        use_slice = prof[:min(k_peak, len(prof))]
                        peak = np.nanmax(use_slice) if np.isfinite(use_slice).any() else np.nanmax(prof)
                        if not (np.isfinite(peak) and peak > 0):
                            return np.nan, sky_med, sky_std
                        half = 0.5 * peak
                        idx = np.where((prof[:-1] >= half) & (prof[1:] < half))[0]
                        if len(idx) == 0:
                            return np.nan, sky_med, sky_std
                        i = int(idx[0])
                        x1_, y1_ = centers[i], prof[i]
                        x2_, y2_ = centers[i + 1], prof[i + 1]
                        if not (np.isfinite(y1_) and np.isfinite(y2_) and (y1_ != y2_)):
                            return np.nan, sky_med, sky_std
                        r_half = x1_ + (half - y1_) * (x2_ - x1_) / (y2_ - y1_)
                        return 2.0 * float(r_half), sky_med, sky_std

                    def copy_if_different(src: Path, dst: Path) -> None:
                        try:
                            if Path(src).resolve() == Path(dst).resolve():
                                return
                        except Exception:
                            if str(src) == str(dst):
                                return
                        shutil.copy2(src, dst)

                    def peak_pass_candidates(data_det, med, seed_fwhm_px):
                        h, w = data_det.shape
                        fill = finite_nanmedian(data_det, default=0.0)
                        det_safe = np.where(np.isfinite(data_det), data_det, fill)
                        out = []

                        for s in (peak_scales or [0.9, 1.3]):
                            if self._stop_requested:
                                return []
                            sig = max(0.7, (s * seed_fwhm_px) / 2.355)
                            fim_local = gaussian_filter(det_safe - med, sig, mode="nearest")
                            _, mF, stdF = sigma_clipped_stats(fim_local, sigma=3.0, maxiters=5)
                            thr = mF + peak_nsigma * max(stdF, 1e-9)
                            tbl = find_peaks(fim_local, threshold=thr, box_size=9)
                            if tbl is None or len(tbl) == 0:
                                continue
                            for r in tbl:
                                x = float(r["x_peak"])
                                y = float(r["y_peak"])
                                rc = refine_centroid(data_det, x, y, seed_fwhm_px)
                                if rc is None:
                                    continue
                                xc, yc, _, _ = rc
                                if not (2 <= xc < w - 2 and 2 <= yc < h - 2):
                                    continue
                                elong, sharp = shape_metrics(data_det, xc, yc, seed_fwhm_px)
                                if (elong <= peak_max_elong) and (sharp >= peak_sharp_lo):
                                    out.append((xc, yc))

                        if not out:
                            return []
                        pts = np.array(out, float)
                        keep = np.ones(len(pts), bool)
                        tree = KDTree(pts)
                        for i in range(len(pts)):
                            if self._stop_requested:
                                return []
                            if not keep[i]:
                                continue
                            j = tree.query_ball_point(pts[i], r=peak_min_sep)
                            for k in j:
                                if k <= i:
                                    continue
                                keep[k] = False
                        pts = pts[keep]
                        if len(pts) > peak_max_add:
                            vals = data_det[pts[:, 1].astype(int), pts[:, 0].astype(int)]
                            order = np.argsort(vals)[::-1][:peak_max_add]
                            pts = pts[order]
                        return [tuple(xy) for xy in pts]

                    n_sources = 0
                    positions = []
                    fwhm_values = []
                    fwhm_seed = float(getattr(P, 'fwhm_seed_px', getattr(P, 'fwhm_pix_guess', 6.0) or 6.0))
                    fwhm_prelim = fwhm_seed
                    detect_engine = _normalize_detect_engine(getattr(P, 'detect_engine', 'segm'))
                    detect_method = "segm" if detect_engine == "segm" else detect_engine
                    median_elongation = np.nan
                    median_roundness = np.nan
                    sat_star_count = 0
                    elong_map = {}
                    fwhm_by_xy = {}
                    source_dao_info = {}
                    segm = None

                    # Background estimation and primary detection.
                    if self._stop_requested:
                        return filename, None

                    if detect_engine == "sep":
                        try:
                            import sep
                        except Exception as exc:
                            raise RuntimeError(
                                "SEP engine selected but the 'sep' package is not installed"
                            ) from exc

                        self.worker_status.emit(worker_id, short_name, "SEP Background", 25)
                        data_sep = np.ascontiguousarray(data.astype(np.float32, copy=False))
                        sep_mask = ~np.isfinite(data_sep)
                        mask_arg = sep_mask if bool(sep_mask.any()) else None
                        data_sep = np.where(np.isfinite(data_sep), data_sep, 0.0).astype(np.float32, copy=False)

                        if bkg2d_enable:
                            fw = max(1, int(getattr(P, 'bkg2d_filter_size', 3)))
                            sep_bkg = sep.Background(
                                data_sep,
                                mask=mask_arg,
                                bw=max(4, int(bkg2d_box)),
                                bh=max(4, int(bkg2d_box)),
                                fw=fw,
                                fh=fw,
                            )
                            bkg_median = float(sep_bkg.globalback)
                            bkg_rms = float(sep_bkg.globalrms)
                            data_sub = np.asarray(data_sep - sep_bkg.back(), dtype=float)
                        else:
                            data_filled = np.where(np.isfinite(data), data, 0.0)
                            _ds = max(1, int(round(data_filled.shape[0] ** 0.5)) // 8)
                            _, bkg_median, bkg_rms = sigma_clipped_stats(
                                data_filled[::_ds, ::_ds], sigma=3.0, maxiters=3
                            )
                            bkg_median = float(bkg_median)
                            bkg_rms = float(bkg_rms)
                            data_sub = data - bkg_median

                        work = np.where(np.isfinite(data_sub), data_sub, 0.0)
                        threshold = nsig * max(bkg_rms, 1e-9)

                        self.worker_status.emit(worker_id, short_name, "SEP Detecting", 40)
                        sep_nlevels = int(deblend_nthresh) if deblend_enable else 1
                        sep_cont = float(deblend_cont) if deblend_enable else 1.0
                        objects = sep.extract(
                            np.ascontiguousarray(work.astype(np.float32, copy=False)),
                            nsig,
                            err=max(bkg_rms, 1e-9),
                            mask=mask_arg,
                            minarea=max(1, int(minarea_pix)),
                            deblend_nthresh=max(1, sep_nlevels),
                            deblend_cont=sep_cont,
                            clean=True,
                            segmentation_map=False,
                        )

                        if objects is not None and len(objects) > 0:
                            x_arr = np.asarray(objects["x"], float)
                            y_arr = np.asarray(objects["y"], float)
                            flux_arr = np.asarray(objects["flux"], float)
                            peak_arr = np.asarray(objects["peak"], float)
                            a_arr = np.asarray(objects["a"], float)
                            b_arr = np.asarray(objects["b"], float)
                            flag_arr = np.asarray(objects["flag"], int)
                            elong_arr = np.divide(
                                a_arr, b_arr,
                                out=np.full_like(a_arr, np.inf, dtype=float),
                                where=b_arr > 0,
                            )
                            round_arr = np.divide(
                                a_arr - b_arr, a_arr + b_arr,
                                out=np.full_like(a_arr, np.nan, dtype=float),
                                where=(a_arr + b_arr) > 0,
                            )
                            median_elongation = finite_nanmedian(elong_arr, default=np.nan)
                            median_roundness = finite_nanmedian(np.abs(round_arr), default=np.nan)

                            order = np.argsort(flux_arr)[::-1]
                            x_arr, y_arr = x_arr[order], y_arr[order]
                            flux_arr, peak_arr = flux_arr[order], peak_arr[order]
                            elong_arr, round_arr = elong_arr[order], round_arr[order]
                            flag_arr = flag_arr[order]

                            elong_max = max(
                                float(getattr(P, 'fwhm_elong_max', 1.3)),
                                float(getattr(P, 'peak_max_elong', 1.6)),
                            )
                            keep = np.isfinite(x_arr) & np.isfinite(y_arr) & (elong_arr <= elong_max)
                            x_arr, y_arr = x_arr[keep], y_arr[keep]
                            flux_arr, peak_arr = flux_arr[keep], peak_arr[keep]
                            elong_arr, round_arr = elong_arr[keep], round_arr[keep]
                            flag_arr = flag_arr[keep]

                            detect_keep_max = int(getattr(P, 'detect_keep_max', 6000))
                            if len(x_arr) > detect_keep_max:
                                sl = slice(0, detect_keep_max)
                                x_arr, y_arr = x_arr[sl], y_arr[sl]
                                flux_arr, peak_arr = flux_arr[sl], peak_arr[sl]
                                elong_arr, round_arr = elong_arr[sl], round_arr[sl]
                                flag_arr = flag_arr[sl]

                            cand = np.vstack([x_arr, y_arr]).T if len(x_arr) else np.zeros((0, 2), float)
                            iso_min_sep = float(getattr(P, 'iso_min_sep_pix', 18.0))
                            if len(cand) > 1:
                                tree = KDTree(cand)
                                d, _ = tree.query(cand, k=2)
                                keep_iso = d[:, 1] >= iso_min_sep
                                cand = cand[keep_iso]
                                flux_arr, peak_arr = flux_arr[keep_iso], peak_arr[keep_iso]
                                elong_arr, round_arr = elong_arr[keep_iso], round_arr[keep_iso]
                                flag_arr = flag_arr[keep_iso]

                            if len(cand) > 0:
                                H, W = data.shape
                                ix = np.clip(np.round(cand[:, 0]).astype(int), 0, W - 1)
                                iy = np.clip(np.round(cand[:, 1]).astype(int), 0, H - 1)
                                sat_adu = float(getattr(P, 'saturation_adu', 60000.0))
                                ok = data[iy, ix] < sat_adu
                                sat_star_count = int((~ok).sum())
                                cand = cand[ok]
                                flux_arr, peak_arr = flux_arr[ok], peak_arr[ok]
                                elong_arr, round_arr = elong_arr[ok], round_arr[ok]
                                flag_arr = flag_arr[ok]

                            positions = [tuple(map(float, p)) for p in cand]
                            n_sources = len(positions)
                            for (x, y), e, r, peak, flux, flag in zip(
                                positions, elong_arr, round_arr, peak_arr, flux_arr, flag_arr
                            ):
                                xf, yf = float(x), float(y)
                                elong_val = float(e) if np.isfinite(e) else np.nan
                                elong_map[(xf, yf)] = elong_val
                                sharp_proxy = 1.0 / max(elong_val, 1e-6) if np.isfinite(elong_val) else np.nan
                                source_dao_info[(xf, yf)] = {
                                    'sharpness': float(np.clip(sharp_proxy, 0.0, 1.0)) if np.isfinite(sharp_proxy) else np.nan,
                                    'roundness1': float(r) if np.isfinite(r) else np.nan,
                                    'roundness2': float(r) if np.isfinite(r) else np.nan,
                                    'peak': float(peak) if np.isfinite(peak) else np.nan,
                                    'flux': float(flux) if np.isfinite(flux) else np.nan,
                                    'sep_flag': int(flag),
                                }
                    else:
                        # Background estimation (Jupyter-style downsample median)
                        data_filled = np.where(np.isfinite(data), data, 0.0)
                        if bkg2d_enable:
                            if self._stop_requested:
                                return filename, None
                            ds = max(1, int(getattr(P, 'bkg2d_downsample', 4)))
                            k = max(3, int(round(bkg2d_box / ds)))
                            small = data_filled[::ds, ::ds]
                            bkg_small = median_filter(small, size=k, mode="nearest")
                            if self._stop_requested:
                                return filename, None
                            bkg = np.repeat(np.repeat(bkg_small, ds, axis=0), ds, axis=1)[:data.shape[0], :data.shape[1]]
                            if self._stop_requested:
                                return filename, None
                            data_sub = data - bkg
                        else:
                            data_sub = data.copy()

                        work = np.where(np.isfinite(data_sub), data_sub, 0.0)
                        _ds = max(1, int(round(work.shape[0] ** 0.5)) // 8)
                        _, bkg_median, _ = sigma_clipped_stats(work[::_ds, ::_ds], sigma=3.0, maxiters=3)
                        if self._stop_requested:
                            return filename, None

                        # Smoothing (Jupyter-style)
                        sig = max(0.8, fwhm_seed / 2.355)
                        fim = gaussian_filter(work - bkg_median, sig, mode="nearest")
                        if self._stop_requested:
                            return filename, None

                        # Threshold (median + nsig * std)
                        _, mF, stdF = sigma_clipped_stats(fim[::_ds, ::_ds], sigma=3.0, maxiters=3)
                        bkg_rms = float(stdF)
                        threshold = mF + nsig * max(stdF, 1e-9)

                        # Stage 3: Detection
                        self.worker_status.emit(worker_id, short_name, "Detecting", 40)

                        # Segmentation detection
                        if detect_engine == "segm":
                            segm = detect_sources(
                                fim, threshold=threshold,
                                npixels=max(3, minarea_pix), connectivity=8
                            )
                        if self._stop_requested:
                            return filename, None

                    if segm is not None and segm.nlabels > 0:
                        # Stage 4: Deblending
                        self.worker_status.emit(worker_id, short_name, f"Deblend ({segm.nlabels})", 55)

                        # Deblending with soft/hard limits
                        if deblend_enable:
                            nlabels0 = int(segm.nlabels)
                            nlevels_soft = int(getattr(P, 'deblend_nlevels_soft', 32))
                            cont_soft = float(getattr(P, 'deblend_contrast_soft', 0.005))
                            if nlabels0 < deblend_label_hard_max:
                                nlevels = deblend_nthresh
                                cont = deblend_cont
                                if nlabels0 >= deblend_max_labels:
                                    nlevels = min(nlevels, nlevels_soft)
                                    cont = max(cont, cont_soft)
                                if self._stop_requested:
                                    return filename, None
                                try:
                                    # mode='exponential' (photutils>=1.0.1) is ~2x faster
                                    # than 'linear'; falls back gracefully on older versions
                                    segm = deblend_sources(
                                        fim, segm,
                                        npixels=minarea_pix,
                                        nlevels=nlevels,
                                        contrast=cont,
                                        mode=deblend_mode,
                                        progress_bar=False
                                    )
                                except TypeError:
                                    # older photutils: no mode/progress_bar kwarg
                                    try:
                                        segm = deblend_sources(
                                            fim, segm,
                                            npixels=minarea_pix,
                                            nlevels=nlevels,
                                            contrast=cont
                                        )
                                    except Exception:
                                        pass  # Keep original segmentation
                                except Exception:
                                    pass  # Keep original segmentation
                            else:
                                self.worker_status.emit(worker_id, short_name, f"Deblend skip ({nlabels0})", 58)

                        # Stage 5: Catalog
                        self.worker_status.emit(worker_id, short_name, "Catalog", 75)

                        cat = SourceCatalog(work, segm)
                        tab = cat.to_table(columns=("xcentroid", "ycentroid", "elongation", "segment_flux"))
                        xcen = np.asarray(tab["xcentroid"], float)
                        ycen = np.asarray(tab["ycentroid"], float)
                        flux = np.asarray(tab["segment_flux"], float)
                        elong = np.asarray(tab["elongation"], float) if "elongation" in tab.colnames else np.ones_like(xcen)
                        if len(elong):
                            try:
                                median_elongation = float(np.nanmedian(elong))
                            except Exception:
                                median_elongation = np.nan

                        order = np.argsort(flux)[::-1]
                        x0 = xcen[order]
                        y0 = ycen[order]
                        e0 = elong[order]

                        elong_max = float(getattr(P, 'fwhm_elong_max', 1.3))
                        keep = np.isfinite(x0) & np.isfinite(y0) & (e0 <= elong_max)
                        x0, y0, e0 = x0[keep], y0[keep], e0[keep]

                        # Estimate frame FWHM from bright sources before DAO refinement.
                        try:
                            _n_prelim = min(10, len(x0))
                            _fwhm_prelim_vals = []
                            for _xi, _yi in zip(x0[:_n_prelim], y0[:_n_prelim]):
                                _rc = refine_centroid(data, _xi, _yi, fwhm_seed)
                                if _rc is None:
                                    continue
                                _xc, _yc, _, _ = _rc
                                _fpx, _, _ = radial_fwhm(
                                    data, _xc, _yc, fwhm_seed,
                                    float(getattr(P, 'fitsky_annulus_scale', 4.0)),
                                    float(getattr(P, 'fitsky_dannulus_scale', 2.0)),
                                    float(getattr(P, 'annulus_min_gap_px', 6.0)),
                                    float(getattr(P, 'annulus_min_width_px', 12.0)),
                                    sigma=float(getattr(P, 'annulus_sigma_clip', 3.0)),
                                    maxiters=int(getattr(P, 'fitsky_max_iter', 5)),
                                    dr=float(getattr(P, 'fwhm_dr', 0.5)),
                                )
                                if np.isfinite(_fpx) and _fpx > 0:
                                    _fpx = float(_fpx)
                                    _fwhm_prelim_vals.append(_fpx)
                                    fwhm_by_xy[(float(_xi), float(_yi))] = _fpx
                            if len(_fwhm_prelim_vals) >= 3:
                                fwhm_prelim = finite_nanmedian(_fwhm_prelim_vals, default=0.0)
                                self.worker_status.emit(
                                    worker_id, short_name,
                                    f"FWHM={fwhm_prelim:.1f}px", 77
                                )
                        except Exception:
                            pass

                        cand = np.vstack([x0, y0]).T
                        e_cand = e0.copy()
                        detect_keep_max = int(getattr(P, 'detect_keep_max', 6000))
                        if len(cand) > detect_keep_max:
                            cand = cand[:detect_keep_max]
                            e_cand = e_cand[:detect_keep_max]

                        iso_min_sep = float(getattr(P, 'iso_min_sep_pix', 18.0))
                        if len(cand) > 1:
                            tree = KDTree(cand)
                            d, _ = tree.query(cand, k=2)
                            sep = d[:, 1]
                            keep_iso = sep >= iso_min_sep
                            cand = cand[keep_iso]
                            e_cand = e_cand[keep_iso]

                        if len(cand) > 0:
                            H, W = data.shape
                            ix = np.clip(np.round(cand[:, 0]).astype(int), 0, W - 1)
                            iy = np.clip(np.round(cand[:, 1]).astype(int), 0, H - 1)
                            sat_adu = float(getattr(P, 'saturation_adu', 60000.0))
                            ok = data[iy, ix] < sat_adu
                            sat_star_count = int((~ok).sum())
                            cand = cand[ok]
                            e_cand = e_cand[ok]

                        positions = [tuple(map(float, p)) for p in cand]
                        if len(cand):
                            elong_map = {
                                (float(x), float(y)): float(e)
                                for (x, y), e in zip(cand, e_cand)
                            }
                        n_sources = len(positions)

                        for src in cat:
                            try:
                                area = float(src.area.value)
                                fwhm_est = 2.0 * np.sqrt(area / np.pi)
                                fwhm_values.append(fwhm_est)
                            except Exception:
                                pass

                    # Optional DAO refine to reject hot pixels (cutout-based for speed)
                    # Also collect DAO statistics for each source
                    dao_primary = detect_engine == "dao"
                    dao_refine_active = bool(dao_refine_enable) and detect_engine != "sep"
                    if dao_refine_active or dao_primary:
                        if self._stop_requested:
                            return filename, None
                        self.worker_status.emit(worker_id, short_name, "DAO refine", 65)

                        if dao_primary:
                            # DAO as primary detector: full image scan
                            try:
                                daofind = DAOStarFinder(
                                    fwhm=fwhm_prelim,
                                    threshold=threshold,
                                    sharplo=dao_sharp_lo,
                                    sharphi=dao_sharp_hi,
                                    roundlo=dao_round_lo,
                                    roundhi=dao_round_hi
                                )
                                dao_cat = daofind(data_sub)
                                if dao_cat is not None and len(dao_cat) > 0:
                                    positions = []
                                    try:
                                        round1 = np.asarray(dao_cat["roundness1"], float)
                                        round2 = np.asarray(dao_cat["roundness2"], float)
                                        round_vals = np.maximum(np.abs(round1), np.abs(round2))
                                        if len(round_vals):
                                            median_roundness = float(np.nanmedian(round_vals))
                                    except Exception:
                                        median_roundness = np.nan
                                    try:
                                        sat_adu = float(getattr(P, 'saturation_adu', 60000.0))
                                        if "peak" in dao_cat.colnames:
                                            peaks = np.asarray(dao_cat["peak"], float)
                                            sat_star_count = int(np.sum(peaks >= sat_adu))
                                    except Exception:
                                        sat_star_count = 0
                                    for i, (x, y) in enumerate(zip(dao_cat['xcentroid'], dao_cat['ycentroid'])):
                                        xf, yf = float(x), float(y)
                                        positions.append((xf, yf))
                                        try:
                                            source_dao_info[(xf, yf)] = {
                                                'sharpness': float(dao_cat['sharpness'][i]),
                                                'roundness1': float(dao_cat['roundness1'][i]),
                                                'roundness2': float(dao_cat['roundness2'][i]),
                                                'peak': float(dao_cat['peak'][i]),
                                                'flux': float(dao_cat['flux'][i]),
                                            }
                                        except Exception:
                                            pass
                                    n_sources = len(positions)
                                    detect_method = "dao"
                                else:
                                    detect_method = "none"
                            except Exception:
                                detect_method = "none"
                        elif positions:
                            # DAO refine: single full-image pass + batched KDTree match.
                            filtered = []

                            daofind = DAOStarFinder(
                                fwhm=fwhm_prelim,
                                threshold=threshold * 0.8,
                                sharplo=dao_sharp_lo,
                                sharphi=dao_sharp_hi,
                                roundlo=dao_round_lo,
                                roundhi=dao_round_hi
                            )

                            if self._stop_requested:
                                return filename, None
                            try:
                                dao_full = daofind(data_sub)
                            except Exception:
                                dao_full = None

                            if dao_full is not None and len(dao_full) > 0:
                                dao_xy = np.c_[
                                    np.asarray(dao_full['xcentroid'], float),
                                    np.asarray(dao_full['ycentroid'], float)
                                ]
                                src_xy = np.array(positions, float)
                                dists, _idxs = KDTree(dao_xy).query(src_xy, k=1)
                                keep_dao = dists <= dao_match_tol
                                for _si in range(len(positions)):
                                    if not keep_dao[_si]:
                                        continue
                                    xf, yf = float(positions[_si][0]), float(positions[_si][1])
                                    filtered.append((xf, yf))
                                    _di = int(_idxs[_si])
                                    try:
                                        source_dao_info[(xf, yf)] = {
                                            'sharpness': float(dao_full['sharpness'][_di]),
                                            'roundness1': float(dao_full['roundness1'][_di]),
                                            'roundness2': float(dao_full['roundness2'][_di]),
                                            'peak': float(dao_full['peak'][_di]),
                                            'flux': float(dao_full['flux'][_di]),
                                        }
                                    except Exception:
                                        pass

                            positions = filtered
                            n_sources = len(positions)
                            detect_method = "segm+dao"

                    # Peak assist (Jupyter-style)
                    added_peak = 0
                    peak_positions = []
                    if peak_enable and detect_engine != "sep":
                        if len(positions) < peak_skip_if_nsrc_ge:
                            if self._stop_requested:
                                return filename, None
                            self.worker_status.emit(worker_id, short_name, "Peak assist", 70)
                            try:
                                peak_xy = peak_pass_candidates(data_sub, bkg_median, fwhm_seed)
                                if positions:
                                    base = np.array(positions, float)
                                    tree_base = KDTree(base)
                                    for (x, y) in peak_xy:
                                        d, _ = tree_base.query([x, y], k=1)
                                        if d >= iso_min_sep * 0.6:
                                            positions.append((float(x), float(y)))
                                            added_peak += 1
                                            peak_positions.append((float(x), float(y)))
                                else:
                                    positions = [(float(x), float(y)) for (x, y) in peak_xy]
                                    added_peak = len(positions)
                                    peak_positions = [(float(x), float(y)) for (x, y) in peak_xy]
                                if added_peak > 0:
                                    if detect_method == "segm+dao":
                                        detect_method = "segm+dao+peak"
                                    elif detect_method.startswith("segm"):
                                        detect_method = "segm+peak"
                                    else:
                                        detect_method = "peak"
                            except Exception:
                                pass
                        n_sources = len(positions)

                    # Calculate median FWHM (radial)
                    fwhm_median = finite_nanmedian(fwhm_values, default=0.0) if fwhm_values else 0.0
                    fwhm_qc_max = int(getattr(P, 'fwhm_qc_max_sources', 40))
                    fwhm_measure_max = int(getattr(P, 'fwhm_measure_max', 25))
                    fwhm_dr = float(getattr(P, 'fwhm_dr', 0.5))
                    ann_sigma = float(getattr(P, 'annulus_sigma_clip', 3.0))
                    ann_maxiter = int(getattr(P, 'fitsky_max_iter', 5))
                    ann_in_scale = float(getattr(P, 'fitsky_annulus_scale', 4.0))
                    ann_out_scale = float(getattr(P, 'fitsky_dannulus_scale', 2.0))
                    ann_min_gap_px = float(getattr(P, 'annulus_min_gap_px', 6.0))
                    ann_min_width_px = float(getattr(P, 'annulus_min_width_px', 12.0))

                    try:
                        if positions:
                            pts = np.array(positions, float)
                            h, w = data.shape
                            ix = np.clip(np.round(pts[:, 0]).astype(int), 0, w - 1)
                            iy = np.clip(np.round(pts[:, 1]).astype(int), 0, h - 1)
                            vals = np.where(np.isfinite(data_sub[iy, ix]), data_sub[iy, ix], data[iy, ix])
                            order = np.argsort(vals)[::-1]
                            pts = pts[order]
                            n_use = min(int(fwhm_qc_max), int(fwhm_measure_max), len(pts))
                            fwhm_values = []
                            for (x, y) in pts[:n_use]:
                                if self._stop_requested:
                                    return filename, None
                                xy_key = (float(x), float(y))
                                cached_fwhm = fwhm_by_xy.get(xy_key)
                                if cached_fwhm is not None and np.isfinite(cached_fwhm):
                                    fwhm_values.append(float(cached_fwhm))
                                    continue
                                rc = refine_centroid(data, x, y, fwhm_seed)
                                if rc is None:
                                    continue
                                xc, yc, _, _ = rc
                                fpx, _, _ = radial_fwhm(
                                    data, xc, yc, fwhm_seed,
                                    ann_in_scale, ann_out_scale,
                                    ann_min_gap_px, ann_min_width_px,
                                    sigma=ann_sigma, maxiters=ann_maxiter, dr=fwhm_dr
                                )
                                if np.isfinite(fpx):
                                    fpx = float(fpx)
                                    fwhm_values.append(fpx)
                                    fwhm_by_xy[xy_key] = fpx
                            if fwhm_values:
                                fwhm_median = finite_nanmedian(fwhm_values, default=0.0)
                    except Exception:
                        fwhm_median = finite_nanmedian(fwhm_values, default=0.0) if fwhm_values else 0.0
                    try:
                        pixscale = float(getattr(P, 'pixel_scale_arcsec', 0.4))
                    except Exception:
                        pixscale = 0.4
                    if not np.isfinite(pixscale) or pixscale <= 0:
                        pixscale = 0.4
                    fwhm_arcsec = fwhm_median * pixscale

                    result = {
                        'n_sources': n_sources,
                        'positions': positions,
                        'peak_positions': peak_positions,
                        'fwhm_px': fwhm_median,
                        'fwhm_arcsec': fwhm_arcsec,
                        'bkg_median': bkg_median,
                        'bkg_rms': bkg_rms,
                        'filter': filt,  # Preserve original case
                        'threshold': threshold,
                        'sigma_used': nsig,
                        'detect_method': detect_method,
                        'median_elongation': median_elongation,
                        'median_roundness': median_roundness,
                        'sat_star_count': sat_star_count,
                        **source_sig,
                    }

                    # Stage 6: Saving
                    self.worker_status.emit(worker_id, short_name, "Saving", 90)
                    if self._stop_requested:
                        return filename, None

                    # Save to cache
                    cache_file = self.cache_dir / f"detect_{filename}.json"
                    step4_file = step4_out / f"detect_{filename}.json"
                    payload = {
                        'n_sources': n_sources,
                        'fwhm_px': fwhm_median,
                        'fwhm_arcsec': fwhm_arcsec,
                        'bkg_median': bkg_median,
                        'bkg_rms': bkg_rms,
                        'filter': filt,
                        'threshold': threshold,
                        'sigma_used': nsig,
                        'detect_method': detect_method,
                        'peak_added': added_peak,
                        'median_elongation': median_elongation,
                        'median_roundness': median_roundness,
                        'sat_star_count': sat_star_count,
                        **source_sig,
                    }
                    with open(cache_file, 'w') as f:
                        json.dump(payload, f)
                    copy_if_different(cache_file, step4_file)

                    # Save positions with extended info (DAO stats, FWHM, peak value)
                    pos_file = self.cache_dir / f"detect_{filename}.csv"
                    step4_pos = step4_out / f"detect_{filename}.csv"
                    if positions:
                        source_records = []
                        peak_set = {(float(px), float(py)) for px, py in peak_positions}
                        measure_all_source_fwhm = bool(getattr(P, 'fwhm_measure_all_sources', False))
                        for idx, (x, y) in enumerate(positions):
                            if self._stop_requested:
                                return filename, None
                            record = {
                                'det_uid': idx,       # 0-indexed, per-frame detection key
                                'id': idx + 1,
                                'x': x,
                                'y': y,
                            }
                            if elong_map:
                                record['elongation'] = elong_map.get((x, y), np.nan)
                            dao_info = source_dao_info.get((x, y), {})
                            record['sharpness'] = dao_info.get('sharpness', np.nan)
                            record['roundness1'] = dao_info.get('roundness1', np.nan)
                            record['roundness2'] = dao_info.get('roundness2', np.nan)
                            try:
                                r1 = float(record['roundness1'])
                                r2 = float(record['roundness2'])
                                if np.isfinite(r1) or np.isfinite(r2):
                                    r1 = r1 if np.isfinite(r1) else 0.0
                                    r2 = r2 if np.isfinite(r2) else 0.0
                                    record['roundness'] = max(abs(r1), abs(r2))
                            except Exception:
                                pass
                            record['dao_peak'] = dao_info.get('peak', np.nan)
                            record['dao_flux'] = dao_info.get('flux', np.nan)
                            if 'sep_flag' in dao_info:
                                record['sep_flag'] = dao_info.get('sep_flag', 0)

                            try:
                                ix, iy = int(round(x)), int(round(y))
                                h_img, w_img = data.shape
                                in_bounds = 0 <= iy < h_img and 0 <= ix < w_img
                                record['peak_adu'] = float(data[iy, ix]) if in_bounds else np.nan

                                r_in_flag = max(ann_in_scale * fwhm_seed, ann_min_gap_px, 12.0)
                                r_out_flag = max(r_in_flag + ann_out_scale * fwhm_seed, r_in_flag + ann_min_width_px, 20.0)
                                near_edge = (
                                    (float(x) - r_out_flag) < 0
                                    or (float(y) - r_out_flag) < 0
                                    or (float(x) + r_out_flag) >= w_img
                                    or (float(y) + r_out_flag) >= h_img
                                )
                                saturated = np.isfinite(record['peak_adu']) and record['peak_adu'] >= float(getattr(P, 'saturation_adu', 60000.0))

                                xy_key = (float(x), float(y))
                                fwhm_src = fwhm_by_xy.get(xy_key, np.nan)
                                if not np.isfinite(fwhm_src) and measure_all_source_fwhm:
                                    rc = refine_centroid(data, x, y, fwhm_seed)
                                    if rc is not None:
                                        xc, yc, _, _ = rc
                                        fwhm_src, _, _ = radial_fwhm(
                                            data, xc, yc, fwhm_seed,
                                            ann_in_scale, ann_out_scale,
                                            ann_min_gap_px, ann_min_width_px,
                                            sigma=ann_sigma, maxiters=ann_maxiter, dr=fwhm_dr
                                        )
                                        if fwhm_src is not None and np.isfinite(fwhm_src):
                                            fwhm_by_xy[xy_key] = float(fwhm_src)
                                record['fwhm_px'] = float(fwhm_src) if np.isfinite(fwhm_src) else np.nan
                                if not np.isfinite(record['fwhm_px']):
                                    record['fwhm_status'] = "not_measured"
                                elif saturated:
                                    record['fwhm_status'] = "saturated"
                                elif near_edge:
                                    record['fwhm_status'] = "edge"
                                else:
                                    record['fwhm_status'] = "ok"
                            except Exception:
                                record['peak_adu'] = np.nan
                                record['fwhm_px'] = np.nan
                                record['fwhm_status'] = "failed"

                            base_source_type = 'sep' if detect_engine == 'sep' else 'segm'
                            record['source_type'] = 'peak' if (x, y) in peak_set else base_source_type
                            source_records.append(record)

                        df_sources = compute_source_quality(
                            pd.DataFrame(source_records),
                            frame_fwhm_px=fwhm_median,
                            image_shape=data.shape,
                            params=P,
                        )
                        quality_summary = {
                            'n_anchor_candidates': int(df_sources.get('anchor_candidate', pd.Series(dtype=bool)).sum()),
                            'n_apcorr_candidates': int(df_sources.get('apcorr_candidate', pd.Series(dtype=bool)).sum()),
                            'n_epsf_candidates': int(df_sources.get('epsf_candidate', pd.Series(dtype=bool)).sum()),
                            'n_psf_seed_candidates': int(df_sources.get('psf_seed_candidate', pd.Series(dtype=bool)).sum()),
                            'quality_score_median': finite_nanmedian(
                                pd.to_numeric(df_sources.get('quality_score', pd.Series(dtype=float)), errors='coerce').to_numpy(float),
                                default=np.nan,
                            ),
                        }
                        result.update(quality_summary)
                        payload.update(quality_summary)
                        with open(cache_file, 'w') as f:
                            json.dump(payload, f)
                        copy_if_different(cache_file, step4_file)
                        df_sources.to_csv(pos_file, index=False)
                        copy_if_different(pos_file, step4_pos)
                    if peak_positions:
                        peak_file = self.cache_dir / f"detect_peak_{filename}.csv"
                        step4_peak = step4_out / f"detect_peak_{filename}.csv"
                        np.savetxt(peak_file, peak_positions, delimiter=',',
                                  header='x,y', comments='')
                        copy_if_different(peak_file, step4_peak)

                    self.worker_status.emit(worker_id, short_name, "Done", 100)

                    return filename, result

                except Exception as e:
                    print(f"[DetectionWorker] Error in detect_single({filename}): {e}")
                    return filename, {'error': str(e)}

            # Run detection in parallel
            print(f"[DetectionWorker] Starting ThreadPoolExecutor with {max_workers} workers")
            self._executor = ThreadPoolExecutor(max_workers=max_workers)
            futures = {self._executor.submit(detect_single, f): f for f in self.file_list}
            try:
                completed_count = 0

                for future in as_completed(futures):
                    if self._stop_requested:
                        for f in futures:
                            f.cancel()
                        break

                    filename = futures[future]
                    completed_count += 1

                    # Count active workers (pending futures)
                    active = sum(1 for f in futures if not f.done())

                    try:
                        fname, result = future.result()
                        if result is not None:
                            if 'error' in result:
                                self.error.emit(fname, result['error'])
                            else:
                                results[fname] = result
                                self.file_done.emit(fname, result)
                    except Exception as e:
                        self.error.emit(filename, str(e))

                    self.progress.emit(completed_count, total, filename, active)
            finally:
                self._executor.shutdown(wait=True, cancel_futures=True)
                self._executor = None

            # Summary
            print(f"[DetectionWorker] Completed. {len(results)} results")
            if self._stop_requested:
                summary = {
                    'stopped': True,
                    'total_files': len(results),
                    'total_sources': sum(r['n_sources'] for r in results.values()) if results else 0,
                }
            elif results:
                total_sources = sum(r['n_sources'] for r in results.values())
                avg_sources = total_sources / len(results)
                fwhm_values_arc = [r['fwhm_arcsec'] for r in results.values() if r['fwhm_arcsec'] > 0]
                fwhm_values_px = [r['fwhm_px'] for r in results.values() if r.get('fwhm_px', 0.0) > 0]
                median_fwhm = float(np.median(fwhm_values_arc)) if fwhm_values_arc else 0.0
                median_fwhm_px = float(np.median(fwhm_values_px)) if fwhm_values_px else float("nan")

                summary = {
                    'total_files': len(results),
                    'total_sources': total_sources,
                    'avg_sources': avg_sources,
                    'median_fwhm_arcsec': median_fwhm,
                    'median_fwhm_px': median_fwhm_px,
                }
            else:
                summary = {}

            self.finished.emit(summary)

        except Exception as e:
            print(f"[DetectionWorker] FATAL ERROR: {e}")
            import traceback
            traceback.print_exc()
            self.error.emit("WORKER", str(e))
            self.finished.emit({})


class QCInspectionPanel(QWidget):
    """QC inspection panel for per-frame quality checks."""
    ALL_FILTER_LABEL = "All"

    def __init__(self, parent_window):
        super().__init__(parent_window)
        self.parent_window = parent_window
        self.params = parent_window.params
        self.file_manager = parent_window.file_manager
        self.file_list = []
        self.frame_df = pd.DataFrame()
        self.exclude_reasons = {}
        self.pending_candidates = {}
        self._header_cache = {}
        self._scatter_map = {}
        self._pending_state = None
        self.current_filter = None
        self._selected_fname = None
        self._build_ui()

    def _build_ui(self):
        layout = QHBoxLayout(self)
        self.setFocusPolicy(Qt.StrongFocus)

        # Left: controls
        control_box = QGroupBox("QC Controls")
        control_layout = QVBoxLayout(control_box)

        filter_row = QHBoxLayout()
        filter_row.addWidget(QLabel("Filter:"))
        self.filter_combo = QComboBox()
        self.filter_combo.currentIndexChanged.connect(self._on_filter_changed)
        filter_row.addWidget(self.filter_combo)
        control_layout.addLayout(filter_row)

        xmode_row = QHBoxLayout()
        xmode_row.addWidget(QLabel("X axis:"))
        self.xmode_combo = QComboBox()
        self.xmode_combo.addItems(["Auto", "Airmass", "Time", "Index"])
        self.xmode_combo.currentIndexChanged.connect(self.update_plots)
        xmode_row.addWidget(self.xmode_combo)
        control_layout.addLayout(xmode_row)

        z_group = QGroupBox("Auto QC (robust z)")
        z_layout = QFormLayout(z_group)
        self.sky_z_spin = QDoubleSpinBox()
        self.sky_z_spin.setRange(1.0, 10.0)
        self.sky_z_spin.setSingleStep(0.5)
        self.sky_z_spin.setValue(4.0)
        z_layout.addRow("Sky z (high):", self.sky_z_spin)

        self.fwhm_z_spin = QDoubleSpinBox()
        self.fwhm_z_spin.setRange(1.0, 10.0)
        self.fwhm_z_spin.setSingleStep(0.5)
        self.fwhm_z_spin.setValue(4.0)
        z_layout.addRow("FWHM z (high):", self.fwhm_z_spin)

        self.nsrc_z_spin = QDoubleSpinBox()
        self.nsrc_z_spin.setRange(1.0, 10.0)
        self.nsrc_z_spin.setSingleStep(0.5)
        self.nsrc_z_spin.setValue(4.0)
        z_layout.addRow("Nsrc z (low):", self.nsrc_z_spin)

        control_layout.addWidget(z_group)

        elong_group = QGroupBox("Elongation Filter (absolute)")
        elong_layout = QFormLayout(elong_group)
        self.elong_thresh_spin = QDoubleSpinBox()
        self.elong_thresh_spin.setRange(1.0, 5.0)
        self.elong_thresh_spin.setSingleStep(0.05)
        self.elong_thresh_spin.setDecimals(2)
        self.elong_thresh_spin.setValue(1.3)
        elong_layout.addRow("elong > (exclude):", self.elong_thresh_spin)
        self.btn_elong_filter = QPushButton("Auto-filter elong")
        self.btn_elong_filter.setToolTip(
            "Exclude all frames whose median elongation exceeds the threshold.\n"
            "Works per filter band. Useful for catching tracking-trailed frames."
        )
        self.btn_elong_filter.clicked.connect(self.auto_filter_elong)
        elong_layout.addRow(self.btn_elong_filter)
        control_layout.addWidget(elong_group)

        btn_row = QHBoxLayout()
        self.btn_find = QPushButton("Find Outliers (z)")
        self.btn_find.clicked.connect(self.find_outliers)
        btn_row.addWidget(self.btn_find)
        self.btn_apply = QPushButton("Exclude Candidates")
        self.btn_apply.clicked.connect(self.apply_candidates)
        btn_row.addWidget(self.btn_apply)
        control_layout.addLayout(btn_row)

        btn_row2 = QHBoxLayout()
        self.btn_reset = QPushButton("Clear Exclusions")
        self.btn_reset.clicked.connect(self.reset_filter_exclusions)
        btn_row2.addWidget(self.btn_reset)
        self.btn_save = QPushButton("Save")
        self.btn_save.clicked.connect(self.save_frame_quality)
        btn_row2.addWidget(self.btn_save)
        control_layout.addLayout(btn_row2)

        self.warning_label = QLabel("")
        self.warning_label.setStyleSheet("QLabel { color: #D32F2F; }")
        self.warning_label.setWordWrap(True)
        control_layout.addWidget(self.warning_label)

        self.hotkey_label = QLabel("Click a point to select. D = exclude, A = include (undo)")
        self.hotkey_label.setStyleSheet("QLabel { color: #455A64; }")
        self.hotkey_label.setWordWrap(True)
        control_layout.addWidget(self.hotkey_label)

        info_group = QGroupBox("Selected Frame")
        info_layout = QVBoxLayout(info_group)
        self.selected_label = QLabel("Click a point to inspect frame details.")
        self.selected_label.setWordWrap(True)
        info_layout.addWidget(self.selected_label)
        self.btn_open_frame = QPushButton("Open in Detection Tab")
        self.btn_open_frame.clicked.connect(self._open_selected_frame)
        info_layout.addWidget(self.btn_open_frame)
        control_layout.addWidget(info_group)

        cand_group = QGroupBox("Outlier Candidates")
        cand_layout = QVBoxLayout(cand_group)
        self.cand_table = QTableWidget()
        self.cand_table.setColumnCount(5)
        self.cand_table.setHorizontalHeaderLabels(["File", "Sky z", "FWHM z", "Nsrc z", "Reasons"])
        self.cand_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.cand_table.horizontalHeader().setStretchLastSection(True)
        self.cand_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.cand_table.setMinimumHeight(160)
        self.cand_table.cellClicked.connect(self._on_candidate_clicked)
        cand_layout.addWidget(self.cand_table)
        control_layout.addWidget(cand_group)

        summary_group = QGroupBox("QC Summary")
        summary_layout = QVBoxLayout(summary_group)
        self.summary_text = QTextEdit()
        self.summary_text.setReadOnly(True)
        self.summary_text.setStyleSheet("QTextEdit { font-family: monospace; font-size: 9pt; }")
        summary_layout.addWidget(self.summary_text)
        control_layout.addWidget(summary_group)

        layout.addWidget(control_box)

        # Right: plots
        plot_box = QGroupBox("Inspection Plots")
        plot_layout = QVBoxLayout(plot_box)
        self.plot_status = QLabel("No data loaded.")
        plot_layout.addWidget(self.plot_status)
        self.fig = Figure(figsize=(6, 6))
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setFocusPolicy(Qt.StrongFocus)
        self.ax_sky = self.fig.add_subplot(2, 1, 1)
        self.ax_fwhm = self.fig.add_subplot(2, 1, 2)
        self.canvas.mpl_connect("pick_event", self._on_pick)
        self.canvas.mpl_connect("key_press_event", self._on_keypress)
        plot_layout.addWidget(self.canvas)
        layout.addWidget(plot_box, stretch=1)

    def _safe_float(self, value, default=np.nan):
        try:
            if value is None:
                return default
            return float(value)
        except Exception:
            return default

    def _resolve_fits_path(self, fname: str, use_cropped: bool) -> Optional[Path]:
        if use_cropped:
            cropped_dir = step2_cropped_dir(self.params.P.result_dir)
            cand = cropped_dir / fname
            if cand.exists():
                return cand
        try:
            return Path(self.params.get_file_path(fname))
        except Exception:
            return None

    def _parse_time_value(self, header: fits.Header) -> tuple[float, str]:
        for key in ("JD", "JULIAN", "BJD", "HJD", "MJD-OBS", "MJD"):
            if key in header:
                val = self._safe_float(header.get(key), np.nan)
                if np.isfinite(val):
                    return float(val), key
        date_obs = header.get("DATE-OBS") or header.get("DATE")
        time_obs = header.get("TIME-OBS") or header.get("UTC") or header.get("UT")
        if date_obs:
            dt_str = str(date_obs).strip()
            if "T" not in dt_str and time_obs:
                dt_str = f"{dt_str}T{str(time_obs).strip()}"
            try:
                t = Time(dt_str, format="isot", scale="utc")
                return float(t.jd), "JD"
            except Exception:
                pass
            try:
                t = Time(dt_str, scale="utc")
                return float(t.jd), "JD"
            except Exception:
                pass
        return np.nan, "index"

    def _load_header_meta(self, fname: str, use_cropped: bool) -> dict:
        if fname in self._header_cache:
            return self._header_cache[fname]
        meta = {"airmass": np.nan, "time_val": np.nan, "time_src": "index"}
        path = self._resolve_fits_path(fname, use_cropped)
        if path and path.exists():
            try:
                with fits.open(path) as hdul:
                    h = hdul[0].header
                info = compute_airmass_from_header(
                    h,
                    float(getattr(self.params.P, "site_lat_deg", np.nan)),
                    float(getattr(self.params.P, "site_lon_deg", np.nan)),
                    float(getattr(self.params.P, "site_alt_m", 0.0)),
                    float(getattr(self.params.P, "site_tz_offset_hours", 0.0)),
                    formula=getattr(self.params.P, "airmass_formula", None),
                )
                meta["airmass"] = self._safe_float(info.get("airmass"), np.nan)
                tval, tsrc = self._parse_time_value(h)
                meta["time_val"] = tval
                meta["time_src"] = tsrc
            except Exception:
                pass
        self._header_cache[fname] = meta
        return meta

    def _load_detect_meta(self, fname: str) -> dict:
        if hasattr(self.parent_window, "_pick_detection_cache"):
            try:
                payload, _, _, _ = self.parent_window._pick_detection_cache(fname, previous=False)
                if isinstance(payload, dict):
                    return payload
            except Exception:
                pass
        cache_dir = self.params.P.cache_dir
        step4_out = step4_dir(self.params.P.result_dir)
        cache_file = cache_dir / f"detect_{fname}.json"
        if not cache_file.exists():
            alt = step4_out / f"detect_{fname}.json"
            if alt.exists():
                cache_file = alt
        if not cache_file.exists():
            return {}
        try:
            return json.loads(cache_file.read_text(encoding="utf-8"))
        except Exception:
            return {}

    def load_frames(self, detection_results: dict, file_list: list, use_cropped: bool) -> None:
        self.exclude_reasons = {}
        self.pending_candidates = {}
        rows = []
        self.file_list = list(file_list)
        for idx, fname in enumerate(self.file_list):
            meta = detection_results.get(fname)
            if meta is None:
                meta = self._load_detect_meta(fname)
            elif "median_elongation" not in meta or "median_roundness" not in meta:
                disk_meta = self._load_detect_meta(fname)
                if disk_meta:
                    merged = dict(disk_meta)
                    merged.update(meta)
                    meta = merged
            if not meta:
                continue
            hmeta = self._load_header_meta(fname, use_cropped)
            rows.append({
                "file": fname,
                "filter": str(meta.get("filter", "") or "").strip(),
                "time_index": idx,
                "time_val": hmeta.get("time_val", np.nan),
                "time_src": hmeta.get("time_src", "index"),
                "airmass": hmeta.get("airmass", np.nan),
                "sky_med": self._safe_float(meta.get("bkg_median"), np.nan),
                "sky_sigma": self._safe_float(meta.get("bkg_rms"), np.nan),
                "fwhm_med": self._safe_float(meta.get("fwhm_px"), np.nan),
                "n_sources": int(meta.get("n_sources", 0) or 0),
                "elong_med": self._safe_float(meta.get("median_elongation"), np.nan),
                "round_med": self._safe_float(meta.get("median_roundness"), np.nan),
            })

        self.frame_df = pd.DataFrame(rows)
        self._refresh_filter_list()
        self._apply_exclusions_from_file()
        self._apply_pending_state()
        self.update_plots()
        self.update_summary()


    def _refresh_filter_list(self):
        self.filter_combo.blockSignals(True)
        self.filter_combo.clear()
        filters = sorted(self.frame_df.get("filter", pd.Series([""])).fillna("").astype(str).unique().tolist())
        if not filters:
            filters = [""]
        self.filter_combo.addItem(self.ALL_FILTER_LABEL)
        for f in filters:
            if not f:
                continue
            self.filter_combo.addItem(f)
        self.filter_combo.setCurrentIndex(0)
        self.filter_combo.blockSignals(False)
        self.current_filter = None

    def _on_filter_changed(self, index: int) -> None:
        if index < 0:
            return
        selected = self.filter_combo.currentText()
        if selected == self.ALL_FILTER_LABEL:
            self.current_filter = None
        else:
            self.current_filter = selected
        self.pending_candidates = {}
        self.cand_table.setRowCount(0)
        self.warning_label.setText("")
        self._refresh_qc_views(force_draw=True)

    def _toggle_exclusion(self, fname: str) -> None:
        reasons = set(self.exclude_reasons.get(fname, set()))
        if reasons:
            self.exclude_reasons[fname] = set()
        else:
            self.exclude_reasons[fname] = {"manual"}
        self._refresh_qc_views(force_draw=True)

    def _on_pick(self, event):
        artist = event.artist
        if artist not in self._scatter_map:
            return
        indices = getattr(event, "ind", None)
        if indices is None or len(indices) == 0:
            return
        fname = self._scatter_map[artist][int(indices[0])]
        self._show_frame_info(fname)
        self.update_plots()
        self.setFocus()
        self.canvas.setFocus()

    def _on_keypress(self, event) -> None:
        key = (getattr(event, "key", "") or "").lower()
        if key not in ("a", "d"):
            return
        fname = getattr(self, "_selected_fname", None)
        if not fname:
            return
        # D = exclude (제외), A = include (되돌리기)
        if key == "d":
            self.exclude_reasons[fname] = {"manual"}
            self.warning_label.setText(f"Excluded: {fname}")
        else:  # key == "a"
            self.exclude_reasons[fname] = set()
            self.warning_label.setText(f"Included: {fname}")
        self._refresh_qc_views(force_draw=True)

    def _show_frame_info(self, fname: str) -> None:
        if self.frame_df.empty:
            return
        row = self.frame_df[self.frame_df["file"] == fname]
        if row.empty:
            return
        r = row.iloc[0]
        elong_val = r.get("elong_med", np.nan)
        elong_str = f"{float(elong_val):.3f}" if np.isfinite(float(elong_val)) else "N/A"
        self.selected_label.setText(
            f"{fname}\n"
            f"filter={r.get('filter','')}, airmass={r.get('airmass', np.nan):.3f}\n"
            f"sky={r.get('sky_med', np.nan):.2f}, fwhm={r.get('fwhm_med', np.nan):.2f}, "
            f"elong={elong_str}, n_sources={int(r.get('n_sources', 0))}"
        )
        self.selected_label.repaint()
        self._selected_fname = fname

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_A:
            self._on_keypress(type("evt", (), {"key": "a"})())
            return
        if event.key() == Qt.Key_D:
            self._on_keypress(type("evt", (), {"key": "d"})())
            return
        super().keyPressEvent(event)

    def _open_selected_frame(self):
        fname = getattr(self, "_selected_fname", None)
        if not fname:
            return
        if hasattr(self.parent_window, "show_frame_in_detection_tab"):
            self.parent_window.show_frame_in_detection_tab(fname)

    def _robust_z(self, values: np.ndarray) -> np.ndarray:
        med, mad = robust_median_mad(values)
        if not np.isfinite(mad) or mad == 0:
            return np.zeros_like(values, dtype=float)
        return 0.6745 * (values - med) / mad

    def _subset_df(self) -> pd.DataFrame:
        if self.frame_df.empty:
            return self.frame_df
        if self.current_filter:
            return self.frame_df[self.frame_df["filter"] == self.current_filter].copy()
        return self.frame_df.copy()

    def find_outliers(self):
        self.pending_candidates = {}
        self.cand_table.setRowCount(0)
        df = self._subset_df()
        if df.empty:
            return
        self.warning_label.setText("")
        sky_th = float(self.sky_z_spin.value())
        fwhm_th = float(self.fwhm_z_spin.value())
        nsrc_th = float(self.nsrc_z_spin.value())

        def _accumulate(df_in: pd.DataFrame) -> None:
            sky = df_in["sky_med"].to_numpy(float)
            fwhm = df_in["fwhm_med"].to_numpy(float)
            nsrc = df_in["n_sources"].to_numpy(float)
            z_sky = self._robust_z(sky)
            z_fwhm = self._robust_z(fwhm)
            z_nsrc = self._robust_z(nsrc)
            for idx, row in enumerate(df_in.itertuples(index=False)):
                reasons = []
                if np.isfinite(z_sky[idx]) and z_sky[idx] > sky_th:
                    reasons.append("sky_outlier")
                if np.isfinite(z_fwhm[idx]) and z_fwhm[idx] > fwhm_th:
                    reasons.append("fwhm_outlier")
                if np.isfinite(z_nsrc[idx]) and z_nsrc[idx] < -nsrc_th:
                    reasons.append("low_nsrc")
                if reasons:
                    fname = row.file
                    self.pending_candidates[fname] = {
                        "sky_z": z_sky[idx],
                        "fwhm_z": z_fwhm[idx],
                        "nsrc_z": z_nsrc[idx],
                        "reasons": reasons,
                    }

        if self.current_filter:
            if len(df) < 10:
                self.warning_label.setText("Too few frames for auto QC (need >=10).")
                return
            _accumulate(df)
        else:
            warn_filters = []
            for filt, grp in df.groupby("filter"):
                if len(grp) < 10:
                    warn_filters.append(filt or "(none)")
                    continue
                _accumulate(grp)
            if warn_filters:
                self.warning_label.setText(
                    "Auto QC skipped (too few frames): " + ", ".join(warn_filters)
                )

        for fname, info in self.pending_candidates.items():
            row_idx = self.cand_table.rowCount()
            self.cand_table.insertRow(row_idx)
            self.cand_table.setItem(row_idx, 0, QTableWidgetItem(str(fname)))
            self.cand_table.setItem(row_idx, 1, QTableWidgetItem(f"{info['sky_z']:.2f}"))
            self.cand_table.setItem(row_idx, 2, QTableWidgetItem(f"{info['fwhm_z']:.2f}"))
            self.cand_table.setItem(row_idx, 3, QTableWidgetItem(f"{info['nsrc_z']:.2f}"))
            self.cand_table.setItem(row_idx, 4, QTableWidgetItem(",".join(info["reasons"])))

        if not self.pending_candidates:
            self.warning_label.setText("No outlier candidates found.")
        else:
            first_candidate = next(iter(self.pending_candidates.keys()))
            self._show_frame_info(first_candidate)
        self._refresh_qc_views(force_draw=True)

    def auto_filter_elong(self):
        """Exclude frames whose median elongation exceeds the absolute threshold."""
        df = self._subset_df()
        if df.empty or "elong_med" not in df.columns:
            self.warning_label.setText("No elongation data available.")
            return
        thresh = float(self.elong_thresh_spin.value())
        elong_vals = df["elong_med"].to_numpy(float)
        n_finite = int(np.isfinite(elong_vals).sum())
        if n_finite == 0:
            self.warning_label.setText(
                "elong_med is N/A for all frames — re-run Step 4 detection to refresh cache."
            )
            return
        n_excluded = 0
        for _, row in df.iterrows():
            elong = float(row.get("elong_med", np.nan))
            fname = str(row["file"])
            if np.isfinite(elong) and elong > thresh:
                self.exclude_reasons.setdefault(fname, set()).add("high_elong")
                n_excluded += 1
        if n_excluded:
            self.warning_label.setText(
                f"Auto-filter: excluded {n_excluded} frame(s) with elong > {thresh:.2f} "
                f"(checked {n_finite}/{len(df)} frames). Click Save to persist."
            )
        else:
            self.warning_label.setText(
                f"No frames with elong > {thresh:.2f} "
                f"(max elong = {np.nanmax(elong_vals):.3f}, n={n_finite}/{len(df)})."
            )
        self._refresh_qc_views(force_draw=True)

    def apply_candidates(self):
        if not self.pending_candidates:
            return
        n_candidates = len(self.pending_candidates)
        for fname, info in self.pending_candidates.items():
            self.exclude_reasons.setdefault(fname, set()).update(info.get("reasons", []))
        self.pending_candidates = {}
        self.cand_table.setRowCount(0)
        self.warning_label.setText(f"Excluded {n_candidates} candidate frame(s). Click Save to persist.")
        self._refresh_qc_views(force_draw=True)

    def _on_candidate_clicked(self, row: int, col: int) -> None:
        if row < 0:
            return
        item = self.cand_table.item(row, 0)
        if not item:
            return
        fname = item.text().strip()
        if not fname:
            return
        self._show_frame_info(fname)
        self._ensure_visible_x(fname)
        self.update_plots()

    def reset_filter_exclusions(self):
        df = self._subset_df()
        for fname in df["file"].tolist():
            self.exclude_reasons[fname] = set()
        self._refresh_qc_views(force_draw=True)

    def _apply_exclusions_from_file(self):
        _rd = self.params.P.result_dir
        fq_path = step4_dir(_rd) / "frame_quality.csv"
        if not fq_path.exists():
            return
        try:
            dfq = pd.read_csv(fq_path)
        except Exception:
            return
        if "file" not in dfq.columns:
            return
        for _, row in dfq.iterrows():
            fname = str(row.get("file", ""))
            passed = is_passed_value(row.get("passed", True), default=True)
            if not fname:
                continue
            if passed:
                self.exclude_reasons[fname] = set()
                continue
            reasons = set()
            reason_str = row.get("exclude_reason", "")
            if isinstance(reason_str, str) and reason_str.strip():
                reasons.update([r.strip() for r in reason_str.split(",") if r.strip()])
            if not reasons:
                reasons.add("manual")
            self.exclude_reasons[fname] = reasons

    def _apply_pending_state(self):
        if not self._pending_state:
            return
        state = self._pending_state
        self._pending_state = None
        self.exclude_reasons = {
            str(k): set(v) for k, v in (state.get("exclude_reasons", {}) or {}).items()
        }
        if "qc_filter" in state and state["qc_filter"]:
            filt = state["qc_filter"]
            idx = self.filter_combo.findText(filt)
            if idx >= 0:
                self.filter_combo.setCurrentIndex(idx)
        if "sky_z" in state:
            self.sky_z_spin.setValue(float(state["sky_z"]))
        if "fwhm_z" in state:
            self.fwhm_z_spin.setValue(float(state["fwhm_z"]))
        if "nsrc_z" in state:
            self.nsrc_z_spin.setValue(float(state["nsrc_z"]))
        if "x_mode" in state:
            idx = self.xmode_combo.findText(state["x_mode"])
            if idx >= 0:
                self.xmode_combo.setCurrentIndex(idx)

    def export_state(self) -> dict:
        return {
            "exclude_reasons": {k: sorted(list(v)) for k, v in self.exclude_reasons.items()},
            "qc_filter": self.current_filter,
            "sky_z": self.sky_z_spin.value(),
            "fwhm_z": self.fwhm_z_spin.value(),
            "nsrc_z": self.nsrc_z_spin.value(),
            "x_mode": self.xmode_combo.currentText(),
        }

    def restore_state(self, state: Optional[dict]) -> None:
        if not state:
            return
        if self.frame_df.empty:
            self._pending_state = state
        else:
            self._pending_state = state
            self._apply_pending_state()
            self.update_plots()
            self.update_summary()

    def _build_quality_df(self) -> pd.DataFrame:
        if self.frame_df.empty:
            return pd.DataFrame()
        df = self.frame_df.copy()
        reasons = []
        passed = []
        for fname in df["file"].tolist():
            r = self.exclude_reasons.get(fname, set())
            reasons.append(",".join(sorted(r)) if r else "")
            passed.append(len(r) == 0)
        df["exclude_reason"] = reasons
        df["passed"] = passed
        return df

    def write_frame_quality_csv(self):
        df = self._build_quality_df()
        if df.empty:
            return None
        out_dir = step4_dir(self.params.P.result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        gain = self._safe_float(getattr(self.params.P, "gain_e_per_adu", np.nan), np.nan)
        if np.isfinite(gain) and "sky_sigma" in df.columns:
            df["sky_sigma_med_e"] = df["sky_sigma"] * gain
        if "sky_sigma" in df.columns:
            df.rename(columns={
                "sky_sigma": "sky_sigma_med_adu",
            }, inplace=True)
        out_path = out_dir / "frame_quality.csv"
        df.to_csv(out_path, index=False)
        return out_path

    def save_frame_quality(self):
        if self.write_frame_quality_csv() is None:
            return
        self._apply_pipeline_flags()
        self.warning_label.setText("QC saved and applied.")
        self.update_summary()

    def apply_to_pipeline(self):
        if self.write_frame_quality_csv() is not None:
            self.warning_label.setText("QC saved and applied.")
            self.update_summary()
        self._apply_pipeline_flags()

    def _apply_pipeline_flags(self):
        self.params.P.wcs_require_qc_pass = True
        self.params.P.phot_use_qc_pass_only = True
        if hasattr(self.params.P, "idmatch_use_qc_pass_only"):
            self.params.P.idmatch_use_qc_pass_only = True
        if hasattr(self.parent_window, "persist_params"):
            self.parent_window.persist_params()
        if hasattr(self.parent_window, "save_state"):
            self.parent_window.save_state()

    def _refresh_qc_views(self, force_draw: bool = False) -> None:
        """Refresh all QC views after in-memory inclusion/exclusion changes."""
        self.update_summary()
        self.update_plots(force_draw=force_draw)
        self.cand_table.viewport().update()
        self.summary_text.viewport().update()
        self.plot_status.repaint()
        if force_draw:
            self.canvas.flush_events()
            self.canvas.repaint()
            QApplication.processEvents()

    def update_plots(self, force_draw: bool = False):
        if self.frame_df.empty:
            self.plot_status.setText("No data loaded.")
            self.fig.clear()
            self.ax_sky = self.fig.add_subplot(2, 1, 1)
            self.ax_fwhm = self.fig.add_subplot(2, 1, 2)
            if force_draw:
                self.canvas.draw()
            else:
                self.canvas.draw_idle()
            return
        df = self._subset_df()
        if df.empty:
            return
        x_mode = self.xmode_combo.currentText()
        if x_mode == "Airmass":
            x_vals = df["airmass"].to_numpy(float)
            x_label = "Airmass"
        elif x_mode == "Time":
            x_vals = df["time_val"].to_numpy(float)
            x_label = "Time"
        elif x_mode == "Index":
            x_vals = df["time_index"].to_numpy(float)
            x_label = "Index"
        else:
            airmass = df["airmass"].to_numpy(float)
            time_vals = df["time_val"].to_numpy(float)
            n_air = int(np.isfinite(airmass).sum())
            n_time = int(np.isfinite(time_vals).sum())
            if n_air == 0 and n_time == 0:
                x_vals = df["time_index"].to_numpy(float)
                x_label = "Index"
            elif n_air >= n_time:
                x_vals = airmass
                x_label = "Airmass"
            else:
                x_vals = time_vals
                x_label = "Time"

        excluded = np.array([len(self.exclude_reasons.get(f, set())) > 0 for f in df["file"].tolist()])
        pending = np.array([f in self.pending_candidates for f in df["file"].tolist()])
        included = ~(excluded | pending)

        self.fig.clear()
        self.ax_sky = self.fig.add_subplot(2, 1, 1)
        self.ax_fwhm = self.fig.add_subplot(2, 1, 2)
        self._scatter_map = {}

        def _scatter(ax, x, y, mask, color, marker, label, size=28, alpha=0.8, edge=None):
            finite = mask & np.isfinite(x) & np.isfinite(y)
            xs = x[finite]
            ys = y[finite]
            if len(xs) == 0:
                return None
            sc = ax.scatter(
                xs, ys, s=size, color=color, marker=marker, alpha=alpha, picker=5,
                label=label, edgecolors=edge
            )
            files = [f for f, m in zip(df["file"].tolist(), finite) if m]
            self._scatter_map[sc] = files
            return sc

        finite_x = np.isfinite(x_vals)
        finite_sky = np.isfinite(df["sky_med"].to_numpy(float))
        finite_fwhm = np.isfinite(df["fwhm_med"].to_numpy(float))
        hidden_sky = int(np.sum(pending & ~(finite_x & finite_sky)))
        hidden_fwhm = int(np.sum(pending & ~(finite_x & finite_fwhm)))

        sky_handles = []
        sc = _scatter(self.ax_sky, x_vals, df["sky_med"].to_numpy(float), included,
                      "#000000", "o", "included", size=22, alpha=0.9)
        if sc:
            sky_handles.append(sc)
        sc = _scatter(self.ax_sky, x_vals, df["sky_med"].to_numpy(float), pending,
                      "#E53935", "o", "outlier", size=58, alpha=0.9, edge="#212121")
        if sc:
            sky_handles.append(sc)
        sc = _scatter(self.ax_sky, x_vals, df["sky_med"].to_numpy(float), excluded,
                      "#9E9E9E", "x", "excluded", size=40, alpha=0.9)
        if sc:
            sky_handles.append(sc)
        self.ax_sky.set_ylabel("sky_med")
        self.ax_sky.set_xlabel(x_label)
        self.ax_sky.grid(True, alpha=0.2)
        if sky_handles:
            self.ax_sky.legend(handles=sky_handles, loc="best", fontsize=8, frameon=False)

        fwhm_handles = []
        sc = _scatter(self.ax_fwhm, x_vals, df["fwhm_med"].to_numpy(float), included,
                      "#000000", "o", "included", size=22, alpha=0.9)
        if sc:
            fwhm_handles.append(sc)
        sc = _scatter(self.ax_fwhm, x_vals, df["fwhm_med"].to_numpy(float), pending,
                      "#E53935", "o", "outlier", size=58, alpha=0.9, edge="#212121")
        if sc:
            fwhm_handles.append(sc)
        sc = _scatter(self.ax_fwhm, x_vals, df["fwhm_med"].to_numpy(float), excluded,
                      "#9E9E9E", "x", "excluded", size=40, alpha=0.9)
        if sc:
            fwhm_handles.append(sc)
        self.ax_fwhm.set_ylabel("fwhm_med")
        self.ax_fwhm.set_xlabel(x_label)
        self.ax_fwhm.grid(True, alpha=0.2)
        fwhm_vals = df["fwhm_med"].to_numpy(float)
        fwhm_med = np.nanmedian(fwhm_vals)
        fwhm_mad = np.nanmedian(np.abs(fwhm_vals - fwhm_med))
        if np.isfinite(fwhm_med) and np.isfinite(fwhm_mad) and fwhm_mad > 0:
            fwhm_cut = fwhm_med + (self.fwhm_z_spin.value() * fwhm_mad / 0.6745)
            self.ax_fwhm.axhline(fwhm_cut, color="#E53935", linestyle="--", linewidth=1.2, alpha=0.8)
        if fwhm_handles:
            self.ax_fwhm.legend(handles=fwhm_handles, loc="best", fontsize=8, frameon=False)

        sel = getattr(self, "_selected_fname", None)
        if sel:
            row = df[df["file"] == sel]
            if not row.empty:
                r = row.iloc[0]
                x_sel = float(x_vals[df["file"].tolist().index(sel)])
                self.ax_sky.scatter(
                    [x_sel], [r["sky_med"]], s=220, facecolors="none",
                    edgecolors="red", linewidths=2.0, zorder=20
                )
                self.ax_fwhm.scatter(
                    [x_sel], [r["fwhm_med"]], s=220, facecolors="none",
                    edgecolors="red", linewidths=2.0, zorder=20
                )

        n_total = len(df)
        n_exc = int(excluded.sum())
        rate = (n_exc / n_total * 100.0) if n_total else 0.0
        filter_label = self.current_filter or "all"
        hidden_note = ""
        if hidden_sky or hidden_fwhm:
            hidden_note = f" | hidden(outlier) sky={hidden_sky} fwhm={hidden_fwhm}"
        self.plot_status.setText(
            f"Filter={filter_label} | frames={n_total} | excluded={n_exc} ({rate:.1f}%) | "
            f"outlier=red dot, selected=red circle, excluded=gray x{hidden_note}"
        )
        self.fig.tight_layout()
        if force_draw:
            self.canvas.draw()
            QApplication.processEvents()
        else:
            self.canvas.draw_idle()

    def _ensure_visible_x(self, fname: str) -> None:
        if self.frame_df.empty or not fname:
            return
        row = self.frame_df[self.frame_df["file"] == fname]
        if row.empty:
            return
        r = row.iloc[0]
        x_mode = self.xmode_combo.currentText()
        airmass = self._safe_float(r.get("airmass"), np.nan)
        time_val = self._safe_float(r.get("time_val"), np.nan)
        if x_mode == "Airmass":
            if not np.isfinite(airmass):
                if np.isfinite(time_val):
                    self.xmode_combo.setCurrentText("Time")
                    self.warning_label.setText("X axis switched to Time (missing AIRMASS).")
                else:
                    self.xmode_combo.setCurrentText("Index")
                    self.warning_label.setText("X axis switched to Index (missing AIRMASS/Time).")
        elif x_mode == "Time":
            if not np.isfinite(time_val):
                if np.isfinite(airmass):
                    self.xmode_combo.setCurrentText("Airmass")
                    self.warning_label.setText("X axis switched to Airmass (missing Time).")
                else:
                    self.xmode_combo.setCurrentText("Index")
                    self.warning_label.setText("X axis switched to Index (missing AIRMASS/Time).")
        elif x_mode == "Auto":
            if not (np.isfinite(airmass) or np.isfinite(time_val)):
                self.xmode_combo.setCurrentText("Index")
                self.warning_label.setText("X axis switched to Index (missing AIRMASS/Time).")

    def update_summary(self):
        df = self._subset_df()
        if df.empty:
            self.summary_text.setText("No data.")
            return
        reasons_count = {"sky_outlier": 0, "fwhm_outlier": 0, "low_nsrc": 0, "high_elong": 0, "manual": 0}
        excluded_files = []
        for fname in df["file"].tolist():
            r = self.exclude_reasons.get(fname, set())
            if r:
                excluded_files.append(fname)
            for key in reasons_count:
                if key in r:
                    reasons_count[key] += 1
        n_total = len(df)
        n_exc = len(excluded_files)
        rate = (n_exc / n_total * 100.0) if n_total else 0.0

        sky_top = df.sort_values("sky_med", ascending=False).head(10)
        fwhm_top = df.sort_values("fwhm_med", ascending=False).head(10)
        nsrc_low = df.sort_values("n_sources", ascending=True).head(10)

        lines = [
            f"Excluded: {n_exc}/{n_total} ({rate:.1f}%)",
            f"Reasons: sky={reasons_count['sky_outlier']} "
            f"fwhm={reasons_count['fwhm_outlier']} "
            f"nsrc={reasons_count['low_nsrc']} "
            f"elong={reasons_count['high_elong']} "
            f"manual={reasons_count['manual']}",
            "",
            "Top sky_med:",
        ]
        for _, r in sky_top.iterrows():
            lines.append(f"  {r['file']}  {r['sky_med']:.2f}")
        lines.append("")
        lines.append("Top fwhm_med:")
        for _, r in fwhm_top.iterrows():
            lines.append(f"  {r['file']}  {r['fwhm_med']:.2f}")
        lines.append("")
        lines.append("Low n_sources:")
        for _, r in nsrc_low.iterrows():
            lines.append(f"  {r['file']}  {int(r['n_sources'])}")

        self.summary_text.setText("\n".join(lines))


class SourceDetectionWindow(StepWindowBase):
    """
    Step 4: Source Detection
    Parallel detection with segmentation/deblending
    """

    def __init__(self, params, file_manager, project_state, main_window):
        """Initialize source detection window"""
        self.file_manager = file_manager
        self.detection_worker = None
        self.detection_results = {}
        self.stop_requested = False
        self._detect_cache_mgr: StepCacheManager | None = None

        # Image data
        self.current_filename = None
        self.image_data = None
        self.header = None

        # Image viewer
        self._viewer: FITSViewerWidget | None = None
        self._selected_star_pos = None
        self._fits_cache: OrderedDict = OrderedDict()  # filename -> (image_data, header), LRU order
        self._FITS_CACHE_SIZE = max(3, int(getattr(params.P, "step4_fits_cache_size", 6)))

        # File list
        self.file_list = []
        self.use_cropped = False

        # Filter-sigma mapping (flexible, user-defined)
        self.filter_sigma_map = {}
        self.log_window = None
        self.worker_progress_bars = {}
        self.worker_last_status = {}
        self._resume_cache_active = False

        # Stretch plot window (2D Plot)
        self.stretch_plot_dialog = None
        self.stretch_plot_canvas = None
        self.stretch_plot_ax = None
        self.stretch_plot_fig = None
        self.stretch_plot_info_label = None
        self._stretch_vmin = None
        self._stretch_vmax = None
        self._stretch_data_range = None
        self._stretch_dragging = False
        self._stretch_drag_target = None
        self._stretch_marker_min_line = None
        self._stretch_marker_max_line = None

        # Initialize base class
        super().__init__(
            step_index=3,  # 0-based index
            step_name="Source Detection",
            params=params,
            project_state=project_state,
            main_window=main_window
        )

        # Setup step-specific UI
        self.setup_step_ui()

        # Load filter sigma map from parameters
        self.load_filter_sigma_map()
        self.params.P.detect_mode = _get_detect_mode_from_params(self.params.P)

        # Restore state
        self.restore_state()

    def load_filter_sigma_map(self):
        """Load filter-specific sigma values from parameters"""
        P = self.params.P
        raw = getattr(P, '_raw', {})

        sigma_by_filter = getattr(P, "detect_sigma_by_filter", None)
        if not isinstance(sigma_by_filter, dict):
            sigma_by_filter = raw.get("detect_sigma_by_filter", {})
        if isinstance(sigma_by_filter, dict):
            for filt, val in sigma_by_filter.items():
                try:
                    self.filter_sigma_map[str(filt)] = float(val)
                except Exception:
                    pass

        # Read all detect_sigma_* patterns
        for key, val in raw.items():
            if key.startswith('detect_sigma_') and key != 'detect_sigma':
                filt = key.replace('detect_sigma_', '')
                try:
                    self.filter_sigma_map[filt] = float(val)
                except Exception:
                    pass

        # Also check uppercase parameter names.
        if hasattr(P, 'detect_sigma_g') and P.detect_sigma_g:
            self.filter_sigma_map['g'] = float(P.detect_sigma_g)
        if hasattr(P, 'detect_sigma_r') and P.detect_sigma_r:
            self.filter_sigma_map['r'] = float(P.detect_sigma_r)
        if hasattr(P, 'detect_sigma_i') and P.detect_sigma_i:
            self.filter_sigma_map['i'] = float(P.detect_sigma_i)

    # === Detection Mode Preset System ===

    def _selected_detect_mode(self) -> str:
        combo = getattr(self, "param_detect_mode", None)
        if combo is None:
            return _get_detect_mode_from_params(self.params.P)
        idx = combo.currentIndex()
        data = combo.itemData(idx)
        return _normalize_detect_mode(data if data not in (None, "") else combo.currentText())

    def _selected_detect_engine(self) -> str:
        combo = getattr(self, "param_engine", None)
        if combo is None:
            return _normalize_detect_engine(getattr(self.params.P, "detect_engine", "segm"))
        idx = combo.currentIndex()
        data = combo.itemData(idx)
        return _normalize_detect_engine(data if data not in (None, "") else combo.currentText())

    def _sync_engine_dependent_dialog_state(self):
        engine = self._selected_detect_engine()
        is_sep = engine == "sep"

        # SEP already performs fast C-level background, extraction, deblending,
        # and cheap morphology measurements. Full-frame DAO/peak passes would
        # defeat the point of selecting the fast engine.
        if hasattr(self, "param_dao_enable"):
            if is_sep:
                self.param_dao_enable.setChecked(False)
            self.param_dao_enable.setEnabled(not is_sep)
        for name in (
            "param_dao_fwhm",
            "param_dao_sharp_lo",
            "param_dao_sharp_hi",
            "param_dao_round_lo",
            "param_dao_round_hi",
            "param_dao_match_tol",
        ):
            widget = getattr(self, name, None)
            if widget is not None:
                widget.setEnabled(not is_sep)
        if hasattr(self, "param_peak_enable"):
            if is_sep:
                self.param_peak_enable.setChecked(False)
            self.param_peak_enable.setEnabled(not is_sep)

    def _apply_detect_mode_preset_to_dialog(self, mode: str) -> bool:
        mode_key = _normalize_detect_mode(mode)
        if mode_key == "custom":
            return False
        preset = _get_detect_mode_preset(mode_key)
        if hasattr(self, "param_detect_sigma"):
            self.param_detect_sigma.setValue(float(preset["detect_sigma"]))
        if hasattr(self, "param_minarea"):
            self.param_minarea.setValue(int(preset["minarea_pix"]))
        if hasattr(self, "param_deblend"):
            self.param_deblend.setChecked(bool(preset["deblend_enable"]))
        if hasattr(self, "param_deblend_nthresh"):
            self.param_deblend_nthresh.setValue(int(preset["deblend_nthresh"]))
        if hasattr(self, "param_deblend_cont"):
            self.param_deblend_cont.setValue(float(preset["deblend_cont"]))
        if hasattr(self, "param_deblend_max_labels"):
            self.param_deblend_max_labels.setValue(int(preset["deblend_max_labels"]))
        if hasattr(self, "param_deblend_label_hard_max"):
            self.param_deblend_label_hard_max.setValue(int(preset["deblend_label_hard_max"]))
        if hasattr(self, "param_dao_enable"):
            self.param_dao_enable.setChecked(bool(preset["dao_refine_enable"]))
        if hasattr(self, "param_peak_enable"):
            self.param_peak_enable.setChecked(bool(preset["peak_pass_enable"]))
        if hasattr(self, "param_peak_nsigma"):
            self.param_peak_nsigma.setValue(float(preset["peak_nsigma"]))
        if hasattr(self, "param_peak_max_add"):
            self.param_peak_max_add.setValue(int(preset["peak_max_add"]))
        self._sync_engine_dependent_dialog_state()
        return True

    def _update_detect_mode_ui_state(self):
        mode_key = self._selected_detect_mode()
        if hasattr(self, "param_detect_mode_apply"):
            self.param_detect_mode_apply.setEnabled(mode_key != "custom")

    def _on_apply_detect_mode_clicked(self):
        mode_key = self._selected_detect_mode()
        self._apply_detect_mode_preset_to_dialog(mode_key)
        self._update_detect_mode_ui_state()

    # === Detection Cache Validation ===

    def _resolve_source_path_for_file(self, filename: str) -> Optional[Path]:
        if self.use_cropped:
            cropped_dir = step2_cropped_dir(self.params.P.result_dir)
            cand = cropped_dir / filename
            if cand.exists():
                return cand
        try:
            cand = Path(self.params.get_file_path(filename))
            if cand.exists():
                return cand
        except Exception:
            pass
        return None

    def _source_signature_for_file(self, filename: str) -> Optional[dict]:
        src_path = self._resolve_source_path_for_file(filename)
        if src_path is None or not src_path.exists():
            return None
        try:
            st = src_path.stat()
        except Exception:
            return None
        return {
            "source_path": norm_path_key(src_path),
            "source_use_cropped": bool(self.use_cropped),
            "source_size": int(st.st_size),
            "source_mtime_ns": int(st.st_mtime_ns),
        }

    def _is_detection_cache_compatible(self, filename: str, payload: dict, meta_path: Path) -> bool:
        if not isinstance(payload, dict):
            return False
        sig = self._source_signature_for_file(filename)
        if sig is None:
            return False
        return detection_cache_signature_matches(
            payload,
            sig,
            min_schema=DETECTION_CACHE_SCHEMA_VERSION,
            current_engine=getattr(self.params.P, "detect_engine", "segm"),
            allow_mtime_drift=True,
        )

    def _pick_detection_cache(self, filename: str, previous: bool = False):
        cache_dir = self.params.P.cache_dir
        step4_out = step4_dir(self.params.P.result_dir)
        if previous:
            prefix = "detect_prev_"
            peak_prefix = "detect_prev_peak_"
            candidates = [cache_dir / f"{prefix}{filename}.json"]
        else:
            prefix = "detect_"
            peak_prefix = "detect_peak_"
            candidates = [cache_dir / f"{prefix}{filename}.json", step4_out / f"{prefix}{filename}.json"]
            candidates = [p for p in candidates if p.exists()]
            candidates.sort(key=lambda p: p.stat().st_mtime_ns, reverse=True)
        for meta_path in candidates:
            if not meta_path.exists():
                continue
            try:
                payload = json.loads(meta_path.read_text(encoding="utf-8"))
            except Exception:
                continue
            if not self._is_detection_cache_compatible(filename, payload, meta_path):
                continue
            base_dir = meta_path.parent
            pos_path = base_dir / f"{prefix}{filename}.csv"
            peak_path = base_dir / f"{peak_prefix}{filename}.csv"
            return payload, pos_path, peak_path, meta_path
        return None, None, None, None

    # === FITS Cache ===

    def _load_fits_cached(self, filename: str):
        if filename in self._fits_cache:
            self._fits_cache.move_to_end(filename)
            return self._fits_cache[filename]

        file_path = self._resolve_fits_path(filename, self.use_cropped)
        if file_path is None:
            raise FileNotFoundError(f"FITS not found: {filename}")
        with fits.open(file_path, memmap=False) as hdul:
            data = np.asarray(hdul[0].data, dtype=np.float32)
            header = hdul[0].header.copy()

        self._fits_cache[filename] = (data, header)
        while len(self._fits_cache) > self._FITS_CACHE_SIZE:
            self._fits_cache.popitem(last=False)
        return self._fits_cache[filename]

    def scan_filters_from_files(self):
        """Scan FITS files to detect which filters are actually present"""
        filters_found = set()

        # Determine data directory
        if self.use_cropped:
            data_dir = step2_cropped_dir(self.params.P.result_dir)
        else:
            data_dir = self.params.P.data_dir

        # Sample files (first 20 or all if less)
        sample_files = self.file_list[:min(20, len(self.file_list))]

        for filename in sample_files:
            try:
                if self.use_cropped:
                    file_path = step2_cropped_dir(self.params.P.result_dir) / filename
                else:
                    file_path = self.params.get_file_path(filename)
                with fits.open(file_path) as hdul:
                    filt = hdul[0].header.get('FILTER', '').strip()
                    if filt:
                        filters_found.add(filt)
            except Exception:
                pass

        return sorted(filters_found)

    def setup_step_ui(self):
        """Setup step-specific UI components"""

        # Tabs
        self.tabs = QTabWidget()
        self.detect_tab = QWidget()
        self.detect_layout = QVBoxLayout(self.detect_tab)
        self.qc_tab = QWidget()
        self.qc_layout = QVBoxLayout(self.qc_tab)
        self.tabs.addTab(self.detect_tab, "Detection")
        self.tabs.addTab(self.qc_tab, "QC")
        self.content_layout.addWidget(self.tabs)

        # === Info Label ===
        info_label = QLabel(
            "Detect sources in all images using segmentation algorithm.\n"
            "Results are cached for subsequent steps. Mouse: Wheel to zoom | Right-click drag to pan"
        )
        info_label.setStyleSheet("QLabel { background-color: #E8F5E9; padding: 10px; border-radius: 5px; }")
        self.detect_layout.addWidget(info_label)

        # === Control Bar ===
        control_layout = QHBoxLayout()

        btn_params = create_parameter_button("Detection Parameters")
        btn_params.clicked.connect(self.open_parameters_dialog)
        control_layout.addWidget(btn_params)

        btn_clear_cache = create_cache_action_button("Clear Detection Cache")
        btn_clear_cache.clicked.connect(self.clear_detection_cache)
        control_layout.addWidget(btn_clear_cache)

        self.chk_resume_cache = create_detection_cache_checkbox(
            bool(getattr(self.params.P, "resume_mode", True))
            and not bool(getattr(self.params.P, "force_redetect", False)),
            "When enabled, skip frames with compatible Step 4 detect_*.json cache. "
            "Disable to force source detection for every frame.",
        )
        control_layout.addWidget(self.chk_resume_cache)

        control_layout.addStretch()

        self.run_bar = RunControlBar(
            "Run Detection", "Log & Workers",
            run_cb=self.run_detection,
            stop_cb=self.stop_detection,
            log_cb=self.show_log_window,
        )
        control_layout.addWidget(self.run_bar)
        self.btn_run = self.run_bar.btn_run
        self.btn_stop = self.run_bar.btn_stop

        self.detect_layout.addLayout(control_layout)

        # === Progress Bar ===
        progress_layout = QHBoxLayout()
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        progress_layout.addWidget(self.progress_bar)

        self.progress_label = QLabel("Ready")
        self.progress_label.setMinimumWidth(350)
        progress_layout.addWidget(self.progress_label)

        self.detect_layout.addLayout(progress_layout)

        # === Main Splitter ===
        main_splitter = QSplitter(Qt.Horizontal)

        # Left: Image Viewer
        viewer_group = QGroupBox("Preview")
        viewer_layout = QVBoxLayout(viewer_group)

        # File selector row
        file_layout = QHBoxLayout()
        file_layout.addWidget(QLabel("File:"))
        self.file_combo = QComboBox()
        self.file_combo.currentIndexChanged.connect(self.on_file_changed)
        file_layout.addWidget(self.file_combo)

        btn_load = QPushButton("Load")
        btn_load.clicked.connect(self.load_and_display)
        file_layout.addWidget(btn_load)

        self.chk_overlay = QCheckBox("Show Sources")
        self.chk_overlay.setChecked(True)
        self.chk_overlay.stateChanged.connect(self.update_overlay)
        file_layout.addWidget(self.chk_overlay)

        file_layout.addStretch()
        viewer_layout.addLayout(file_layout)

        # 2D Plot button (top bar)
        btn_2d_plot = QPushButton("2D Plot")
        btn_2d_plot.setStyleSheet("QPushButton { background-color: #FF9800; color: white; font-weight: bold; }")
        btn_2d_plot.clicked.connect(self.open_stretch_plot)
        file_layout.addWidget(btn_2d_plot)

        # GPU-accelerated FITS viewer
        self._viewer = FITSViewerWidget(self)
        self._viewer.gl.mouse_pressed.connect(self._on_viewer_click)
        viewer_layout.addWidget(self._viewer)

        main_splitter.addWidget(viewer_group)

        # Right: Results
        results_group = QGroupBox("Detection Results")
        results_layout = QVBoxLayout(results_group)

        # Summary stats
        self.summary_label = QLabel("No detection run yet")
        self.summary_label.setStyleSheet("QLabel { font-family: monospace; padding: 10px; background-color: #f5f5f5; }")
        self.summary_label.setWordWrap(True)
        results_layout.addWidget(self.summary_label)

        # Results table - updated columns
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(6)
        self.results_table.setHorizontalHeaderLabels(['File', 'N', 'FWHM', 'Bkg', 'Filt', 'Sig'])
        results_header = self.results_table.horizontalHeader()
        results_header.setStretchLastSection(False)
        results_header.setSectionResizeMode(0, QHeaderView.Stretch)
        for col, width in ((1, 42), (2, 92), (3, 52), (4, 42), (5, 42)):
            results_header.setSectionResizeMode(col, QHeaderView.Fixed)
            self.results_table.setColumnWidth(col, width)
        self.results_table.setWordWrap(False)
        self.results_table.setTextElideMode(Qt.ElideMiddle)
        self.results_table.verticalHeader().setDefaultSectionSize(22)
        self.results_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.results_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.results_table.cellClicked.connect(self.on_table_cell_clicked)
        results_layout.addWidget(self.results_table)

        # === Selected Star Info Panel ===
        star_info_group = QGroupBox("Selected Star Info")
        star_info_group.setStyleSheet("QGroupBox { font-weight: bold; }")
        star_info_layout = QVBoxLayout(star_info_group)

        self.star_info_label = QLabel("Right-click on a star in the image to see details\n(Right-click near a detected source)")
        self.star_info_label.setStyleSheet("""
            QLabel {
                font-family: monospace;
                font-size: 11px;
                padding: 8px;
                background-color: #FFFDE7;
                border: 1px solid #FBC02D;
                border-radius: 4px;
            }
        """)
        self.star_info_label.setWordWrap(True)
        self.star_info_label.setMinimumHeight(180)
        star_info_layout.addWidget(self.star_info_label)

        results_layout.addWidget(star_info_group)

        main_splitter.addWidget(results_group)
        main_splitter.setStretchFactor(0, 2)
        main_splitter.setStretchFactor(1, 1)

        self.detect_layout.addWidget(main_splitter)

        self.setup_log_window()

        self.qc_panel = QCInspectionPanel(self)
        self.qc_layout.addWidget(self.qc_panel)

        # Populate file list
        self.populate_file_list()

    def populate_file_list(self):
        """Populate file combo box"""
        crop_active = crop_is_active(self.params.P.result_dir)
        cropped_dir = step2_cropped_dir(self.params.P.result_dir)

        excluded = getattr(self.file_manager, "excluded_files", set())
        if crop_active and cropped_dir.exists() and list(cropped_dir.glob("*.fit*")):
            files = sorted([f.name for f in cropped_dir.glob("*.fit*")
                            if f.name not in excluded])
            self.use_cropped = True
        else:
            if not self.file_manager.filenames:
                try:
                    self.file_manager.scan_files()
                except Exception:
                    pass
            files = [f for f in self.file_manager.filenames if f not in excluded]
            self.use_cropped = False

        self.file_list = list(files)
        self.file_combo.clear()
        self.file_combo.addItems(files)
        self._fits_cache.clear()
        self.load_cached_results()

    def _refresh_qc_panel(self):
        if hasattr(self, "qc_panel") and self.qc_panel is not None:
            self.qc_panel.load_frames(self.detection_results, self.file_list, self.use_cropped)

    def show_frame_in_detection_tab(self, filename: str) -> None:
        idx = self.file_combo.findText(filename)
        if idx >= 0:
            self.file_combo.setCurrentIndex(idx)
            self.load_and_display(quick_switch=True)
        if hasattr(self, "tabs"):
            self.tabs.setCurrentIndex(0)

    def load_cached_results(self):
        """Load cached detection results from disk"""
        cache_dir = self.params.P.cache_dir
        step4_out = step4_dir(self.params.P.result_dir)
        if not cache_dir.exists() and not step4_out.exists():
            return

        results = {}
        skipped_incompatible = 0
        mgr = self._get_detect_cache_mgr()
        for filename in self.file_list:
            manifest_invalid = False
            # manifest validation (supplemental — missing manifest is OK for backward compat)
            try:
                vr = mgr.validate_key(filename, required_payloads=["detect_json"])
                if vr.manifest is not None and not vr.valid:
                    manifest_invalid = True
            except Exception:
                pass

            data, pos_file, peak_file, cache_file = self._pick_detection_cache(filename, previous=False)
            if data is None:
                # meta exists but all incompatible -> count for log
                maybe_meta = [
                    cache_dir / f"detect_{filename}.json",
                    step4_out / f"detect_{filename}.json",
                ]
                if manifest_invalid or any(p.exists() for p in maybe_meta):
                    skipped_incompatible += 1
                continue
            try:
                positions = []
                peak_positions = []
                if pos_file is not None and pos_file.exists():
                    try:
                        df_pos = pd.read_csv(pos_file)
                        if 'x' in df_pos.columns and 'y' in df_pos.columns:
                            positions = [(row['x'], row['y']) for _, row in df_pos.iterrows()]
                        else:
                            pos_data = np.loadtxt(pos_file, delimiter=',', skiprows=1)
                            if pos_data.ndim == 1:
                                positions = [(pos_data[0], pos_data[1])]
                            else:
                                positions = [(row[0], row[1]) for row in pos_data]
                    except Exception:
                        positions = []
                if peak_file is not None and peak_file.exists():
                    try:
                        peak_positions = np.loadtxt(peak_file, delimiter=',', skiprows=1).tolist()
                        if peak_positions and isinstance(peak_positions[0], float):
                            peak_positions = [peak_positions]
                    except Exception:
                        peak_positions = []

                result = dict(data)
                result.update({
                    'n_sources': data.get('n_sources', 0),
                    'positions': positions,
                    'peak_positions': peak_positions,
                    'fwhm_px': data.get('fwhm_px', 0.0),
                    'fwhm_arcsec': data.get('fwhm_arcsec', 0.0),
                    'bkg_median': data.get('bkg_median', 0.0),
                    'bkg_rms': data.get('bkg_rms', 0.0),
                    'filter': data.get('filter', ''),
                    'threshold': data.get('threshold', 0.0),
                    'sigma_used': data.get('sigma_used', 0.0),
                    'detect_method': data.get('detect_method', 'segm'),
                })
                results[filename] = result
            except Exception:
                continue

        if results:
            self.detection_results = results
            self.populate_results_table()
            self.update_summary_from_results()
            self.log(f"Loaded cached results: {len(results)} files")
            if skipped_incompatible > 0:
                self.log(f"Skipped incompatible detection cache: {skipped_incompatible} files")
            self.update_navigation_buttons()
            self._refresh_qc_panel()
        elif skipped_incompatible > 0:
            self.log(f"Ignored incompatible detection cache: {skipped_incompatible} files")


    def setup_log_window(self):
        """Create the log/workers window"""
        if self.log_window is not None:
            return

        worker_group = QGroupBox("Workers")
        worker_group.setMinimumWidth(430)
        worker_layout = QVBoxLayout(worker_group)
        worker_layout.setContentsMargins(5, 5, 5, 5)
        self.worker_panel = WorkerStatusPanel(worker_group)
        worker_layout.addWidget(self.worker_panel)

        self.log_window = WorkflowLogWindow(
            self,
            "Detection Log & Workers",
            width=900,
            height=500,
            side_widget=worker_group,
        )
        self.log_text = self.log_window.log_text
        self.log_window.show()

    def show_log_window(self):
        """Show log/workers window"""
        if self.log_window is None:
            self.setup_log_window()
        show_raised(self.log_window)

    def clear_worker_status(self):
        """Clear worker status UI"""
        self.worker_progress_bars = {}
        self.worker_last_status = {}
        if hasattr(self, "worker_panel") and self.worker_panel is not None:
            self.worker_panel.clear()

    def on_file_changed(self, index):
        """Handle file selection change"""
        pass

    def _resolve_fits_path(self, fname: str, use_cropped: bool) -> Optional[Path]:
        """Resolve display FITS path for Step4 viewer."""
        if use_cropped:
            cropped_dir = step2_cropped_dir(self.params.P.result_dir)
            cand = cropped_dir / fname
            if cand.exists():
                return cand
        try:
            cand = Path(self.params.get_file_path(fname))
            if cand.exists():
                return cand
        except Exception:
            pass
        cand = self.params.P.data_dir / fname
        if cand.exists():
            return cand
        return None

    def load_and_display(self, quick_switch=False):
        """Load and display selected image"""
        filename = self.file_combo.currentText()
        if not filename:
            return

        try:
            new_data, new_header = self._load_fits_cached(filename)
            self.image_data = new_data
            self.header = new_header
            self.current_filename = filename
            self._stretch_vmin = None
            self._stretch_vmax = None
            self._selected_star_pos = None
            self._update_viewer_image(keep_view=quick_switch)
            self.update_overlay()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load: {str(e)}")

    # === Viewer methods ===

    def _update_viewer_image(self, keep_view=False):
        if self._viewer is None or self.image_data is None:
            return
        self._viewer.set_data(self.image_data)
        self._viewer.auto_stf()
        if not keep_view:
            self._viewer.fit_in_view()

    def _on_viewer_click(self, x, y, btn):
        if btn == int(Qt.RightButton):
            self.select_nearest_star(int(round(x)), int(round(y)))

    def open_stretch_plot(self):
        """Open stretch plot window showing histogram with draggable min/max markers"""
        if self.image_data is None:
            QMessageBox.warning(self, "Warning", "Load an image first")
            return

        if self.stretch_plot_dialog is not None and self.stretch_plot_dialog.isVisible():
            self.stretch_plot_dialog.raise_()
            self.stretch_plot_dialog.activateWindow()
            self._update_stretch_plot()
            return

        self.stretch_plot_dialog = QDialog(self)
        self.stretch_plot_dialog.setWindowTitle("2D Plot - Stretch Control")
        self.stretch_plot_dialog.resize(500, 250)

        layout = QVBoxLayout(self.stretch_plot_dialog)

        self.stretch_plot_info_label = QLabel("Drag min/max markers to adjust stretch")
        self.stretch_plot_info_label.setStyleSheet(
            "QLabel { padding: 5px; background-color: #E3F2FD; border-radius: 3px; }"
        )
        layout.addWidget(self.stretch_plot_info_label)

        self.stretch_plot_fig = Figure(figsize=(6, 2.5))
        self.stretch_plot_canvas = FigureCanvas(self.stretch_plot_fig)
        self.stretch_plot_ax = self.stretch_plot_fig.add_subplot(111)
        self.stretch_plot_fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15, top=0.9)

        self.stretch_plot_canvas.mpl_connect('button_press_event', self._on_stretch_plot_press)
        self.stretch_plot_canvas.mpl_connect('motion_notify_event', self._on_stretch_plot_motion)
        self.stretch_plot_canvas.mpl_connect('button_release_event', self._on_stretch_plot_release)

        layout.addWidget(self.stretch_plot_canvas)

        hint_label = QLabel("Click and drag < > markers to adjust min/max | Changes apply in real-time")
        hint_label.setStyleSheet("QLabel { color: #666; font-size: 10px; }")
        layout.addWidget(hint_label)

        self.stretch_plot_dialog.show()
        self._update_stretch_plot()

    def _update_stretch_plot(self):
        """Update the stretch plot histogram and markers"""
        if self.stretch_plot_ax is None or self.image_data is None:
            return

        ax = self.stretch_plot_ax
        ax.clear()

        data = self.image_data.copy()
        finite_mask = np.isfinite(data)
        if not finite_mask.any():
            return

        flat = data[finite_mask].flatten()

        p_low, p_high = np.percentile(flat, [1, 99])
        display_data = flat[(flat >= p_low) & (flat <= p_high)]
        if len(display_data) == 0:
            display_data = flat

        self._stretch_data_range = (float(p_low), float(p_high))

        if self._stretch_vmin is None or self._stretch_vmax is None:
            if self._viewer is not None:
                shadow, highlight, _ = self._viewer.get_stf_params()
                vmin, vmax = float(shadow), float(highlight)
            else:
                _, median_val, std_val = sigma_clipped_stats(flat, sigma=3.0, maxiters=5)
                vmin = max(float(np.min(flat)), median_val - 2.8 * std_val)
                vmax = min(float(np.max(flat)), np.percentile(flat, 99.9))

            if vmax <= vmin:
                vmin = float(np.min(flat))
                vmax = float(np.max(flat))

            self._stretch_vmin = float(vmin)
            self._stretch_vmax = float(vmax)

        ax.hist(display_data, bins=128, color='#3a6ea5', edgecolor='none', alpha=0.7)
        ax.set_xlim(p_low, p_high)

        vmin = self._stretch_vmin
        vmax = self._stretch_vmax

        vmin_display = max(p_low, min(p_high, vmin))
        vmax_display = max(p_low, min(p_high, vmax))

        self._stretch_marker_min_line = ax.axvline(
            vmin_display, color='#FF5722', linewidth=2, linestyle='-', label=f"Min: {vmin:.1f}"
        )
        self._stretch_marker_max_line = ax.axvline(
            vmax_display, color='#4CAF50', linewidth=2, linestyle='-', label=f"Max: {vmax:.1f}"
        )

        y_max = ax.get_ylim()[1]
        ax.text(vmin_display, y_max * 0.95, '<', color='#FF5722', fontsize=14,
                ha='center', va='top', fontweight='bold')
        ax.text(vmax_display, y_max * 0.95, '>', color='#4CAF50', fontsize=14,
                ha='center', va='top', fontweight='bold')

        ax.set_xlabel('Pixel Value')
        ax.set_ylabel('Count')
        ax.set_title('Image Histogram')
        ax.legend(loc='upper right', fontsize=8)

        if self.stretch_plot_info_label:
            mode_name = self._viewer._mode_combo.currentText() if self._viewer else "STF"
            self.stretch_plot_info_label.setText(
                f"Stretch: {mode_name} | Min: {vmin:.2f} | Max: {vmax:.2f}"
            )

        self.stretch_plot_canvas.draw_idle()

    def _on_stretch_plot_press(self, event):
        """Handle mouse press on stretch plot"""
        if event.inaxes != self.stretch_plot_ax or event.xdata is None:
            return
        if self._stretch_vmin is None or self._stretch_vmax is None:
            return

        x = event.xdata
        dist_to_min = abs(x - self._stretch_vmin)
        dist_to_max = abs(x - self._stretch_vmax)
        self._stretch_drag_target = "min" if dist_to_min < dist_to_max else "max"
        self._stretch_dragging = True

    def _on_stretch_plot_motion(self, event):
        """Handle mouse motion on stretch plot (dragging)"""
        if not self._stretch_dragging or event.xdata is None:
            return

        x = event.xdata
        if self._stretch_drag_target == "min":
            new_val = min(x, self._stretch_vmax - 1)
            self._stretch_vmin = new_val
        else:
            new_val = max(x, self._stretch_vmin + 1)
            self._stretch_vmax = new_val

        self._update_stretch_plot()
        self._apply_custom_stretch()

    def _on_stretch_plot_release(self, event):
        """Handle mouse release on stretch plot"""
        self._stretch_dragging = False
        self._stretch_drag_target = None

    def _apply_custom_stretch(self):
        if self.image_data is None or self._viewer is None:
            return
        if self._stretch_vmin is None or self._stretch_vmax is None:
            return
        vmin, vmax = self._stretch_vmin, self._stretch_vmax
        if vmax <= vmin:
            vmax = vmin + 1
        self._viewer.set_shadow_highlight(vmin, vmax)

    def select_nearest_star(self, click_x, click_y):
        """Find and display info for the nearest detected star"""
        if click_x is None or click_y is None:
            return
        if not self.current_filename:
            return

        cache_dir = self.params.P.cache_dir
        pos_file = cache_dir / f"detect_{self.current_filename}.csv"

        if not pos_file.exists():
            self.star_info_label.setText("No detection data available for this frame")
            return

        try:
            df = pd.read_csv(pos_file)
            if df.empty or 'x' not in df.columns or 'y' not in df.columns:
                self.star_info_label.setText("No sources in detection file")
                return

            distances = np.sqrt((df['x'] - click_x) ** 2 + (df['y'] - click_y) ** 2)
            min_idx = distances.idxmin()
            min_dist = distances[min_idx]

            if min_dist > 20:
                self.star_info_label.setText(f"No star found near click position\n(nearest is {min_dist:.1f} px away)")
                return

            src = df.iloc[min_idx]
            det_id = int(src.get('id', min_idx + 1))
            info_lines = [
                f"{'═' * 32}",
                f"  Detection #{det_id}",
                f"{'═' * 32}",
                f"",
                f"Position:",
                f"  X: {src['x']:.2f} px",
                f"  Y: {src['y']:.2f} px",
                f"",
            ]

            if 'fwhm_px' in src and pd.notna(src['fwhm_px']):
                fwhm_px = float(src['fwhm_px'])
                pixscale = getattr(self.params.P, 'pixel_scale_arcsec', 0.4)
                fwhm_arcsec = fwhm_px * pixscale
                info_lines.append(f"FWHM: {fwhm_px:.2f} px ({fwhm_arcsec:.2f}\")")
            else:
                info_lines.append("FWHM: (measurement failed)")

            if 'peak_adu' in src and pd.notna(src['peak_adu']):
                info_lines.append(f"Peak: {src['peak_adu']:.1f} ADU")

            info_lines.append("")

            if 'sharpness' in src and pd.notna(src['sharpness']):
                info_lines.append("DAO Statistics:")
                info_lines.append(f"  Sharpness:  {src['sharpness']:.4f}")
                if 'roundness1' in src and pd.notna(src['roundness1']):
                    info_lines.append(f"  Roundness1: {src['roundness1']:.4f}")
                if 'roundness2' in src and pd.notna(src['roundness2']):
                    info_lines.append(f"  Roundness2: {src['roundness2']:.4f}")
                if 'dao_peak' in src and pd.notna(src['dao_peak']):
                    info_lines.append(f"  DAO Peak:   {src['dao_peak']:.1f}")
                if 'dao_flux' in src and pd.notna(src['dao_flux']):
                    info_lines.append(f"  DAO Flux:   {src['dao_flux']:.1f}")
            else:
                info_lines.append("DAO Statistics: N/A")
                info_lines.append("  (DAO refine disabled or source")
                info_lines.append("   added via peak-assist)")

            if 'source_type' in src:
                info_lines.append("")
                info_lines.append(f"Source Type: {src['source_type']}")

            self.star_info_label.setText("\n".join(info_lines))
            self.highlight_selected_star(src['x'], src['y'])
        except Exception as e:
            self.star_info_label.setText(f"Error reading source data:\n{str(e)}")

    def highlight_selected_star(self, x, y):
        self._selected_star_pos = (float(x), float(y))
        self.update_overlay()

    def update_overlay(self):
        """Update source overlay on image"""
        if self._viewer is None or self.image_data is None or self.current_filename is None:
            return

        markers = []
        lime = QColor(0, 255, 0, 180)
        cyan = QColor(0, 220, 255, 180)

        if self.chk_overlay.isChecked() and self.current_filename in self.detection_results:
            result = self.detection_results[self.current_filename]
            for (x, y) in result.get('positions', []):
                markers.append(OverlayMarker(col=float(x), row=float(y), radius=6.0, color=lime))
            for (x, y) in result.get('peak_positions', []):
                markers.append(OverlayMarker(col=float(x), row=float(y), radius=7.5, color=cyan))

        if self._selected_star_pos is not None:
            sx, sy = self._selected_star_pos
            markers.append(OverlayMarker(col=sx, row=sy, radius=15.0, color=QColor(255, 220, 0, 220)))

        self._viewer.set_overlay_markers(markers)

    def run_detection(self):
        """Start detection process"""
        if not self.file_list:
            QMessageBox.warning(self, "Warning", "No files to process")
            return

        if self.detection_worker and self.detection_worker.isRunning():
            QMessageBox.information(self, "Detection Running", "Detection is already running.")
            return

        # Prepare
        use_cache = bool(getattr(self, "chk_resume_cache", None) and self.chk_resume_cache.isChecked())
        self._resume_cache_active = use_cache
        self.params.P.resume_mode = use_cache
        self.params.P.force_redetect = not use_cache
        if hasattr(self, "persist_params"):
            self.persist_params()

        if use_cache:
            self.load_cached_results()
            cached = set(self.detection_results.keys())
            pending_files = [f for f in self.file_list if f not in cached]
            if not pending_files:
                self.update_summary_from_results(title="Detection Complete (cache)")
                self.progress_label.setText("Done")
                self.log("All frames already cached. Nothing to do.")
                return
        else:
            pending_files = list(self.file_list)

        if not use_cache:
            self.detection_results = {}
            self._refresh_qc_panel()
            self.results_table.setRowCount(0)
        else:
            self.populate_results_table()
            self._refresh_qc_panel()

        self.log_text.clear()
        self.clear_worker_status()
        self.stop_requested = False
        if use_cache:
            self.log(f"Starting detection (resume). Cached: {len(self.file_list) - len(pending_files)}, Pending: {len(pending_files)}")
        self.log(f"Starting detection on {len(pending_files)} files...")
        self.log(f"Use cropped: {self.use_cropped}")
        self.log(f"Data dir: {self.params.P.data_dir}")
        self.log(f"Filter sigma map: {self.filter_sigma_map}")

        # Ensure cache directory exists
        cache_dir = self.params.P.cache_dir
        cache_dir.mkdir(parents=True, exist_ok=True)
        self.log(f"Cache dir: {cache_dir}")

        # Create worker with filter sigma map
        self.detection_worker = DetectionWorker(
            pending_files,
            self.params,
            self.params.P.data_dir,
            cache_dir,
            self.use_cropped,
            self.filter_sigma_map
        )

        # Connect signals
        self.detection_worker.progress.connect(self.on_progress)
        self.detection_worker.file_done.connect(self.on_file_done)
        self.detection_worker.finished.connect(self.on_detection_finished)
        self.detection_worker.error.connect(self.on_detection_error)
        self.detection_worker.worker_status.connect(self.on_worker_status)

        # Update UI
        self.run_bar.set_running(True)
        self.progress_bar.setValue(0)
        self.progress_bar.setMaximum(len(pending_files))
        self._detect_start_time = time.monotonic()
        self.progress_label.setText(
            progress_status_text(0, len(pending_files), self._detect_start_time, message="Starting...")
        )

        # Start
        self.log("Starting worker thread...")
        self.detection_worker.start()
        self.log("Worker thread started")
        self.show_log_window()

    def stop_detection(self):
        """Stop detection process"""
        if self.detection_worker and self.detection_worker.isRunning():
            if self.stop_requested:
                return
            self.stop_requested = True
            self.run_bar.set_stopping()
            self.progress_label.setText("Stopping...")
            self.log("Stopping...")
            self.detection_worker.stop()

    def on_progress(self, current, total, filename, active_workers):
        """Handle progress update"""
        self.progress_bar.setValue(current)
        self.progress_label.setText(
            progress_status_text(
                current, total, getattr(self, "_detect_start_time", None),
                workers=active_workers, message=filename,
            )
        )

    def on_worker_status(self, worker_id, filename, status, progress):
        """Update worker status panel + log meaningful state changes."""
        if not hasattr(self, "worker_panel") or self.worker_panel is None:
            self.setup_log_window()
        self.worker_panel.update_worker(worker_id, filename, status, progress)

        last = self.worker_last_status.get(worker_id)
        current = (filename, status)
        if last != current:
            short_fname = Path(str(filename)).name
            self.log(f"W{int(worker_id):02d} {short_fname}  {status}")
            self.worker_last_status[worker_id] = current

    def _get_detect_cache_mgr(self) -> StepCacheManager:
        if self._detect_cache_mgr is None:
            self._detect_cache_mgr = StepCacheManager(
                self.params.P.cache_dir, "source_detection", cache_schema_version=1
            )
        return self._detect_cache_mgr

    def _write_detect_manifest(self, filename: str) -> None:
        try:
            from apex.utils.step_paths import step2_cropped_dir, crop_is_active
            cache_dir = Path(self.params.P.cache_dir)
            fits_path = (
                step2_cropped_dir(self.params.P.result_dir) / filename
                if getattr(self, "use_cropped", False)
                else Path(self.params.P.data_dir) / filename
            )
            mgr = self._get_detect_cache_mgr()
            manifest = mgr.build_manifest(
                input_paths=[fits_path] if fits_path.exists() else [],
                payload_paths={
                    "detect_json": cache_dir / f"detect_{filename}.json",
                    "detect_csv": cache_dir / f"detect_{filename}.csv",
                },
            )
            mgr.write_manifest(filename, manifest)
        except Exception:
            pass

    def on_file_done(self, filename, result):
        """Handle single file completion"""
        self.detection_results[filename] = result
        self._write_detect_manifest(filename)

        # Periodic state save every 20 frames (protects against overnight crash)
        if len(self.detection_results) % 20 == 0:
            self.save_state()

        # Add to table - with detection method in FWHM column
        row = self.results_table.rowCount()
        self.results_table.insertRow(row)
        filename_item = QTableWidgetItem(filename)
        filename_item.setToolTip(filename)
        self.results_table.setItem(row, 0, filename_item)
        self.results_table.setItem(row, 1, QTableWidgetItem(str(result['n_sources'])))
        # FWHM with detection method
        fwhm_px = float(result.get("fwhm_px", 0.0))
        fwhm_arcsec = float(result.get("fwhm_arcsec", 0.0))
        fwhm_str = f'{fwhm_arcsec:.2f}" ({fwhm_px:.2f} px; {result.get("detect_method", "segm")})'
        self.results_table.setItem(row, 2, QTableWidgetItem(fwhm_str))
        self.results_table.setItem(row, 3, QTableWidgetItem(f"{result['bkg_median']:.1f}"))
        # Filter - preserve original case from header
        self.results_table.setItem(row, 4, QTableWidgetItem(result['filter']))
        # Sigma used
        self.results_table.setItem(row, 5, QTableWidgetItem(f"{result.get('sigma_used', 3.2):.1f}"))

        self.log(
            f"{filename}: {result['n_sources']} sources, "
            f"FWHM={fwhm_arcsec:.2f}\" ({fwhm_px:.2f} px; {result.get('detect_method', 'segm')}), "
            f"sigma={result.get('sigma_used', 3.2):.1f}"
        )

    def on_detection_error(self, filename, error):
        """Handle detection error"""
        self.log(f"ERROR {filename}: {error}")

    def on_detection_finished(self, summary):
        """Handle detection completion"""
        self.run_bar.set_running(False)
        self.stop_requested = False
        resume_mode = self._resume_cache_active
        self._resume_cache_active = False
        if self.detection_worker:
            self.detection_worker.wait(2000)
            self.detection_worker.deleteLater()
            self.detection_worker = None

        if summary and summary.get('stopped'):
            elapsed_txt = ""
            if hasattr(self, "_detect_start_time"):
                elapsed_txt = f" | elapsed {format_duration(time.monotonic() - self._detect_start_time)}"
            self.summary_label.setText(
                f"Detection Stopped\n"
                f"{'─' * 30}\n"
                f"Files processed: {summary.get('total_files', 0)}\n"
                f"Total sources: {summary.get('total_sources', 0)}"
            )
            self.log("Detection stopped by user")
            self.progress_label.setText(f"Stopped{elapsed_txt}")
        elif summary and not resume_mode:
            median_arc = float(summary.get("median_fwhm_arcsec", np.nan))
            median_px = float(summary.get("median_fwhm_px", np.nan))
            if np.isfinite(median_arc) and np.isfinite(median_px):
                fwhm_note = f'Median FWHM: {median_arc:.2f}" ({median_px:.2f} px)'
            elif np.isfinite(median_arc):
                fwhm_note = f'Median FWHM: {median_arc:.2f}"'
            elif np.isfinite(median_px):
                fwhm_note = f"Median FWHM: {median_px:.2f} px"
            else:
                fwhm_note = "Median FWHM: N/A"
            self.summary_label.setText(
                f"Detection Complete\n"
                f"{'─' * 30}\n"
                f"Files processed: {summary['total_files']}\n"
                f"Total sources: {summary['total_sources']}\n"
                f"Average per frame: {summary['avg_sources']:.1f}\n"
                f"{fwhm_note}"
            )
            self.log(f"Detection complete: {summary['total_files']} files, {summary['total_sources']} total sources")

            # Save state
            self.save_state()
        elif summary and resume_mode:
            self.update_summary_from_results(title="Detection Complete (cache+new)")
            self.log(f"Detection complete (cache+new): {len(self.detection_results)} files")
            self.save_state()
        else:
            self.summary_label.setText("Detection stopped or failed")
            if self.detection_results:
                self.save_state()

        if self.detection_results:
            self.update_navigation_buttons()
            self._refresh_qc_panel()
            # Auto-switch to QC tab after detection completes
            if not (summary and summary.get('stopped')):
                if hasattr(self, "tabs"):
                    self.tabs.setCurrentIndex(1)

        if not summary or not summary.get('stopped'):
            elapsed_txt = ""
            if hasattr(self, "_detect_start_time"):
                elapsed_txt = f" | elapsed {format_duration(time.monotonic() - self._detect_start_time)}"
            self.progress_label.setText(f"Done{elapsed_txt}")

    def closeEvent(self, event):
        """Ensure worker thread is stopped before closing window"""
        if self.detection_worker and self.detection_worker.isRunning():
            self.stop_detection()
            self.detection_worker.wait(5000)
        super().closeEvent(event)

    def populate_results_table(self):
        """Populate results table from detection_results"""
        self.results_table.setRowCount(0)
        for filename, result in self.detection_results.items():
            row = self.results_table.rowCount()
            self.results_table.insertRow(row)
            filename_item = QTableWidgetItem(filename)
            filename_item.setToolTip(filename)
            self.results_table.setItem(row, 0, filename_item)
            self.results_table.setItem(row, 1, QTableWidgetItem(str(result.get('n_sources', 0))))
            fwhm_arcsec = float(result.get('fwhm_arcsec', 0.0))
            fwhm_px = float(result.get('fwhm_px', 0.0))
            method = result.get('detect_method', 'segm')
            fwhm_str = f'{fwhm_arcsec:.2f}" ({fwhm_px:.2f} px; {method})'
            self.results_table.setItem(row, 2, QTableWidgetItem(fwhm_str))
            self.results_table.setItem(row, 3, QTableWidgetItem(f"{float(result.get('bkg_median', 0.0)):.1f}"))
            self.results_table.setItem(row, 4, QTableWidgetItem(result.get('filter', '')))
            self.results_table.setItem(row, 5, QTableWidgetItem(f"{float(result.get('sigma_used', 3.2)):.1f}"))

    def update_summary_from_results(self, title: str = "Detection Loaded"):
        """Update summary label from detection_results"""
        if not self.detection_results:
            self.summary_label.setText("No detection run yet")
            return
        total_sources = sum(r.get('n_sources', 0) for r in self.detection_results.values())
        avg_sources = total_sources / len(self.detection_results)
        fwhm_values = [r.get('fwhm_arcsec', 0.0) for r in self.detection_results.values() if r.get('fwhm_arcsec', 0.0) > 0]
        median_fwhm = float(np.median(fwhm_values)) if fwhm_values else 0.0
        self.summary_label.setText(
            f"{title}\n"
            f"{'─' * 30}\n"
            f"Files processed: {len(self.detection_results)}\n"
            f"Total sources: {total_sources}\n"
            f"Average per frame: {avg_sources:.1f}\n"
            f"Median FWHM: {median_fwhm:.2f}\""
        )


    def on_table_cell_clicked(self, row, col):
        """Handle table cell click - load that file"""
        filename = self.results_table.item(row, 0).text()
        idx = self.file_combo.findText(filename)
        if idx >= 0:
            self.file_combo.setCurrentIndex(idx)
            self.load_and_display()

    def log(self, message):
        """Add message to log"""
        append_timestamped_log(self.log_text, message)

    def clear_detection_cache(self):
        cache_dir = self.params.P.cache_dir
        if not cache_dir.exists():
            self.log("Cache directory not found.")
            return
        step4_out = step4_dir(self.params.P.result_dir)
        patterns = [
            "detect_*.csv",
            "detect_*.json",
            "detect_peak_*.csv",
            "detect_prev_peak_*.csv",
            "detect_tmp_peak_*.csv",
        ]
        removed = 0
        for pattern in patterns:
            for path in cache_dir.glob(pattern):
                try:
                    path.unlink()
                    removed += 1
                except Exception:
                    pass
            if step4_out.exists():
                for path in step4_out.glob(pattern):
                    try:
                        path.unlink()
                        removed += 1
                    except Exception:
                        pass
        self.detection_results = {}
        self.populate_results_table()
        self.update_summary_from_results()
        self.update_overlay()
        self.save_state()
        self.update_navigation_buttons()
        self._refresh_qc_panel()
        self.log(f"Detection cache cleared: {removed} files removed.")

    def open_parameters_dialog(self):
        """Open detection parameters dialog"""
        dialog = QDialog(self)
        configure_parameter_dialog(dialog, "Detection Parameters", 540, 680)

        outer_layout = QVBoxLayout(dialog)

        # Scrollable content area
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setContentsMargins(4, 4, 4, 4)
        scroll.setWidget(content)
        outer_layout.addWidget(scroll, 1)

        # Info
        info = QLabel("Adjust source detection parameters.\nChanges apply to next detection run.")
        info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 10px; margin-bottom: 10px; }")
        layout.addWidget(info)

        # Main form – always visible (mode, engine, base sigma)
        form = QFormLayout()

        self.param_detect_mode = QComboBox()
        self.param_detect_mode.addItem("Normal (기본)", "normal")
        self.param_detect_mode.addItem("Crowded (혼잡장)", "crowded")
        self.param_detect_mode.addItem("Faint (희미한 장)", "faint")
        self.param_detect_mode.addItem("Custom (수동)", "custom")
        current_mode = _get_detect_mode_from_params(self.params.P)
        for i in range(self.param_detect_mode.count()):
            if _normalize_detect_mode(self.param_detect_mode.itemData(i)) == current_mode:
                self.param_detect_mode.setCurrentIndex(i)
                break
        mode_row = QWidget()
        mode_row_layout = QHBoxLayout(mode_row)
        mode_row_layout.setContentsMargins(0, 0, 0, 0)
        mode_row_layout.setSpacing(8)
        mode_row_layout.addWidget(self.param_detect_mode, 1)
        self.param_detect_mode_apply = QPushButton("Apply Preset")
        mode_row_layout.addWidget(self.param_detect_mode_apply)
        form.addRow("Detection Mode:", mode_row)

        self.param_engine = QComboBox()
        self.param_engine.addItem("SEP (fast catalog)", "sep")
        self.param_engine.addItem("Photutils segmentation", "segm")
        self.param_engine.addItem("DAOStarFinder", "dao")
        self.param_engine.addItem("Peak finder", "peak")
        current_engine = _normalize_detect_engine(getattr(self.params.P, "detect_engine", "segm"))
        idx = self.param_engine.findData(current_engine)
        if idx >= 0:
            self.param_engine.setCurrentIndex(idx)
        form.addRow("Detection Engine:", self.param_engine)

        self.param_detect_sigma = QDoubleSpinBox()
        self.param_detect_sigma.setRange(1.0, 10.0)
        self.param_detect_sigma.setSingleStep(0.1)
        self.param_detect_sigma.setValue(float(getattr(self.params.P, 'detect_sigma', 3.2)))
        form.addRow("Detection Sigma (base):", self.param_detect_sigma)

        layout.addLayout(form)

        # === Per-Filter Sigma (collapsible) ===
        filter_group, filter_container = create_collapsible_section("Per-Filter Sigma (overrides base)")
        filter_layout = QGridLayout(filter_container)

        detected_filters = self.scan_filters_from_files()
        filter_info = (
            f"Detected filters from data: {', '.join(detected_filters)}"
            if detected_filters
            else "No filters detected (will use base sigma for all)"
        )
        filter_layout.addWidget(QLabel(filter_info), 0, 0, 1, 2)

        self.filter_sigma_edits = {}
        current_filters = set(detected_filters)
        current_filters.update(self.filter_sigma_map.keys())
        if not current_filters:
            current_filters = {'g', 'r', 'i'}

        row = 1
        for filt in sorted(current_filters):
            spin = QDoubleSpinBox()
            spin.setRange(1.0, 10.0)
            spin.setSingleStep(0.1)
            spin.setDecimals(2)
            spin.setValue(self.filter_sigma_map.get(filt, self.param_detect_sigma.value()))
            spin.setSpecialValueText("use base")
            filter_layout.addWidget(QLabel(f"{filt}:"), row, 0)
            filter_layout.addWidget(spin, row, 1)
            self.filter_sigma_edits[filt] = spin
            row += 1

        filter_layout.addWidget(QLabel("Add filter:"), row, 0)
        custom_layout = QHBoxLayout()
        self.custom_filter_name = QLineEdit()
        self.custom_filter_name.setPlaceholderText("e.g., u or z")
        self.custom_filter_name.setMaximumWidth(60)
        custom_layout.addWidget(self.custom_filter_name)
        self.custom_filter_sigma = QDoubleSpinBox()
        self.custom_filter_sigma.setRange(1.0, 10.0)
        self.custom_filter_sigma.setSingleStep(0.1)
        self.custom_filter_sigma.setValue(3.2)
        custom_layout.addWidget(self.custom_filter_sigma)
        btn_add = QPushButton("Add")
        btn_add.clicked.connect(lambda: self.add_custom_filter(filter_layout, row + 1))
        custom_layout.addWidget(btn_add)
        filter_layout.addLayout(custom_layout, row, 1)

        layout.addWidget(filter_group)

        # === Detection Options (collapsible) ===
        detect_opts_group, detect_opts_container = create_collapsible_section(
            "Detection Options",
            initial_expanded=True,
        )
        detect_opts_form = QFormLayout(detect_opts_container)
        detect_opts_form.setContentsMargins(0, 0, 0, 0)

        self.param_minarea = QSpinBox()
        self.param_minarea.setRange(1, 50)
        self.param_minarea.setValue(int(getattr(self.params.P, 'minarea_pix', 3)))
        detect_opts_form.addRow("Min Area (pixels):", self.param_minarea)

        self.param_fwhm_all_sources = QCheckBox("Enable")
        self.param_fwhm_all_sources.setChecked(
            bool(getattr(self.params.P, 'fwhm_measure_all_sources', False))
        )
        detect_opts_form.addRow("All-source radial FWHM:", self.param_fwhm_all_sources)

        self.param_deblend = QCheckBox("Enable")
        self.param_deblend.setChecked(getattr(self.params.P, 'deblend_enable', True))
        detect_opts_form.addRow("Deblending:", self.param_deblend)

        self.param_deblend_nthresh = QSpinBox()
        self.param_deblend_nthresh.setRange(8, 128)
        self.param_deblend_nthresh.setValue(int(getattr(self.params.P, 'deblend_nthresh', 64)))
        detect_opts_form.addRow("Deblend Levels:", self.param_deblend_nthresh)

        self.param_deblend_cont = QDoubleSpinBox()
        self.param_deblend_cont.setRange(0.001, 0.1)
        self.param_deblend_cont.setDecimals(4)
        self.param_deblend_cont.setSingleStep(0.001)
        self.param_deblend_cont.setValue(float(getattr(self.params.P, 'deblend_cont', 0.004)))
        detect_opts_form.addRow("Deblend Contrast:", self.param_deblend_cont)

        self.param_deblend_max_labels = QSpinBox()
        self.param_deblend_max_labels.setRange(500, 20000)
        self.param_deblend_max_labels.setValue(int(getattr(self.params.P, 'deblend_max_labels', 4000)))
        detect_opts_form.addRow("Deblend Soft Max Labels:", self.param_deblend_max_labels)

        self.param_deblend_label_hard_max = QSpinBox()
        self.param_deblend_label_hard_max.setRange(500, 50000)
        self.param_deblend_label_hard_max.setValue(int(getattr(self.params.P, 'deblend_label_hard_max', 7000)))
        detect_opts_form.addRow("Deblend Hard Max Labels:", self.param_deblend_label_hard_max)

        self.param_bkg2d = QCheckBox("Enable")
        self.param_bkg2d.setChecked(getattr(self.params.P, 'bkg2d_in_detect', True))
        detect_opts_form.addRow("2D Background:", self.param_bkg2d)

        self.param_bkg_box = QSpinBox()
        self.param_bkg_box.setRange(16, 256)
        self.param_bkg_box.setValue(int(getattr(self.params.P, 'bkg2d_box', 64)))
        detect_opts_form.addRow("Background Box:", self.param_bkg_box)

        layout.addWidget(detect_opts_group)

        # === DAO Refine (collapsible) ===
        dao_group, dao_container = create_collapsible_section("DAO Refine (hot pixel filter)")
        dao_layout = QFormLayout(dao_container)
        dao_layout.setContentsMargins(0, 0, 0, 0)

        self.param_dao_enable = QCheckBox("Enable")
        self.param_dao_enable.setChecked(getattr(self.params.P, 'dao_refine_enable', False))
        dao_layout.addRow("DAO refine:", self.param_dao_enable)

        self.param_dao_fwhm = QDoubleSpinBox()
        self.param_dao_fwhm.setRange(0.5, 20.0)
        self.param_dao_fwhm.setSingleStep(0.1)
        self.param_dao_fwhm.setValue(float(getattr(self.params.P, 'dao_fwhm_px', getattr(self.params.P, 'fwhm_seed_px', 6.0))))
        dao_layout.addRow("DAO FWHM (px):", self.param_dao_fwhm)

        self.param_dao_sharp_lo = QDoubleSpinBox()
        self.param_dao_sharp_lo.setRange(0.0, 2.0)
        self.param_dao_sharp_lo.setSingleStep(0.05)
        self.param_dao_sharp_lo.setValue(float(getattr(self.params.P, 'dao_sharp_lo', 0.2)))
        dao_layout.addRow("Sharpness min:", self.param_dao_sharp_lo)

        self.param_dao_sharp_hi = QDoubleSpinBox()
        self.param_dao_sharp_hi.setRange(0.0, 2.0)
        self.param_dao_sharp_hi.setSingleStep(0.05)
        self.param_dao_sharp_hi.setValue(float(getattr(self.params.P, 'dao_sharp_hi', 1.0)))
        dao_layout.addRow("Sharpness max:", self.param_dao_sharp_hi)

        self.param_dao_round_lo = QDoubleSpinBox()
        self.param_dao_round_lo.setRange(-2.0, 2.0)
        self.param_dao_round_lo.setSingleStep(0.05)
        self.param_dao_round_lo.setValue(float(getattr(self.params.P, 'dao_round_lo', -0.5)))
        dao_layout.addRow("Roundness min:", self.param_dao_round_lo)

        self.param_dao_round_hi = QDoubleSpinBox()
        self.param_dao_round_hi.setRange(-2.0, 2.0)
        self.param_dao_round_hi.setSingleStep(0.05)
        self.param_dao_round_hi.setValue(float(getattr(self.params.P, 'dao_round_hi', 0.5)))
        dao_layout.addRow("Roundness max:", self.param_dao_round_hi)

        self.param_dao_match_tol = QDoubleSpinBox()
        self.param_dao_match_tol.setRange(0.5, 10.0)
        self.param_dao_match_tol.setSingleStep(0.5)
        self.param_dao_match_tol.setValue(float(getattr(self.params.P, 'dao_match_tol_px', 2.0)))
        dao_layout.addRow("Match tolerance (px):", self.param_dao_match_tol)

        layout.addWidget(dao_group)

        # === Peak Assist (collapsible) ===
        peak_group, peak_container = create_collapsible_section("Peak Assist (segm supplement)")
        peak_layout = QFormLayout(peak_container)
        peak_layout.setContentsMargins(0, 0, 0, 0)
        peak_layout.setSpacing(4)

        self.param_peak_enable = QCheckBox("Enable")
        self.param_peak_enable.setChecked(getattr(self.params.P, 'peak_pass_enable', False))
        peak_layout.addRow("Peak assist:", self.param_peak_enable)

        self.param_peak_nsigma = QDoubleSpinBox()
        self.param_peak_nsigma.setRange(1.0, 10.0)
        self.param_peak_nsigma.setSingleStep(0.1)
        self.param_peak_nsigma.setValue(float(getattr(self.params.P, 'peak_nsigma', 3.2)))
        peak_layout.addRow("Peak n-sigma:", self.param_peak_nsigma)

        self.param_peak_scales = QLineEdit()
        self.param_peak_scales.setText(str(getattr(self.params.P, 'peak_kernel_scales', "0.9,1.3")))
        peak_layout.addRow("Kernel scales:", self.param_peak_scales)

        self.param_peak_min_sep = QDoubleSpinBox()
        self.param_peak_min_sep.setRange(0.5, 20.0)
        self.param_peak_min_sep.setSingleStep(0.5)
        self.param_peak_min_sep.setValue(float(getattr(self.params.P, 'peak_min_sep_px', 4.0)))
        peak_layout.addRow("Min separation (px):", self.param_peak_min_sep)

        self.param_peak_max_add = QSpinBox()
        self.param_peak_max_add.setRange(0, 10000)
        self.param_peak_max_add.setValue(int(getattr(self.params.P, 'peak_max_add', 600)))
        peak_layout.addRow("Max add:", self.param_peak_max_add)

        self.param_peak_max_elong = QDoubleSpinBox()
        self.param_peak_max_elong.setRange(1.0, 5.0)
        self.param_peak_max_elong.setSingleStep(0.1)
        self.param_peak_max_elong.setValue(float(getattr(self.params.P, 'peak_max_elong', 1.6)))
        peak_layout.addRow("Max elongation:", self.param_peak_max_elong)

        self.param_peak_sharp_lo = QDoubleSpinBox()
        self.param_peak_sharp_lo.setRange(0.0, 2.0)
        self.param_peak_sharp_lo.setSingleStep(0.05)
        self.param_peak_sharp_lo.setValue(float(getattr(self.params.P, 'peak_sharp_lo', 0.12)))
        peak_layout.addRow("Sharpness min:", self.param_peak_sharp_lo)

        self.param_peak_skip_if_nsrc = QSpinBox()
        self.param_peak_skip_if_nsrc.setRange(0, 20000)
        self.param_peak_skip_if_nsrc.setValue(int(getattr(self.params.P, 'peak_skip_if_nsrc_ge', 4500)))
        peak_layout.addRow("Skip if Nsrc >=:", self.param_peak_skip_if_nsrc)

        layout.addWidget(peak_group)
        layout.addStretch(1)

        # Connect detect mode signals
        self.param_detect_mode.currentIndexChanged.connect(self._update_detect_mode_ui_state)
        self.param_detect_mode_apply.clicked.connect(self._on_apply_detect_mode_clicked)
        self.param_engine.currentIndexChanged.connect(self._sync_engine_dependent_dialog_state)
        self._update_detect_mode_ui_state()
        self._sync_engine_dependent_dialog_state()

        # Buttons outside the scroll area
        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        buttons.accepted.connect(lambda: self.save_parameters(dialog))
        buttons.rejected.connect(dialog.reject)
        outer_layout.addWidget(buttons)

        dialog.exec_()

    def add_custom_filter(self, layout, row):
        """Add custom filter to the dialog"""
        name = self.custom_filter_name.text().strip()
        if not name:
            return

        sigma_val = self.custom_filter_sigma.value()

        if name not in self.filter_sigma_edits:
            lbl = QLabel(f"{name}:")
            spin = QDoubleSpinBox()
            spin.setRange(1.0, 10.0)
            spin.setSingleStep(0.1)
            spin.setDecimals(2)
            spin.setValue(sigma_val)
            layout.addWidget(lbl, row, 0)
            layout.addWidget(spin, row, 1)
            self.filter_sigma_edits[name] = spin

        self.custom_filter_name.clear()

    def save_parameters(self, dialog):
        """Save detection parameters"""
        selected_mode = self._selected_detect_mode()
        if selected_mode != "custom":
            self._apply_detect_mode_preset_to_dialog(selected_mode)
        self._sync_engine_dependent_dialog_state()
        self.params.P.detect_mode = selected_mode
        self.params.P.detect_engine = self._selected_detect_engine()
        self.params.P.detect_sigma = self.param_detect_sigma.value()
        self.params.P.minarea_pix = self.param_minarea.value()
        self.params.P.fwhm_measure_all_sources = self.param_fwhm_all_sources.isChecked()
        self.params.P.deblend_enable = self.param_deblend.isChecked()
        self.params.P.deblend_nthresh = self.param_deblend_nthresh.value()
        self.params.P.deblend_cont = self.param_deblend_cont.value()
        self.params.P.deblend_max_labels = self.param_deblend_max_labels.value()
        self.params.P.deblend_label_hard_max = self.param_deblend_label_hard_max.value()
        self.params.P.bkg2d_in_detect = self.param_bkg2d.isChecked()
        self.params.P.bkg2d_box = self.param_bkg_box.value()
        self.params.P.dao_refine_enable = self.param_dao_enable.isChecked()
        self.params.P.dao_fwhm_px = self.param_dao_fwhm.value()
        self.params.P.dao_sharp_lo = self.param_dao_sharp_lo.value()
        self.params.P.dao_sharp_hi = self.param_dao_sharp_hi.value()
        self.params.P.dao_round_lo = self.param_dao_round_lo.value()
        self.params.P.dao_round_hi = self.param_dao_round_hi.value()
        self.params.P.dao_match_tol_px = self.param_dao_match_tol.value()
        self.params.P.peak_pass_enable = self.param_peak_enable.isChecked()
        self.params.P.peak_nsigma = self.param_peak_nsigma.value()
        self.params.P.peak_kernel_scales = self.param_peak_scales.text().strip()
        self.params.P.peak_min_sep_px = self.param_peak_min_sep.value()
        self.params.P.peak_max_add = self.param_peak_max_add.value()
        self.params.P.peak_max_elong = self.param_peak_max_elong.value()
        self.params.P.peak_sharp_lo = self.param_peak_sharp_lo.value()
        self.params.P.peak_skip_if_nsrc_ge = self.param_peak_skip_if_nsrc.value()

        # Save filter-sigma mappings
        self.filter_sigma_map = {}
        for filt, spin in self.filter_sigma_edits.items():
            val = spin.value()
            self.filter_sigma_map[filt] = val
        self.params.P.detect_sigma_by_filter = dict(self.filter_sigma_map)
        self.params.P.detect_sigma_g = self.filter_sigma_map.get("g")
        self.params.P.detect_sigma_r = self.filter_sigma_map.get("r")
        self.params.P.detect_sigma_i = self.filter_sigma_map.get("i")

        saved = self.persist_params()
        self.save_state()

        msg = "Parameters saved!" if saved else "Parameters updated, but TOML save failed."
        QMessageBox.information(dialog, "Success", msg)
        dialog.accept()

    def validate_step(self) -> bool:
        """Validate if step can be completed"""
        return len(self.detection_results) > 0

    def save_state(self):
        """Save step state"""
        state_data = {
            "detection_complete": len(self.detection_results) > 0,
            "n_files": len(self.detection_results),
            "use_cropped": self.use_cropped,
            "filter_sigma_map": self.filter_sigma_map,
            "detect_mode": _get_detect_mode_from_params(self.params.P),
            "detect_engine": _normalize_detect_engine(getattr(self.params.P, "detect_engine", "segm")),
            "detect_sigma": self.params.P.detect_sigma,
            "minarea_pix": self.params.P.minarea_pix,
            "fwhm_measure_all_sources": getattr(self.params.P, "fwhm_measure_all_sources", False),
            "deblend_enable": self.params.P.deblend_enable,
            "deblend_nthresh": self.params.P.deblend_nthresh,
            "deblend_cont": self.params.P.deblend_cont,
            "deblend_max_labels": getattr(self.params.P, "deblend_max_labels", 4000),
            "deblend_label_hard_max": getattr(self.params.P, "deblend_label_hard_max", 7000),
            "bkg2d_in_detect": self.params.P.bkg2d_in_detect,
            "bkg2d_box": self.params.P.bkg2d_box,
            "dao_refine_enable": self.params.P.dao_refine_enable,
            "dao_fwhm_px": self.params.P.dao_fwhm_px,
            "dao_sharp_lo": self.params.P.dao_sharp_lo,
            "dao_sharp_hi": self.params.P.dao_sharp_hi,
            "dao_round_lo": self.params.P.dao_round_lo,
            "dao_round_hi": self.params.P.dao_round_hi,
            "dao_match_tol_px": self.params.P.dao_match_tol_px,
            "peak_pass_enable": getattr(self.params.P, "peak_pass_enable", True),
            "peak_nsigma": getattr(self.params.P, "peak_nsigma", 3.2),
            "peak_kernel_scales": getattr(self.params.P, "peak_kernel_scales", "0.9,1.3"),
            "peak_min_sep_px": getattr(self.params.P, "peak_min_sep_px", 4.0),
            "peak_max_add": getattr(self.params.P, "peak_max_add", 600),
            "peak_max_elong": getattr(self.params.P, "peak_max_elong", 1.6),
            "peak_sharp_lo": getattr(self.params.P, "peak_sharp_lo", 0.12),
            "peak_skip_if_nsrc_ge": getattr(self.params.P, "peak_skip_if_nsrc_ge", 4500),
        }
        if hasattr(self, "qc_panel") and self.qc_panel is not None:
            state_data["qc_state"] = self.qc_panel.export_state()
            self.qc_panel.write_frame_quality_csv()
        self.project_state.store_step_data("source_detection", state_data)

    def restore_state(self):
        """Restore step state"""
        state_data = self.project_state.get_step_data("source_detection")
        if state_data:
            # Restore filter sigma map
            if 'filter_sigma_map' in state_data:
                self.filter_sigma_map = state_data['filter_sigma_map']
            for key in [
                "detect_mode",
                "detect_engine",
                "detect_sigma",
                "minarea_pix",
                "fwhm_measure_all_sources",
                "deblend_enable",
                "deblend_nthresh",
                "deblend_cont",
                "deblend_max_labels",
                "deblend_label_hard_max",
                "bkg2d_in_detect",
                "bkg2d_box",
                "dao_refine_enable",
                "dao_fwhm_px",
                "dao_sharp_lo",
                "dao_sharp_hi",
                "dao_round_lo",
                "dao_round_hi",
                "dao_match_tol_px",
                "peak_pass_enable",
                "peak_nsigma",
                "peak_kernel_scales",
                "peak_min_sep_px",
                "peak_max_add",
                "peak_max_elong",
                "peak_sharp_lo",
                "peak_skip_if_nsrc_ge",
            ]:
                if key in state_data:
                    setattr(self.params.P, key, state_data[key])
            if "qc_state" in state_data and hasattr(self, "qc_panel") and self.qc_panel is not None:
                self.qc_panel.restore_state(state_data["qc_state"])
