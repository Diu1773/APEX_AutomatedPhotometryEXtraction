"""Qt-free source detection (headless port of DetectionWorker.run()).

This module is a VERBATIM RELOCATION of the compute body of
``apex.gui.workflow.step4_source_detection.DetectionWorker.run`` with only the
Qt coupling substituted by optional callbacks. The numerical algorithm, parameter
lookups, rounding, and file writes are unchanged.

LAYER RULE: this module must NOT import apex.gui or PyQt5.
"""

from __future__ import annotations

from pathlib import Path
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
from scipy.ndimage import gaussian_filter, median_filter
from scipy.spatial import cKDTree as KDTree

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
from apex.utils.astro_utils import compute_airmass_from_header, normalize_filter_name
from apex.utils.filter_parameters import (
    build_active_filter_float_map,
    normalize_filter_float_map,
)
from apex.utils.noise_params import resolve_effective_noise_params
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


def run_detection(file_list, params, data_dir, cache_dir, use_cropped=False,
                  filter_sigma_map=None, *, progress_cb=None, worker_status_cb=None,
                  file_done_cb=None, error_cb=None, should_stop=None, logger=None) -> dict:
    """Headless source detection over ``file_list``.

    Verbatim relocation of ``DetectionWorker.run``: same algorithm, parameter
    lookups, rounding, and on-disk outputs. Qt signals are replaced by optional
    callbacks; the summary dict is returned instead of emitted.
    """
    # Bind call/state shims so the relocated body keeps using the same names.
    data_dir = Path(data_dir)
    cache_dir = Path(cache_dir)
    result_dir = Path(getattr(params.P, "result_dir", data_dir))
    filter_sigma_map = normalize_filter_float_map(filter_sigma_map)

    def _should_stop():
        return should_stop() if should_stop else False

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

        print(f"[DetectionWorker] Starting with {len(file_list)} files")
        print(f"[DetectionWorker] Data dir: {data_dir}, use_cropped: {use_cropped}")

        results = {}
        total = len(file_list)
        P = params.P
        step4_out = step4_dir(result_dir)
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

        max_workers = get_parallel_workers(params)

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
            if _should_stop():
                return filename, None

            worker_id = get_worker_id()
            short_name = filename[:25] + "..." if len(filename) > 28 else filename

            try:
                # Stage 1: Loading
                (worker_status_cb(worker_id, short_name, "Loading", 10) if worker_status_cb else None)

                # Determine path
                if use_cropped:
                    cropped_dir = step2_cropped_dir(result_dir)
                    file_path = cropped_dir / filename
                else:
                    file_path = params.get_file_path(filename)
                source_sig = build_detection_cache_signature(
                    Path(file_path),
                    use_cropped=bool(use_cropped),
                    detect_engine=getattr(P, 'detect_engine', 'segm'),
                )

                # Load FITS. float32, not float64: halves per-frame RAM (244 MB
                # vs 488 MB for a 61 MP frame) and the SEP path downcasts to
                # float32 anyway, so this also drops the float64 intermediate.
                with fits.open(file_path) as hdul:
                    data = hdul[0].data.astype(np.float32)
                    header = hdul[0].header

                # Get filter - preserve original case from header
                filt = header.get('FILTER', '').strip()
                filt_key = normalize_filter_name(filt)

                # Canonical keys preserve Johnson R/I versus Sloan r/i.
                nsig = float(filter_sigma_map.get(filt_key, detect_sigma_base))

                # Stage 2: Background
                (worker_status_cb(worker_id, short_name, "Background", 25) if worker_status_cb else None)
                if _should_stop():
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
                        if _should_stop():
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
                        if _should_stop():
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
                if _should_stop():
                    return filename, None

                if detect_engine == "sep":
                    try:
                        import sep
                    except Exception as exc:
                        raise RuntimeError(
                            "SEP engine selected but the 'sep' package is not installed"
                        ) from exc

                    (worker_status_cb(worker_id, short_name, "SEP Background", 25) if worker_status_cb else None)
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

                    (worker_status_cb(worker_id, short_name, "SEP Detecting", 40) if worker_status_cb else None)
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
                        if _should_stop():
                            return filename, None
                        ds = max(1, int(getattr(P, 'bkg2d_downsample', 4)))
                        k = max(3, int(round(bkg2d_box / ds)))
                        small = data_filled[::ds, ::ds]
                        bkg_small = median_filter(small, size=k, mode="nearest")
                        if _should_stop():
                            return filename, None
                        bkg = np.repeat(np.repeat(bkg_small, ds, axis=0), ds, axis=1)[:data.shape[0], :data.shape[1]]
                        if _should_stop():
                            return filename, None
                        data_sub = data - bkg
                    else:
                        data_sub = data.copy()

                    work = np.where(np.isfinite(data_sub), data_sub, 0.0)
                    _ds = max(1, int(round(work.shape[0] ** 0.5)) // 8)
                    _, bkg_median, _ = sigma_clipped_stats(work[::_ds, ::_ds], sigma=3.0, maxiters=3)
                    if _should_stop():
                        return filename, None

                    # Smoothing (Jupyter-style)
                    sig = max(0.8, fwhm_seed / 2.355)
                    fim = gaussian_filter(work - bkg_median, sig, mode="nearest")
                    if _should_stop():
                        return filename, None

                    # Threshold (median + nsig * std)
                    _, mF, stdF = sigma_clipped_stats(fim[::_ds, ::_ds], sigma=3.0, maxiters=3)
                    bkg_rms = float(stdF)
                    threshold = mF + nsig * max(stdF, 1e-9)

                    # Stage 3: Detection
                    (worker_status_cb(worker_id, short_name, "Detecting", 40) if worker_status_cb else None)

                    # Segmentation detection
                    if detect_engine == "segm":
                        segm = detect_sources(
                            fim, threshold=threshold,
                            npixels=max(3, minarea_pix), connectivity=8
                        )
                    if _should_stop():
                        return filename, None

                if segm is not None and segm.nlabels > 0:
                    # Stage 4: Deblending
                    (worker_status_cb(worker_id, short_name, f"Deblend ({segm.nlabels})", 55) if worker_status_cb else None)

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
                            if _should_stop():
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
                            (worker_status_cb(worker_id, short_name, f"Deblend skip ({nlabels0})", 58) if worker_status_cb else None)

                    # Stage 5: Catalog
                    (worker_status_cb(worker_id, short_name, "Catalog", 75) if worker_status_cb else None)

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
                            (worker_status_cb(
                                worker_id, short_name,
                                f"FWHM={fwhm_prelim:.1f}px", 77
                            ) if worker_status_cb else None)
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
                    if _should_stop():
                        return filename, None
                    (worker_status_cb(worker_id, short_name, "DAO refine", 65) if worker_status_cb else None)

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

                        if _should_stop():
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
                        if _should_stop():
                            return filename, None
                        (worker_status_cb(worker_id, short_name, "Peak assist", 70) if worker_status_cb else None)
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
                            if _should_stop():
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
                (worker_status_cb(worker_id, short_name, "Saving", 90) if worker_status_cb else None)
                if _should_stop():
                    return filename, None

                # Save to cache
                cache_file = cache_dir / f"detect_{filename}.json"
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
                pos_file = cache_dir / f"detect_{filename}.csv"
                step4_pos = step4_out / f"detect_{filename}.csv"
                if positions:
                    source_records = []
                    peak_set = {(float(px), float(py)) for px, py in peak_positions}
                    measure_all_source_fwhm = bool(getattr(P, 'fwhm_measure_all_sources', False))
                    for idx, (x, y) in enumerate(positions):
                        if _should_stop():
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
                    peak_file = cache_dir / f"detect_peak_{filename}.csv"
                    step4_peak = step4_out / f"detect_peak_{filename}.csv"
                    np.savetxt(peak_file, peak_positions, delimiter=',',
                              header='x,y', comments='')
                    copy_if_different(peak_file, step4_peak)

                (worker_status_cb(worker_id, short_name, "Done", 100) if worker_status_cb else None)

                return filename, result

            except Exception as e:
                print(f"[DetectionWorker] Error in detect_single({filename}): {e}")
                return filename, {'error': str(e)}

        # Run detection in parallel
        print(f"[DetectionWorker] Starting ThreadPoolExecutor with {max_workers} workers")
        _executor = ThreadPoolExecutor(max_workers=max_workers)
        futures = {_executor.submit(detect_single, f): f for f in file_list}
        try:
            completed_count = 0

            for future in as_completed(futures):
                if _should_stop():
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
                            (error_cb(fname, result['error']) if error_cb else None)
                        else:
                            results[fname] = result
                            (file_done_cb(fname, result) if file_done_cb else None)
                except Exception as e:
                    (error_cb(filename, str(e)) if error_cb else None)

                (progress_cb(completed_count, total, filename, active) if progress_cb else None)
        finally:
            # Running frames can still write cache products. Join them
            # before emitting finished so a rerun cannot overlap writes.
            _executor.shutdown(wait=True, cancel_futures=True)
            _executor = None

        # Summary
        print(f"[DetectionWorker] Completed. {len(results)} results")
        if _should_stop():
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

        return summary

    except Exception as e:
        print(f"[DetectionWorker] FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
        (error_cb("WORKER", str(e)) if error_cb else None)
        return {}
