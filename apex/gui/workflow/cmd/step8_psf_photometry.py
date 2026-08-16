"""
Step 8 - PSF Photometry  (photutils EPSFBuilder + PSFPhotometry, skippable)

Reads:
    step4_detection/detect_{fname}.csv   (det_uid, x, y)
    step7_forced_phot/photometry_{fname}.tsv (initial flux estimate, optional)
    <FITS image>

Writes to cmd_psf/:
    photometry_{fname}.tsv   – det_uid, x_fit, y_fit, mag_psf, mag_psf_err,
                                chi2, iter_found, flags_psf
    epsf_model_{filter}_{frame_stem}.fits – oversampled ePSF model (per-frame)
    residual_iter{N}_{fname}.fits – sky-subtracted residual image (per iteration)
    starsub_iter{N}_{fname}.fits  – raw image with fitted stars removed (per iteration)
    photometry_index.csv     – per-frame summary
    psf_output_signature.json – selected frames, input mtimes, and PSF params

Step can be SKIPPED: clicking "Skip PSF" marks step as complete and
passes control to Step 9 (Master ID Editor). Downstream steps use Step 7 forced aperture
photometry results when PSF outputs are unavailable.
"""
from __future__ import annotations

import json
import hashlib
import traceback
import time
import copy
import threading
from dataclasses import replace
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from threading import Lock

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.nddata import NDData
from astropy.stats import sigma_clipped_stats, mad_std as _mad_std


from scipy.spatial import cKDTree

from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QGroupBox, QMessageBox,
    QTextEdit, QDialog, QFormLayout, QDialogButtonBox, QProgressBar,
    QCheckBox, QSpinBox, QDoubleSpinBox, QWidget, QTabWidget,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView,
    QSplitter, QComboBox, QListWidget, QScrollArea,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle, Patch

from apex.gui.layout_rules import FittedDialog, prevent_collapse, scroll_wrap, tame_canvas

from apex.analysis.cmd.psf_photometry_runner import (  # noqa: F401  (re-exported)
    PsfPhotometryRunner,
    _allstar_apply_model_inplace,
    _allstar_build_model,
    _allstar_fit,
    _allstar_newton_group,
    _allstar_newton_one,
    _build_groups,
    _build_psf_frame_qc_table,
    _build_psf_qc_summary,
    _fit_variance,
    _float32_difference,
    _load_detect_positions,
    build_ap_psf_comparison,
    build_psf_output_signature,
    export_psf_qc_products,
    load_psf_qc_inputs,
    render_psf_final_diagnostics,
    write_psf_output_signature,
)
from apex.gui.workflow.step_window_base import StepWindowBase
from apex.gui.workflow.log_panel import WorkflowLogWindow, WorkerStatusPanel, append_timestamped_log, show_raised
from apex.gui.workflow.ui_helpers import (
    add_parameter_reset_button,
    create_collapsible_section,
    create_output_reuse_checkbox,
    create_parameter_button,
    configure_parameter_dialog,
    set_table_row_background,
    status_row_background,
)
from apex.utils.step_paths_cmd import (
    step2_cropped_dir, step4_dir, step8_psf_dir,
    crop_is_active,
)
from apex.utils.step_paths import step7_forced_phot_dir
from apex.utils.astro_utils import normalize_filter_name
from apex.utils.constants import get_parallel_workers
from apex.utils.noise_params import resolve_effective_noise_params
from apex.utils.qc_utils import filter_files_by_qc, should_use_frame_quality_qc
from apex.utils.psf_core import (
    PSFCoreCut,
    estimate_psf_core_cut,
    psf_core_keep_mask,
    target_pixel_from_wcs,
)
from apex.analysis.psf_policy import (
    estimate_psf_flux_seeds,
    local_group_policy,
    merge_forced_catalog_seeds,
    plan_epsf_stars,
    plan_psf_fit_window,
    psf_symmetric_mask,
    select_epsf_reference_stars,
    select_spatially_balanced,
)
from apex.analysis.psf_iteration import (
    IterationSnapshot,
    PSFFitFlag,
    assess_psf_frame_quality,
    decide_residual_iteration,
    fit_parameters_changed,
    measure_psf_fit_quality,
    qfit_noise_diagnostics,
)
from apex.analysis.psf_diagnostics import (
    draw_psf_final_diagnostics,
    load_psf_final_diagnostic_data,
)
from apex.analysis.psf_flux_scale import (
    PSFApertureScale,
    apply_psf_aperture_scale,
    estimate_psf_aperture_scale,
)


# ── Scalar helpers ────────────────────────────────────────────────────────────











































# Korean → ASCII translation for matplotlib (matplotlib default font lacks Korean glyphs)
_KO_TO_ASCII = {
    "신규검출 (step4 미검출)": "New (not in step4)",
    "재검출 (step4 기검출)": "Re-detected (step4)",
    "경계소스": "Edge",
}














_PSF_MODE_PRESETS = {
    "normal": {
        "psf_n_stars_max": 0,
        "psf_isolation_fwhm_mult": 3.0,
        "psf_epsf_contamination_filter": True,
        "psf_flux_scale_correction": False,
        "psf_fit_shape_fwhm_mult": 2.4,
        "psf_fit_window_mode": "auto",
        "psf_fit_encircled_energy": 0.90,
        "psf_max_iter": 2,
        "psf_fitter_max_iter": 6,
        "psf_redetect_sigma": 4.0,
        "psf_duplicate_radius_fwhm_mult": 0.8,
        "psf_new_sources_cap_per_iter": 70,
        "psf_new_sources_cap_frac": 0.02,
        "psf_fit_init_max_sources": 0,
        "psf_core_cut_enable": False,
        "psf_core_cut_radius_px": 0.0,
        "psf_core_cut_radius_fwhm_mult": 20.0,
        "psf_core_cut_auto_min_density_ratio": 1.5,
        "psf_substar_neighbor_r_fwhm_mult": 8.0,
        "psf_substar_max_sources": 1500,
        "psf_substar_iters": 1,
        "psf_conv_new_frac": 0.02,
        "psf_postfit_snr_min": 3.0,
        "psf_postfit_qfit_max": 3.0,
        "psf_postfit_reduced_chi2_max": 25.0,
        "psf_blend_residual_ratio": 0.3,
        "psf_flux_conv_threshold": 0.01,
        "psf_use_grouper": False,
        "psf_grouper_radius_fwhm": 1.5,
        "psf_grouper_budget_frac": 0.10,
        "psf_grouper_budget_cap": 200,
        "psf_profile_error_frac": 0.0,
        "psf_final_pass_max_iter": 2,
        "psf_forced_position_lock": "always",
        "psf_forced_match_radius_fwhm": 1.25,
        "psf_redetect_sharp_lo": 0.15,
        "psf_redetect_sharp_hi": 0.95,
        "psf_redetect_round_abs_max": 0.8,
    },
    "crowded": {
        "psf_n_stars_max": 0,
        "psf_isolation_fwhm_mult": 2.0,
        "psf_epsf_contamination_filter": True,
        "psf_flux_scale_correction": False,
        "psf_fit_shape_fwhm_mult": 2.4,
        "psf_fit_window_mode": "auto",
        "psf_fit_encircled_energy": 0.90,
        "psf_max_iter": 2,
        "psf_fitter_max_iter": 8,
        "psf_redetect_sigma": 4.5,
        "psf_duplicate_radius_fwhm_mult": 0.4,
        "psf_new_sources_cap_per_iter": 50,
        "psf_new_sources_cap_frac": 0.015,
        "psf_fit_init_max_sources": 3000,
        "psf_core_cut_enable": False,
        "psf_core_cut_radius_px": 0.0,
        "psf_core_cut_radius_fwhm_mult": 20.0,
        "psf_core_cut_auto_min_density_ratio": 1.5,
        "psf_substar_neighbor_r_fwhm_mult": 5.0,
        "psf_substar_max_sources": 1000,
        "psf_substar_iters": 1,
        "psf_conv_new_frac": 0.02,
        "psf_postfit_snr_min": 3.0,
        "psf_postfit_qfit_max": 3.0,
        "psf_postfit_reduced_chi2_max": 25.0,
        "psf_blend_residual_ratio": 0.3,
        "psf_flux_conv_threshold": 0.01,
        "psf_use_grouper": False,
        "psf_grouper_radius_fwhm": 1.5,
        "psf_grouper_budget_frac": 0.10,
        "psf_grouper_budget_cap": 200,
        "psf_profile_error_frac": 0.0,
        "psf_final_pass_max_iter": 2,
        "psf_forced_position_lock": "always",
        "psf_forced_match_radius_fwhm": 1.25,
        "psf_redetect_sharp_lo": 0.2,
        "psf_redetect_sharp_hi": 0.9,
        "psf_redetect_round_abs_max": 0.6,
    },
    "faint": {
        "psf_n_stars_max": 0,
        "psf_isolation_fwhm_mult": 2.5,
        "psf_epsf_contamination_filter": True,
        "psf_flux_scale_correction": False,
        "psf_fit_shape_fwhm_mult": 2.4,
        "psf_fit_window_mode": "auto",
        "psf_fit_encircled_energy": 0.90,
        "psf_max_iter": 3,
        "psf_fitter_max_iter": 8,
        "psf_redetect_sigma": 3.0,
        "psf_duplicate_radius_fwhm_mult": 1.0,
        "psf_new_sources_cap_per_iter": 100,
        "psf_new_sources_cap_frac": 0.05,
        "psf_fit_init_max_sources": 0,
        "psf_core_cut_enable": False,
        "psf_core_cut_radius_px": 0.0,
        "psf_core_cut_radius_fwhm_mult": 20.0,
        "psf_core_cut_auto_min_density_ratio": 1.5,
        "psf_substar_neighbor_r_fwhm_mult": 8.0,
        "psf_substar_max_sources": 1500,
        "psf_substar_iters": 1,
        "psf_conv_new_frac": 0.03,
        "psf_postfit_snr_min": 3.0,
        "psf_postfit_qfit_max": 3.0,
        "psf_postfit_reduced_chi2_max": 25.0,
        "psf_blend_residual_ratio": 0.25,
        "psf_flux_conv_threshold": 0.01,
        "psf_use_grouper": False,
        "psf_grouper_radius_fwhm": 1.5,
        "psf_grouper_budget_frac": 0.10,
        "psf_grouper_budget_cap": 200,
        "psf_profile_error_frac": 0.0,
        "psf_final_pass_max_iter": 2,
        "psf_forced_position_lock": "always",
        "psf_forced_match_radius_fwhm": 1.25,
        "psf_redetect_sharp_lo": 0.1,
        "psf_redetect_sharp_hi": 0.95,
        "psf_redetect_round_abs_max": 0.9,
    },
}




# ── FITS helpers ───────────────────────────────────────────────────────────────





# ── Detect helpers ────────────────────────────────────────────────────────────





# ── Moffat PSF builder ────────────────────────────────────────────────────────







# ── Unified PSF evaluator (EPSF, Moffat, or the hybrid) ──────────────────────





# ── APEX iterative engine (ALLSTAR-inspired) ──────────────────────────────────

# Independent real-image injections in M13 and M3 show that single-source fits
# inside 1.5 FWHM can retain acceptable qfit/chi2 while absorbing neighbour flux.









def _allstar_fit_group(cleaned_patch: np.ndarray,
                       group_info: list,
                       patch_y0: int, patch_x0: int,
                       eval_psf, max_shift: float):
    """Simultaneously fit N stars on a pre-neighbour-subtracted patch.

    group_info : list of (x, y, flux) — absolute image coordinates.
    Returns    : list of (x_fit, y_fit, flux_fit, chi2, ok) per star.

    Parameters are subpixel offsets (dx, dy) + log-flux for numerical stability.
    Each dx/dy is relative to the integer-rounded star centre, keeping values near 0.
    log-flux prevents negative-flux blow-up during LM iterations.
    """
    from scipy.optimize import least_squares

    N = len(group_info)
    if N == 0:
        return []
    if N == 1:
        x0, y0, fl0 = group_info[0]
        return [_allstar_fit_one(cleaned_patch, x0, y0, patch_y0, patch_x0,
                                 fl0, eval_psf, max_shift)]

    ny, nx = cleaned_patch.shape
    if ny < 3 or nx < 3:
        return [(x, y, max(fl, 1.0), np.nan, False) for x, y, fl in group_info]

    yy_abs = np.arange(ny, dtype=float) + patch_y0
    xx_abs = np.arange(nx, dtype=float) + patch_x0
    YY, XX = np.meshgrid(yy_abs, xx_abs, indexing='ij')

    # Integer reference centres (keep dx/dy near zero for good LM conditioning)
    xi_refs = np.array([int(round(float(x))) for x, y, fl in group_info], dtype=int)
    yi_refs = np.array([int(round(float(y))) for x, y, fl in group_info], dtype=int)

    # p = [dx1, dy1, log_fl1, dx2, dy2, log_fl2, ...]
    # Using log-flux so flux stays positive and scale is comparable to dx/dy
    p0 = []
    fl_refs = []
    for n, (x, y, fl) in enumerate(group_info):
        fl_safe = max(float(fl), 1.0)
        fl_refs.append(fl_safe)
        p0.extend([float(x) - xi_refs[n], float(y) - yi_refs[n],
                   float(np.log(fl_safe))])
    p0 = np.array(p0, dtype=float)

    def _res(p):
        model = np.zeros((ny, nx), dtype=float)
        for n in range(N):
            dx, dy = p[3 * n], p[3 * n + 1]
            fl = float(np.exp(min(p[3 * n + 2], 30.0)))  # cap prevents overflow
            xc = xi_refs[n] + dx
            yc = yi_refs[n] + dy
            model += eval_psf(XX - xc, YY - yc) * fl
        diff = (cleaned_patch - model).ravel()
        return np.where(np.isfinite(diff), diff, 0.0)

    try:
        r = least_squares(_res, p0, method='lm',
                          ftol=1e-4, xtol=0.01, gtol=1e-6,
                          max_nfev=40 * N)
        chi2_base = float(np.mean(r.fun ** 2))
        out = []
        for n in range(N):
            dx_f, dy_f = r.x[3 * n], r.x[3 * n + 1]
            fl_f = float(np.exp(min(r.x[3 * n + 2], 30.0)))
            x0, y0, fl0 = group_info[n]
            xc_f = xi_refs[n] + dx_f
            yc_f = yi_refs[n] + dy_f
            if abs(dx_f) > max_shift or abs(dy_f) > max_shift or fl_f <= 0:
                out.append((x0, y0, max(fl0, 1.0), np.nan, False))
            else:
                out.append((xc_f, yc_f, fl_f,
                             chi2_base / max(fl_f ** 2 * 1e-8, 1e-20), True))
        return out
    except Exception:
        return [(x, y, max(fl, 1.0), np.nan, False) for x, y, fl in group_info]









def _allstar_fit_one(cleaned_patch: np.ndarray,
                     x0: float, y0: float,
                     patch_y0: int, patch_x0: int,
                     flux0: float, eval_psf,
                     max_shift: float = 2.0):
    """Fit one source on a pre-cleaned local patch.

    cleaned_patch : 2D array, already neighbour-subtracted, positioned at
                    image coords [patch_y0:patch_y0+ny, patch_x0:patch_x0+nx].
    Returns (x_fit, y_fit, flux_fit, chi2, ok).
    """
    from scipy.optimize import least_squares
    ny, nx = cleaned_patch.shape
    if ny < 3 or nx < 3:
        return x0, y0, flux0, np.nan, False

    xi, yi = int(round(x0)), int(round(y0))
    # Pixel offset grids (relative to integer star centre)
    yy = np.arange(ny, dtype=float) + patch_y0 - yi
    xx = np.arange(nx, dtype=float) + patch_x0 - xi
    YY, XX = np.meshgrid(yy, xx, indexing='ij')

    dx0, dy0 = x0 - xi, y0 - yi
    flux_safe = max(float(flux0) if np.isfinite(flux0) else 1.0, 1.0)

    # Cluster cores retain unresolved stellar light after global sky removal.
    # Fit a local constant term so that diffuse core light is not forced into
    # the target star's PSF flux.
    edge_mask = np.zeros(cleaned_patch.shape, dtype=bool)
    edge_mask[0, :] = True
    edge_mask[-1, :] = True
    edge_mask[:, 0] = True
    edge_mask[:, -1] = True
    edge_vals = cleaned_patch[edge_mask]
    edge_vals = edge_vals[np.isfinite(edge_vals)]
    if edge_vals.size:
        bg0 = float(np.nanmedian(edge_vals))
        bg_scale = float(_mad_std(edge_vals)) if edge_vals.size > 3 else float(np.nanstd(edge_vals))
    else:
        bg0 = 0.0
        bg_scale = 0.0
    if not np.isfinite(bg0):
        bg0 = 0.0
    if not np.isfinite(bg_scale) or bg_scale <= 0:
        bg_scale = max(1.0, abs(bg0) * 0.1)

    def _res(p):
        dx, dy, fl, bg = p
        diff = (cleaned_patch - (eval_psf(XX - dx, YY - dy) * fl + bg)).ravel()
        return np.where(np.isfinite(diff), diff, 0.0)

    try:
        bg_pad = max(5.0 * bg_scale, abs(bg0) + 10.0)
        r = least_squares(
            _res,
            [dx0, dy0, flux_safe, bg0],
            bounds=(
                [-max_shift, -max_shift, 1e-12, bg0 - bg_pad],
                [ max_shift,  max_shift, np.inf, bg0 + bg_pad],
            ),
            method='trf',
            ftol=1e-4,
            xtol=0.01,
            gtol=1e-6,
            max_nfev=80,
        )
        dx_f, dy_f, fl_f, _bg_f = r.x
        if abs(dx_f) > max_shift or abs(dy_f) > max_shift or fl_f <= 0:
            return x0, y0, flux_safe, np.nan, False
        chi2 = float(np.mean(r.fun ** 2)) / max(fl_f ** 2 * 1e-8, 1e-20)
        return xi + dx_f, yi + dy_f, fl_f, chi2, True
    except Exception:
        return x0, y0, flux_safe, np.nan, False




# ── PSF Worker ────────────────────────────────────────────────────────────────



# ── PSF Photometry Window ─────────────────────────────────────────────────────

# The calculation lives in `apex.analysis.cmd.psf_photometry_runner` so a script
# can run Step 8 without PyQt5. This adds the thread and the Qt signals, and
# subscribes each of the runner's channels to one — the window and the headless
# pipeline drive the very same object (2026-08-16).
# The runner comes FIRST in the bases — see the note on
# `ZeropointCalibrationWorker` in step10_zeropoint_calibration.py. In short:
# `QThread` first would shadow `run` with its own empty one, and sip's
# `QThread.__init__` would walk the cooperative chain into the runner's
# `__init__` with no arguments and raise.
class Step6PSFWorker(PsfPhotometryRunner, QThread):
    progress = pyqtSignal(int, int, str)
    worker_status = pyqtSignal(int, str, str, int)
    frame_done = pyqtSignal(str, dict)
    epsf_ready = pyqtSignal(str, str, object)
    residual_ready = pyqtSignal(str, object, object)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str, str)
    log = pyqtSignal(str)

    def __init__(self, file_list, params, data_dir, result_dir, cache_dir,
                 use_cropped=False):
        QThread.__init__(self)
        PsfPhotometryRunner.__init__(
            self, file_list, params, data_dir, result_dir, cache_dir, use_cropped)
        for name in self._CHANNELS:
            getattr(self, f"on_{name}").subscribe(getattr(self, name).emit)


class PSFPhotometryWindow(StepWindowBase):
    """Step 8 - PSF Photometry (skippable).

    If skipped, Step 9 Master ID Editor falls back to Step 7 forced photometry results.
    """

    def __init__(self, params, file_manager, project_state, main_window):
        self.file_manager = file_manager
        self.workflow_mode = str(getattr(main_window, "mode", "cmd")).lower()
        self.downstream_name = (
            "Target/Comparison Selection"
            if self.workflow_mode == "lc"
            else "Master ID Editor"
        )
        self.worker = None
        self.file_list = []
        self.use_cropped = False
        self.log_window = None
        self._skip_psf = False
        # In-memory cache (lost on window close → reloaded from disk in restore_state)
        self._last_epsf: dict[str, np.ndarray] = {}          # display_key → epsf array
        self._residual_meta: dict[str, dict] = {}            # fname  → residual metadata + iter records
        self._last_new_xy: dict[str, np.ndarray | None] = {} # fname  → new-detect XY or None
        self._cutout_idx: int = 0  # current star index in cutout viewer
        self._run_started_ts: float | None = None
        self._log_worker_frame: dict[int, str] = {}    # worker_id → current frame name
        self._current_psf_run_frames: list[str] = []
        self._current_psf_run_signature: dict | None = None
        self._psf_cache_validation_key: str | None = None
        self._psf_cache_validation_result: tuple[bool, str] = (False, "not checked")
        self._final_diag_data = pd.DataFrame()
        self._final_diag_summary: dict[str, object] = {}

        super().__init__(
            step_index=7,
            step_name="PSF Photometry",
            params=params,
            project_state=project_state,
            main_window=main_window,
        )
        self.setup_step_ui()
        self.restore_state()

    def setup_step_ui(self):
        info = QLabel(
            "Optional PSF photometry using photutils EPSFBuilder.\n"
            f"Click Skip PSF to continue to {self.downstream_name}; downstream "
            "steps will use Step 7 forced aperture photometry."
        )
        info.setWordWrap(True)
        info.setStyleSheet("QLabel { background-color: #E8F5E9; padding: 8px; margin-bottom: 6px; }")
        self.content_layout.addWidget(info)

        # ── Controls ──────────────────────────────────────────────────────────
        ctrl = QHBoxLayout()
        self.btn_params = create_parameter_button("PSF Parameters")
        self.btn_params.clicked.connect(self.open_parameters_dialog)
        ctrl.addWidget(self.btn_params)

        self.btn_skip = QPushButton("Skip PSF →")
        self.btn_skip.setStyleSheet(
            "QPushButton { background-color: #FF7043; color: white; font-weight: bold; padding: 8px 20px; }"
        )
        self.btn_skip.setToolTip(
            f"Skip PSF photometry; {self.downstream_name} will use Step 7 forced aperture results."
        )
        self.btn_skip.clicked.connect(self.skip_psf)
        ctrl.addWidget(self.btn_skip)

        self.chk_use_existing_output = create_output_reuse_checkbox(
            True,
            "Load Step 8 PSF outputs when photometry_index.csv, per-frame TSVs, "
            "ePSF/residual files, and the saved PSF signature all match the current run. "
            "Disable to force recomputation.",
        )
        ctrl.addWidget(self.chk_use_existing_output)

        ctrl.addStretch()

        self.btn_run = QPushButton("Run PSF")
        self.btn_run.setStyleSheet(
            "QPushButton { background-color: #388E3C; color: white; font-weight: bold; padding: 8px 24px; }"
        )
        self.btn_run.clicked.connect(self.run_psf)
        ctrl.addWidget(self.btn_run)

        self.btn_stop = QPushButton("Stop")
        self.btn_stop.setEnabled(False)
        self.btn_stop.clicked.connect(self.stop_psf)
        ctrl.addWidget(self.btn_stop)

        self.btn_log = QPushButton("Log")
        self.btn_log.clicked.connect(self.show_log_window)
        ctrl.addWidget(self.btn_log)

        self.content_layout.addLayout(ctrl)

        # ── Progress ──────────────────────────────────────────────────────────
        prog = QHBoxLayout()
        self.progress_bar = QProgressBar()
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        prog.addWidget(self.progress_bar, 1)
        self.progress_label = QLabel("Ready")
        self.progress_label.setMinimumWidth(420)
        prog.addWidget(self.progress_label)
        self.content_layout.addLayout(prog)


        # ── Skip status label ─────────────────────────────────────────────────
        self.skip_label = QLabel("")
        self.skip_label.setStyleSheet("QLabel { color: #FF7043; font-weight: bold; padding: 4px; }")
        self.content_layout.addWidget(self.skip_label)

        # ── Diagnostic tabs ───────────────────────────────────────────────────
        self.main_tabs = QTabWidget()
        self.content_layout.addWidget(self.main_tabs, 1)

        # Tab 0: EPSF Model
        epsf_tab = QWidget()
        epsf_layout = QVBoxLayout(epsf_tab)
        epsf_top = QHBoxLayout()
        epsf_top.addWidget(QLabel("Model:"))
        self.epsf_filter_combo = QComboBox()
        self.epsf_filter_combo.currentTextChanged.connect(self._on_epsf_filter_changed)
        epsf_top.addWidget(self.epsf_filter_combo)
        epsf_top.addStretch()
        epsf_layout.addLayout(epsf_top)

        self.epsf_fig = Figure(figsize=(8, 4))
        self.epsf_canvas = FigureCanvas(self.epsf_fig)
        self.epsf_toolbar = NavigationToolbar(self.epsf_canvas, self)
        epsf_layout.addWidget(self.epsf_toolbar)
        epsf_layout.addWidget(tame_canvas(self.epsf_canvas), 1)
        self.main_tabs.addTab(epsf_tab, "PSF Model")

        # Tab 1: Cutout viewer – Raw | Star-subtracted per star, per iter
        res_tab = QWidget()
        res_layout = QVBoxLayout(res_tab)

        # Row 1: frame / iter / mode selectors
        res_top = QHBoxLayout()
        res_top.addWidget(QLabel("Frame:"))
        self.res_file_combo = QComboBox()
        self.res_file_combo.currentTextChanged.connect(self._on_residual_frame_changed)
        res_top.addWidget(self.res_file_combo)
        res_top.addWidget(QLabel("Iter:"))
        self.res_iter_combo = QComboBox()
        self.res_iter_combo.currentTextChanged.connect(self._on_residual_iter_changed)
        res_top.addWidget(self.res_iter_combo)
        res_top.addStretch()
        res_layout.addLayout(res_top)

        # Row 2: star navigation (◀ / label / ▶)
        res_nav = QHBoxLayout()
        self.res_prev_btn = QPushButton("◀")
        self.res_prev_btn.setFixedWidth(36)
        self.res_prev_btn.clicked.connect(self._on_cutout_prev)
        res_nav.addWidget(self.res_prev_btn)
        self.res_star_label = QLabel("—")
        self.res_star_label.setMinimumWidth(70)
        self.res_star_label.setAlignment(Qt.AlignCenter)
        res_nav.addWidget(self.res_star_label)
        self.res_next_btn = QPushButton("▶")
        self.res_next_btn.setFixedWidth(36)
        self.res_next_btn.clicked.connect(self._on_cutout_next)
        res_nav.addWidget(self.res_next_btn)
        res_nav.addStretch()
        self.res_info_label = QLabel("")
        res_nav.addWidget(self.res_info_label)
        res_layout.addLayout(res_nav)

        self.res_fig = Figure(figsize=(8, 4))
        self.res_canvas = FigureCanvas(self.res_fig)
        self.res_toolbar = NavigationToolbar(self.res_canvas, self)
        res_layout.addWidget(self.res_toolbar)
        res_layout.addWidget(tame_canvas(self.res_canvas), 1)
        self.main_tabs.addTab(res_tab, "Residuals")

        # Tab 2: Photometry Table
        phot_tab = QWidget()
        phot_layout = QVBoxLayout(phot_tab)
        self.frame_table = QTableWidget()
        self.frame_table.setColumnCount(9)
        self.frame_table.setHorizontalHeaderLabels(
            ["Frame", "Filter", "N_psf", "N_goodmag", "N_fail", "N_new_iter",
             "Forced %", "PSF QC", "Time"]
        )
        # Forced % 는 PSF 결과를 읽을 때 반드시 함께 봐야 하는 값이다. 강제 측광
        # 위치(그 프레임에서 검출되지 않은 어두운 별)는 구경이 거의 0 을 재므로
        # PSF/구경 비교가 그쪽에서 폭주한다. 이 값이 높은 프레임의 PSF vs 구경
        # 통계는 검출된 별로만 걸러서 봐야 한다.
        _ft_tips = [
            "FITS 파일명", "필터명",
            "PSF 적합 대상 소스 수",
            "유효 등급을 얻은 수",
            "적합 실패 수",
            "잔차 재검출로 추가된 수",
            "강제 측광 비율 = n_forced / N_psf.\n"
            "그 프레임에서 검출되지 않아 마스터 위치로 강제 측광한 별의 비율.\n"
            "노출이 짧거나 청색 필터일수록 높다(실측: 같은 밤 B 75% vs R 50%).\n"
            "높으면 PSF/구경 플럭스 비교는 검출된 별로만 해야 한다.",
            "PSF QC 판정", "프레임 처리 시간",
        ]
        for _c, _tip in enumerate(_ft_tips):
            _it = self.frame_table.horizontalHeaderItem(_c)
            if _it is not None:
                _it.setToolTip(_tip)
        self.frame_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        for c in range(1, 9):
            self.frame_table.horizontalHeader().setSectionResizeMode(c, QHeaderView.ResizeToContents)
        self.frame_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        phot_layout.addWidget(self.frame_table)
        self.main_tabs.addTab(phot_tab, "Photometry")

        # Tab 3: QC Report (PSF statistics + Ap vs PSF comparison)
        qc_tab = QWidget()
        qc_outer = QVBoxLayout(qc_tab)

        qc_top_bar = QHBoxLayout()
        qc_refresh_btn = QPushButton("Refresh QC")
        qc_refresh_btn.clicked.connect(self._refresh_qc)
        qc_top_bar.addWidget(qc_refresh_btn)
        qc_top_bar.addStretch()
        qc_outer.addLayout(qc_top_bar)

        qc_splitter = QSplitter(Qt.Vertical)
        prevent_collapse(qc_splitter)
        qc_outer.addWidget(qc_splitter, 1)

        # Top: PSF statistics summary text
        self.qc_text = QTextEdit()
        self.qc_text.setReadOnly(True)
        self.qc_text.setStyleSheet("QTextEdit { font-family: monospace; font-size: 9pt; }")
        self.qc_text.setMaximumHeight(300)
        qc_splitter.addWidget(self.qc_text)

        # Bottom: Ap vs PSF matplotlib plot
        cmp_widget = QWidget()
        cmp_layout = QVBoxLayout(cmp_widget)
        cmp_top = QHBoxLayout()
        cmp_refresh_btn = QPushButton("Refresh Plot")
        cmp_refresh_btn.clicked.connect(self._plot_mag_comparison)
        cmp_top.addWidget(cmp_refresh_btn)
        cmp_top.addWidget(QLabel("Filter:"))
        self.cmp_filter_combo = QComboBox()
        self.cmp_filter_combo.addItem("all")
        self.cmp_filter_combo.currentTextChanged.connect(self._plot_mag_comparison)
        cmp_top.addWidget(self.cmp_filter_combo)
        cmp_top.addWidget(QLabel("Frame:"))
        self.cmp_frame_combo = QComboBox()
        self.cmp_frame_combo.addItem("all")
        self.cmp_frame_combo.currentTextChanged.connect(self._plot_mag_comparison)
        cmp_top.addWidget(self.cmp_frame_combo)
        self.cmp_flags0_only = QCheckBox("flags=0 only")
        self.cmp_flags0_only.setChecked(False)
        self.cmp_flags0_only.toggled.connect(self._plot_mag_comparison)
        cmp_top.addWidget(self.cmp_flags0_only)
        cmp_top.addWidget(QLabel("SNR ≥"))
        self.cmp_snr_min = QDoubleSpinBox()
        self.cmp_snr_min.setRange(0.0, 200.0)
        self.cmp_snr_min.setSingleStep(1.0)
        self.cmp_snr_min.setValue(0.0)
        self.cmp_snr_min.setDecimals(1)
        self.cmp_snr_min.setToolTip("0 = off")
        self.cmp_snr_min.valueChanged.connect(self._plot_mag_comparison)
        cmp_top.addWidget(self.cmp_snr_min)
        cmp_top.addWidget(QLabel("qfit/noise ≤"))
        self.cmp_qfit_max = QDoubleSpinBox()
        self.cmp_qfit_max.setRange(0.0, 10.0)
        self.cmp_qfit_max.setSingleStep(0.05)
        self.cmp_qfit_max.setValue(0.0)
        self.cmp_qfit_max.setDecimals(3)
        self.cmp_qfit_max.setToolTip("0 = off")
        self.cmp_qfit_max.valueChanged.connect(self._plot_mag_comparison)
        cmp_top.addWidget(self.cmp_qfit_max)
        cmp_top.addWidget(QLabel("|Δmag| ≤"))
        self.cmp_dmag_clip = QDoubleSpinBox()
        self.cmp_dmag_clip.setRange(0.0, 5.0)
        self.cmp_dmag_clip.setSingleStep(0.05)
        self.cmp_dmag_clip.setValue(0.0)
        self.cmp_dmag_clip.setDecimals(3)
        self.cmp_dmag_clip.setToolTip("0 = off")
        self.cmp_dmag_clip.valueChanged.connect(self._plot_mag_comparison)
        cmp_top.addWidget(self.cmp_dmag_clip)
        self.cmp_stats_label = QLabel("")
        self.cmp_stats_label.setWordWrap(True)
        cmp_top.addWidget(self.cmp_stats_label, 1)
        cmp_layout.addLayout(cmp_top)
        self.cmp_fig = Figure(figsize=(10, 4))
        self.cmp_canvas = FigureCanvas(self.cmp_fig)
        self.cmp_toolbar = NavigationToolbar(self.cmp_canvas, self)
        cmp_layout.addWidget(self.cmp_toolbar)
        cmp_layout.addWidget(tame_canvas(self.cmp_canvas), 1)
        qc_splitter.addWidget(cmp_widget)

        # Tallest page (478 px): scroll it so it does not set the window's
        # minimum height — see layout_rules.scroll_wrap.
        self.main_tabs.addTab(scroll_wrap(qc_tab), "QC")

        # Tab 4: per-frame final PSF diagnostics.
        final_tab = QWidget()
        final_layout = QVBoxLayout(final_tab)
        final_top = QHBoxLayout()
        final_top.addWidget(QLabel("Frame:"))
        self.final_diag_frame_combo = QComboBox()
        self.final_diag_frame_combo.setMinimumWidth(280)
        self.final_diag_frame_combo.currentTextChanged.connect(self._plot_final_diagnostics)
        final_top.addWidget(self.final_diag_frame_combo)
        final_refresh_btn = QPushButton("Refresh Diagnostics")
        final_refresh_btn.clicked.connect(self._refresh_final_diagnostics)
        final_top.addWidget(final_refresh_btn)
        self.final_diag_status = QLabel("Run Step 8 to generate diagnostics.")
        self.final_diag_status.setWordWrap(True)
        final_top.addWidget(self.final_diag_status, 1)
        final_layout.addLayout(final_top)

        self.final_diag_fig = Figure(figsize=(12, 7.2))
        self.final_diag_canvas = FigureCanvas(self.final_diag_fig)
        self.final_diag_toolbar = NavigationToolbar(self.final_diag_canvas, self)
        final_layout.addWidget(self.final_diag_toolbar)
        final_layout.addWidget(tame_canvas(self.final_diag_canvas), 1)
        self.main_tabs.addTab(final_tab, "Final Diagnostics")

        self.main_tabs.setCurrentIndex(0)

        # ── Log window ────────────────────────────────────────────────────────
        _log_worker_group = QGroupBox("Workers")
        _log_worker_group_layout = QVBoxLayout(_log_worker_group)
        _log_worker_group_layout.setContentsMargins(5, 5, 5, 5)
        self._worker_panel = WorkerStatusPanel(_log_worker_group)
        _log_worker_group_layout.addWidget(self._worker_panel)

        self.log_window = WorkflowLogWindow(
            self, "PSF Photometry Log & Workers",
            width=900, height=500,
            side_widget=_log_worker_group,
        )
        self.log_text = self.log_window.log_text

        # Keyboard shortcuts: ← → navigate cutout stars
        from PyQt5.QtWidgets import QShortcut
        from PyQt5.QtGui import QKeySequence
        QShortcut(QKeySequence(Qt.Key_Left),  self).activated.connect(self._on_cutout_prev)
        QShortcut(QKeySequence(Qt.Key_Right), self).activated.connect(self._on_cutout_next)

        self.populate_file_list()
        self.update_frame_table()
        self._update_skip_label()
        self._refresh_final_diagnostics()

    # ── File list ─────────────────────────────────────────────────────────────

    def populate_file_list(self):
        crop_active = crop_is_active(self.params.P.result_dir)
        cropped_dir = step2_cropped_dir(self.params.P.result_dir)
        if crop_active and cropped_dir.exists() and list(cropped_dir.glob("*.fit*")):
            files = sorted([f.name for f in cropped_dir.glob("*.fit*")])
            self.use_cropped = True
        else:
            if not self.file_manager.filenames:
                try:
                    self.file_manager.scan_files()
                except Exception:
                    pass
            files = self.file_manager.filenames
            self.use_cropped = False
        files = list(files)

        # Hard gate: downstream should skip frames where apcorr was not applied.
        apcorr_sum = step7_forced_phot_dir(self.params.P.result_dir) / "apcorr_summary.csv"
        if apcorr_sum.exists():
            try:
                df_apc = pd.read_csv(apcorr_sum)
                if (not df_apc.empty) and {"file", "apply"} <= set(df_apc.columns):
                    ok_vals = df_apc["apply"].astype(str).str.strip().str.lower().isin(
                        {"true", "1", "yes", "y", "on"}
                    )
                    ok_set = set(df_apc.loc[ok_vals, "file"].astype(str).map(lambda s: Path(str(s)).name))
                    before_n = len(files)
                    files = [f for f in files if Path(str(f)).name in ok_set]
                    self.log(f"[APCORR] apply=True frame filter: {len(files)}/{before_n} kept")
            except Exception as e:
                self.log(f"[APCORR] frame filter skipped ({e})")

        self.file_list = list(files)

    def _cache_dir_path(self) -> Path:
        return Path(getattr(self.params.P, "cache_dir", self.params.P.result_dir))

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
            return [PSFPhotometryWindow._signature_value(v) for v in value]
        if isinstance(value, dict):
            return {
                str(k): PSFPhotometryWindow._signature_value(v)
                for k, v in sorted(value.items(), key=lambda item: str(item[0]))
            }
        try:
            if pd.isna(value):
                return None
        except Exception:
            pass
        return str(value)

    @staticmethod
    def _first_existing_path(candidates: list[Path]) -> Path | None:
        for p in candidates:
            try:
                if p.exists() and p.stat().st_size > 0:
                    return p
            except Exception:
                continue
        return None

    @staticmethod
    def _newest_existing_path(candidates: list[Path]) -> Path | None:
        found = []
        for p in candidates:
            try:
                if p.exists() and p.stat().st_size > 0:
                    found.append(p)
            except Exception:
                continue
        if not found:
            return None
        return max(found, key=lambda p: p.stat().st_mtime_ns)

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

    def _psf_frames_after_qc(self) -> tuple[list[str], dict]:
        use_qc = should_use_frame_quality_qc(
            Path(self.params.P.result_dir),
            self.params.P,
            "phot_use_qc_pass_only",
            default=False,
        )
        frames, qc_info = filter_files_by_qc(
            Path(self.params.P.result_dir),
            list(self.file_list),
            require_qc=use_qc,
        )
        return list(frames), dict(qc_info or {})

    def _log_step4_qc_selection(self, qc_info: dict):
        if not should_use_frame_quality_qc(
            Path(self.params.P.result_dir),
            self.params.P,
            "phot_use_qc_pass_only",
            default=False,
        ):
            return
        if qc_info.get("applied"):
            self.log(f"Step4 QC: {qc_info.get('kept', 0)}/{qc_info.get('total', 0)} frame(s) kept.")
        elif qc_info.get("path") is None:
            self.log("Step4 QC: frame_quality.csv not found; using all frames.")
        else:
            self.log(f"Step4 QC: frame_quality.csv ignored ({qc_info.get('reason', 'unknown')}); using all frames.")

    def _build_psf_output_signature(self, frames: list[str]) -> dict:
        return build_psf_output_signature(
            self.params,
            frames,
            use_cropped=self.use_cropped,
            cache_dir=self._cache_dir_path(),
        )

    def _stored_psf_signature(self) -> dict | None:
        sig_path = step8_psf_dir(self.params.P.result_dir) / _PSF_SIGNATURE_FILE
        if not sig_path.exists():
            return None
        try:
            data = json.loads(sig_path.read_text(encoding="utf-8"))
            return data if isinstance(data, dict) else None
        except Exception:
            return None

    def _write_psf_output_signature(self, signature: dict):
        write_psf_output_signature(self.params.P.result_dir, signature)

    def _psf_signature_matches(self, signature: dict) -> tuple[bool, str]:
        stored = self._stored_psf_signature()
        if not stored:
            return False, "missing signature"
        if stored.get("signature_version") != _PSF_SIGNATURE_VERSION:
            return False, "signature version mismatch"
        if stored.get("signature_hash") != signature.get("signature_hash"):
            return False, "signature hash mismatch"
        return True, "ok"

    @staticmethod
    def _path_from_meta(out_dir: Path, name) -> Path | None:
        if name is None:
            return None
        text = str(name).strip()
        if not text:
            return None
        return out_dir / text

    def _existing_psf_output_covers(self, frames: list[str], signature: dict) -> tuple[bool, str]:
        ok, reason = self._psf_signature_matches(signature)
        if not ok:
            return False, reason

        out_dir = step8_psf_dir(self.params.P.result_dir)
        idx_path = out_dir / "photometry_index.csv"
        if not idx_path.exists():
            return False, "missing photometry_index.csv"
        try:
            idx = pd.read_csv(idx_path)
        except Exception as exc:
            return False, f"cannot read photometry_index.csv: {exc}"
        required_idx_cols = {"file", "filter"}
        if not required_idx_cols <= set(idx.columns):
            return False, "photometry_index.csv missing file/filter columns"

        expected_frames = [Path(str(f)).name for f in frames]
        expected_set = set(expected_frames)
        idx_files = [Path(str(f)).name for f in idx["file"].astype(str).tolist()]
        if len(idx_files) != len(expected_frames) or set(idx_files) != expected_set:
            return False, "photometry_index.csv frame set mismatch"

        expected_tsv_names = {f"photometry_{fname}.tsv" for fname in expected_frames}
        actual_tsv_names = {p.name for p in out_dir.glob("photometry_*.tsv")}
        if actual_tsv_names != expected_tsv_names:
            return False, "per-frame photometry TSV set mismatch"

        expected_epsf_paths: set[Path] = set()
        shared_epsf = bool(getattr(self.params.P, "psf_shared_filter_epsf", False))
        required_tsv_cols = {"det_uid", "x_fit", "y_fit", "mag_psf", "flags_psf"}
        if bool(getattr(self.params.P, "psf_flux_scale_correction", False)):
            required_tsv_cols.update({
                "flux_psf_raw_e",
                "psf_aperture_scale",
                "psf_aperture_scale_applied",
            })

        for fname in expected_frames:
            rows = idx[idx["file"].astype(str).map(lambda s: Path(s).name) == fname]
            if rows.empty:
                return False, f"missing index row for {fname}"
            filt = str(rows.iloc[0].get("filter", "")).strip()
            if not filt:
                return False, f"missing filter for {fname}"

            tsv_path = out_dir / f"photometry_{fname}.tsv"
            if not tsv_path.exists() or tsv_path.stat().st_size <= 0:
                return False, f"missing TSV for {fname}"
            try:
                tsv_head = pd.read_csv(tsv_path, sep="\t", nrows=5)
            except Exception as exc:
                return False, f"cannot read TSV for {fname}: {exc}"
            if not required_tsv_cols <= set(tsv_head.columns):
                return False, f"TSV columns incomplete for {fname}"

            for product_name in (f"residual_{fname}", f"starsub_{fname}"):
                product_path = out_dir / product_name
                if not product_path.exists() or product_path.stat().st_size <= 0:
                    return False, f"missing {product_name}"

            meta_path = out_dir / f"residual_meta_{fname}.json"
            if not meta_path.exists() or meta_path.stat().st_size <= 0:
                return False, f"missing residual metadata for {fname}"
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
            except Exception as exc:
                return False, f"cannot read residual metadata for {fname}: {exc}"
            if not isinstance(meta, dict):
                return False, f"invalid residual metadata for {fname}"
            epsf_reference = meta.get("epsf_reference", {})
            if (
                isinstance(epsf_reference, dict)
                and bool(epsf_reference.get("contamination_aware", False))
            ):
                reference_path = self._path_from_meta(
                    out_dir,
                    epsf_reference.get("catalog_path"),
                )
                if (
                    reference_path is None
                    or not reference_path.exists()
                    or reference_path.stat().st_size <= 0
                ):
                    return False, f"missing ePSF reference catalogue for {fname}"
            flux_scale = meta.get("flux_scale", {})
            if bool(getattr(self.params.P, "psf_flux_scale_correction", False)):
                if not isinstance(flux_scale, dict) or not bool(flux_scale.get("enabled", False)):
                    return False, f"missing PSF aperture-scale metadata for {fname}"
                scale_catalog = self._path_from_meta(out_dir, flux_scale.get("catalog_path"))
                if scale_catalog is not None and (
                    not scale_catalog.exists() or scale_catalog.stat().st_size <= 0
                ):
                    return False, f"missing PSF aperture-scale catalogue for {fname}"
            iters = meta.get("iters", [])
            if not isinstance(iters, list) or len(iters) == 0:
                return False, f"missing iteration metadata for {fname}"
            for key in ("seedxy_path", "rawxy_iter2_path"):
                p = self._path_from_meta(out_dir, meta.get(key))
                if p is None or not p.exists() or p.stat().st_size <= 0:
                    return False, f"missing {key} for {fname}"
            for rec in iters:
                if not isinstance(rec, dict):
                    return False, f"invalid iteration record for {fname}"
                for key in ("fitxy_path", "modelxy_path", "detxy_path", "appliedxy_path", "boxxy_path"):
                    p = self._path_from_meta(out_dir, rec.get(key))
                    if p is None or not p.exists() or p.stat().st_size <= 0:
                        return False, f"missing {key} for {fname}"
                for key in ("residual_path", "starsub_path"):
                    p = self._path_from_meta(out_dir, rec.get(key))
                    if p is not None and (not p.exists() or p.stat().st_size <= 0):
                        return False, f"missing {key} for {fname}"

            if shared_epsf:
                expected_epsf_paths.add(out_dir / f"epsf_model_{filt}.fits")
            else:
                expected_epsf_paths.add(out_dir / f"epsf_model_{filt}_{Path(fname).stem}.fits")

        for epsf_path in expected_epsf_paths:
            if not epsf_path.exists() or epsf_path.stat().st_size <= 0:
                return False, f"missing {epsf_path.name}"

        return True, "ok"

    def _current_psf_cache_status(self) -> tuple[bool, str]:
        frames, _ = self._psf_frames_after_qc()
        if not frames:
            return False, "no current frames"
        signature = self._build_psf_output_signature(frames)
        key = str(signature.get("signature_hash", ""))
        if key and key == self._psf_cache_validation_key:
            return self._psf_cache_validation_result
        result = self._existing_psf_output_covers(frames, signature)
        self._psf_cache_validation_key = key
        self._psf_cache_validation_result = result
        return result

    def _clear_psf_outputs(self) -> int:
        out_dir = step8_psf_dir(self.params.P.result_dir)
        if not out_dir.exists():
            return 0
        patterns = [
            _PSF_SIGNATURE_FILE,
            "photometry_index.csv",
            "photometry_*.tsv",
            "epsf_reference_*.csv",
            "psf_flux_scale_reference_*.csv",
            "epsf_model_*.fits",
            "residual_*",
            "starsub_*",
            "fitxy_iter*.npy",
            "modelxy_iter*.npy",
            "appliedxy_iter*.npy",
            "detxy_iter*.npy",
            "boxxy_iter*.npy",
            "seed_xy_*.npy",
            "rawxy_iter*.npy",
        ]
        removed = 0
        seen: set[Path] = set()
        for pat in patterns:
            for p in out_dir.glob(pat):
                if p in seen or not p.is_file():
                    continue
                seen.add(p)
                try:
                    p.unlink()
                    removed += 1
                except Exception:
                    pass
        return removed

    # ── Actions ───────────────────────────────────────────────────────────────

    def skip_psf(self):
        self._skip_psf = True
        self.save_state()
        self._update_skip_label()
        self.update_navigation_buttons()
        self.log(
            f"PSF skipped; {self.downstream_name} will use Step 7 forced aperture results."
        )

    def run_psf(self):
        if not self.file_list:
            QMessageBox.warning(self, "Warning", "No files to process")
            return
        if self.worker and self.worker.isRunning():
            return
        if not (step7_forced_phot_dir(self.params.P.result_dir) / "photometry_index.csv").exists():
            QMessageBox.warning(
                self, "Prerequisite",
                "Step 7 Forced Aperture Photometry must be completed first."
            )
            return

        frames_for_run, qc_info = self._psf_frames_after_qc()
        if not frames_for_run:
            QMessageBox.warning(
                self,
                "No frames",
                "No frames remain after Step 4 QC / Step 7 apcorr filtering.",
            )
            return
        signature = self._build_psf_output_signature(frames_for_run)
        self._psf_cache_validation_key = None
        self._current_psf_run_frames = list(frames_for_run)
        self._current_psf_run_signature = signature

        self._skip_psf = False
        self.log_text.clear()
        self.frame_table.setRowCount(0)
        self._last_epsf.clear()
        self._residual_meta.clear()
        self._last_new_xy.clear()
        self.epsf_filter_combo.clear()
        self.res_file_combo.clear()
        self.res_iter_combo.clear()
        # Keep frame selector populated during run so UI is not visually blank
        # before first residual metadata arrives.
        self.res_file_combo.addItems(frames_for_run)
        if self.res_file_combo.count() > 0:
            self.res_file_combo.setCurrentIndex(0)
        self._log_step4_qc_selection(qc_info)

        if getattr(self, "chk_use_existing_output", None) is not None and self.chk_use_existing_output.isChecked():
            cache_ok, cache_reason = self._existing_psf_output_covers(frames_for_run, signature)
            if cache_ok:
                self._psf_cache_validation_key = str(signature.get("signature_hash", ""))
                self._psf_cache_validation_result = (True, "ok")
                self.log(
                    f"[PSF][CACHE] Existing Step 8 output is complete "
                    f"({len(frames_for_run)} frame(s)); loading from disk."
                )
                self.progress_bar.setMaximum(len(frames_for_run))
                self.progress_bar.setValue(len(frames_for_run))
                self.progress_label.setText(
                    f"Cached Step 8 PSF output loaded ({len(frames_for_run)} frame(s))"
                )
                self._load_from_disk()
                self.update_frame_table()
                self._refresh_qc()
                self.save_state()
                self.update_navigation_buttons()
                self._current_psf_run_frames = []
                self._current_psf_run_signature = None
                return
            self.log(f"[PSF][CACHE] Existing output not reusable: {cache_reason}")

        removed = self._clear_psf_outputs()
        if removed:
            self.log(f"[PSF][CACHE] Removed {removed} stale Step 8 output file(s) before recompute.")

        # Clear log window worker bars from previous run
        self._log_worker_frame.clear()
        if hasattr(self, "_worker_panel") and self._worker_panel is not None:
            self._worker_panel.clear()

        self.log(f"Start PSF photometry | {len(frames_for_run)} frames")
        self._run_started_ts = time.time()

        self.worker = Step6PSFWorker(
            frames_for_run, self.params,
            self.params.P.data_dir, self.params.P.result_dir,
            self.params.P.cache_dir, self.use_cropped,
        )
        self.worker.progress.connect(self.on_progress)
        self.worker.worker_status.connect(self.on_worker_status)
        self.worker.frame_done.connect(self.on_frame_done)
        self.worker.epsf_ready.connect(self.on_epsf_ready)
        self.worker.residual_ready.connect(self.on_residual_ready)
        self.worker.finished.connect(self.on_finished)
        self.worker.error.connect(self.on_error)
        self.worker.log.connect(self.log)

        self.btn_run.setEnabled(False)
        self.btn_stop.setEnabled(True)
        self.btn_skip.setEnabled(False)
        self.progress_bar.setMaximum(len(frames_for_run))
        self.progress_bar.setValue(0)
        self.progress_label.setText(f"0/{len(frames_for_run)} | ETA --:-- | Starting...")
        self.worker.start()
        self.show_log_window()

    def stop_psf(self):
        if self.worker and self.worker.isRunning():
            self.worker.stop()
            self.btn_stop.setEnabled(False)
            self.progress_label.setText("Stopping... (running frames will finish current fit)")

    # ── Callbacks ─────────────────────────────────────────────────────────────

    def on_worker_status(self, worker_id: int, frame: str, stage: str, pct: int):
        self._log_worker_frame[int(worker_id)] = frame
        if hasattr(self, "_worker_panel") and self._worker_panel is not None:
            self._worker_panel.update_worker(worker_id, frame, stage, pct)

    def on_progress(self, current, total, filename):
        self.progress_bar.setValue(current)
        eta_txt = "pending"
        if self._run_started_ts is not None and current > 0 and total > 0:
            elapsed = max(0.0, time.time() - float(self._run_started_ts))
            per_frame = elapsed / float(current)
            eta_sec = max(0.0, per_frame * float(max(0, total - current)))
            eta_txt = self._fmt_duration(eta_sec)
        self.progress_label.setText(f"{current}/{total} | ETA {eta_txt} | {filename}")

    def on_frame_done(self, filename, result):
        r = self.frame_table.rowCount()
        self.frame_table.insertRow(r)
        self.frame_table.setItem(r, 0, QTableWidgetItem(filename))
        self.frame_table.setItem(r, 1, QTableWidgetItem(str(result.get("filter", ""))))
        self.frame_table.setItem(r, 2, QTableWidgetItem(str(result.get("n", 0))))
        self.frame_table.setItem(r, 3, QTableWidgetItem(str(result.get("n_goodmag", 0))))
        self.frame_table.setItem(r, 4, QTableWidgetItem(str(result.get("n_fail", 0))))
        self.frame_table.setItem(r, 5, QTableWidgetItem(str(result.get("n_new_iter", 0))))
        qc_status = str(result.get("psf_qc_status", "") or "")
        qc_item = QTableWidgetItem(qc_status)
        qc_item.setToolTip(str(result.get("psf_qc_reasons", "") or ""))
        self.frame_table.setItem(r, 6, qc_item)
        elapsed = _safe_float(result.get("frame_total_elapsed_s", np.nan), np.nan)
        self.frame_table.setItem(
            r, 7, QTableWidgetItem(f"{elapsed:.1f} s" if np.isfinite(elapsed) else "")
        )
        try:
            has_good_phot = int(result.get("n_goodmag", 0) or 0) > 0
        except (TypeError, ValueError):
            has_good_phot = False
        if qc_status == "FAIL":
            row_background = status_row_background(False)
        elif qc_status == "REVIEW":
            row_background = status_row_background(True, warning=True)
        else:
            row_background = status_row_background(has_good_phot)
        set_table_row_background(self.frame_table, r, row_background)
        self.frame_table.scrollToBottom()
        # Mark log window worker bar as done
        for w_key, fname in self._log_worker_frame.items():
            if fname == filename:
                if hasattr(self, "_worker_panel") and self._worker_panel is not None:
                    self._worker_panel.update_worker(w_key, fname, "Done", 100)
                break

    def on_epsf_ready(self, display_key: str, _frame_name: str, epsf_arr: np.ndarray):
        self._last_epsf[display_key] = epsf_arr
        current = self.epsf_filter_combo.currentText()
        self.epsf_filter_combo.blockSignals(True)
        self.epsf_filter_combo.clear()
        self.epsf_filter_combo.addItems(sorted(self._last_epsf.keys()))
        if current in self._last_epsf:
            self.epsf_filter_combo.setCurrentText(current)
        else:
            try:
                self.epsf_filter_combo.setCurrentText(display_key)
            except Exception:
                self.epsf_filter_combo.setCurrentIndex(0)
        self.epsf_filter_combo.blockSignals(False)
        self._plot_epsf(display_key)

    def on_residual_ready(self, fname: str, residual_meta: dict, new_xy):
        if isinstance(residual_meta, dict):
            self._residual_meta[fname] = residual_meta
        self._last_new_xy[fname] = new_xy  # ndarray or None
        current = self.res_file_combo.currentText()
        self.res_file_combo.blockSignals(True)
        self.res_file_combo.clear()
        self.res_file_combo.addItems(sorted(self._residual_meta.keys()))
        if current in self._residual_meta:
            self.res_file_combo.setCurrentText(current)
        else:
            self.res_file_combo.setCurrentIndex(self.res_file_combo.count() - 1)
        self.res_file_combo.blockSignals(False)
        self._cutout_idx = 0
        self._refresh_residual_iter_combo(fname)
        self._plot_cutout(fname)

    def on_error(self, src, err):
        self.log(f"ERROR {src}: {err}")

    def on_finished(self, summary):
        self.btn_run.setEnabled(True)
        self.btn_stop.setEnabled(False)
        self.btn_skip.setEnabled(True)
        elapsed_txt = ""
        if self._run_started_ts is not None:
            elapsed_txt = f" | elapsed {self._fmt_duration(time.time() - float(self._run_started_ts))}"
        self.progress_label.setText(f"Done{elapsed_txt}")
        self._run_started_ts = None
        self.log(f"PSF done: {summary}")
        if isinstance(summary, dict) and self._current_psf_run_signature:
            expected = len(self._current_psf_run_frames)
            processed = _to_int(summary.get("processed", 0), 0)
            stopped = _to_int(summary.get("stopped", 0), 0)
            if expected > 0 and stopped == 0 and processed == expected:
                self._write_psf_output_signature(self._current_psf_run_signature)
                cache_ok, cache_reason = self._existing_psf_output_covers(
                    self._current_psf_run_frames,
                    self._current_psf_run_signature,
                )
                if cache_ok:
                    self.log("[PSF][CACHE] Output signature saved for future reuse.")
                else:
                    sig_path = step8_psf_dir(self.params.P.result_dir) / _PSF_SIGNATURE_FILE
                    try:
                        sig_path.unlink()
                    except Exception:
                        pass
                    self.log(f"[PSF][CACHE] Signature not saved: output incomplete ({cache_reason}).")
            else:
                self.log(
                    "[PSF][CACHE] Output reuse disabled for this run: "
                    f"processed={processed}/{expected}, stopped={stopped}."
                )
        self._cleanup_worker()
        self._update_skip_label()
        self.update_frame_table()  # refresh Photometry tab from disk
        self._refresh_qc()           # refresh QC tab (stats + Ap vs PSF plot)
        self.save_state()
        self._psf_cache_validation_key = None
        self.update_navigation_buttons()
        self._current_psf_run_frames = []
        self._current_psf_run_signature = None

    # ── EPSF plot ─────────────────────────────────────────────────────────────

    def _on_epsf_filter_changed(self, display_key: str):
        self._plot_epsf(display_key)

    @staticmethod
    def _fmt_duration(seconds: float) -> str:
        s = int(max(0, round(float(seconds))))
        h, rem = divmod(s, 3600)
        m, sec = divmod(rem, 60)
        if h > 0:
            return f"{h:d}:{m:02d}:{sec:02d}"
        return f"{m:02d}:{sec:02d}"

    def _plot_epsf(self, display_key: str):
        if display_key not in self._last_epsf:
            return
        epsf_arr = self._last_epsf[display_key]

        is_moffat = display_key.startswith("[MOFFAT]")
        is_epsf   = display_key.startswith("[EPSF]")
        psf_label = "Moffat PSF" if is_moffat else "ePSF"
        px_label  = "px (native)" if is_moffat else "px (oversampled)"

        self.epsf_fig.clf()
        ax2d  = self.epsf_fig.add_subplot(121)
        ax_rad = self.epsf_fig.add_subplot(122)

        vmax = np.nanpercentile(epsf_arr, 99)
        im = ax2d.imshow(epsf_arr, origin="lower", cmap="viridis",
                         norm=mcolors.PowerNorm(gamma=0.5, vmin=0, vmax=max(vmax, 1e-10)))
        self.epsf_fig.colorbar(im, ax=ax2d, fraction=0.046, pad=0.04)
        ax2d.set_title(f"{psf_label} — {display_key}", fontsize=9)
        ax2d.set_xlabel(px_label, fontsize=8)
        ax2d.set_ylabel(px_label, fontsize=8)

        cy, cx = np.array(epsf_arr.shape) / 2.0
        yy, xx = np.mgrid[0:epsf_arr.shape[0], 0:epsf_arr.shape[1]]
        rr = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
        r_flat, v_flat = rr.ravel(), epsf_arr.ravel()
        order = np.argsort(r_flat)
        ax_rad.plot(r_flat[order], v_flat[order], ".", markersize=1, alpha=0.3, color="#1565C0")
        ax_rad.set_xlabel(f"Radius ({px_label})", fontsize=8)
        ax_rad.set_ylabel("PSF value", fontsize=8)
        ax_rad.set_title("Radial profile", fontsize=9)
        ax_rad.grid(True, alpha=0.3)

        self.epsf_fig.tight_layout()
        self.epsf_canvas.draw_idle()

    # ── Residual plot ─────────────────────────────────────────────────────────

    def _on_residual_frame_changed(self, fname: str):
        self._cutout_idx = 0
        self._refresh_residual_iter_combo(fname)
        self._plot_cutout(fname)

    def _on_residual_iter_changed(self, _iter_label: str):
        self._cutout_idx = 0
        fname = self.res_file_combo.currentText().strip()
        if fname:
            self._plot_cutout(fname)

    def _on_cutout_prev(self):
        if self._cutout_idx > 0:
            self._cutout_idx -= 1
            fname = self.res_file_combo.currentText().strip()
            if fname:
                self._plot_cutout(fname)

    def _on_cutout_next(self):
        self._cutout_idx += 1
        fname = self.res_file_combo.currentText().strip()
        if fname:
            self._plot_cutout(fname)

    def _get_iter_records(self, fname: str) -> list[dict]:
        meta = self._residual_meta.get(fname, {})
        recs = meta.get("iters", []) if isinstance(meta, dict) else []
        if not isinstance(recs, list):
            return []
        out = [r for r in recs if isinstance(r, dict)]
        out.sort(key=lambda d: int(d.get("iter", 0)))
        return out

    def _get_iter_record_by_no(self, fname: str, iter_no: int) -> dict | None:
        for rec in self._get_iter_records(fname):
            if int(rec.get("iter", -1)) == int(iter_no):
                return rec
        return None

    def _refresh_residual_iter_combo(self, fname: str):
        recs = self._get_iter_records(fname)
        self.res_iter_combo.blockSignals(True)
        self.res_iter_combo.clear()
        for rec in recs:
            i = int(rec.get("iter", 0))
            phase = str(rec.get("phase", "residual_fit"))
            label = f"{i} final flux" if phase == "final_flux" else f"{i} residual"
            self.res_iter_combo.addItem(label, i)
        if self.res_iter_combo.count() > 0:
            self.res_iter_combo.setCurrentIndex(self.res_iter_combo.count() - 1)
        self.res_iter_combo.blockSignals(False)

    def _load_starsub_for_iter(self, fname: str, rec: dict) -> np.ndarray | None:
        out_dir = step8_psf_dir(self.params.P.result_dir)
        meta = self._residual_meta.get(fname, {})
        bkg_med = float(meta.get("bkg_med", 0.0)) if isinstance(meta, dict) else 0.0

        starsub_name = str(rec.get("starsub_path", "")).strip()
        if starsub_name:
            p = out_dir / starsub_name
            if p.exists():
                try:
                    return fits.getdata(str(p)).astype(float)
                except Exception:
                    pass

        residual_name = str(rec.get("residual_path", "")).strip()
        if residual_name:
            p = out_dir / residual_name
            if p.exists():
                try:
                    res = fits.getdata(str(p)).astype(float)
                    return res + bkg_med
                except Exception:
                    pass

        return None

    def _load_snapshot_image(self, rec: dict, key: str) -> np.ndarray | None:
        out_dir = step8_psf_dir(self.params.P.result_dir)
        image_name = str(rec.get(key, "")).strip()
        if not image_name:
            return None
        path = out_dir / image_name
        if not path.exists():
            return None
        try:
            return fits.getdata(str(path)).astype(float)
        except Exception:
            return None

    def _load_xy_npy_for_iter(self, rec: dict, key: str, max_points: int = 500) -> np.ndarray:
        out_dir = step8_psf_dir(self.params.P.result_dir)
        arr_name = str(rec.get(key, "")).strip()
        if not arr_name:
            return np.zeros((0, 2), dtype=float)
        p = out_dir / arr_name
        if not p.exists():
            return np.zeros((0, 2), dtype=float)
        try:
            arr = np.load(str(p), allow_pickle=False)
            arr = np.asarray(arr, dtype=float)
            if arr.ndim != 2 or arr.shape[1] < 2:
                return np.zeros((0, 2), dtype=float)
            arr = arr[:, :2]
            finite = np.isfinite(arr[:, 0]) & np.isfinite(arr[:, 1])
            arr = arr[finite]
            if int(max_points) > 0:
                return arr[:max(0, int(max_points))]
            return arr
        except Exception:
            return np.zeros((0, 2), dtype=float)

    def _load_boxxy_for_iter(self, rec: dict, max_boxes: int = 500) -> np.ndarray:
        # Preferred: delta boxes (applied-from-previous + detected-this-iter).
        arr = self._load_xy_npy_for_iter(rec, "boxxy_path", max_points=max_boxes)
        if len(arr):
            return arr
        # Fallback 1: compose from separate arrays if present.
        arr_applied = self._load_xy_npy_for_iter(rec, "appliedxy_path", max_points=0)
        arr_detected = self._load_xy_npy_for_iter(rec, "detxy_path", max_points=0)
        if len(arr_applied) and len(arr_detected):
            arr = np.vstack([arr_applied, arr_detected])
        elif len(arr_applied):
            arr = arr_applied
        elif len(arr_detected):
            arr = arr_detected
        else:
            # Fallback 2: all fitted stars in this iteration.
            arr = self._load_xy_npy_for_iter(rec, "fitxy_path", max_points=0)
        if int(max_boxes) > 0:
            return arr[:max(0, int(max_boxes))]
        return arr

    def _load_modelxy_for_iter(self, rec: dict, max_points: int = 500) -> np.ndarray:
        arr = self._load_xy_npy_for_iter(rec, "modelxy_path", max_points=max_points)
        if len(arr):
            return arr
        # Backward compatibility with runs before modelxy_path existed:
        # iter>=2 should map to "new in this iter" rather than cumulative fit list.
        iter_no = int(rec.get("iter", 1))
        if iter_no > 1:
            arr = self._load_xy_npy_for_iter(rec, "detxy_path", max_points=max_points)
            if len(arr):
                return arr
        return self._load_xy_npy_for_iter(rec, "fitxy_path", max_points=max_points)

    def _resolve_fits_path_window(self, fname: str) -> Path | None:
        if self.use_cropped and crop_is_active(self.params.P.result_dir):
            cdir = step2_cropped_dir(self.params.P.result_dir)
            cpath = cdir / fname
            if cpath.exists():
                return cpath
        fpath = Path(self.params.P.data_dir) / fname
        return fpath if fpath.exists() else None

    # ── Cutout viewer ─────────────────────────────────────────────────────────

    def _plot_cutout(self, fname: str):  # noqa: C901
        """Show Raw | Star-subtracted cutouts for the selected star."""
        if fname not in self._residual_meta:
            self.res_fig.clf()
            ax = self.res_fig.add_subplot(111)
            ax.text(
                0.5, 0.5,
                "No residual result yet for this frame.\n(Still processing or frame skipped/failed)",
                transform=ax.transAxes,
                ha="center", va="center", fontsize=10, color="gray",
            )
            ax.set_title(fname, fontsize=9)
            self.res_star_label.setText("0/0")
            self.res_info_label.setText("waiting for residual_meta...")
            self.res_canvas.draw_idle()
            return
        recs = self._get_iter_records(fname)
        if not recs:
            return

        iter_no = self.res_iter_combo.currentData()
        selected = recs[-1]
        if iter_no is not None:
            for rec in recs:
                if int(rec.get("iter", -1)) == int(iter_no):
                    selected = rec
                    break

        # Fixed semantics requested by user:
        # - iter1: detected/fitted stars in iter1
        # - iter>=2: stars detected from residual(iter-1)
        iter_val = int(selected.get("iter", 0))
        phase = str(selected.get("phase", "residual_fit"))
        det_xy = self._load_xy_npy_for_iter(selected, "detxy_path", max_points=0)
        model_xy = self._load_modelxy_for_iter(selected, max_points=0)
        if phase == "final_flux":
            xy_list = self._load_xy_npy_for_iter(selected, "fitxy_path", max_points=0)
            mode_label = "fixed-position final flux"
        elif iter_val <= 1:
            xy_list = model_xy
            mode_label = "iter1 fitted stars"
        else:
            xy_list = det_xy if len(det_xy) > 0 else model_xy
            mode_label = f"iter{iter_val} detected-from-residual"

        res_std = float(selected.get("residual_std", np.nan))
        n_new_raw = int(selected.get("n_new_raw", 0))
        n_new_kept = int(selected.get("n_new_kept", 0))
        n_candidate_raw = int(selected.get("n_candidates_raw", 0))
        n_candidate_accepted = int(selected.get("n_candidates_accepted", 0))
        median_qfit = float(selected.get("median_qfit", np.nan))
        median_redchi = float(selected.get("median_reduced_chi2", np.nan))
        stop_reason = str(selected.get("stop_reason", ""))

        # Cutout half-size: driven by epsf_size_px so the PSF footprint is visible
        epsf_sz = int(selected.get("epsf_size_px", 25))
        half = max(epsf_sz // 2 + 4, 10)

        # Load raw FITS
        raw_img = None
        try:
            p = self._resolve_fits_path_window(fname)
            if p is not None:
                raw_img = fits.getdata(str(p)).astype(float)
        except Exception:
            pass

        model_img_snapshot = self._load_snapshot_image(selected, "model_path")
        residual_img_snapshot = self._load_snapshot_image(selected, "residual_path")

        full_cut_sz = 2 * half + 1

        def _cut_at(img, x_val: float, y_val: float):
            if img is None:
                return None
            nr, nc = img.shape
            cx = int(round(float(x_val)))
            cy = int(round(float(y_val)))
            x0 = max(0, cx - half)
            x1 = min(nc, cx + half + 1)
            y0 = max(0, cy - half)
            y1 = min(nr, cy + half + 1)
            cut = img[y0:y1, x0:x1]
            if cut.size == 0:
                return None
            if cut.shape == (full_cut_sz, full_cut_sz):
                return cut
            # Edge source: pad with NaN so the star stays centred in the panel.
            # NaN renders as the colormap bad-colour (neutral), making the
            # padding region visually distinct from real background.
            padded = np.full((full_cut_sz, full_cut_sz), np.nan, dtype=np.float64)
            dst_y = max(0, half - cy)
            dst_x = max(0, half - cx)
            padded[dst_y:dst_y + (y1 - y0), dst_x:dst_x + (x1 - x0)] = cut
            return padded

        # ── Filter xy_list to photometry-successful sources only ─────────────
        # Load the Step 8 PSF photometry TSV and keep only positions where
        # mag_psf is finite and FLAG_SAT is not set.  Saturated / edge /
        # fit-fail sources have NaN mag_psf or FLAG_SAT=1 — showing their
        # cutouts is misleading (over-subtraction rings, clipped PSF, etc.).
        if len(xy_list) > 0:
            try:
                _psf_tsv = step8_psf_dir(self.params.P.result_dir) / f"photometry_{fname}.tsv"
                if _psf_tsv.exists():
                    _df_phot = pd.read_csv(_psf_tsv, sep="\t")
                    if "saturated_psf" in _df_phot.columns:
                        _not_saturated = ~_df_phot["saturated_psf"].astype(str).str.lower().isin(
                            {"true", "1", "yes"}
                        )
                    else:
                        # Compatibility with catalogs written before standard fit flags.
                        _not_saturated = (
                            pd.to_numeric(
                                _df_phot.get("flags_psf", pd.Series(0, index=_df_phot.index)),
                                errors="coerce",
                            ).fillna(0).astype(int) & 1
                        ) == 0
                    _good = (
                        pd.to_numeric(_df_phot.get("mag_psf", pd.Series(dtype=float)), errors="coerce").notna() &
                        _not_saturated
                    )
                    _good_xy = _df_phot.loc[_good, ["x_fit", "y_fit"]].to_numpy(dtype=float)
                    if len(_good_xy) > 0:
                        _tree_good = cKDTree(_good_xy)
                        _d, _ = _tree_good.query(xy_list, k=1, workers=1)
                        _match_r = 1.5  # px
                        _keep = np.asarray(_d, dtype=float) <= _match_r
                        if np.any(_keep):
                            xy_list = xy_list[_keep]
            except Exception:
                pass  # on any error fall through and show all sources

        # Move edge sources to the end so idx=0 always shows a well-centred star.
        if raw_img is not None and len(xy_list) > 1:
            _nr, _nc = raw_img.shape
            _not_edge = (
                (xy_list[:, 0] >= half) & (xy_list[:, 0] < _nc - half) &
                (xy_list[:, 1] >= half) & (xy_list[:, 1] < _nr - half)
            )
            if np.any(_not_edge) and not np.all(_not_edge):
                _order = np.concatenate([np.where(_not_edge)[0], np.where(~_not_edge)[0]])
                xy_list = xy_list[_order]

        n = int(len(xy_list))
        idx = max(0, min(self._cutout_idx, n - 1)) if n > 0 else 0
        self._cutout_idx = idx
        self.res_star_label.setText(f"{idx + 1}/{n}" if n > 0 else "0/0")

        self.res_fig.clf()
        if n == 0:
            self.res_info_label.setText(
                f"iter {iter_val} | stars={n} | "
                f"new(raw/used)={n_new_raw}/{n_new_kept} | res_std={res_std:.3f}"
            )
            ax = self.res_fig.add_subplot(111)
            ax.text(
                0.5, 0.5,
                f"No sources for {mode_label}.",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=10, color="gray",
            )
            ax.set_title(f"{fname}  iter {iter_val}", fontsize=9)
            self.res_canvas.draw_idle()
            return

        x_c, y_c = float(xy_list[idx, 0]), float(xy_list[idx, 1])

        # ── Classify iter2+ detection: "신규검출" vs "재검출" ──────────────────
        # Compare the star's fitted position against step4 seed positions.
        # Within 2px → the seed existed in step4 (iter1 fit it but subtraction was poor).
        # Beyond 2px → genuinely new source not in step4.
        seed_tag = ""
        if iter_val >= 2:
            out_dir_ui = step8_psf_dir(self.params.P.result_dir)
            meta_ui = self._residual_meta.get(fname, {})
            seedxy_name = meta_ui.get("seedxy_path", "") if isinstance(meta_ui, dict) else ""
            if seedxy_name:
                seed_path_ui = out_dir_ui / seedxy_name
                if seed_path_ui.exists():
                    try:
                        seed_xy_ui = np.load(str(seed_path_ui)).astype(float)
                        if len(seed_xy_ui) > 0:
                            d_seed = np.hypot(seed_xy_ui[:, 0] - x_c, seed_xy_ui[:, 1] - y_c)
                            if np.min(d_seed) <= 2.0:
                                seed_tag = "재검출 (step4 기검출)"
                            else:
                                seed_tag = "신규검출 (step4 미검출)"
                    except Exception:
                        pass

        # Edge tag: shown when the star is within `half` pixels of the image boundary.
        edge_tag = ""
        if raw_img is not None:
            _nr, _nc = raw_img.shape
            if not (half <= x_c < _nc - half and half <= y_c < _nr - half):
                edge_tag = "경계소스"

        tags = "  " + "  ".join(f"[{t}]" for t in [seed_tag, edge_tag] if t) if (seed_tag or edge_tag) else ""
        self.res_info_label.setText(
            f"pass {iter_val} {phase} | stars={n} | new fit={n_new_kept} | "
            f"candidates={n_candidate_raw}/{n_candidate_accepted} | "
            f"res_std={res_std:.3f} qfit50={median_qfit:.3f} "
            f"chi2r50={median_redchi:.2f} | stop={stop_reason or '-'} | "
            f"xy=({x_c:.2f},{y_c:.2f}){tags}"
        )

        def _cut(img):
            return _cut_at(img, x_c, y_c)

        cut_raw = _cut(raw_img)
        cut_model = _cut(model_img_snapshot)
        cut_residual = _cut(residual_img_snapshot)
        panels: list[dict] = [
            {"img": cut_raw, "title": "Raw", "mark_detect": False, "residual": False},
            {"img": cut_model, "title": "PSF model", "mark_detect": False, "residual": False},
            {
                "img": cut_residual,
                "title": "Sky-sub residual",
                "mark_detect": phase != "final_flux" and iter_val > 1,
                "residual": True,
            },
        ]

        for i, p in enumerate(panels):
            cut = p.get("img", None)
            title = str(p.get("title", ""))
            mark_detect = bool(p.get("mark_detect", False))
            ax = self.res_fig.add_subplot(1, len(panels), i + 1)
            if cut is not None and cut.size > 0:
                if bool(p.get("residual", False)):
                    vmax = float(np.nanpercentile(np.abs(cut), 99))
                    panel_vmin, panel_vmax = -max(vmax, 1e-10), max(vmax, 1e-10)
                    panel_cmap = "coolwarm"
                else:
                    panel_vmin, panel_vmax = np.nanpercentile(cut, [1, 99])
                    panel_cmap = "gray"
                im = ax.imshow(
                    cut,
                    origin="lower",
                    cmap=panel_cmap,
                    vmin=panel_vmin,
                    vmax=panel_vmax,
                    interpolation="nearest",
                )
                self.res_fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                # Crosshair at star center
                cy_c = (cut.shape[0] - 1) / 2.0
                cx_c = (cut.shape[1] - 1) / 2.0
                ax.axhline(cy_c, color="#FF4444", lw=0.8, alpha=0.55, ls="--")
                ax.axvline(cx_c, color="#FF4444", lw=0.8, alpha=0.55, ls="--")
                if mark_detect:
                    ax.plot(
                        [cx_c], [cy_c],
                        marker="o", markersize=8,
                        markerfacecolor="none",
                        markeredgecolor="#FF5555",
                        markeredgewidth=1.2,
                    )
            else:
                ax.text(0.5, 0.5, "N/A", transform=ax.transAxes,
                        ha="center", va="center", fontsize=10, color="gray")
            ax.set_title(title, fontsize=8)
            ax.set_xlabel("Δx (px)", fontsize=7)
            ax.set_ylabel("Δy (px)", fontsize=7)

        suptitle_tags = "  " + "  ".join(f"[{_KO_TO_ASCII.get(t, t)}]" for t in [seed_tag, edge_tag] if t) if (seed_tag or edge_tag) else ""
        self.res_fig.suptitle(
            f"{fname}  |  {mode_label} #{idx + 1}/{n}{suptitle_tags}",
            fontsize=8, y=1.01,
        )
        self.res_fig.tight_layout()
        self.res_canvas.draw_idle()

    # ── Frame table refresh ───────────────────────────────────────────────────

    def update_frame_table(self):
        idx_path = step8_psf_dir(self.params.P.result_dir) / "photometry_index.csv"
        if not idx_path.exists() or not hasattr(self, "frame_table"):
            return
        try:
            idx = pd.read_csv(idx_path)
        except Exception:
            return
        self.frame_table.setRowCount(len(idx))
        for r, row in enumerate(idx.itertuples(index=False)):
            self.frame_table.setItem(r, 0, QTableWidgetItem(str(getattr(row, "file", ""))))
            self.frame_table.setItem(r, 1, QTableWidgetItem(str(getattr(row, "filter", ""))))
            self.frame_table.setItem(r, 2, QTableWidgetItem(str(int(_safe_float(getattr(row, "n", 0), 0)))))
            self.frame_table.setItem(r, 3, QTableWidgetItem(str(int(_safe_float(getattr(row, "n_goodmag", 0), 0)))))
            self.frame_table.setItem(r, 4, QTableWidgetItem(str(int(_safe_float(getattr(row, "n_fail", 0), 0)))))
            self.frame_table.setItem(r, 5, QTableWidgetItem(str(int(_safe_float(getattr(row, "n_new_iter", 0), 0)))))
            n_psf = _safe_float(getattr(row, "n", 0), 0)
            n_forced = _safe_float(getattr(row, "n_forced", np.nan), np.nan)
            if np.isfinite(n_forced) and n_psf > 0:
                frac = 100.0 * n_forced / n_psf
                fitem = QTableWidgetItem(f"{frac:.0f}%")
                fitem.setToolTip(
                    f"{int(n_forced)} / {int(n_psf)} — 그 프레임에서 검출되지 않아 "
                    "마스터 위치로 강제 측광한 별.\n높으면 PSF vs 구경 비교를 "
                    "검출된 별로만 해야 한다(강제 위치는 구경이 거의 0)."
                )
            else:
                fitem = QTableWidgetItem("")
            self.frame_table.setItem(r, 6, fitem)
            qc_status = str(getattr(row, "psf_qc_status", "") or "")
            qc_item = QTableWidgetItem(qc_status)
            qc_item.setToolTip(str(getattr(row, "psf_qc_reasons", "") or ""))
            self.frame_table.setItem(r, 7, qc_item)
            elapsed = _safe_float(getattr(row, "frame_total_elapsed_s", np.nan), np.nan)
            self.frame_table.setItem(
                r, 8, QTableWidgetItem(f"{elapsed:.1f} s" if np.isfinite(elapsed) else "")
            )
            if qc_status == "FAIL":
                set_table_row_background(self.frame_table, r, status_row_background(False))
            elif qc_status == "REVIEW":
                set_table_row_background(
                    self.frame_table, r, status_row_background(True, warning=True)
                )
            else:
                has_good_phot = int(
                    _safe_float(getattr(row, "n_goodmag", 0), 0)
                ) > 0
                set_table_row_background(
                    self.frame_table, r, status_row_background(has_good_phot)
                )

    # ── QC Report ─────────────────────────────────────────────────────────────

    def _refresh_qc(self):
        """Compute PSF QC statistics and update the QC tab text + Ap vs PSF plot."""
        if not hasattr(self, "qc_text"):
            return
        psf_dir = step8_psf_dir(self.params.P.result_dir)
        idx_path = psf_dir / "photometry_index.csv"
        export_inputs = None
        if not idx_path.exists():
            self.qc_text.setPlainText("photometry_index.csv not found.\nRun Step 8 first.")
            self._cmp_merged_df = None
            self._plot_mag_comparison()
            self._refresh_final_diagnostics()
            return
        try:
            # 헤드리스 러너와 **같은 함수**로 읽는다 — 창과 배치가 서로 다른
            # 코드로 QC 를 만들면 둘의 산출물이 조용히 갈라진다.
            idx, all_df, meta_df = load_psf_qc_inputs(psf_dir)
            if all_df.empty:
                self.qc_text.setPlainText("No photometry TSV files found.")
                self._refresh_final_diagnostics()
                return

            good = all_df[all_df["flags_psf"] == 0].copy() if "flags_psf" in all_df.columns else all_df.copy()
            filters = sorted(all_df["FILTER"].dropna().unique().tolist()) if "FILTER" in all_df.columns else []
            n_total = len(all_df)
            n_clean = len(good)

            W = 60
            lines = []
            lines.append("─" * W)
            lines.append("  PSF Photometry QC Report")
            lines.append("─" * W)
            filt_counts = "  ".join(
                f"{f}:{(idx['filter'] == f).sum()}" for f in filters if "filter" in idx.columns
            )
            lines.append(f"  총 프레임    : {len(idx)}  ({filt_counts})")
            lines.append(f"  총 검출 소스 : {n_total:,}")
            lines.append(f"  flags=0      : {n_clean:,} / {n_total:,} = {100 * n_clean / max(n_total, 1):.1f}%")
            lines.append("")

            # 1. Filter stats
            lines.append("  1. 필터별 검출 통계")
            lines.append(f"  {'필터':^4}  {'프레임':^6}  {'평균검출':^8}  {'성공률':^7}  {'실패율':^6}  {'iter2추가':^9}")
            lines.append("  " + "─" * 52)
            for filt in filters:
                si = idx[idx["filter"] == filt] if "filter" in idx.columns else pd.DataFrame()
                if si.empty:
                    continue
                avg_n = si["n"].mean() if "n" in si.columns else 0
                avg_g = si["n_goodmag"].mean() if "n_goodmag" in si.columns else avg_n
                ok_pct = 100 * avg_g / max(avg_n, 1)
                avg_new = si["n_new_iter"].mean() if "n_new_iter" in si.columns else 0
                lines.append(
                    f"  {filt:^4}  {len(si):^6}  {avg_n:^8.0f}  {ok_pct:^6.1f}%  "
                    f"{100 - ok_pct:^5.1f}%  avg {avg_new:.1f}"
                )
            lines.append("")

            # 2. Mag range & error
            lines.append("  2. 등급 범위 (mag_psf, flags=0)")
            lines.append(f"  {'필터':^4}  {'범위':^15}  {'평균':^6}  {'중앙값':^6}  {'σ':^5}  {'err중앙값':^9}")
            lines.append("  " + "─" * 58)
            for filt in filters:
                sub = good[good["FILTER"] == filt] if "FILTER" in good.columns else pd.DataFrame()
                mag = sub["mag_psf"].dropna() if "mag_psf" in sub.columns else pd.Series()
                err = sub["mag_psf_err"].dropna() if "mag_psf_err" in sub.columns else pd.Series()
                if mag.empty:
                    continue
                lines.append(
                    f"  {filt:^4}  {mag.min():.2f} ~ {mag.max():.2f}  "
                    f"{mag.mean():^6.2f}  {mag.median():^6.2f}  {mag.std():^5.2f}  "
                    f"{err.median():^9.4f}" if not err.empty else
                    f"  {filt:^4}  {mag.min():.2f} ~ {mag.max():.2f}  "
                    f"{mag.mean():^6.2f}  {mag.median():^6.2f}  {mag.std():^5.2f}  {'N/A':^9}"
                )
            lines.append("")

            # 3. SNR
            if "snr_psf" in good.columns:
                lines.append("  3. SNR 분포 (flags=0)")
                lines.append(f"  {'필터':^4}  {'10%':^8}  {'median':^8}  {'90%':^8}")
                lines.append("  " + "─" * 36)
                for filt in filters:
                    sub = good[good["FILTER"] == filt]["snr_psf"].dropna() if "FILTER" in good.columns else pd.Series()
                    if sub.empty:
                        continue
                    lines.append(
                        f"  {filt:^4}  {np.percentile(sub, 10):^8.1f}  "
                        f"{sub.median():^8.1f}  {np.percentile(sub, 90):^8.1f}"
                    )
                lines.append("")

            # 4. qfit
            if "qfit" in good.columns:
                lines.append("  4. PSF 적합 품질 (qfit, flags=0)")
                lines.append(f"  {'필터':^4}  {'중앙값':^8}  {'>5 비율':^8}")
                lines.append("  " + "─" * 28)
                for filt in filters:
                    sub = good[good["FILTER"] == filt]["qfit"].dropna() if "FILTER" in good.columns else pd.Series()
                    if sub.empty:
                        continue
                    bad_pct = 100 * (sub > 5).sum() / max(len(sub), 1)
                    warn = " ⚠" if bad_pct > 5 else ""
                    lines.append(f"  {filt:^4}  {sub.median():^8.3f}  {bad_pct:^7.1f}%{warn}")
                lines.append("")

            # 5. Residual STD
            if not meta_df.empty and "residual_std" in meta_df.columns:
                lines.append("  5. Residual STD (ADU, per frame mean)")
                lines.append(f"  {'필터':^4}  {'iter1':^10}  {'iter2':^10}")
                lines.append("  " + "─" * 30)
                for filt in filters:
                    i1 = meta_df[(meta_df["filter"] == filt) & (meta_df["iter"] == 1)]["residual_std"]
                    i2 = meta_df[(meta_df["filter"] == filt) & (meta_df["iter"] == 2)]["residual_std"]
                    if i1.empty:
                        continue
                    i2_mean = f"{i2.mean():.2f}" if not i2.empty else "N/A"
                    lines.append(f"  {filt:^4}  {i1.mean():^10.2f}  {i2_mean:^10}")
                lines.append("")

            self.qc_text.setPlainText("\n".join(lines))
            export_inputs = (idx, all_df, meta_df)
        except Exception as e:
            self.qc_text.setPlainText(f"QC 생성 오류: {e}")

        self._cmp_merged_df = None
        self._plot_mag_comparison()
        self._refresh_final_diagnostics()
        if export_inputs is not None:
            self._export_psf_qc_products(*export_inputs)

    def _refresh_final_diagnostics(self, *_args) -> None:
        """Refresh the frame selector and redraw the selected final diagnostic."""
        if not hasattr(self, "final_diag_frame_combo"):
            return
        psf_dir = step8_psf_dir(self.params.P.result_dir)
        frames = [
            path.name[len("photometry_"):-len(".tsv")]
            for path in sorted(psf_dir.glob("photometry_*.tsv"))
        ]
        current = self.final_diag_frame_combo.currentText()
        self.final_diag_frame_combo.blockSignals(True)
        self.final_diag_frame_combo.clear()
        self.final_diag_frame_combo.addItems(frames)
        if current in frames:
            self.final_diag_frame_combo.setCurrentText(current)
        elif frames:
            self.final_diag_frame_combo.setCurrentIndex(0)
        self.final_diag_frame_combo.blockSignals(False)
        self._plot_final_diagnostics(self.final_diag_frame_combo.currentText())

    def _show_final_diagnostic_message(self, message: str) -> None:
        self.final_diag_fig.clear()
        ax = self.final_diag_fig.add_subplot(111)
        ax.set_axis_off()
        ax.text(0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes, color="0.4")
        self.final_diag_status.setText(message)
        self.final_diag_status.setStyleSheet("QLabel { color: #616161; }")
        self.final_diag_canvas.draw_idle()

    def _final_diagnostic_meta(self, fname: str) -> dict:
        meta = self._residual_meta.get(fname, {})
        if isinstance(meta, dict) and meta:
            return meta
        path = step8_psf_dir(self.params.P.result_dir) / f"residual_meta_{fname}.json"
        if not path.exists():
            return {}
        try:
            meta = json.loads(path.read_text(encoding="utf-8"))
            return meta if isinstance(meta, dict) else {}
        except Exception:
            return {}

    def _final_diagnostic_epsf(self, fname: str, meta: dict) -> tuple[np.ndarray | None, Path | None]:
        psf_dir = step8_psf_dir(self.params.P.result_dir)
        filt = str(meta.get("filter", "")).strip()
        stem = Path(fname).stem
        candidates = []
        if filt:
            candidates.extend(
                [
                    psf_dir / f"epsf_model_{filt}_{stem}.fits",
                    psf_dir / f"epsf_model_{filt.lower()}_{stem}.fits",
                    psf_dir / f"epsf_model_{filt.upper()}_{stem}.fits",
                    psf_dir / f"epsf_model_{filt}.fits",
                    psf_dir / f"epsf_model_{filt.lower()}.fits",
                    psf_dir / f"epsf_model_{filt.upper()}.fits",
                ]
            )
        candidates.extend(sorted(psf_dir.glob(f"epsf_model_*_{stem}.fits")))
        seen: set[Path] = set()
        for path in candidates:
            if path in seen or not path.exists():
                continue
            seen.add(path)
            try:
                return np.asarray(fits.getdata(path), dtype=float), path
            except Exception:
                continue
        return None, None

    def _final_diagnostic_reference_catalog(self, fname: str, meta: dict) -> pd.DataFrame:
        reference = meta.get("epsf_reference", {})
        catalog_name = reference.get("catalog_path", "") if isinstance(reference, dict) else ""
        path = step8_psf_dir(self.params.P.result_dir) / (
            str(catalog_name).strip() or f"epsf_reference_{fname}.csv"
        )
        if not path.exists():
            return pd.DataFrame()
        try:
            return pd.read_csv(path)
        except Exception:
            return pd.DataFrame()

    def _final_diagnostic_pixel_scale(self, fname: str) -> float:
        fits_path = self._resolve_fits_path_window(fname)
        if fits_path is not None:
            try:
                from astropy.wcs import WCS
                from astropy.wcs.utils import proj_plane_pixel_scales

                celestial = WCS(fits.getheader(fits_path)).celestial
                scales = np.asarray(proj_plane_pixel_scales(celestial), dtype=float) * 3600.0
                scale = float(np.nanmedian(np.abs(scales)))
                if np.isfinite(scale) and scale > 0:
                    return scale
            except Exception:
                pass

        for path in (
            self._cache_dir_path() / f"detect_{fname}.json",
            step4_dir(self.params.P.result_dir) / f"detect_{fname}.json",
        ):
            if not path.exists():
                continue
            try:
                meta = json.loads(path.read_text(encoding="utf-8"))
                fwhm_arcsec = _safe_float(meta.get("fwhm_arcsec"), np.nan)
                fwhm_px = _safe_float(
                    meta.get("fwhm_px", meta.get("fwhm_med_px", meta.get("fwhm_med"))),
                    np.nan,
                )
                scale = fwhm_arcsec / fwhm_px
                if np.isfinite(scale) and scale > 0:
                    return float(scale)
            except Exception:
                continue
        return np.nan

    def _plot_final_diagnostics(self, fname: str | None = None) -> None:
        """Draw the six-panel final diagnostic for one Step 8 frame."""
        if not hasattr(self, "final_diag_fig"):
            return
        if not isinstance(fname, str) or not fname:
            fname = self.final_diag_frame_combo.currentText()
        if not fname:
            self._final_diag_data = pd.DataFrame()
            self._final_diag_summary = {}
            self._show_final_diagnostic_message("No PSF result is available.")
            return

        try:
            # 헤드리스와 **같은 함수**로 그린다 — 조립을 창이 따로 하면 두
            # 산출물이 갈라진다.
            data, summary = render_psf_final_diagnostics(
                self.final_diag_fig,
                self.params,
                Path(self.params.P.result_dir),
                fname,
                use_cropped=bool(getattr(self, "use_cropped", False)),
            )
            self._final_diag_data = data
            self._final_diag_summary = summary
            status = str(summary.get("status", "CHECK"))
            status_parts = [
                status,
                f"N={int(summary.get('n_matched', 0))}",
                f"offset={_safe_float(summary.get('high_snr_reference_offset_mag')):+.3f} mag",
                f"low-SNR={_safe_float(summary.get('low_snr_5_10_median_centered_mag')):+.3f} mag",
                f"high-SNR scatter={_safe_float(summary.get('high_snr_robust_scatter_mag')):.3f} mag",
                f"e={_safe_float(summary.get('epsf_ellipticity')):.3f}",
                f"A180={_safe_float(summary.get('epsf_rotation_asymmetry')):.3f}",
                f"refs={int(summary.get('epsf_reference_n', 0))}",
            ]
            if bool(summary.get("psf_aperture_scale_applied", False)):
                status_parts.append(
                    f"scale={_safe_float(summary.get('psf_aperture_scale'), 1.0):.4f} "
                    f"(N={_to_int(summary.get('psf_aperture_scale_n', 0), 0)})"
                )
            warnings = summary.get("warnings", [])
            if isinstance(warnings, list) and warnings:
                status_parts.append("; ".join(str(item) for item in warnings))
            self.final_diag_status.setText(" | ".join(status_parts))
            color = "#2E7D32" if status == "OK" else "#E65100"
            self.final_diag_status.setStyleSheet(f"QLabel {{ color: {color}; font-weight: bold; }}")
            self.final_diag_canvas.draw_idle()
        except Exception as exc:
            self._final_diag_data = pd.DataFrame()
            self._final_diag_summary = {"file": fname, "status": "ERROR", "error": str(exc)}
            self._show_final_diagnostic_message(f"Final diagnostics failed for {fname}:\n{exc}")
            try:
                self.log(f"Final diagnostics failed for {fname}: {exc}")
            except Exception:
                pass

    def _export_psf_qc_products(
        self,
        idx: pd.DataFrame,
        all_df: pd.DataFrame,
        meta_df: pd.DataFrame,
    ) -> list[Path]:
        """Export reproducible Step 8 QC products for papers and run audits."""
        psf_dir = step8_psf_dir(self.params.P.result_dir)
        psf_dir.mkdir(parents=True, exist_ok=True)

        saved: list[Path] = []
        summary = _build_psf_qc_summary(
            idx,
            all_df,
            meta_df,
            getattr(self, "_cmp_merged_df", None),
        )
        if not summary.empty:
            summary_path = psf_dir / "psf_qc_summary.csv"
            summary.to_csv(summary_path, index=False)
            saved.append(summary_path)

        frame_qc = _build_psf_frame_qc_table(idx, meta_df)
        if not frame_qc.empty:
            frame_qc_path = psf_dir / "psf_frame_qc.csv"
            frame_qc.to_csv(frame_qc_path, index=False)
            saved.append(frame_qc_path)

            frame_fig = Figure(figsize=(10.5, 6.8), dpi=120)
            if _draw_psf_frame_qc_overview(frame_fig, frame_qc):
                frame_fig_path = psf_dir / "step8_residual_core_qc.png"
                frame_fig.savefig(frame_fig_path, dpi=160, bbox_inches="tight")
                saved.append(frame_fig_path)

        report = self.qc_text.toPlainText() if hasattr(self, "qc_text") else ""
        if report.strip():
            report_path = psf_dir / "psf_qc_report.txt"
            report_path.write_text(report, encoding="utf-8")
            saved.append(report_path)

        if hasattr(self, "cmp_fig"):
            fig_path = psf_dir / "step8_ap_vs_psf_comparison.png"
            self.cmp_fig.savefig(fig_path, dpi=160, bbox_inches="tight")
            saved.append(fig_path)

        final_summary = getattr(self, "_final_diag_summary", {})
        final_data = getattr(self, "_final_diag_data", pd.DataFrame())
        if isinstance(final_summary, dict) and final_summary.get("file"):
            stem = Path(str(final_summary["file"])).stem
            final_fig_path = psf_dir / f"step8_final_diagnostics_{stem}.png"
            self.final_diag_fig.savefig(final_fig_path, dpi=160, bbox_inches="tight")
            saved.append(final_fig_path)

            summary_path = psf_dir / f"psf_final_diagnostics_{stem}.json"
            summary_path.write_text(
                json.dumps(final_summary, ensure_ascii=False, indent=2, allow_nan=True),
                encoding="utf-8",
            )
            saved.append(summary_path)
            if isinstance(final_data, pd.DataFrame) and not final_data.empty:
                data_path = psf_dir / f"psf_final_diagnostics_{stem}.csv"
                final_data.to_csv(data_path, index=False)
                saved.append(data_path)

        if saved:
            try:
                names = ", ".join(path.name for path in saved)
                self.log(f"Step8 QC products exported: {names}")
            except Exception:
                pass
        return saved

    # ── Aperture vs PSF magnitude comparison ──────────────────────────────────

    def _load_or_build_comparison(self) -> tuple[pd.DataFrame, int]:
        """구경 vs PSF 병합표를 디스크 캐시에서 읽고, 없거나 낡았으면 다시 만든다.

        이 병합은 프레임 수에 비례해 무겁다 — M13 15프레임에 **10.9초**이고,
        창을 열 때마다 `_refresh_qc()` 가 캐시를 버리고 처음부터 다시 만든다.
        그래서 Step 8 창 로드가 17.8초였다. Step 8 산출물이 그대로면 결과도
        같으므로, `photometry_index.csv` 보다 새 캐시가 있으면 그것을 쓴다.
        """
        psf_dir = step8_psf_dir(self.params.P.result_dir)
        cache_path = psf_dir / "psf_ap_vs_psf.csv"
        meta_path = psf_dir / "psf_ap_vs_psf_meta.json"
        index_path = psf_dir / "photometry_index.csv"
        try:
            if (
                cache_path.exists()
                and index_path.exists()
                and cache_path.stat().st_mtime >= index_path.stat().st_mtime
            ):
                cached = pd.read_csv(cache_path)
                n_split = 0
                if meta_path.exists():
                    n_split = int(
                        json.loads(meta_path.read_text(encoding="utf-8")).get(
                            "split_excluded_total", 0
                        )
                    )
                return cached, n_split
        except Exception:
            pass

        merged, split_excluded_total = build_ap_psf_comparison(
            self.params, self.params.P.result_dir
        )
        try:
            psf_dir.mkdir(parents=True, exist_ok=True)
            merged.to_csv(cache_path, index=False)
            meta_path.write_text(
                json.dumps({"split_excluded_total": int(split_excluded_total)}),
                encoding="utf-8",
            )
        except Exception:
            pass
        return merged, int(split_excluded_total)

    def _plot_mag_comparison(self):  # noqa: C901
        """Scatter: mag_ap (Step5) vs mag_psf (Step6), merged on det_uid."""
        if not hasattr(self, "cmp_fig"):
            return

        _FILT_COLORS = {
            "u": "#9467bd", "g": "#2ca02c", "r": "#d62728",
            "i": "#ff7f0e", "z": "#8c564b", "b": "#1f77b4",
            "v": "#bcbd22", "ha": "#e377c2",
        }

        # 병합은 헤드리스와 **같은 함수**로 한다 — 창과 배치가 각자 병합하면
        # 두 산출물이 조용히 갈라진다. 캐시(_cmp_merged_df)는 창 쪽 사정이다.
        if not hasattr(self, "_cmp_merged_df") or self._cmp_merged_df is None:
            merged, split_excluded_total = self._load_or_build_comparison()
            self._cmp_merged_df = merged
            self._cmp_split_excluded_total = int(split_excluded_total)

        self.cmp_fig.clf()

        def _empty(msg):
            ax = self.cmp_fig.add_subplot(111)
            ax.text(0.5, 0.5, msg, transform=ax.transAxes,
                    ha="center", va="center", fontsize=10, color="gray")
            self.cmp_canvas.draw_idle()
            self.cmp_stats_label.setText(msg)

        if self._cmp_merged_df.empty:
            _empty("No matched data.\nRun Step 7 and Step 8 first.")
            return

        df = self._cmp_merged_df.copy()

        # Refresh filter/frame selectors from available merged data.
        if hasattr(self, "cmp_filter_combo"):
            try:
                _fvals = sorted(df["FILTER"].dropna().astype(str).unique().tolist()) if "FILTER" in df.columns else []
                _cur = self.cmp_filter_combo.currentText().strip() or "all"
                self.cmp_filter_combo.blockSignals(True)
                self.cmp_filter_combo.clear()
                self.cmp_filter_combo.addItem("all")
                for _v in _fvals:
                    self.cmp_filter_combo.addItem(_v)
                self.cmp_filter_combo.setCurrentText(_cur if _cur in (["all"] + _fvals) else "all")
                self.cmp_filter_combo.blockSignals(False)
            except Exception:
                pass
        if hasattr(self, "cmp_frame_combo"):
            try:
                _frames = sorted(df["FRAME"].dropna().astype(str).unique().tolist()) if "FRAME" in df.columns else []
                _cur = self.cmp_frame_combo.currentText().strip() or "all"
                self.cmp_frame_combo.blockSignals(True)
                self.cmp_frame_combo.clear()
                self.cmp_frame_combo.addItem("all")
                for _v in _frames:
                    self.cmp_frame_combo.addItem(_v)
                self.cmp_frame_combo.setCurrentText(_cur if _cur in (["all"] + _frames) else "all")
                self.cmp_frame_combo.blockSignals(False)
            except Exception:
                pass

        df["mag_ap"] = pd.to_numeric(df.get("mag_ap"), errors="coerce")
        df["mag_psf"] = pd.to_numeric(df.get("mag_psf"), errors="coerce")
        df = df[np.isfinite(df["mag_ap"]) & np.isfinite(df["mag_psf"])].copy()

        if len(df) == 0:
            _empty("All magnitudes are NaN.\nCheck Step 7 forced photometry and Step 8 PSF outputs.")
            return

        df["delta"] = df["mag_ap"] - df["mag_psf"]
        n_before = int(len(df))

        # Pre-convert numeric filter columns once
        if "flags_psf" in df.columns:
            df["flags_psf"] = pd.to_numeric(df["flags_psf"], errors="coerce")
        if "snr_psf" in df.columns:
            df["snr_psf"] = pd.to_numeric(df["snr_psf"], errors="coerce")
        if "qfit" in df.columns:
            df["qfit"] = pd.to_numeric(df["qfit"], errors="coerce")
        if "qfit_noise_ratio" in df.columns:
            df["qfit_noise_ratio"] = pd.to_numeric(
                df["qfit_noise_ratio"], errors="coerce"
            )

        # Selector filters (filter/frame)
        if hasattr(self, "cmp_filter_combo") and "FILTER" in df.columns:
            _fsel = str(self.cmp_filter_combo.currentText()).strip()
            if _fsel and _fsel.lower() != "all":
                _fkey = normalize_filter_name(_fsel)
                df = df[df["FILTER"].astype(str).map(normalize_filter_name) == _fkey].copy()
        if hasattr(self, "cmp_frame_combo") and "FRAME" in df.columns:
            _rsel = str(self.cmp_frame_combo.currentText()).strip()
            if _rsel and _rsel.lower() != "all":
                df = df[df["FRAME"].astype(str) == _rsel].copy()

        # User filters
        if getattr(self, "cmp_flags0_only", None) is not None and self.cmp_flags0_only.isChecked() and "flags_psf" in df.columns:
            df = df[np.isfinite(df["flags_psf"]) & (df["flags_psf"] == 0)].copy()
        if getattr(self, "cmp_snr_min", None) is not None:
            _snr_min = float(self.cmp_snr_min.value())
            if _snr_min > 0 and "snr_psf" in df.columns:
                df = df[np.isfinite(df["snr_psf"]) & (df["snr_psf"] >= _snr_min)].copy()
        if getattr(self, "cmp_qfit_max", None) is not None:
            _qmax = float(self.cmp_qfit_max.value())
            _qfit_column = (
                "qfit_noise_ratio"
                if "qfit_noise_ratio" in df.columns
                and np.isfinite(df["qfit_noise_ratio"]).any()
                else "qfit"
            )
            if _qmax > 0 and _qfit_column in df.columns:
                df = df[
                    np.isfinite(df[_qfit_column]) & (df[_qfit_column] <= _qmax)
                ].copy()
        if getattr(self, "cmp_dmag_clip", None) is not None:
            _dclip = float(self.cmp_dmag_clip.value())
            if _dclip > 0:
                df = df[np.isfinite(df["delta"]) & (np.abs(df["delta"]) <= _dclip)].copy()

        if len(df) == 0:
            _empty("No data after filters.\nRelax SNR/qfit/|Δmag| settings.")
            return

        filt_col = "FILTER" if "FILTER" in df.columns else None

        ax1 = self.cmp_fig.add_subplot(121)
        ax2 = self.cmp_fig.add_subplot(122)

        stats_parts = []
        groups = df.groupby(filt_col) if filt_col else [("all", df)]
        for filt, sub in groups:
            color = _FILT_COLORS.get(str(filt).lower(), "#999999")
            ax1.scatter(sub["mag_ap"], sub["mag_psf"],
                        s=4, alpha=0.35, color=color, label=str(filt), rasterized=True)
            ax2.scatter(sub["mag_ap"], sub["delta"],
                        s=4, alpha=0.35, color=color, label=str(filt), rasterized=True)
            n = len(sub)
            med = float(np.nanmedian(sub["delta"]))
            std = float(np.nanstd(sub["delta"]))
            stats_parts.append(f"{filt}: N={n}  Δmed={med:+.3f}  σ={std:.3f}")

        # 1:1 reference line (ax1)
        all_mag = np.concatenate([df["mag_ap"].values, df["mag_psf"].values])
        lo, hi = np.nanmin(all_mag) - 0.2, np.nanmax(all_mag) + 0.2
        ax1.plot([lo, hi], [lo, hi], "k--", lw=0.8, alpha=0.5, zorder=0)
        ax1.set_xlim(lo, hi)
        ax1.set_ylim(lo, hi)
        ax1.set_xlabel("mag_ap", fontsize=9)
        ax1.set_ylabel("mag_psf", fontsize=9)
        ax1.set_title("Aperture vs PSF magnitude", fontsize=9)
        ax1.legend(fontsize=7, markerscale=2, loc="upper left")

        # Zero ± 0.05 reference lines (ax2) + robust median guide
        dmed_all = float(np.nanmedian(df["delta"]))
        ax2.axhline(
            dmed_all,
            color="#D62728",
            lw=2.0,
            ls="-",
            alpha=0.95,
            zorder=1,
            label=f"Δmag median {dmed_all:+.3f}",
        )
        ax2.axhline(0.0,  color="k",    lw=0.8, ls="--", alpha=0.6, zorder=0)
        ax2.axhline(+0.05, color="gray", lw=0.5, ls=":",  alpha=0.5, zorder=0)
        ax2.axhline(-0.05, color="gray", lw=0.5, ls=":",  alpha=0.5, zorder=0)
        ax2.set_xlabel("mag_ap", fontsize=9)
        ax2.set_ylabel("Δmag  (Ap − PSF)", fontsize=9)
        ax2.set_title("Δmag vs mag_ap", fontsize=9)
        ax2.legend(fontsize=7, markerscale=2, loc="upper left")

        self.cmp_fig.tight_layout()
        self.cmp_canvas.draw_idle()
        self.cmp_stats_label.setText(
            f"N={len(df)}/{n_before}  |  split_excluded={int(getattr(self, '_cmp_split_excluded_total', 0))}  |  "
            + "  |  ".join(stats_parts)
        )

    # ── Load existing results from disk (called on restore_state) ─────────────

    def _load_from_disk(self):
        """Reload EPSF models and residual images from disk into memory caches."""
        out_dir = step8_psf_dir(self.params.P.result_dir)
        if not out_dir.exists():
            return
        self._last_epsf.clear()
        self._residual_meta.clear()
        self._last_new_xy.clear()
        self.epsf_filter_combo.clear()
        self.res_file_combo.clear()
        self.res_iter_combo.clear()

        def _epsf_display_key_from_path(epsf_path: Path) -> str:
            stem = epsf_path.stem  # epsf_model_{filter}_{frame_stem} or epsf_model_{filter}
            body = stem.replace("epsf_model_", "", 1)
            if "_" not in body:
                return body
            filt, frame_stem = body.split("_", 1)
            return f"{frame_stem} | {filt}"

        # Load EPSF FITS files
        for epsf_path in out_dir.glob("epsf_model_*.fits"):
            try:
                display_key = _epsf_display_key_from_path(epsf_path)
                if not display_key:
                    continue
                arr = fits.getdata(str(epsf_path)).astype(float)
                self._last_epsf[display_key] = arr
            except Exception:
                pass

        if self._last_epsf:
            self.epsf_filter_combo.blockSignals(True)
            self.epsf_filter_combo.clear()
            self.epsf_filter_combo.addItems(sorted(self._last_epsf.keys()))
            self.epsf_filter_combo.setCurrentIndex(0)
            self.epsf_filter_combo.blockSignals(False)
            first_filter = self.epsf_filter_combo.currentText()
            if first_filter:
                self._plot_epsf(first_filter)

        # Load residual metadata (preferred, supports iteration-wise view).
        for meta_path in sorted(out_dir.glob("residual_meta_*.json")):
            try:
                name = meta_path.name
                if not name.startswith("residual_meta_") or not name.endswith(".json"):
                    continue
                fname = name[len("residual_meta_"):-len(".json")]
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                if isinstance(meta, dict):
                    self._residual_meta[fname] = meta
                    self._last_new_xy.setdefault(fname, None)
            except Exception:
                pass

        if self._residual_meta:
            self.res_file_combo.blockSignals(True)
            self.res_file_combo.clear()
            self.res_file_combo.addItems(sorted(self._residual_meta.keys()))
            self.res_file_combo.setCurrentIndex(0)
            self.res_file_combo.blockSignals(False)
            first_fname = self.res_file_combo.currentText()
            if first_fname:
                self._refresh_residual_iter_combo(first_fname)
                self._plot_cutout(first_fname)

        self._refresh_qc()  # refresh QC tab (stats + Ap vs PSF plot)

    # ── Parameters dialog ─────────────────────────────────────────────────────

    def open_parameters_dialog(self):
        dialog = FittedDialog(self)
        configure_parameter_dialog(dialog, "Step 8 PSF Parameters", 620, 720)
        layout = QVBoxLayout(dialog)

        scroll = QScrollArea(dialog)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        content = QWidget()
        body = QVBoxLayout(content)
        body.setContentsMargins(4, 4, 4, 4)
        body.setSpacing(8)
        scroll.setWidget(content)
        layout.addWidget(scroll, 1)

        _info = QLabel("Adjust PSF photometry parameters. Changes apply to the next run.")
        _info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 10px; margin-bottom: 10px; }")
        _info.setWordWrap(True)
        body.addWidget(_info)

        def _add_group(title: str, *, expanded: bool = False) -> QFormLayout:
            group, container = create_collapsible_section(title, initial_expanded=expanded)
            form = QFormLayout(container)
            form.setLabelAlignment(Qt.AlignRight)
            form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
            body.addWidget(group)
            return form

        mode_form = _add_group("Mode", expanded=True)
        epsf_form = _add_group("ePSF Model", expanded=True)
        scale_form = _add_group("PSF Flux Scale")
        fit_form = _add_group("PSF Fit")
        core_form = _add_group("Crowded Core Cut")
        redetect_form = _add_group("Residual Re-detection")
        output_form = _add_group("Output")

        # ── Field mode preset ──────────────────────────────────────────────
        mode_combo = QComboBox()
        mode_combo.addItem("Normal (일반)", "normal")
        mode_combo.addItem("Crowded (구상성단/혼잡장)", "crowded")
        mode_combo.addItem("Faint (희미한 필드)", "faint")
        mode_combo.addItem("Custom (수동)", "custom")
        _saved_mode = str(getattr(self.params.P, "psf_mode", "normal"))
        _mi = mode_combo.findData(_saved_mode)
        mode_combo.setCurrentIndex(_mi if _mi >= 0 else 0)
        mode_form.addRow("Field mode:", mode_combo)

        self.p_model_mode = QComboBox()
        self.p_model_mode.addItems(["per_frame"])
        self.p_model_mode.setCurrentText("per_frame")

        self.p_fit_engine = QComboBox()
        self.p_fit_engine.addItem("APEX iterative - CPU, recommended", "apex_iterative")
        self.p_fit_engine.addItem("photutils - validation, slower", "photutils")
        _eng = str(
            getattr(self.params.P, "psf_fit_engine", "apex_iterative")
        ).strip().lower()
        if _eng == "allstar":
            _eng = "apex_iterative"
        _ei = self.p_fit_engine.findData(_eng)
        self.p_fit_engine.setCurrentIndex(_ei if _ei >= 0 else 0)
        mode_form.addRow("Fit engine:", self.p_fit_engine)

        self.p_build_mode = QComboBox()
        self.p_build_mode.addItem("epsf  —  EPSFBuilder 경험적 PSF", "epsf")
        self.p_build_mode.addItem("moffat  —  해석 Moffat (γ·β 적합)", "moffat")
        self.p_build_mode.addItem("moffat+residual  —  해석 + 잔차격자", "moffat_hybrid")
        self.p_build_mode.setToolTip(
            "epsf: 별 이미지를 그대로 평균낸 경험적 모형. 어떤 모양이든 담지만"
            " 가파른 핵심까지 보간해야 한다.\n"
            "moffat: γ·β 를 적합한 해석 모형. 핵심을 정확히 계산하지만"
            " 원형이라 늘어난 별은 못 맞춘다.\n"
            "moffat+residual: 핵심은 해석식으로 정확히, 남은 모양은 격자로."
            " DAOPHOT 과 같은 구성."
        )
        _bm = str(getattr(self.params.P, "psf_build_mode", "epsf")).strip().lower()
        _bi = self.p_build_mode.findData(_bm)
        self.p_build_mode.setCurrentIndex(_bi if _bi >= 0 else 0)

        self.p_workers = QSpinBox()
        self.p_workers.setRange(0, 64)
        self.p_workers.setValue(_to_int(getattr(self.params.P, "psf_parallel_workers", 0), 0))
        self.p_workers.setToolTip("0 = auto/global parallel workers")
        mode_form.addRow("PSF workers (0=auto):", self.p_workers)

        self.p_oversampling = QSpinBox()
        self.p_oversampling.setRange(1, 8)
        self.p_oversampling.setValue(_to_int(getattr(self.params.P, "psf_epsf_oversampling", 2), 2))
        epsf_form.addRow("EPSF oversampling:", self.p_oversampling)

        self.p_epsf_mult = QDoubleSpinBox()
        self.p_epsf_mult.setRange(1.0, 10.0)
        self.p_epsf_mult.setSingleStep(0.1)
        self.p_epsf_mult.setValue(_to_float(getattr(self.params.P, "psf_epsf_size_fwhm_mult", 4.0), 4.0))
        epsf_form.addRow("EPSF cutout (×FWHM):", self.p_epsf_mult)

        self.p_n_stars = QSpinBox()
        self.p_n_stars.setRange(0, 500)
        self.p_n_stars.setSpecialValueText("auto")
        self.p_n_stars.setValue(_to_int(getattr(self.params.P, "psf_n_stars_max", 0), 0))
        self.p_n_stars.setToolTip(
            "0 = automatic per-frame budget; positive values are a CPU/memory cap (minimum 3)"
        )
        epsf_form.addRow("PSF star cap:", self.p_n_stars)

        self.p_isolation = QDoubleSpinBox()
        self.p_isolation.setRange(1.0, 10.0)
        self.p_isolation.setSingleStep(0.5)
        self.p_isolation.setValue(_to_float(getattr(self.params.P, "psf_isolation_fwhm_mult", 3.0), 3.0))
        epsf_form.addRow("Isolation (×FWHM):", self.p_isolation)

        self.p_epsf_contamination_filter = QCheckBox(
            "Reject locally contaminated/core ePSF reference stars"
        )
        self.p_epsf_contamination_filter.setChecked(
            _as_bool(getattr(self.params.P, "psf_epsf_contamination_filter", True), True)
        )
        self.p_epsf_contamination_filter.setToolTip(
            "Uses all Step 4 detections, local annulus residuals, and the automatic "
            "cluster-core estimate only for ePSF reference-star selection."
        )
        epsf_form.addRow("", self.p_epsf_contamination_filter)

        self.p_flux_scale_correction = QCheckBox(
            "Anchor PSF fluxes to clean Step 7 aperture references"
        )
        self.p_flux_scale_correction.setChecked(
            _as_bool(getattr(self.params.P, "psf_flux_scale_correction", False), False)
        )
        self.p_flux_scale_correction.setToolTip(
            "Applies one robust per-frame multiplicative scale after PSF fitting. "
            "Raw PSF fluxes are preserved in separate output columns."
        )
        scale_form.addRow("", self.p_flux_scale_correction)

        self.p_flux_scale_min_snr = QDoubleSpinBox()
        self.p_flux_scale_min_snr.setRange(5.0, 1000.0)
        self.p_flux_scale_min_snr.setSingleStep(10.0)
        self.p_flux_scale_min_snr.setValue(
            _to_float(getattr(self.params.P, "psf_flux_scale_min_snr", 50.0), 50.0)
        )
        scale_form.addRow("Minimum aperture SNR:", self.p_flux_scale_min_snr)

        self.p_flux_scale_min_stars = QSpinBox()
        self.p_flux_scale_min_stars.setRange(3, 500)
        self.p_flux_scale_min_stars.setValue(
            _to_int(getattr(self.params.P, "psf_flux_scale_min_stars", 8), 8)
        )
        scale_form.addRow("Minimum references:", self.p_flux_scale_min_stars)

        self.p_flux_scale_min_neighbor = QDoubleSpinBox()
        self.p_flux_scale_min_neighbor.setRange(0.0, 20.0)
        self.p_flux_scale_min_neighbor.setSingleStep(0.5)
        self.p_flux_scale_min_neighbor.setValue(
            _to_float(
                getattr(self.params.P, "psf_flux_scale_min_neighbor_fwhm", 4.0),
                4.0,
            )
        )
        scale_form.addRow("Minimum neighbor distance (xFWHM):", self.p_flux_scale_min_neighbor)

        self.p_flux_scale_max_scatter = QDoubleSpinBox()
        self.p_flux_scale_max_scatter.setRange(0.01, 1.0)
        self.p_flux_scale_max_scatter.setSingleStep(0.01)
        self.p_flux_scale_max_scatter.setDecimals(3)
        self.p_flux_scale_max_scatter.setValue(
            _to_float(
                getattr(self.params.P, "psf_flux_scale_max_scatter_mag", 0.10),
                0.10,
            )
        )
        scale_form.addRow("Maximum reference scatter (mag):", self.p_flux_scale_max_scatter)

        self.p_fit_window_mode = QComboBox()
        self.p_fit_window_mode.addItem("Auto (PSF energy)", "auto")
        self.p_fit_window_mode.addItem("Manual (FWHM multiplier)", "manual")
        _fit_window_mode = str(
            getattr(self.params.P, "psf_fit_window_mode", "auto")
        ).strip().lower()
        _fit_window_index = self.p_fit_window_mode.findData(_fit_window_mode)
        self.p_fit_window_mode.setCurrentIndex(
            _fit_window_index if _fit_window_index >= 0 else 0
        )
        fit_form.addRow("Fit window mode:", self.p_fit_window_mode)

        self.p_fit_energy = QDoubleSpinBox()
        self.p_fit_energy.setRange(0.50, 0.995)
        self.p_fit_energy.setSingleStep(0.01)
        self.p_fit_energy.setDecimals(3)
        self.p_fit_energy.setValue(
            _to_float(
                getattr(self.params.P, "psf_fit_encircled_energy", 0.90), 0.90
            )
        )
        fit_form.addRow("Target PSF energy:", self.p_fit_energy)

        self.p_fit_mult = QDoubleSpinBox()
        self.p_fit_mult.setRange(0.5, 5.0)
        self.p_fit_mult.setSingleStep(0.1)
        self.p_fit_mult.setValue(
            _to_float(getattr(self.params.P, "psf_fit_shape_fwhm_mult", 2.4), 2.4)
        )
        fit_form.addRow("Manual fit window (xFWHM):", self.p_fit_mult)

        def _sync_fit_window_controls() -> None:
            automatic = self.p_fit_window_mode.currentData() == "auto"
            self.p_fit_energy.setEnabled(automatic)
            self.p_fit_mult.setEnabled(not automatic)

        self.p_fit_window_mode.currentIndexChanged.connect(
            lambda _index: _sync_fit_window_controls()
        )
        _sync_fit_window_controls()

        self.p_max_iter = QSpinBox()
        self.p_max_iter.setRange(1, 3)
        self.p_max_iter.setValue(_to_int(getattr(self.params.P, "psf_max_iter", 2), 2))
        fit_form.addRow("Residual passes:", self.p_max_iter)

        self.p_fitter_max_iter = QSpinBox()
        self.p_fitter_max_iter.setRange(1, 10)
        self.p_fitter_max_iter.setValue(
            _to_int(getattr(self.params.P, "psf_fitter_max_iter", 6), 6)
        )
        self.p_fitter_max_iter.setToolTip(
            "Maximum weighted Newton updates inside each residual pass"
        )
        fit_form.addRow("Fitter updates/pass:", self.p_fitter_max_iter)

        # The pass that sets every published flux. Two steps solve an isolated
        # star and do not solve a blended group; ALLSTAR allows 50 for the same
        # solve. Convergence still exits early, so a high ceiling only costs
        # time where it is actually used.
        self.p_final_pass_max_iter = QSpinBox()
        self.p_final_pass_max_iter.setRange(1, 100)
        self.p_final_pass_max_iter.setValue(
            _to_int(getattr(self.params.P, "psf_final_pass_max_iter", 2), 2)
        )
        self.p_final_pass_max_iter.setToolTip(
            "Newton updates allowed in the final fixed-position pass. Blended "
            "stars need more than a couple; isolated stars stop early anyway."
        )
        fit_form.addRow("Final pass updates:", self.p_final_pass_max_iter)

        # Catalog positions are worth trusting for an isolated star. In a blend
        # they make the fit asymmetric: only the neighbour can move, so it
        # absorbs the pair's misfit and the catalog star comes out faint.
        self.p_forced_position_lock = QComboBox()
        self.p_forced_position_lock.addItem("Always (catalog positions fixed)", "always")
        self.p_forced_position_lock.addItem("Never (fit every position)", "never")
        _lock_now = str(getattr(self.params.P, "psf_forced_position_lock", "always"))
        _lock_idx = self.p_forced_position_lock.findData(_lock_now)
        self.p_forced_position_lock.setCurrentIndex(max(0, _lock_idx))
        self.p_forced_position_lock.setToolTip(
            "Whether catalog (forced) sources keep their positions during the "
            "fit. Releasing them puts both members of a blend on equal terms."
        )
        fit_form.addRow("Catalog position lock:", self.p_forced_position_lock)

        self.p_redetect = QDoubleSpinBox()
        self.p_redetect.setRange(1.0, 10.0)
        self.p_redetect.setSingleStep(0.5)
        self.p_redetect.setValue(_to_float(getattr(self.params.P, "psf_redetect_sigma", 4.0), 4.0))
        redetect_form.addRow("Re-detect sigma (base):", self.p_redetect)

        def _make_filter_sigma_spin(attr, label):
            sp = QDoubleSpinBox()
            sp.setRange(0.0, 10.0)
            sp.setSingleStep(0.5)
            sp.setDecimals(1)
            sp.setSpecialValueText("base")
            _v = _to_float(getattr(self.params.P, attr, float("nan")), float("nan"))
            sp.setValue(0.0 if not np.isfinite(_v) else float(_v))
            sp.setToolTip("0 = use base sigma")
            redetect_form.addRow(label, sp)
            return sp

        self.p_redetect_g = _make_filter_sigma_spin("psf_redetect_sigma_g", "g-band override:")
        self.p_redetect_r = _make_filter_sigma_spin("psf_redetect_sigma_r", "r-band override:")
        self.p_redetect_i = _make_filter_sigma_spin("psf_redetect_sigma_i", "i-band override:")

        self.p_dup_mult = QDoubleSpinBox()
        self.p_dup_mult.setRange(0.0, 5.0)
        self.p_dup_mult.setSingleStep(0.1)
        self.p_dup_mult.setValue(_to_float(getattr(self.params.P, "psf_duplicate_radius_fwhm_mult", 0.8), 0.8))
        redetect_form.addRow("Duplicate radius (×FWHM):", self.p_dup_mult)

        self.p_dup_px = QDoubleSpinBox()
        self.p_dup_px.setRange(0.0, 50.0)
        self.p_dup_px.setSingleStep(0.1)
        self.p_dup_px.setDecimals(2)
        _dup_px = _to_float(getattr(self.params.P, "psf_duplicate_radius_px", np.nan), np.nan)
        self.p_dup_px.setValue(0.0 if not np.isfinite(_dup_px) else float(_dup_px))
        self.p_dup_px.setToolTip("0이면 비활성(×FWHM 값 사용), >0이면 절대 px 반경 사용")
        redetect_form.addRow("Duplicate radius (px override):", self.p_dup_px)

        self.p_cap_per_iter = QSpinBox()
        self.p_cap_per_iter.setRange(0, 50000)
        self.p_cap_per_iter.setSingleStep(50)
        self.p_cap_per_iter.setValue(_to_int(getattr(self.params.P, "psf_new_sources_cap_per_iter", 70), 70))
        redetect_form.addRow("Max new/iter (abs):", self.p_cap_per_iter)

        self.p_cap_frac = QDoubleSpinBox()
        self.p_cap_frac.setRange(0.0, 1.0)
        self.p_cap_frac.setSingleStep(0.01)
        self.p_cap_frac.setValue(_to_float(getattr(self.params.P, "psf_new_sources_cap_frac", 0.02), 0.02))
        redetect_form.addRow("Max new/iter (frac):", self.p_cap_frac)

        self.p_blend_ratio = QDoubleSpinBox()
        self.p_blend_ratio.setRange(0.0, 1.0)
        self.p_blend_ratio.setDecimals(2)
        self.p_blend_ratio.setSingleStep(0.05)
        self.p_blend_ratio.setValue(
            _to_float(getattr(self.params.P, "psf_blend_residual_ratio", 0.3), 0.3)
        )
        self.p_blend_ratio.setToolTip(
            "Reject residual peaks that are too weak relative to the current source model; 0 disables"
        )
        redetect_form.addRow("Residual/model minimum:", self.p_blend_ratio)

        self.p_postfit_snr = QDoubleSpinBox()
        self.p_postfit_snr.setRange(0.0, 100.0)
        self.p_postfit_snr.setDecimals(1)
        self.p_postfit_snr.setSingleStep(0.5)
        self.p_postfit_snr.setValue(
            _to_float(getattr(self.params.P, "psf_postfit_snr_min", 3.0), 3.0)
        )
        self.p_postfit_snr.setToolTip(
            "New residual detections below this fitted S/N are removed; initial Step 4 sources are retained"
        )
        redetect_form.addRow("Post-fit S/N minimum:", self.p_postfit_snr)

        self.p_postfit_qfit = QDoubleSpinBox()
        self.p_postfit_qfit.setRange(0.0, 100.0)
        self.p_postfit_qfit.setDecimals(2)
        self.p_postfit_qfit.setSingleStep(0.1)
        self.p_postfit_qfit.setValue(
            _to_float(getattr(self.params.P, "psf_postfit_qfit_max", 3.0), 3.0)
        )
        self.p_postfit_qfit.setToolTip(
            "Remove new residual sources above qfit / expected-noise qfit; "
            "0 disables. Initial Step 4 sources are retained."
        )
        redetect_form.addRow("Post-fit qfit/noise maximum:", self.p_postfit_qfit)

        self.p_postfit_redchi = QDoubleSpinBox()
        self.p_postfit_redchi.setRange(0.0, 100000.0)
        self.p_postfit_redchi.setDecimals(1)
        self.p_postfit_redchi.setSingleStep(5.0)
        self.p_postfit_redchi.setValue(
            _to_float(
                getattr(self.params.P, "psf_postfit_reduced_chi2_max", 25.0),
                25.0,
            )
        )
        self.p_postfit_redchi.setToolTip(
            "Remove newly detected sources above this reduced chi-square; 0 disables."
        )
        redetect_form.addRow("Post-fit reduced chi2 maximum:", self.p_postfit_redchi)

        self.p_fit_init_max = QSpinBox()
        self.p_fit_init_max.setRange(0, 200000)
        self.p_fit_init_max.setSingleStep(100)
        self.p_fit_init_max.setValue(_to_int(getattr(self.params.P, "psf_fit_init_max_sources", 0), 0))
        self.p_fit_init_max.setToolTip("0이면 초기 피팅 소스 무제한")
        fit_form.addRow("Initial fit source cap (0=off):", self.p_fit_init_max)

        self.p_core_enable = QCheckBox("Hard-exclude crowded core during PSF fit")
        self.p_core_enable.setChecked(bool(getattr(self.params.P, "psf_core_cut_enable", False)))
        self.p_core_enable.setToolTip(
            "Optional CPU/quality safeguard. Leave off to fit the full field; unresolved pairs are retained with a crowding flag."
        )
        core_form.addRow("", self.p_core_enable)

        self.p_core_center_mode = QComboBox()
        self.p_core_center_mode.addItem("Auto density peak", "auto")
        self.p_core_center_mode.addItem("Image center", "image")
        self.p_core_center_mode.addItem("Manual x/y", "manual")
        _core_mode = str(getattr(self.params.P, "psf_core_cut_center_mode", "auto")).strip().lower()
        _core_mode_i = self.p_core_center_mode.findData(_core_mode)
        self.p_core_center_mode.setCurrentIndex(_core_mode_i if _core_mode_i >= 0 else 0)
        core_form.addRow("Center:", self.p_core_center_mode)

        self.p_core_x = QDoubleSpinBox()
        self.p_core_x.setRange(0.0, 200000.0)
        self.p_core_x.setDecimals(1)
        self.p_core_x.setSingleStep(10.0)
        _cx = _to_float(getattr(self.params.P, "psf_core_cut_x_px", 0.0), 0.0)
        self.p_core_x.setValue(0.0 if not np.isfinite(_cx) else float(_cx))
        core_form.addRow("Manual center x (px):", self.p_core_x)

        self.p_core_y = QDoubleSpinBox()
        self.p_core_y.setRange(0.0, 200000.0)
        self.p_core_y.setDecimals(1)
        self.p_core_y.setSingleStep(10.0)
        _cy = _to_float(getattr(self.params.P, "psf_core_cut_y_px", 0.0), 0.0)
        self.p_core_y.setValue(0.0 if not np.isfinite(_cy) else float(_cy))
        core_form.addRow("Manual center y (px):", self.p_core_y)

        self.p_core_radius_px = QDoubleSpinBox()
        self.p_core_radius_px.setRange(0.0, 200000.0)
        self.p_core_radius_px.setDecimals(1)
        self.p_core_radius_px.setSingleStep(10.0)
        self.p_core_radius_px.setSpecialValueText("auto")
        self.p_core_radius_px.setValue(_to_float(getattr(self.params.P, "psf_core_cut_radius_px", 0.0), 0.0))
        self.p_core_radius_px.setToolTip("0 = estimate radius from the detection-density profile")
        core_form.addRow("Cut radius (px):", self.p_core_radius_px)

        self.p_core_radius_mult = QDoubleSpinBox()
        self.p_core_radius_mult.setRange(1.0, 200.0)
        self.p_core_radius_mult.setDecimals(1)
        self.p_core_radius_mult.setSingleStep(1.0)
        self.p_core_radius_mult.setValue(_to_float(getattr(self.params.P, "psf_core_cut_radius_fwhm_mult", 20.0), 20.0))
        self.p_core_radius_mult.setToolTip("Fallback and safety cap for the automatic core radius")
        core_form.addRow("Auto radius cap (xFWHM):", self.p_core_radius_mult)

        self.p_core_density_ratio = QDoubleSpinBox()
        self.p_core_density_ratio.setRange(1.0, 20.0)
        self.p_core_density_ratio.setDecimals(2)
        self.p_core_density_ratio.setSingleStep(0.1)
        self.p_core_density_ratio.setValue(
            _to_float(getattr(self.params.P, "psf_core_cut_auto_min_density_ratio", 1.5), 1.5)
        )
        self.p_core_density_ratio.setToolTip("Auto center/cut is disabled when the density peak is below this contrast")
        core_form.addRow("Min density contrast:", self.p_core_density_ratio)

        self.p_substar_nei_mult = QDoubleSpinBox()
        self.p_substar_nei_mult.setRange(2.0, 30.0)
        self.p_substar_nei_mult.setSingleStep(0.5)
        self.p_substar_nei_mult.setValue(_to_float(getattr(self.params.P, "psf_substar_neighbor_r_fwhm_mult", 8.0), 8.0))
        fit_form.addRow("Substar neighbor radius (×FWHM):", self.p_substar_nei_mult)

        self.p_substar_iters = QSpinBox()
        self.p_substar_iters.setRange(0, 2)
        self.p_substar_iters.setValue(_to_int(getattr(self.params.P, "psf_substar_iters", 1), 1))
        self.p_substar_iters.setToolTip("0 disables neighbour cleaning; 1 is the recommended CPU default")
        fit_form.addRow("Substar passes:", self.p_substar_iters)

        self.p_substar_max_src = QSpinBox()
        self.p_substar_max_src.setRange(0, 200000)
        self.p_substar_max_src.setSingleStep(100)
        self.p_substar_max_src.setValue(_to_int(getattr(self.params.P, "psf_substar_max_sources", 1500), 1500))
        self.p_substar_max_src.setToolTip("0이면 substar 이웃 소스 캡 무제한")
        fit_form.addRow("Substar max neighbor sources:", self.p_substar_max_src)

        self.p_conv_new = QDoubleSpinBox()
        self.p_conv_new.setRange(0.0, 1.0)
        self.p_conv_new.setSingleStep(0.005)
        self.p_conv_new.setValue(_to_float(getattr(self.params.P, "psf_conv_new_frac", 0.02), 0.02))
        self.p_conv_new.setToolTip("Uses unique candidates before the CPU source cap")
        redetect_form.addRow("Converge candidate frac <", self.p_conv_new)

        self.p_conv_flux = QDoubleSpinBox()
        self.p_conv_flux.setRange(0.0, 1.0)
        self.p_conv_flux.setSingleStep(0.001)
        self.p_conv_flux.setValue(_to_float(getattr(self.params.P, "psf_flux_conv_threshold", 0.01), 0.01))
        fit_form.addRow("Flux convergence fraction:", self.p_conv_flux)

        self.p_use_grouper = QCheckBox("Fit close neighbours together (CPU-limited)")
        self.p_use_grouper.setChecked(bool(getattr(self.params.P, "psf_use_grouper", False)))
        fit_form.addRow("", self.p_use_grouper)

        self.p_grouper_max_size = QSpinBox()
        self.p_grouper_max_size.setRange(1, 200)
        self.p_grouper_max_size.setValue(_to_int(getattr(self.params.P, "psf_grouper_max_size", 3), 3))
        self.p_grouper_max_size.setToolTip(
            "1 disables grouping; 2-3 is the CPU default. Groups above 4 use sparse LSQR."
        )
        fit_form.addRow("Group max size:", self.p_grouper_max_size)

        self.p_grouper_radius = QDoubleSpinBox()
        self.p_grouper_radius.setRange(0.5, 5.0)
        self.p_grouper_radius.setSingleStep(0.25)
        self.p_grouper_radius.setSuffix(" FWHM")
        self.p_grouper_radius.setValue(
            _to_float(getattr(self.params.P, "psf_grouper_radius_fwhm", 1.5), 1.5)
        )
        self.p_grouper_radius.setToolTip(
            "Neighbours inside this separation are fit simultaneously. Larger values cost more CPU."
        )
        fit_form.addRow("Group radius:", self.p_grouper_radius)

        # How much of a frame may be solved jointly. This was a literal 0.10 in
        # the fitting code, so only a tenth of each frame was ever grouped and
        # no setting could say otherwise — which made an ALLSTAR-style
        # comparison impossible to arrange (2026-08-14). 100 % is the
        # ALLSTAR-like end of the dial.
        self.p_grouper_budget_frac = QDoubleSpinBox()
        self.p_grouper_budget_frac.setRange(0.0, 100.0)
        self.p_grouper_budget_frac.setSingleStep(5.0)
        self.p_grouper_budget_frac.setSuffix(" %")
        self.p_grouper_budget_frac.setValue(
            _to_float(getattr(self.params.P, "psf_grouper_budget_frac", 0.10), 0.10) * 100.0
        )
        self.p_grouper_budget_frac.setToolTip(
            "Share of a frame's sources eligible for simultaneous fitting. "
            "100 % fits every crowded star jointly, as ALLSTAR does; the "
            "default 10 % keeps the cost of a dense field bounded."
        )
        fit_form.addRow("Group budget:", self.p_grouper_budget_frac)

        self.p_grouper_budget_cap = QSpinBox()
        self.p_grouper_budget_cap.setRange(0, 100000)
        self.p_grouper_budget_cap.setSpecialValueText("no cap")
        self.p_grouper_budget_cap.setValue(
            _to_int(getattr(self.params.P, "psf_grouper_budget_cap", 200), 200)
        )
        self.p_grouper_budget_cap.setToolTip(
            "Hard ceiling on grouped sources per frame, applied after the "
            "share above. 0 removes the ceiling."
        )
        fit_form.addRow("Group budget cap:", self.p_grouper_budget_cap)

        # DAOPHOT's `proferr`. Off by default and staying that way: measured
        # across three instruments on 2026-08-15 it bought 2-4 mmag in the most
        # blended bin on two of them and cost 8.5 mmag everywhere on the third,
        # where everything fainter than the bright anchor came back too bright.
        self.p_profile_error_frac = QDoubleSpinBox()
        self.p_profile_error_frac.setRange(0.0, 50.0)
        self.p_profile_error_frac.setSingleStep(0.5)
        self.p_profile_error_frac.setSuffix(" %")
        self.p_profile_error_frac.setSpecialValueText("off")
        self.p_profile_error_frac.setValue(
            _to_float(getattr(self.params.P, "psf_profile_error_frac", 0.0), 0.0) * 100.0
        )
        self.p_profile_error_frac.setToolTip(
            "Assumed PSF model error, added to the fit variance so a bright "
            "core stops dominating — DAOPHOT uses 5 %. Off by default: it "
            "helps blended stars on some instruments and reports faint stars "
            "too bright on others. Check the faint end before adopting it."
        )
        fit_form.addRow("Profile error:", self.p_profile_error_frac)

        self.p_forced_match_radius = QDoubleSpinBox()
        self.p_forced_match_radius.setRange(0.1, 3.0)
        self.p_forced_match_radius.setSingleStep(0.05)
        self.p_forced_match_radius.setSuffix(" FWHM")
        self.p_forced_match_radius.setValue(
            _to_float(
                getattr(self.params.P, "psf_forced_match_radius_fwhm", 1.25),
                1.25,
            )
        )
        self.p_forced_match_radius.setToolTip(
            "Step 4 detections inside this radius are anchored to their Step 7 catalog positions."
        )
        fit_form.addRow("Forced-catalog match:", self.p_forced_match_radius)

        self.p_use_error_img = QCheckBox("Use error image (slower, higher RAM)")
        self.p_use_error_img.setChecked(bool(getattr(self.params.P, "psf_use_error_image", True)))
        fit_form.addRow("", self.p_use_error_img)

        self.p_shared_filter_epsf = QCheckBox(
            "Share EPSF per filter (faster; disable if seeing varies >1px across frames)"
        )
        self.p_shared_filter_epsf.setChecked(
            bool(getattr(self.params.P, "psf_shared_filter_epsf", False))
        )
        epsf_form.addRow("", self.p_shared_filter_epsf)

        self.p_min_epsf_stars = QSpinBox()
        self.p_min_epsf_stars.setRange(1, 200)
        self.p_min_epsf_stars.setSingleStep(1)
        self.p_min_epsf_stars.setValue(_to_int(getattr(self.params.P, "psf_min_epsf_stars", 10), 10))
        self.p_min_epsf_stars.setToolTip(
            "Min isolated PSF stars required to build/cache a new EPSF.\n"
            "With 'Share EPSF' ON: frames below this threshold reuse the cached filter EPSF.\n"
            "Raise to avoid bad EPSF from crowded/trailed frames (e.g. 10–20)."
        )
        epsf_form.addRow("Min isolated PSF stars:", self.p_min_epsf_stars)

        self.p_sharp_lo = QDoubleSpinBox()
        self.p_sharp_lo.setRange(0.0, 1.0)
        self.p_sharp_lo.setSingleStep(0.05)
        self.p_sharp_lo.setDecimals(2)
        self.p_sharp_lo.setValue(_to_float(getattr(self.params.P, "psf_redetect_sharp_lo", 0.15), 0.15))
        redetect_form.addRow("Re-detect sharpness min:", self.p_sharp_lo)

        self.p_sharp_hi = QDoubleSpinBox()
        self.p_sharp_hi.setRange(0.0, 1.0)
        self.p_sharp_hi.setSingleStep(0.05)
        self.p_sharp_hi.setDecimals(2)
        self.p_sharp_hi.setValue(_to_float(getattr(self.params.P, "psf_redetect_sharp_hi", 0.95), 0.95))
        redetect_form.addRow("Re-detect sharpness max:", self.p_sharp_hi)

        self.p_round_max = QDoubleSpinBox()
        self.p_round_max.setRange(0.0, 2.0)
        self.p_round_max.setSingleStep(0.05)
        self.p_round_max.setDecimals(2)
        self.p_round_max.setValue(_to_float(getattr(self.params.P, "psf_redetect_round_abs_max", 0.8), 0.8))
        redetect_form.addRow("Re-detect |roundness| max:", self.p_round_max)

        self.p_save_residuals = QCheckBox("Save residual FITS (required for iter viewer)")
        self.p_save_residuals.setChecked(True)
        self.p_save_residuals.setEnabled(False)
        output_form.addRow("", self.p_save_residuals)

        self.p_save_all_iter_residuals = QCheckBox(
            "Also save background-restored star-subtracted images for every pass"
        )
        self.p_save_all_iter_residuals.setChecked(
            bool(getattr(self.params.P, "psf_save_all_iter_residuals", False))
        )
        output_form.addRow("", self.p_save_all_iter_residuals)

        # ── mode logic ────────────────────────────────────────────────────
        _manual_widgets = [
            self.p_n_stars, self.p_isolation, self.p_epsf_contamination_filter,
            self.p_flux_scale_correction, self.p_flux_scale_min_snr,
            self.p_flux_scale_min_stars, self.p_flux_scale_min_neighbor,
            self.p_flux_scale_max_scatter,
            self.p_fit_window_mode, self.p_fit_energy, self.p_fit_mult,
            self.p_max_iter,
            self.p_fitter_max_iter,
            self.p_redetect, self.p_dup_mult, self.p_dup_px,
            self.p_cap_per_iter, self.p_cap_frac, self.p_blend_ratio,
            self.p_postfit_snr, self.p_postfit_qfit, self.p_postfit_redchi,
            self.p_fit_init_max,
            self.p_core_enable, self.p_core_center_mode, self.p_core_x, self.p_core_y,
            self.p_core_radius_px, self.p_core_radius_mult, self.p_core_density_ratio,
            self.p_substar_iters, self.p_substar_nei_mult, self.p_substar_max_src,
            self.p_conv_new, self.p_conv_flux, self.p_use_grouper,
            self.p_grouper_max_size, self.p_grouper_radius,
            self.p_grouper_budget_frac, self.p_grouper_budget_cap,
            self.p_final_pass_max_iter, self.p_forced_position_lock,
            self.p_forced_match_radius,
            self.p_sharp_lo, self.p_sharp_hi, self.p_round_max,
        ]

        def _apply_mode_to_widgets(mode_key):
            p = _PSF_MODE_PRESETS.get(mode_key, _PSF_MODE_PRESETS["normal"])
            self.p_n_stars.setValue(p["psf_n_stars_max"])
            self.p_isolation.setValue(p["psf_isolation_fwhm_mult"])
            self.p_epsf_contamination_filter.setChecked(
                bool(p["psf_epsf_contamination_filter"])
            )
            self.p_flux_scale_correction.setChecked(
                bool(p["psf_flux_scale_correction"])
            )
            self.p_fit_window_mode.setCurrentIndex(
                max(0, self.p_fit_window_mode.findData(p["psf_fit_window_mode"]))
            )
            self.p_fit_energy.setValue(p["psf_fit_encircled_energy"])
            self.p_fit_mult.setValue(p["psf_fit_shape_fwhm_mult"])
            self.p_max_iter.setValue(p["psf_max_iter"])
            self.p_fitter_max_iter.setValue(p["psf_fitter_max_iter"])
            self.p_redetect.setValue(p["psf_redetect_sigma"])
            self.p_dup_mult.setValue(p["psf_duplicate_radius_fwhm_mult"])
            self.p_dup_px.setValue(0.0)
            self.p_cap_per_iter.setValue(p["psf_new_sources_cap_per_iter"])
            self.p_cap_frac.setValue(p["psf_new_sources_cap_frac"])
            self.p_blend_ratio.setValue(p["psf_blend_residual_ratio"])
            self.p_postfit_snr.setValue(p["psf_postfit_snr_min"])
            self.p_postfit_qfit.setValue(p["psf_postfit_qfit_max"])
            self.p_postfit_redchi.setValue(p["psf_postfit_reduced_chi2_max"])
            self.p_fit_init_max.setValue(p["psf_fit_init_max_sources"])
            self.p_core_enable.setChecked(bool(p["psf_core_cut_enable"]))
            self.p_core_center_mode.setCurrentIndex(max(0, self.p_core_center_mode.findData("auto")))
            self.p_core_x.setValue(0.0)
            self.p_core_y.setValue(0.0)
            self.p_core_radius_px.setValue(p["psf_core_cut_radius_px"])
            self.p_core_radius_mult.setValue(p["psf_core_cut_radius_fwhm_mult"])
            self.p_core_density_ratio.setValue(p["psf_core_cut_auto_min_density_ratio"])
            self.p_substar_iters.setValue(p["psf_substar_iters"])
            self.p_substar_nei_mult.setValue(p["psf_substar_neighbor_r_fwhm_mult"])
            self.p_substar_max_src.setValue(p["psf_substar_max_sources"])
            self.p_conv_new.setValue(p["psf_conv_new_frac"])
            self.p_conv_flux.setValue(p["psf_flux_conv_threshold"])
            self.p_use_grouper.setChecked(p["psf_use_grouper"])
            self.p_grouper_radius.setValue(p["psf_grouper_radius_fwhm"])
            self.p_grouper_budget_frac.setValue(p["psf_grouper_budget_frac"] * 100.0)
            self.p_grouper_budget_cap.setValue(p["psf_grouper_budget_cap"])
            self.p_profile_error_frac.setValue(p["psf_profile_error_frac"] * 100.0)
            self.p_final_pass_max_iter.setValue(p["psf_final_pass_max_iter"])
            _li = self.p_forced_position_lock.findData(p["psf_forced_position_lock"])
            self.p_forced_position_lock.setCurrentIndex(max(0, _li))
            self.p_forced_match_radius.setValue(p["psf_forced_match_radius_fwhm"])
            self.p_sharp_lo.setValue(p["psf_redetect_sharp_lo"])
            self.p_sharp_hi.setValue(p["psf_redetect_sharp_hi"])
            self.p_round_max.setValue(p["psf_redetect_round_abs_max"])

        _epsf_only_widgets = [
            self.p_oversampling, self.p_shared_filter_epsf, self.p_min_epsf_stars,
        ]

        def _refresh_controls():
            engine = self.p_fit_engine.currentData()
            is_moffat  = (self.p_build_mode.currentData() == "moffat")

            for w in _manual_widgets:
                w.setEnabled(True)
            for w in _epsf_only_widgets:
                w.setEnabled(not is_moffat)
            self.p_use_error_img.setEnabled(engine == "photutils")
            self.p_grouper_max_size.setEnabled(self.p_use_grouper.isChecked())
            self.p_grouper_radius.setEnabled(self.p_use_grouper.isChecked())
            self.p_grouper_budget_frac.setEnabled(self.p_use_grouper.isChecked())
            self.p_grouper_budget_cap.setEnabled(self.p_use_grouper.isChecked())
            scale_enabled = self.p_flux_scale_correction.isChecked()
            for widget in (
                self.p_flux_scale_min_snr,
                self.p_flux_scale_min_stars,
                self.p_flux_scale_min_neighbor,
                self.p_flux_scale_max_scatter,
            ):
                widget.setEnabled(scale_enabled)
            _sync_fit_window_controls()

        def _on_mode_changed():
            mode_key = mode_combo.currentData()
            if mode_key != "custom":
                _apply_mode_to_widgets(mode_key)
            _refresh_controls()

        mode_combo.currentIndexChanged.connect(lambda *_: _on_mode_changed())
        self.p_fit_engine.currentIndexChanged.connect(lambda *_: _refresh_controls())
        self.p_build_mode.currentIndexChanged.connect(lambda *_: _refresh_controls())
        self.p_use_grouper.toggled.connect(lambda *_: _refresh_controls())
        self.p_flux_scale_correction.toggled.connect(lambda *_: _refresh_controls())
        _refresh_controls()

        btns = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        def _after_psf_reset():
            _refresh_controls()

        add_parameter_reset_button(
            btns,
            [
                (mode_combo, "crowded"),
                (self.p_fit_engine, "apex_iterative"),
                (self.p_build_mode, "epsf"),
                (self.p_workers, 2),
                (self.p_oversampling, 2),
                (self.p_epsf_mult, 4.0),
                (self.p_n_stars, 0),
                (self.p_isolation, 2.0),
                (self.p_epsf_contamination_filter, True),
                (self.p_flux_scale_correction, False),
                (self.p_flux_scale_min_snr, 50.0),
                (self.p_flux_scale_min_stars, 8),
                (self.p_flux_scale_min_neighbor, 4.0),
                (self.p_flux_scale_max_scatter, 0.10),
                (self.p_fit_window_mode, "auto"),
                (self.p_fit_energy, 0.90),
                (self.p_fit_mult, 2.4),
                (self.p_max_iter, 2),
                (self.p_fitter_max_iter, 8),
                (self.p_redetect, 4.5),
                (self.p_redetect_g, 4.0),
                (self.p_redetect_r, 4.0),
                (self.p_redetect_i, 4.5),
                (self.p_dup_mult, 0.4),
                (self.p_dup_px, 0.0),
                (self.p_cap_per_iter, 50),
                (self.p_cap_frac, 0.01),
                (self.p_blend_ratio, 0.3),
                (self.p_postfit_snr, 3.0),
                (self.p_postfit_qfit, 3.0),
                (self.p_postfit_redchi, 25.0),
                (self.p_fit_init_max, 3000),
                (self.p_core_enable, False),
                (self.p_core_center_mode, "auto"),
                (self.p_core_x, 0.0),
                (self.p_core_y, 0.0),
                (self.p_core_radius_px, 0.0),
                (self.p_core_radius_mult, 20.0),
                (self.p_core_density_ratio, 1.5),
                (self.p_substar_iters, 1),
                (self.p_substar_nei_mult, 5.0),
                (self.p_substar_max_src, 1000),
                (self.p_conv_new, 0.02),
                (self.p_conv_flux, 0.01),
                (self.p_use_grouper, False),
                (self.p_grouper_max_size, 3),
                (self.p_grouper_radius, 1.5),
                (self.p_grouper_budget_frac, 10.0),
                (self.p_grouper_budget_cap, 200),
                (self.p_final_pass_max_iter, 2),
                (self.p_forced_match_radius, 1.25),
                (self.p_use_error_img, False),
                (self.p_shared_filter_epsf, False),
                (self.p_min_epsf_stars, 10),
                (self.p_sharp_lo, 0.2),
                (self.p_sharp_hi, 0.9),
                (self.p_round_max, 0.6),
                (self.p_save_residuals, True),
                (self.p_save_all_iter_residuals, False),
            ],
            on_reset=_after_psf_reset,
        )
        btns.accepted.connect(lambda: self._save_params(dialog, mode_combo.currentData()))
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)
        dialog.exec_()


    def _save_params(self, dialog, mode_key="normal"):
        self.params.P.psf_mode = mode_key
        self.params.P.psf_model_mode = "per_frame"
        self.params.P.psf_fit_engine = self.p_fit_engine.currentData()
        self.params.P.psf_build_mode = self.p_build_mode.currentData()
        self.params.P.psf_parallel_workers = self.p_workers.value()
        self.params.P.psf_epsf_oversampling = self.p_oversampling.value()
        self.params.P.psf_epsf_size_fwhm_mult = self.p_epsf_mult.value()
        self.params.P.psf_n_stars_max = self.p_n_stars.value()
        self.params.P.psf_isolation_fwhm_mult = self.p_isolation.value()
        self.params.P.psf_epsf_contamination_filter = (
            self.p_epsf_contamination_filter.isChecked()
        )
        self.params.P.psf_flux_scale_correction = self.p_flux_scale_correction.isChecked()
        self.params.P.psf_flux_scale_min_snr = self.p_flux_scale_min_snr.value()
        self.params.P.psf_flux_scale_min_stars = self.p_flux_scale_min_stars.value()
        self.params.P.psf_flux_scale_min_neighbor_fwhm = self.p_flux_scale_min_neighbor.value()
        self.params.P.psf_flux_scale_max_scatter_mag = self.p_flux_scale_max_scatter.value()
        self.params.P.psf_fit_window_mode = self.p_fit_window_mode.currentData()
        self.params.P.psf_fit_encircled_energy = self.p_fit_energy.value()
        self.params.P.psf_fit_shape_fwhm_mult = self.p_fit_mult.value()
        self.params.P.psf_max_iter = self.p_max_iter.value()
        self.params.P.psf_fitter_max_iter = self.p_fitter_max_iter.value()
        self.params.P.psf_redetect_sigma = self.p_redetect.value()
        def _spin_to_sigma(sp):
            v = sp.value()
            return float("nan") if v <= 0.0 else v
        self.params.P.psf_redetect_sigma_g = _spin_to_sigma(self.p_redetect_g)
        self.params.P.psf_redetect_sigma_r = _spin_to_sigma(self.p_redetect_r)
        self.params.P.psf_redetect_sigma_i = _spin_to_sigma(self.p_redetect_i)
        self.params.P.psf_duplicate_radius_fwhm_mult = self.p_dup_mult.value()
        self.params.P.psf_duplicate_radius_px = self.p_dup_px.value() if self.p_dup_px.value() > 0 else np.nan
        self.params.P.psf_new_sources_cap_per_iter = self.p_cap_per_iter.value()
        self.params.P.psf_new_sources_cap_frac = self.p_cap_frac.value()
        self.params.P.psf_blend_residual_ratio = self.p_blend_ratio.value()
        self.params.P.psf_postfit_snr_min = self.p_postfit_snr.value()
        self.params.P.psf_postfit_qfit_max = self.p_postfit_qfit.value()
        self.params.P.psf_postfit_reduced_chi2_max = self.p_postfit_redchi.value()
        self.params.P.psf_fit_init_max_sources = self.p_fit_init_max.value()
        self.params.P.psf_core_cut_enable = self.p_core_enable.isChecked()
        self.params.P.psf_core_cut_center_mode = self.p_core_center_mode.currentData()
        self.params.P.psf_core_cut_x_px = self.p_core_x.value()
        self.params.P.psf_core_cut_y_px = self.p_core_y.value()
        self.params.P.psf_core_cut_radius_px = self.p_core_radius_px.value()
        self.params.P.psf_core_cut_radius_fwhm_mult = self.p_core_radius_mult.value()
        self.params.P.psf_core_cut_auto_min_density_ratio = self.p_core_density_ratio.value()
        self.params.P.psf_substar_iters = self.p_substar_iters.value()
        self.params.P.psf_substar_neighbor_r_fwhm_mult = self.p_substar_nei_mult.value()
        self.params.P.psf_substar_max_sources = self.p_substar_max_src.value()
        self.params.P.psf_conv_new_frac = self.p_conv_new.value()
        self.params.P.psf_flux_conv_threshold = self.p_conv_flux.value()
        self.params.P.psf_use_grouper = self.p_use_grouper.isChecked()
        self.params.P.psf_grouper_max_size = self.p_grouper_max_size.value()
        self.params.P.psf_grouper_radius_fwhm = self.p_grouper_radius.value()
        self.params.P.psf_grouper_budget_frac = self.p_grouper_budget_frac.value() / 100.0
        self.params.P.psf_grouper_budget_cap = self.p_grouper_budget_cap.value()
        self.params.P.psf_profile_error_frac = self.p_profile_error_frac.value() / 100.0
        self.params.P.psf_final_pass_max_iter = self.p_final_pass_max_iter.value()
        self.params.P.psf_forced_position_lock = self.p_forced_position_lock.currentData()
        self.params.P.psf_forced_match_radius_fwhm = self.p_forced_match_radius.value()
        self.params.P.psf_use_error_image = self.p_use_error_img.isChecked()
        self.params.P.psf_shared_filter_epsf = self.p_shared_filter_epsf.isChecked()
        self.params.P.psf_min_epsf_stars = self.p_min_epsf_stars.value()
        self.params.P.psf_save_all_iter_residuals = self.p_save_all_iter_residuals.isChecked()
        self.params.P.psf_redetect_sharp_lo = self.p_sharp_lo.value()
        self.params.P.psf_redetect_sharp_hi = self.p_sharp_hi.value()
        self.params.P.psf_redetect_round_abs_max = self.p_round_max.value()
        self.params.P.psf_save_residuals = self.p_save_residuals.isChecked()
        self.save_state()
        self.persist_params()
        QMessageBox.information(dialog, "Saved", "Parameters saved.")
        dialog.accept()

    # ── Thread cleanup ────────────────────────────────────────────────────────

    def _cleanup_worker(self, timeout_ms=5000):
        if not self.worker:
            return
        if self.worker.isRunning():
            self.worker.stop()
            self.worker.quit()
            self.worker.wait(timeout_ms)
        try:
            self.worker.deleteLater()
        except Exception:
            pass
        self.worker = None

    # ── Log ───────────────────────────────────────────────────────────────────

    def log(self, message: str):
        append_timestamped_log(self.log_text, message)

    def show_log_window(self):
        show_raised(self.log_window)

    # ── Skip label ────────────────────────────────────────────────────────────

    def _update_skip_label(self):
        if not hasattr(self, "skip_label"):
            return
        if self._skip_psf:
            self.skip_label.setText(
                f"PSF SKIPPED — {self.downstream_name} will use Step 7 forced aperture results."
            )
        else:
            psf_idx = step8_psf_dir(self.params.P.result_dir) / "photometry_index.csv"
            if psf_idx.exists():
                self.skip_label.setText("PSF photometry results available.")
            else:
                self.skip_label.setText("")

    # ── Validation / State ────────────────────────────────────────────────────

    def validate_step(self) -> bool:
        """Step 8 is always valid: either PSF was run or it was skipped."""
        if self._skip_psf:
            return True
        valid, _ = self._current_psf_cache_status()
        return valid

    def save_state(self):
        self.project_state.store_step_data("psf_photometry", {
            "skip_psf": self._skip_psf,
            "use_existing_psf_output": (
                self.chk_use_existing_output.isChecked()
                if hasattr(self, "chk_use_existing_output")
                else True
            ),
            "psf_mode": getattr(self.params.P, "psf_mode", "normal"),
            "psf_model_mode": getattr(self.params.P, "psf_model_mode", "per_frame"),
            "psf_fit_engine": getattr(self.params.P, "psf_fit_engine", "apex_iterative"),
            "psf_build_mode": getattr(self.params.P, "psf_build_mode", "epsf"),
            "psf_parallel_workers": getattr(self.params.P, "psf_parallel_workers", 0),
            "psf_epsf_oversampling": getattr(self.params.P, "psf_epsf_oversampling", 2),
            "psf_epsf_size_px": getattr(self.params.P, "psf_epsf_size_px", 25),
            "psf_epsf_size_fwhm_mult": getattr(self.params.P, "psf_epsf_size_fwhm_mult", 4.0),
            "psf_n_stars_max": getattr(self.params.P, "psf_n_stars_max", 0),
            "psf_isolation_fwhm_mult": getattr(self.params.P, "psf_isolation_fwhm_mult", 3.0),
            "psf_epsf_contamination_filter": getattr(
                self.params.P,
                "psf_epsf_contamination_filter",
                True,
            ),
            "psf_flux_scale_correction": getattr(
                self.params.P, "psf_flux_scale_correction", False
            ),
            "psf_flux_scale_min_snr": getattr(
                self.params.P, "psf_flux_scale_min_snr", 50.0
            ),
            "psf_flux_scale_min_stars": getattr(
                self.params.P, "psf_flux_scale_min_stars", 8
            ),
            "psf_flux_scale_min_neighbor_fwhm": getattr(
                self.params.P, "psf_flux_scale_min_neighbor_fwhm", 4.0
            ),
            "psf_flux_scale_max_scatter_mag": getattr(
                self.params.P, "psf_flux_scale_max_scatter_mag", 0.10
            ),
            "psf_fit_shape_px": getattr(self.params.P, "psf_fit_shape_px", 5),
            "psf_fit_shape_fwhm_mult": getattr(self.params.P, "psf_fit_shape_fwhm_mult", 2.4),
            "psf_fit_window_mode": getattr(
                self.params.P, "psf_fit_window_mode", "auto"
            ),
            "psf_fit_encircled_energy": getattr(
                self.params.P, "psf_fit_encircled_energy", 0.90
            ),
            "psf_use_grouper": getattr(self.params.P, "psf_use_grouper", False),
            "psf_max_iter": getattr(self.params.P, "psf_max_iter", 2),
            "psf_fitter_max_iter": getattr(self.params.P, "psf_fitter_max_iter", 6),
            "psf_redetect_sigma": getattr(self.params.P, "psf_redetect_sigma", 4.0),
            "psf_redetect_sigma_g": getattr(self.params.P, "psf_redetect_sigma_g", float("nan")),
            "psf_redetect_sigma_r": getattr(self.params.P, "psf_redetect_sigma_r", float("nan")),
            "psf_redetect_sigma_i": getattr(self.params.P, "psf_redetect_sigma_i", float("nan")),
            "psf_duplicate_radius_fwhm_mult": getattr(self.params.P, "psf_duplicate_radius_fwhm_mult", 0.8),
            "psf_duplicate_radius_px": getattr(self.params.P, "psf_duplicate_radius_px", np.nan),
            "psf_new_sources_cap_per_iter": getattr(self.params.P, "psf_new_sources_cap_per_iter", 70),
            "psf_new_sources_cap_frac": getattr(self.params.P, "psf_new_sources_cap_frac", 0.02),
            "psf_blend_residual_ratio": getattr(self.params.P, "psf_blend_residual_ratio", 0.3),
            "psf_postfit_snr_min": getattr(self.params.P, "psf_postfit_snr_min", 3.0),
            "psf_postfit_qfit_max": getattr(self.params.P, "psf_postfit_qfit_max", 3.0),
            "psf_postfit_reduced_chi2_max": getattr(
                self.params.P, "psf_postfit_reduced_chi2_max", 25.0
            ),
            "psf_fit_init_max_sources": getattr(self.params.P, "psf_fit_init_max_sources", 0),
            "psf_core_cut_enable": getattr(self.params.P, "psf_core_cut_enable", False),
            "psf_core_cut_center_mode": getattr(self.params.P, "psf_core_cut_center_mode", "auto"),
            "psf_core_cut_x_px": getattr(self.params.P, "psf_core_cut_x_px", 0.0),
            "psf_core_cut_y_px": getattr(self.params.P, "psf_core_cut_y_px", 0.0),
            "psf_core_cut_radius_px": getattr(self.params.P, "psf_core_cut_radius_px", 0.0),
            "psf_core_cut_radius_fwhm_mult": getattr(self.params.P, "psf_core_cut_radius_fwhm_mult", 20.0),
            "psf_core_cut_auto_min_density_ratio": getattr(self.params.P, "psf_core_cut_auto_min_density_ratio", 1.5),
            "psf_substar_iters": getattr(self.params.P, "psf_substar_iters", 1),
            "psf_substar_neighbor_r_fwhm_mult": getattr(self.params.P, "psf_substar_neighbor_r_fwhm_mult", 8.0),
            "psf_substar_max_sources": getattr(self.params.P, "psf_substar_max_sources", 1500),
            "psf_conv_new_frac": getattr(self.params.P, "psf_conv_new_frac", 0.02),
            "psf_flux_conv_threshold": getattr(self.params.P, "psf_flux_conv_threshold", 0.01),
            "psf_use_error_image": getattr(self.params.P, "psf_use_error_image", False),
            "psf_shared_filter_epsf": getattr(self.params.P, "psf_shared_filter_epsf", False),
            "psf_grouper_max_size": getattr(self.params.P, "psf_grouper_max_size", 3),
            "psf_grouper_radius_fwhm": getattr(
                self.params.P, "psf_grouper_radius_fwhm", 1.5
            ),
            "psf_forced_match_radius_fwhm": getattr(
                self.params.P, "psf_forced_match_radius_fwhm", 1.25
            ),
            "psf_min_epsf_stars": getattr(self.params.P, "psf_min_epsf_stars", 10),
            "psf_save_all_iter_residuals": getattr(self.params.P, "psf_save_all_iter_residuals", False),
            "psf_redetect_sharp_lo": getattr(self.params.P, "psf_redetect_sharp_lo", 0.15),
            "psf_redetect_sharp_hi": getattr(self.params.P, "psf_redetect_sharp_hi", 0.95),
            "psf_redetect_round_abs_max": getattr(self.params.P, "psf_redetect_round_abs_max", 0.8),
            "psf_save_residuals": getattr(self.params.P, "psf_save_residuals", True),
        })

    def restore_state(self):
        state = self.project_state.get_step_data("psf_photometry")
        if state:
            self._skip_psf = bool(state.get("skip_psf", False))
            if hasattr(self, "chk_use_existing_output"):
                self.chk_use_existing_output.setChecked(
                    bool(state.get("use_existing_psf_output", True))
                )
            for k, v in state.items():
                if k not in {"skip_psf", "use_existing_psf_output"} and hasattr(self.params.P, k):
                    setattr(self.params.P, k, v)
        if str(getattr(self.params.P, "psf_model_mode", "per_frame")).strip().lower() != "per_frame":
            self.params.P.psf_model_mode = "per_frame"
        if str(getattr(self.params.P, "psf_fit_engine", "apex_iterative")).strip().lower() == "allstar":
            self.params.P.psf_fit_engine = "apex_iterative"
        _mode = str(getattr(self.params.P, "psf_mode", "normal")).strip().lower()
        _star_cap = _to_int(getattr(self.params.P, "psf_n_stars_max", 0), 0)
        if _mode != "custom" and _star_cap in {30, 40, 50}:
            self.params.P.psf_n_stars_max = 0
        self.params.P.psf_grouper_max_size = min(
            25,
            max(1, _to_int(getattr(self.params.P, "psf_grouper_max_size", 3), 3)),
        )
        self.params.P.psf_grouper_radius_fwhm = min(
            5.0,
            max(
                0.5,
                _to_float(getattr(self.params.P, "psf_grouper_radius_fwhm", 1.5), 1.5),
            ),
        )
        self.params.P.psf_forced_match_radius_fwhm = min(
            3.0,
            max(
                0.1,
                _to_float(
                    getattr(self.params.P, "psf_forced_match_radius_fwhm", 1.25),
                    1.25,
                ),
            ),
        )
        # Clamp fit_shape_fwhm_mult to a sensible minimum (< 1.0 is unusable).
        _fmult = _to_float(getattr(self.params.P, "psf_fit_shape_fwhm_mult", 1.5), 1.5)
        if _fmult < 1.0:
            self.params.P.psf_fit_shape_fwhm_mult = 1.5
        # Migrate broad defaults to tuned defaults unless user explicitly changed them.
        _rsig = _to_float(getattr(self.params.P, "psf_redetect_sigma", 4.0), 4.0)
        if abs(_rsig - 6.0) < 1e-6 or abs(_rsig - 7.5) < 1e-6:
            # 6.0 and 7.5 were old defaults; migrate to current default.
            self.params.P.psf_redetect_sigma = 4.0
        _cap_abs = _to_int(getattr(self.params.P, "psf_new_sources_cap_per_iter", 70), 70)
        if _cap_abs == 100:
            self.params.P.psf_new_sources_cap_per_iter = 70
        _cap_frac = _to_float(getattr(self.params.P, "psf_new_sources_cap_frac", 0.02), 0.02)
        if abs(_cap_frac - 0.04) < 1e-6:
            self.params.P.psf_new_sources_cap_frac = 0.02
        _slo = _to_float(getattr(self.params.P, "psf_redetect_sharp_lo", 0.15), 0.15)
        _shi = _to_float(getattr(self.params.P, "psf_redetect_sharp_hi", 0.95), 0.95)
        _rnd = _to_float(getattr(self.params.P, "psf_redetect_round_abs_max", 0.8), 0.8)
        if _slo <= -900.0 and _shi >= 900.0:
            self.params.P.psf_redetect_sharp_lo = 0.15
            self.params.P.psf_redetect_sharp_hi = 0.95
        if _rnd >= 9.0:
            self.params.P.psf_redetect_round_abs_max = 0.8
        # Migrate overly-loose UI values.
        if _slo <= 0.01 and _shi >= 0.99 and _rnd >= 1.5:
            self.params.P.psf_redetect_sharp_lo = 0.15
            self.params.P.psf_redetect_sharp_hi = 0.95
            self.params.P.psf_redetect_round_abs_max = 0.8
        self._update_skip_label()
        self.update_frame_table()
        idx_path = step8_psf_dir(self.params.P.result_dir) / "photometry_index.csv"
        if (not self._skip_psf) and idx_path.exists():
            valid, reason = self._current_psf_cache_status()
            if valid:
                try:
                    idx = pd.read_csv(idx_path)
                    n_frames = len(idx)
                    self.progress_bar.setMaximum(max(1, n_frames))
                    self.progress_bar.setValue(n_frames)
                    self.progress_label.setText(f"Loaded previous PSF output ({n_frames} frames)")
                    self.log(f"[PSF][CACHE] Loaded previous Step 8 PSF output from disk ({n_frames} frames).")
                    self._load_from_disk()
                except Exception as exc:
                    self.log(f"[PSF][CACHE] Previous output could not be loaded: {exc}")
            else:
                self.log(f"[PSF][CACHE] Previous output not restored ({reason}).")

    def closeEvent(self, event):
        if self.worker and self.worker.isRunning():
            self.stop_psf()
        self._cleanup_worker(timeout_ms=10000)
        super().closeEvent(event)
