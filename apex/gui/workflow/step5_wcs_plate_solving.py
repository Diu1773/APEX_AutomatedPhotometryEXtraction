"""
Step 5: WCS Plate Solving (ASTAP)
WCS solving window and workers.
"""

from __future__ import annotations

import json
import sys
import time
import subprocess
import threading
import tempfile
import warnings
import shlex
import shutil
from pathlib import Path, PureWindowsPath
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.coordinates import SkyCoord
import astropy.units as u

from scipy.spatial import cKDTree as KDTree

from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QGroupBox, QMessageBox,
    QTextEdit, QComboBox, QDialog, QFormLayout, QLineEdit, QDialogButtonBox,
    QProgressBar, QCheckBox, QSpinBox, QDoubleSpinBox, QTableWidget,
    QTableWidgetItem, QHeaderView, QAbstractItemView, QWidget, QTabWidget,
    QScrollArea, QFileDialog
)

from PyQt5.QtCore import Qt, QThread, pyqtSignal

from .step_window_base import StepWindowBase
from .run_control import RunControlBar, format_duration, progress_status_text
from .log_panel import WorkflowLogWindow, WorkerStatusPanel, append_timestamped_log, show_raised
from .ui_helpers import (
    add_parameter_reset_button,
    create_cache_checkbox,
    create_collapsible_section,
    create_parameter_button,
    configure_parameter_dialog,
    set_table_row_background,
    status_row_background,
)
from apex.core.cache_manager import StepCacheManager
from apex.utils.step_paths import (
    step2_cropped_dir,
    crop_is_active,
    crop_rect_path,
    step4_dir,
    step5_wcs_dir,
    step1_dir,
)
from apex.utils.constants import get_parallel_workers, MAD_TO_SIGMA
from apex.utils.gaia_catalog_service import (
    GaiaCatalogService,
    gaia_runtime_available,
)
from apex.utils.qc_utils import filter_files_by_qc, should_use_frame_quality_qc
from apex.utils.cache_utils import (
    norm_path_key,
    build_file_signature,
    detection_cache_signature_matches,
    file_signature_matches_relaxed,
    astap_wcs_candidates,
    parse_astap_wcs_file,
)


# ── WCS solving workers (Qt delegators) ──────────────────────────────────────
# The compute lives Qt-free in apex.analysis.wcs_solve (relocated VERBATIM from
# the previous in-file workers). Here we keep thin QThread subclasses that add
# the real pyqtSignals and re-emit the relocated bodies' shim-signal events, so
# the GUI keeps its cross-thread (queued-connection) signal delivery and the
# step window code is unchanged. Module-level solver helpers (preflight checks,
# QC metrics, etc.) are re-exported from the analysis module so the rest of this
# file (e.g. _log_preflight) keeps using them by the same names.
from apex.analysis.wcs_solve import (
    _SUBPROCESS_TEXT_KWARGS,
    _tail_text,
    _exc_brief,
    _is_explicit_stopped_message,
    _coord_sep_deg,
    _format_coord_hint,
    _header_pointing_coord,
    _strip_outer_quotes,
    _split_command,
    _local_executable_available,
    _APEX_WCS_META_KEYS,
    _load_simbad_target_coord,
    _reset_apex_wcs_meta,
    _shared_wcs_center_coords,
    _shared_empty_wcs_qc_metrics,
    _shared_compute_wcs_qc_metrics,
    _shared_evaluate_wcs_qc_pass,
    _check_astap_available,
    _classify_astap_failure,
    _check_astnet_available,
    _check_gaia_runtime_available,
    _astnet_center_candidates,
    _wsl_path_exists_probe,
    _wsl_ensure_writable_dir_probe,
    _preferred_astnet_solution_path,
    _astnet_solution_artifacts_ready,
    _cleanup_redundant_astnet_new,
    _solution_header_shape,
    WcsWorkerBase,
    AstrometryNetWorkerBase,
    InternalWcsWorkerBase,
)


class WcsWorker(WcsWorkerBase, QThread):
    """QThread wrapper around the Qt-free WcsWorkerBase (ASTAP solving)."""
    progress = pyqtSignal(int, int, str)
    file_done = pyqtSignal(str, dict)
    worker_status = pyqtSignal(int, str, str, int)
    log_message = pyqtSignal(str)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str, str)

    def _init_signals(self):
        # Real pyqtSignals are class attributes; nothing to set per-instance.
        pass

    def __init__(self, file_list, params, data_dir, result_dir, cache_dir,
                 use_cropped=False, target_coord=None):
        QThread.__init__(self)
        WcsWorkerBase.__init__(
            self, file_list, params, data_dir, result_dir, cache_dir,
            use_cropped=use_cropped, target_coord=target_coord,
        )


class AstrometryNetWorker(AstrometryNetWorkerBase, QThread):
    """QThread wrapper around the Qt-free AstrometryNetWorkerBase (solve-field)."""
    progress = pyqtSignal(int, int, str)
    file_done = pyqtSignal(str, dict)
    refine_done = pyqtSignal(str, dict)
    log_message = pyqtSignal(str)
    worker_status = pyqtSignal(int, str, str, int)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str, str)

    def _init_signals(self):
        pass

    def __init__(self, file_list, params, data_dir, result_dir, cache_dir,
                 use_cropped=False, target_coord=None):
        QThread.__init__(self)
        AstrometryNetWorkerBase.__init__(
            self, file_list, params, data_dir, result_dir, cache_dir,
            use_cropped=use_cropped, target_coord=target_coord,
        )


class InternalWcsWorker(InternalWcsWorkerBase, QThread):
    """QThread wrapper around the Qt-free InternalWcsWorkerBase (Python engine)."""
    progress      = pyqtSignal(int, int, str)
    log_message   = pyqtSignal(str)
    frame_done    = pyqtSignal(str, dict)
    finished      = pyqtSignal(dict)
    error         = pyqtSignal(str)
    worker_status = pyqtSignal(int, str, str, int)

    def _init_signals(self):
        pass

    def __init__(self, file_list, params, data_dir, result_dir,
                 use_cropped=False,
                 min_matches=8,
                 rms_max_arcsec=2.0,
                 advanced_params=None):
        QThread.__init__(self)
        InternalWcsWorkerBase.__init__(
            self, file_list, params, data_dir, result_dir,
            use_cropped=use_cropped,
            min_matches=min_matches,
            rms_max_arcsec=rms_max_arcsec,
            advanced_params=advanced_params,
        )



class WcsPlateSolvingWindow(StepWindowBase):
    """Step 5: WCS Plate Solving"""

    def __init__(self, params, file_manager, project_state, main_window):
        self.file_manager = file_manager
        self.worker = None
        self.results = {}
        self.stop_requested = False
        self.log_window = None
        self._wcs_cache_mgr: StepCacheManager | None = None

        self.file_list = []
        self.use_cropped = False

        super().__init__(
            step_index=4,
            step_name="WCS Plate Solving",
            params=params,
            project_state=project_state,
            main_window=main_window
        )

        self.setup_step_ui()
        self.restore_state()

    def setup_step_ui(self):
        # ── Persistent "which solver ran each file?" banner ─────────────────
        # The three tabs each only render their own solver's runs, so after a
        # project re-open it is otherwise invisible which engine produced the
        # current WCS for each frame. This banner reads self.results (loaded
        # from wcs_solve_summary.csv) and prints a per-solver tally above the
        # tabs so the answer is one glance away, regardless of which tab is
        # active.
        self.solver_breakdown_label = QLabel("")
        self.solver_breakdown_label.setWordWrap(True)
        self.solver_breakdown_label.setStyleSheet(
            "QLabel { background-color: #ECEFF1; padding: 6px 10px;"
            " border-radius: 4px; color: #37474F; }"
        )
        self.content_layout.addWidget(self.solver_breakdown_label)

        # Create tab widget
        self.tab_widget = QTabWidget()
        self.content_layout.addWidget(self.tab_widget)

        # Internal (Python) Tab — first, no external tools required
        self.internal_tab = QWidget()
        self.setup_internal_tab()
        self.tab_widget.addTab(self.internal_tab, "Internal (Python)")

        # ASTAP Tab
        self.astap_tab = QWidget()
        self.setup_astap_tab()
        self.tab_widget.addTab(self.astap_tab, "ASTAP (Local)")

        # Astrometry.net Tab
        self.astrometrynet_tab = QWidget()
        self.setup_astrometrynet_tab()
        self.tab_widget.addTab(self.astrometrynet_tab, "Astrometry.net (Local)")

        self.setup_log_window()
        self.populate_file_list()

    def setup_internal_tab(self):
        """Internal (Python) WCS solver — no external executables needed."""
        layout = QVBoxLayout(self.internal_tab)

        info = QLabel(
            "Solve WCS using the built-in astnet-style quad-hash matcher and "
            "Gaia refinement (numpy/scipy/astropy only — no external binary). "
            "Needs Step 1 target RA/Dec and the Instrument pixel scale."
        )
        info.setWordWrap(True)
        info.setStyleSheet("QLabel { background-color: #E8F5E9; padding: 10px; border-radius: 5px; }")
        layout.addWidget(info)

        # All Internal-solver parameters live in one dict edited by the
        # parameter dialog — same pattern as ASTAP/astnet so the three tabs
        # behave consistently.
        self._internal_params = {
            "min_matches": 8,
            "rms_max_arcsec": 2.0,
            "n_brightest_src": 350,
            "n_brightest_cat": 900,
            "quad_k_neighbor": 10,
            "quad_neighbor_pool_factor": 3,
            "quad_code_tol": 0.020,
            "quad_scale_ratio_tol": 0.25,
            "quad_max_per_side": 3000,
            "ransac_inlier_radius_px": 4.0,
            "ransac_max_trials": 4000,
            "ransac_keep_candidates": 8,
            "allow_reflection": True,
            "local_blind_fallback": True,
            "local_blind_radius_factor": 2.5,
            "sip_degree": 3,
            "sip_min_pairs": 30,
            "sip_holdout_fraction": 0.25,
            "sip_min_improvement": 0.10,
        }

        control_layout = QHBoxLayout()
        btn_int_params = create_parameter_button("Internal Parameters")
        btn_int_params.clicked.connect(self._open_internal_parameters_dialog)
        control_layout.addWidget(btn_int_params)

        self.run_bar_internal = RunControlBar(
            "Run Internal Solver", "Log",
            run_cb=self.run_wcs_internal_solver,
            stop_cb=self.stop_wcs_internal_solver,
            log_cb=self.show_log_window,
        )
        control_layout.addWidget(self.run_bar_internal)
        layout.addLayout(control_layout)

        progress_layout = QHBoxLayout()
        self.internal_progress = QProgressBar()
        self.internal_progress.setMinimum(0)
        self.internal_progress.setMaximum(100)
        self.internal_progress.setValue(0)
        progress_layout.addWidget(self.internal_progress)

        self.internal_status = QLabel("Ready")
        self.internal_status.setMinimumWidth(350)
        progress_layout.addWidget(self.internal_status)
        layout.addLayout(progress_layout)

        results_group = QGroupBox("Internal Solver Results")
        results_layout = QVBoxLayout(results_group)

        self.internal_results_table = QTableWidget()
        self.internal_results_table.setColumnCount(9)
        self.internal_results_table.setHorizontalHeaderLabels([
            "File", "Status", "Hint", "RA", "Dec", "PixScale", "Refine", "Resid(\")", "Elapsed (s)"
        ])
        self.internal_results_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.internal_results_table.horizontalHeader().setStretchLastSection(True)
        self.internal_results_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.internal_results_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        results_layout.addWidget(self.internal_results_table)

        layout.addWidget(results_group)

    # ── Shared run front-matter (all three solver tabs) ──────────────────
    def _qc_filter_and_log(self, file_list=None):
        """Apply the Step 4 QC frame filter and log the outcome.

        Shared by the ASTAP / astrometry.net / Internal tabs so the QC line is
        computed and worded identically everywhere. Returns the filtered file
        list, or ``None`` when nothing remains (the caller should abort — a
        warning dialog has already been shown).
        """
        src = list(self.file_list if file_list is None else file_list)
        require_qc = should_use_frame_quality_qc(
            Path(self.params.P.result_dir),
            self.params.P,
            "wcs_require_qc_pass",
            default=True,
        )
        filtered, qc_info = filter_files_by_qc(
            Path(self.params.P.result_dir), src, require_qc=require_qc,
        )
        if require_qc:
            if qc_info.get("applied"):
                self.log(f"[WCS][QC] Frame QC filter: {qc_info['kept']}/{qc_info['total']} kept.")
            elif qc_info.get("path") is None:
                self.log("[WCS][QC] frame_quality.csv not found; using all frames.")
            else:
                self.log(f"[WCS][QC] frame_quality.csv ignored ({qc_info['reason']}); using all frames.")
        if not filtered:
            QMessageBox.warning(self, "Warning", "No frames remain after Step 4 QC filtering.")
            return None
        return filtered

    def _log_preflight(self, *, astap=False, astnet=False, gaia=False):
        """Log a uniform ``[WCS][Preflight]`` header and return availability.

        Each tab requests only the backends it actually uses, so no slow probe
        runs unnecessarily (e.g. the Internal tab never spins up the WSL
        astrometry.net check). The wording/format is shared so the preflight
        header reads the same regardless of which tab launched the run.

        Returns a dict like ``{"astap": (ok, detail), ...}`` for requested
        backends only.
        """
        out: dict[str, tuple[bool, str]] = {}
        if astap:
            ok, detail = _check_astap_available(self.params)
            self.log(f"[WCS][Preflight] ASTAP: {detail}")
            out["astap"] = (ok, detail)
        if astnet:
            ok, detail = _check_astnet_available(self.params)
            self.log(f"[WCS][Preflight] astrometry.net: {detail}")
            out["astnet"] = (ok, detail)
        if gaia:
            ok, detail = _check_gaia_runtime_available()
            self.log(f"[WCS][Preflight] Gaia: {detail}")
            out["gaia"] = (ok, detail)
        return out

    def run_wcs_internal_solver(self):
        if not self.file_list:
            QMessageBox.warning(self, "Warning", "No files to process")
            return
        if hasattr(self, "_internal_worker") and self._internal_worker and \
                self._internal_worker.isRunning():
            return

        self.log_text.clear()
        file_list = self._qc_filter_and_log()
        if file_list is None:
            return

        # Internal can still solve from a cached gaia_fov.ecsv even when the
        # runtime probe is unavailable, so this is informational — no early
        # return.
        self._log_preflight(gaia=True)

        self.internal_results_table.setRowCount(0)
        # Drop any prior run's rows so a fresh Internal run starts from a
        # clean slate. ASTAP / astnet handlers reset implicitly through their
        # own paths; we do it explicitly here so cache writers below see only
        # the current solver's output and never accidentally pick up stale
        # ASTAP/astnet rows that share filenames.
        self.results = {}
        self._internal_run_results = {}

        p = dict(self._internal_params)
        min_matches = int(p.pop("min_matches"))
        rms_max = float(p.pop("rms_max_arcsec"))

        use_cropped = bool(getattr(self, "use_cropped", False))
        self._internal_worker = InternalWcsWorker(
            file_list,
            self.params,
            self.params.P.data_dir,
            self.params.P.result_dir,
            use_cropped=use_cropped,
            min_matches=min_matches,
            rms_max_arcsec=rms_max,
            advanced_params=p,
        )
        self._internal_worker.progress.connect(self._on_internal_progress)
        self._internal_worker.log_message.connect(self.log)
        self._internal_worker.frame_done.connect(self._on_internal_frame_done)
        self._internal_worker.finished.connect(self._on_internal_finished)
        self._internal_worker.error.connect(self._on_internal_error)

        # Wire the shared Worker Status panel so the Internal tab shows one row
        # per parallel worker, just like the ASTAP / astnet tabs.
        self.setup_log_window()
        if hasattr(self, "_worker_panel") and self._worker_panel is not None:
            self._worker_panel.clear()
            self._internal_worker.worker_status.connect(self._worker_panel.update_worker)

        self.run_bar_internal.set_running(True)
        self.internal_progress.setValue(0)
        self.internal_progress.setMaximum(len(file_list))
        self.internal_status.setText("Solving…")
        self._internal_worker.start()

    def _on_internal_progress(self, i: int, total: int, fname: str):
        self.internal_progress.setMaximum(max(total, 1))
        self.internal_progress.setValue(int(i))
        self.internal_status.setText(f"{i}/{total}  {fname}")

    def stop_wcs_internal_solver(self):
        if hasattr(self, "_internal_worker") and self._internal_worker:
            self._internal_worker.stop()

    def _open_internal_parameters_dialog(self):
        """Edit Internal Solver parameters in a single dialog.

        Same styling pattern as ASTAP / Astrometry.net parameter dialogs
        (configure_parameter_dialog + grouped form).
        """
        from PyQt5.QtWidgets import (
            QDialog, QFormLayout, QVBoxLayout, QSpinBox, QDoubleSpinBox,
            QCheckBox, QDialogButtonBox, QGroupBox,
        )

        d = FittedDialog(self)
        configure_parameter_dialog(d, "Internal Parameters", 560, 780)
        root = QVBoxLayout(d)
        p = self._internal_params

        basic = QGroupBox("Basic")
        bf = QFormLayout(basic)
        sp_min = QSpinBox(); sp_min.setRange(4, 100); sp_min.setValue(int(p["min_matches"]))
        sp_min.setToolTip("Reject solutions with fewer than this many matched stars")
        bf.addRow("min_matches:", sp_min)
        sp_rms = QDoubleSpinBox(); sp_rms.setRange(0.1, 10.0); sp_rms.setDecimals(2)
        sp_rms.setSingleStep(0.1); sp_rms.setSuffix(" \"")
        sp_rms.setValue(float(p["rms_max_arcsec"]))
        sp_rms.setToolTip("Reject solutions with residual RMS above this")
        bf.addRow("rms_max:", sp_rms)
        root.addWidget(basic)

        quad = QGroupBox("Quad matching")
        qf = QFormLayout(quad)
        sp_src = QSpinBox(); sp_src.setRange(20, 2000); sp_src.setValue(int(p["n_brightest_src"]))
        sp_src.setToolTip("Brightest detected sources fed into quad construction")
        qf.addRow("n_brightest_src:", sp_src)
        sp_cat = QSpinBox(); sp_cat.setRange(20, 5000); sp_cat.setValue(int(p["n_brightest_cat"]))
        sp_cat.setToolTip("Brightest in-field Gaia stars fed into quad construction. "
                          "Crowded fields may need ≥500 so that halo stars enter the pool.")
        qf.addRow("n_brightest_cat:", sp_cat)
        sp_k = QSpinBox(); sp_k.setRange(4, 30); sp_k.setValue(int(p["quad_k_neighbor"]))
        sp_k.setToolTip("Neighbour sample size per star for quad construction")
        qf.addRow("quad_k_neighbor:", sp_k)
        sp_pool = QSpinBox(); sp_pool.setRange(1, 8); sp_pool.setValue(int(p["quad_neighbor_pool_factor"]))
        sp_pool.setToolTip("Sample neighbours from a wider pool; higher values add long baselines in crowded fields")
        qf.addRow("quad_neighbor_pool_factor:", sp_pool)
        sp_tol = QDoubleSpinBox(); sp_tol.setRange(0.001, 0.20); sp_tol.setDecimals(4)
        sp_tol.setSingleStep(0.005); sp_tol.setValue(float(p["quad_code_tol"]))
        sp_tol.setToolTip("L2 tolerance in 4-D code space for two quads to be a match")
        qf.addRow("quad_code_tol:", sp_tol)
        sp_sr = QDoubleSpinBox(); sp_sr.setRange(0.02, 1.0); sp_sr.setDecimals(3)
        sp_sr.setSingleStep(0.05); sp_sr.setValue(float(p["quad_scale_ratio_tol"]))
        sp_sr.setToolTip("Allowed |log(src_side / cat_side)| between matched quads")
        qf.addRow("quad_scale_ratio_tol:", sp_sr)
        sp_mq = QSpinBox(); sp_mq.setRange(100, 10000); sp_mq.setValue(int(p["quad_max_per_side"]))
        sp_mq.setSingleStep(100)
        sp_mq.setToolTip("Cap on quads kept per side (larger = slower but more thorough)")
        qf.addRow("quad_max_per_side:", sp_mq)
        root.addWidget(quad)

        ransac = QGroupBox("RANSAC verification")
        rf = QFormLayout(ransac)
        sp_ir = QDoubleSpinBox(); sp_ir.setRange(0.5, 20.0); sp_ir.setDecimals(2)
        sp_ir.setSingleStep(0.5); sp_ir.setValue(float(p["ransac_inlier_radius_px"]))
        sp_ir.setSuffix(" px"); sp_ir.setToolTip("Inlier radius for the candidate similarity")
        rf.addRow("ransac_inlier_radius_px:", sp_ir)
        sp_mt = QSpinBox(); sp_mt.setRange(50, 20000); sp_mt.setValue(int(p["ransac_max_trials"]))
        sp_mt.setSingleStep(100)
        sp_mt.setToolTip("Maximum quad-pair candidates to test (runtime cap)")
        rf.addRow("ransac_max_trials:", sp_mt)
        sp_keep = QSpinBox(); sp_keep.setRange(1, 20); sp_keep.setValue(int(p["ransac_keep_candidates"]))
        sp_keep.setToolTip("Distinct RANSAC candidates carried into final WCS fitting")
        rf.addRow("ransac_keep_candidates:", sp_keep)
        cb_refl = QCheckBox("Allow mirrored transform (image/sky parity flip)")
        cb_refl.setChecked(bool(p["allow_reflection"]))
        rf.addRow(cb_refl)
        root.addWidget(ransac)

        fallback = QGroupBox("Fallback")
        lf = QFormLayout(fallback)
        cb_lb = QCheckBox("Retry local-blind when normal solve fails")
        cb_lb.setChecked(bool(p.get("local_blind_fallback", True)))
        cb_lb.setToolTip(
            "Retry once with a wider Gaia tangent-plane window. This is not "
            "all-sky blind solving; Step 1 target RA/Dec and pixel scale are "
            "still required."
        )
        lf.addRow(cb_lb)
        sp_lb = QDoubleSpinBox(); sp_lb.setRange(1.0, 8.0); sp_lb.setDecimals(2)
        sp_lb.setSingleStep(0.25); sp_lb.setValue(float(p.get("local_blind_radius_factor", 2.5)))
        sp_lb.setToolTip("Gaia window radius in half-image-diagonal units for local-blind retry")
        lf.addRow("local_blind_radius_factor:", sp_lb)
        root.addWidget(fallback)

        refine = QGroupBox("Refine")
        ff = QFormLayout(refine)
        sp_sip = QSpinBox(); sp_sip.setRange(0, 5); sp_sip.setValue(int(p["sip_degree"]))
        sp_sip.setToolTip("SIP polynomial degree for distortion correction. "
                          "0 = pure TAN. 2–3 typical. Needs ≥30 matched pairs to engage.")
        ff.addRow("sip_degree:", sp_sip)
        sp_sip_min = QSpinBox(); sp_sip_min.setRange(12, 300)
        sp_sip_min.setValue(int(p.get("sip_min_pairs", 30)))
        sp_sip_min.setToolTip("Minimum matched pairs required before SIP validation is attempted")
        ff.addRow("sip_min_pairs:", sp_sip_min)
        sp_sip_hold = QDoubleSpinBox(); sp_sip_hold.setRange(0.05, 0.50)
        sp_sip_hold.setDecimals(2); sp_sip_hold.setSingleStep(0.05)
        sp_sip_hold.setValue(float(p.get("sip_holdout_fraction", 0.25)))
        sp_sip_hold.setToolTip("Fraction of matched pairs held out for SIP-vs-TAN validation")
        ff.addRow("sip_holdout_fraction:", sp_sip_hold)
        sp_sip_imp = QDoubleSpinBox(); sp_sip_imp.setRange(0.00, 0.50)
        sp_sip_imp.setDecimals(2); sp_sip_imp.setSingleStep(0.05)
        sp_sip_imp.setValue(float(p.get("sip_min_improvement", 0.10)))
        sp_sip_imp.setToolTip("Minimum holdout RMS improvement required to adopt SIP")
        ff.addRow("sip_min_improvement:", sp_sip_imp)
        root.addWidget(refine)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        add_parameter_reset_button(
            bb,
            [
                (sp_min, 8),
                (sp_rms, 2.0),
                (sp_src, 350),
                (sp_cat, 900),
                (sp_k, 10),
                (sp_pool, 3),
                (sp_tol, 0.020),
                (sp_sr, 0.25),
                (sp_mq, 3000),
                (sp_ir, 4.0),
                (sp_mt, 4000),
                (sp_keep, 8),
                (cb_refl, True),
                (cb_lb, True),
                (sp_lb, 2.5),
                (sp_sip, 3),
                (sp_sip_min, 30),
                (sp_sip_hold, 0.25),
                (sp_sip_imp, 0.10),
            ],
        )
        bb.accepted.connect(d.accept)
        bb.rejected.connect(d.reject)
        root.addWidget(bb)

        if d.exec_() != QDialog.Accepted:
            return

        self._internal_params.update({
            "min_matches": int(sp_min.value()),
            "rms_max_arcsec": float(sp_rms.value()),
            "n_brightest_src": int(sp_src.value()),
            "n_brightest_cat": int(sp_cat.value()),
            "quad_k_neighbor": int(sp_k.value()),
            "quad_neighbor_pool_factor": int(sp_pool.value()),
            "quad_code_tol": float(sp_tol.value()),
            "quad_scale_ratio_tol": float(sp_sr.value()),
            "quad_max_per_side": int(sp_mq.value()),
            "ransac_inlier_radius_px": float(sp_ir.value()),
            "ransac_max_trials": int(sp_mt.value()),
            "ransac_keep_candidates": int(sp_keep.value()),
            "allow_reflection": bool(cb_refl.isChecked()),
            "local_blind_fallback": bool(cb_lb.isChecked()),
            "local_blind_radius_factor": float(sp_lb.value()),
            "sip_degree": int(sp_sip.value()),
            "sip_min_pairs": int(sp_sip_min.value()),
            "sip_holdout_fraction": float(sp_sip_hold.value()),
            "sip_min_improvement": float(sp_sip_imp.value()),
        })

    def _on_internal_frame_done(self, fname: str, info: dict):
        """Add one row to the results table as each frame completes.

        Column layout mirrors the Astrometry.net tab so the three solver tabs
        present comparable data side-by-side. Also records the result into
        ``self.results`` so that ``validate_step()`` recognises the step as
        complete and the Next Step button enables — the ASTAP and astnet
        worker callbacks already do the same.
        """
        from PyQt5.QtWidgets import QTableWidgetItem

        if not hasattr(self, "_internal_run_results"):
            self._internal_run_results = {}
        self._internal_run_results[fname] = info
        if bool(info.get("ok", False)):
            self.results[fname] = info
        else:
            self.results.pop(fname, None)
        # Cache writes (parity with ASTAP/astnet on_file_done):
        #   * StepCacheManager manifest per successful frame
        #   * Navigation refresh so the Next Step button enables as soon as
        #     at least one frame is solved.
        # Both are wrapped in broad try/except — a cache hiccup must never
        # stop the per-frame loop from displaying results.
        try:
            if bool(info.get("ok", False)):
                self._write_wcs_manifest(fname, info)
        except Exception:
            pass
        try:
            self.update_navigation_buttons()
        except Exception:
            pass
        try:
            self._refresh_solver_breakdown_label()
        except Exception:
            pass

        t = self.internal_results_table
        row = t.rowCount()
        t.insertRow(row)
        ok = bool(info.get("ok", False))

        def _fmt_num(key, fmt):
            v = info.get(key)
            if v is None or not np.isfinite(float(v)):
                return ""
            return fmt.format(float(v))

        if ok:
            ra_s = _fmt_num("center_ra_deg", "{:.4f}")
            dec_s = _fmt_num("center_dec_deg", "{:.4f}")
            pix_s = _fmt_num("scale_arcsec_per_px", "{:.4f}\"/px")
            model_s = str(info.get("model", "") or "").strip()
            refine_s = f"{model_s} m={int(info.get('n_matches', 0) or 0)}".strip()
            resid_s = _fmt_num("rms_arcsec", "{:.3f}")
            status_s = "✓ QC" if self._boolish(info.get("wcs_qc_pass")) else "✓ QC?"
        else:
            ra_s = dec_s = pix_s = refine_s = resid_s = ""
            status_s = "✗ " + str(info.get("reason", "fail"))

        elapsed_s = _fmt_num("elapsed_s", "{:.2f}")
        hint_s = str(info.get("hint_source", "") or "")

        values = [fname, status_s, hint_s, ra_s, dec_s, pix_s, refine_s, resid_s, elapsed_s]
        for col, val in enumerate(values):
            t.setItem(row, col, QTableWidgetItem(val))
        set_table_row_background(t, row, self._wcs_result_row_background(info))
        t.scrollToBottom()

    def _on_internal_finished(self, summary: dict):
        self.run_bar_internal.set_running(False)
        ok = int(summary.get("ok", 0) or 0)
        total = int(summary.get("total", 0) or 0)
        stopped = bool(summary.get("stopped", False))
        msg = f"Internal solver: {ok}/{total} solved"
        if stopped:
            msg += " (stopped)"
        self.internal_status.setText(msg)
        qc_not_eval = int(summary.get("wcs_qc_not_evaluated", 0) or 0)
        qc_text = (
            f"WCS-QC not evaluated: {qc_not_eval} (Gaia unavailable)"
            if qc_not_eval
            else f"WCS-QC pass: {summary.get('wcs_qc_pass', 0)}"
        )
        self.log(f"[Internal WCS] {msg} | {qc_text}")

        # Final cache layer: write the solver-agnostic summary CSV. Without
        # this, _restore_success_results_from_summary() on the next project
        # open finds no rows for Internal-solved frames, the Next Step button
        # stays disabled, and downstream QC tools that read this CSV miss the
        # frames entirely. Wrapped so a write failure is logged but does not
        # break the user-visible "N/N solved" completion.
        try:
            self._write_internal_summary_csv()
        except Exception as exc:
            self.log(f"[Internal WCS] summary CSV write failed: {exc}")

    def _refresh_solver_breakdown_label(self) -> None:
        """Re-tally ``self.results`` by WCS solver and update the banner.

        Source of truth for the per-frame tag is, in order:
          1. ``info["wcssrc"]``  (the FITS-header keyword written by the
             worker — survives across project re-opens because it's also in
             the summary CSV)
          2. ``info["solver"]`` (legacy / Internal-native key)

        Frames whose result row doesn't carry either key are bucketed under
        ``"unknown"`` rather than silently dropped, so a missing tag is
        visible instead of pretending everything was solved.
        """
        label = getattr(self, "solver_breakdown_label", None)
        if label is None:
            return
        results = getattr(self, "results", {}) or {}
        if not results:
            label.setText("No WCS results yet — run a solver tab below.")
            return
        counts: dict[str, int] = {}
        for info in results.values():
            if not isinstance(info, dict):
                continue
            tag = info.get("wcssrc") or info.get("solver") or "unknown"
            tag = str(tag).strip() or "unknown"
            counts[tag] = counts.get(tag, 0) + 1
        total = sum(counts.values())
        parts = [f"<b>{tag}</b>: {n}" for tag, n in sorted(counts.items())]
        label.setText(
            f"Solved frames ({total}) — " + ",&nbsp; ".join(parts)
        )

    @staticmethod
    def _finite_result_float(value, default=np.nan) -> float:
        try:
            out = float(value)
        except Exception:
            return float(default)
        return out if np.isfinite(out) else float(default)

    def _append_restored_internal_row(self, fname: str, info: dict) -> None:
        if not hasattr(self, "internal_results_table"):
            return
        t = self.internal_results_table
        row = t.rowCount()
        t.insertRow(row)
        ok = self._is_successful_wcs_result(info)

        def _fmt(keys, fmt, default=""):
            for key in keys:
                val = self._finite_result_float(info.get(key))
                if np.isfinite(val):
                    return fmt.format(val)
            return default

        if ok:
            ra_s = _fmt(("center_ra_deg", "ra"), "{:.4f}")
            dec_s = _fmt(("center_dec_deg", "dec"), "{:.4f}")
            pix_s = _fmt(("scale_arcsec_per_px", "pix_fit", "pixscale"), "{:.4f}\"/px")
            model_s = str(info.get("model", "") or "").strip()
            refine_s = str(info.get("refine", "") or "").strip()
            if not refine_s and info.get("n_matches") is not None:
                try:
                    refine_s = f"{model_s} m={int(info.get('n_matches', 0) or 0)}".strip()
                except Exception:
                    refine_s = model_s
            resid_s = _fmt(("rms_arcsec", "resid_med", "resid_med_px"), "{:.3f}")
            status_s = str(info.get("status", "") or "ok")
        else:
            ra_s = dec_s = pix_s = refine_s = resid_s = ""
            status_s = str(info.get("status", "") or info.get("reason", "") or "fail")

        elapsed_s = _fmt(("elapsed_s", "elapsed"), "{:.2f}")
        hint_s = str(info.get("hint_source", "") or "")
        values = [fname, status_s, hint_s, ra_s, dec_s, pix_s, refine_s, resid_s, elapsed_s]
        for col, val in enumerate(values):
            t.setItem(row, col, QTableWidgetItem(val))
        set_table_row_background(t, row, self._wcs_result_row_background(info))

    def _append_restored_astnet_row(self, fname: str, info: dict) -> None:
        if not hasattr(self, "astrometrynet_results_table"):
            return
        t = self.astrometrynet_results_table
        row = t.rowCount()
        t.insertRow(row)
        pixscale = self._finite_result_float(info.get("pixscale", info.get("pix_fit")))
        resid_med = self._finite_result_float(info.get("resid_med", info.get("resid_med_px")))
        elapsed = self._finite_result_float(info.get("elapsed_s", info.get("elapsed")))
        values = [
            fname,
            str(info.get("status", "") or "solved"),
            f"{pixscale:.4f}" if np.isfinite(pixscale) and pixscale > 0 else "-",
            str(info.get("refine", "") or "-"),
            f"{resid_med:.2f}" if np.isfinite(resid_med) else "-",
            f"{elapsed:.1f}" if np.isfinite(elapsed) and elapsed > 0 else "-",
        ]
        for col, val in enumerate(values):
            t.setItem(row, col, QTableWidgetItem(val))
        set_table_row_background(t, row, self._wcs_result_row_background(info))

    def _append_restored_astap_row(self, fname: str, info: dict) -> None:
        if not hasattr(self, "results_table"):
            return
        row = self.results_table.rowCount()
        self.results_table.insertRow(row)
        self.results_table.setItem(row, 0, QTableWidgetItem(fname))
        self.results_table.setItem(row, 1, QTableWidgetItem(str(info.get("status", "") or "ok")))
        pix_fit = self._finite_result_float(info.get("pix_fit", info.get("pixscale")))
        self.results_table.setItem(
            row, 2,
            QTableWidgetItem(f"{pix_fit:.4f}" if np.isfinite(pix_fit) and pix_fit > 0 else "-"),
        )
        self.results_table.setItem(row, 3, QTableWidgetItem(str(info.get("refine", "") or "-")))
        resid_med = self._finite_result_float(info.get("resid_med"))
        if np.isfinite(resid_med):
            resid_str = f"{resid_med:.3f}\""
        else:
            resid_px = self._finite_result_float(info.get("resid_med_px"))
            resid_str = f"{resid_px:.3f}px" if np.isfinite(resid_px) else "-"
        self.results_table.setItem(row, 4, QTableWidgetItem(resid_str))
        elapsed = self._finite_result_float(info.get("elapsed", info.get("elapsed_s")), 0.0)
        self.results_table.setItem(row, 5, QTableWidgetItem(f"{elapsed:.1f}" if np.isfinite(elapsed) else "0.0"))
        set_table_row_background(self.results_table, row, self._wcs_result_row_background(info))

    def _populate_restored_wcs_tables(self, restored: dict[str, dict]) -> dict[str, int]:
        counts = {"internal": 0, "astnet": 0, "astap": 0}
        for table_name in ("internal_results_table", "results_table", "astrometrynet_results_table"):
            table = getattr(self, table_name, None)
            if table is not None:
                table.setRowCount(0)

        for fname, info in restored.items():
            tag = str(info.get("wcssrc") or info.get("solver") or "").strip().lower()
            if "internal" in tag or "apex_internal" in tag:
                self._append_restored_internal_row(fname, info)
                counts["internal"] += 1
            elif "astnet" in tag or "astrometry" in tag or "solve-field" in tag:
                self._append_restored_astnet_row(fname, info)
                counts["astnet"] += 1
            else:
                self._append_restored_astap_row(fname, info)
                counts["astap"] += 1
        return counts

    def _write_internal_summary_csv(self):
        """Persist ``self.results`` to ``step5/wcs_solve_summary.csv``.

        Schema-aligned with the ASTAP and astrometry.net writers so the same
        downstream loaders (project re-open, Step 6 QC, diagnostics) work
        unchanged. The Internal solver's native result dict uses slightly
        different field names than ASTAP's, so we map the common columns:

        ============================ ==========================================
        Internal key                  Aliased to (for cross-solver consumers)
        ============================ ==========================================
        ``rms_arcsec``                ``resid_med``
        ``scale_arcsec_per_px``       ``pix_fit``
        ``elapsed_s``                 ``elapsed``
        ``n_matches``                 ``refine`` (formatted as ``m={N}``)
        ============================ ==========================================

        We never drop the original field — the alias is added alongside so
        Internal-specific consumers and the cross-solver loader both see what
        they expect. ``file`` and ``fname`` are both written for the same
        reason (``_restore_success_results_from_summary`` accepts either).
        """
        source_results = getattr(self, "_internal_run_results", None) or self.results
        if not source_results:
            return
        rows: list[dict] = []
        for fname, info in source_results.items():
            row = dict(info) if isinstance(info, dict) else {}
            row.setdefault("file", str(fname))
            row.setdefault("fname", str(fname))
            ok = bool(row.get("ok", False))
            row["ok"] = ok
            row.setdefault("status", "ok" if ok else "fail")
            row.setdefault("solver", "apex_internal")
            row.setdefault("wcssrc", "APEX_INTERNAL")
            row.setdefault("hint_source", "")
            if "rms_arcsec" in row and "resid_med" not in row:
                row["resid_med"] = row["rms_arcsec"]
            if "scale_arcsec_per_px" in row and "pix_fit" not in row:
                row["pix_fit"] = row["scale_arcsec_per_px"]
            if "elapsed_s" in row and "elapsed" not in row:
                row["elapsed"] = row["elapsed_s"]
            if "n_matches" in row and "refine" not in row:
                try:
                    row["refine"] = f"m={int(row.get('n_matches', 0) or 0)}"
                except Exception:
                    row["refine"] = ""
            rows.append(row)
        step5_out = step5_wcs_dir(self.params.P.result_dir)
        step5_out.mkdir(parents=True, exist_ok=True)
        df = pd.DataFrame(rows)
        out_path = step5_out / "wcs_solve_summary.csv"
        df.to_csv(out_path, index=False)
        self.log(f"[Internal WCS] Summary CSV: {len(df)} rows → {out_path.name}")

    def _on_internal_error(self, msg: str):
        self.run_bar_internal.set_running(False)
        self.log(f"[Internal WCS] ERROR: {msg}")
        QMessageBox.warning(self, "Internal WCS Solver Error", msg)

    def setup_astap_tab(self):
        """Setup ASTAP tab UI"""
        layout = QVBoxLayout(self.astap_tab)

        info = QLabel(
            "Solve WCS for all frames using ASTAP (local). "
            "ASTAP and its D80/D50 star database are installed separately; "
            "see Help > WCS Solver Installation Help."
        )
        info.setWordWrap(True)
        info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 10px; border-radius: 5px; }")
        layout.addWidget(info)

        control_layout = QHBoxLayout()
        btn_params = create_parameter_button("ASTAP Parameters")
        btn_params.clicked.connect(self.open_parameters_dialog)
        control_layout.addWidget(btn_params)

        self.run_bar_astap = RunControlBar(
            "Run ASTAP", "Log",
            run_cb=self.run_wcs,
            stop_cb=self.stop_wcs,
            log_cb=self.show_log_window,
        )
        control_layout.addWidget(self.run_bar_astap)
        self.btn_run = self.run_bar_astap.btn_run
        self.btn_stop = self.run_bar_astap.btn_stop

        layout.addLayout(control_layout)

        progress_layout = QHBoxLayout()
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        progress_layout.addWidget(self.progress_bar)

        self.progress_label = QLabel("Ready")
        self.progress_label.setMinimumWidth(350)
        progress_layout.addWidget(self.progress_label)
        layout.addLayout(progress_layout)

        results_group = QGroupBox("WCS Results")
        results_layout = QVBoxLayout(results_group)

        self.results_table = QTableWidget()
        self.results_table.setColumnCount(6)
        self.results_table.setHorizontalHeaderLabels([
            "File", "Status", "PixScale Fit", "Refine", "Resid Med", "Elapsed (s)"
        ])
        self.results_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.results_table.horizontalHeader().setStretchLastSection(True)
        self.results_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.results_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        results_layout.addWidget(self.results_table)

        layout.addWidget(results_group)

    def setup_astrometrynet_tab(self):
        """Setup Astrometry.net tab UI"""
        layout = QVBoxLayout(self.astrometrynet_tab)

        info_text = (
            "Solve WCS for all frames using local astrometry.net (solve-field). "
            "This optional fallback needs a local solve-field installation and "
            "matching astrometry.net index files; see Help > WCS Solver Installation Help."
        )
        info_style = "QLabel { background-color: #E8F5E9; padding: 10px; border-radius: 5px; }"
        info = QLabel(info_text)
        info.setWordWrap(True)
        info.setStyleSheet(info_style)
        layout.addWidget(info)

        control_layout = QHBoxLayout()

        btn_astnet_params = create_parameter_button("Astrometry.net Parameters")
        btn_astnet_params.clicked.connect(self.open_astrometrynet_parameters_dialog)
        control_layout.addWidget(btn_astnet_params)

        self.run_bar_astnet = RunControlBar(
            "Solve All Frames", "Log",
            run_cb=self.run_astrometrynet_solve,
            stop_cb=self.stop_astrometrynet_solve,
            log_cb=self.show_log_window,
        )
        control_layout.addWidget(self.run_bar_astnet)
        self.btn_solve_astrometrynet = self.run_bar_astnet.btn_run
        self.btn_stop_astrometrynet = self.run_bar_astnet.btn_stop

        layout.addLayout(control_layout)

        progress_layout = QHBoxLayout()
        self.astrometrynet_progress = QProgressBar()
        self.astrometrynet_progress.setMinimum(0)
        self.astrometrynet_progress.setMaximum(100)
        self.astrometrynet_progress.setValue(0)
        progress_layout.addWidget(self.astrometrynet_progress)

        self.astrometrynet_status = QLabel("Ready")
        self.astrometrynet_status.setMinimumWidth(350)
        progress_layout.addWidget(self.astrometrynet_status)
        layout.addLayout(progress_layout)

        results_group = QGroupBox("Astrometry.net Results")
        results_layout = QVBoxLayout(results_group)

        self.astrometrynet_results_table = QTableWidget()
        self.astrometrynet_results_table.setColumnCount(6)
        self.astrometrynet_results_table.setHorizontalHeaderLabels([
            "File", "Status", "PixScale", "Refine", "Resid(\")", "Elapsed (s)"
        ])
        self.astrometrynet_results_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.astrometrynet_results_table.horizontalHeader().setStretchLastSection(True)
        self.astrometrynet_results_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.astrometrynet_results_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        results_layout.addWidget(self.astrometrynet_results_table)

        layout.addWidget(results_group)

    def run_astrometrynet_solve(self):
        if not self.file_list:
            QMessageBox.warning(self, "Warning", "No frames found to solve")
            return

        self.log_text.clear()
        file_list = self._qc_filter_and_log()
        if file_list is None:
            return

        pf = self._log_preflight(astnet=True, gaia=True)
        astnet_ok, astnet_detail = pf["astnet"]
        gaia_ok, gaia_detail = pf["gaia"]
        if not astnet_ok:
            self.astrometrynet_status.setText("solve-field not found or not configured")
            QMessageBox.warning(
                self,
                "Astrometry.net Not Found",
                f"{astnet_detail}\n\n"
                "Install local astrometry.net with matching index files, or update Step 5 > "
                "Astrometry.net Parameters. See Help > WCS Solver Installation Help.",
            )
            return

        if not gaia_ok:
            QMessageBox.warning(
                self,
                "Gaia Runtime Not Available",
                f"{gaia_detail}\n\n"
                "Local astrometry.net can still solve WCS, but Gaia attach, "
                "refine, residual medians, and Step 6 Gaia stats will be absent "
                "unless a compatible gaia_fov.ecsv cache already exists.",
            )

        target_coord = None
        ra = getattr(self.params.P, "target_ra_deg", None)
        dec = getattr(self.params.P, "target_dec_deg", None)
        if ra is not None and dec is not None:
            try:
                target_coord = SkyCoord(float(ra) * u.deg, float(dec) * u.deg)
            except Exception:
                target_coord = None

        self.astrometrynet_results_table.setRowCount(0)
        self.results = {}

        # Start solving in thread
        self.astrometrynet_worker = AstrometryNetWorker(
            file_list,
            self.params,
            self.params.P.data_dir,
            self.params.P.result_dir,
            self.params.P.cache_dir,
            self.use_cropped,
            target_coord=target_coord,
        )
        self.astrometrynet_worker.progress.connect(self.on_astrometrynet_progress)
        self.astrometrynet_worker.file_done.connect(self.on_astrometrynet_file_done)
        self.astrometrynet_worker.refine_done.connect(self.on_astrometrynet_refine_done)
        self.astrometrynet_worker.log_message.connect(self.log)
        self.astrometrynet_worker.finished.connect(self.on_astrometrynet_finished)
        self.astrometrynet_worker.error.connect(self.on_astrometrynet_error)
        self.setup_log_window()
        if hasattr(self, "_worker_panel") and self._worker_panel is not None:
            self._worker_panel.clear()
            self.astrometrynet_worker.worker_status.connect(self._worker_panel.update_worker)

        self.run_bar_astnet.set_running(True)
        self.astrometrynet_progress.setValue(0)
        self._astnet_start_time = time.monotonic()
        self.astrometrynet_status.setText(
            progress_status_text(0, len(file_list), self._astnet_start_time, message="Starting local astrometry.net...")
        )
        self.log("=" * 50)
        self.log("Starting local astrometry.net (solve-field) plate solving...")
        self.log(f"Frames: {len(file_list)}")
        self.astrometrynet_worker.start()
        self.show_log_window()

    def stop_astrometrynet_solve(self):
        if self.astrometrynet_worker and self.astrometrynet_worker.isRunning():
            self.run_bar_astnet.set_stopping()
            self.astrometrynet_status.setText("Stopping...")
            self.log("Astrometry.net stop requested...")
            self.astrometrynet_worker.stop()

    def on_astrometrynet_progress(self, current, total, status):
        pct = int(100 * current / max(1, total))
        self.astrometrynet_progress.setValue(pct)
        self.astrometrynet_status.setText(
            progress_status_text(current, total, getattr(self, "_astnet_start_time", None), message=status)
        )

    def on_astrometrynet_file_done(self, filename, result):
        row = self.astrometrynet_results_table.rowCount()
        self.astrometrynet_results_table.insertRow(row)
        self.astrometrynet_results_table.setItem(row, 0, QTableWidgetItem(filename))
        self.astrometrynet_results_table.setItem(row, 1, QTableWidgetItem(result.get("status", "")))
        pixscale = float(result.get("pixscale", 0.0))
        refine = result.get("refine", "-")
        resid_med = result.get("resid_med", np.nan)
        elapsed = float(result.get("elapsed_s", 0.0))
        self.astrometrynet_results_table.setItem(row, 2, QTableWidgetItem(f"{pixscale:.4f}" if np.isfinite(pixscale) and pixscale > 0 else "-"))
        self.astrometrynet_results_table.setItem(row, 3, QTableWidgetItem(str(refine) if refine else "-"))
        self.astrometrynet_results_table.setItem(row, 4, QTableWidgetItem(f"{resid_med:.2f}" if np.isfinite(resid_med) else "-"))
        self.astrometrynet_results_table.setItem(row, 5, QTableWidgetItem(f"{elapsed:.1f}" if np.isfinite(elapsed) and elapsed > 0 else "-"))
        set_table_row_background(
            self.astrometrynet_results_table,
            row,
            self._wcs_result_row_background(result),
        )

        if result.get("ok"):
            self.results[filename] = result
            self.log(f"Astrometry.net solved: {filename} (RA={result.get('ra', 0):.4f}, Dec={result.get('dec', 0):.4f})")
        try:
            self._refresh_solver_breakdown_label()
        except Exception:
            pass

    def on_astrometrynet_error(self, filename, error):
        self.log(f"Astrometry.net ERROR {filename}: {error}")

    def on_astrometrynet_refine_done(self, filename, result):
        """Refine 결과로 테이블 업데이트"""
        # 파일명으로 해당 행 찾기
        for row in range(self.astrometrynet_results_table.rowCount()):
            item = self.astrometrynet_results_table.item(row, 0)
            if item and item.text() == filename:
                refine = result.get("refine", "-")
                resid_med = result.get("resid_med", np.nan)
                self.astrometrynet_results_table.setItem(row, 3, QTableWidgetItem(str(refine) if refine else "-"))
                self.astrometrynet_results_table.setItem(row, 4, QTableWidgetItem(f"{resid_med:.2f}" if np.isfinite(resid_med) else "-"))
                # results에도 업데이트
                if filename in self.results:
                    self.results[filename].update(result)
                    row_result = self.results[filename]
                else:
                    row_result = result
                set_table_row_background(
                    self.astrometrynet_results_table,
                    row,
                    self._wcs_result_row_background(row_result),
                )
                break

    def on_astrometrynet_finished(self, summary):
        self.run_bar_astnet.set_running(False)
        stopped = bool(summary.get("stopped")) if isinstance(summary, dict) else False
        n_ok = summary.get("ok", 0)
        n_qc = summary.get("wcs_qc_pass", 0)
        n_qc_not_eval = int(summary.get("wcs_qc_not_evaluated", 0) or 0)
        if not stopped:
            self.astrometrynet_progress.setValue(100)
            self.astrometrynet_status.setText(f"Done: {n_ok}/{summary.get('total', 0)} solved")
        else:
            self.astrometrynet_status.setText(f"Stopped: {n_ok}/{summary.get('total', 0)} solved")
        if n_ok > 0:
            qc_text = (
                f"WCS-QC not evaluated: {n_qc_not_eval} (Gaia unavailable)"
                if n_qc_not_eval
                else f"WCS-QC pass: {n_qc}"
            )
            self.log(f"Astrometry.net: {n_ok} frames solved successfully | {qc_text}")
        if stopped:
            self.log("Astrometry.net solve stopped by user")
        self.save_state()
        self.update_navigation_buttons()

    def setup_log_window(self):
        if self.log_window is not None:
            return

        worker_group = QGroupBox("Workers")
        worker_group.setMinimumWidth(430)
        wg_layout = QVBoxLayout(worker_group)
        wg_layout.setContentsMargins(5, 5, 5, 5)
        self._worker_panel = WorkerStatusPanel(worker_group)
        wg_layout.addWidget(self._worker_panel)

        self.log_window = WorkflowLogWindow(
            self, "WCS Log & Workers", width=900, height=500,
            side_widget=worker_group,
        )
        self.log_text = self.log_window.log_text

    def show_log_window(self):
        if self.log_window is None:
            self.setup_log_window()
        show_raised(self.log_window)

    def log(self, message: str):
        append_timestamped_log(self.log_text, message)

    def populate_file_list(self):
        crop_active = crop_is_active(self.params.P.result_dir)
        cropped_dir = step2_cropped_dir(self.params.P.result_dir)
        excluded = getattr(self.file_manager, "excluded_files", set()) if self.file_manager else set()

        if crop_active and cropped_dir.exists() and list(cropped_dir.glob("*.fit*")):
            files = sorted([f.name for f in cropped_dir.glob("*.fit*")
                            if f.name not in excluded])
            self.use_cropped = True
        else:
            files = self.file_manager.get_file_list() if self.file_manager else []
            self.use_cropped = False

        self.file_list = list(files)

    def open_parameters_dialog(self):
        dialog = FittedDialog(self)
        configure_parameter_dialog(dialog, "WCS Parameters", 560, 720)

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

        _info = QLabel(
            "Adjust WCS plate-solving parameters. ASTAP and the selected D80/D50 "
            "star database must be installed outside APEX. See Help > WCS Solver "
            "Installation Help."
        )
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

        astap_form = _add_group("ASTAP Solver", expanded=True)
        refine_form = _add_group("WCS Refinement", expanded=True)
        gaia_form = _add_group("Gaia / Hybrid ID")

        self.param_astap_exe = QLineEdit(str(getattr(self.params.P, "astap_exe", "astap_cli.exe")))
        astap_exe_row = QWidget()
        astap_exe_layout = QHBoxLayout(astap_exe_row)
        astap_exe_layout.setContentsMargins(0, 0, 0, 0)
        astap_exe_layout.setSpacing(6)
        astap_exe_layout.addWidget(self.param_astap_exe, 1)
        btn_browse_astap = QPushButton("Browse...")
        btn_browse_astap.clicked.connect(self.browse_astap_cli)
        astap_exe_layout.addWidget(btn_browse_astap)
        astap_form.addRow("ASTAP CLI Path:", astap_exe_row)

        self.param_timeout = QDoubleSpinBox()
        self.param_timeout.setRange(10, 1000)
        self.param_timeout.setValue(float(getattr(self.params.P, "astap_timeout_s", 120.0)))
        astap_form.addRow("Timeout (s):", self.param_timeout)

        self.param_radius = QDoubleSpinBox()
        self.param_radius.setRange(0.5, 30.0)
        self.param_radius.setValue(float(getattr(self.params.P, "astap_search_radius_deg", 8.0)))
        astap_form.addRow("Search Radius (deg):", self.param_radius)

        self.param_astap_db = QComboBox()
        self.param_astap_db.addItems(["D80", "D50"])
        current_db = str(getattr(self.params.P, "astap_database", "D80"))
        idx = self.param_astap_db.findText(current_db)
        if idx >= 0:
            self.param_astap_db.setCurrentIndex(idx)
        astap_form.addRow("ASTAP Star DB:", self.param_astap_db)

        self.param_annotate_variables = QCheckBox("Enable")
        self.param_annotate_variables.setChecked(bool(getattr(self.params.P, "astap_annotate_variables", False)))
        self.param_annotate_variables.setToolTip("ASTAP 변광성 데이터베이스로 변광성 주석 표시 (별도 설치 필요)")
        astap_form.addRow("Annotate Variable Stars:", self.param_annotate_variables)

        self.param_fov_fudge = QDoubleSpinBox()
        self.param_fov_fudge.setRange(0.5, 2.0)
        self.param_fov_fudge.setSingleStep(0.05)
        self.param_fov_fudge.setValue(float(getattr(self.params.P, "astap_fov_fudge", 1.0)))
        astap_form.addRow("FOV Fudge:", self.param_fov_fudge)

        self.param_downsample = QSpinBox()
        self.param_downsample.setRange(1, 8)
        self.param_downsample.setValue(int(getattr(self.params.P, "astap_downsample_z", 2)))
        astap_form.addRow("Downsample Z:", self.param_downsample)

        self.param_max_stars = QSpinBox()
        self.param_max_stars.setRange(50, 5000)
        self.param_max_stars.setValue(int(getattr(self.params.P, "astap_max_stars_s", 500)))
        astap_form.addRow("Max Stars (S):", self.param_max_stars)

        self.param_max_workers = QSpinBox()
        self.param_max_workers.setRange(1, 16)
        self.param_max_workers.setValue(int(getattr(self.params.P, "wcs_max_workers", 1)))
        astap_form.addRow("Max Workers:", self.param_max_workers)

        self.param_require_qc = QCheckBox("Enable")
        self.param_require_qc.setChecked(bool(getattr(self.params.P, "wcs_require_qc_pass", True)))
        refine_form.addRow("QC Pass Only:", self.param_require_qc)

        self.param_refine_enable = QCheckBox("Enable")
        self.param_refine_enable.setChecked(bool(getattr(self.params.P, "wcs_refine_enable", True)))
        refine_form.addRow("Refine CRPIX:", self.param_refine_enable)

        self.param_refine_max_match = QSpinBox()
        self.param_refine_max_match.setRange(50, 5000)
        self.param_refine_max_match.setValue(int(getattr(self.params.P, "wcs_refine_max_match", 600)))
        refine_form.addRow("Refine Max Match:", self.param_refine_max_match)

        self.param_refine_match_r = QDoubleSpinBox()
        self.param_refine_match_r.setRange(0.5, 5.0)
        self.param_refine_match_r.setSingleStep(0.1)
        self.param_refine_match_r.setValue(float(getattr(self.params.P, "wcs_refine_match_r_fwhm", 1.6)))
        refine_form.addRow("Refine Match R (×FWHM):", self.param_refine_match_r)

        self.param_refine_min_match = QSpinBox()
        self.param_refine_min_match.setRange(5, 500)
        self.param_refine_min_match.setValue(int(getattr(self.params.P, "wcs_refine_min_match", 50)))
        refine_form.addRow("Refine Min Match:", self.param_refine_min_match)

        self.param_gaia_fudge = QDoubleSpinBox()
        self.param_gaia_fudge.setRange(0.5, 3.0)
        self.param_gaia_fudge.setSingleStep(0.05)
        self.param_gaia_fudge.setValue(float(getattr(self.params.P, "gaia_radius_fudge", 1.35)))
        gaia_form.addRow("Gaia Radius Fudge:", self.param_gaia_fudge)

        self.param_gaia_mag_max = QDoubleSpinBox()
        self.param_gaia_mag_max.setRange(10.0, 25.0)
        self.param_gaia_mag_max.setSingleStep(0.5)
        self.param_gaia_mag_max.setValue(float(getattr(self.params.P, "gaia_mag_max", 18.0)))
        gaia_form.addRow("Gaia Mag Max:", self.param_gaia_mag_max)

        self.param_gaia_wcs_mag_max = QDoubleSpinBox()
        self.param_gaia_wcs_mag_max.setRange(10.0, 21.0)
        self.param_gaia_wcs_mag_max.setSingleStep(0.5)
        self.param_gaia_wcs_mag_max.setValue(float(getattr(self.params.P, "gaia_wcs_mag_max", 18.0)))
        self.param_gaia_wcs_mag_max.setToolTip(
            "Server-side Gaia G magnitude cap used by Step 5 WCS refinement/QC. "
            "Keep this near 18-20 to avoid TAP timeouts."
        )
        gaia_form.addRow("Gaia WCS Mag Max:", self.param_gaia_wcs_mag_max)

        self.param_ref_gaia_match_tol = QDoubleSpinBox()
        self.param_ref_gaia_match_tol.setRange(0.1, 30.0)
        self.param_ref_gaia_match_tol.setDecimals(2)
        self.param_ref_gaia_match_tol.setSingleStep(0.1)
        self.param_ref_gaia_match_tol.setValue(float(getattr(self.params.P, "ref_wcs_match_radius_arcsec", 2.0)))
        gaia_form.addRow("Gaia Match Tol (Ref, arcsec):", self.param_ref_gaia_match_tol)

        self.param_gaia_g_limit = QDoubleSpinBox()
        self.param_gaia_g_limit.setRange(10.0, 25.0)
        self.param_gaia_g_limit.setDecimals(2)
        self.param_gaia_g_limit.setSingleStep(0.5)
        self.param_gaia_g_limit.setValue(
            float(
                getattr(
                    self.params.P,
                    "idmatch_gaia_g_limit",
                    getattr(self.params.P, "gaia_mag_max", 18.0),
                )
            )
        )
        gaia_form.addRow("Gaia G limit (Hybrid ID):", self.param_gaia_g_limit)

        self.param_gaia_retry = QSpinBox()
        self.param_gaia_retry.setRange(0, 10)
        self.param_gaia_retry.setValue(int(getattr(self.params.P, "gaia_retry", 2)))
        gaia_form.addRow("Gaia Retry:", self.param_gaia_retry)

        self.param_gaia_timeout = QDoubleSpinBox()
        self.param_gaia_timeout.setRange(5.0, 300.0)
        self.param_gaia_timeout.setSingleStep(5.0)
        self.param_gaia_timeout.setValue(float(getattr(self.params.P, "gaia_timeout_s", 30.0)))
        self.param_gaia_timeout.setSuffix(" s")
        gaia_form.addRow("Gaia Timeout:", self.param_gaia_timeout)

        self.param_gaia_backoff = QDoubleSpinBox()
        self.param_gaia_backoff.setRange(0.0, 30.0)
        self.param_gaia_backoff.setSingleStep(1.0)
        self.param_gaia_backoff.setValue(float(getattr(self.params.P, "gaia_backoff_s", 6.0)))
        gaia_form.addRow("Gaia Backoff (s):", self.param_gaia_backoff)

        self.param_gaia_allow_no_cache = QCheckBox("Allow query when cache missing")
        self.param_gaia_allow_no_cache.setChecked(bool(getattr(self.params.P, "gaia_allow_no_cache", True)))
        self.param_gaia_allow_no_cache.setToolTip(
            "This is not a reuse-cache toggle. It controls whether Step 5 may query Gaia online "
            "when no local Gaia cache exists."
        )
        gaia_form.addRow("Gaia cache miss:", self.param_gaia_allow_no_cache)

        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        add_parameter_reset_button(
            buttons,
            [
                (self.param_astap_exe, "C:/Program Files/astap/astap_cli.exe"),
                (self.param_timeout, 120.0),
                (self.param_radius, 8.0),
                (self.param_astap_db, "D50"),
                (self.param_annotate_variables, False),
                (self.param_fov_fudge, 1.0),
                (self.param_downsample, 2),
                (self.param_max_stars, 500),
                (self.param_max_workers, 4),
                (self.param_require_qc, True),
                (self.param_refine_enable, True),
                (self.param_refine_max_match, 600),
                (self.param_refine_match_r, 2.0),
                (self.param_refine_min_match, 25),
                (self.param_gaia_fudge, 1.35),
                (self.param_gaia_mag_max, 25.0),
                (self.param_gaia_wcs_mag_max, 18.0),
                (self.param_ref_gaia_match_tol, 2.0),
                (self.param_gaia_g_limit, 25.0),
                (self.param_gaia_retry, 2),
                (self.param_gaia_timeout, 30.0),
                (self.param_gaia_backoff, 6.0),
                (self.param_gaia_allow_no_cache, True),
            ],
        )
        buttons.accepted.connect(lambda: self.save_parameters(dialog))
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        dialog.exec_()

    def open_astrometrynet_parameters_dialog(self):
        dialog = FittedDialog(self)
        configure_parameter_dialog(dialog, "Astrometry.net Parameters", 540, 620)

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

        _info = QLabel(
            "Configure Astrometry.net local solver as a fallback for ASTAP failures. "
            "Local solve-field and matching index files must be installed outside APEX. "
            "See Help > WCS Solver Installation Help."
        )
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

        form = _add_group("Local Solve", expanded=True)
        retry_form = _add_group("Blind Retry")

        self.param_astnet_enable = QCheckBox("Enable")
        self.param_astnet_enable.setChecked(bool(getattr(self.params.P, "astnet_local_enable", False)))
        form.addRow("Enable Local Solve:", self.param_astnet_enable)

        self.param_astnet_use_wsl = QCheckBox("Use WSL")
        self.param_astnet_use_wsl.setChecked(bool(getattr(self.params.P, "astnet_local_use_wsl", True)))
        form.addRow("Use WSL:", self.param_astnet_use_wsl)

        self.param_astnet_command = QLineEdit(str(getattr(self.params.P, "astnet_local_command", "solve-field")))
        form.addRow("solve-field Command:", self.param_astnet_command)

        self.param_astnet_timeout = QDoubleSpinBox()
        self.param_astnet_timeout.setRange(30, 3600)
        self.param_astnet_timeout.setValue(float(getattr(self.params.P, "astnet_local_timeout_s", 300.0)))
        form.addRow("Timeout (s):", self.param_astnet_timeout)

        self.param_astnet_downsample = QSpinBox()
        self.param_astnet_downsample.setRange(1, 8)
        self.param_astnet_downsample.setValue(int(getattr(self.params.P, "astnet_local_downsample", 2)))
        form.addRow("Downsample:", self.param_astnet_downsample)

        self.param_astnet_scale_low = QDoubleSpinBox()
        self.param_astnet_scale_low.setRange(0.0, 10.0)
        self.param_astnet_scale_low.setDecimals(5)
        self.param_astnet_scale_low.setValue(float(getattr(self.params.P, "astnet_local_scale_low", 0.0)))
        form.addRow("Scale Low (arcsec/pix):", self.param_astnet_scale_low)

        self.param_astnet_scale_high = QDoubleSpinBox()
        self.param_astnet_scale_high.setRange(0.0, 10.0)
        self.param_astnet_scale_high.setDecimals(5)
        self.param_astnet_scale_high.setValue(float(getattr(self.params.P, "astnet_local_scale_high", 0.0)))
        form.addRow("Scale High (arcsec/pix):", self.param_astnet_scale_high)

        self.param_astnet_radius = QDoubleSpinBox()
        self.param_astnet_radius.setRange(0.1, 30.0)
        self.param_astnet_radius.setValue(float(getattr(self.params.P, "astnet_local_radius_deg", 8.0)))
        form.addRow("Radius (deg):", self.param_astnet_radius)

        self.param_astnet_keep_outputs = QCheckBox("Keep .wcs/debug files")
        self.param_astnet_keep_outputs.setChecked(bool(getattr(self.params.P, "astnet_local_keep_outputs", True)))
        self.param_astnet_keep_outputs.setToolTip(
            "Keep .wcs and debug sidecars; redundant .new FITS copies are removed automatically."
        )
        form.addRow("Keep Debug Outputs:", self.param_astnet_keep_outputs)

        self.param_astnet_use_cache = create_cache_checkbox(
            "Use Cache",
            bool(getattr(self.params.P, "astnet_local_use_cache", True)),
            "Reuse compatible local astrometry.net sidecar outputs instead of running solve-field again.",
        )
        form.addRow("Use Cached Outputs:", self.param_astnet_use_cache)

        self.param_astnet_max_objs = QSpinBox()
        self.param_astnet_max_objs.setRange(100, 20000)
        self.param_astnet_max_objs.setValue(int(getattr(self.params.P, "astnet_local_max_objs", 2000)))
        form.addRow("Max Objects:", self.param_astnet_max_objs)

        self.param_astnet_cpulimit = QDoubleSpinBox()
        self.param_astnet_cpulimit.setRange(1, 300)
        self.param_astnet_cpulimit.setValue(float(getattr(self.params.P, "astnet_local_cpulimit_s", 30.0)))
        form.addRow("CPU Limit (s):", self.param_astnet_cpulimit)

        self.param_astnet_blind_retry = QCheckBox("Retry blind when hint-based solve fails")
        self.param_astnet_blind_retry.setChecked(bool(getattr(self.params.P, "astnet_blind_retry_on_fail", True)))
        retry_form.addRow("Blind Retry:", self.param_astnet_blind_retry)

        self.param_astnet_blind_cpulimit = QDoubleSpinBox()
        self.param_astnet_blind_cpulimit.setRange(10, 600)
        self.param_astnet_blind_cpulimit.setValue(float(getattr(self.params.P, "astnet_blind_cpulimit_s", 120.0)))
        self.param_astnet_blind_cpulimit.setSuffix(" s")
        retry_form.addRow("Blind CPU Limit:", self.param_astnet_blind_cpulimit)

        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        add_parameter_reset_button(
            buttons,
            [
                (self.param_astnet_enable, True),
                (self.param_astnet_use_wsl, True),
                (self.param_astnet_command, "solve-field"),
                (self.param_astnet_timeout, 300.0),
                (self.param_astnet_downsample, 2),
                (self.param_astnet_scale_low, 0.3),
                (self.param_astnet_scale_high, 0.5),
                (self.param_astnet_radius, 30.0),
                (self.param_astnet_keep_outputs, True),
                (self.param_astnet_use_cache, True),
                (self.param_astnet_max_objs, 2000),
                (self.param_astnet_cpulimit, 30.0),
                (self.param_astnet_blind_retry, True),
                (self.param_astnet_blind_cpulimit, 120.0),
            ],
        )
        buttons.accepted.connect(lambda: self.save_astrometrynet_parameters(dialog))
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        dialog.exec_()

    def browse_astap_cli(self):
        current = _strip_outer_quotes(self.param_astap_exe.text())
        start_dir = Path(current).parent if current and Path(current).parent.exists() else Path.home()
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select ASTAP CLI executable",
            str(start_dir),
            "ASTAP CLI (astap_cli.exe);;Executables (*.exe);;All Files (*)",
        )
        if path:
            self.param_astap_exe.setText(path)

    def save_parameters(self, dialog):
        self.params.P.astap_exe = self.param_astap_exe.text().strip()
        self.params.P.astap_timeout_s = self.param_timeout.value()
        self.params.P.astap_search_radius_deg = self.param_radius.value()
        self.params.P.astap_database = self.param_astap_db.currentText()
        self.params.P.astap_annotate_variables = self.param_annotate_variables.isChecked()
        self.params.P.astap_fov_fudge = self.param_fov_fudge.value()
        self.params.P.astap_downsample_z = self.param_downsample.value()
        self.params.P.astap_max_stars_s = self.param_max_stars.value()
        self.params.P.wcs_max_workers = self.param_max_workers.value()
        self.params.P.wcs_require_qc_pass = self.param_require_qc.isChecked()
        self.params.P.wcs_refine_enable = self.param_refine_enable.isChecked()
        self.params.P.wcs_refine_max_match = self.param_refine_max_match.value()
        self.params.P.wcs_refine_match_r_fwhm = self.param_refine_match_r.value()
        self.params.P.wcs_refine_min_match = self.param_refine_min_match.value()
        self.params.P.gaia_radius_fudge = self.param_gaia_fudge.value()
        self.params.P.gaia_mag_max = self.param_gaia_mag_max.value()
        self.params.P.gaia_wcs_mag_max = self.param_gaia_wcs_mag_max.value()
        self.params.P.ref_wcs_match_radius_arcsec = self.param_ref_gaia_match_tol.value()
        self.params.P.idmatch_gaia_g_limit = self.param_gaia_g_limit.value()
        self.params.P.gaia_retry = self.param_gaia_retry.value()
        self.params.P.gaia_timeout_s = self.param_gaia_timeout.value()
        self.params.P.gaia_backoff_s = self.param_gaia_backoff.value()
        self.params.P.gaia_allow_no_cache = self.param_gaia_allow_no_cache.isChecked()
        self.persist_params()
        self.save_state()
        QMessageBox.information(dialog, "Success", "Parameters saved!")
        dialog.accept()

    def save_astrometrynet_parameters(self, dialog):
        self.params.P.astnet_local_enable = self.param_astnet_enable.isChecked()
        self.params.P.astnet_local_use_wsl = self.param_astnet_use_wsl.isChecked()
        self.params.P.astnet_local_command = self.param_astnet_command.text().strip()
        self.params.P.astnet_local_timeout_s = self.param_astnet_timeout.value()
        self.params.P.astnet_local_downsample = self.param_astnet_downsample.value()
        self.params.P.astnet_local_scale_low = self.param_astnet_scale_low.value()
        self.params.P.astnet_local_scale_high = self.param_astnet_scale_high.value()
        self.params.P.astnet_local_radius_deg = self.param_astnet_radius.value()
        self.params.P.astnet_local_keep_outputs = self.param_astnet_keep_outputs.isChecked()
        self.params.P.astnet_local_use_cache = self.param_astnet_use_cache.isChecked()
        self.params.P.astnet_local_max_objs = self.param_astnet_max_objs.value()
        self.params.P.astnet_local_cpulimit_s = self.param_astnet_cpulimit.value()
        self.params.P.astnet_blind_retry_on_fail = self.param_astnet_blind_retry.isChecked()
        self.params.P.astnet_blind_cpulimit_s = self.param_astnet_blind_cpulimit.value()
        self.persist_params()
        self.save_state()
        QMessageBox.information(dialog, "Success", "Astrometry.net parameters saved!")
        dialog.accept()

    def run_wcs(self):
        if not self.file_list:
            QMessageBox.warning(self, "Warning", "No files to process")
            return
        if self.worker and self.worker.isRunning():
            return

        self.results = {}
        self.results_table.setRowCount(0)
        self.log_text.clear()
        file_list = self._qc_filter_and_log()
        if file_list is None:
            return
        self.stop_requested = False

        astap_ok, astap_detail = self._log_preflight(astap=True)["astap"]
        astnet_enabled = bool(getattr(self.params.P, "astnet_local_enable", False))
        # astrometry.net is a fallback for ASTAP: only probe it (slow WSL spin-up)
        # when ASTAP is unavailable, otherwise log a deferred note instead.
        if astnet_enabled and not astap_ok:
            astnet_ok, astnet_detail = self._log_preflight(astnet=True)["astnet"]
        else:
            astnet_ok = False
            astnet_detail = (
                "local astrometry.net fallback deferred until ASTAP fails"
                if astnet_enabled
                else "local astrometry.net fallback disabled"
            )
            if not astap_ok:
                self.log(f"[WCS][Preflight] astrometry.net: {astnet_detail}")

        gaia_ok, gaia_detail = self._log_preflight(gaia=True)["gaia"]
        if not gaia_ok:
            QMessageBox.warning(
                self,
                "Gaia Runtime Not Available",
                f"{gaia_detail}\n\n"
                "ASTAP can still solve WCS, but Gaia attach, refine, residual "
                "medians, and Step 6 Gaia stats will be absent unless a compatible "
                "gaia_fov.ecsv cache already exists.",
            )

        if not astap_ok and not astnet_ok:
            self.progress_label.setText("WCS solver not found or not configured")
            QMessageBox.warning(
                self,
                "WCS Solver Not Found",
                f"{astap_detail}\n\n{astnet_detail}\n\n"
                "Install ASTAP with a D80/D50 star database, or enable/install local "
                "astrometry.net with matching index files. See Help > WCS Solver "
                "Installation Help.",
            )
            return

        if not astap_ok and astnet_ok:
            self.progress_label.setText("ASTAP not found; using local astrometry.net fallback")
            QMessageBox.information(
                self,
                "ASTAP Not Found",
                f"{astap_detail}\n\n"
                "Local astrometry.net is enabled and available, so APEX will continue "
                "with solve-field fallback. See Help > WCS Solver Installation Help "
                "for ASTAP/D80/D50 setup.",
            )

        # Get target from params.P (loaded from TOML)
        target_coord = None
        ra = getattr(self.params.P, "target_ra_deg", None)
        dec = getattr(self.params.P, "target_dec_deg", None)
        if ra is not None and dec is not None:
            try:
                target_coord = SkyCoord(float(ra) * u.deg, float(dec) * u.deg)
            except Exception:
                target_coord = None

        self.worker = WcsWorker(
            file_list,
            self.params,
            self.params.P.data_dir,
            self.params.P.result_dir,
            self.params.P.cache_dir,
            self.use_cropped,
            target_coord=target_coord
        )
        self.worker.progress.connect(self.on_progress)
        self.worker.file_done.connect(self.on_file_done)
        self.worker.log_message.connect(self.log)
        self.worker.finished.connect(self.on_finished)
        self.worker.error.connect(self.on_error)
        self.setup_log_window()
        if hasattr(self, "_worker_panel") and self._worker_panel is not None:
            self._worker_panel.clear()
            self.worker.worker_status.connect(self._worker_panel.update_worker)

        self.run_bar_astap.set_running(True)
        self.progress_bar.setValue(0)
        self.progress_bar.setMaximum(len(file_list))
        self._wcs_start_time = time.monotonic()
        self.progress_label.setText(
            progress_status_text(0, len(file_list), self._wcs_start_time, message="Starting...")
        )
        self.worker.start()
        self.show_log_window()

    def stop_wcs(self):
        if self.worker and self.worker.isRunning():
            self.stop_requested = True
            self.run_bar_astap.set_stopping()
            self.progress_label.setText("Stopping...")
            self.log("Stop requested...")
            self.worker.stop()

    def on_progress(self, current, total, filename):
        self.progress_bar.setValue(current)
        self.progress_label.setText(
            progress_status_text(current, total, getattr(self, "_wcs_start_time", None), message=filename)
        )

    @staticmethod
    def _boolish(value) -> bool:
        if isinstance(value, (bool, np.bool_)):
            return bool(value)
        if value is None:
            return False
        if isinstance(value, (int, float, np.integer, np.floating)):
            return bool(value) and np.isfinite(value)
        return str(value).strip().lower() in {"1", "true", "t", "yes", "y", "ok", "solved"}

    @classmethod
    def _is_successful_wcs_result(cls, result: dict) -> bool:
        if not isinstance(result, dict):
            return False
        if cls._boolish(result.get("ok")):
            return True
        status = str(result.get("status", "")).strip().lower()
        return status in {"ok", "ok_astnet_wsl", "solved"}

    @classmethod
    def _wcs_result_row_background(cls, result: dict) -> str:
        ok = cls._is_successful_wcs_result(result)
        warning = False
        if ok and isinstance(result, dict):
            if "wcs_qc_pass" in result and not cls._boolish(result.get("wcs_qc_pass")):
                warning = True
            refine = str(result.get("refine", "") or "").strip().lower()
            if refine.startswith("match_too_small") or refine.startswith("error"):
                warning = True
        return status_row_background(ok, warning=warning)

    def _restore_success_results_from_summary(self):
        summary_path = step5_wcs_dir(self.params.P.result_dir) / "wcs_solve_summary.csv"
        if not summary_path.exists():
            return
        try:
            df = pd.read_csv(summary_path)
        except Exception:
            return
        if df.empty:
            return

        file_col = "file" if "file" in df.columns else "fname" if "fname" in df.columns else None
        if not file_col:
            return

        current_files = {Path(str(f)).name for f in getattr(self, "file_list", []) or []}
        restored = {}
        skipped_foreign = 0
        for row in df.to_dict("records"):
            if not self._is_successful_wcs_result(row):
                continue
            filename = row.get(file_col)
            if filename is None or pd.isna(filename):
                continue
            fname = Path(str(filename)).name
            if current_files and fname not in current_files:
                skipped_foreign += 1
                continue
            row["file"] = fname
            row.setdefault("fname", fname)
            restored[fname] = row

        if restored:
            self.results = restored
            counts = self._populate_restored_wcs_tables(restored)
            total = len(restored)
            parts = []
            if counts.get("internal"):
                parts.append(f"Internal {counts['internal']}")
                if hasattr(self, "internal_status"):
                    self.internal_status.setText(f"Loaded previous results: {counts['internal']} solved")
            if counts.get("astap"):
                parts.append(f"ASTAP {counts['astap']}")
            if counts.get("astnet"):
                parts.append(f"Astrometry.net {counts['astnet']}")
                if hasattr(self, "astrometrynet_status"):
                    self.astrometrynet_status.setText(f"Loaded previous results: {counts['astnet']} solved")
            if hasattr(self, "progress_label"):
                detail = ", ".join(parts) if parts else f"{total} solved"
                self.progress_label.setText(f"Loaded previous WCS summary: {detail}")
            try:
                self.log(f"[WCS][CACHE] Loaded previous WCS summary: {total} solved frame(s).")
            except Exception:
                pass
        elif skipped_foreign:
            try:
                self.log(f"[WCS][CACHE] Ignored {skipped_foreign} WCS row(s) not in the current file list.")
            except Exception:
                pass
        try:
            self._refresh_solver_breakdown_label()
        except Exception:
            pass

    def _get_wcs_cache_mgr(self) -> StepCacheManager:
        if self._wcs_cache_mgr is None:
            self._wcs_cache_mgr = StepCacheManager(
                self.params.P.cache_dir, "wcs_plate_solving", cache_schema_version=1
            )
        return self._wcs_cache_mgr

    def _write_wcs_manifest(self, filename: str, result: dict) -> None:
        try:
            from apex.utils.step_paths import step5_wcs_dir
            fits_path = (
                step2_cropped_dir(self.params.P.result_dir) / filename
                if self.use_cropped
                else Path(self.params.P.data_dir) / filename
            )
            wcs_out = step5_wcs_dir(self.params.P.result_dir)
            stem = Path(filename).stem
            mgr = self._get_wcs_cache_mgr()
            manifest = mgr.build_manifest(
                input_paths=[fits_path] if fits_path.exists() else [],
                payload_paths={"wcs_summary": wcs_out / "wcs_solve_summary.csv"},
                extra={"status": result.get("status", ""), "filename": filename},
            )
            mgr.write_manifest(filename, manifest)
        except Exception:
            pass

    @staticmethod
    def _wcs_failure_detail(result: dict, *, limit: int = 260) -> str:
        reason = str(result.get("fail_reason", "") or "").strip()
        for key in (
            "astap_stderr",
            "astap_stdout",
            "astnet_wsl_stderr",
            "astnet_wsl_stdout",
        ):
            tail = _tail_text(result.get(key), limit=limit, max_lines=2)
            if tail:
                return f"{reason} | {tail}" if reason else tail
        return reason

    @staticmethod
    def _wcs_failure_tooltip(result: dict) -> str:
        parts: list[str] = []
        reason = str(result.get("fail_reason", "") or "").strip()
        if reason:
            parts.append(f"reason: {reason}")
        for label, key in (
            ("ASTAP stderr", "astap_stderr"),
            ("ASTAP stdout", "astap_stdout"),
            ("astrometry.net stderr", "astnet_wsl_stderr"),
            ("astrometry.net stdout", "astnet_wsl_stdout"),
        ):
            tail = _tail_text(result.get(key), limit=900, max_lines=5)
            if tail:
                parts.append(f"{label}: {tail}")
        return "\n".join(parts)

    def on_file_done(self, filename, result):
        if self._is_successful_wcs_result(result):
            self.results[filename] = result
        else:
            self.results.pop(filename, None)
        try:
            self._refresh_solver_breakdown_label()
        except Exception:
            pass
        row = self.results_table.rowCount()
        self.results_table.insertRow(row)
        self.results_table.setItem(row, 0, QTableWidgetItem(filename))
        status_item = QTableWidgetItem(str(result.get("status", "")))
        if not self._is_successful_wcs_result(result):
            tooltip = self._wcs_failure_tooltip(result)
            if tooltip:
                status_item.setToolTip(tooltip)
        self.results_table.setItem(row, 1, status_item)
        pix_fit = result.get("pix_fit")
        pix_str = f"{pix_fit:.4f}" if isinstance(pix_fit, float) and np.isfinite(pix_fit) else "-"
        self.results_table.setItem(row, 2, QTableWidgetItem(pix_str))
        refine = result.get("refine", "")
        self.results_table.setItem(row, 3, QTableWidgetItem(str(refine)))
        resid_med = result.get("resid_med")
        if isinstance(resid_med, float) and np.isfinite(resid_med):
            resid_str = f"{resid_med:.3f}\""
        else:
            resid_px = result.get("resid_med_px")
            resid_str = f"{resid_px:.3f}px" if isinstance(resid_px, float) and np.isfinite(resid_px) else "-"
        self.results_table.setItem(row, 4, QTableWidgetItem(resid_str))
        elapsed = result.get("elapsed", 0.0)
        self.results_table.setItem(row, 5, QTableWidgetItem(f"{elapsed:.1f}"))
        set_table_row_background(self.results_table, row, self._wcs_result_row_background(result))
        self.results_table.scrollToBottom()
        pix_fit_log = result.get("pix_fit")
        pix_log = f"{pix_fit_log:.4f}" if isinstance(pix_fit_log, float) and np.isfinite(pix_fit_log) else "-"
        resid_log = result.get("resid_med")
        if isinstance(resid_log, float) and np.isfinite(resid_log):
            resid_str = f"{resid_log:.3f}\""
        else:
            resid_px = result.get("resid_med_px")
            resid_str = f"{resid_px:.3f}px" if isinstance(resid_px, float) and np.isfinite(resid_px) else "-"
        if self._is_successful_wcs_result(result):
            self.log(f"{filename}: {result.get('status', '')} pix={pix_log} refine={refine or '-'} resid_med={resid_str}")
            self._write_wcs_manifest(filename, result)
        else:
            detail = self._wcs_failure_detail(result)
            detail_txt = f" reason={detail}" if detail else ""
            self.log(
                f"{filename}: {result.get('status', '')} pix={pix_log} "
                f"refine={refine or '-'} resid_med={resid_str}{detail_txt}"
            )

    def on_error(self, filename, error):
        self.log(f"ERROR {filename}: {error}")

    def on_finished(self, summary):
        self.run_bar_astap.set_running(False)
        self.stop_requested = False
        stopped = bool(summary.get("stopped")) if isinstance(summary, dict) else False
        elapsed_txt = ""
        if hasattr(self, "_wcs_start_time"):
            elapsed_txt = f" | elapsed {format_duration(time.monotonic() - self._wcs_start_time)}"
        self.progress_label.setText(("Stopped" if stopped else "Done") + elapsed_txt)
        if summary:
            qc_not_eval = int(summary.get("wcs_qc_not_evaluated", 0) or 0)
            qc_text = (
                f"WCS-QC not evaluated: {qc_not_eval} (Gaia unavailable)"
                if qc_not_eval
                else f"WCS-QC pass: {summary.get('wcs_qc_pass', 0)}"
            )
            self.log(f"WCS done: {summary.get('ok', 0)}/{summary.get('total', 0)} OK | {qc_text}")
        self.save_state()
        self.update_navigation_buttons()

    def validate_step(self) -> bool:
        return len(self.results) > 0

    def save_state(self):
        state_data = {
            "wcs_complete": len(self.results) > 0,
            "n_files": len(self.results),
            "use_cropped": self.use_cropped,
            "astap_exe": getattr(self.params.P, "astap_exe", "astap_cli.exe"),
            "astap_timeout_s": getattr(self.params.P, "astap_timeout_s", 120.0),
            "astap_search_radius_deg": getattr(self.params.P, "astap_search_radius_deg", 8.0),
            "astap_database": getattr(self.params.P, "astap_database", "D80"),
            "astap_annotate_variables": getattr(self.params.P, "astap_annotate_variables", False),
            "astap_fov_fudge": getattr(self.params.P, "astap_fov_fudge", 1.0),
            "astap_downsample_z": getattr(self.params.P, "astap_downsample_z", 2),
            "astap_max_stars_s": getattr(self.params.P, "astap_max_stars_s", 500),
            "wcs_max_workers": getattr(self.params.P, "wcs_max_workers", 1),
            "wcs_require_qc_pass": getattr(self.params.P, "wcs_require_qc_pass", True),
            "wcs_refine_enable": getattr(self.params.P, "wcs_refine_enable", True),
            "wcs_refine_max_match": getattr(self.params.P, "wcs_refine_max_match", 600),
            "wcs_refine_match_r_fwhm": getattr(self.params.P, "wcs_refine_match_r_fwhm", 1.6),
            "wcs_refine_min_match": getattr(self.params.P, "wcs_refine_min_match", 50),
            "gaia_radius_fudge": getattr(self.params.P, "gaia_radius_fudge", 1.35),
            "gaia_mag_max": getattr(self.params.P, "gaia_mag_max", 18.0),
            "gaia_wcs_mag_max": getattr(self.params.P, "gaia_wcs_mag_max", 18.0),
            "gaia_retry": getattr(self.params.P, "gaia_retry", 2),
            "gaia_timeout_s": getattr(self.params.P, "gaia_timeout_s", 30.0),
            "gaia_backoff_s": getattr(self.params.P, "gaia_backoff_s", 6.0),
            "gaia_allow_no_cache": getattr(self.params.P, "gaia_allow_no_cache", True),
            "internal_params": dict(getattr(self, "_internal_params", {}) or {}),
        }
        self.project_state.store_step_data("wcs_plate_solve", state_data)

    def restore_state(self):
        state_data = self.project_state.get_step_data("wcs_plate_solve")
        if state_data:
            for key, val in state_data.items():
                if hasattr(self.params.P, key):
                    setattr(self.params.P, key, val)
            internal_params = state_data.get("internal_params")
            if isinstance(internal_params, dict) and hasattr(self, "_internal_params"):
                for key, val in internal_params.items():
                    if key in self._internal_params:
                        self._internal_params[key] = val
        self._restore_success_results_from_summary()
        self.update_navigation_buttons()

    def closeEvent(self, event):
        workers = [
            (getattr(self, "worker", None), self.stop_wcs),
            (getattr(self, "astrometrynet_worker", None), self.stop_astrometrynet_solve),
            (getattr(self, "_internal_worker", None), self.stop_wcs_internal_solver),
        ]
        for worker, stop_fn in workers:
            if worker is not None and worker.isRunning():
                stop_fn()

        for worker, _ in workers:
            if worker is not None and worker.isRunning() and not worker.wait(10000):
                QMessageBox.warning(
                    self,
                    "Background Task Running",
                    "A WCS solver is still stopping. Please wait and close again.",
                )
                event.ignore()
                return
        super().closeEvent(event)
