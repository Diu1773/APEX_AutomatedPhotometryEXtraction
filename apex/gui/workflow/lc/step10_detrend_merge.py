"""
Step 11: Detrend & Night Merge
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime
import json
import re

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from PyQt5.QtWidgets import (
    QApplication,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QLabel,
    QDialog,
    QGroupBox,
    QLineEdit,
    QMessageBox,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QTextEdit,
    QWidget,
    QFormLayout,
    QDoubleSpinBox,
    QSpinBox,
    QCheckBox,
    QComboBox,
    QListWidget,
    QListWidgetItem,
    QRadioButton,
    QButtonGroup,
    QTabWidget,
    QSplitter,
    QFrame,
    QTextBrowser,
    QProgressBar,
    QScrollArea,
)
from PyQt5.QtCore import Qt, QSignalBlocker, QThread, QTimer, pyqtSignal
from PyQt5.QtGui import QColor, QFont

from apex.gui.layout_rules import FittedDialog, scroll_wrap
from apex.gui.theme import Tokens, readable_on, refresh, style_button
from apex.gui.workflow.step_window_base import StepWindowBase
from apex.analysis.light_curve.global_ensemble import solve_global_ensemble
from apex.analysis.light_curve.photometry_source_service import (
    load_lightcurve_frame_photometry,
    resolve_lightcurve_photometry_source,
)
from apex.analysis.light_curve.detrend_output_service import (
    annotate_step10_output,
    build_detrend_summary_report_text,
    build_global_summary_text,
    write_step10_current_meta,
)
from apex.analysis.light_curve.detrend_runner import (
    DetrendRunner,
    _load_headers_table,
    _load_night_assignments_from_disk,
    _load_step8_source_to_id_map,
    _load_step9_comparison_ids_by_filter,
    _load_step9_selection_ids,
    _parse_color_expr,
)
from apex.utils.step_paths_lc import (
    step1_dir,
    step8_selection_dir,
    step9_lc_dir,
    step10_detrend_dir,
    step10_current_lc_path,
    step10_current_params_path,
    step10_current_summary_path,
    step10_current_plot_path,
    step10_current_meta_path,
    step10_current_global_zp_path,
    step10_current_global_mean_path,
    step10_current_global_diag_path,
    step10_history_dir,
    load_detrend_preference,
    list_lightcurve_csvs,
)
from apex.utils.step_paths import forced_phot_input_dir
from apex.utils.common_helpers import safe_float as _safe_float, normalize_filter_key as _normalize_filter_key, parse_jd as _parse_jd
from apex.utils.photometry_provenance import (
    build_photometry_provenance,
    format_photometry_provenance,
    summarize_photometry_table,
)
from apex.utils.io_utils import (
    read_csv_int64_source_id,
    coerce_int64_source_id,
    load_night_assignments as _load_night_assignments_util,
    load_headers_table as _load_headers_table_util,
)
from apex.utils.qc_utils import (
    filter_frame_df_by_qc,
    load_frame_excludes,
    should_use_frame_quality_qc,
)
from apex.utils.astro_utils import (
    compute_airmass_from_jd_array,
    compute_airmass_from_header,
    compute_bjd_tdb_array,
    is_reasonable_airmass,
)


class _DetrendTaskWorker(QThread):
    result_ready = pyqtSignal(str, object)
    error = pyqtSignal(str, str)

    def __init__(self, owner, kind: str, fn):
        super().__init__(owner)
        self.kind = str(kind)
        self._fn = fn

    def stop(self):
        self.requestInterruption()

    def run(self):
        try:
            result = self._fn()
            if not self.isInterruptionRequested():
                self.result_ready.emit(self.kind, result)
        except Exception as exc:
            if not self.isInterruptionRequested():
                self.error.emit(self.kind, str(exc))


def _wrap_scroll(widget: QWidget) -> QScrollArea:
    """Frameless scroll wrapper for a tab page.

    QTabWidget's minimum height is its tallest page's — the mode/options
    stacks alone pushed the whole window past small screens. Wrapped pages
    scroll instead of dictating window height.

    Horizontal scrolling stays enabled: these pages sit in a 384 px column but
    their option rows want ~455 px, and with the horizontal bar switched off
    that overflow was clipped with no way to reach it.
    """
    return scroll_wrap(widget, horizontal=True)


def _fmt_float(value, default: str = "") -> str:
    try:
        if value is None:
            return default
        v = float(value)
        if not np.isfinite(v):
            return default
        return f"{v:.5f}"
    except Exception:
        return default


def _load_check_star_for_plot(result_dir: Path, filt: str | None = None):
    """Load check star CSV from lc_lightcurve/ for plotting. Returns (check_id, df_or_None)."""
    try:
        from apex.analysis.light_curve.check_star_io import load_check_star_csv
        check_id, df = load_check_star_csv(result_dir, filt=filt)
        return check_id, (df if not df.empty else None)
    except Exception:
        return None, None


class DetrendNightMergeWindow(StepWindowBase, DetrendRunner):
    """Step 11: Nightly detrend + merge"""

    background_log = pyqtSignal(str)

    # -- the seam: what the batch path supplies itself -----------------------

    def _target_id_text(self) -> str:
        return self.target_edit.text().strip()

    def _filter_selection(self) -> str:
        return self.filter_combo.currentText() or "All"

    def _use_global_k2(self) -> bool:
        return bool(self.chk_global_k2.isChecked())

    def _ensure_plot_drawn(self) -> None:
        """The canvas already holds the current view; saving must not redraw it."""

    def _plot_figure(self):
        return self.plot_canvas.figure

    def _plot_redraw(self) -> None:
        self.plot_canvas.draw()

    def _refresh_style(self, widget=None) -> None:
        if widget is not None:
            refresh(widget)

    def _tell_user(self, level: str, title: str, message: str) -> None:
        if level == "warn":
            QMessageBox.warning(self, title, message)
        elif level == "error":
            QMessageBox.critical(self, title, message)
        else:
            QMessageBox.information(self, title, message)

    def __init__(self, params, file_manager, project_state, main_window, runtime_mode: bool = False):
        self.file_manager = file_manager
        self.runtime_mode = bool(runtime_mode)
        self.datasets: list[tuple[str, Path]] = []
        self.raw_df = pd.DataFrame()
        self.corrected_df = pd.DataFrame()
        self.params_df = pd.DataFrame()
        self.global_mean_df = pd.DataFrame()
        self.global_diagnostics: dict = {}
        self.global_input_df = pd.DataFrame()
        self._busy_active = False
        self._busy_log_buffer: list[str] | None = None
        self._busy_message = ""
        self._detrend_worker: _DetrendTaskWorker | None = None
        self._pending_fit_request: dict | None = None

        self.comp_active_ids: list[int] = []
        self.comp_candidate_ids: list[int] = []
        self.color_map_by_filter: dict[str, str] = {}
        self.color_by = "Date"
        self.mode = "offset"  # "offset" | "color" | "global" | "sysrem"
        self.sigma_clip = True
        self.clip_sigma = 3.0
        self.clip_iters = 2
        self.x_axis_mode = "time"
        self.phase_cycles = 2.0  # 기본 2주기 표시
        # Plot view: one big plot by default; "all" = classic 3-stack.
        self._plot_view_mode = "corr"

        self.delta_c_map: dict[str, float] = {}

        self.time_masks: list[tuple[float, float]] = []
        self.phase_masks: list[tuple[float, float]] = []
        self.phase_period = 0.0
        self.phase_t0 = 0.0

        # Global ensemble (method C) params
        self.global_min_comps = 3
        self.global_sigma = 3.0
        self.global_iters = 3
        self.global_rms_pct = 10.0
        self.global_rms_threshold = 0.0
        self.global_frame_sigma = 3.0
        self.global_gauge = "meanZ0"
        self.global_robust = True
        self.global_interp_missing = False
        self.global_normalize = False
        self.global_rescale_errors = True
        self.sysrem_iter = 5
        self.sysrem_apply = 3

        super().__init__(
            step_index=10,
            step_name="Detrend & Night Merge",
            params=params,
            project_state=project_state,
            main_window=main_window,
        )
        if self.runtime_mode:
            self._setup_runtime_ui()
            self.restore_state()
        else:
            self.setup_step_ui()
            self.restore_state()
            self._auto_load_ids()
            self.load_raw_data(
                silent=True,
                rebuild_missing=False,
                allow_fits_airmass=False,
            )
        self.background_log.connect(self._append_background_log)
        if not self.runtime_mode:
            QTimer.singleShot(0, self._start_airmass_backfill_if_needed)

    def _setup_runtime_ui(self):
        """Build a lightweight runtime-only UI for merger inline execution."""
        note = QLabel("Runtime mode: merger inline Step11 execution")
        note.setProperty("role", "caption")
        self.content_layout.addWidget(note)

        rd = Path(self.params.P.result_dir)
        self.datasets = [(rd.name, rd)]

        self.target_edit = QLineEdit()
        self.id_info_label = QLabel("Target / Comp: (runtime)")
        self.id_info_label.setWordWrap(True)
        self.date_list = QListWidget()
        self.filter_combo = QComboBox()
        self.filter_combo.addItem("All")
        self.color_by_combo = QComboBox()
        self.color_by_combo.addItems(["Date", "Filter"])
        self.analysis_text = QLabel("Runtime analysis")
        self.analysis_text.setWordWrap(True)
        self.recommendation_label = QLabel("")
        self.recommendation_label.setWordWrap(True)
        self.color_status_label = QLabel("")
        self.color_map_group = QGroupBox("Color Index")
        self.left_tabs = QTabWidget()

        self.mode_offset = QRadioButton("Offset")
        self.mode_color = QRadioButton("Color")
        self.mode_global = QRadioButton("Global")
        self.mode_offset.setChecked(True)
        self.chk_global_k2 = QCheckBox("Global k''")
        self.chk_global_k2.setChecked(True)

        self.phase_mode_combo = QComboBox()
        self.phase_mode_combo.addItems(["Time (JD)", "Phase"])
        self.spin_period = QDoubleSpinBox()
        self.spin_period.setDecimals(6)
        self.spin_period.setRange(0.0, 1000.0)
        self.spin_t0 = QDoubleSpinBox()
        self.spin_t0.setDecimals(6)
        self.spin_t0.setRange(0.0, 3000000.0)
        self.spin_cycles = QDoubleSpinBox()
        self.spin_cycles.setDecimals(2)
        self.spin_cycles.setRange(1.0, 5.0)
        self.spin_cycles.setValue(self.phase_cycles)

        self.chk_clip = QCheckBox()
        self.chk_clip.setChecked(self.sigma_clip)
        self.spin_clip = QDoubleSpinBox()
        self.spin_clip.setDecimals(1)
        self.spin_clip.setRange(1.0, 10.0)
        self.spin_clip.setValue(self.clip_sigma)
        self.spin_iters = QSpinBox()
        self.spin_iters.setRange(1, 5)
        self.spin_iters.setValue(self.clip_iters)

        self.spin_global_min_comps = QSpinBox()
        self.spin_global_min_comps.setRange(1, 50)
        self.spin_global_min_comps.setValue(self.global_min_comps)
        self.spin_global_sigma = QDoubleSpinBox()
        self.spin_global_sigma.setRange(1.0, 10.0)
        self.spin_global_sigma.setDecimals(1)
        self.spin_global_sigma.setValue(self.global_sigma)
        self.spin_global_iters = QSpinBox()
        self.spin_global_iters.setRange(1, 5)
        self.spin_global_iters.setValue(self.global_iters)
        self.spin_global_rms_pct = QDoubleSpinBox()
        self.spin_global_rms_pct.setRange(0.0, 80.0)
        self.spin_global_rms_pct.setDecimals(1)
        self.spin_global_rms_pct.setValue(self.global_rms_pct)
        self.spin_global_rms_pct.setToolTip(
            "Maximum fraction of statistically high-RMS comparison stars removed per iteration."
        )
        self.spin_global_rms_thr = QDoubleSpinBox()
        self.spin_global_rms_thr.setRange(0.0, 1.0)
        self.spin_global_rms_thr.setDecimals(4)
        self.spin_global_rms_thr.setValue(self.global_rms_threshold)
        self.spin_global_frame_sigma = QDoubleSpinBox()
        self.spin_global_frame_sigma.setRange(1.0, 10.0)
        self.spin_global_frame_sigma.setDecimals(1)
        self.spin_global_frame_sigma.setValue(self.global_frame_sigma)
        self.combo_global_gauge = QComboBox()
        self.combo_global_gauge.addItems(["meanZ0", "ref"])
        self.combo_global_gauge.setCurrentText(self.global_gauge)
        self.chk_global_robust = QCheckBox("Robust (MAD)")
        self.chk_global_robust.setChecked(self.global_robust)
        self.chk_global_interp = QCheckBox("Z_t 보간")
        self.chk_global_interp.setChecked(self.global_interp_missing)
        self.chk_global_normalize = QCheckBox("Target 중앙값 0")
        self.chk_global_normalize.setChecked(self.global_normalize)
        self.chk_global_rescale_err = QCheckBox("Chi²_red 오차 보정")
        self.chk_global_rescale_err.setChecked(self.global_rescale_errors)

        self.btn_apply = QPushButton("Fit && Apply")
        self.btn_revert = QPushButton("Revert")
        self.busy_status_label = QLabel("")
        self.busy_progress_bar = QProgressBar()
        self.busy_progress_bar.setRange(0, 0)
        self.busy_progress_bar.hide()

        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)

        self.plot_canvas = FigureCanvas(Figure(figsize=(9, 7)))
        self.ax_raw = self.plot_canvas.figure.add_subplot(311)
        self.ax_corr = self.plot_canvas.figure.add_subplot(312)
        self.ax_diag = self.plot_canvas.figure.add_subplot(313)

    def setup_step_ui(self):
        # Main horizontal splitter
        main_splitter = QSplitter(Qt.Horizontal)
        self.content_layout.addWidget(main_splitter, 1)

        # ===== LEFT PANEL =====
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(6)

        # Info banner (compact)
        info = QLabel("여러 밤의 차등측광 데이터를 보정하여 병합합니다.")
        info.setProperty("role", "info")
        info.setWordWrap(True)
        left_layout.addWidget(info)

        # Left tabs — pane border comes from the themed QSS
        self.left_tabs = QTabWidget()
        left_layout.addWidget(self.left_tabs, 1)

        # ----- Tab 1: Data & Target -----
        data_tab = QWidget()
        data_layout = QVBoxLayout(data_tab)
        data_layout.setSpacing(8)

        # Hidden QLineEdit (used internally by load logic)
        self.target_edit = QLineEdit()

        # Fix base dataset to current result_dir
        rd = Path(self.params.P.result_dir)
        self.datasets = [(rd.name, rd)]

        # Read-only info bar (auto-loaded from Step 10 / Step 9)
        self.id_info_label = QLabel("Target / Comp: (loading...)")
        self.id_info_label.setProperty("banner", "warn")
        self.id_info_label.setWordWrap(True)
        data_layout.addWidget(self.id_info_label)

        source_info = resolve_lightcurve_photometry_source(
            rd, self.project_state
        )
        self.photometry_source_label = QLabel(
            format_photometry_provenance(source_info)
        )
        self.photometry_source_label.setProperty("role", "caption")
        self.photometry_source_label.setToolTip(str(source_info.get("reason", "")))
        data_layout.addWidget(self.photometry_source_label)

        # Date & Filter selection
        selection_row = QHBoxLayout()

        date_group = QGroupBox("날짜")
        date_layout = QVBoxLayout(date_group)
        date_layout.setContentsMargins(4, 4, 4, 4)
        self.date_list = QListWidget()
        # No max height: the night list fills the tab (the 120px cap left the
        # bottom half of the panel as dead space with many nights hidden).
        self.date_list.itemChanged.connect(self._on_date_selection_changed)
        date_layout.addWidget(self.date_list)
        selection_row.addWidget(date_group)

        filter_group = QGroupBox("필터 / 색상")
        filter_layout = QFormLayout(filter_group)
        filter_layout.setSpacing(4)
        self.filter_combo = QComboBox()
        self.filter_combo.addItem("All")
        self.filter_combo.currentIndexChanged.connect(self._on_filter_changed)
        filter_layout.addRow("필터:", self.filter_combo)
        self.color_by_combo = QComboBox()
        self.color_by_combo.addItems(["Date", "Filter"])
        self.color_by_combo.currentIndexChanged.connect(self._on_color_by_changed)
        filter_layout.addRow("Color by:", self.color_by_combo)
        selection_row.addWidget(filter_group)

        data_layout.addLayout(selection_row, 1)
        self.left_tabs.addTab(data_tab, "데이터")

        # ----- Tab 2: Correction Mode -----
        mode_tab = QWidget()
        mode_layout = QVBoxLayout(mode_tab)
        mode_layout.setSpacing(8)

        # Analysis feedback panel
        self.analysis_panel = QGroupBox("데이터 분석 결과")
        self.analysis_panel.setStyleSheet(
            "QGroupBox { font-weight: bold; } "
            f"QGroupBox::title {{ color: {Tokens.ACCENT}; }}"
        )
        analysis_layout = QVBoxLayout(self.analysis_panel)
        analysis_layout.setSpacing(4)

        self.analysis_text = QLabel("Step 10 결과가 있으면 자동 로드되어 분석 결과가 표시됩니다.")
        self.analysis_text.setWordWrap(True)
        self.analysis_text.setStyleSheet(
            f"QLabel {{ background-color: {Tokens.SURFACE_ALT}; padding: 8px; "
            "border-radius: 4px; font-size: 9pt; }"
        )
        analysis_layout.addWidget(self.analysis_text)

        self.recommendation_label = QLabel("")
        self.recommendation_label.setWordWrap(True)
        self.recommendation_label.setProperty("banner", "ok")
        analysis_layout.addWidget(self.recommendation_label)
        mode_layout.addWidget(self.analysis_panel)

        # Mode selection
        mode_group = QGroupBox("보정 모드 선택")
        mode_group_layout = QVBoxLayout(mode_group)
        mode_group_layout.setSpacing(6)

        mode_header = QHBoxLayout()
        mode_header.setSpacing(6)
        mode_help_label = QLabel("각 모드의 입력 데이터, 보정식, 결과 설명")
        mode_help_label.setStyleSheet(f"QLabel {{ color: {Tokens.TEXT_SUB}; font-size: 8pt; }}")
        mode_header.addWidget(mode_help_label)
        mode_header.addStretch()
        btn_mode_help = QPushButton("?")
        btn_mode_help.setFixedSize(24, 24)
        btn_mode_help.setToolTip("보정 모드 도움말")
        btn_mode_help.setStyleSheet(
            f"QPushButton {{ font-weight: bold; border: 1px solid {Tokens.BORDER_STRONG}; "
            f"border-radius: 12px; background: {Tokens.SURFACE}; }}"
        )
        btn_mode_help.clicked.connect(self._show_mode_help_dialog)
        mode_header.addWidget(btn_mode_help)
        mode_group_layout.addLayout(mode_header)

        self.mode_group = QButtonGroup(self)

        # Offset mode
        offset_frame = QFrame()
        offset_frame.setStyleSheet(
            f"QFrame {{ border: 1px solid {Tokens.BORDER}; "
            "border-radius: 4px; padding: 4px; }"
        )
        offset_layout = QVBoxLayout(offset_frame)
        offset_layout.setSpacing(2)
        self.mode_offset = QRadioButton("Offset Only (ZP₀)")
        _f = self.mode_offset.font(); _f.setBold(True); self.mode_offset.setFont(_f)
        offset_layout.addWidget(self.mode_offset)
        offset_desc = QLabel(
            "Δm_corr = Δm_raw - ZP₀\n"
            "• 밤별 영점 오프셋만 보정\n"
            "• Target-Comp 색차가 작을 때 (|ΔC| < 0.3)"
        )
        offset_desc.setStyleSheet(f"QLabel {{ color: {Tokens.TEXT_SUB}; font-size: 8pt; margin-left: 16px; }}")
        offset_layout.addWidget(offset_desc)
        mode_group_layout.addWidget(offset_frame)

        # Color mode
        color_frame = QFrame()
        color_frame.setStyleSheet(
            f"QFrame {{ border: 1px solid {Tokens.BORDER}; "
            "border-radius: 4px; padding: 4px; }"
        )
        color_layout = QVBoxLayout(color_frame)
        color_layout.setSpacing(2)
        self.mode_color = QRadioButton("Color-dependent (ZP₀ + k''·ΔC·X)")
        _f = self.mode_color.font(); _f.setBold(True); self.mode_color.setFont(_f)
        color_layout.addWidget(self.mode_color)
        color_desc = QLabel(
            "Δm_corr = Δm_raw - ZP₀ - k''·ΔC·X\n"
            "• 2차 소광계수(k'')로 색차 보정\n"
            "• Target-Comp 색차가 클 때 (|ΔC| ≥ 0.3)"
        )
        color_desc.setStyleSheet(f"QLabel {{ color: {Tokens.TEXT_SUB}; font-size: 8pt; margin-left: 16px; }}")
        color_layout.addWidget(color_desc)

        self.chk_global_k2 = QCheckBox("Global k'' (전체 데이터로 k'' 한 번 피팅)")
        self.chk_global_k2.setStyleSheet(f"QCheckBox {{ margin-left: {Tokens.MARGIN}px; }}")
        self.chk_global_k2.setChecked(True)
        color_layout.addWidget(self.chk_global_k2)
        mode_group_layout.addWidget(color_frame)

        # Global ensemble mode (method C)
        global_frame = QFrame()
        global_frame.setStyleSheet(
            f"QFrame {{ border: 1px solid {Tokens.BORDER}; "
            "border-radius: 4px; padding: 4px; }"
        )
        global_layout = QVBoxLayout(global_frame)
        global_layout.setSpacing(2)
        self.mode_global = QRadioButton("Global Ensemble (Method C)")
        _f = self.mode_global.font(); _f.setBold(True); self.mode_global.setFont(_f)
        global_layout.addWidget(self.mode_global)
        global_desc = QLabel(
            "Δm_corr = (m_target - Z_t) - <M_comp>\n"
            "• 프레임별 Z_t를 전역 최소제곱으로 동시 추정\n"
            "• comparison ensemble 기준 차등광도 스케일 유지"
        )
        global_desc.setStyleSheet(f"QLabel {{ color: {Tokens.TEXT_SUB}; font-size: 8pt; margin-left: 16px; }}")
        global_layout.addWidget(global_desc)
        mode_group_layout.addWidget(global_frame)

        # SYSREM mode
        sysrem_frame = QFrame()
        sysrem_frame.setStyleSheet(
            f"QFrame {{ border: 1px solid {Tokens.BORDER}; "
            "border-radius: 4px; padding: 4px; }"
        )
        sysrem_layout = QVBoxLayout(sysrem_frame)
        sysrem_layout.setSpacing(2)
        self.mode_sysrem = QRadioButton("SYSREM (Tamuz+2005)")
        _f = self.mode_sysrem.font(); _f.setBold(True); self.mode_sysrem.setFont(_f)
        sysrem_layout.addWidget(self.mode_sysrem)
        sysrem_desc = QLabel(
            "비교성 집합에서 공통 체계오차 벡터를 반복 추출하여 제거\n"
            "• 외부 영향(대기, 측기 변화)에 강인\n"
            "• Step 7 강제측광 데이터 직접 사용"
        )
        sysrem_desc.setStyleSheet(f"QLabel {{ color: {Tokens.TEXT_SUB}; font-size: 8pt; margin-left: 16px; }}")
        sysrem_layout.addWidget(sysrem_desc)
        sysrem_params_layout = QHBoxLayout()
        sysrem_params_layout.addWidget(QLabel("반복 수:"))
        self.spin_sysrem_iter = QSpinBox()
        self.spin_sysrem_iter.setRange(1, 20)
        self.spin_sysrem_iter.setValue(5)
        self.spin_sysrem_iter.setToolTip("추출할 체계오차 성분 수 (보통 3–7)")
        sysrem_params_layout.addWidget(self.spin_sysrem_iter)
        sysrem_params_layout.addWidget(QLabel("적용 수:"))
        self.spin_sysrem_apply = QSpinBox()
        self.spin_sysrem_apply.setRange(1, 20)
        self.spin_sysrem_apply.setValue(3)
        self.spin_sysrem_apply.setToolTip("타겟에 적용할 성분 수 (≤ 반복 수, 너무 많으면 신호 제거됨)")
        sysrem_params_layout.addWidget(self.spin_sysrem_apply)
        sysrem_params_layout.addStretch()
        sysrem_layout.addLayout(sysrem_params_layout)
        mode_group_layout.addWidget(sysrem_frame)

        self.mode_group.addButton(self.mode_offset)
        self.mode_group.addButton(self.mode_color)
        self.mode_group.addButton(self.mode_global)
        self.mode_group.addButton(self.mode_sysrem)
        self.mode_offset.setChecked(True)
        self.mode_offset.toggled.connect(lambda checked: checked and self._set_mode("offset"))
        self.mode_color.toggled.connect(lambda checked: checked and self._set_mode("color"))
        self.mode_global.toggled.connect(lambda checked: checked and self._set_mode("global"))
        self.mode_sysrem.toggled.connect(lambda checked: checked and self._set_mode("sysrem"))

        self.color_status_label = QLabel("")
        self.color_status_label.setStyleSheet(f"QLabel {{ color: {Tokens.ERROR}; font-size: 9pt; }}")
        self.color_status_label.setWordWrap(True)
        mode_group_layout.addWidget(self.color_status_label)

        mode_layout.addWidget(mode_group)

        # Color index mapping (compact)
        self.color_map_group = QGroupBox("Color Index (ΔC) 설정")
        self.color_map_layout = QFormLayout(self.color_map_group)
        self.color_map_layout.setSpacing(4)
        self.color_map_combos = {}
        mode_layout.addWidget(self.color_map_group)

        mode_layout.addStretch()
        # Scroll-wrapped: QTabWidget's minimum height is the TALLEST tab's —
        # this stack of mode groups alone pushed the window past 1200px.
        self.left_tabs.addTab(_wrap_scroll(mode_tab), "보정 모드")

        # ----- Tab 3: Phase & Options -----
        options_tab = QWidget()
        options_layout = QVBoxLayout(options_tab)
        options_layout.setSpacing(8)

        # Phase folding
        phase_group = QGroupBox("Phase Folding")
        phase_layout = QFormLayout(phase_group)
        phase_layout.setSpacing(4)

        self.phase_mode_combo = QComboBox()
        self.phase_mode_combo.addItems(["Time (JD)", "Phase"])
        self.phase_mode_combo.currentIndexChanged.connect(self._on_xaxis_changed)
        phase_layout.addRow("X축:", self.phase_mode_combo)

        self.spin_period = QDoubleSpinBox()
        self.spin_period.setDecimals(6)
        self.spin_period.setRange(0.0, 1000.0)
        self.spin_period.valueChanged.connect(self._on_phase_params_changed)
        phase_layout.addRow("주기 (days):", self.spin_period)

        t0_row = QHBoxLayout()
        self.spin_t0 = QDoubleSpinBox()
        self.spin_t0.setDecimals(6)
        self.spin_t0.setRange(0.0, 3000000.0)
        self.spin_t0.setToolTip("기준 시점 (예: 주극소 JD). 0이면 min(JD) 자동 사용")
        self.spin_t0.valueChanged.connect(self._on_phase_params_changed)
        t0_row.addWidget(self.spin_t0)
        btn_t0_auto = QPushButton("Auto")
        btn_t0_auto.setMaximumWidth(50)
        btn_t0_auto.setToolTip("데이터의 min(JD)로 설정")
        btn_t0_auto.clicked.connect(self._auto_set_t0)
        t0_row.addWidget(btn_t0_auto)
        phase_layout.addRow("T₀ (JD):", t0_row)

        t0_note = QLabel("※ T₀ = 기준시점 (주극소/극대). 0이면 min(JD) 사용")
        t0_note.setStyleSheet(f"QLabel {{ color: {Tokens.TEXT_MUTED}; font-size: 8pt; }}")
        phase_layout.addRow("", t0_note)

        self.spin_cycles = QDoubleSpinBox()
        self.spin_cycles.setDecimals(2)
        self.spin_cycles.setRange(1.0, 5.0)
        self.spin_cycles.setSingleStep(0.5)
        self.spin_cycles.setValue(2.0)  # 기본 2주기
        self.spin_cycles.valueChanged.connect(self._on_phase_params_changed)
        phase_layout.addRow("표시 주기:", self.spin_cycles)
        options_layout.addWidget(phase_group)

        # Sigma clipping
        clip_group = QGroupBox("Sigma Clipping (이상치 제거)")
        clip_layout = QVBoxLayout(clip_group)
        clip_layout.setSpacing(4)

        clip_desc = QLabel(
            "피팅 시 잔차가 Nσ 이상인 이상치를 반복 제거합니다.\n"
            "(구름, 장비 오류 등으로 인한 outlier 제거용)"
        )
        clip_desc.setStyleSheet(f"QLabel {{ color: {Tokens.TEXT_SUB}; font-size: 8pt; }}")
        clip_desc.setWordWrap(True)
        clip_layout.addWidget(clip_desc)

        clip_form = QFormLayout()
        clip_form.setSpacing(4)

        self.chk_clip = QCheckBox()
        self.chk_clip.setChecked(True)
        self.chk_clip.stateChanged.connect(self._on_clip_changed)
        clip_form.addRow("활성화:", self.chk_clip)

        self.spin_clip = QDoubleSpinBox()
        self.spin_clip.setDecimals(1)
        self.spin_clip.setRange(1.0, 10.0)
        self.spin_clip.setValue(self.clip_sigma)
        self.spin_clip.setToolTip("이 값 × σ 이상의 잔차를 가진 점을 제거")
        self.spin_clip.valueChanged.connect(self._on_clip_changed)
        clip_form.addRow("σ 임계값:", self.spin_clip)

        self.spin_iters = QSpinBox()
        self.spin_iters.setRange(1, 5)
        self.spin_iters.setValue(self.clip_iters)
        self.spin_iters.setToolTip("클리핑 후 재피팅 반복 횟수")
        self.spin_iters.valueChanged.connect(self._on_clip_changed)
        clip_form.addRow("반복 횟수:", self.spin_iters)

        clip_layout.addLayout(clip_form)
        options_layout.addWidget(clip_group)

        # Global ensemble parameters
        global_group = QGroupBox("Global Ensemble (Method C)")
        global_layout = QFormLayout(global_group)
        global_layout.setSpacing(4)

        self.spin_global_min_comps = QSpinBox()
        self.spin_global_min_comps.setRange(1, 50)
        self.spin_global_min_comps.setValue(self.global_min_comps)
        global_layout.addRow("최소 비교성 수:", self.spin_global_min_comps)

        self.spin_global_sigma = QDoubleSpinBox()
        self.spin_global_sigma.setRange(1.0, 10.0)
        self.spin_global_sigma.setDecimals(1)
        self.spin_global_sigma.setValue(self.global_sigma)
        global_layout.addRow("σ 클립:", self.spin_global_sigma)

        self.spin_global_iters = QSpinBox()
        self.spin_global_iters.setRange(1, 5)
        self.spin_global_iters.setValue(self.global_iters)
        global_layout.addRow("반복 횟수:", self.spin_global_iters)

        self.spin_global_rms_pct = QDoubleSpinBox()
        self.spin_global_rms_pct.setRange(0.0, 80.0)
        self.spin_global_rms_pct.setDecimals(1)
        self.spin_global_rms_pct.setValue(self.global_rms_pct)
        self.spin_global_rms_pct.setToolTip(
            "Maximum fraction of statistically high-RMS comparison stars removed per iteration."
        )
        global_layout.addRow("RMS 상위% 제거:", self.spin_global_rms_pct)

        self.spin_global_rms_thr = QDoubleSpinBox()
        self.spin_global_rms_thr.setRange(0.0, 1.0)
        self.spin_global_rms_thr.setDecimals(4)
        self.spin_global_rms_thr.setValue(self.global_rms_threshold)
        self.spin_global_rms_thr.setToolTip("0이면 비활성")
        global_layout.addRow("RMS 절대컷:", self.spin_global_rms_thr)

        self.spin_global_frame_sigma = QDoubleSpinBox()
        self.spin_global_frame_sigma.setRange(1.0, 10.0)
        self.spin_global_frame_sigma.setDecimals(1)
        self.spin_global_frame_sigma.setValue(self.global_frame_sigma)
        global_layout.addRow("프레임 σ:", self.spin_global_frame_sigma)

        self.combo_global_gauge = QComboBox()
        self.combo_global_gauge.addItems(["meanZ0", "ref"])
        self.combo_global_gauge.setCurrentText(self.global_gauge)
        global_layout.addRow("Gauge:", self.combo_global_gauge)

        self.chk_global_robust = QCheckBox("Robust (MAD)")
        self.chk_global_robust.setChecked(self.global_robust)
        global_layout.addRow(self.chk_global_robust)

        self.chk_global_interp = QCheckBox("Z_t 보간 (부족 프레임)")
        self.chk_global_interp.setChecked(self.global_interp_missing)
        global_layout.addRow(self.chk_global_interp)

        self.chk_global_normalize = QCheckBox("Target 중앙값 0으로 정규화")
        self.chk_global_normalize.setChecked(self.global_normalize)
        global_layout.addRow(self.chk_global_normalize)

        self.chk_global_rescale_err = QCheckBox("Chi²_red 기반 오차 보정")
        self.chk_global_rescale_err.setChecked(self.global_rescale_errors)
        self.chk_global_rescale_err.setToolTip(
            "프레임별 chi²_red에 따라 오차를 보정합니다 (1.0 이상일 때만 확대, 축소하지 않음)"
        )
        global_layout.addRow(self.chk_global_rescale_err)

        options_layout.addWidget(global_group)

        options_layout.addStretch()
        self.left_tabs.addTab(_wrap_scroll(options_tab), "옵션")

        # ----- Tab 4: Log -----
        log_tab = QWidget()
        log_layout = QVBoxLayout(log_tab)
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setObjectName("Log")
        log_layout.addWidget(self.log_text)
        self.left_tabs.addTab(log_tab, "로그")

        # Action buttons (at bottom of left panel)
        btn_group = QFrame()
        btn_group.setObjectName("Card")  # themed flat surface
        btn_layout = QHBoxLayout(btn_group)
        btn_layout.setSpacing(8)

        # The step's single primary action (was a hand-painted green fill).
        self.btn_apply = QPushButton("Fit && Apply (저장)")
        style_button(self.btn_apply, "primary", height=Tokens.H_ACTION)
        self.btn_apply.setToolTip("피팅 수행 후 종합 결과 파일 자동 저장")
        self.btn_apply.clicked.connect(self.fit_and_apply)
        btn_layout.addWidget(self.btn_apply)

        self.btn_revert = QPushButton("Revert")
        style_button(self.btn_revert, height=Tokens.H_ACTION)
        self.btn_revert.setToolTip("보정 결과 초기화")
        self.btn_revert.clicked.connect(self.revert_raw)
        btn_layout.addWidget(self.btn_revert)

        self.busy_status_label = QLabel("")
        self.busy_status_label.setProperty("role", "subtitle")
        self.busy_status_label.hide()
        btn_layout.addWidget(self.busy_status_label)

        self.busy_progress_bar = QProgressBar()
        self.busy_progress_bar.setRange(0, 0)
        self.busy_progress_bar.setMaximumWidth(140)
        self.busy_progress_bar.setMaximumHeight(10)
        self.busy_progress_bar.setTextVisible(False)
        self.busy_progress_bar.setStyleSheet(
            f"QProgressBar {{ border: none; background-color: {Tokens.SURFACE_ALT}; border-radius: 4px; }}"
            f"QProgressBar::chunk {{ background-color: {Tokens.ACCENT}; border-radius: 4px; }}"
        )
        self.busy_progress_bar.hide()
        btn_layout.addWidget(self.busy_progress_bar)

        btn_layout.addStretch()
        left_layout.addWidget(btn_group)

        main_splitter.addWidget(left_widget)

        # ===== RIGHT PANEL =====
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(4)

        # View selector — one big plot with a segment switch; "All" restores
        # the classic 3-stack (the old always-on stack needed 1283px minimum
        # height and clipped on laptops).
        view_row = QHBoxLayout()
        view_row.setSpacing(Tokens.GAP)
        view_row.addWidget(QLabel("보기:"))
        self._plot_view_group = QButtonGroup(self)
        self._plot_view_group.setExclusive(True)
        self._plot_view_buttons = {}
        for key, label in (("raw", "Raw"), ("corr", "Corrected"),
                           ("diag", "Diagnostics"), ("all", "모두 보기")):
            btn = QPushButton(label)
            btn.setCheckable(True)
            style_button(btn, height=Tokens.H_COMPACT)
            btn.clicked.connect(lambda _=False, k=key: self._set_plot_view(k))
            self._plot_view_group.addButton(btn)
            self._plot_view_buttons[key] = btn
            view_row.addWidget(btn)
        view_row.addStretch()
        right_layout.addLayout(view_row)

        # Plot canvas
        self.plot_canvas = FigureCanvas(Figure(figsize=(9, 7)))
        self.ax_raw = self.plot_canvas.figure.add_subplot(311)
        self.ax_corr = self.plot_canvas.figure.add_subplot(312)
        self.ax_diag = self.plot_canvas.figure.add_subplot(313)
        right_layout.addWidget(self.plot_canvas, 3)

        # 피팅 결과 테이블은 좌측 "결과" 탭으로 이동 — 우측은 플롯 전용.
        self.result_table = QTableWidget()
        self.result_table.setColumnCount(9)
        self.result_table.setHorizontalHeaderLabels(
            ["Date", "Filter", "N", "ZP₀", "±σ(ZP)", "k''", "±σ(k'')", "RMS전", "RMS후"]
        )
        self.result_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        result_tab = QWidget()
        result_tab_layout = QVBoxLayout(result_tab)
        result_tab_layout.setContentsMargins(4, 4, 4, 4)
        result_tab_layout.addWidget(self.result_table)
        # Insert before the trailing 로그 tab.
        self.left_tabs.insertTab(self.left_tabs.count() - 1, result_tab, "결과")

        self._plot_view_buttons[self._plot_view_mode].setChecked(True)
        self._apply_plot_view(redraw=False)

        main_splitter.addWidget(right_widget)
        main_splitter.setSizes([320, 680])
        main_splitter.setChildrenCollapsible(False)

    def _set_plot_view(self, mode: str) -> None:
        if mode == getattr(self, "_plot_view_mode", "corr"):
            return
        self._plot_view_mode = mode
        self._apply_plot_view()
        self.save_state()

    def log(self, msg: str):
        if self._busy_log_buffer is not None:
            self._busy_log_buffer.append(str(msg))
            return
        if QThread.currentThread() != self.thread():
            self.background_log.emit(str(msg))
            return
        self._append_background_log(str(msg))

    def _set_busy_state(self, busy: bool, message: str = "") -> None:
        if self.runtime_mode:
            if busy:
                if self._busy_active:
                    self._busy_message = message or self._busy_message
                    return
                self._busy_active = True
                self._busy_log_buffer = []
                self._busy_message = message or "Working..."
                return
            if not self._busy_active:
                return
            self._flush_busy_log_buffer()
            self._busy_log_buffer = None
            self._busy_active = False
            self._busy_message = ""
            return

        if busy:
            if self._busy_active:
                self._set_busy_message(message)
                return
            self._busy_active = True
            self._busy_log_buffer = []
            self._busy_message = message or "Working..."
            self.busy_status_label.setText(self._busy_message)
            self.busy_status_label.show()
            self.busy_progress_bar.show()
            self.btn_apply.setEnabled(False)
            self.btn_revert.setEnabled(False)
            self.left_tabs.setEnabled(False)
            QApplication.setOverrideCursor(Qt.WaitCursor)
            QApplication.processEvents()
            return

        if not self._busy_active:
            return
        self._flush_busy_log_buffer()
        self._busy_log_buffer = None
        self._busy_active = False
        self._busy_message = ""
        self.busy_status_label.clear()
        self.busy_status_label.hide()
        self.busy_progress_bar.hide()
        self.btn_apply.setEnabled(True)
        self.btn_revert.setEnabled(True)
        self.left_tabs.setEnabled(True)
        QApplication.restoreOverrideCursor()
        QApplication.processEvents()

    def _set_busy_message(self, message: str) -> None:
        if not self._busy_active:
            return
        self._busy_message = message or self._busy_message
        if self.runtime_mode:
            return
        self.busy_status_label.setText(self._busy_message)
        QApplication.processEvents()

    def _launch_detrend_worker(self, kind: str, fn) -> None:
        worker = _DetrendTaskWorker(self, kind, fn)
        self._detrend_worker = worker
        worker.result_ready.connect(self._on_detrend_task_ready)
        worker.error.connect(self._on_detrend_task_error)
        worker.finished.connect(
            lambda current=worker: self._on_detrend_task_finished(current)
        )
        worker.start()

    def _start_airmass_backfill_if_needed(self) -> None:
        if self.raw_df.empty or "airmass" not in self.raw_df.columns:
            return
        airmass = pd.to_numeric(self.raw_df["airmass"], errors="coerce")
        missing_before = int(
            (airmass.isna() | ~airmass.map(is_reasonable_airmass)).sum()
        )
        if missing_before == 0 or (
            self._detrend_worker is not None and self._detrend_worker.isRunning()
        ):
            return

        def _backfill():
            self._fill_airmass_from_headers(allow_fits_fallback=True)
            updated = pd.to_numeric(self.raw_df["airmass"], errors="coerce")
            missing_after = int(
                (updated.isna() | ~updated.map(is_reasonable_airmass)).sum()
            )
            return {
                "filled": max(0, missing_before - missing_after),
                "remaining": missing_after,
            }

        self.log(f"[LOAD] Background airmass backfill: {missing_before} points")
        self._set_busy_state(True, "Completing airmass metadata...")
        self._launch_detrend_worker("airmass", _backfill)

    def _start_pending_background_fit(self) -> None:
        request = self._pending_fit_request
        if request is None:
            return
        self._pending_fit_request = None

        def _fit():
            self.fit_and_apply(
                update_ui=False,
                save_outputs=False,
                selected_dates=request["selected_dates"],
                use_global_k2=request["use_global_k2"],
                target_id_override=request["target_id"],
                comp_ids_override=request["comp_ids"],
                sync_controls=False,
            )
            return {
                "mode": self.mode,
                "save_outputs": request["save_outputs"],
            }

        self._set_busy_message("Computing correction...")
        self._launch_detrend_worker("fit", _fit)

    def _on_detrend_task_ready(self, kind: str, payload: object) -> None:
        if kind == "airmass":
            result = payload if isinstance(payload, dict) else {}
            self.log(
                f"[LOAD] Background airmass complete: "
                f"filled={result.get('filled', 0)}, remaining={result.get('remaining', 0)}"
            )
            return
        result = payload if isinstance(payload, dict) else {}
        try:
            self._sync_mode_controls_from_state()
            if self.mode == "global":
                self._populate_date_list()
                self._refresh_filter_combo(
                    self.raw_df.get("filter", pd.Series([], dtype=str))
                    .astype(str)
                    .tolist()
                )
            self._update_results_table()
            self._update_plots()
            self._update_analysis_panel()
            if self.mode not in {"global", "sysrem"}:
                self._log_fit_summary()
            if bool(result.get("save_outputs", True)):
                self._set_busy_message("Saving outputs...")
                self._save_comprehensive_results()
        finally:
            self._set_busy_state(False)

    def _on_detrend_task_error(self, kind: str, message: str) -> None:
        if kind == "fit":
            self._set_busy_state(False)
            QMessageBox.warning(self, "Detrend", message)
        else:
            self.log(f"[LOAD] Background airmass failed: {message}")

    def _on_detrend_task_finished(self, worker: _DetrendTaskWorker) -> None:
        if self._detrend_worker is worker:
            self._detrend_worker = None
        worker.deleteLater()
        if self._pending_fit_request is not None:
            self._start_pending_background_fit()
        elif self._busy_active:
            self._set_busy_state(False)

    def _show_mode_help_dialog(self):
        dialog = FittedDialog(self)
        dialog.setWindowTitle("Step 11 보정 모드 도움말")
        dialog.resize(780, 680)

        layout = QVBoxLayout(dialog)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(8)

        browser = QTextBrowser(dialog)
        browser.setReadOnly(True)
        browser.setOpenExternalLinks(False)
        browser.setStyleSheet(
            f"QTextBrowser {{ background-color: {Tokens.SURFACE_ALT}; "
            f"border: 1px solid {Tokens.BORDER}; padding: 8px; font-size: 9pt; }}"
        )
        browser.setHtml(self._build_mode_help_html())
        layout.addWidget(browser, 1)

        btn_row = QHBoxLayout()
        btn_row.addStretch()
        btn_close = QPushButton("닫기")
        btn_close.clicked.connect(dialog.accept)
        btn_row.addWidget(btn_close)
        layout.addLayout(btn_row)

        dialog.exec_()

    def _build_mode_help_html(self) -> str:
        return """
<html>
<body style="font-family:'Malgun Gothic'; font-size:9pt; line-height:1.45;">
<h3 style="margin-top:0;">Step 11 도움말: 보정 모드 설명</h3>
<p>
Step 11은 여러 밤의 관측을 합칠 때 기준선을 맞추는 단계입니다.
모드에 따라 <b>Step 10 차등측광 결과</b>를 그대로 후처리하거나,
<b>Step 7 forced photometry</b>을 다시 읽어 전역 ensemble 해를 구합니다.
</p>

<h4>1. Offset Only (ZP₀)</h4>
<p><b>입력 데이터</b></p>
<ul>
  <li><code>lc_lightcurve/lightcurve_ID*_raw.csv</code></li>
  <li>주요 컬럼: <code>diff_mag_raw</code>, <code>diff_err</code>, <code>airmass</code>, <code>date/night</code></li>
</ul>
<p><b>무엇을 하나</b></p>
<ul>
  <li>각 밤/필터 그룹에서 <code>diff_mag_raw</code>의 가중평균을 <code>ZP₀</code>로 추정합니다.</li>
  <li>기울기 없이 상수항만 제거합니다.</li>
  <li>적용식: <code>Δm_corr = Δm_raw - ZP₀</code></li>
</ul>
<p><b>결과</b></p>
<ul>
  <li>출력은 baseline이 0 근처로 맞춰진 차등광도입니다.</li>
  <li><code>ZP₀</code>는 절대측광 zero point가 아니라, 그 밤 차등광도의 baseline입니다.</li>
  <li>저장 파일: <code>lightcurve_ID*_current.csv</code>, <code>params_ID*_current.csv</code>, <code>summary_ID*_current.txt</code> (모드는 CSV 내부 <code>correction_mode</code>와 <code>result_ID*_current.json</code>에 기록)</li>
</ul>
<p><b>권장 사용처 / 주의</b></p>
<ul>
  <li>빠른 위상 접기, exploratory plotting, baseline centering</li>
  <li>target 자체의 밤평균 변화까지 같이 제거할 수 있으므로, 장기/저주파 변광 해석에는 주의</li>
</ul>

<h4>2. Color-dependent (ZP₀ + k''·ΔC·X)</h4>
<p><b>입력 데이터</b></p>
<ul>
  <li><code>lc_lightcurve/lightcurve_ID*_raw.csv</code></li>
  <li>주요 컬럼: <code>diff_mag_raw</code>, <code>diff_err</code>, <code>airmass</code></li>
  <li>추가 정보: 필터별 <code>ΔC = target color - comparison ensemble color</code></li>
</ul>
<p><b>무엇을 하나</b></p>
<ul>
  <li><code>x = X · ΔC</code>를 만들고 차등광도를 선형 fit합니다.</li>
  <li>적용식: <code>Δm_corr = Δm_raw - ZP₀ - k''·ΔC·X</code></li>
  <li><code>Global k''</code>를 켜면 필터별 <code>k''</code>를 전체 선택 밤 데이터로 먼저 한 번 fit하고, 각 밤에서는 <code>ZP₀</code>만 다시 구합니다.</li>
</ul>
<p><b>결과</b></p>
<ul>
  <li>색차에 의한 잔여 소광 경향까지 제거한 차등광도입니다.</li>
  <li>저장 파일: <code>lightcurve_ID*_current.csv</code>, <code>params_ID*_current.csv</code>, <code>summary_ID*_current.txt</code> (기존 top-level 결과는 <code>_history/</code>로 이동)</li>
</ul>
<p><b>권장 사용처 / 주의</b></p>
<ul>
  <li>target과 comparison ensemble의 색차가 크고, residual이 <code>X·ΔC</code>와 실제로 상관될 때</li>
  <li>여전히 target differential curve를 직접 fit하므로, ΔX가 좁거나 데이터가 적으면 astrophysical trend를 흡수할 수 있음</li>
  <li>check star가 있다면 보정 후 check curve가 더 flat해지는지 같이 확인 권장</li>
</ul>

<h4>3. Global Ensemble (Method C)</h4>
<p><b>입력 데이터</b></p>
<ul>
  <li><code>step7_forced_phot/photometry_index.csv</code>와 각 프레임의 <code>photometry_*.tsv</code></li>
  <li>사용 별: Step 9/10에서 선택한 target + comparison ensemble</li>
  <li>solver는 comparison stars만으로 frame-wise zero-point를 풉니다. target은 기준선 해를 구할 때 제외됩니다.</li>
</ul>
<p><b>무엇을 하나</b></p>
<ul>
  <li>모델: <code>mag_inst(i,t,f) = M_i,f + Z_t,f + ε</code></li>
  <li>comparison ensemble 평균광도 <code>M_i</code>와 프레임별 zero-point <code>Z_t</code>를 전역 최소제곱으로 동시에 추정합니다.</li>
  <li>outlier rejection, RMS가 큰 comparison 제거, frame sigma clipping을 같이 수행합니다.</li>
</ul>
<p><b>결과</b></p>
<ul>
  <li><code>diff_mag_raw</code>: 기존 방식의 <code>m_target - mean(comp)</code></li>
  <li><code>mag_ensemble_corr</code>: <code>m_target - Z_t</code> (ensemble-corrected instrumental magnitude)</li>
  <li><code>diff_mag_corr</code>: <code>(m_target - Z_t) - &lt;M_comp&gt;</code> 로 다시 차등광도 스케일에 맞춘 결과</li>
  <li>추가 산출물: <code>global_zp_ID*_current.csv</code>, <code>global_mean_ID*_current.csv</code>, <code>global_diagnostics_ID*_current.json</code></li>
  <li>저장 파일: <code>lightcurve_ID*_current.csv</code>, <code>params_ID*_current.csv</code>, <code>summary_ID*_current.txt</code></li>
</ul>
<p><b>권장 사용처 / 주의</b></p>
<ul>
  <li>멀티나잇 기준선을 comparison ensemble로 정의하고 싶을 때</li>
  <li>target 신호를 기준선 정의에 직접 넣고 싶지 않을 때</li>
  <li>현재 <code>Z_t</code>는 nightly가 아니라 <b>frame-wise</b>입니다.</li>
  <li>comparison에 변광성이 섞이면 결과가 바로 무너집니다.</li>
</ul>

<h4>모드 선택 요약</h4>
<ul>
  <li><b>Offset</b>: 가장 단순한 중심 맞춤. 빠른 검사/위상 확인용.</li>
  <li><b>Color</b>: 색차 소광이 실제로 보일 때만 선택.</li>
  <li><b>Global</b>: comparison ensemble 기준의 문헌적인 접근이 필요할 때 가장 권장.</li>
</ul>
</body>
</html>
"""

    def _auto_load_ids(self):
        """Step 10 comp_selection.json → Step 9 merged selection 순서로 자동 로드."""
        rd = Path(self.params.P.result_dir)
        # 1) Step 10 output (comp_selection.json)
        sel_path = step9_lc_dir(rd) / "comp_selection.json"
        if sel_path.exists():
            try:
                data = json.loads(sel_path.read_text(encoding="utf-8"))
                target_id = data.get("target_id")
                if target_id is not None:
                    self.target_edit.setText(str(target_id))
                self.comp_active_ids = [int(x) for x in data.get("comp_active_ids", []) if str(x).strip()]
                self.comp_candidate_ids = [int(x) for x in data.get("comp_candidate_ids", []) if str(x).strip()]
                self._update_id_info_label()
                return
            except Exception:
                pass
        # 2) Step 9 merged selection
        tid, cids = _load_step9_selection_ids(rd)
        if tid is not None:
            self.target_edit.setText(str(int(tid)))
        if cids:
            self.comp_active_ids = list(cids)
            self.comp_candidate_ids = list(cids)
        if tid is not None or cids:
            self._update_id_info_label()
            return
        # 3) Fallback: parse target ID from raw filenames
        raw_paths = list(step9_lc_dir(rd).glob("lightcurve_ID*_raw.csv"))
        if raw_paths:
            target_id = self._parse_target_id_from_name(raw_paths[0].name)
            if target_id is not None:
                self.target_edit.setText(str(target_id))
        self._update_id_info_label()

    def _update_id_info_label(self):
        """Read-only info label 갱신."""
        t = self.target_edit.text().strip()
        comp_ids = self.comp_active_ids or self.comp_candidate_ids
        c = ", ".join(str(i) for i in comp_ids[:8]) if comp_ids else ""
        if len(comp_ids) > 8:
            c += f" +{len(comp_ids) - 8}"
        if t and c:
            self.id_info_label.setText(f"Target: ID {t}  |  Comp: {c}")
            banner = "ok"
        elif t:
            self.id_info_label.setText(f"Target: ID {t}  |  Comp: (없음)")
            banner = "warn"
        else:
            self.id_info_label.setText("Target / Comp: Step 9를 먼저 실행해주세요")
            banner = "error"
        self.id_info_label.setProperty("banner", banner)
        refresh(self.id_info_label)

    def _sync_mode_controls_from_state(self) -> None:
        is_global = self.mode in ("global", "sysrem")
        self.color_map_group.setEnabled(not is_global)
        self.chk_global_k2.setEnabled(not is_global)
        blockers = [
            QSignalBlocker(self.mode_offset),
            QSignalBlocker(self.mode_color),
            QSignalBlocker(self.mode_global),
        ]
        if hasattr(self, "mode_sysrem"):
            blockers.append(QSignalBlocker(self.mode_sysrem))
        try:
            if self.mode == "color":
                self.mode_color.setChecked(True)
            elif self.mode == "global":
                self.mode_global.setChecked(True)
            elif self.mode == "sysrem" and hasattr(self, "mode_sysrem"):
                self.mode_sysrem.setChecked(True)
            else:
                self.mode = "offset"
                self.mode_offset.setChecked(True)
        finally:
            del blockers

    def _available_color_expressions(self, result_dir: Path, bands_hint: list[str] | None = None) -> list[str]:
        """Return color-index expressions built from filters present in the median table.

        Uses :func:`filter_bands_from_columns` so case (Johnson uppercase vs SDSS
        lowercase) is preserved and derived columns like ``mag_inst_err_B`` are skipped.
        """
        from apex.utils.gaia_transforms import filter_bands_from_columns, build_color_pairs
        df = self._load_color_median_table(result_dir)
        bands: list[str] = []
        if not df.empty:
            bands = (
                filter_bands_from_columns(df.columns, "mag_cal_")
                or filter_bands_from_columns(df.columns, "mag_std_")
                or filter_bands_from_columns(df.columns, "mag_inst_")
            )
        if not bands and bands_hint:
            bands = [_normalize_filter_key(b) for b in bands_hint if str(b).strip()]
            bands = sorted({b for b in bands if b})
        # All adjacent pairs in both directions (blue-red and red-blue) for UI choice
        pairs = build_color_pairs(bands, adjacent_only=True)
        exprs = []
        for a, b in pairs:
            exprs.append(f"{a}-{b}")
            exprs.append(f"{b}-{a}")
        return exprs

    def _normalize_color_map(self, mapping) -> dict[str, str]:
        if not isinstance(mapping, dict):
            return {}
        out: dict[str, str] = {}
        for key, value in mapping.items():
            fkey = _normalize_filter_key(key)
            expr = str(value).strip()
            if fkey and expr:
                out[fkey] = expr
        return out

    def _on_color_map_changed(self) -> None:
        if not hasattr(self, "color_map_combos"):
            return
        mapping: dict[str, str] = {}
        for fkey, combo in self.color_map_combos.items():
            expr = combo.currentText().strip()
            if expr and expr.lower() != "none":
                mapping[fkey] = expr
        self.color_map_by_filter = mapping
        self._refresh_delta_c_map()
        self._log_color_index_info()
        self._update_color_mode_enabled()
        self._update_analysis_panel()

    def _on_date_selection_changed(self, _item: QListWidgetItem):
        self._update_plots()

    def _on_filter_changed(self, _idx: int) -> None:
        self._update_plots()

    def _on_color_by_changed(self, _idx: int) -> None:
        self.color_by = self.color_by_combo.currentText() if hasattr(self, "color_by_combo") else "Date"
        self._update_plots()

    def _set_mode(self, mode: str):
        if self.mode == mode:
            return
        self.mode = mode
        is_global = mode in ("global", "sysrem")
        self.color_map_group.setEnabled(not is_global)
        self.chk_global_k2.setEnabled(not is_global)
        self.mode_color.setEnabled(True)
        if not self.raw_df.empty:
            self.corrected_df = pd.DataFrame()
            self.params_df = pd.DataFrame()
            self.global_mean_df = pd.DataFrame()
            self.global_diagnostics = {}
            self._update_results_table()
            self._update_plots()
            self.log(f"[MODE] Switched to {mode}. Run Fit && Apply to update.")

    def _on_clip_changed(self):
        self.sigma_clip = bool(self.chk_clip.isChecked())
        self.clip_sigma = float(self.spin_clip.value())
        self.clip_iters = int(self.spin_iters.value())

    def _on_xaxis_changed(self, idx: int) -> None:
        self.x_axis_mode = "phase" if idx == 1 else "time"
        self._update_plots()

    def _on_phase_params_changed(self):
        self.phase_period = float(self.spin_period.value())
        self.phase_t0 = float(self.spin_t0.value())
        if hasattr(self, "spin_cycles"):
            self.phase_cycles = float(self.spin_cycles.value())
        if self.x_axis_mode == "phase":
            self._update_plots()

    def _auto_set_t0(self):
        """Set T0 to min(JD) from loaded data."""
        if self.raw_df.empty:
            QMessageBox.information(self, "T₀ Auto", "먼저 데이터를 로드하세요.")
            return
        if "JD" not in self.raw_df.columns:
            QMessageBox.information(self, "T₀ Auto", "JD 컬럼이 없습니다.")
            return
        jd = pd.to_numeric(self.raw_df["JD"], errors="coerce").to_numpy(float)
        jd = jd[np.isfinite(jd)]
        if jd.size == 0:
            return
        t0_auto = float(np.min(jd))
        self.spin_t0.setValue(t0_auto)
        self.log(f"[PHASE] T₀ auto-set to min(JD) = {t0_auto:.6f}")

    def revert_raw(self):
        self.corrected_df = pd.DataFrame()
        self.params_df = pd.DataFrame()
        self.global_mean_df = pd.DataFrame()
        self.global_diagnostics = {}
        self._update_results_table()
        self._update_plots()

    def validate_step(self) -> bool:
        target_id = None
        try:
            # AttributeError: base __init__ validates once before
            # setup_step_ui() has created target_edit.
            text = self.target_edit.text().strip()
            if text:
                target_id = int(text)
        except (TypeError, ValueError, AttributeError):
            pass
        detrend_dir = step10_detrend_dir(self.params.P.result_dir).resolve()
        for path in list_lightcurve_csvs(self.params.P.result_dir, target_id):
            try:
                if path.resolve().parent != detrend_dir:
                    continue
                if not pd.read_csv(path, nrows=1).empty:
                    return True
            except Exception:
                continue
        return False

    def save_state(self):
        state_data = {
            "mode": self.mode,
            "clip_sigma": self.clip_sigma,
            "clip_iters": self.clip_iters,
            "sigma_clip": self.sigma_clip,
            "plot_view_mode": self._plot_view_mode,
            "x_axis_mode": self.x_axis_mode,
            "phase_period": self.phase_period,
            "phase_t0": self.phase_t0,
            "phase_cycles": self.phase_cycles,
            "color_map_by_filter": self.color_map_by_filter,
            "color_by": self.color_by,
            "use_global_k2": self.chk_global_k2.isChecked(),
            "global_min_comps": self.global_min_comps,
            "global_sigma": self.global_sigma,
            "global_iters": self.global_iters,
            "global_rms_pct": self.global_rms_pct,
            "global_rms_threshold": self.global_rms_threshold,
            "global_frame_sigma": self.global_frame_sigma,
            "global_gauge": self.global_gauge,
            "global_robust": self.global_robust,
            "global_interp_missing": self.global_interp_missing,
            "global_normalize": self.global_normalize,
            "global_rescale_errors": self.global_rescale_errors,
        }
        if hasattr(self, "spin_sysrem_iter"):
            state_data["sysrem_iter"] = int(self.spin_sysrem_iter.value())
        if hasattr(self, "spin_sysrem_apply"):
            state_data["sysrem_apply"] = int(self.spin_sysrem_apply.value())
        self.project_state.store_step_data("detrend_merge", state_data)

    def restore_state(self):
        state_data = self.project_state.get_step_data("detrend_merge")
        if state_data:
            self.mode = state_data.get("mode", "offset")
            self.clip_sigma = float(state_data.get("clip_sigma", 3.0))
            self.clip_iters = int(state_data.get("clip_iters", 2))
            self.sigma_clip = bool(state_data.get("sigma_clip", True))
            view = state_data.get("plot_view_mode", self._plot_view_mode)
            if view in ("raw", "corr", "diag", "all"):
                self._plot_view_mode = view
                if hasattr(self, "_plot_view_buttons"):
                    self._plot_view_buttons[view].setChecked(True)
                    self._apply_plot_view(redraw=False)
            self.x_axis_mode = state_data.get("x_axis_mode", self.x_axis_mode)
            self.phase_period = float(state_data.get("phase_period", 0.0))
            self.phase_t0 = float(state_data.get("phase_t0", 0.0))
            self.phase_cycles = float(state_data.get("phase_cycles", self.phase_cycles))
            self.color_map_by_filter = self._normalize_color_map(state_data.get("color_map_by_filter", {}))
            self.color_by = state_data.get("color_by", self.color_by)
            use_global_k2 = state_data.get("use_global_k2", True)
            self.chk_global_k2.setChecked(bool(use_global_k2))
            self.global_min_comps = int(state_data.get("global_min_comps", self.global_min_comps))
            self.global_sigma = float(state_data.get("global_sigma", self.global_sigma))
            self.global_iters = int(state_data.get("global_iters", self.global_iters))
            self.global_rms_pct = float(state_data.get("global_rms_pct", self.global_rms_pct))
            self.global_rms_threshold = float(state_data.get("global_rms_threshold", self.global_rms_threshold))
            self.global_frame_sigma = float(state_data.get("global_frame_sigma", self.global_frame_sigma))
            self.global_gauge = str(state_data.get("global_gauge", self.global_gauge))
            self.global_robust = bool(state_data.get("global_robust", self.global_robust))
            self.global_interp_missing = bool(state_data.get("global_interp_missing", self.global_interp_missing))
            self.global_normalize = bool(state_data.get("global_normalize", self.global_normalize))
            self.global_rescale_errors = bool(state_data.get("global_rescale_errors", self.global_rescale_errors))
            if hasattr(self, "spin_sysrem_iter"):
                self.spin_sysrem_iter.setValue(int(state_data.get("sysrem_iter", self.spin_sysrem_iter.value())))
            if hasattr(self, "spin_sysrem_apply"):
                self.spin_sysrem_apply.setValue(int(state_data.get("sysrem_apply", self.spin_sysrem_apply.value())))

        if not self.color_map_by_filter:
            self.color_map_by_filter = self._normalize_color_map(
                getattr(self.params.P, "lightcurve_color_index_by_filter", {})
            )

        # Sync UI
        self._sync_mode_controls_from_state()

        self.chk_clip.setChecked(self.sigma_clip)
        self.spin_clip.setValue(self.clip_sigma)
        self.spin_iters.setValue(self.clip_iters)
        self.spin_period.setValue(self.phase_period)
        self.spin_t0.setValue(self.phase_t0)
        if hasattr(self, "spin_cycles"):
            self.spin_cycles.setValue(self.phase_cycles)
        if hasattr(self, "phase_mode_combo"):
            self.phase_mode_combo.setCurrentIndex(1 if self.x_axis_mode == "phase" else 0)
        if hasattr(self, "color_by_combo"):
            self.color_by_combo.setCurrentText(self.color_by)
        if hasattr(self, "spin_global_min_comps"):
            self.spin_global_min_comps.setValue(self.global_min_comps)
        if hasattr(self, "spin_global_sigma"):
            self.spin_global_sigma.setValue(self.global_sigma)
        if hasattr(self, "spin_global_iters"):
            self.spin_global_iters.setValue(self.global_iters)
        if hasattr(self, "spin_global_rms_pct"):
            self.spin_global_rms_pct.setValue(self.global_rms_pct)
        if hasattr(self, "spin_global_rms_thr"):
            self.spin_global_rms_thr.setValue(self.global_rms_threshold)
        if hasattr(self, "spin_global_frame_sigma"):
            self.spin_global_frame_sigma.setValue(self.global_frame_sigma)
        if hasattr(self, "combo_global_gauge"):
            self.combo_global_gauge.setCurrentText(self.global_gauge)
        if hasattr(self, "chk_global_robust"):
            self.chk_global_robust.setChecked(self.global_robust)
        if hasattr(self, "chk_global_interp"):
            self.chk_global_interp.setChecked(self.global_interp_missing)
        if hasattr(self, "chk_global_normalize"):
            self.chk_global_normalize.setChecked(self.global_normalize)
        if hasattr(self, "chk_global_rescale_err"):
            self.chk_global_rescale_err.setChecked(self.global_rescale_errors)
