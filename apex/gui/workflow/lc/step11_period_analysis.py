"""
Step 12: Period Analysis (Lomb-Scargle / PDM / BLS)

Quick period scan.  Detailed analysis (refine, bootstrap, T₀, O-C,
variable-star / transit / EB fitting) lives in the Tools menu.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional, List

import numpy as np
import pandas as pd

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from PyQt5.QtWidgets import (
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QLabel,
    QGroupBox,
    QCheckBox,
    QMessageBox,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QTextEdit,
    QWidget,
    QComboBox,
    QFormLayout,
    QDoubleSpinBox,
    QSpinBox,
    QTabWidget,
    QScrollArea,
    QFrame,
)
from PyQt5.QtCore import Qt, QThread, QTimer, pyqtSignal

from apex.gui.layout_rules import tame_canvas
from apex.gui.theme import style_button, Tokens
from apex.gui.workflow.step_window_base import StepWindowBase
from apex.analysis.light_curve.target_config import read_step8_selection
from apex.analysis.light_curve.period_analysis_service import (
    compute_ls,
    run_period_analysis,
)
from apex.analysis.light_curve.period_alias_service import (
    analyze_period_aliases,
    diagnose_multimode_suspicion,
    periods_are_window_aliases,
)
from apex.analysis.light_curve.period_io_service import (
    ALL_FILTER_KEY,
    corr_mode_label,
    detect_corr_mode_from_df,
    load_period_lightcurve_csv,
    save_period_analysis_outputs,
)
from apex.analysis.light_curve.period_plot import (
    PeriodSummaryPlotter,
    _FILTER_PLOT_COLORS,
    _load_check_star_for_plot,
)
from apex.analysis.light_curve.variable_analysis_contract import (
    ValidatedLightCurveBundle,
    compute_file_fingerprint,
)
from apex.utils.common_helpers import normalize_filter_key as _normalize_filter_key
from apex.utils.io_utils import read_csv_int64_source_id, coerce_int64_source_id
from apex.utils.photometry_provenance import (
    format_photometry_provenance,
    summarize_photometry_table,
)
from apex.utils.step_paths_lc import (
    step8_selection_dir,
    step9_lc_dir,
    step11_period_dir,
    find_best_lightcurve_csv,
    load_detrend_preference,
)


_ALL_FILTER_LABEL = "All filters (median aligned)"
def _load_step9_source_to_id_map(result_dir: Path, flt: str | None = None) -> dict[int, int]:
    step9_out = step8_selection_dir(result_dir)
    if not step9_out.exists():
        return {}

    key = _normalize_filter_key(flt)
    candidates: list[tuple[Path, str]] = []
    if key:
        candidates.extend(
            [
                (step9_out / f"master_catalog_{key}.tsv", "\t"),
                (step9_out / f"id_mapping_{key}.csv", ","),
            ]
        )
    candidates.extend((p, "\t") for p in sorted(step9_out.glob("master_catalog_*.tsv")))
    candidates.extend((p, ",") for p in sorted(step9_out.glob("id_mapping_*.csv")))

    mapping: dict[int, int] = {}
    for path, sep in candidates:
        if not path.exists():
            continue
        try:
            df = read_csv_int64_source_id(path, sep=sep)
        except Exception:
            continue
        if not {"source_id", "ID"} <= set(df.columns):
            continue
        sid_vals = coerce_int64_source_id(df["source_id"])
        id_vals = pd.to_numeric(df["ID"], errors="coerce").astype("Int64")
        for sid_val, id_val in zip(sid_vals, id_vals):
            if pd.isna(sid_val) or pd.isna(id_val):
                continue
            sid_int = int(sid_val)
            if sid_int not in mapping:
                mapping[sid_int] = int(id_val)
    return mapping


def _load_step9_target_id(result_dir: Path) -> int | None:
    step9_out = step8_selection_dir(result_dir)
    if not step9_out.exists():
        return None

    target_ids: set[int] = set()
    for sp in sorted(step9_out.glob("selection_*.json")):
        try:
            data = json.loads(sp.read_text(encoding="utf-8"))
        except Exception:
            continue
        tid = data.get("target_id")
        if tid is None:
            flt = sp.stem.replace("selection_", "")
            sid_map = _load_step9_source_to_id_map(result_dir, flt)
            target_sid = data.get("target_source_id")
            if target_sid is not None and int(target_sid) in sid_map:
                tid = int(sid_map[int(target_sid)])
        if tid is not None:
            target_ids.add(int(tid))

    if len(target_ids) == 1:
        return next(iter(target_ids))
    return None


def _detect_corr_mode(filename: str) -> tuple[str, str]:
    return detect_corr_mode_from_df(pd.DataFrame(), filename)


class PeriodAnalysisWorker(QThread):
    """Worker thread for period analysis (Lomb-Scargle, PDM, BLS)."""
    progress = pyqtSignal(str)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)

    def __init__(
        self,
        time: np.ndarray,
        mag_raw: np.ndarray,
        mag_corr: np.ndarray,
        mag_err: Optional[np.ndarray],
        min_period: float,
        max_period: float,
        samples_per_peak: int = 10,
        methods: Optional[List[str]] = None,
        pdm_n_bins: int = 10,
        night_id: Optional[np.ndarray] = None,
        correction_preserves_nightly_baseline: bool = True,
        include_alias_diagnostics: bool = False,
        include_multimode_diagnostic: bool = False,
    ):
        super().__init__()
        self.time = time
        self.mag_raw = mag_raw
        self.mag_corr = mag_corr
        self.mag_err = mag_err
        self.min_period = min_period
        self.max_period = max_period
        self.samples_per_peak = samples_per_peak
        self.methods = methods or ["ls"]
        self.pdm_n_bins = pdm_n_bins
        self.night_id = night_id
        self.correction_preserves_nightly_baseline = bool(
            correction_preserves_nightly_baseline
        )
        self.include_alias_diagnostics = bool(include_alias_diagnostics)
        self.include_multimode_diagnostic = bool(include_multimode_diagnostic)

    def run(self):
        try:
            results = run_period_analysis(
                time=self.time,
                mag_raw=self.mag_raw,
                mag_corr=self.mag_corr,
                mag_err=self.mag_err,
                min_period=self.min_period,
                max_period=self.max_period,
                samples_per_peak=self.samples_per_peak,
                methods=self.methods,
                pdm_n_bins=self.pdm_n_bins,
                progress_cb=self.progress.emit,
            )
            use_corrected = self.correction_preserves_nightly_baseline and (
                self.mag_corr is not None and np.any(np.isfinite(self.mag_corr))
            )
            mag = self.mag_corr if use_corrected else self.mag_raw
            needs_diagnostics = self.include_alias_diagnostics or self.include_multimode_diagnostic
            ls_result = (
                results.get("corr_ls") if use_corrected else results.get("raw_ls")
            ) or results.get("raw_ls")
            if needs_diagnostics and not ls_result:
                ls_result = compute_ls(
                    self.time,
                    mag,
                    self.mag_err,
                    "alias-scan",
                    self.min_period,
                    self.max_period,
                    self.samples_per_peak,
                )
            if needs_diagnostics and ls_result and "error" not in ls_result and "frequency" in ls_result:
                self.progress.emit("Evaluating sampling-window aliases...")
                alias_analysis = analyze_period_aliases(
                    self.time,
                    mag,
                    self.mag_err,
                    self.night_id,
                    ls_result["frequency"],
                    ls_result["power"],
                    self.min_period,
                    self.max_period,
                    harmonics=2,
                )
                alias_analysis["input_series"] = "corrected" if use_corrected else "raw"
                alias_analysis["nightly_baseline_preserved"] = bool(
                    self.correction_preserves_nightly_baseline
                )
                if self.include_alias_diagnostics:
                    results["alias_analysis"] = alias_analysis
                if self.include_multimode_diagnostic:
                    self.progress.emit("Testing for independent residual frequencies...")
                    results["multimode_diagnostic"] = diagnose_multimode_suspicion(
                        self.time,
                        mag,
                        self.mag_err,
                        self.night_id,
                        alias_analysis,
                        self.min_period,
                        self.max_period,
                        harmonics=2,
                        samples_per_peak=max(self.samples_per_peak, 10),
                    )
            self.finished.emit(results)
        except Exception as e:
            import traceback
            self.error.emit(f"{e}\n{traceback.format_exc()}")


class PeriodAnalysisWindow(StepWindowBase, PeriodSummaryPlotter):
    """Step 12: Period Analysis — quick scan with Lomb-Scargle / PDM / BLS."""

    # -- the seam: what the batch path is handed, the window reads ----------

    def _search_window(self) -> tuple[float, float, int]:
        return (float(self.min_period_spin.value()),
                float(self.max_period_spin.value()),
                int(self.samples_spin.value()))

    def _show_alias_marks(self) -> bool:
        return bool(self.chk_alias.isChecked())

    def _load_check_star(self, requested_filter):
        return _load_check_star_for_plot(
            Path(self.params.P.result_dir), requested_filter)

    def _summary_figure(self):
        return self.summary_canvas.figure

    def _summary_canvas_width(self) -> int:
        return self.summary_canvas.width()

    def _summary_redraw(self) -> None:
        self.summary_canvas.setMinimumHeight(260)
        self.summary_canvas.draw_idle()

    def __init__(self, params, file_manager, project_state, main_window):
        self.file_manager = file_manager
        self.worker = None
        self.results = {}
        self.multinight = None
        self.alias_analysis = None
        self.multimode_diagnostic = None
        self.variable_analysis_bundle = None
        self.lc_data = None
        self.current_filter = None
        self._ui_ready = False
        self._summary_layout_compact = None
        self._summary_layout_stacked = None
        self._summary_resize_pending = False
        self._periodogram_layout_compact = None
        self._periodogram_resize_pending = False
        self._shell_compact = None
        self._responsive_shell_pending = False
        self._check_star_plot_cache_key = None
        self._check_star_plot_cache = None
        self._check_star_ls_cache = {}

        super().__init__(
            step_index=11,
            step_name="Period Analysis",
            params=params,
            project_state=project_state,
            main_window=main_window,
        )

        self.setup_step_ui()
        self.restore_state()
        self._auto_load_target_id()
        self._ui_ready = True
        self._load_light_curve(silent=True)

    def _auto_load_target_id(self):
        """Step 10 light curve output -> Step 9 selection 순서로 target ID 자동 로드."""
        rd = Path(self.params.P.result_dir)
        sel_path = step9_lc_dir(rd) / "comp_selection.json"
        if sel_path.exists():
            try:
                data = json.loads(sel_path.read_text(encoding="utf-8"))
                target_id = data.get("target_id")
                if target_id is not None:
                    self.target_id_spin.setValue(int(target_id))
                    self.target_hint.setText("(Step 10)")
                    return
            except Exception:
                pass
        tid = _load_step9_target_id(rd)
        if tid is not None:
            self.target_id_spin.setValue(int(tid))
            self.target_hint.setText("(Step 9)")
            return
        self.target_hint.setText("")

    def setup_step_ui(self):
        # 2-column "tool" layout: a fixed-width, scrollable control column on
        # the left; the plot/result tabs take every remaining pixel on the
        # right. The old vertical stack needed 1309px minimum height (clipped
        # even on a 1528px monitor, and Qt squeezed the plot canvas first).
        # The intro banner is gone — the same text lives in the ⓘ 가이드 popup.
        t = Tokens
        self.content_layout.setContentsMargins(t.MARGIN, t.MARGIN, t.MARGIN, t.MARGIN)
        self.content_layout.setSpacing(t.S3)

        root = QHBoxLayout()
        root.setSpacing(t.S3)
        self.content_layout.addLayout(root, 1)

        self.btn_controls = self.add_header_action(
            "Controls",
            self._toggle_control_column,
            tooltip="Show or hide data and period-search controls.",
            min_width=86,
        )
        self.btn_controls.setCheckable(True)

        # ── Left: control column (scrolls; never drives window height) ──
        left_host = QWidget()
        lv = QVBoxLayout(left_host)
        lv.setContentsMargins(0, 0, t.S2, 0)  # right pad so the scrollbar doesn't kiss the fields
        lv.setSpacing(t.S3)

        # Prominent result callout — the trustworthy period, surfaced out of the
        # log so the user's eye lands on the answer, not the scroll-back.
        self.result_callout = QLabel("")
        self.result_callout.setWordWrap(True)
        self.result_callout.setVisible(False)
        lv.addWidget(self.result_callout)

        self.mode_callout = QLabel("")
        self.mode_callout.setWordWrap(True)
        self.mode_callout.setVisible(False)
        lv.addWidget(self.mode_callout)

        # Data selection
        data_group = QGroupBox("Data Selection")
        data_layout = QFormLayout(data_group)

        target_row = QHBoxLayout()
        self.target_id_spin = QSpinBox()
        self.target_id_spin.setRange(1, 99999)
        self.target_id_spin.setValue(1)
        self.target_id_spin.valueChanged.connect(self._on_target_id_changed)
        target_row.addWidget(self.target_id_spin)
        self.target_hint = QLabel("")
        self.target_hint.setStyleSheet(f"QLabel {{ color: {Tokens.OK}; }}")
        target_row.addWidget(self.target_hint)
        target_row.addStretch()
        data_layout.addRow("Target ID:", target_row)

        self.source_label = QLabel("—")
        self.source_label.setStyleSheet(f"QLabel {{ color: {Tokens.TEXT_SUB}; }}")
        data_layout.addRow("Data source:", self.source_label)

        self.photometry_source_label = QLabel(
            "Photometry: Unknown | MAG: unknown"
        )
        self.photometry_source_label.setStyleSheet(
            f"QLabel {{ color: {Tokens.TEXT_SUB}; font-weight: bold; }}"
        )
        data_layout.addRow(self.photometry_source_label)

        self.filter_combo = QComboBox()
        self.filter_combo.currentIndexChanged.connect(self._on_filter_changed)
        data_layout.addRow("Filter:", self.filter_combo)

        self.data_status = QLabel("Loading light curve data...")
        self.data_status.setWordWrap(True)
        data_layout.addRow("Status:", self.data_status)

        lv.addWidget(data_group)

        # Period search parameters
        param_group = QGroupBox("Period Search Parameters")
        param_layout = QFormLayout(param_group)

        self.min_period_spin = QDoubleSpinBox()
        self.min_period_spin.setRange(0.001, 100.0)
        self.min_period_spin.setDecimals(4)
        self.min_period_spin.setValue(0.01)
        self.min_period_spin.setSuffix(" days")
        param_layout.addRow("Min Period:", self.min_period_spin)

        self.max_period_spin = QDoubleSpinBox()
        self.max_period_spin.setRange(0.01, 1000.0)
        self.max_period_spin.setDecimals(4)
        self.max_period_spin.setValue(10.0)
        self.max_period_spin.setSuffix(" days")
        param_layout.addRow("Max Period:", self.max_period_spin)

        self.samples_spin = QSpinBox()
        self.samples_spin.setRange(5, 100)
        self.samples_spin.setValue(10)
        param_layout.addRow("Samples per peak:", self.samples_spin)

        # Spanning row (no label column): with the form's label column the
        # three checkboxes need ~430px and clip the 400px control column.
        method_row = QHBoxLayout()
        method_row.addWidget(QLabel("Methods:"))
        self.chk_ls = QCheckBox("Lomb-Scargle")
        self.chk_ls.setChecked(True)
        method_row.addWidget(self.chk_ls)
        self.chk_pdm = QCheckBox("PDM")
        self.chk_pdm.setChecked(True)
        method_row.addWidget(self.chk_pdm)
        self.chk_bls = QCheckBox("BLS")
        self.chk_bls.setChecked(False)
        method_row.addWidget(self.chk_bls)
        method_row.addStretch()
        param_layout.addRow(method_row)

        self.pdm_bins_spin = QSpinBox()
        self.pdm_bins_spin.setRange(5, 50)
        self.pdm_bins_spin.setValue(10)
        param_layout.addRow("PDM bins:", self.pdm_bins_spin)

        lv.addWidget(param_group)

        # Run — single primary action, full column width.
        self.btn_run = QPushButton("Compute Periodogram")
        style_button(self.btn_run, "primary", height=Tokens.H_ACTION)
        self.btn_run.clicked.connect(self._run_analysis)
        self.btn_run.setEnabled(False)
        lv.addWidget(self.btn_run)

        self.progress_label = QLabel("")
        self.progress_label.setStyleSheet(f"QLabel {{ color: {Tokens.TEXT_SUB}; }}")
        lv.addWidget(self.progress_label)

        self.btn_open_varstar = QPushButton("Open in Variable Star Tool")
        style_button(self.btn_open_varstar, height=Tokens.H_BUTTON)
        self.btn_open_varstar.setToolTip(
            "Open the current light curve, adopted period, alias candidates, "
            "and mode diagnostic in the Variable Star Analysis Tool."
        )
        self.btn_open_varstar.setEnabled(False)
        self.btn_open_varstar.clicked.connect(self._open_variable_star_tool)
        lv.addWidget(self.btn_open_varstar)

        lv.addStretch(1)

        ctrl_scroll = QScrollArea()
        ctrl_scroll.setWidgetResizable(True)
        ctrl_scroll.setFrameShape(QFrame.NoFrame)
        ctrl_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        ctrl_scroll.setWidget(left_host)

        left_col = QWidget()
        # Measured on the real (windows) font stack: the control stack's
        # minimumSizeHint is 453px wide — the Methods checkbox row dominates.
        # 480 = 453 + vertical-scrollbar allowance; narrower clips the right
        # edge because the horizontal scrollbar is off.
        left_col.setFixedWidth(480)
        left_v = QVBoxLayout(left_col)
        left_v.setContentsMargins(0, 0, 0, 0)
        left_v.setSpacing(t.S2)
        left_v.addWidget(ctrl_scroll, 1)
        # Log stays pinned below the scroll so it's always visible.
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMaximumHeight(100)
        self.log_text.setStyleSheet(
            f"QTextEdit {{ background: {Tokens.SURFACE_ALT}; color: {Tokens.TEXT_SUB}; }}"
        )
        left_v.addWidget(self.log_text)
        self.control_column = left_col
        root.addWidget(left_col)

        # Results tabs
        self.tabs = QTabWidget()

        # Summary tab: the complete quick-look result on one canvas. Detailed
        # method and interactive phase controls remain available in later tabs.
        summary_tab = QWidget()
        summary_layout = QVBoxLayout(summary_tab)
        self.summary_canvas = FigureCanvas(
            Figure(figsize=(11, 8), constrained_layout=True)
        )
        self.summary_canvas.mpl_connect(
            "resize_event", self._on_summary_canvas_resize
        )
        summary_layout.addWidget(tame_canvas(self.summary_canvas), 1)
        self.tabs.addTab(summary_tab, "Summary")

        # Periodogram tab
        periodogram_tab = QWidget()
        periodogram_layout = QVBoxLayout(periodogram_tab)

        alias_row = QHBoxLayout()
        self.chk_alias = QCheckBox("Show sampling-window candidates")
        self.chk_alias.setChecked(True)
        self.chk_alias.toggled.connect(self._update_periodogram_plot)
        self.chk_alias.toggled.connect(self._update_summary_plot)
        alias_row.addWidget(self.chk_alias)
        alias_row.addStretch()
        periodogram_layout.addLayout(alias_row)

        self.periodogram_canvas = FigureCanvas(Figure(figsize=(10, 5)))
        self.periodogram_canvas.mpl_connect(
            "resize_event", self._on_periodogram_canvas_resize
        )
        tame_canvas(self.periodogram_canvas, min_h=400)
        self.periodogram_scroll = QScrollArea()
        self.periodogram_scroll.setWidgetResizable(True)
        self.periodogram_scroll.setFrameShape(QFrame.NoFrame)
        self.periodogram_scroll.setWidget(self.periodogram_canvas)
        periodogram_layout.addWidget(self.periodogram_scroll, 1)

        self.tabs.addTab(periodogram_tab, "Periodogram")

        # Phase plot tab
        phase_tab = QWidget()
        phase_layout = QVBoxLayout(phase_tab)

        phase_control = QHBoxLayout()
        phase_control.addWidget(QLabel("Period for phase plot:"))

        self.phase_period_combo = QComboBox()
        self.phase_period_combo.currentIndexChanged.connect(self._update_phase_plot)
        phase_control.addWidget(self.phase_period_combo)

        self.phase_period_edit = QDoubleSpinBox()
        self.phase_period_edit.setRange(0.0001, 10000.0)
        self.phase_period_edit.setDecimals(6)
        self.phase_period_edit.setSuffix(" days")
        self.phase_period_edit.valueChanged.connect(self._update_phase_plot_custom)
        phase_control.addWidget(self.phase_period_edit)

        phase_control.addStretch()
        phase_layout.addLayout(phase_control)

        self.phase_canvas = FigureCanvas(Figure(figsize=(10, 6)))
        phase_layout.addWidget(tame_canvas(self.phase_canvas), 1)

        self.tabs.addTab(phase_tab, "Phase Plot")

        # Results table tab
        results_tab = QWidget()
        results_layout = QVBoxLayout(results_tab)

        self.results_table = QTableWidget()
        self.results_table.setColumnCount(7)
        self.results_table.setHorizontalHeaderLabels([
            "Method", "Data", "Best Period (days)", "Power", "FAP", "Alias?", "Top 3 Periods"
        ])
        self.results_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.results_table.horizontalHeader().setStretchLastSection(True)
        results_layout.addWidget(self.results_table)

        self.alias_table = QTableWidget()
        self.alias_table.setColumnCount(7)
        self.alias_table.setHorizontalHeaderLabels([
            "Rank", "Period (days)", "Freq (d^-1)", "Rel. power",
            "Delta BIC", "LOO votes", "Relation",
        ])
        self.alias_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.alias_table.horizontalHeader().setStretchLastSection(True)
        self.alias_table.setMaximumHeight(190)
        results_layout.addWidget(self.alias_table)

        # Bootstrap FAP controls
        bsfap_row = QHBoxLayout()
        bsfap_row.addWidget(QLabel("Bootstrap FAP (LS only):"))
        self.spin_bootstrap_n = QSpinBox()
        self.spin_bootstrap_n.setRange(100, 10000)
        self.spin_bootstrap_n.setSingleStep(100)
        self.spin_bootstrap_n.setValue(1000)
        self.spin_bootstrap_n.setToolTip("Monte-Carlo 반복 수 (많을수록 정확, 느림)")
        bsfap_row.addWidget(self.spin_bootstrap_n)
        self.btn_bootstrap = QPushButton("Bootstrap FAP 계산")
        style_button(self.btn_bootstrap, height=Tokens.H_BUTTON)
        self.btn_bootstrap.setToolTip("분석 완료 후 LS 피크에 대한 Bootstrap FAP를 계산합니다")
        self.btn_bootstrap.setEnabled(False)
        self.btn_bootstrap.clicked.connect(self._run_bootstrap_fap)
        bsfap_row.addWidget(self.btn_bootstrap)
        self.bootstrap_progress = QLabel("")
        bsfap_row.addWidget(self.bootstrap_progress)
        bsfap_row.addStretch()
        results_layout.addLayout(bsfap_row)

        self.tabs.addTab(results_tab, "Results")

        # ── Right: plots/results take everything that's left ──
        root.addWidget(self.tabs, 1)

        self._scan_available_data()
        QTimer.singleShot(0, self._apply_responsive_shell)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if not hasattr(self, "control_column") or self._responsive_shell_pending:
            return
        self._responsive_shell_pending = True
        QTimer.singleShot(0, self._apply_responsive_shell)

    def _apply_responsive_shell(self) -> None:
        self._responsive_shell_pending = False
        if not hasattr(self, "control_column"):
            return
        compact = self._uses_compact_shell(self.width())
        if compact == self._shell_compact:
            return
        self._shell_compact = compact
        self._set_control_column_visible(not compact)

    @staticmethod
    def _uses_compact_shell(window_width: int) -> bool:
        return int(window_width) < 1350

    def _set_control_column_visible(self, visible: bool) -> None:
        self.control_column.setVisible(bool(visible))
        self.btn_controls.setChecked(bool(visible))
        self.btn_controls.setToolTip(
            "Hide data and period-search controls."
            if visible
            else "Show data and period-search controls."
        )

    def _toggle_control_column(self) -> None:
        self._set_control_column_visible(not self.control_column.isVisible())

    # ------------------------------------------------------------------
    # Data scanning / loading
    # ------------------------------------------------------------------

    def _scan_available_data(self):
        """Auto-select best lightcurve and populate filter combo."""
        result_dir = Path(self.params.P.result_dir)
        target_id = self.target_id_spin.value()

        lc_path = find_best_lightcurve_csv(result_dir, target_id)
        self._auto_lc_path = lc_path

        if lc_path and lc_path.exists():
            mode_key, mode_label = _detect_corr_mode(lc_path.name)
            pref = load_detrend_preference(result_dir, target_id)
            if pref:
                mode_label = corr_mode_label(pref, mode_label)
            self.source_label.setText(mode_label)

            # Scan filters
            filters_found: set[str] = set()
            try:
                df_head = pd.read_csv(lc_path, nrows=500)
                self.photometry_source_label.setText(
                    format_photometry_provenance(
                        summarize_photometry_table(df_head)
                    )
                )
                if "filter" in df_head.columns:
                    for flt in df_head["filter"].dropna().astype(str).unique():
                        fkey = _normalize_filter_key(flt)
                        if fkey and fkey.lower() != "nan":
                            filters_found.add(fkey)
            except Exception:
                self.photometry_source_label.setText(
                    "Photometry: Unknown | MAG: unknown"
                )
        else:
            self.source_label.setText("No data")
            self.photometry_source_label.setText(
                "Photometry: Unknown | MAG: unknown"
            )
            filters_found = set()

        self.filter_combo.blockSignals(True)
        self.filter_combo.clear()
        if filters_found:
            ordered_filters = sorted(filters_found)
            if len(ordered_filters) > 1:
                self.filter_combo.addItem(_ALL_FILTER_LABEL, ALL_FILTER_KEY)
            for flt in ordered_filters:
                self.filter_combo.addItem(flt, flt)
        else:
            self.filter_combo.addItem("(no data)")
        self.filter_combo.blockSignals(False)
        self.current_filter = self._selected_filter_key()

    def _selected_filter_key(self) -> str:
        value = self.filter_combo.currentData()
        if value is None:
            value = self.filter_combo.currentText()
        return str(value)

    def _on_filter_changed(self, index: int):
        self.current_filter = self._selected_filter_key()
        if self._ui_ready:
            self._load_light_curve(silent=True)

    def _on_target_id_changed(self, value: int):
        if self._ui_ready:
            # Re-scan sources for this target ID
            self._scan_available_data()
            self._load_light_curve(silent=True)

    def _set_data_status(self, text: str, ok: bool = False):
        color = Tokens.OK if ok else Tokens.ERROR
        self.data_status.setText(text)
        self.data_status.setStyleSheet(f"QLabel {{ color: {color}; }}")

    def _clear_analysis_results(self):
        self.results = {}
        self.alias_analysis = None
        self.multimode_diagnostic = None
        self.variable_analysis_bundle = None
        self.multinight = None
        self.results_table.setRowCount(0)
        self.alias_table.setRowCount(0)
        self.phase_period_combo.blockSignals(True)
        self.phase_period_combo.clear()
        self.phase_period_combo.blockSignals(False)
        self.summary_canvas.figure.clear()
        self.summary_canvas.draw_idle()
        self.periodogram_canvas.figure.clear()
        self.periodogram_canvas.draw_idle()
        self.phase_canvas.figure.clear()
        self.phase_canvas.draw_idle()
        self.progress_label.setText("")
        self.result_callout.setVisible(False)
        self.mode_callout.setVisible(False)
        self.btn_open_varstar.setEnabled(False)
        self._check_star_plot_cache_key = None
        self._check_star_plot_cache = None
        self._check_star_ls_cache.clear()

    def _load_light_curve(self, silent: bool = False):
        target_id = self.target_id_spin.value()
        flt = self._selected_filter_key()

        if not flt or flt == "(no data)":
            self.lc_data = None
            self.btn_run.setEnabled(False)
            self._clear_analysis_results()
            self._set_data_status("No light curve data available.", ok=False)
            if not silent:
                QMessageBox.warning(self, "No Data", "No light curve data available.")
            return

        lc_file = getattr(self, "_auto_lc_path", None)
        if not lc_file or not lc_file.exists():
            self.lc_data = None
            self.btn_run.setEnabled(False)
            self._clear_analysis_results()
            self._set_data_status("No data source found.", ok=False)
            return

        self.log(f"Loading: {lc_file}")

        try:
            self.lc_data = load_period_lightcurve_csv(lc_file, flt, target_id)
            self.current_filter = flt
            self._clear_analysis_results()
            n_valid = np.sum(np.isfinite(self.lc_data["time"]) & np.isfinite(self.lc_data["mag_raw"]))
            corr_info = f"corr={self.lc_data.get('col_corr')}" if self.lc_data.get("col_corr") else "corr=없음"
            err_info = f"err={self.lc_data.get('col_err')}" if self.lc_data.get("col_err") else "err=없음"
            corr_mode_label = self.lc_data.get("corr_mode_label", "Unknown")
            self._set_data_status(
                f"{n_valid}점  [{lc_file.name}]\n"
                f"Detrend: {corr_mode_label}  |  raw={self.lc_data.get('col_raw')}  {corr_info}\n"
                f"{err_info}  time={self.lc_data.get('col_time')}",
                ok=True,
            )
            self.btn_run.setEnabled(True)

            self.log(f"Loaded {self.lc_data.get('n_rows', 0)} rows, columns: {self.lc_data.get('columns', [])}")
            self.log(
                f"Time: {self.lc_data.get('col_time')}, Raw: {self.lc_data.get('col_raw')}, "
                f"Corr: {self.lc_data.get('col_corr')}, Err: {self.lc_data.get('col_err')}"
            )
            self.log(f"Filter: {flt}, Target ID: {target_id}, Valid points: {n_valid}, Detrend: {corr_mode_label}")

        except Exception as e:
            self.lc_data = None
            self.btn_run.setEnabled(False)
            self._clear_analysis_results()
            self._set_data_status(f"Load error: {e}", ok=False)
            if not silent:
                QMessageBox.warning(self, "Load Error", str(e))
            self.log(f"[ERROR] {e}")

    # ------------------------------------------------------------------
    # Analysis
    # ------------------------------------------------------------------

    def _run_analysis(self):
        if self.lc_data is None:
            QMessageBox.warning(self, "No Data", "Load light curve data first.")
            return

        if self.worker is not None and self.worker.isRunning():
            return

        min_period = self.min_period_spin.value()
        max_period = self.max_period_spin.value()
        samples = self.samples_spin.value()

        if min_period >= max_period:
            QMessageBox.warning(self, "Invalid Range", "Min period must be less than max period.")
            return

        methods = []
        if self.chk_ls.isChecked():
            methods.append("ls")
        if self.chk_pdm.isChecked():
            methods.append("pdm")
        if self.chk_bls.isChecked():
            methods.append("bls")
        if not methods:
            QMessageBox.warning(self, "No Method", "Select at least one method.")
            return

        if self._shell_compact and self.control_column.isVisible():
            self._set_control_column_visible(False)

        self.btn_run.setEnabled(False)
        self.btn_open_varstar.setEnabled(False)
        self.progress_label.setText("Computing...")

        self.worker = PeriodAnalysisWorker(
            time=self.lc_data["time"],
            mag_raw=self.lc_data["mag_raw"],
            mag_corr=self.lc_data["mag_corr"],
            mag_err=self.lc_data["mag_err"],
            min_period=min_period,
            max_period=max_period,
            samples_per_peak=samples,
            methods=methods,
            pdm_n_bins=self.pdm_bins_spin.value(),
            night_id=self.lc_data.get("night_id"),
            correction_preserves_nightly_baseline=self.lc_data.get(
                "correction_preserves_nightly_baseline", True
            ),
            include_alias_diagnostics=True,
            include_multimode_diagnostic=True,
        )
        self.worker.progress.connect(self._on_progress)
        self.worker.finished.connect(self._on_finished)
        self.worker.error.connect(self._on_error)
        self.worker.start()

    def _on_progress(self, msg: str):
        self.progress_label.setText(msg)
        self.log(msg)

    def _on_error(self, msg: str):
        self.btn_run.setEnabled(True)
        self.progress_label.setText("Error")
        QMessageBox.warning(self, "Error", msg)
        self.log(f"[ERROR] {msg}")

    def _on_finished(self, results: dict):
        self.btn_run.setEnabled(True)
        self.progress_label.setText("Rendering and saving...")
        # Keep diagnostics separate from per-method periodogram result dicts.
        self.alias_analysis = results.pop("alias_analysis", None)
        self.multimode_diagnostic = results.pop("multimode_diagnostic", None)
        self.multinight = results.pop("multinight", None)
        self.results = results

        self._show_result_callout()
        self._show_multimode_callout()
        self._update_summary_plot()
        self._update_periodogram_plot()
        self._update_results_table()
        self._update_alias_candidate_table()
        self._populate_phase_periods()
        self._update_phase_plot()
        self._save_results()
        self._log_alias_warnings()
        self.log("Analysis complete")
        self.progress_label.setText("Done")
        self.btn_open_varstar.setEnabled(
            bool(self.results)
            and self.variable_analysis_bundle is not None
            and self.variable_analysis_bundle.can_launch
        )

        # Enable Bootstrap FAP if LS result is available
        has_ls = any("ls" in k for k in results if "error" not in results.get(k, {}))
        self.btn_bootstrap.setEnabled(has_ls)

    def _bundle_filter_keys(self) -> list[str]:
        if self.lc_data is None:
            return []
        selected = str(self.lc_data.get("filter", ""))
        if selected and selected != ALL_FILTER_KEY:
            key = _normalize_filter_key(selected) or selected
            return [key]
        values = np.asarray(self.lc_data.get("filter_values", []), dtype=str)
        return sorted({
            _normalize_filter_key(value) or str(value).strip()
            for value in values
            if str(value).strip() and str(value).strip().lower() != "nan"
        })

    def _build_comparison_provenance(self, filter_keys: list[str]) -> tuple[dict, list[str]]:
        selection_dir = step8_selection_dir(Path(self.params.P.result_dir))
        filters: dict[str, dict] = {}
        release_reasons: list[str] = []
        for filt in filter_keys:
            stability_path = selection_dir / f"comparison_stability_{filt}.json"
            stability: dict = {}
            # Either half's record of the same decision. The window writes
            # `selection_<filter>.json`, the pipeline step writes
            # `lc_target_selection.json`; reading only the first held the
            # release with "comparison selection metadata is missing" on every
            # batch-built workspace — beside a stability report that was there.
            selection = read_step8_selection(self.params.P.result_dir, filt) or {}
            selection_path = Path(selection.get("source_file", "")) if selection                 else selection_dir / f"selection_{filt}.json"
            try:
                if stability_path.exists():
                    stability = json.loads(stability_path.read_text(encoding="utf-8"))
            except Exception:
                stability = {}

            comparison_ids = [int(value) for value in selection.get("comparison_ids", [])]
            comparison_source_ids = [
                int(value) for value in selection.get("comparison_source_ids", [])
            ]
            check_id = selection.get("check_id")
            check_source_id = selection.get("check_source_id")
            check_metrics = dict(stability.get("check_metrics") or {})
            filters[filt] = {
                "selection_file": str(selection_path) if selection else "",
                "stability_file": str(stability_path) if stability_path.exists() else "",
                "photometry_source": str(stability.get("photometry_source", "unknown")),
                "comparison_ids": comparison_ids,
                "comparison_source_ids": comparison_source_ids,
                "check_id": int(check_id) if check_id is not None else None,
                "check_source_id": int(check_source_id) if check_source_id is not None else None,
                "check_metrics": check_metrics,
                "selection_timestamp": str(selection.get("timestamp", "")),
                "stability_timestamp": str(stability.get("timestamp", "")),
            }

            if not selection:
                release_reasons.append(f"{filt}: comparison selection metadata is missing")
            elif not comparison_ids:
                release_reasons.append(f"{filt}: no comparison ensemble is selected")
            if check_id is None and check_source_id is None:
                release_reasons.append(f"{filt}: no check star is selected")
            if not stability:
                release_reasons.append(f"{filt}: comparison stability metadata is missing")
            elif not check_metrics:
                release_reasons.append(f"{filt}: check-star stability metrics are missing")
            elif int(check_metrics.get("n", 0) or 0) < 10:
                release_reasons.append(f"{filt}: check-star stability has fewer than 10 points")

        if not filter_keys:
            release_reasons.append("No analysis filter is available for Main QC release")
        return {
            "selection_dir": str(selection_dir),
            "requested_filters": list(filter_keys),
            "filters": filters,
        }, release_reasons

    def _build_check_star_qc(self, comparison_provenance: dict) -> tuple[dict, list[str]]:
        check = self._check_star_plot_data()
        reasons: list[str] = []
        ids_by_filter = {
            filt: entry.get("check_id")
            for filt, entry in dict(comparison_provenance.get("filters") or {}).items()
            if entry.get("check_id") is not None
        }
        if check is None:
            reasons.append("No check-star light curve is available")
            return {
                "status": "BLOCKED",
                "available": False,
                "check_ids_by_filter": ids_by_filter,
                "n_points": 0,
                "period_test": {},
            }, reasons

        period_test: dict[str, object] = {}
        if self._summary_method_result("ls") is not None:
            check_ls = self._check_star_ls_result(check)
            if "error" in check_ls:
                reasons.append(f"Check-star period test failed: {check_ls['error']}")
                period_test = {"error": str(check_ls["error"])}
            else:
                frequency = np.asarray(check_ls.get("frequency", []), dtype=float)
                power = np.asarray(check_ls.get("power", []), dtype=float)
                adopted = float(self._summary_period())
                candidate_power = None
                if (
                    np.isfinite(adopted)
                    and adopted > 0
                    and len(frequency) == len(power)
                    and len(frequency) > 0
                ):
                    candidate_index = int(np.nanargmin(np.abs(frequency - 1.0 / adopted)))
                    candidate_power = float(power[candidate_index])
                period_test = {
                    "best_period": float(check_ls.get("best_period", np.nan)),
                    "best_power": float(check_ls.get("best_power", np.nan)),
                    "fap": float(check_ls.get("fap", np.nan)),
                    "power_at_target_candidate": candidate_power,
                    "target_candidate_period": adopted,
                }

        filters = sorted({str(value) for value in check.get("filters", [])})
        return {
            "status": "PASS" if not reasons else "BLOCKED",
            "available": True,
            "check_id": check.get("check_id"),
            "check_ids_by_filter": ids_by_filter,
            "n_points": int(len(check["time"])),
            "filters": filters,
            "period_test": period_test,
        }, reasons

    def _build_variable_star_handoff(self) -> ValidatedLightCurveBundle:
        if self.lc_data is None or not self.results:
            raise ValueError("Run period analysis before opening the Variable Star Tool.")

        dtype = self._summary_data_type()
        scan_results: dict[str, dict] = {}
        for method in ("ls", "pdm", "bls"):
            result = self._summary_method_result(method)
            if result is not None:
                compact_result = dict(result)
                for transient_key in ("time", "mag", "mag_err"):
                    compact_result.pop(transient_key, None)
                scan_results[f"raw_{method}"] = compact_result

        if not scan_results:
            raise ValueError("No usable periodogram result is available for handoff.")

        filter_key = str(self.lc_data.get("filter", ""))
        analysis_filter = "__all__" if filter_key == ALL_FILTER_KEY else filter_key
        mag_col = (
            self.lc_data.get("col_corr")
            if dtype == "corr"
            else self.lc_data.get("col_raw")
        ) or self.lc_data.get("col_raw")
        input_mag = self.lc_data.get("mag_corr") if dtype == "corr" else self.lc_data.get("mag_raw")
        input_time = np.asarray(self.lc_data.get("time", []), dtype=float)
        input_mag = np.asarray(input_mag, dtype=float)
        input_valid = np.isfinite(input_time) & np.isfinite(input_mag)
        input_err = self.lc_data.get("mag_err")
        if input_err is not None:
            input_err = np.asarray(input_err, dtype=float)
            input_valid &= np.isfinite(input_err) & (input_err > 0)

        filter_keys = self._bundle_filter_keys()
        comparison_provenance, comparison_reasons = self._build_comparison_provenance(
            filter_keys
        )
        check_star_qc, check_reasons = self._build_check_star_qc(comparison_provenance)
        release_reasons = comparison_reasons + check_reasons
        main_qc = {
            "status": "PASS" if not release_reasons else "BLOCKED",
            "owner": "main_workflow",
            "check_star": check_star_qc,
            "reasons": list(release_reasons),
        }
        photometry_provenance = {
            "source": str(self.lc_data.get("photometry_source", "unknown")),
            "mag_input_column": str(self.lc_data.get("mag_input_column", "unknown")),
            "mag_error_input_column": str(
                self.lc_data.get("mag_error_input_column", "unknown")
            ),
            "series_mode": dtype,
            "correction_mode": str(self.lc_data.get("corr_mode", "")),
            "correction_mode_label": str(self.lc_data.get("corr_mode_label", "")),
        }
        source_file = Path(str(self.lc_data.get("source_file", "")))
        input_signature = {
            "n_points": int(np.count_nonzero(input_valid)),
            "time_min": float(np.nanmin(input_time[input_valid])) if np.any(input_valid) else None,
            "time_max": float(np.nanmax(input_time[input_valid])) if np.any(input_valid) else None,
        }
        if source_file.is_file():
            input_signature.update(compute_file_fingerprint(source_file))

        return ValidatedLightCurveBundle(
            workspace_dir=str(Path(self.params.P.result_dir)),
            source_file=str(source_file),
            target_id=int(self.lc_data.get("target_id", 0)),
            analysis_filter=analysis_filter,
            series_mode=dtype,
            mag_col=str(mag_col or ""),
            correction_mode=str(self.lc_data.get("corr_mode", "")),
            correction_preserves_nightly_baseline=bool(
                self.lc_data.get("correction_preserves_nightly_baseline", True)
            ),
            input_signature=input_signature,
            adopted_period=float(self._summary_period()),
            scan_results=scan_results,
            alias_analysis=dict(self.alias_analysis or {}),
            multimode_diagnostic=dict(self.multimode_diagnostic or {}),
            search={
                "min_period": float(self.min_period_spin.value()),
                "max_period": float(self.max_period_spin.value()),
                "samples_per_peak": int(self.samples_spin.value()),
                "pdm_bins": int(self.pdm_bins_spin.value()),
                "methods": [key.split("_", 1)[1] for key in scan_results],
            },
            release_status="APPROVED" if not release_reasons else "BLOCKED",
            release_reasons=release_reasons,
            main_qc=main_qc,
            comparison_provenance=comparison_provenance,
            photometry_provenance=photometry_provenance,
        )

    def _open_variable_star_tool(self):
        try:
            handoff = self.variable_analysis_bundle or self._build_variable_star_handoff()
            if not handoff.can_launch:
                raise ValueError(
                    "Main workflow did not release this light curve for advanced analysis: "
                    f"{handoff.release_message}"
                )
            launcher = getattr(self.main_window, "open_variable_star_tool", None)
            if not callable(launcher):
                raise RuntimeError("Variable Star Tool launcher is unavailable.")
            launcher(handoff=handoff)
            self.log(
                "Opened Variable Star Tool with "
                f"P={handoff.adopted_period:.8f} d, "
                f"filter={handoff.analysis_filter}, release={handoff.release_status}."
            )
        except Exception as exc:
            QMessageBox.warning(self, "Variable Star Tool", str(exc))
            self.log(f"[ERROR] Variable Star Tool handoff failed: {exc}")

    def _show_result_callout(self):
        """Surface the adopted candidate without hiding unresolved aliases."""
        t = Tokens
        analysis = self.alias_analysis
        if analysis and np.isfinite(analysis.get("adopted_period", np.nan)):
            status = str(analysis.get("status", "AMBIGUOUS")).upper()
            period = float(analysis["adopted_period"])
            reason = str(analysis.get("reason", "")).strip()
            text = f"{status}  |  candidate period {period:.6f} d"
            if reason:
                text += f"\n{reason}"
            if status == "RESOLVED":
                bg, fg = "#E8F5E9", t.OK
            elif status == "INSUFFICIENT":
                bg, fg = "#FDECEC", t.ERROR
            else:
                bg, fg = "#FFF4E5", t.WARN
            self.log(
                f"[ALIAS {status}] P={period:.8f} d; "
                f"window peaks={len(analysis.get('window_peaks', []))}; {reason}"
            )
        else:
            # single-night: show the strongest LS/PDM period
            best = None
            for key, data in self.results.items():
                bp = data.get("best_period", np.nan) if isinstance(data, dict) else np.nan
                if np.isfinite(bp):
                    best = bp
                    break
            if best is None:
                self.result_callout.setVisible(False)
                return
            text = f"Best period  {best:.6f} d   (single night)"
            bg, fg = t.SURFACE_ALT, t.TEXT
        self.result_callout.setText(text)
        self.result_callout.setStyleSheet(
            f"QLabel {{ background: {bg}; color: {fg}; font-weight: 600; "
            f"padding: {t.S2}px {t.S3}px; border-radius: {t.RADIUS_SM}px; }}"
        )
        self.result_callout.setVisible(True)

    def _show_multimode_callout(self):
        diagnostic = self.multimode_diagnostic or {}
        status = str(diagnostic.get("status", "")).upper()
        if not status:
            self.mode_callout.setVisible(False)
            return

        reason = str(diagnostic.get("reason", "")).strip()
        text = f"MODE DIAGNOSTIC  |  {status}"
        if reason:
            text += f"\n{reason}"
        if status == "MULTIMODE-SUSPECT":
            bg, fg = "#FCE8E6", "#B3261E"
        elif status == "SINGLE-COMPATIBLE":
            bg, fg = "#E8F5E9", Tokens.OK
        else:
            bg, fg = "#FFF4E5", Tokens.WARN
        self.mode_callout.setText(text)
        self.mode_callout.setStyleSheet(
            f"QLabel {{ background: {bg}; color: {fg}; font-weight: 600; "
            f"padding: {Tokens.S2}px {Tokens.S3}px; border-radius: {Tokens.RADIUS_SM}px; }}"
        )
        self.mode_callout.setVisible(True)
        self.log(f"[MODE {status}] {reason}")

    def _run_bootstrap_fap(self):
        """Run Bootstrap FAP for the best LS result in a background thread."""
        from apex.analysis.light_curve.period_analysis_service import bootstrap_fap

        # Find best LS result
        ls_key = next(
            (k for k in self.results if "ls" in k and "error" not in self.results[k]),
            None,
        )
        if ls_key is None:
            QMessageBox.information(self, "Bootstrap FAP", "LS 결과가 없습니다.")
            return

        data = self.results[ls_key]
        best_power = float(data["best_power"])
        time_arr   = data.get("time")
        mag_arr    = data.get("mag")
        mag_err    = data.get("mag_err")

        if time_arr is None or mag_arr is None:
            QMessageBox.warning(self, "Bootstrap FAP", "LS 데이터가 없습니다. 분석을 다시 실행하세요.")
            return

        n_bootstrap = int(self.spin_bootstrap_n.value())
        min_period  = float(self.min_period_spin.value())
        max_period  = float(self.max_period_spin.value())
        spp         = int(self.samples_spin.value())

        self.btn_bootstrap.setEnabled(False)
        self.bootstrap_progress.setText("계산 중…")

        class _BootstrapWorker(QThread):
            finished = pyqtSignal(float)
            progress = pyqtSignal(int, int)
            error    = pyqtSignal(str)

            def __init__(self, t, m, e, bp, minp, maxp, spp_, n):
                super().__init__()
                self._t, self._m, self._e = t, m, e
                self._bp = bp
                self._minp, self._maxp, self._spp = minp, maxp, spp_
                self._n = n

            def run(self):
                try:
                    fap = bootstrap_fap(
                        self._t, self._m, self._e,
                        self._bp,
                        self._minp, self._maxp,
                        samples_per_peak=self._spp,
                        n_bootstrap=self._n,
                        progress_cb=self.progress.emit,
                    )
                    self.finished.emit(fap)
                except Exception as exc:
                    self.error.emit(str(exc))

        self._bsworker = _BootstrapWorker(
            time_arr, mag_arr, mag_err,
            best_power, min_period, max_period, spp, n_bootstrap,
        )
        self._bsworker.progress.connect(
            lambda cur, tot: self.bootstrap_progress.setText(f"{cur}/{tot}")
        )
        self._bsworker.finished.connect(self._on_bootstrap_done)
        self._bsworker.error.connect(
            lambda e: (self.bootstrap_progress.setText("오류"), self.btn_bootstrap.setEnabled(True))
        )
        self._bsworker.start()

    def _on_bootstrap_done(self, fap: float):
        self.btn_bootstrap.setEnabled(True)
        fap_str = f"{fap:.4f}" if np.isfinite(fap) else "N/A"
        self.bootstrap_progress.setText(f"Bootstrap FAP = {fap_str}")
        self.log(f"[Bootstrap FAP] {fap_str}  (LS best peak)")

        # Update FAP column in results table for LS rows
        for row in range(self.results_table.rowCount()):
            method_item = self.results_table.item(row, 0)
            if method_item and "LS" in method_item.text():
                self.results_table.setItem(
                    row, 4,
                    QTableWidgetItem(f"{fap:.2e} (BS)" if np.isfinite(fap) else "-"),
                )

    # ------------------------------------------------------------------
    # Plots
    # ------------------------------------------------------------------

    @staticmethod
    @staticmethod
    @staticmethod
    @staticmethod
    def _on_summary_canvas_resize(self, event) -> None:
        compact = self._summary_uses_compact_layout(event.width)
        stacked = self._summary_uses_stacked_layout(event.width)
        unchanged = (
            compact == self._summary_layout_compact
            and stacked == self._summary_layout_stacked
        )
        if unchanged or self._summary_resize_pending:
            return
        self._summary_resize_pending = True
        QTimer.singleShot(0, self._refresh_summary_layout_after_resize)

    def _refresh_summary_layout_after_resize(self) -> None:
        self._summary_resize_pending = False
        if not self.results or self.lc_data is None:
            return
        width = self.summary_canvas.width()
        compact = self._summary_uses_compact_layout(width)
        stacked = self._summary_uses_stacked_layout(width)
        if (
            compact != self._summary_layout_compact
            or stacked != self._summary_layout_stacked
        ):
            self._update_summary_plot()

    @staticmethod
    def _periodogram_uses_compact_layout(canvas_width: int) -> bool:
        """Stack method/source panels when a two-column grid would collapse."""
        return int(canvas_width) < 620

    def _on_periodogram_canvas_resize(self, event) -> None:
        compact = self._periodogram_uses_compact_layout(event.width)
        if (
            compact == self._periodogram_layout_compact
            or self._periodogram_resize_pending
        ):
            return
        self._periodogram_resize_pending = True
        QTimer.singleShot(0, self._refresh_periodogram_layout_after_resize)

    def _refresh_periodogram_layout_after_resize(self) -> None:
        self._periodogram_resize_pending = False
        if not self.results:
            return
        compact = self._periodogram_uses_compact_layout(
            self.periodogram_canvas.width()
        )
        if compact != self._periodogram_layout_compact:
            self._update_periodogram_plot()

    def _update_periodogram_plot(self):
        fig = self.periodogram_canvas.figure
        fig.clear()

        if not self.results:
            self.periodogram_canvas.draw_idle()
            return

        method_labels = {"ls": "Lomb-Scargle", "pdm": "PDM (1-\u0398)", "bls": "BLS"}
        data_labels = {"raw": "Raw", "corr": "Corrected"}
        method_colors = {"ls": "#1E88E5", "pdm": "#E53935", "bls": "#FF9800"}
        y_labels = {"ls": "LS Power", "pdm": "1 - \u0398", "bls": "BLS Power"}

        methods_present = []
        data_types_present = []
        for key in self.results:
            parts = key.split("_", 1)
            if len(parts) == 2:
                dt, mt = parts
                if mt not in methods_present:
                    methods_present.append(mt)
                if dt not in data_types_present:
                    data_types_present.append(dt)

        panel_specs = [
            (method, dtype)
            for method in methods_present
            for dtype in data_types_present
        ]
        compact = self._periodogram_uses_compact_layout(
            self.periodogram_canvas.width()
        )
        self._periodogram_layout_compact = compact
        if compact:
            n_rows = len(panel_specs) or 1
            n_cols = 1
            self.periodogram_canvas.setMinimumHeight(max(520, n_rows * 190))
        else:
            n_rows = len(methods_present) or 1
            n_cols = len(data_types_present) or 1
            self.periodogram_canvas.setMinimumHeight(400)
        axes = fig.subplots(n_rows, n_cols, squeeze=False)

        for panel_index, (method, dtype) in enumerate(panel_specs):
            if compact:
                ax = axes[panel_index][0]
            else:
                ri = methods_present.index(method)
                ci = data_types_present.index(dtype)
                ax = axes[ri][ci]
            key = f"{dtype}_{method}"
            data = self.results.get(key)
            if data is None or "error" in (data or {}):
                err_msg = data.get("error", "No data") if data else "No data"
                ax.text(
                    0.5, 0.5, err_msg, ha="center", va="center",
                    transform=ax.transAxes, fontsize=9,
                )
                ax.set_title(
                    f"{data_labels.get(dtype, dtype)} / "
                    f"{method_labels.get(method, method)}"
                )
                continue

            power = data["power"]
            best_period = data["best_period"]
            best_power = data["best_power"]

            if "frequency" in data:
                periods = 1.0 / data["frequency"]
            elif "trial_periods" in data:
                periods = data["trial_periods"]
            else:
                ax.text(
                    0.5, 0.5, "No period axis", ha="center", va="center",
                    transform=ax.transAxes,
                )
                continue

            color = method_colors.get(method, "#666")
            ax.plot(periods, power, color=color, lw=0.8, alpha=0.8)
            ax.axvline(
                best_period, color="red", ls="--", lw=1.5, alpha=0.8,
                label=f"P={best_period:.6f}d",
            )
            ax.scatter([best_period], [best_power], color="red", s=50, zorder=5)

            ax.set_xlabel("Period (days)")
            ax.set_ylabel(y_labels.get(method, "Power"))
            ax.set_title(
                f"{data_labels.get(dtype, dtype)} / "
                f"{method_labels.get(method, method)}\n"
                f"P = {best_period:.6f} d"
            )
            ax.set_xscale("log")

            # Candidate family generated from the actual sampling window.
            if hasattr(self, "chk_alias") and self.chk_alias.isChecked():
                p_min = self.min_period_spin.value()
                p_max = self.max_period_spin.value()
                analysis = self.alias_analysis or {}
                for candidate in analysis.get("candidates", [])[1:6]:
                    ap = float(candidate.get("period", np.nan))
                    if np.isfinite(ap) and p_min <= ap <= p_max:
                        ax.axvline(
                            ap,
                            color="orange",
                            ls="--",
                            lw=1.0,
                            alpha=0.65,
                            label=(
                                f"window candidate {ap:.4f}d"
                                if int(candidate.get("rank", 0)) == 2
                                else None
                            ),
                        )

            ax.legend(loc="upper right", fontsize=7)
            ax.grid(True, alpha=0.3)

        if compact:
            fig.subplots_adjust(
                left=0.14, right=0.97, bottom=0.06, top=0.95, hspace=1.05
            )
        else:
            short_canvas = self.periodogram_canvas.height() < 520
            fig.subplots_adjust(
                left=0.08, right=0.98, bottom=0.08,
                top=0.84 if short_canvas else 0.90,
                wspace=0.35,
                hspace=0.90 if short_canvas else 0.60,
            )
        self.periodogram_canvas.draw_idle()

    def _update_results_table(self):
        self.results_table.setRowCount(0)

        if not self.results:
            return

        method_labels = {"ls": "Lomb-Scargle", "pdm": "PDM", "bls": "BLS"}
        data_labels = {"raw": "Raw", "corr": "Corrected"}

        # Collect all best periods per data type for cross-method alias detection
        best_periods_by_dtype: dict[str, dict[str, float]] = {}
        for key, data in self.results.items():
            if "error" in data:
                continue
            parts = key.split("_", 1)
            if len(parts) == 2:
                dt, mt = parts
                best_periods_by_dtype.setdefault(dt, {})[mt] = data["best_period"]

        for key, data in self.results.items():
            row = self.results_table.rowCount()
            self.results_table.insertRow(row)

            parts = key.split("_", 1)
            dtype = parts[0] if len(parts) == 2 else key
            method = parts[1] if len(parts) == 2 else ""

            self.results_table.setItem(row, 0, QTableWidgetItem(method_labels.get(method, method)))
            self.results_table.setItem(row, 1, QTableWidgetItem(data_labels.get(dtype, dtype)))

            if "error" in data:
                self.results_table.setItem(row, 2, QTableWidgetItem(data["error"]))
                continue

            self.results_table.setItem(row, 2, QTableWidgetItem(f"{data['best_period']:.6f}"))
            self.results_table.setItem(row, 3, QTableWidgetItem(f"{data['best_power']:.4f}"))

            fap = data.get("fap", np.nan)
            fap_str = f"{fap:.2e}" if np.isfinite(fap) else "-"
            self.results_table.setItem(row, 4, QTableWidgetItem(fap_str))

            # Compare method peaks against the measured sampling-window family.
            alias_tag = ""
            bp = data["best_period"]
            others = best_periods_by_dtype.get(dtype, {})
            for other_method, other_p in others.items():
                if other_method == method:
                    continue
                if self._periods_are_sampling_aliases(bp, other_p):
                    oml = method_labels.get(other_method, other_method.upper())
                    alias_tag = f"window alias of {oml} ({other_p:.4f}d)"
                    break
            item_alias = QTableWidgetItem(alias_tag)
            if alias_tag:
                item_alias.setForeground(Qt.red)
            self.results_table.setItem(row, 5, item_alias)

            top_periods = data.get("top_periods", [])[:3]
            top_str = ", ".join(f"{p:.4f}" for p in top_periods)
            self.results_table.setItem(row, 6, QTableWidgetItem(top_str))

    def _periods_are_sampling_aliases(self, period1: float, period2: float) -> bool:
        analysis = self.alias_analysis or {}
        return periods_are_window_aliases(
            float(period1),
            float(period2),
            analysis.get("window_peaks", []),
            float(analysis.get("baseline_days", 1.0)),
        )

    def _update_alias_candidate_table(self):
        self.alias_table.setRowCount(0)
        analysis = self.alias_analysis or {}
        for row_idx, candidate in enumerate(analysis.get("candidates", [])):
            self.alias_table.insertRow(row_idx)
            values = [
                str(candidate.get("rank", row_idx + 1)),
                f"{float(candidate.get('period', np.nan)):.8f}",
                f"{float(candidate.get('freq_cd', np.nan)):.5f}",
                f"{float(candidate.get('relative_power', np.nan)):.3f}",
                f"{float(candidate.get('delta_bic', np.nan)):.2f}",
                str(candidate.get("leave_one_out_votes", 0)),
                str(candidate.get("relation_to_best", "")),
            ]
            relation = str(candidate.get("relation_to_best", ""))
            for col_idx, value in enumerate(values):
                item = QTableWidgetItem(value)
                if relation == "adopted":
                    item.setForeground(Qt.darkGreen)
                elif relation == "window-alias":
                    item.setForeground(Qt.darkYellow)
                self.alias_table.setItem(row_idx, col_idx, item)

    def _populate_phase_periods(self):
        self.phase_period_combo.blockSignals(True)
        self.phase_period_combo.clear()

        method_labels = {"ls": "LS", "pdm": "PDM", "bls": "BLS"}
        data_labels = {"raw": "Raw", "corr": "Corr"}

        periods = []
        seen = set()
        analysis = self.alias_analysis or {}
        for candidate in analysis.get("candidates", [])[:5]:
            period = float(candidate.get("period", np.nan))
            if not np.isfinite(period) or period <= 0:
                continue
            relation = str(candidate.get("relation_to_best", "candidate"))
            periods.append((f"Alias rank {candidate.get('rank', '?')} ({relation}): {period:.6f} d", period))
            seen.add(round(period, 8))
        for key, data in self.results.items():
            if "error" in data:
                continue
            parts = key.split("_", 1)
            dtype = parts[0] if len(parts) == 2 else key
            method = parts[1] if len(parts) == 2 else ""

            ml = method_labels.get(method, method.upper())
            dl = data_labels.get(dtype, dtype)
            best_p = data.get("best_period", np.nan)
            if np.isfinite(best_p):
                tag = f"{ml}/{dl}"
                periods.append((f"{tag}: {best_p:.6f} d", best_p))
                p2 = round(best_p * 2, 8)
                ph = round(best_p / 2, 8)
                if p2 not in seen:
                    periods.append((f"{tag} x2: {p2:.6f} d", p2))
                    seen.add(p2)
                if ph not in seen:
                    periods.append((f"{tag} /2: {ph:.6f} d", ph))
                    seen.add(ph)

        for label, p in periods:
            self.phase_period_combo.addItem(label, p)

        if periods:
            self.phase_period_edit.setValue(periods[0][1])

        self.phase_period_combo.blockSignals(False)

    def _update_phase_plot(self, index: int = 0):
        if self.phase_period_combo.count() == 0:
            return
        period = self.phase_period_combo.currentData()
        if period is None or not np.isfinite(period) or period <= 0:
            return
        self.phase_period_edit.blockSignals(True)
        self.phase_period_edit.setValue(period)
        self.phase_period_edit.blockSignals(False)
        self._draw_phase_plot(period)

    def _update_phase_plot_custom(self):
        period = self.phase_period_edit.value()
        if period <= 0:
            return
        self._draw_phase_plot(period)

    def _draw_phase_plot(self, period: float):
        fig = self.phase_canvas.figure
        fig.clear()

        if not self.results:
            self.phase_canvas.draw_idle()
            return

        ax = fig.add_subplot(111)

        colors = {"raw": "#1E88E5", "corr": "#43A047"}
        markers = {"raw": "o", "corr": "s"}
        col_raw = self.lc_data.get("col_raw", "") if self.lc_data else ""
        col_corr = self.lc_data.get("col_corr", "") if self.lc_data else ""
        labels_map = {
            "raw": f"Raw ({col_raw})" if col_raw else "Raw",
            "corr": f"Corrected ({col_corr})" if col_corr else "Corrected",
        }
        phase_t0 = self._phase_plot_reference_t0()

        plotted_dtypes = set()
        for key, data in self.results.items():
            if "error" in data:
                continue
            parts = key.split("_", 1)
            dtype = parts[0] if len(parts) == 2 else key
            if dtype in plotted_dtypes:
                continue
            plotted_dtypes.add(dtype)

            t = data["time"]
            mag = data["mag"]
            mag_err = data.get("mag_err")

            if not np.isfinite(phase_t0):
                continue
            phase = ((t - phase_t0) / period) % 1.0
            phase_ext = np.concatenate([phase, phase + 1.0])
            mag_ext = np.concatenate([mag, mag])

            color = colors.get(dtype, "#666")
            marker = markers.get(dtype, "o")
            label = labels_map.get(dtype, dtype)

            if mag_err is not None and np.any(np.isfinite(mag_err)):
                err_ext = np.concatenate([mag_err, mag_err])
                ax.errorbar(
                    phase_ext, mag_ext, yerr=err_ext,
                    fmt=marker, color=color, markersize=4,
                    elinewidth=0.5, capsize=0, alpha=0.7,
                    label=label
                )
            else:
                ax.scatter(
                    phase_ext, mag_ext, c=color, marker=marker,
                    s=20, alpha=0.7, label=label
                )

        ax.invert_yaxis()
        ax.set_xlabel("Phase")
        ax.set_ylabel("Magnitude")
        src_name = Path(self.lc_data.get("source_file", "")).name if self.lc_data else ""
        ax.set_title(f"Phase Folded Light Curve  P = {period:.6f} d\n{src_name}", fontsize=9)
        ax.set_xlim(0, 2)
        ax.legend(loc="upper right")
        ax.grid(True, alpha=0.3)
        ax.axvline(0, color="gray", ls=":", alpha=0.5)
        ax.axvline(1, color="gray", ls=":", alpha=0.5)

        # Check star phase-folded overlay
        try:
            result_dir = Path(self.params.P.result_dir)
            _flt = self.lc_data.get("filter", "") if self.lc_data else ""
            _check_filter = None if _flt == ALL_FILTER_KEY else _flt
            _ck_id, _ck_df = _load_check_star_for_plot(result_dir, _check_filter)
            if _ck_df is not None and not _ck_df.empty:
                _t_col = self._phase_plot_check_time_column(_ck_df)
                _y_col = next((c for c in ["diff_mag_raw", "diff_mag", "mag"] if c in _ck_df.columns), None)
                if _t_col and _y_col:
                    if _check_filter and "filter" in _ck_df.columns:
                        _ck_df = _ck_df[
                            _ck_df["filter"].astype(str).map(_normalize_filter_key)
                            == _normalize_filter_key(_check_filter)
                        ]
                    _ct = pd.to_numeric(_ck_df[_t_col], errors="coerce").to_numpy(float)
                    if _t_col == "rel_time_hr":
                        _ct = _ct / 24.0
                    _cm = pd.to_numeric(_ck_df[_y_col], errors="coerce").to_numpy(float)
                    if _flt == ALL_FILTER_KEY and "filter" in _ck_df.columns:
                        _check_labels = (
                            _ck_df["filter"]
                            .astype(str)
                            .map(lambda value: _normalize_filter_key(value) or str(value).strip())
                            .to_numpy(dtype=str)
                        )
                        _cm = self._median_center_for_plot(_cm, _check_labels)
                    _mask = np.isfinite(_ct) & np.isfinite(_cm)
                    if _mask.any() and np.isfinite(phase_t0):
                        _ck_label = f"Check ID {_ck_id}" if _ck_id is not None else "Check"
                        _phase = ((_ct[_mask] - phase_t0) / period) % 1.0
                        _phase_ext = np.concatenate([_phase, _phase + 1.0])
                        _mag_ext = np.concatenate([_cm[_mask], _cm[_mask]])
                        ax.scatter(_phase_ext, _mag_ext, s=8, color="#FFD700", alpha=0.4,
                                   zorder=2, label=_ck_label, marker="^")
                        ax.legend(loc="upper right", fontsize=8)
        except Exception:
            pass

        fig.tight_layout()
        self.phase_canvas.draw_idle()

    def _phase_plot_reference_t0(self) -> float:
        if self.lc_data is not None:
            t = pd.to_numeric(pd.Series(self.lc_data.get("time", np.array([]))), errors="coerce").to_numpy(float)
            finite = t[np.isfinite(t)]
            if finite.size:
                return float(np.nanmin(finite))
        finite_chunks: list[np.ndarray] = []
        for data in self.results.values():
            if not isinstance(data, dict) or "error" in data or "time" not in data:
                continue
            t = pd.to_numeric(pd.Series(data["time"]), errors="coerce").to_numpy(float)
            finite = t[np.isfinite(t)]
            if finite.size:
                finite_chunks.append(finite)
        if not finite_chunks:
            return np.nan
        return float(np.nanmin(np.concatenate(finite_chunks)))

    # ------------------------------------------------------------------
    # Save / validate
    # ------------------------------------------------------------------

    def _log_alias_warnings(self):
        """Log warnings when methods select different sampling-window aliases."""
        if not self.results:
            return
        # Group best periods by data type
        by_dtype: dict[str, dict[str, float]] = {}
        for key, data in self.results.items():
            if "error" in data:
                continue
            parts = key.split("_", 1)
            if len(parts) == 2:
                dt, mt = parts
                by_dtype.setdefault(dt, {})[mt] = data["best_period"]

        method_labels = {"ls": "Lomb-Scargle", "pdm": "PDM", "bls": "BLS"}
        for dtype, methods in by_dtype.items():
            if len(methods) < 2:
                continue
            keys = list(methods.keys())
            for i in range(len(keys)):
                for j in range(i + 1, len(keys)):
                    m1, m2 = keys[i], keys[j]
                    p1, p2 = methods[m1], methods[m2]
                    if abs(p1 - p2) / max(p1, p2) < 0.005:
                        continue  # same period, no alias issue
                    if self._periods_are_sampling_aliases(p1, p2):
                        ml1 = method_labels.get(m1, m1.upper())
                        ml2 = method_labels.get(m2, m2.upper())
                        self.log(
                            f"[ALIAS WARNING] {dtype}: {ml1}={p1:.6f}d ↔ "
                            f"{ml2}={p2:.6f}d match the measured sampling window. "
                            "Method agreement alone does not resolve the alias."
                        )

    def _save_results(self):
        if not self.results or self.lc_data is None:
            return

        result_dir = Path(self.params.P.result_dir)
        out_dir = step11_period_dir(result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        flt = self.lc_data.get("filter", "unknown")
        target_id = self.lc_data.get("target_id", 0)

        summary_path = save_period_analysis_outputs(
            result_dir=result_dir,
            lc_data=self.lc_data,
            results=self.results,
            min_period=self.min_period_spin.value(),
            max_period=self.max_period_spin.value(),
            alias_analysis=self.alias_analysis,
            multimode_diagnostic=self.multimode_diagnostic,
        )
        self.log(f"Saved: {summary_path}")
        if hasattr(self, "summary_canvas") and self.summary_canvas.figure.axes:
            plot_path = out_dir / f"period_summary_{flt}_ID{target_id}.png"
            self.summary_canvas.figure.savefig(plot_path, dpi=160, bbox_inches="tight")
            self.log(f"Saved: {plot_path}")
        try:
            bundle = self._build_variable_star_handoff()
            bundle_path = bundle.write_json(out_dir / "variable_analysis_bundle.json")
            self.variable_analysis_bundle = bundle
            self.log(f"Saved: {bundle_path}")
            self.log(
                f"[MAIN QC] Variable analysis release={bundle.release_status}: "
                f"{bundle.release_message}"
            )
        except Exception as exc:
            self.variable_analysis_bundle = None
            self.log(f"[ERROR] Variable analysis bundle failed: {exc}")

    def log(self, msg: str):
        if self.log_text is not None:
            self.log_text.append(msg)

    def validate_step(self) -> bool:
        result_dir = Path(self.params.P.result_dir)
        out_dir = step11_period_dir(result_dir)
        return out_dir.exists() and any(out_dir.glob("period_analysis_*.json"))

    def save_state(self):
        state = {
            "min_period": self.min_period_spin.value(),
            "max_period": self.max_period_spin.value(),
            "samples_per_peak": self.samples_spin.value(),
            "pdm_bins": self.pdm_bins_spin.value(),
            "use_ls": self.chk_ls.isChecked(),
            "use_pdm": self.chk_pdm.isChecked(),
            "use_bls": self.chk_bls.isChecked(),
            "show_alias": self.chk_alias.isChecked(),
        }
        self.project_state.store_step_data("period_analysis", state)

    def restore_state(self):
        state = self.project_state.get_step_data("period_analysis")
        if not state:
            return
        if "min_period" in state:
            self.min_period_spin.setValue(float(state["min_period"]))
        if "max_period" in state:
            self.max_period_spin.setValue(float(state["max_period"]))
        if "samples_per_peak" in state:
            self.samples_spin.setValue(int(state["samples_per_peak"]))
        if "pdm_bins" in state:
            self.pdm_bins_spin.setValue(int(state["pdm_bins"]))
        if "use_ls" in state:
            self.chk_ls.setChecked(bool(state["use_ls"]))
        if "use_pdm" in state:
            self.chk_pdm.setChecked(bool(state["use_pdm"]))
        if "use_bls" in state:
            self.chk_bls.setChecked(bool(state["use_bls"]))
        if "show_alias" in state:
            self.chk_alias.setChecked(bool(state["show_alias"]))
