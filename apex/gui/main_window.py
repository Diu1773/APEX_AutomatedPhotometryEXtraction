"""
APEX Main Window — unified CMD + LC workflow.

Usage:
    MainWindowWorkflow(mode="cmd")  ← cluster photometry
    MainWindowWorkflow(mode="lc")   ← light curve analysis
"""

from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QGroupBox, QMessageBox, QFileDialog,
    QAction, QApplication, QLineEdit, QTextEdit,
    QPlainTextEdit, QComboBox, QSpinBox, QDoubleSpinBox
)
from PyQt5.QtCore import Qt, pyqtSignal, QObject, QEvent
from PyQt5.QtGui import QFont, QIcon, QPixmap, QPainter
from pathlib import Path
from typing import Optional, List

_RESOURCES = Path(__file__).resolve().parent.parent / "resources"


def _svg_to_pixmap(svg_path: Path, size: int = 256) -> QPixmap:
    pixmap = QPixmap(size, size)
    pixmap.fill(Qt.transparent)
    try:
        from PyQt5.QtSvg import QSvgRenderer
        renderer = QSvgRenderer(str(svg_path))
        painter = QPainter(pixmap)
        renderer.render(painter)
        painter.end()
    except Exception:
        pass
    return pixmap


def _apply_opacity(pixmap: QPixmap, opacity: float) -> QPixmap:
    result = QPixmap(pixmap.size())
    result.fill(Qt.transparent)
    painter = QPainter(result)
    painter.setOpacity(opacity)
    painter.drawPixmap(0, 0, pixmap)
    painter.end()
    return result


def _load_icon(mode: str) -> QIcon:
    svg_path = _RESOURCES / f"logo_{mode}.svg"
    if not svg_path.exists():
        svg_path = _RESOURCES / "logo_base.svg"
    if not svg_path.exists():
        return QIcon()
    px = _svg_to_pixmap(svg_path, 256)
    icon = QIcon()
    for size in (16, 24, 32, 48, 64, 128, 256):
        icon.addPixmap(px.scaled(size, size, Qt.KeepAspectRatio, Qt.SmoothTransformation))
    return icon


class StepButton(QPushButton):
    """Step button with completion/accessibility status indication."""

    def __init__(self, step_number: int, step_name: str, parent=None):
        super().__init__(parent)
        self.step_number = step_number
        self.step_name = step_name
        self.completed = False
        self.accessible = False
        self.setText(f"Step {step_number}: {step_name}")
        self.setMinimumHeight(50)
        self.setMinimumWidth(300)
        self.update_appearance()

    def set_completed(self, completed: bool):
        self.completed = completed
        self.update_appearance()

    def set_accessible(self, accessible: bool):
        self.accessible = accessible
        self.update_appearance()

    def update_appearance(self):
        if self.completed:
            self.setStyleSheet("""
                QPushButton {
                    background-color: #4CAF50; color: white;
                    font-size: 14px; font-weight: bold;
                    border: 2px solid #45a049; border-radius: 5px;
                    text-align: left; padding: 10px;
                }
                QPushButton:hover { background-color: #45a049; }
            """)
            self.setText(f"\u2713 Step {self.step_number}: {self.step_name}")
            self.setEnabled(True)
        elif self.accessible:
            self.setStyleSheet("""
                QPushButton {
                    background-color: #2196F3; color: white;
                    font-size: 14px; font-weight: bold;
                    border: 2px solid #1976D2; border-radius: 5px;
                    text-align: left; padding: 10px;
                }
                QPushButton:hover { background-color: #1976D2; }
            """)
            self.setText(f"\u25cb Step {self.step_number}: {self.step_name}")
            self.setEnabled(True)
        else:
            self.setStyleSheet("""
                QPushButton {
                    background-color: #E0E0E0; color: #999999;
                    font-size: 14px;
                    border: 2px solid #CCCCCC; border-radius: 5px;
                    text-align: left; padding: 10px;
                }
                QPushButton:disabled { background-color: #E0E0E0; color: #999999; }
            """)
            self.setText(f"\U0001f512 Step {self.step_number}: {self.step_name} (Locked)")
            self.setEnabled(False)


class ShortcutRouter(QObject):
    def __init__(self, main_window: "MainWindowWorkflow"):
        super().__init__(main_window)
        self.main_window = main_window

    @staticmethod
    def _is_text_input(widget) -> bool:
        if isinstance(widget, (QLineEdit, QTextEdit, QPlainTextEdit)):
            return True
        if isinstance(widget, QComboBox) and widget.isEditable():
            return True
        return False

    def eventFilter(self, obj, event):
        if event.type() != QEvent.KeyPress:
            return False
        if self._is_text_input(QApplication.focusWidget()):
            return False
        key = event.key()
        if key not in (Qt.Key_Period, Qt.Key_BracketLeft, Qt.Key_BracketRight):
            return False
        target = self.main_window.current_step_window
        if target is None:
            target = QApplication.activeWindow()
        if target is None:
            return False
        if key == Qt.Key_Period and hasattr(target, "cycle_filter"):
            target.cycle_filter()
            return True
        if key == Qt.Key_BracketLeft:
            if hasattr(target, "navigate_frame"):
                target.navigate_frame(-1)
                return True
            if hasattr(target, "step_frame"):
                target.step_frame(-1)
                return True
        if key == Qt.Key_BracketRight:
            if hasattr(target, "navigate_frame"):
                target.navigate_frame(1)
                return True
            if hasattr(target, "step_frame"):
                target.step_frame(1)
                return True
        return False


class MainWindowWorkflow(QMainWindow):
    """Unified APEX main window for CMD (cluster) and LC (light curve) modes."""

    step_requested = pyqtSignal(int)
    log_message = pyqtSignal(str)

    def __init__(self, mode: str = "lc", param_file: Optional[str] = None):
        super().__init__()
        if mode not in ("cmd", "lc"):
            raise ValueError(f"mode must be \'cmd\' or \'lc\', got {mode!r}")
        self.mode = mode

        try:
            if mode == "cmd":
                from apex.config.parameters_cmd import Parameters
            else:
                from apex.config.parameters_lc import Parameters
            param_path = Path(param_file) if param_file else Path("parameters.toml")
            self.params = Parameters(param_path)

            from apex.core import InstrumentConfig, FileManager, ProjectState
            self.instrument = InstrumentConfig(self.params)
            self.file_manager = FileManager(self.params)

            project_root = Path(__file__).parent.parent
            state_dir = project_root / ".state" / mode
            state_dir.mkdir(parents=True, exist_ok=True)
            self.project_state = ProjectState(state_dir)

            if mode == "lc":
                self._bootstrap_file_selection_state()

        except Exception as e:
            QMessageBox.critical(self, "Initialization Error",
                                 f"Failed to load parameters:\n{e}")
            raise

        if mode == "cmd":
            self.step_names = [
                "File Selection",
                "Image Crop",
                "Sky Preview & QC",
                "Source Detection",
                "Aperture Photometry",
                "PSF Photometry",
                "WCS Plate Solving",
                "Reference Catalog Build",
                "Star ID Matching",
                "Master ID Editor",
                "Zeropoint Calibration",
                "CMD Plot",
                "Isochrone Model",
            ]
        else:
            self.step_names = [
                "File Selection",
                "Image Crop",
                "Sky Preview & QC",
                "Source Detection",
                "Aperture Photometry",
                "WCS Plate Solving",
                "Reference Build",
                "Star ID Matching",
                "Target/Comparison Selection",
                "Light Curve Builder",
                "Detrend & Night Merge",
                "Period Analysis",
            ]

        self.project_state.assign_steps(self.step_names)
        self.step_buttons: List[StepButton] = []
        self.current_step_window = None
        self.varstar_window = None

        self.setup_ui()
        self.setup_menu()
        self.update_step_buttons()

        mode_label = "CMD Cluster Photometry" if mode == "cmd" else "LC Light Curve Analysis"
        self.append_log(f"APEX {mode_label} initialized")
        self.append_log(f"Project: {self.project_state.state['project_name']}")

        self._shortcut_router = ShortcutRouter(self)
        QApplication.instance().installEventFilter(self._shortcut_router)

        if mode == "lc" and hasattr(self, "_offline_data_dir"):
            QMessageBox.warning(
                self, "Previous Data Path Unavailable",
                f"The last-used data path is inaccessible:\n\n"
                f"  {self._offline_data_dir}\n\n"
                "External drive may be disconnected.\n"
                "Please set a data path in Step 1."
            )

    # ── LC file-selection bootstrap ──────────────────────────────────────────

    def _bootstrap_file_selection_state(self) -> None:
        state_data = self.project_state.get_step_data("file_selection")
        if not state_data:
            return
        data_dir = state_data.get("data_dir")
        if data_dir:
            self.params.P.data_dir = Path(data_dir)
            saved_result_dir = state_data.get("result_dir")
            self.params.P.result_dir = (
                Path(saved_result_dir) if saved_result_dir
                else self.params.P.data_dir / "result"
            )
            self.params.P.cache_dir = self.params.P.result_dir / "cache"
            try:
                self.params.P.result_dir.mkdir(parents=True, exist_ok=True)
                self.params.P.cache_dir.mkdir(parents=True, exist_ok=True)
            except OSError:
                self._offline_data_dir = str(data_dir)
        prefix = state_data.get("filename_prefix")
        if prefix:
            self.params.P.filename_prefix = prefix
        ref_frame = state_data.get("reference_frame")
        if ref_frame:
            self.file_manager.ref_filename = ref_frame
        multi_night = bool(state_data.get("multi_night"))
        root_dir = state_data.get("root_dir") or data_dir
        night_dirs = [Path(p) for p in state_data.get("night_dirs", []) if p]
        if multi_night and night_dirs:
            root_path = Path(root_dir) if root_dir else self.params.P.data_dir
            self.file_manager.set_multi_night_dirs(root_path, night_dirs)
        else:
            self.file_manager.clear_multi_night_dirs()

    # ── UI setup ─────────────────────────────────────────────────────────────

    def setup_ui(self):
        mode_title = "CMD Cluster Photometry" if self.mode == "cmd" else "Light Curve Analysis"
        self.setWindowTitle(f"APEX — {mode_title}")
        icon = _load_icon(self.mode)
        self.setWindowIcon(icon)
        QApplication.instance().setWindowIcon(icon)
        self.setMinimumSize(800, 700)

        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        # Watermark logo — bottom-right corner of the central widget
        wm_svg = _RESOURCES / f"logo_{self.mode}.svg"
        if not wm_svg.exists():
            wm_svg = _RESOURCES / "logo_base.svg"
        if wm_svg.exists():
            wm_px = _svg_to_pixmap(wm_svg, 180)
            wm_px = _apply_opacity(wm_px, 0.07)
            wm_label = QLabel(central)
            wm_label.setPixmap(wm_px)
            wm_label.setAttribute(Qt.WA_TransparentForMouseEvents)
            wm_label.setStyleSheet("background: transparent;")
            wm_label.resize(180, 180)
            central.resizeEvent = lambda e, lbl=wm_label: (
                lbl.move(e.size().width() - 190, e.size().height() - 190)
            )

        title = QLabel("Aperture Photometry Toolkit")
        title.setFont(QFont("Arial", 18, QFont.Bold))
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)

        subtitle = QLabel(f"KNUEMAO Observatory — {mode_title}")
        subtitle.setFont(QFont("Arial", 10))
        subtitle.setAlignment(Qt.AlignCenter)
        layout.addWidget(subtitle)

        settings_layout = QHBoxLayout()
        settings_layout.addStretch()
        btn_settings = QPushButton("\u2699 Instrument Settings")
        btn_settings.setFont(QFont("Arial", 11, QFont.Bold))
        btn_settings.setMinimumHeight(40)
        btn_settings.setStyleSheet("""
            QPushButton {
                background-color: #9C27B0; color: white;
                border: 2px solid #7B1FA2; border-radius: 5px; padding: 5px 15px;
            }
            QPushButton:hover { background-color: #7B1FA2; }
        """)
        btn_settings.clicked.connect(self.open_settings)
        settings_layout.addWidget(btn_settings)
        settings_layout.addStretch()
        layout.addLayout(settings_layout)

        progress_group = QGroupBox("Workflow Progress")
        progress_layout = QVBoxLayout(progress_group)
        self.progress_label = QLabel(f"Progress: 0/{len(self.step_names)} steps finished")
        self.progress_label.setFont(QFont("Arial", 10, QFont.Bold))
        progress_layout.addWidget(self.progress_label)
        layout.addWidget(progress_group)

        steps_group = QGroupBox("Processing Steps")
        steps_layout = QVBoxLayout(steps_group)
        for i, step_name in enumerate(self.step_names):
            btn = StepButton(i + 1, step_name)
            btn.clicked.connect(lambda checked, idx=i: self.open_step(idx))
            self.step_buttons.append(btn)
            steps_layout.addWidget(btn)
        layout.addWidget(steps_group)

        action_layout = QHBoxLayout()
        btn_resume = QPushButton("Resume Next Step")
        btn_resume.setFont(QFont("Arial", 11, QFont.Bold))
        btn_resume.setMinimumHeight(40)
        btn_resume.clicked.connect(self.resume_next_step)
        action_layout.addWidget(btn_resume)
        btn_reset = QPushButton("Reset Progress")
        btn_reset.setMinimumHeight(40)
        btn_reset.clicked.connect(self.reset_progress)
        action_layout.addWidget(btn_reset)
        layout.addLayout(action_layout)

        log_group = QGroupBox("Activity Log")
        log_layout = QVBoxLayout(log_group)
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMaximumHeight(120)
        self.log_text.setFont(QFont("Courier", 8))
        log_layout.addWidget(self.log_text)
        layout.addWidget(log_group)

        self.status_bar = self.statusBar()
        self.status_bar.showMessage("Ready")

    # ── Menu setup ───────────────────────────────────────────────────────────

    def setup_menu(self):
        menubar = self.menuBar()

        file_menu = menubar.addMenu("&File")
        action_save = QAction("&Save Project State", self)
        action_save.setShortcut("Ctrl+S")
        action_save.triggered.connect(self.save_project_state)
        file_menu.addAction(action_save)
        action_export = QAction("&Export Summary...", self)
        action_export.triggered.connect(self.export_summary)
        file_menu.addAction(action_export)
        file_menu.addSeparator()
        action_exit = QAction("E&xit", self)
        action_exit.setShortcut("Ctrl+Q")
        action_exit.triggered.connect(self.close)
        file_menu.addAction(action_exit)

        tools_menu = menubar.addMenu("&Tools")

        action_qa = QAction("QA Reports", self)
        action_qa.setShortcut("Ctrl+R")
        action_qa.triggered.connect(self.open_qa_report)
        tools_menu.addAction(action_qa)

        action_iraf = QAction("IRAF/DAOPHOT Tool", self)
        action_iraf.setShortcut("Ctrl+I")
        action_iraf.triggered.connect(self.open_iraf_tool)
        tools_menu.addAction(action_iraf)

        if self.mode == "lc":
            action_ext = QAction("Extinction Fit Tool", self)
            action_ext.setShortcut("Ctrl+E")
            action_ext.triggered.connect(self.open_extinction_tool)
            tools_menu.addAction(action_ext)

            action_airmass = QAction("Airmass Header Debug", self)
            action_airmass.triggered.connect(self.open_airmass_debug_tool)
            tools_menu.addAction(action_airmass)

            tools_menu.addSeparator()

            action_merger = QAction("Multi-Night Light Curve Merger", self)
            action_merger.setShortcut("Ctrl+M")
            action_merger.triggered.connect(self.open_multi_night_merger)
            tools_menu.addAction(action_merger)

            tools_menu.addSeparator()

            action_varstar = QAction("Variable Star Analysis", self)
            action_varstar.setShortcut("Ctrl+Shift+V")
            action_varstar.triggered.connect(self.open_variable_star_tool)
            tools_menu.addAction(action_varstar)

            action_transit = QAction("Exoplanet Transit Analysis", self)
            action_transit.setShortcut("Ctrl+Shift+T")
            action_transit.triggered.connect(self.open_transit_tool)
            tools_menu.addAction(action_transit)

            action_eb = QAction("Eclipsing Binary Analysis", self)
            action_eb.setShortcut("Ctrl+Shift+B")
            action_eb.triggered.connect(self.open_eb_tool)
            tools_menu.addAction(action_eb)

        elif self.mode == "cmd":
            action_ext_fit = QAction("Extinction (Airmass Fit)", self)
            action_ext_fit.triggered.connect(self.open_extinction_fit_cmd)
            tools_menu.addAction(action_ext_fit)

            action_cmd_prev = QAction("CMD + Isochrone (From Results)...", self)
            action_cmd_prev.triggered.connect(self.open_cmd_iso_tool)
            tools_menu.addAction(action_cmd_prev)

            action_gaia_3d = QAction("Gaia 3D Cluster Viewer", self)
            action_gaia_3d.triggered.connect(self.open_gaia_3d_viewer)
            tools_menu.addAction(action_gaia_3d)

            action_cluster = QAction("Analyze Cluster Structure", self)
            action_cluster.triggered.connect(self.open_cluster_structure_tool)
            tools_menu.addAction(action_cluster)

        help_menu = menubar.addMenu("&Help")
        action_about = QAction("&About", self)
        action_about.triggered.connect(self.show_about)
        help_menu.addAction(action_about)

    # ── Step button state ────────────────────────────────────────────────────

    def update_step_buttons(self):
        completed_count = len(self.project_state.state["completed_steps"])
        self.progress_label.setText(
            f"Progress: {completed_count}/{len(self.step_names)} steps finished"
        )
        for i, btn in enumerate(self.step_buttons):
            completed = self.project_state.is_step_completed(i)
            accessible = self.project_state.is_step_accessible(i)
            btn.set_completed(completed)
            btn.set_accessible(accessible)
            btn.setEnabled(accessible)

    # ── Step dispatch ────────────────────────────────────────────────────────

    def open_step(self, step_index: int):
        if not self.project_state.is_step_accessible(step_index):
            prev_idx = step_index - 1
            prev_name = (self.step_names[prev_idx]
                         if 0 <= prev_idx < len(self.step_names) else "previous step")
            QMessageBox.warning(self, "Step Not Accessible",
                                f"Please finish Step {step_index}: {prev_name} first.")
            return

        self.project_state.set_current_step(step_index)
        self.append_log(f"Opening Step {step_index + 1}: {self.step_names[step_index]}")

        if self.current_step_window:
            self.current_step_window.close()

        win = self._open_step_window(step_index)
        if win is None:
            QMessageBox.information(self, "Step Not Implemented",
                                    f"Step {step_index + 1} is not yet implemented.")
            return

        self.current_step_window = win
        win.show()

    def _open_step_window(self, step_index: int):  # noqa: C901 (complexity ok)
        p, fm, ps = self.params, self.file_manager, self.project_state

        # ── Step 0: File selection (mode-specific) ──
        if step_index == 0:
            if self.mode == "cmd":
                from apex.gui.workflow.cmd.step1_file_selection import FileSelectionWindow
            else:
                from apex.gui.workflow.lc.step1_file_selection import FileSelectionWindow
            return FileSelectionWindow(p, fm, ps, self)

        # ── Steps 1-4: shared ──
        elif step_index == 1:
            from apex.gui.workflow.step2_crop_selector import CropSelectorWindow
            return CropSelectorWindow(p, fm, ps, self)
        elif step_index == 2:
            from apex.gui.workflow.step3_sky_preview import SkyPreviewWindow
            return SkyPreviewWindow(p, fm, ps, self)
        elif step_index == 3:
            from apex.gui.workflow.step4_source_detection import SourceDetectionWindow
            return SourceDetectionWindow(p, fm, ps, self)
        elif step_index == 4:
            from apex.gui.workflow.step5_aperture_photometry import AperturePhotometryWindow
            return AperturePhotometryWindow(p, fm, ps, self)

        # ── Step 5: PSF (CMD only) / WCS (LC) ──
        elif step_index == 5:
            if self.mode == "cmd":
                from apex.gui.workflow.cmd.step6_psf_photometry import PSFPhotometryWindow
                return PSFPhotometryWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.step6_wcs_plate_solving import WcsPlateSolvingWindow
                return WcsPlateSolvingWindow(p, fm, ps, self)

        # ── Steps 6-7: WCS/Ref in CMD mode; Ref/IDMatch in LC ──
        elif step_index == 6:
            if self.mode == "cmd":
                from apex.gui.workflow.step6_wcs_plate_solving import WcsPlateSolvingWindow
                return WcsPlateSolvingWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.step7_ref_build import RefBuildWindow
                return RefBuildWindow(p, fm, ps, self)

        elif step_index == 7:
            if self.mode == "cmd":
                from apex.gui.workflow.step7_ref_build import RefBuildWindow
                return RefBuildWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.step8_star_id_matching import StarIdMatchingWindow
                return StarIdMatchingWindow(p, fm, ps, self)

        # ── Step 8: IDMatch (CMD) / Target selection (LC) ──
        elif step_index == 8:
            if self.mode == "cmd":
                from apex.gui.workflow.step8_star_id_matching import StarIdMatchingWindow
                return StarIdMatchingWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.lc.step9_target_selection import TargetComparisonSelectionWindow
                return TargetComparisonSelectionWindow(p, fm, ps, self)

        # ── CMD steps 9-12 ──
        elif step_index == 9:
            if self.mode == "cmd":
                from apex.gui.workflow.cmd.step10_master_id_editor import MasterIdEditorWindow
                return MasterIdEditorWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.lc.step10_lightcurve_builder import LightCurveBuilderWindow
                return LightCurveBuilderWindow(p, fm, ps, self)

        elif step_index == 10:
            if self.mode == "cmd":
                from apex.gui.workflow.cmd.step11_zeropoint_calibration import ZeropointCalibrationWindow
                return ZeropointCalibrationWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.lc.step11_detrend_merge import DetrendNightMergeWindow
                return DetrendNightMergeWindow(p, fm, ps, self)

        elif step_index == 11:
            if self.mode == "cmd":
                from apex.gui.workflow.cmd.step12_cmd_plot import CmdPlotWindow
                return CmdPlotWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.lc.step12_period_analysis import PeriodAnalysisWindow
                return PeriodAnalysisWindow(p, fm, ps, self)

        elif step_index == 12 and self.mode == "cmd":
            from apex.gui.workflow.cmd.step13_isochrone_model import IsochroneModelWindow
            return IsochroneModelWindow(p, fm, ps, self)

        return None

    def on_step_completed(self, step_index: int):
        self.project_state.mark_step_completed(step_index)
        self.update_step_buttons()
        self.append_log(f"\u2713 Step {step_index + 1} finished: {self.step_names[step_index]}")

    def resume_next_step(self):
        next_step = self.project_state.get_next_incomplete_step()
        if next_step is not None:
            if not self.project_state.is_step_accessible(next_step):
                prev_idx = next_step - 1
                prev_name = (self.step_names[prev_idx]
                             if 0 <= prev_idx < len(self.step_names) else "previous step")
                QMessageBox.warning(self, "Step Not Accessible",
                                    f"Please finish earlier steps first.\n"
                                    f"Required now: Step {next_step}: {prev_name}")
                return
            self.open_step(next_step)
        else:
            QMessageBox.information(self, "Workflow Finished",
                                    "All workflow steps are finished.")

    def reset_progress(self):
        reply = QMessageBox.question(
            self, "Reset Progress",
            "Are you sure you want to reset all progress?\n"
            "This will clear completion status but keep your data files.",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No
        )
        if reply == QMessageBox.Yes:
            self.project_state.reset()
            self.update_step_buttons()
            self.append_log("Progress reset")

    def save_project_state(self):
        self.project_state.save()
        self.append_log("Project state saved")
        QMessageBox.information(self, "Saved", "Project state saved successfully.")

    def export_summary(self):
        summary = self.project_state.export_summary()
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export Summary",
            str(self.params.P.result_dir / "project_summary.txt"),
            "Text Files (*.txt)"
        )
        if file_path:
            Path(file_path).write_text(summary, encoding="utf-8")
            self.append_log(f"Summary exported to {file_path}")

    def append_log(self, message: str):
        from datetime import datetime
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.log_text.append(f"[{timestamp}] {message}")

    def show_about(self):
        mode_label = "CMD Cluster Photometry" if self.mode == "cmd" else "LC Light Curve Analysis"
        QMessageBox.about(
            self, "About APEX",
            f"<h2>APEX — Automated Photometry EXtraction</h2>"
            f"<p><b>Mode: {mode_label}</b></p>"
            "<p>KNUEMAO Observatory — CDK500 + Moravian C3-61000</p>"
            "<p>Version 2.0.0</p>"
        )

    # ── Tool launchers ───────────────────────────────────────────────────────

    def open_qa_report(self, tab: int = 0):
        from apex.gui.tools.qa_report import QAReportWindow
        self.qa_window = QAReportWindow(self.params, self.params.P.result_dir, parent=None)
        self.qa_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        if hasattr(self.qa_window, "tabs") and tab >= 0:
            self.qa_window.tabs.setCurrentIndex(
                min(tab, self.qa_window.tabs.count() - 1))
        self.qa_window.show()
        self.qa_window.raise_()
        self.qa_window.activateWindow()
        self.append_log("Opened QA Report window")

    def open_iraf_tool(self):
        from apex.gui.tools.iraf_photometry import IRAFPhotometryWindow
        self.iraf_window = IRAFPhotometryWindow(
            self.params, self.params.P.data_dir, self.params.P.result_dir,
            self.project_state, parent=None
        )
        self.iraf_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.iraf_window.show()
        self.iraf_window.raise_()
        self.iraf_window.activateWindow()
        self.append_log("Opened IRAF/DAOPHOT Tool")

    # ── LC tools ─────────────────────────────────────────────────────────────

    def open_extinction_tool(self):
        from apex.gui.tools.extinction_fit import ExtinctionFitWindow
        self.extinction_window = ExtinctionFitWindow(
            self.params, self.params.P.data_dir, self.params.P.result_dir, parent=None
        )
        self.extinction_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.extinction_window.show()
        self.extinction_window.raise_()
        self.extinction_window.activateWindow()
        self.append_log("Opened Extinction (Airmass Fit) Tool")

    def open_airmass_debug_tool(self):
        from apex.gui.tools.airmass_debug import AirmassHeaderDebugToolWindow
        self.airmass_debug_window = AirmassHeaderDebugToolWindow(
            self.params, self.project_state, parent=None, file_manager=self.file_manager
        )
        self.airmass_debug_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.airmass_debug_window.show()
        self.airmass_debug_window.raise_()
        self.airmass_debug_window.activateWindow()
        self.append_log("Opened Airmass Header Debug Tool")

    def open_multi_night_merger(self):
        from apex.gui.tools.multi_night_merger import MultiNightMergerWindow
        self.merger_window = MultiNightMergerWindow(
            self.params, self.project_state, main_window=self
        )
        self.merger_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.hide()
        self.merger_window.show()
        self.merger_window.raise_()
        self.merger_window.activateWindow()
        self.append_log("Opened Multi-Night Light Curve Merger")

    def open_variable_star_tool(self):
        from apex.gui.tools.variable_star import VariableStarToolWindow
        if self.varstar_window is None:
            self.varstar_window = VariableStarToolWindow(
                self.params, self.project_state, parent=None
            )
            self.varstar_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.varstar_window.show()
        self.varstar_window.raise_()
        self.varstar_window.activateWindow()
        self.append_log("Opened Variable Star Analysis Tool")

    def open_transit_tool(self):
        from apex.gui.tools.transit_tool import TransitToolWindow
        self.transit_window = TransitToolWindow(
            self.params, self.project_state, parent=None
        )
        self.transit_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.transit_window.show()
        self.transit_window.raise_()
        self.transit_window.activateWindow()
        self.append_log("Opened Exoplanet Transit Analysis Tool")

    def open_eb_tool(self):
        from apex.gui.tools.eb_tool import EclipsingBinaryToolWindow
        self.eb_window = EclipsingBinaryToolWindow(
            self.params, self.project_state, parent=None
        )
        self.eb_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.eb_window.show()
        self.eb_window.raise_()
        self.eb_window.activateWindow()
        self.append_log("Opened Eclipsing Binary Analysis Tool")

    # ── CMD tools ────────────────────────────────────────────────────────────

    def open_extinction_fit_cmd(self):
        from apex.gui.tools.extinction_fit import ExtinctionFitWindow
        self.ext_window = ExtinctionFitWindow(
            self.params, self.params.P.data_dir, self.params.P.result_dir, parent=None
        )
        self.ext_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.ext_window.show()
        self.ext_window.raise_()
        self.ext_window.activateWindow()
        self.append_log("Opened Extinction Fit window")

    def open_cmd_iso_tool(self):
        start_dir = str(getattr(self.params.P, "result_dir", Path.cwd()))
        selected = QFileDialog.getExistingDirectory(self, "Select Result Folder", start_dir)
        if not selected:
            return
        from apex.gui.tools.cmd_iso_tool import CmdIsoToolWindow
        self.cmd_tool_window = CmdIsoToolWindow(
            self.params, self.file_manager, self.project_state, self,
            initial_result_dir=Path(selected), parent=None
        )
        self.cmd_tool_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.cmd_tool_window.show()
        self.cmd_tool_window.raise_()
        self.cmd_tool_window.activateWindow()
        self.append_log(f"Opened CMD + Isochrone tool: {selected}")

    def open_gaia_3d_viewer(self):
        from apex.gui.tools.gaia_3d_viewer import Gaia3DViewerWindow
        self.gaia_3d_window = Gaia3DViewerWindow(
            self.params, self.params.P.result_dir, parent=None
        )
        self.gaia_3d_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.gaia_3d_window.show()
        self.gaia_3d_window.raise_()
        self.gaia_3d_window.activateWindow()
        self.append_log(f"Opened Gaia 3D Viewer: {self.params.P.result_dir}")

    def open_cluster_structure_tool(self):
        from apex.gui.tools.cluster_structure import ClusterStructureWindow
        self.cluster_structure_window = ClusterStructureWindow(
            self.params, self.params.P.result_dir, parent=None
        )
        self.cluster_structure_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.cluster_structure_window.show()
        self.cluster_structure_window.raise_()
        self.cluster_structure_window.activateWindow()
        self.append_log(f"Opened Cluster Structure Tool: {self.params.P.result_dir}")

    # ── Instrument settings ──────────────────────────────────────────────────

    def open_settings(self):
        from PyQt5.QtWidgets import (
            QDialog, QVBoxLayout, QLabel, QLineEdit,
            QDialogButtonBox, QGroupBox, QFormLayout
        )
        dialog = QDialog(self)
        dialog.setWindowTitle("Instrument Settings")
        dialog.setMinimumWidth(650)
        layout = QVBoxLayout(dialog)
        title = QLabel("Instrument Configuration")
        title.setFont(QFont("Arial", 12, QFont.Bold))
        layout.addWidget(title)

        tel_group = QGroupBox("Telescope")
        tel_layout = QFormLayout(tel_group)
        tel_name_edit = QLineEdit(self.instrument.telescope_name)
        tel_layout.addRow("Name:", tel_name_edit)
        tel_aperture_edit = QLineEdit(str(self.instrument.aperture_mm))
        tel_aperture_edit.setMaximumWidth(150)
        tel_layout.addRow("Aperture (mm):", tel_aperture_edit)
        tel_focal_edit = QLineEdit(str(self.instrument.focal_length_mm))
        tel_focal_edit.setMaximumWidth(150)
        tel_layout.addRow("Focal Length (mm):", tel_focal_edit)
        layout.addWidget(tel_group)

        cam_group = QGroupBox("Camera")
        cam_layout = QFormLayout(cam_group)
        cam_name_edit = QLineEdit(self.instrument.camera_name)
        cam_layout.addRow("Name:", cam_name_edit)
        cam_pixsize_edit = QLineEdit(str(self.instrument.pix_size_um))
        cam_pixsize_edit.setMaximumWidth(150)
        cam_layout.addRow("Pixel Size (\u03bcm):", cam_pixsize_edit)
        cam_nx_edit = QLineEdit(str(self.instrument.sensor_nx_1x))
        cam_nx_edit.setMaximumWidth(150)
        cam_layout.addRow("Sensor Width (px):", cam_nx_edit)
        cam_ny_edit = QLineEdit(str(self.instrument.sensor_ny_1x))
        cam_ny_edit.setMaximumWidth(150)
        cam_layout.addRow("Sensor Height (px):", cam_ny_edit)
        cam_binning_edit = QLineEdit(str(self.instrument.binning))
        cam_binning_edit.setMaximumWidth(150)
        cam_layout.addRow("Binning:", cam_binning_edit)
        layout.addWidget(cam_group)

        params_group = QGroupBox("Camera Parameters")
        params_layout = QFormLayout(params_group)
        gain_edit = QLineEdit(str(self.params.P.gain_e_per_adu))
        gain_edit.setMaximumWidth(150)
        params_layout.addRow("Gain (e-/ADU):", gain_edit)
        rdnoise_edit = QLineEdit(str(self.params.P.rdnoise_e))
        rdnoise_edit.setMaximumWidth(150)
        params_layout.addRow("Read Noise (e-):", rdnoise_edit)
        saturation_edit = QLineEdit(str(self.params.P.saturation_adu))
        saturation_edit.setMaximumWidth(150)
        params_layout.addRow("Saturation (ADU):", saturation_edit)
        layout.addWidget(params_group)

        parallel_group = QGroupBox("Parallel Processing")
        parallel_layout = QFormLayout(parallel_group)
        parallel_workers_spin = QSpinBox()
        parallel_workers_spin.setRange(0, 16)
        parallel_workers_spin.setValue(int(getattr(self.params.P, "max_workers", 0)))
        parallel_workers_spin.setToolTip("0 = auto (use ~75% of CPU cores)")
        parallel_layout.addRow("Max Workers (0=auto):", parallel_workers_spin)
        layout.addWidget(parallel_group)

        site_group = QGroupBox("Observatory Location")
        site_layout = QFormLayout(site_group)
        site_lat_edit = QLineEdit(str(getattr(self.params.P, "site_lat_deg", 0.0)))
        site_lat_edit.setMaximumWidth(150)
        site_layout.addRow("Latitude (deg):", site_lat_edit)
        site_lon_edit = QLineEdit(str(getattr(self.params.P, "site_lon_deg", 0.0)))
        site_lon_edit.setMaximumWidth(150)
        site_layout.addRow("Longitude (deg):", site_lon_edit)
        site_alt_edit = QLineEdit(str(getattr(self.params.P, "site_alt_m", 0.0)))
        site_alt_edit.setMaximumWidth(150)
        site_layout.addRow("Altitude (m):", site_alt_edit)
        site_tz_edit = QLineEdit(str(getattr(self.params.P, "site_tz_offset_hours", 0.0)))
        site_tz_edit.setMaximumWidth(150)
        site_layout.addRow("UTC Offset (hours):", site_tz_edit)
        layout.addWidget(site_group)

        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)

        if dialog.exec_() == QDialog.Accepted:
            try:
                self.instrument.telescope_name = tel_name_edit.text().strip()
                self.instrument.aperture_mm = float(tel_aperture_edit.text())
                self.instrument.focal_length_mm = float(tel_focal_edit.text())
                self.instrument.focal_ratio = self.instrument.focal_length_mm / self.instrument.aperture_mm
                self.instrument.camera_name = cam_name_edit.text().strip()
                self.instrument.pix_size_um = float(cam_pixsize_edit.text())
                self.instrument.sensor_nx_1x = int(cam_nx_edit.text())
                self.instrument.sensor_ny_1x = int(cam_ny_edit.text())
                self.instrument.binning = int(cam_binning_edit.text())
                self.params.P.gain_e_per_adu = float(gain_edit.text())
                self.params.P.rdnoise_e = float(rdnoise_edit.text())
                self.params.P.saturation_adu = float(saturation_edit.text())
                self.params.P.site_lat_deg = float(site_lat_edit.text())
                self.params.P.site_lon_deg = float(site_lon_edit.text())
                self.params.P.site_alt_m = float(site_alt_edit.text())
                self.params.P.site_tz_offset_hours = float(site_tz_edit.text())
                self.params.P.max_workers = int(parallel_workers_spin.value())
                if not self.params.save_toml():
                    self.append_log("Warning: could not save settings to TOML")
                QMessageBox.information(self, "Settings Saved",
                                        "Instrument settings updated.")
                self.append_log("Instrument settings updated")
            except ValueError as e:
                QMessageBox.warning(self, "Invalid Input",
                                    f"Please enter valid numeric values.\n{e}")

    # ── Window close ─────────────────────────────────────────────────────────

    def closeEvent(self, event):
        self.project_state.save()
        if self.current_step_window:
            self.current_step_window.close()
        event.accept()
