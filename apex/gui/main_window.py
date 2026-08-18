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
from PyQt5.QtCore import Qt, pyqtSignal, QObject, QEvent, QEventLoop
from PyQt5.QtGui import QFont, QIcon, QPixmap, QPainter
import shutil
from pathlib import Path
from typing import Optional, List

from apex.gui.layout_rules import AutoFitMixin, FittedDialog
from apex.gui.theme import Tokens, style_button
from apex.gui.tools.registry import iter_tools_for_mode
from apex.config.config_io import (
    check_workspace_identity, load_config_data, save_config_data,
)
from apex.config.recent_workspaces import (
    display_label as recent_label,
    forget as forget_workspace,
    load_recent as load_recent_workspaces,
    remember as remember_workspace,
)

_RESOURCES = Path(__file__).resolve().parent.parent / "resources"


def _migrate_lc_optional_psf_state(project_state) -> bool:
    """Insert the optional PSF step without losing legacy LC progress."""
    step_data = project_state.state.setdefault("step_data", {})
    migrations = step_data.setdefault("workflow_migrations", {})
    if migrations.get("lc_optional_psf_v1"):
        return False
    completed = [int(value) for value in project_state.state.get("completed_steps", [])]
    current = int(project_state.state.get("current_step", 0))
    legacy_downstream = current >= 7 or any(value >= 7 for value in completed)
    if legacy_downstream:
        shifted = {value if value < 7 else value + 1 for value in completed}
        shifted.add(7)
        project_state.state["completed_steps"] = sorted(shifted)
        project_state.state["current_step"] = current if current < 7 else current + 1
        psf_state = step_data.setdefault("psf_photometry", {})
        psf_state.setdefault("skip_psf", True)
        psf_state.setdefault("migration_reason", "legacy LC session")
    migrations["lc_optional_psf_v1"] = True
    project_state.save()
    return legacy_downstream


def _svg_to_pixmap(svg_path: Path, size: int = 256) -> QPixmap:
    # Rasterize at the display's device-pixel ratio so the vector logo stays
    # crisp under HiDPI scaling instead of being upscaled from 1×.
    app = QApplication.instance()
    screen = app.primaryScreen() if app is not None else None
    dpr = float(screen.devicePixelRatio()) if screen is not None else 1.0
    side = max(1, int(round(size * dpr)))
    pixmap = QPixmap(side, side)
    pixmap.fill(Qt.transparent)
    try:
        from PyQt5.QtSvg import QSvgRenderer
        renderer = QSvgRenderer(str(svg_path))
        painter = QPainter(pixmap)
        renderer.render(painter)
        painter.end()
    except Exception:
        pass
    pixmap.setDevicePixelRatio(dpr)
    return pixmap


def _apply_opacity(pixmap: QPixmap, opacity: float) -> QPixmap:
    result = QPixmap(pixmap.size())
    result.fill(Qt.transparent)
    result.setDevicePixelRatio(pixmap.devicePixelRatio())  # preserve HiDPI scale
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
        self.emphasized = False  # the single "do this next" step
        self.setText(f"Step {step_number}: {step_name}")
        # Fixed-vertical size policy so a cramped window cannot crush the
        # 12 step buttons into illegible 25px slivers on small displays.
        self.setMinimumHeight(46)
        self.setMinimumWidth(300)
        from PyQt5.QtWidgets import QSizePolicy
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.update_appearance()

    def set_completed(self, completed: bool):
        self.completed = completed
        self.update_appearance()

    def set_accessible(self, accessible: bool):
        self.accessible = accessible
        self.update_appearance()

    def set_emphasized(self, emphasized: bool):
        if emphasized != self.emphasized:
            self.emphasized = emphasized
            self.update_appearance()

    def update_appearance(self):
        # Quiet status rows instead of a wall of saturated green/blue fills:
        # every row is a neutral surface; state is carried by a 3px left bar
        # plus the leading glyph. Only the *next actionable* step gets an
        # accent-tinted background so the eye lands on "what to do next".
        # Styles are rebuilt from Tokens on every call, so a theme switch
        # restyles the list via update_step_buttons().
        t = Tokens
        base = (
            f"font-size: 14px; border-radius: {t.RADIUS_SM}px;"
            f" text-align: left; padding: 10px 14px;"
            f" border: 1px solid {t.BORDER};"
        )
        if self.completed:
            self.setStyleSheet(
                f"QPushButton {{ {base} background: {t.SURFACE}; color: {t.TEXT};"
                f" border-left: 3px solid {t.OK}; }}"
                f"QPushButton:hover {{ background: {t.SURFACE_ALT}; }}"
            )
            self.setText(f"\u2713 Step {self.step_number}: {self.step_name}")
            self.setEnabled(True)
        elif self.accessible and self.emphasized:
            self.setStyleSheet(
                f"QPushButton {{ {base} background: {t.ACCENT_SOFT}; color: {t.TEXT};"
                f" border: 1px solid {t.ACCENT}; border-left: 3px solid {t.ACCENT};"
                f" font-weight: 600; }}"
                f"QPushButton:hover {{ border-color: {t.ACCENT_HOVER}; }}"
            )
            self.setText(f"\u25cb Step {self.step_number}: {self.step_name}")
            self.setEnabled(True)
        elif self.accessible:
            # Reachable but not the next action \u2014 quiet row, accent bar only.
            self.setStyleSheet(
                f"QPushButton {{ {base} background: {t.SURFACE}; color: {t.TEXT};"
                f" border-left: 3px solid {t.ACCENT}; }}"
                f"QPushButton:hover {{ background: {t.SURFACE_ALT}; }}"
            )
            self.setText(f"\u25cb Step {self.step_number}: {self.step_name}")
            self.setEnabled(True)
        else:
            # No glyph: the \ud83d\udd12 emoji stays multicolour on Windows even with
            # U+FE0E, and a muted row reads as locked on its own.
            self.setStyleSheet(
                f"QPushButton {{ {base} background: transparent; color: {t.TEXT_MUTED};"
                f" border-left: 3px solid {t.BORDER}; }}"
                f"QPushButton:disabled {{ color: {t.TEXT_MUTED}; }}"
            )
            self.setText(f"Step {self.step_number}: {self.step_name}")
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


class MainWindowWorkflow(AutoFitMixin, QMainWindow):
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
            param_path = Path(param_file) if param_file else Path("apex_config.json")
            self.params = Parameters(param_path)
            # 파라미터 파일을 명시해서 연 경우에는 그 파일이 기준이다.
            # 기본 parameters.toml 로 연 평소 실행에서는 예전처럼 마지막으로
            # 고른 폴더를 복원한다(_bootstrap_file_selection_state).
            self._explicit_param_file = bool(param_file)

            from apex.core import InstrumentConfig, FileManager, ProjectState
            self.instrument = InstrumentConfig(self.params)
            self.file_manager = FileManager(self.params)

            self.project_state = self._make_project_state(migrate_legacy=True)
            self._register_state_mirror()

            self._bootstrap_file_selection_state(respect_explicit_param=True)
            # Without this the workspace opened at launch is missing from
            # Recent, so switching away from a `--params` workspace leaves no
            # way back to it.
            try:
                remember_workspace(self.params.param_file)
            except Exception:
                pass

        except Exception as e:
            QMessageBox.critical(self, "Initialization Error",
                                 f"Failed to load parameters:\n{e}")
            raise

        if mode == "cmd":
            self.step_names = [
                "File Selection",           # 0
                "Image Crop",               # 1
                "Sky Preview & QC",         # 2
                "Source Detection",         # 3
                "WCS Plate Solving",        # 4
                "Master Catalog Build",     # 5
                "Forced Aperture Phot",     # 6
                "PSF Photometry",           # 7  (skippable)
                "Master ID Editor",         # 8
                "Zeropoint Calibration",    # 9
                "CMD Plot",                 # 10
                "Isochrone Model",          # 11
            ]
        else:
            self.step_names = [
                "File Selection",           # 0
                "Image Crop",               # 1
                "Sky Preview & QC",         # 2
                "Source Detection",         # 3
                "WCS Plate Solving",        # 4
                "Master Catalog Build",     # 5
                "Forced Aperture Phot",     # 6
                "PSF Photometry",           # 7  (optional)
                "Target/Comparison Selection",  # 8
                "Light Curve Builder",      # 9
                "Detrend & Night Merge",    # 10
                "Period Analysis",          # 11
            ]

        if mode == "lc":
            _migrate_lc_optional_psf_state(self.project_state)
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

        if hasattr(self, "_offline_data_dir"):
            QMessageBox.warning(
                self, "Previous Data Path Unavailable",
                f"The last-used data path is inaccessible:\n\n"
                f"  {self._offline_data_dir}\n\n"
                "External drive may be disconnected.\n"
                "Please set a data path in Step 1."
            )

    # ── LC file-selection bootstrap ──────────────────────────────────────────

    def _bootstrap_file_selection_state(self, respect_explicit_param: bool = False) -> None:
        """저장된 file_selection 으로 params 경로를 되살린다.

        `respect_explicit_param` 은 창을 처음 만들 때만 켠다. 프로젝트를
        불러오는 경로(_load_project_state)는 일부러 params 를 그 프로젝트로
        옮기는 것이므로 막으면 안 된다.
        """
        state_data = self.project_state.get_step_data("file_selection")
        if not state_data:
            return
        data_dir = state_data.get("data_dir")
        if data_dir and respect_explicit_param and getattr(self, "_explicit_param_file", False):
            # 명시한 파라미터 파일과 다른 프로젝트의 상태면 덮어쓰지 않는다.
            # 이 상태는 모드당 하나(apex/.state/<mode>)라 마지막에 GUI 로 연
            # 프로젝트가 남아 있고, 그대로 두면 지정한 result_dir 이 조용히
            # 바뀌어 엉뚱한 자료를 읽는다.
            if Path(data_dir) != Path(self.params.P.data_dir):
                data_dir = None
        if data_dir:
            self.params.P.data_dir = Path(data_dir)
            saved_result_dir = state_data.get("result_dir")
            if saved_result_dir:
                self.params.P.result_dir = Path(saved_result_dir)
            elif self.mode == "lc":
                # Legacy LC sessions defaulted result_dir to data_dir/result.
                self.params.P.result_dir = self.params.P.data_dir / "result"
            # CMD without saved result_dir keeps params.toml's value.
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
        self._register_state_mirror()

    def _make_project_state(self, migrate_legacy: bool = False):
        """Progress state belonging to the workspace, not to the application.

        It used to live in one shared ``apex/.state/<mode>/`` regardless of
        which target was open, which was invisible while changing target meant
        relaunching. With a workspace menu it is not: opening an untouched
        cluster showed "9/12 steps finished" inherited from the previous one.
        The per-workspace copy that ``_register_state_mirror`` has been writing
        into ``result_dir`` all along is the right home, so that copy becomes
        the state itself.

        ``migrate_legacy`` seeds a workspace that has no copy yet from the old
        shared file — at launch that file belongs to the last workspace used,
        which is the one being opened.
        """
        from apex.core import ProjectState

        result_dir = getattr(self.params.P, "result_dir", None)
        if not result_dir:
            legacy = Path(__file__).parent.parent / ".state" / self.mode
            legacy.mkdir(parents=True, exist_ok=True)
            return ProjectState(legacy)

        state_dir = Path(result_dir)
        state_dir.mkdir(parents=True, exist_ok=True)
        if migrate_legacy and not (state_dir / "project_state.json").exists():
            legacy_file = (Path(__file__).parent.parent / ".state" / self.mode
                           / "project_state.json")
            if legacy_file.is_file():
                try:
                    shutil.copy2(legacy_file, state_dir / "project_state.json")
                except OSError:
                    pass
        return ProjectState(state_dir)

    def _register_state_mirror(self) -> None:
        """Mirror project_state.json into the current result_dir so previous
        sessions are discoverable by Load Previous Session."""
        result_dir = getattr(self.params.P, "result_dir", None)
        if not result_dir:
            return
        try:
            mirror = Path(result_dir) / "project_state.json"
            if mirror.resolve() == self.project_state.state_file.resolve():
                return
        except OSError:
            return
        self.project_state.add_mirror_path(mirror)

    # ── UI setup ─────────────────────────────────────────────────────────────

    def setup_ui(self):
        mode_title = "CMD Cluster Photometry" if self.mode == "cmd" else "Light Curve Analysis"
        # Same title format the switch produces, so a workspace opened at
        # launch is labelled like one opened from the menu.
        self.setWindowTitle(
            self._workspace_title(getattr(self.params, "param_file", "")))
        icon = _load_icon(self.mode)
        self.setWindowIcon(icon)
        QApplication.instance().setWindowIcon(icon)
        # Portrait ~3:4 proportion (narrow width): the workflow is a vertical
        # list of step buttons, so a tall/narrow window suits it and leaves the
        # desktop uncluttered. Steps scroll inside the Processing Steps group on
        # short screens; the Activity Log is collapsed by default.
        self.setMinimumSize(520, 680)
        self.resize(600, 800)

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

        title = QLabel("APEX — Automated Photometry EXtraction")
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
        style_button(btn_settings, height=Tokens.H_ACTION)
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
        steps_outer = QVBoxLayout(steps_group)
        steps_outer.setContentsMargins(8, 12, 8, 8)
        # The button list goes inside a QScrollArea so cramped windows can
        # scroll instead of crushing 12 buttons into 25px slivers.
        from PyQt5.QtWidgets import QScrollArea, QFrame, QSizePolicy
        steps_scroll = QScrollArea()
        steps_scroll.setWidgetResizable(True)
        steps_scroll.setFrameShape(QFrame.NoFrame)
        steps_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        steps_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        steps_inner = QWidget()
        steps_layout = QVBoxLayout(steps_inner)
        steps_layout.setContentsMargins(0, 0, 0, 0)
        steps_layout.setSpacing(4)
        # Optional off-chain "Step 0": detector calibration. Rendered at the top
        # of the workflow but tracked separately (never in completed_steps), so
        # it shifts no numbered index. Always accessible / re-runnable, and kept
        # OUT of self.step_buttons (that list stays 1:1 with the linear chain).
        self.calib_button = StepButton(0, "Detector Calibration")
        self.calib_button.clicked.connect(lambda checked=False: self.open_calibration())
        steps_layout.addWidget(self.calib_button)
        for i, step_name in enumerate(self.step_names):
            btn = StepButton(i + 1, step_name)
            btn.clicked.connect(lambda checked, idx=i: self.open_step(idx))
            self.step_buttons.append(btn)
            steps_layout.addWidget(btn)
        steps_layout.addStretch(1)
        steps_scroll.setWidget(steps_inner)
        steps_outer.addWidget(steps_scroll)
        steps_group.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        layout.addWidget(steps_group, stretch=1)

        action_layout = QHBoxLayout()
        btn_load = QPushButton("Load Previous Session...")
        btn_load.setFont(QFont("Arial", 11, QFont.Bold))
        btn_load.setMinimumHeight(40)
        btn_load.clicked.connect(self.load_previous_session)
        action_layout.addWidget(btn_load)
        btn_reset = QPushButton("Reset Progress")
        btn_reset.setMinimumHeight(40)
        btn_reset.clicked.connect(self.reset_progress)
        action_layout.addWidget(btn_reset)
        layout.addLayout(action_layout)

        # Activity Log as a collapsible section: a flat header row with a
        # disclosure arrow, folded by default so it doesn't eat vertical space.
        from PyQt5.QtWidgets import QToolButton
        log_container = QWidget()
        log_v = QVBoxLayout(log_container)
        log_v.setContentsMargins(0, 0, 0, 0)
        log_v.setSpacing(4)
        self.log_toggle = QToolButton()
        self.log_toggle.setText("Activity Log")
        self.log_toggle.setCheckable(True)
        self.log_toggle.setChecked(False)
        self.log_toggle.setArrowType(Qt.RightArrow)
        self.log_toggle.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.log_toggle.setAutoRaise(True)
        self.log_toggle.setStyleSheet("QToolButton { font-weight: bold; border: none; }")
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMaximumHeight(110)
        self.log_text.setFont(QFont("Courier", 8))
        self.log_text.setVisible(False)

        def _toggle_log(checked):
            self.log_text.setVisible(checked)
            self.log_toggle.setArrowType(Qt.DownArrow if checked else Qt.RightArrow)

        self.log_toggle.toggled.connect(_toggle_log)
        log_v.addWidget(self.log_toggle)
        log_v.addWidget(self.log_text)
        layout.addWidget(log_container)

        self.status_bar = self.statusBar()
        self.status_bar.showMessage("Ready")

    # ── Workspace switching ──────────────────────────────────────────────────

    def _close_child_windows(self) -> None:
        """Close every step/tool window before the parameters change under them.

        Child windows receive ``params``, ``file_manager`` and ``project_state``
        when they open and keep their own references, so one left open across a
        switch would go on reading and writing the previous target's paths —
        the quiet kind of wrong. They are found by attribute name rather than
        from a hand-written list, so a tool added later is covered without
        anyone remembering to update this.
        """
        for name in [n for n in vars(self) if n.endswith("_window")]:
            window = getattr(self, name, None)
            if window is self or not isinstance(window, QWidget):
                continue
            try:
                window.close()
            except Exception:
                pass
            setattr(self, name, None)

    def _workspace_title(self, config_path) -> str:
        """Title that names the open workspace, not just the mode.

        The accident this helps with: two workspaces open on two days look
        identical in the title bar, and the one being edited is not the one the
        user thinks.
        """
        mode_title = ("CMD Cluster Photometry" if self.mode == "cmd"
                      else "Light Curve Analysis")
        target = ""
        try:
            data, _ = load_config_data(config_path)
            target = str((data.get("target") or {}).get("name") or "").strip()
        except Exception:
            pass
        workspace = Path(config_path).parent.name
        label = " · ".join(x for x in (target, workspace) if x)
        return f"APEX — {mode_title}" + (f" — {label}" if label else "")

    def switch_workspace(self, config_path) -> bool:
        """Load a different workspace without restarting the app."""
        config_path = Path(config_path)
        if not config_path.is_file():
            QMessageBox.warning(self, "워크스페이스",
                                f"파일이 없습니다:\n{config_path}")
            forget_workspace(config_path)
            return False

        if self.mode == "cmd":
            from apex.config.parameters_cmd import Parameters
        else:
            from apex.config.parameters_lc import Parameters
        from apex.core import InstrumentConfig, FileManager

        try:
            params = Parameters(config_path)
        except Exception as exc:
            QMessageBox.critical(self, "워크스페이스",
                                 f"파라미터를 읽지 못했습니다:\n{exc}")
            return False

        self._close_child_windows()
        self.current_step_window = None
        self.params = params
        self.instrument = InstrumentConfig(params)
        self.file_manager = FileManager(params)
        # Progress belongs to the workspace; without this the new target
        # inherits the previous one's completed steps.
        self.project_state = self._make_project_state()
        self._register_state_mirror()
        # The workspace names its own directories, so the remembered
        # file-selection folder from the previous target must not win.
        self._explicit_param_file = True
        try:
            self._bootstrap_file_selection_state(respect_explicit_param=True)
        except Exception as exc:
            self.append_log(f"[WARN] file_selection 복원 실패: {exc}")

        self.setWindowTitle(self._workspace_title(config_path))
        try:
            self.update_step_buttons()
        except Exception as exc:
            self.append_log(f"[WARN] 스텝 목록 갱신 실패: {exc}")
        remember_workspace(config_path)
        self.append_log(f"워크스페이스 전환: {config_path}")

        # The historical accident this guards: io paths pointing at one cluster
        # while [target] name still says another. Warn-only by design.
        try:
            data, resolved = load_config_data(config_path)
            issues = check_workspace_identity(data, resolved)
        except Exception:
            issues = []
        if issues:
            QMessageBox.warning(
                self, "워크스페이스 확인",
                "이 워크스페이스의 설정이 서로 어긋나 보입니다:\n\n• "
                + "\n• ".join(issues)
                + "\n\n실행은 막지 않습니다. 설정을 확인하세요.")
        return True

    def open_workspace_dialog(self) -> None:
        start = str(Path(getattr(self.params, "param_file", ".")).parent)
        path, _ = QFileDialog.getOpenFileName(
            self, "워크스페이스 열기", start,
            "APEX 설정 (apex_config*.json);;JSON (*.json);;모든 파일 (*)")
        if path:
            self.switch_workspace(path)

    def new_workspace_dialog(self) -> None:
        """Create a workspace from the current one's settings.

        Instrument, detection and fitting settings are what a user tunes once
        and reuses; the paths are what changes per target. So a new workspace
        inherits the former and clears the latter, rather than starting from a
        template that has neither.
        """
        directory = QFileDialog.getExistingDirectory(
            self, "새 워크스페이스 폴더 선택",
            str(Path(getattr(self.params, "param_file", ".")).parent.parent))
        if not directory:
            return
        target_dir = Path(directory)
        config_path = target_dir / "apex_config.json"
        if config_path.exists():
            answer = QMessageBox.question(
                self, "새 워크스페이스",
                f"이미 설정이 있습니다:\n{config_path}\n\n그것을 열까요?",
                QMessageBox.Yes | QMessageBox.No)
            if answer == QMessageBox.Yes:
                self.switch_workspace(config_path)
            return
        try:
            data, _ = load_config_data(getattr(self.params, "param_file", ""))
            data = dict(data)
            io_block = dict(data.get("io") or {})
            io_block["data_dir"] = ""
            io_block["result_dir"] = str(target_dir / "result")
            data["io"] = io_block
            data["target"] = dict(data.get("target") or {}, name=target_dir.name)
            save_config_data(config_path, data)
        except Exception as exc:
            QMessageBox.critical(self, "새 워크스페이스",
                                 f"설정을 만들지 못했습니다:\n{exc}")
            return
        if self.switch_workspace(config_path):
            QMessageBox.information(
                self, "새 워크스페이스",
                f"만들었습니다:\n{config_path}\n\n"
                "Step 1 에서 데이터 폴더를 지정하세요.")

    def _rebuild_recent_menu(self) -> None:
        self._recent_menu.clear()
        entries = load_recent_workspaces()
        if not entries:
            action = QAction("(없음)", self)
            action.setEnabled(False)
            self._recent_menu.addAction(action)
            return
        current = str(Path(getattr(self.params, "param_file", "")).absolute()).lower()
        for path in entries:
            action = QAction(recent_label(path), self)
            if str(path).lower() == current:
                action.setEnabled(False)      # already open
            action.triggered.connect(
                lambda _checked=False, p=str(path): self.switch_workspace(p))
            self._recent_menu.addAction(action)

    # ── Menu setup ───────────────────────────────────────────────────────────

    def setup_menu(self):
        menubar = self.menuBar()

        file_menu = menubar.addMenu("&File")
        # Changing target used to mean editing a config by hand or relaunching
        # with --params; these three make it a menu action.
        action_open_ws = QAction("&Open Workspace...", self)
        action_open_ws.setShortcut("Ctrl+O")
        action_open_ws.triggered.connect(self.open_workspace_dialog)
        file_menu.addAction(action_open_ws)
        action_new_ws = QAction("&New Workspace...", self)
        action_new_ws.triggered.connect(self.new_workspace_dialog)
        file_menu.addAction(action_new_ws)
        self._recent_menu = file_menu.addMenu("Open &Recent")
        # Rebuilt on open so a workspace that vanished since launch is not
        # offered, and the one already open is greyed out.
        self._recent_menu.aboutToShow.connect(self._rebuild_recent_menu)
        file_menu.addSeparator()
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
        for spec in iter_tools_for_mode(self.mode):
            if spec.separator_before:
                tools_menu.addSeparator()
            action = QAction(spec.label, self)
            if spec.shortcut:
                action.setShortcut(spec.shortcut)
            action.triggered.connect(getattr(self, spec.launcher))
            tools_menu.addAction(action)

        # Settings — app-wide configuration lives here (not under File).
        settings_menu = menubar.addMenu("&Settings")
        action_instrument = QAction("&Instrument Settings...", self)
        action_instrument.triggered.connect(self.open_settings)
        settings_menu.addAction(action_instrument)
        settings_menu.addSeparator()

        # Theme presets — applied app-wide immediately and persisted to
        # ~/.apex/theme.txt (picked up by all entry points on next launch).
        from apex.gui.theme import STANDARD_THEMES, THEME_PRESETS, current_theme
        theme_menu = settings_menu.addMenu("&Theme")
        self._theme_actions = {}
        active = current_theme()
        separated = False
        for key, label in THEME_PRESETS:
            # Split APEX's own presets from the published third-party designs.
            if key in STANDARD_THEMES and not separated:
                theme_menu.addSeparator()
                separated = True
            act = QAction(label, self)
            act.setCheckable(True)
            act.setChecked(key == active)
            act.triggered.connect(lambda _=False, k=key: self._switch_theme(k))
            theme_menu.addAction(act)
            self._theme_actions[key] = act

        help_menu = menubar.addMenu("&Help")
        action_wcs_help = QAction("WCS Solver Installation Help...", self)
        action_wcs_help.triggered.connect(self.show_wcs_solver_help)
        action_about = QAction("&About", self)
        action_about.triggered.connect(self.show_about)
        help_menu.addAction(action_wcs_help)
        help_menu.addSeparator()
        help_menu.addAction(action_about)

    def _switch_theme(self, key: str):
        """Apply + persist a theme preset and refresh hand-styled widgets."""
        from apex.gui.theme import set_theme
        applied = set_theme(QApplication.instance(), key)
        for k, act in self._theme_actions.items():
            act.setChecked(k == applied)
        # Step buttons build their stylesheet from Tokens at set-time — rerun.
        self.update_step_buttons()
        self.append_log(f"Theme: {applied}")

    # ── Step button state ────────────────────────────────────────────────────

    def update_step_buttons(self):
        completed_count = len(self.project_state.state["completed_steps"])
        self.progress_label.setText(
            f"Progress: {completed_count}/{len(self.step_names)} steps finished"
        )
        next_emphasized = False  # only the FIRST reachable-incomplete step pops
        for i, btn in enumerate(self.step_buttons):
            completed = self.project_state.is_step_completed(i)
            accessible = self.project_state.is_step_accessible(i)
            btn.set_completed(completed)
            btn.set_accessible(accessible)
            emphasize = accessible and not completed and not next_emphasized
            if emphasize:
                next_emphasized = True
            btn.set_emphasized(emphasize)
            btn.setEnabled(accessible)

        # Off-chain calibration button: driven by its own status, always usable.
        calib_button = getattr(self, "calib_button", None)
        if calib_button is not None:
            status = self.project_state.calibration_status()
            calib_button.set_completed(status == "done")
            calib_button.set_accessible(True)
            calib_button.setEnabled(True)

    # ── Step dispatch ────────────────────────────────────────────────────────

    def open_step(self, step_index: int):
        if getattr(self, "_step_open_in_progress", False):
            return
        if not self.project_state.is_step_accessible(step_index):
            prev_idx = step_index - 1
            prev_name = (self.step_names[prev_idx]
                         if 0 <= prev_idx < len(self.step_names) else "previous step")
            QMessageBox.warning(self, "Step Not Accessible",
                                f"Please finish Step {step_index}: {prev_name} first.")
            return

        self.append_log(f"Opening Step {step_index + 1}: {self.step_names[step_index]}")
        previous_window = self.current_step_window
        status = self.statusBar()
        status.showMessage(f"Opening Step {step_index + 1}...")
        self._step_open_in_progress = True
        QApplication.setOverrideCursor(Qt.WaitCursor)
        QApplication.processEvents(QEventLoop.ExcludeUserInputEvents)
        try:
            win = self._open_step_window(step_index)
        except Exception as exc:
            import traceback

            detail = traceback.format_exc()
            self.append_log(
                f"[ERROR] Failed to open Step {step_index + 1}: {exc}\n{detail}"
            )
            QMessageBox.critical(
                self,
                "Step Load Failed",
                f"Could not open Step {step_index + 1}:\n{exc}",
            )
            return
        finally:
            QApplication.restoreOverrideCursor()
            status.clearMessage()
            self._step_open_in_progress = False

        if win is None:
            QMessageBox.information(self, "Step Not Implemented",
                                    f"Step {step_index + 1} is not yet implemented.")
            return

        self.project_state.set_current_step(step_index)
        if previous_window is not None and previous_window is not win:
            previous_window.close()
        self.current_step_window = win
        win.show()

    def open_calibration(self):
        """Open the optional off-chain Step 0 (detector calibration) window.

        Always accessible; never sets current_step or touches completed_steps.
        The window records its outcome via ``project_state.mark_calibration``.
        """
        try:
            from apex.gui.workflow.step0_detector_calibration import (
                DetectorCalibrationWindow,
            )
        except Exception as exc:  # pragma: no cover - defensive
            QMessageBox.warning(self, "Detector Calibration",
                                f"Could not open the calibration window:\n{exc}")
            return

        self.append_log("Opening Step 0: Detector Calibration")
        if self.current_step_window:
            self.current_step_window.close()
        win = DetectorCalibrationWindow(self.params, self.project_state, self)
        self.current_step_window = win
        win.show()

    def _open_step_window(self, step_index: int):  # noqa: C901 (complexity ok)
        p, fm, ps = self.params, self.file_manager, self.project_state

        # ── Step index constants (0-based) ───────────────────────────────────
        # Shared steps (identical index in both modes)
        _FILE_SEL   = 0
        _CROP       = 1
        _SKY        = 2
        _DETECT     = 3
        _WCS        = 4
        _MASTER     = 5
        _FORCED     = 6
        # CMD-only
        _CMD_PSF    = 7
        _CMD_EDITOR = 8
        _CMD_ZP     = 9
        _CMD_PLOT   = 10
        _CMD_ISO    = 11
        # LC-only (after shared steps)
        _LC_PSF     = 7
        _LC_TARGET  = 8
        _LC_BUILDER = 9
        _LC_DETREND = 10
        _LC_PERIOD  = 11

        # ── Step 0: File selection (mode-specific) ──
        if step_index == _FILE_SEL:
            if self.mode == "cmd":
                from apex.gui.workflow.cmd.step1_file_selection import FileSelectionWindow
            else:
                from apex.gui.workflow.lc.step1_file_selection import FileSelectionWindow
            return FileSelectionWindow(p, fm, ps, self)

        # ── Steps 1-3: shared ──
        elif step_index == _CROP:
            from apex.gui.workflow.step2_crop_selector import CropSelectorWindow
            return CropSelectorWindow(p, fm, ps, self)
        elif step_index == _SKY:
            from apex.gui.workflow.step3_sky_preview import SkyPreviewWindow
            return SkyPreviewWindow(p, fm, ps, self)
        elif step_index == _DETECT:
            from apex.gui.workflow.step4_source_detection import SourceDetectionWindow
            return SourceDetectionWindow(p, fm, ps, self)

        # ── Step 5: WCS (shared) ──
        elif step_index == _WCS:
            from apex.gui.workflow.step5_wcs_plate_solving import WcsPlateSolvingWindow
            return WcsPlateSolvingWindow(p, fm, ps, self)

        # ── Step 6: MasterBuild (shared) ──
        elif step_index == _MASTER:
            from apex.gui.workflow.step6_ref_build import RefBuildWindow
            return RefBuildWindow(p, fm, ps, self)

        # ── Step 7: Forced Aperture Phot (shared) ──
        elif step_index == _FORCED:
            from apex.gui.workflow.step7_forced_aperture_phot import ForcedPhotWindow
            return ForcedPhotWindow(p, fm, ps, self)

        # ── Step 8+: mode-specific ──
        elif step_index == 7:
            from apex.gui.workflow.cmd.step8_psf_photometry import PSFPhotometryWindow
            return PSFPhotometryWindow(p, fm, ps, self)

        elif step_index == 8:
            if self.mode == "cmd":
                from apex.gui.workflow.cmd.step9_master_id_editor import MasterIdEditorWindow
                return MasterIdEditorWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.lc.step8_target_selection import TargetComparisonSelectionWindow
                return TargetComparisonSelectionWindow(p, fm, ps, self)

        elif step_index == 9:
            if self.mode == "cmd":
                from apex.gui.workflow.cmd.step10_zeropoint_calibration import ZeropointCalibrationWindow
                return ZeropointCalibrationWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.lc.step9_lightcurve_builder import LightCurveBuilderWindow
                return LightCurveBuilderWindow(p, fm, ps, self)

        elif step_index == 10:
            if self.mode == "cmd":
                from apex.gui.workflow.cmd.step11_cmd_plot import CmdPlotWindow
                return CmdPlotWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.lc.step10_detrend_merge import DetrendNightMergeWindow
                return DetrendNightMergeWindow(p, fm, ps, self)

        elif step_index == 11:
            if self.mode == "cmd":
                from apex.gui.workflow.cmd.step12_isochrone_model import IsochroneModelWindow
                return IsochroneModelWindow(p, fm, ps, self)
            from apex.gui.workflow.lc.step11_period_analysis import PeriodAnalysisWindow
            return PeriodAnalysisWindow(p, fm, ps, self)

        return None

    def on_step_completed(self, step_index: int):
        self.project_state.mark_step_completed(step_index)
        self.update_step_buttons()
        self.append_log(f"\u2713 Step {step_index + 1} finished: {self.step_names[step_index]}")

    def load_previous_session(self):
        """Pick a saved project_state.json (typically mirrored under a previous
        run's result_dir) and rebind this window to that session."""
        result_dir = getattr(self.params.P, "result_dir", None)
        start_dir = str(Path(result_dir).parent) if result_dir else str(Path.home())
        chosen, _ = QFileDialog.getOpenFileName(
            self,
            "Load Previous Session — select project_state.json",
            start_dir,
            "APEX session (project_state.json);;JSON files (*.json);;All files (*)",
        )
        if not chosen:
            return

        chosen_path = Path(chosen)
        try:
            import json as _json
            with open(chosen_path, "r", encoding="utf-8") as f:
                loaded = _json.load(f)
            if not isinstance(loaded, dict) or "step_data" not in loaded:
                raise ValueError("File does not look like an APEX project_state.json")
        except Exception as exc:
            QMessageBox.critical(
                self, "Load Session Failed",
                f"Could not read session file:\n  {chosen_path}\n\n{exc}",
            )
            return

        if self.current_step_window is not None:
            try:
                self.current_step_window.close()
            except Exception:
                pass
            self.current_step_window = None

        # Replace state contents in place, then persist (so the current
        # workspace state file mirrors the loaded session) and re-bootstrap.
        self.project_state.state.update(loaded)
        self.project_state._normalize_state()
        self.project_state.clear_mirror_paths()
        if chosen_path.resolve() != self.project_state.state_file.resolve():
            self.project_state.add_mirror_path(chosen_path)
        self.project_state.save()

        # Re-derive params.data_dir/result_dir from the loaded file_selection
        # state, then rebuild file_manager so downstream steps see the swap.
        if hasattr(self, "_offline_data_dir"):
            del self._offline_data_dir
        self._bootstrap_file_selection_state()
        try:
            from apex.core import FileManager
            self.file_manager = FileManager(self.params)
            # _bootstrap re-applies multi-night/ref_filename to the new manager
            self._bootstrap_file_selection_state()
        except Exception:
            pass

        self.update_step_buttons()
        self.append_log(f"Loaded session: {chosen_path}")
        result_dir_now = getattr(self.params.P, "result_dir", "(unknown)")
        completed = len(self.project_state.state.get("completed_steps", []))
        QMessageBox.information(
            self, "Session Loaded",
            f"Result dir: {result_dir_now}\n"
            f"Completed steps: {completed}/{len(self.step_names)}",
        )

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

    def _write_tool_error_log(self, tool_name: str, traceback_text: str) -> Path | None:
        """Persist tool-launch tracebacks for GUI sessions without a console."""
        try:
            from datetime import datetime
            log_dir = Path(getattr(self.params.P, "result_dir", Path.cwd())) / "logs"
            log_dir.mkdir(parents=True, exist_ok=True)
            log_path = log_dir / "apex_tool_errors.log"
            stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            with log_path.open("a", encoding="utf-8") as fh:
                fh.write(f"\n[{stamp}] {tool_name}\n")
                fh.write(traceback_text.rstrip())
                fh.write("\n")
            return log_path
        except Exception:
            return None

    def show_about(self):
        mode_label = "CMD Cluster Photometry" if self.mode == "cmd" else "LC Light Curve Analysis"
        QMessageBox.about(
            self, "About APEX",
            f"<h2>APEX — Automated Photometry EXtraction</h2>"
            f"<p><b>Mode: {mode_label}</b></p>"
            "<p>KNUEMAO Observatory — CDK500 + Moravian C3-61000</p>"
            "<p>Version 2.0.0</p>"
        )

    def show_wcs_solver_help(self):
        from PyQt5.QtWidgets import QDialog, QDialogButtonBox, QTextBrowser

        dialog = FittedDialog(self)
        dialog.setWindowTitle("WCS Solver Installation Help")
        dialog.resize(760, 620)

        layout = QVBoxLayout(dialog)
        browser = QTextBrowser(dialog)
        browser.setOpenExternalLinks(True)
        browser.setHtml(
            """
            <h2>WCS Solver Installation Help</h2>
            <p>
            APEX does not bundle external WCS solvers or their star/index
            databases. Install at least ASTAP for normal Step 5 operation;
            local astrometry.net is an optional fallback when ASTAP fails.
            </p>

            <h3>ASTAP (recommended first solver)</h3>
            <ol>
              <li>Install ASTAP for Windows.</li>
              <li>
                Install one ASTAP star database matching the APEX parameter.
                APEX currently passes <b>D80</b> or <b>D50</b> to ASTAP.
                Use D80 as the default; D50 is smaller and is normally useful
                when the field of view is comfortably above about 0.2 deg.
              </li>
              <li>
                In Step 5 &gt; ASTAP Parameters, set <b>ASTAP CLI Path</b>
                to <code>astap_cli.exe</code> if it is not already on PATH,
                and set <b>ASTAP Star DB</b> to the installed database.
              </li>
            </ol>
            <p>
              ASTAP:
              <a href="https://www.hnsky.org/astap.htm">https://www.hnsky.org/astap.htm</a><br>
              ASTAP star databases:
              <a href="https://sourceforge.net/projects/astap-program/files/star_databases/">
              https://sourceforge.net/projects/astap-program/files/star_databases/</a>
            </p>

            <h3>Local astrometry.net / solve-field (optional fallback)</h3>
            <ol>
              <li>Install WSL/Ubuntu or another local Linux environment.</li>
              <li>
                Install astrometry.net so that <code>solve-field</code> runs
                from the shell APEX will call.
              </li>
              <li>
                Install astrometry.net index files that match the image field
                of view. Missing or mismatched indexes are the most common
                cause of no-solution failures.
              </li>
              <li>
                In Step 5 &gt; Astrometry.net Parameters, enable
                <b>Use WSL</b> for Windows/WSL setups and leave the command as
                <code>solve-field</code> unless your installation needs a
                different wrapper.
              </li>
            </ol>
            <p>
              Astrometry.net README:
              <a href="https://astrometry.net/doc/readme.html">
              https://astrometry.net/doc/readme.html</a><br>
              Astrometry.net index files:
              <a href="https://data.astrometry.net/">
              https://data.astrometry.net/</a>
            </p>

            <h3>APEX solve flow</h3>
            <p>
            Step 5 tries ASTAP first. If ASTAP fails or does not leave a valid
            WCS header, and local astrometry.net is enabled, APEX calls
            <code>solve-field</code> as a fallback. Both solvers then pass
            through the same WCS header validation and QC summary path.
            </p>

            <h3>Startup checks</h3>
            <p>
            Before solving, APEX checks whether <code>astap_cli.exe</code> is
            reachable and whether <code>solve-field</code> is reachable in the
            configured local/WSL environment. APEX also checks whether
            <code>astroquery.gaia</code> is importable for Gaia attach, WCS
            refine, and residual statistics. Star databases and astrometry.net
            index coverage are still external data requirements; if the solver
            executable is found but solving fails, check the selected D80/D50
            database or the installed index-file scales.
            </p>
            """
        )
        layout.addWidget(browser, 1)

        buttons = QDialogButtonBox(QDialogButtonBox.Close)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        dialog.exec_()

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

    def open_detector_ptc_tool(self):
        from apex.gui.tools.detector_ptc import DetectorPTCWindow
        self.detector_ptc_window = DetectorPTCWindow(
            self.params, self.project_state, parent=None
        )
        self.detector_ptc_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.detector_ptc_window.show()
        self.detector_ptc_window.raise_()
        self.detector_ptc_window.activateWindow()
        self.append_log("Opened Detector Characterisation (PTC) Tool")

    # ── LC tools ─────────────────────────────────────────────────────────────

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

    def open_variable_star_tool(self, handoff: dict | None = None):
        from apex.gui.tools.variable_star import VariableStarToolWindow
        if self.varstar_window is None:
            self.varstar_window = VariableStarToolWindow(
                self.params, self.project_state, parent=None
            )
            self.varstar_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        if handoff:
            self.varstar_window.apply_period_handoff(handoff)
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
        self.open_extinction_tool()

    def open_gaia_3d_viewer(self):
        try:
            from apex.gui.tools.gaia_3d_viewer import Gaia3DViewerWindow
            self.gaia_3d_window = Gaia3DViewerWindow(
                self.params, self.params.P.result_dir, parent=None
            )
            self.gaia_3d_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
            self.gaia_3d_window.show()
            self.gaia_3d_window.raise_()
            self.gaia_3d_window.activateWindow()
            self.append_log(f"Opened Gaia 3D Viewer: {self.params.P.result_dir}")
        except Exception as e:
            import traceback
            tb = traceback.format_exc()
            log_path = self._write_tool_error_log("Gaia 3D Cluster Viewer", tb)
            short = f"{type(e).__name__}: {e}"
            self.append_log(f"Failed to open Gaia 3D Viewer: {short}")
            if log_path is not None:
                self.append_log(f"Tool error log: {log_path}")
            msg = QMessageBox(self)
            msg.setIcon(QMessageBox.Critical)
            msg.setWindowTitle("Gaia 3D Viewer Error")
            msg.setText("Failed to open Gaia 3D Cluster Viewer.")
            detail = f"{short}\n\n"
            if log_path is not None:
                detail += f"Log written to:\n{log_path}"
            msg.setInformativeText(detail)
            msg.setDetailedText(tb)
            msg.exec_()

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
            QDialogButtonBox, QGroupBox, QFormLayout, QCheckBox,
            QScrollArea, QFrame,
        )
        dialog = FittedDialog(self)
        dialog.setWindowTitle("Instrument Settings")
        dialog.setMinimumWidth(650)
        layout = QVBoxLayout(dialog)
        title = QLabel("Instrument Configuration")
        title.setFont(QFont("Arial", 12, QFont.Bold))
        layout.addWidget(title)

        # Group boxes live inside a scroll area so a screen too short for all
        # six groups scrolls instead of clipping the Save/Cancel row below.
        _settings_body = QWidget()
        _settings_layout = QVBoxLayout(_settings_body)
        _settings_layout.setContentsMargins(0, 0, 0, 0)
        _settings_scroll = QScrollArea()
        _settings_scroll.setWidgetResizable(True)
        _settings_scroll.setFrameShape(QFrame.NoFrame)
        _settings_scroll.setWidget(_settings_body)
        layout.addWidget(_settings_scroll, 1)
        _settings_outer = layout      # holds the footer button row
        layout = _settings_layout     # subsequent group adds go into the scroll body

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
        gain_value = getattr(self.params.P, "gain_e_per_adu", None)
        gain_edit = QLineEdit("" if gain_value is None else str(gain_value))
        gain_edit.setMaximumWidth(150)
        params_layout.addRow("Gain (e-/ADU):", gain_edit)
        rdnoise_value = getattr(self.params.P, "rdnoise_e", None)
        rdnoise_edit = QLineEdit("" if rdnoise_value is None else str(rdnoise_value))
        rdnoise_edit.setMaximumWidth(150)
        params_layout.addRow("Read Noise (e-):", rdnoise_edit)
        noise_header_check = QCheckBox()
        noise_header_check.setChecked(bool(getattr(self.params.P, "noise_use_fits_header", False)))
        noise_header_check.setToolTip("Trust FITS EGAIN/RDNOISE when present. Leave off when measured PTC values should override bad headers.")
        params_layout.addRow("Trust FITS EGAIN/RDNOISE:", noise_header_check)
        noise_ref_bin_spin = QSpinBox()
        noise_ref_bin_spin.setRange(1, 8)
        noise_ref_bin_spin.setValue(int(getattr(
            self.params.P,
            "noise_reference_binning",
            getattr(self.params.P, "binning_default", self.instrument.binning),
        ) or getattr(self.params.P, "binning_default", self.instrument.binning)))
        noise_ref_bin_spin.setToolTip("Binning of the manual gain/read-noise values. Use 1 for BIN1 camera specs.")
        params_layout.addRow("Manual Noise Ref. Binning:", noise_ref_bin_spin)
        noise_scale_check = QCheckBox()
        noise_scale_check.setChecked(bool(getattr(self.params.P, "noise_scale_by_binning", False)))
        noise_scale_check.setToolTip("Scale manual gain/read-noise from the reference binning to each FITS frame binning. Leave off for PTC values measured from the current FITS data.")
        params_layout.addRow("Scale Manual Noise by Binning:", noise_scale_check)
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
        _settings_outer.addWidget(buttons)

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
                self.params.P.telescope_focal_mm = float(self.instrument.focal_length_mm)
                self.params.P.camera_pixel_um = float(self.instrument.pix_size_um)
                self.params.P.binning_default = int(self.instrument.binning)
                self.params.P.pixel_scale_arcsec = (
                    206.265
                    * float(self.instrument.pix_size_um)
                    * float(self.instrument.binning)
                    / float(self.instrument.focal_length_mm)
                )
                gain_text = gain_edit.text().strip()
                rdnoise_text = rdnoise_edit.text().strip()
                self.params.P.gain_e_per_adu = float(gain_text) if gain_text else None
                self.params.P.rdnoise_e = float(rdnoise_text) if rdnoise_text else None
                self.params.P.noise_use_fits_header = bool(noise_header_check.isChecked())
                self.params.P.noise_reference_binning = int(noise_ref_bin_spin.value())
                self.params.P.noise_scale_by_binning = bool(noise_scale_check.isChecked())
                self.params.P.saturation_adu = float(saturation_edit.text())
                self.params.P.site_lat_deg = float(site_lat_edit.text())
                self.params.P.site_lon_deg = float(site_lon_edit.text())
                self.params.P.site_alt_m = float(site_alt_edit.text())
                self.params.P.site_tz_offset_hours = float(site_tz_edit.text())
                self.params.P.max_workers = int(parallel_workers_spin.value())
                self.params.P.parallel_max_workers = int(parallel_workers_spin.value())
                saved = self.params.save_toml()
                if not saved:
                    self.append_log("Warning: could not save settings to TOML")
                QMessageBox.information(self, "Settings Saved",
                                        "Instrument settings updated." if saved else
                                        "Instrument settings updated, but TOML save failed.")
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
