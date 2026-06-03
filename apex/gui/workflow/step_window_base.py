"""
Base class for step windows
Each step gets its own window with Previous/Next navigation
"""

from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QMessageBox, QApplication
)
from PyQt5.QtCore import Qt, pyqtSignal, QThread
from PyQt5.QtGui import QFont

from apex.gui.window_chrome import WindowChromeMixin


class StepWindowBase(WindowChromeMixin, QMainWindow):
    """Base class for workflow step windows with navigation buttons"""

    # Signals
    step_completed = pyqtSignal(int)  # Emit step index when completed
    go_to_step = pyqtSignal(int)  # Request to go to another step

    def __init__(self, step_index: int, step_name: str, params, project_state, main_window, parent=None):
        """
        Initialize step window

        Args:
            step_index: Index of this step (0-based)
            step_name: Display name of this step
            params: Parameters object
            project_state: ProjectState object
            main_window: Reference to main window
            parent: Parent widget
        """
        super().__init__(parent)

        self.step_index = step_index
        self.step_name = step_name
        self.params = params
        self.project_state = project_state
        self.main_window = main_window

        self.completed = False

        # Connect to main window
        self.step_completed.connect(main_window.on_step_completed)

        # Setup base UI
        self.setWindowTitle(f"Step {step_index + 1}: {step_name}")
        self.setMinimumSize(900, 700)
        # Default size = 85% of available screen, capped at 1200×950 so
        # the window opens comfortably without being annoying on very
        # large monitors.  Slightly narrower default than the screen so
        # CMD/isochrone plots don't get stretched horizontally past
        # usefulness.
        try:
            screen = QApplication.primaryScreen()
            if screen is not None:
                avail = screen.availableGeometry()
                w = max(900, min(1200, int(avail.width() * 0.85)))
                h = max(700, min(950, int(avail.height() * 0.85)))
                self.resize(w, h)
        except Exception:
            pass

        # Central widget
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_layout = QVBoxLayout(self.central_widget)

        # Title row: centered title + a right-aligned action cluster.
        # Subclasses register buttons via add_header_action()/add_param_button()
        # instead of hand-rolling their own toolbar in each window.
        title_row = QHBoxLayout()
        self.title_label = QLabel(f"Step {step_index + 1}: {step_name}")
        self.title_label.setFont(QFont("Arial", 16, QFont.Bold))
        self.title_label.setAlignment(Qt.AlignCenter)
        # A spacer matching the action cluster keeps the title visually centered.
        title_row.addSpacing(110)
        title_row.addWidget(self.title_label, stretch=1)

        # Action cluster: subclass-added actions appear left of the guide button.
        self.header_actions = QHBoxLayout()
        self.header_actions.setSpacing(6)
        title_row.addLayout(self.header_actions)

        self.btn_guide = QPushButton("❔ 가이드")
        self.btn_guide.setFixedWidth(100)
        self.btn_guide.setMinimumHeight(30)
        self.btn_guide.setToolTip("이 단계에서 무엇을 하는지 간단히 안내합니다")
        self.btn_guide.setStyleSheet(
            "QPushButton { background-color: #1976D2; color: white;"
            " border-radius: 5px; font-weight: bold; padding: 4px 8px; }"
            "QPushButton:hover { background-color: #1565C0; }"
        )
        self.btn_guide.clicked.connect(self.show_step_guide)
        title_row.addWidget(self.btn_guide)
        self.main_layout.addLayout(title_row)

        # Optional skeleton slots (created lazily by helper methods below).
        self._log_text = None
        self._log_window = None

        # Content area (to be filled by subclasses).  stretch=1 so the
        # content fills the window whenever the user resizes — without
        # this, subclasses whose inner widgets have small sizeHint() values
        # (e.g. LC step8) render in a tiny corner of an otherwise empty
        # window.
        self.content_widget = QWidget()
        self.content_layout = QVBoxLayout(self.content_widget)
        self.main_layout.addWidget(self.content_widget, stretch=1)

        # Navigation buttons
        self.setup_navigation_buttons()

        # NOTE: restore_state() should be called by subclass after UI setup

    def show_step_guide(self):
        """Show the short per-step guide in a popup.

        Steps may override ``guide_text()`` to customise; by default the
        text comes from the central :mod:`step_guides` table keyed by
        (mode, step_index).
        """
        text = self.guide_text()
        box = QMessageBox(self)
        box.setWindowTitle(f"가이드 — Step {self.step_index + 1}: {self.step_name}")
        box.setTextFormat(Qt.RichText)
        box.setText(text)
        box.setStandardButtons(QMessageBox.Ok)
        box.setIcon(QMessageBox.Information)
        box.exec_()

    def guide_text(self) -> str:
        """Return the HTML guide for this step (override to customise)."""
        from apex.gui.workflow.step_guides import step_guide
        mode = getattr(self.main_window, "mode", "cmd")
        return step_guide(mode, self.step_index, self.step_name)

    # NOTE: header actions, the workspace bar, the log dock and log() come
    # from WindowChromeMixin (apex/gui/window_chrome.py) — shared with tools.

    def setup_navigation_buttons(self):
        """Setup Previous/Next navigation buttons at bottom."""
        nav_layout = QHBoxLayout()

        # Previous button (GREEN - always enabled if not first step)
        self.btn_previous = QPushButton("← Previous Step")
        self.btn_previous.setMinimumHeight(40)
        self.btn_previous.setFont(QFont("Arial", 10, QFont.Bold))
        self.btn_previous.clicked.connect(self.go_previous)

        if self.step_index > 0:
            # Green style for Previous
            self.btn_previous.setStyleSheet("""
                QPushButton {
                    background-color: #4CAF50;
                    color: white;
                    border: 2px solid #45a049;
                    border-radius: 5px;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #45a049;
                }
            """)
            self.btn_previous.setEnabled(True)
        else:
            # Disabled for Step 1
            self.btn_previous.setEnabled(False)

        nav_layout.addWidget(self.btn_previous)

        nav_layout.addStretch()

        # Keep a hidden placeholder for compatibility with embedded tools that
        # still expect this attribute to exist, but remove it from the UI.
        self.btn_complete = QPushButton("Mark as Complete")
        self.btn_complete.setMinimumHeight(40)
        self.btn_complete.setFont(QFont("Arial", 10, QFont.Bold))
        self.btn_complete.setEnabled(False)
        self.btn_complete.hide()
        self.btn_complete.clicked.connect(self.mark_complete)

        # Next / Exit button (RED when not ready, GREEN when ready)
        self._is_last_step = (self.step_index >= len(self.main_window.step_names) - 1)
        self.btn_next = QPushButton("Exit ✕" if self._is_last_step else "Next Step →")
        self.btn_next.setMinimumHeight(40)
        self.btn_next.setFont(QFont("Arial", 10, QFont.Bold))
        self.btn_next.setEnabled(False)
        self.btn_next.clicked.connect(self.go_next)

        # Initial red style (not complete)
        self.btn_next.setStyleSheet("""
            QPushButton {
                background-color: #f44336;
                color: white;
                border: 2px solid #da190b;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:disabled {
                background-color: #CCCCCC;
                color: #666666;
                border: 2px solid #AAAAAA;
            }
        """)

        nav_layout.addWidget(self.btn_next)

        self.main_layout.addLayout(nav_layout)

        # Update button states
        self.update_navigation_buttons()

    def update_navigation_buttons(self):
        """Update navigation button states."""
        can_proceed = bool(self.validate_step())
        self.btn_complete.setEnabled(False)
        self.btn_complete.hide()
        self.btn_next.setEnabled(can_proceed)

        if can_proceed:
            self.btn_next.setStyleSheet("""
                QPushButton {
                    background-color: #4CAF50;
                    color: white;
                    border: 2px solid #45a049;
                    border-radius: 5px;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #45a049;
                }
            """)
        else:
            self.btn_next.setStyleSheet("""
                QPushButton {
                    background-color: #f44336;
                    color: white;
                    border: 2px solid #da190b;
                    border-radius: 5px;
                    font-weight: bold;
                }
                QPushButton:disabled {
                    background-color: #CCCCCC;
                    color: #666666;
                    border: 2px solid #AAAAAA;
                }
            """)

    def _mark_step_complete(self, notify_main: bool = True) -> bool:
        """Mark the current step complete if needed."""
        if self.project_state.is_step_completed(self.step_index):
            self.completed = True
            return False

        self.completed = True
        if notify_main:
            self.step_completed.emit(self.step_index)
        else:
            self.project_state.mark_step_completed(self.step_index)
            if hasattr(self.main_window, "update_step_buttons"):
                self.main_window.update_step_buttons()
        return True

    def _mark_step_incomplete(self) -> bool:
        """Clear completion state for the current step if present."""
        if not self.project_state.is_step_completed(self.step_index):
            self.completed = False
            return False

        self.project_state.mark_step_incomplete(self.step_index)
        self.completed = False
        if hasattr(self.main_window, "update_step_buttons"):
            self.main_window.update_step_buttons()
        return True

    def _finalize_valid_step(self, notify_main: bool = True) -> bool:
        """Persist state/params and mark the step completed when valid."""
        if not self.validate_step():
            return False

        self.save_state()
        self.persist_params()
        self._mark_step_complete(notify_main=notify_main)
        return True

    def validate_step(self) -> bool:
        """
        Validate if step can be marked as complete
        Override in subclasses

        Returns:
            True if step can be completed
        """
        return False

    def persist_params(self) -> bool:
        """Persist current parameters to TOML, if supported."""
        if hasattr(self.params, "save_toml"):
            try:
                return bool(self.params.save_toml())
            except Exception:
                return False
        return False

    def _iter_qthreads(self):
        """Yield unique QThread instances attached to this window."""
        seen = set()
        for name, value in vars(self).items():
            if not isinstance(value, QThread):
                continue
            ident = id(value)
            if ident in seen:
                continue
            seen.add(ident)
            yield name, value

    def _cleanup_running_threads(self, timeout_ms: int = 5000) -> bool:
        """Request stop for running worker threads and release finished ones."""
        all_stopped = True

        for name, thread in self._iter_qthreads():
            if thread.isRunning():
                stop_fn = getattr(thread, "stop", None)
                if callable(stop_fn):
                    try:
                        stop_fn()
                    except Exception:
                        pass
                thread.quit()
                if not thread.wait(int(timeout_ms)):
                    all_stopped = False
                    continue

            try:
                thread.deleteLater()
            except Exception:
                pass

            try:
                setattr(self, name, None)
            except Exception:
                pass

        return all_stopped

    def mark_complete(self):
        """Backward-compatible alias for the removed manual-complete action."""
        if not self._finalize_valid_step(notify_main=True):
            QMessageBox.warning(
                self, "Step Not Ready",
                "Please finish all required tasks before proceeding."
            )
            return
        self.update_navigation_buttons()

    def go_previous(self):
        """Go to previous step"""
        if self.step_index > 0:
            self.close()
            self.main_window.open_step(self.step_index - 1)

    def go_next(self):
        """Go to next step, or close (exit) if this is the last step."""
        if not self._finalize_valid_step(notify_main=True):
            QMessageBox.warning(
                self, "Step Not Ready",
                "Please finish all required tasks before proceeding to the next step."
            )
            return

        self.close()
        if not self._is_last_step:
            self.main_window.open_step(self.step_index + 1)

    def save_state(self):
        """
        Save step-specific state
        Override in subclasses to save data
        """
        pass

    def restore_state(self):
        """
        Restore step-specific state
        Override in subclasses to restore data
        """
        pass

    def closeEvent(self, event):
        """Handle window close"""
        if not self._cleanup_running_threads(timeout_ms=5000):
            QMessageBox.warning(
                self,
                "Background Task Running",
                "A background worker is still stopping. Please wait a few seconds and close again."
            )
            event.ignore()
            return

        was_completed = self.project_state.is_step_completed(self.step_index)
        if self.validate_step():
            self._finalize_valid_step(notify_main=not was_completed)
        else:
            self._mark_step_incomplete()
            self.persist_params()

        event.accept()
