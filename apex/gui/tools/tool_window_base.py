"""Base window for standalone analysis tools.

Tools historically each rolled their own QWidget/QMainWindow with bespoke
headers, directory pickers and log panels — the source of the "tools 중구난방"
feel. ToolWindowBase gives them the same standard chrome the step windows get
(via WindowChromeMixin): a centered title + header action cluster, an optional
workspace/IO bar, an optional log dock, and a standard Parameters button.

Tools do NOT have Previous/Next navigation or step-completion — they are
standalone, so this base is deliberately lighter than StepWindowBase.

Usage:
    class MyToolWindow(ToolWindowBase):
        def __init__(self, params, project_state, parent=None):
            super().__init__("My Tool", params=params,
                             project_state=project_state, parent=parent)
            self.add_workspace_bar()
            self.add_log_dock(popup=False)
            self.content_layout.addWidget(self._build_body())
"""

from __future__ import annotations

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import (
    QHBoxLayout, QLabel, QMainWindow, QVBoxLayout, QWidget,
)

from apex.gui.window_chrome import WindowChromeMixin


class ToolWindowBase(WindowChromeMixin, QMainWindow):
    """Standalone tool window with standard chrome (no step navigation)."""

    def __init__(self, title, *, params=None, project_state=None, parent=None,
                 min_size=(820, 600)):
        super().__init__(parent)
        self.params = params
        self.project_state = project_state
        self._chrome_title = title

        self.setWindowTitle(title)
        if min_size:
            # Clamp the minimum to the monitor so a large min_size (e.g. a
            # 1180×820 plot tool) can never exceed a small screen and leave the
            # bottom controls permanently clipped. First-show auto-fit
            # (WindowChromeMixin.showEvent) then grows toward the content.
            from apex.gui.layout_rules import clamp_to_screen
            min_w, min_h = clamp_to_screen(int(min_size[0]), int(min_size[1]), self)
            self.setMinimumSize(min_w, min_h)

        from apex.gui.theme import Tokens

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_layout = QVBoxLayout(self.central_widget)
        # Same 8px outer rhythm as step windows so tools feel like one app.
        self.main_layout.setContentsMargins(
            Tokens.MARGIN, Tokens.MARGIN, Tokens.MARGIN, Tokens.S3
        )
        self.main_layout.setSpacing(Tokens.S3)

        # Title row: centered title + right-aligned action cluster.
        title_row = QHBoxLayout()
        self.title_label = QLabel(title)
        self.title_label.setFont(QFont("Arial", 16, QFont.Bold))
        self.title_label.setAlignment(Qt.AlignCenter)
        title_row.addSpacing(110)
        title_row.addWidget(self.title_label, stretch=1)

        self.header_actions = QHBoxLayout()
        self.header_actions.setSpacing(Tokens.GAP)
        title_row.addLayout(self.header_actions)
        self.main_layout.addLayout(title_row)

        # Content area (subclasses fill this).
        self.content_widget = QWidget()
        self.content_layout = QVBoxLayout(self.content_widget)
        self.main_layout.addWidget(self.content_widget, stretch=1)

        # Lazy chrome slots (created by mixin helpers when called).
        self._log_text = None
        self._log_window = None

    def persist_params(self) -> bool:
        """Persist current parameters to TOML, if supported (used by the
        standard Parameters button)."""
        params = getattr(self, "params", None)
        if params is not None and hasattr(params, "save_toml"):
            try:
                return bool(params.save_toml())
            except Exception:
                return False
        return False
