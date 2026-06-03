"""Reusable window chrome shared by step windows and tool windows.

This is the single home for the three things every APEX window used to
re-implement by hand: a header action cluster (incl. the Parameters button),
a workspace / input-output context bar, and a log panel.

A host class mixes this in and must provide three layouts before calling any
helper:
    self.main_layout      (QVBoxLayout)  — the window's outer vertical layout
    self.content_layout   (QVBoxLayout)  — where the page body goes
    self.header_actions   (QHBoxLayout)  — the title-row action cluster

Optional host attributes the helpers will use if present:
    self.params           — Parameters object (needed for params button / workspace bar)
    self._chrome_title    — str used as a default dialog/log title
    self.persist_params   — callable used as the default param-save callback
"""

from __future__ import annotations

from pathlib import Path

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QFileDialog, QFrame, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget,
)

from apex.gui.workflow.log_panel import (
    append_timestamped_log, create_log_text, show_raised,
)


class WindowChromeMixin:
    """Standard header / workspace bar / log dock for step & tool windows."""

    _HEADER_BTN_QSS = (
        "QPushButton { background-color: #ECEFF1; color: #263238;"
        " border: 1px solid #B0BEC5; border-radius: 5px;"
        " font-weight: bold; padding: 4px 10px; }"
        "QPushButton:hover { background-color: #CFD8DC; }"
    )

    # ── internal helpers ──────────────────────────────────────────────────
    def _chrome_default_title(self) -> str:
        return (
            getattr(self, "_chrome_title", "")
            or getattr(self, "step_name", "")
            or self.windowTitle()
            or "Window"
        )

    # ── header actions ────────────────────────────────────────────────────
    def add_header_action(self, label, on_click, *, tooltip="", min_width=None):
        """Add a button to the title-row action cluster. Returns the button."""
        btn = QPushButton(label)
        btn.setMinimumHeight(30)
        if min_width:
            btn.setFixedWidth(int(min_width))
        if tooltip:
            btn.setToolTip(tooltip)
        btn.setStyleSheet(self._HEADER_BTN_QSS)
        if callable(on_click):
            btn.clicked.connect(on_click)
        self.header_actions.addWidget(btn)
        return btn

    def add_param_button(self, specs, *, title=None, info_text="",
                         label="⚙ 파라미터", on_saved=None, resize=(460, 400)):
        """Add a standard Parameters button wired to the shared param dialog."""
        from apex.gui.workflow.param_dialog import run_param_dialog

        dialog_title = title or f"{self._chrome_default_title()} — 파라미터"
        default_save = getattr(self, "persist_params", None)

        def _open():
            run_param_dialog(
                self, dialog_title, specs,
                on_save=(on_saved or default_save or (lambda: None)),
                info_text=info_text, resize=resize,
            )

        return self.add_header_action(label, _open, tooltip="이 창의 파라미터 설정")

    # ── info banner ───────────────────────────────────────────────────────
    def add_info_banner(self, text):
        """Standard one-line info banner at the top of the content area."""
        info = QLabel(text)
        info.setWordWrap(True)
        info.setStyleSheet(
            "QLabel { background-color: #E3F2FD; color: #0D47A1;"
            " padding: 8px 10px; border-radius: 5px; }"
        )
        self.content_layout.addWidget(info)
        return info

    # ── workspace / IO bar ────────────────────────────────────────────────
    def add_workspace_bar(self, *, allow_change=True, insert_index=1):
        """Compact IO context row: shows data/result dirs + an optional change.

        Replaces the per-window directory pickers scattered across windows.
        Override on_workspace_changed() to react to a folder change.
        """
        bar = QFrame()
        bar.setFrameShape(QFrame.StyledPanel)
        bar.setStyleSheet(
            "QFrame { background-color: #FAFAFA; border: 1px solid #E0E0E0;"
            " border-radius: 5px; }"
        )
        row = QHBoxLayout(bar)
        row.setContentsMargins(10, 5, 10, 5)
        row.setSpacing(8)

        self._workspace_label = QLabel()
        self._workspace_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        row.addWidget(self._workspace_label, stretch=1)

        if allow_change:
            btn = QPushButton("출력 폴더 변경…")
            btn.setStyleSheet(self._HEADER_BTN_QSS)
            btn.clicked.connect(self._change_workspace)
            row.addWidget(btn)

        self.main_layout.insertWidget(int(insert_index), bar)
        self._refresh_workspace_label()
        return bar

    def _refresh_workspace_label(self):
        if getattr(self, "_workspace_label", None) is None:
            return
        P = getattr(getattr(self, "params", None), "P", None)
        data_dir = getattr(P, "data_dir", "—") if P is not None else "—"
        result_dir = getattr(P, "result_dir", "—") if P is not None else "—"
        self._workspace_label.setText(
            f"📂 입력  {data_dir}      →      💾 출력  {result_dir}"
        )

    def _change_workspace(self):
        P = getattr(getattr(self, "params", None), "P", None)
        if P is None:
            return
        start = str(getattr(P, "result_dir", "") or getattr(P, "data_dir", "") or Path.home())
        chosen = QFileDialog.getExistingDirectory(self, "출력 폴더 선택", start)
        if not chosen:
            return
        P.result_dir = Path(chosen)
        P.cache_dir = P.result_dir / "cache"
        try:
            P.result_dir.mkdir(parents=True, exist_ok=True)
            P.cache_dir.mkdir(parents=True, exist_ok=True)
        except OSError:
            pass
        self._refresh_workspace_label()
        self.on_workspace_changed()

    def on_workspace_changed(self):
        """Hook called after the output folder changes (override as needed)."""
        pass

    # ── log dock ──────────────────────────────────────────────────────────
    def add_log_dock(self, *, popup=False, height=140, title=None):
        """Standard log panel. popup=True puts it in a separate window with a
        '로그' header action; otherwise it docks at the bottom.

        Provides self.log(message); replaces hand-rolled QTextEdit logs.
        """
        self._log_text = create_log_text(self)
        if popup:
            self._log_window = QWidget(self, Qt.Window)
            self._log_window.setWindowTitle(title or f"{self._chrome_default_title()} — Log")
            self._log_window.resize(700, 360)
            lay = QVBoxLayout(self._log_window)
            lay.addWidget(self._log_text)
            self.add_header_action(
                "📜 로그", lambda: show_raised(self._log_window),
                tooltip="처리 로그 보기",
            )
        else:
            box = QFrame()
            box.setFrameShape(QFrame.StyledPanel)
            lay = QVBoxLayout(box)
            lay.setContentsMargins(6, 6, 6, 6)
            hdr = QLabel(title or "로그")
            hdr.setStyleSheet("QLabel { color: #555; font-weight: bold; }")
            lay.addWidget(hdr)
            self._log_text.setFixedHeight(int(height))
            lay.addWidget(self._log_text)
            # Insert just above the last item (nav row for steps).
            self.main_layout.insertWidget(self.main_layout.count() - 1, box)
        return self._log_text

    def log(self, message: str):
        """Append a timestamped line to the standard log dock, if present."""
        append_timestamped_log(getattr(self, "_log_text", None), message)
