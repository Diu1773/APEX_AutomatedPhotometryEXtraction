"""Shared workflow log window widgets."""

from __future__ import annotations

import time

from pathlib import Path

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QHBoxLayout, QLabel, QProgressBar, QSplitter,
    QTextEdit, QVBoxLayout, QWidget,
)


def create_log_text(parent: QWidget | None = None) -> QTextEdit:
    text = QTextEdit(parent)
    text.setReadOnly(True)
    text.setStyleSheet("QTextEdit { font-family: monospace; font-size: 9pt; }")
    return text


def append_timestamped_log(log_text: QTextEdit | None, message: str) -> None:
    if log_text is None:
        return
    timestamp = time.strftime("%H:%M:%S")
    log_text.append(f"[{timestamp}] {message}")


def show_raised(window: QWidget | None) -> None:
    if window is None:
        return
    window.show()
    window.raise_()
    window.activateWindow()


class WorkflowLogWindow(QWidget):
    """Reusable workflow log window with an optional right-side panel."""

    def __init__(
        self,
        parent: QWidget | None,
        title: str,
        *,
        width: int = 800,
        height: int = 400,
        side_widget: QWidget | None = None,
        log_stretch: int = 2,
        side_stretch: int = 1,
    ):
        super().__init__(parent, Qt.Window)
        self.setWindowTitle(title)
        self.resize(width, height)

        layout = QVBoxLayout(self)
        self.log_text = create_log_text(self)

        if side_widget is None:
            layout.addWidget(self.log_text)
            return

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(self.log_text)
        splitter.addWidget(side_widget)
        splitter.setStretchFactor(0, int(log_stretch))
        splitter.setStretchFactor(1, int(side_stretch))
        splitter.setChildrenCollapsible(False)
        layout.addWidget(splitter)

    def append(self, message: str) -> None:
        append_timestamped_log(self.log_text, message)

    def clear(self) -> None:
        self.log_text.clear()

    def show_raised(self) -> None:
        show_raised(self)


class WorkerStatusPanel(QWidget):
    """Reusable per-worker status/progress panel.

    Renders one row per worker_id: ``W## <fname>  <status>  [progress bar]``.
    Use ``update_worker(worker_id, filename, status, progress)`` from a slot
    connected to a worker signal of shape ``(int, str, str, int)``.
    """

    _MONO_QSS = "QLabel { font-family: Consolas, 'Courier New', monospace; font-size: 9pt; }"

    def __init__(self, parent: QWidget | None = None, *, min_width: int = 430):
        super().__init__(parent)
        self.setMinimumWidth(min_width)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(5, 5, 5, 5)
        outer.setSpacing(2)
        self._rows_layout = QVBoxLayout()
        outer.addLayout(self._rows_layout)
        outer.addStretch()
        self._rows: dict[int, tuple[QLabel, QLabel, QLabel, QProgressBar]] = {}

    def update_worker(self, worker_id: int, filename: str, status: str, progress: int) -> None:
        if worker_id not in self._rows:
            row = QWidget()
            row_layout = QHBoxLayout(row)
            row_layout.setContentsMargins(0, 0, 0, 0)
            row_layout.setSpacing(4)

            id_lbl = QLabel()
            id_lbl.setFixedWidth(32)
            id_lbl.setStyleSheet(self._MONO_QSS)

            fname_lbl = QLabel()
            fname_lbl.setStyleSheet(self._MONO_QSS)
            fname_lbl.setMinimumWidth(160)

            status_lbl = QLabel()
            status_lbl.setFixedWidth(110)
            status_lbl.setStyleSheet(self._MONO_QSS)

            bar = QProgressBar()
            bar.setRange(0, 100)
            bar.setValue(0)
            bar.setTextVisible(True)
            bar.setFixedWidth(80)

            row_layout.addWidget(id_lbl)
            row_layout.addWidget(fname_lbl, stretch=1)
            row_layout.addWidget(status_lbl)
            row_layout.addWidget(bar)
            self._rows_layout.addWidget(row)
            self._rows[worker_id] = (id_lbl, fname_lbl, status_lbl, bar)

        id_lbl, fname_lbl, status_lbl, bar = self._rows[worker_id]
        short_fname = Path(str(filename)).name
        id_lbl.setText(f"W{int(worker_id):02d}")
        fm = fname_lbl.fontMetrics()
        elided = fm.elidedText(short_fname, Qt.ElideMiddle, fname_lbl.width() or 160)
        fname_lbl.setText(elided)
        fname_lbl.setToolTip(short_fname)
        status_lbl.setText(str(status))
        bar.setValue(int(progress))

    def clear(self) -> None:
        while self._rows_layout.count():
            item = self._rows_layout.takeAt(0)
            w = item.widget()
            if w is not None:
                w.setParent(None)
        self._rows.clear()
