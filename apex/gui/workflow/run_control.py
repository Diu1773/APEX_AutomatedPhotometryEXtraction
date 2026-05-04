"""Shared Run / Stop / Log button bar for workflow step windows.

Usage:
    self.run_bar = RunControlBar(
        "Run Detection", "Log & Workers",
        run_cb=self.run_detection,
        stop_cb=self.stop_detection,
        log_cb=self.show_log_window,
    )
    layout.addWidget(self.run_bar)

    # start job
    self.run_bar.set_running(True)

    # mid-stop (waiting for worker to finish)
    self.run_bar.set_stopping()

    # job done
    self.run_bar.set_running(False)
"""
from __future__ import annotations

from PyQt5.QtWidgets import QHBoxLayout, QPushButton, QSizePolicy, QWidget

_RUN_SS = (
    "QPushButton { background-color: #4CAF50; color: white; font-weight: bold; padding: 8px 20px; }"
)
_STOP_SS = (
    "QPushButton { background-color: #f44336; color: white; font-weight: bold; padding: 8px 15px; }"
)
_LOG_SS = (
    "QPushButton { background-color: #607D8B; color: white; font-weight: bold; padding: 8px 15px; }"
)


class RunControlBar(QWidget):
    """[Run]──[Stop]──────────────────────[Log]

    btn_run and btn_stop are public for direct access when needed.
    Prefer set_running(bool) over touching them individually.
    """

    def __init__(
        self,
        run_label: str = "Run",
        log_label: str = "Log",
        *,
        run_cb,
        stop_cb,
        log_cb,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)

        self.btn_run = QPushButton(run_label)
        self.btn_run.setStyleSheet(_RUN_SS)
        self.btn_run.clicked.connect(run_cb)
        layout.addWidget(self.btn_run)

        self.btn_stop = QPushButton("Stop")
        self.btn_stop.setStyleSheet(_STOP_SS)
        self.btn_stop.setEnabled(False)
        self.btn_stop.clicked.connect(stop_cb)
        layout.addWidget(self.btn_stop)

        layout.addStretch()

        self.btn_log = QPushButton(log_label)
        self.btn_log.setStyleSheet(_LOG_SS)
        self.btn_log.clicked.connect(log_cb)
        layout.addWidget(self.btn_log)

    def set_running(self, running: bool) -> None:
        """Toggle run/stop enabled states."""
        self.btn_run.setEnabled(not running)
        self.btn_stop.setEnabled(running)
        if running:
            self.btn_stop.setText("Stop")

    def set_stopping(self) -> None:
        """Intermediate state while waiting for the worker to cancel."""
        self.btn_stop.setEnabled(False)
        self.btn_stop.setText("Stopping…")
