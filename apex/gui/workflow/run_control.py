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

import math
import time

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


def format_duration(seconds: float) -> str:
    """Format elapsed/remaining seconds as MM:SS or H:MM:SS."""
    try:
        value = float(seconds)
    except Exception:
        value = 0.0
    if not math.isfinite(value):
        value = 0.0
    total = int(max(0, round(value)))
    hours, rem = divmod(total, 3600)
    minutes, secs = divmod(rem, 60)
    if hours > 0:
        return f"{hours:d}:{minutes:02d}:{secs:02d}"
    return f"{minutes:02d}:{secs:02d}"


def progress_status_text(
    current: int,
    total: int,
    start_time: float | None,
    *,
    message: str = "",
    workers: int | None = None,
    now: float | None = None,
) -> str:
    """Build a consistent progress label with elapsed time and ETA."""
    current_i = int(max(0, current))
    total_i = int(max(0, total))
    parts = [f"{current_i}/{total_i}"]

    if start_time is not None:
        clock_now = time.monotonic() if now is None else float(now)
        elapsed = max(0.0, clock_now - float(start_time))
        parts.append(f"elapsed {format_duration(elapsed)}")
        if current_i > 0 and total_i > current_i:
            remaining = elapsed / float(current_i) * float(total_i - current_i)
            parts.append(f"ETA {format_duration(remaining)}")
        else:
            parts.append("ETA --:--")

    if workers is not None:
        parts.append(f"W:{int(workers)}")
    if message:
        parts.append(str(message))
    return " | ".join(parts)


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
