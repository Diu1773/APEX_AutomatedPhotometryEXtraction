"""
Step 11: CMD Plot
Loads calibrated CMD tables and opens the CMD viewer.
"""

from __future__ import annotations

import time
from pathlib import Path

import pandas as pd

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QGroupBox, QTextEdit, QMessageBox,
    QScrollArea, QFrame, QPushButton, QWidget,
)

from apex.gui.workflow.step_window_base import StepWindowBase
from .step10_zeropoint_calibration import (
    CmdViewerWindow,
    resolve_cmd_photometry_provenance,
)
from apex.utils.step_paths import step7_forced_phot_dir
from apex.utils.step_paths_cmd import step8_psf_dir, step10_zp_dir, step11_cmd_dir
from apex.utils.photometry_provenance import (
    format_photometry_provenance,
    summarize_photometry_table,
)


class CmdPlotWindow(StepWindowBase):
    """Step 11: CMD Plot"""

    def __init__(self, params, file_manager, project_state, main_window):
        self.file_manager = file_manager
        self.viewer = None
        self._viewer_scroll = None
        super().__init__(
            step_index=10,
            step_name="CMD Plot",
            params=params,
            project_state=project_state,
            main_window=main_window
        )
        self.setup_step_ui()
        self.restore_state()

    def setup_step_ui(self):
        info = QLabel("Open CMD viewer using calibrated ZP products.")
        info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 10px; border-radius: 5px; }")
        self.content_layout.addWidget(info)

        self.photometry_source_label = QLabel(
            format_photometry_provenance(
                resolve_cmd_photometry_provenance(self.params.P.result_dir)
            )
        )
        self.photometry_source_label.setStyleSheet(
            "QLabel { color: #455A64; font-weight: bold; padding: 2px 4px; }"
        )
        self.content_layout.addWidget(self.photometry_source_label)

        self.viewer_group = QGroupBox("CMD Viewer")
        self.viewer_layout = QVBoxLayout(self.viewer_group)
        self.viewer_placeholder = QLabel("CMD viewer will appear here after ZP calibration.")
        self.viewer_placeholder.setAlignment(Qt.AlignCenter)
        self.viewer_placeholder.setStyleSheet("QLabel { color: #666666; padding: 12px; }")
        self.viewer_layout.addWidget(self.viewer_placeholder)
        # stretch=1 so the viewer takes priority over the log when the
        # window resizes — otherwise the QTextEdit greedily grabs vertical
        # space and squeezes the CMD plot to the top.
        self.content_layout.addWidget(self.viewer_group, stretch=1)

        # Log lives in a popup window so it doesn't waste vertical space
        # below the CMD plot.  The main page just exposes a "Show Log"
        # button.
        self.log_window = QWidget(self, Qt.Window)
        self.log_window.setWindowTitle("Step 11 — CMD Plot Log")
        self.log_window.resize(700, 360)
        _log_layout = QVBoxLayout(self.log_window)
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setStyleSheet("QTextEdit { font-family: monospace; font-size: 9pt; }")
        _log_layout.addWidget(self.log_text)

        log_row = QHBoxLayout()
        log_row.addStretch()
        self.btn_show_log = QPushButton("Show Log")
        self.btn_show_log.setStyleSheet(
            "QPushButton { background-color: #607D8B; color: white; "
            "font-weight: bold; padding: 6px 14px; }"
            "QPushButton:hover { background-color: #455A64; }"
        )
        self.btn_show_log.clicked.connect(self._show_log_window)
        log_row.addWidget(self.btn_show_log)
        self.content_layout.addLayout(log_row)

    def _show_log_window(self) -> None:
        self.log_window.show()
        self.log_window.raise_()
        self.log_window.activateWindow()

    def log(self, message: str):
        timestamp = time.strftime("%H:%M:%S")
        self.log_text.append(f"[{timestamp}] {message}")

    def open_cmd_viewer(self):
        input_dir = step10_zp_dir(self.params.P.result_dir)
        if not input_dir.exists():
            input_dir = self.params.P.result_dir
        output_dir = step11_cmd_dir(self.params.P.result_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        wide_path = input_dir / "median_by_ID_filter_wide_cmd.csv"
        if not wide_path.exists():
            wide_path = input_dir / "median_by_ID_filter_wide.csv"
        if not wide_path.exists():
            self._reset_viewer("CMD wide CSV not found")
            QMessageBox.warning(self, "Missing Data", "CMD wide CSV not found")
            return

        idx_candidates = [
            step7_forced_phot_dir(self.params.P.result_dir) / "photometry_index.csv",
            step8_psf_dir(self.params.P.result_dir) / "photometry_index.csv",
            self.params.P.result_dir / "photometry_index.csv",
            self.params.P.result_dir / "phot_index.csv",
            self.params.P.result_dir / "phot" / "photometry_index.csv",
            self.params.P.result_dir / "phot" / "phot_index.csv",
        ]
        idx_path = next((p for p in idx_candidates if p.exists()), None)
        if idx_path and wide_path.stat().st_mtime < idx_path.stat().st_mtime:
            self.log(f"WARNING: {wide_path.name} is older than {idx_path.name}; rerun ZP calibration.")

        df = pd.read_csv(wide_path, dtype={"source_id": str, "gaia_source_id": str})
        self.photometry_source_label.setText(
            format_photometry_provenance(summarize_photometry_table(df))
        )
        self._reset_viewer()
        viewer = CmdViewerWindow(df, self.params.P.result_dir, self, embedded=True, params=self.params)
        # Wrap in a scroll area so cramped windows can still reach the
        # Parallax/ROI filter row + Prev/Next View buttons; otherwise the
        # 800×400-min canvas pushes those controls off-screen.
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll.setWidget(viewer)
        self.viewer_layout.addWidget(scroll)
        self.viewer = viewer
        self._viewer_scroll = scroll
        self.viewer_placeholder.setVisible(False)

    def _reset_viewer(self, placeholder_text: str | None = None):
        if self.viewer is not None:
            self.viewer.setParent(None)
            self.viewer.deleteLater()
            self.viewer = None
        scroll = getattr(self, "_viewer_scroll", None)
        if scroll is not None:
            scroll.setParent(None)
            scroll.deleteLater()
            self._viewer_scroll = None
        if placeholder_text:
            self.viewer_placeholder.setText(placeholder_text)
        self.viewer_placeholder.setVisible(self.viewer is None)

    def validate_step(self) -> bool:
        input_dir = step10_zp_dir(self.params.P.result_dir)
        if not input_dir.exists():
            input_dir = self.params.P.result_dir
        return (input_dir / "median_by_ID_filter_wide_cmd.csv").exists() or (input_dir / "median_by_ID_filter_wide.csv").exists()

    def save_state(self):
        self.project_state.store_step_data("cmd_plot", {})

    def restore_state(self):
        if self.validate_step():
            try:
                self.open_cmd_viewer()
            except Exception:
                pass
        return
