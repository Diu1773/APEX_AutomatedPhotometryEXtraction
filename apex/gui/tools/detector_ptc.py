"""Detector characterisation tool — measure gain and read noise from calibration frames.

The gain and read noise the error model needs are not printed on the camera's
spec sheet: manufacturers quote a read-noise ceiling and a full-well figure, but
the conversion gain depends on the readout mode.  Worse, the FITS ``EGAIN``
keyword is often the sensor's nominal register value rather than the realised
gain — on the camera used for APEX's own validation it is wrong by a factor of
14, which would scale every reported magnitude error with it.

This tool measures both from the user's own bias and flat frames using the
photon-transfer method, compares the result against the header, and can write
the measured values straight into the workspace configuration.

Runs the measurement on a worker thread; only the main thread touches widgets.
"""

from __future__ import annotations

import glob
from pathlib import Path

from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtWidgets import (
    QCheckBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox,
    QHBoxLayout, QLabel, QLineEdit, QMessageBox, QProgressBar, QPushButton,
    QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from apex.gui.theme import ICON, Tokens, style_button
from apex.gui.tools.tool_window_base import ToolWindowBase

_BIAS_PATTERNS = ("bias*.fit", "bias*.fits", "Bias*.fit", "BIAS*.fit")
_FLAT_PATTERNS = ("*.fit", "*.fits")


def _gather(folder: str, patterns) -> list[str]:
    found: set[str] = set()
    for pattern in patterns:
        found.update(glob.glob(str(Path(folder) / "**" / pattern), recursive=True))
    return sorted(found)


class _PTCWorker(QThread):
    """Runs the measurement off the GUI thread."""

    progressed = pyqtSignal(str, int, int)
    logged = pyqtSignal(str)
    finished_ok = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, bias_paths, flat_paths, *, signal_floor, binning, parent=None):
        super().__init__(parent)
        self._bias = list(bias_paths)
        self._flats = list(flat_paths)
        self._floor = float(signal_floor)
        self._binning = int(binning) if binning else None

    def run(self):  # noqa: D102 — QThread entry point
        try:
            from apex.analysis.detector_ptc import characterize_detector
            result = characterize_detector(
                self._bias, self._flats,
                binning=self._binning,
                signal_floor=self._floor,
                fix_intercept=True,
                progress=lambda stage, i, n: self.progressed.emit(stage, i, n),
                log_fn=self.logged.emit,
            )
        except Exception as exc:  # noqa: BLE001 — surfaced in the UI
            self.failed.emit(str(exc))
            return
        self.finished_ok.emit(result)


class DetectorPTCWindow(ToolWindowBase):
    """Measure gain / read noise from bias and flat frames."""

    def __init__(self, params=None, project_state=None, parent=None):
        super().__init__("Detector Characterisation (PTC)", params=params,
                         project_state=project_state, parent=parent,
                         min_size=(900, 680))
        self._worker: _PTCWorker | None = None
        self._result = None

        self.add_log_dock(popup=False)
        self.content_layout.addWidget(self._build_inputs())
        self.content_layout.addWidget(self._build_run_bar())
        self.content_layout.addWidget(self._build_results(), stretch=1)
        self._prefill_from_params()

    # ── layout ───────────────────────────────────────────────────────────
    def _build_inputs(self) -> QWidget:
        box = QGroupBox("Calibration frames")
        form = QFormLayout(box)
        form.setSpacing(Tokens.GAP)

        self.bias_edit = QLineEdit()
        self.flat_edit = QLineEdit()
        form.addRow("Bias folder:", self._with_browse(self.bias_edit, "bias"))
        form.addRow("Flat folder:", self._with_browse(self.flat_edit, "flat"))

        self.floor_spin = QDoubleSpinBox()
        self.floor_spin.setRange(0.0, 60000.0)
        self.floor_spin.setSingleStep(1000.0)
        self.floor_spin.setValue(20000.0)
        self.floor_spin.setSuffix(" ADU")
        self.floor_spin.setToolTip(
            "Flats below this signal are ignored. Keep some spread in exposure "
            "level — a narrow range leaves the slope poorly constrained."
        )
        form.addRow("Minimum flat signal:", self.floor_spin)

        self.apply_check = QCheckBox("Write the measured values into the workspace config")
        self.apply_check.setChecked(True)
        form.addRow("", self.apply_check)
        return box

    def _with_browse(self, edit: QLineEdit, kind: str) -> QWidget:
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(Tokens.GAP)
        layout.addWidget(edit, stretch=1)
        button = QPushButton(f"{ICON['input']} Browse")
        style_button(button, height=Tokens.H_BUTTON)
        button.clicked.connect(lambda: self._pick_folder(edit, kind))
        layout.addWidget(button)
        return row

    def _build_run_bar(self) -> QWidget:
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(Tokens.GAP)

        self.run_button = QPushButton("Measure")
        style_button(self.run_button, "primary", height=Tokens.H_ACTION)
        self.run_button.clicked.connect(self._start)
        layout.addWidget(self.run_button)

        self.progress = QProgressBar()
        self.progress.setTextVisible(True)
        self.progress.setFormat("%p%")
        layout.addWidget(self.progress, stretch=1)
        return row

    def _build_results(self) -> QWidget:
        box = QGroupBox("Measured detector constants")
        layout = QVBoxLayout(box)

        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["Quantity", "Value"])
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setVisible(False)
        layout.addWidget(self.table, stretch=1)

        self.verdict = QLabel("")
        self.verdict.setWordWrap(True)
        self.verdict.setTextFormat(Qt.TextFormat.RichText)
        layout.addWidget(self.verdict)
        return box

    # ── behaviour ────────────────────────────────────────────────────────
    def _prefill_from_params(self) -> None:
        data_dir = getattr(getattr(self.params, "P", None), "data_dir", None)
        if data_dir:
            self.bias_edit.setText(str(data_dir))
            self.flat_edit.setText(str(data_dir))

    def _pick_folder(self, edit: QLineEdit, kind: str) -> None:
        folder = QFileDialog.getExistingDirectory(
            self, f"Select the {kind} folder", edit.text() or "")
        if folder:
            edit.setText(folder)

    def _start(self) -> None:
        bias = _gather(self.bias_edit.text(), _BIAS_PATTERNS)
        flats = _gather(self.flat_edit.text(), _FLAT_PATTERNS)
        if len(bias) < 2:
            QMessageBox.warning(self, "Not enough bias frames",
                                "At least two bias frames are needed to measure "
                                "the read noise.")
            return
        if len(flats) < 6:
            QMessageBox.warning(self, "Not enough flat frames",
                                "At least six flats are needed, ideally spread "
                                "over several exposure levels.")
            return

        self.append_log(f"bias {len(bias)} frames, flat candidates {len(flats)}")
        self.run_button.setEnabled(False)
        self.progress.setValue(0)

        binning = getattr(getattr(self.params, "P", None), "binning", None)
        self._worker = _PTCWorker(bias, flats,
                                  signal_floor=self.floor_spin.value(),
                                  binning=binning, parent=self)
        self._worker.progressed.connect(self._on_progress)
        self._worker.logged.connect(self.append_log)
        self._worker.finished_ok.connect(self._on_done)
        self._worker.failed.connect(self._on_failed)
        self._worker.start()

    def _on_progress(self, stage: str, done: int, total: int) -> None:
        self.progress.setFormat(f"{stage} %p%")
        self.progress.setValue(int(100 * done / total) if total else 0)

    def _on_failed(self, message: str) -> None:
        self.run_button.setEnabled(True)
        self.progress.setValue(0)
        self.append_log(f"measurement failed: {message}")
        QMessageBox.critical(self, "Measurement failed", message)

    def _on_done(self, result) -> None:
        self._result = result
        self.run_button.setEnabled(True)
        self.progress.setValue(100)
        self._fill_table(result)
        self._fill_verdict(result)
        if self.apply_check.isChecked():
            self._write_config(result)

    def _fill_table(self, r) -> None:
        rows = [
            ("Gain (stored pixel)", f"{r.gain_eff:.4f} ± {r.gain_eff_err:.4f} e⁻/ADU"),
            ("Read noise (stored pixel)", f"{r.read_noise_eff:.3f} e⁻"),
            ("Gain (per photosite)", f"{r.gain_pixel:.4f} e⁻/ADU"),
            ("Read noise (per photosite)", f"{r.read_noise_pixel:.3f} e⁻"),
            ("Binning", f"{r.binning}×{r.binning} ({r.binning_mode})"),
            ("Flat pairs used", f"{r.n_pairs}"),
            ("Signal range", f"{r.signal_min:,.0f} – {r.signal_max:,.0f} ADU"),
            ("Lever arm", f"{100 * r.lever_arm:.0f} % of the top signal"),
            ("Fit R²", f"{r.r_squared:.5f}"),
            ("Largest residual", f"{100 * r.max_residual_frac:.2f} % of variance"),
        ]
        if r.header_egain:
            rows.append(("Header EGAIN", f"{r.header_egain:.5g} e⁻/ADU"))
            rows.append(("Measured / header", f"{r.header_ratio:.2f}×"))

        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(rows))
        for i, (name, value) in enumerate(rows):
            self.table.setItem(i, 0, QTableWidgetItem(name))
            self.table.setItem(i, 1, QTableWidgetItem(value))
        self.table.resizeColumnsToContents()

    def _fill_verdict(self, r) -> None:
        notes: list[str] = []
        if r.header_ratio is not None and not (0.9 <= r.header_ratio <= 1.1):
            notes.append(
                f"The FITS header reports {r.header_egain:.4g} e⁻/ADU, "
                f"<b>{r.header_ratio:.1f}× away</b> from the measured gain. "
                "Leave <code>noise_use_fits_header</code> off so the error model "
                "uses the measured value."
            )
        if r.message and r.message != "ok":
            notes.append(f"<b>Caution:</b> {r.message}.")
        if r.max_residual_frac > 0.15:
            notes.append(
                "Some points sit well off the fitted line — check for saturated "
                "flats or frames taken under a changing sky."
            )
        if not notes:
            notes.append("The photon-transfer fit looks well behaved.")
        self.verdict.setText("<br><br>".join(notes))

    def _write_config(self, r) -> None:
        target = getattr(self.params, "P", None)
        if target is None:
            self.append_log("no workspace config to write to — values not applied")
            return
        try:
            target.gain_e_per_adu = float(r.gain_eff)
            target.rdnoise_e = float(r.read_noise_eff)
            target.noise_use_fits_header = False
            save = getattr(self.params, "save", None)
            if callable(save):
                save()
            self.append_log(
                f"config updated: gain_e_per_adu = {r.gain_eff:.4f}, "
                f"rdnoise_e = {r.read_noise_eff:.3f}, noise_use_fits_header = False"
            )
        except Exception as exc:  # noqa: BLE001
            self.append_log(f"could not write the config: {exc}")

    def closeEvent(self, event):  # noqa: N802 — Qt naming
        if self._worker is not None and self._worker.isRunning():
            self._worker.wait(3000)
        super().closeEvent(event)
