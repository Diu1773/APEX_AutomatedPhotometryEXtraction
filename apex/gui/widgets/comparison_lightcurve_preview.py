"""Compact non-modal light-curve preview for comparison-star review."""

from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (
    QButtonGroup,
    QDialog,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
)

from apex.gui.theme import Tokens, style_button


class ComparisonLightCurvePreview(QDialog):
    use_as_comparison = pyqtSignal(int)
    reject_source = pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__(parent, Qt.Tool)
        self.setWindowTitle("Comparison Light Curve")
        self.setAttribute(Qt.WA_DeleteOnClose, False)
        self.resize(920, 310)
        self.setMinimumSize(720, 260)
        self._source_id: int | None = None
        self._candidate = pd.DataFrame()
        self._loo = pd.DataFrame()
        self._target_delta = pd.DataFrame()
        self._metrics: dict = {}

        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 8, 10, 8)
        layout.setSpacing(6)

        header = QHBoxLayout()
        self.title_label = QLabel("Select a source and press P")
        self.title_label.setProperty("role", "subtitle")
        header.addWidget(self.title_label)
        header.addStretch()

        self.mode_group = QButtonGroup(self)
        self.mode_group.setExclusive(True)
        self.mode_buttons: dict[str, QPushButton] = {}
        for mode, text in (
            ("loo", "LOO Residual"),
            ("mag", "Inst. Mag"),
            ("target", "Target Delta"),
        ):
            button = QPushButton(text)
            button.setCheckable(True)
            button.setFixedHeight(Tokens.H_COMPACT)
            self.mode_group.addButton(button)
            self.mode_buttons[mode] = button
            header.addWidget(button)
            button.clicked.connect(self._draw)
        self.mode_buttons["loo"].setChecked(True)
        layout.addLayout(header)

        self.metrics_label = QLabel("")
        self.metrics_label.setProperty("role", "caption")
        self.metrics_label.setWordWrap(True)
        layout.addWidget(self.metrics_label)

        self.figure = Figure(figsize=(9.0, 2.2), constrained_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.axes = self.figure.add_subplot(111)
        self.canvas.setMinimumHeight(155)
        layout.addWidget(self.canvas, 1)

        actions = QHBoxLayout()
        self.status_label = QLabel("")
        self.status_label.setProperty("role", "caption")
        actions.addWidget(self.status_label)
        actions.addStretch()
        self.reject_button = QPushButton("Reject")
        style_button(self.reject_button, "danger", height=Tokens.H_COMPACT)
        self.reject_button.clicked.connect(self._emit_reject)
        actions.addWidget(self.reject_button)
        self.use_button = QPushButton("Use as Comp")
        style_button(self.use_button, "primary", height=Tokens.H_COMPACT)
        self.use_button.clicked.connect(self._emit_use)
        actions.addWidget(self.use_button)
        layout.addLayout(actions)
        self._set_actions_enabled(False)

    def _set_actions_enabled(self, enabled: bool) -> None:
        self.reject_button.setEnabled(enabled)
        self.use_button.setEnabled(enabled)

    def set_loading(self, source_id: int, display_id: int | None, filter_name: str) -> None:
        self._source_id = int(source_id)
        display = display_id if display_id is not None else source_id
        self.title_label.setText(f"ID {display} | source {source_id} | {filter_name}")
        self.metrics_label.setText("Loading time series...")
        self.status_label.setText("")
        self._candidate = pd.DataFrame()
        self._loo = pd.DataFrame()
        self._target_delta = pd.DataFrame()
        self._set_actions_enabled(False)
        self.axes.clear()
        self.axes.text(0.5, 0.5, "Loading...", ha="center", va="center", transform=self.axes.transAxes)
        self.canvas.draw_idle()

    def set_payload(
        self,
        *,
        source_id: int,
        display_id: int | None,
        filter_name: str,
        photometry_source: str,
        candidate: pd.DataFrame,
        loo: pd.DataFrame,
        target_delta: pd.DataFrame,
        metrics: dict | None,
        is_comparison: bool,
        is_rejected: bool = False,
        can_assign: bool = True,
    ) -> None:
        self._source_id = int(source_id)
        self._candidate = candidate.copy()
        self._loo = loo.copy()
        self._target_delta = target_delta.copy()
        self._metrics = dict(metrics or {})
        display = display_id if display_id is not None else source_id
        self.title_label.setText(
            f"ID {display} | source {source_id} | {filter_name} | {photometry_source.upper()}"
        )
        self.use_button.setText("Remove Comp" if is_comparison else "Use as Comp")
        if is_comparison:
            self.status_label.setText("Current comparison")
        elif is_rejected:
            self.status_label.setText("Manually rejected")
        else:
            self.status_label.setText("")
        self._set_actions_enabled(can_assign)
        self.reject_button.setEnabled(can_assign and not is_rejected)
        self._update_metrics_label()
        self._draw()

    def set_error(self, message: str) -> None:
        self.metrics_label.setText(str(message))
        self._set_actions_enabled(False)
        self.axes.clear()
        self.axes.text(0.5, 0.5, str(message), ha="center", va="center", transform=self.axes.transAxes)
        self.canvas.draw_idle()

    def _selected_mode(self) -> str:
        for mode, button in self.mode_buttons.items():
            if button.isChecked():
                return mode
        return "loo"

    def _update_metrics_label(self) -> None:
        def _fmt(key: str, digits: int = 4, suffix: str = "") -> str:
            try:
                value = float(self._metrics.get(key, np.nan))
            except (TypeError, ValueError):
                value = np.nan
            return f"{value:.{digits}f}{suffix}" if np.isfinite(value) else "-"

        try:
            coverage = 100.0 * float(self._metrics.get("coverage", np.nan))
        except (TypeError, ValueError):
            coverage = np.nan
        coverage_text = f"{coverage:.1f}%" if np.isfinite(coverage) else "-"

        status = str(self._metrics.get("status", "unscored"))
        reason = str(self._metrics.get("reasons", "") or "")
        self.metrics_label.setText(
            f"{status.upper()} | N {_fmt('n', 0)} | coverage {coverage_text} | "
            f"RMS {_fmt('rms')} | robust sigma {_fmt('robust_sigma')} | "
            f"eta {_fmt('eta', 3)} | night sigma {_fmt('night_scatter')} | {reason}"
        )

    @staticmethod
    def _plot_coordinates(data: pd.DataFrame, value_column: str) -> tuple[np.ndarray, np.ndarray, str]:
        values = pd.to_numeric(data.get(value_column), errors="coerce").to_numpy(float)
        if "time" in data.columns:
            times = pd.to_numeric(data["time"], errors="coerce").to_numpy(float)
        else:
            times = np.full(len(data), np.nan)
        finite_time = times[np.isfinite(times)]
        if finite_time.size >= 2 and float(np.nanmax(finite_time) - np.nanmin(finite_time)) > 0:
            origin = float(np.floor(np.nanmin(finite_time)))
            span = float(np.nanmax(finite_time) - np.nanmin(finite_time))
            if span <= 2.0:
                x = (times - float(np.nanmin(finite_time))) * 24.0
                label = "Hours from first exposure"
            else:
                x = times - origin
                label = f"JD - {origin:.0f}"
        else:
            x = np.arange(len(data), dtype=float)
            label = "Frame order"
        return x, values, label

    def _draw(self) -> None:
        mode = self._selected_mode()
        if mode == "loo":
            data, value_column, y_label = self._loo, "residual", "LOO residual (mag)"
        elif mode == "target":
            data, value_column, y_label = self._target_delta, "value", "Target - candidate (mag)"
        else:
            data, value_column, y_label = self._candidate, "mag", "Instrumental magnitude"

        self.axes.clear()
        if data is None or data.empty or value_column not in data.columns:
            self.axes.text(0.5, 0.5, "No usable series", ha="center", va="center", transform=self.axes.transAxes)
            self.canvas.draw_idle()
            return
        x, y, x_label = self._plot_coordinates(data, value_column)
        if "night_id" in data.columns:
            night_values = pd.to_numeric(data["night_id"], errors="coerce").fillna(0).astype(int).to_numpy()
        else:
            night_values = np.zeros(len(data), dtype=int)
        unique_nights = sorted(set(int(value) for value in night_values))
        colors = ["#3A66DB", "#C44949", "#2F8A5B", "#9A6A22", "#6F5AA8", "#287D8E"]
        for color_index, night_id in enumerate(unique_nights):
            mask = (night_values == night_id) & np.isfinite(x) & np.isfinite(y)
            if not np.any(mask):
                continue
            label = f"Night {night_id}" if night_id > 0 and len(unique_nights) > 1 else None
            self.axes.scatter(
                x[mask], y[mask], s=15, color=colors[color_index % len(colors)],
                alpha=0.85, linewidths=0, label=label,
            )
        if mode in {"loo", "target"}:
            self.axes.axhline(0.0, color="#98A2B3", linewidth=0.8, linestyle="--")
        self.axes.set_xlabel(x_label, fontsize=8)
        self.axes.set_ylabel(y_label, fontsize=8)
        self.axes.tick_params(labelsize=8)
        self.axes.grid(True, alpha=0.18, linewidth=0.6)
        self.axes.invert_yaxis()
        if len(unique_nights) > 1:
            self.axes.legend(loc="best", fontsize=7, frameon=False, ncol=min(4, len(unique_nights)))
        self.canvas.draw_idle()

    def _emit_use(self) -> None:
        if self._source_id is not None:
            self.use_as_comparison.emit(int(self._source_id))

    def _emit_reject(self) -> None:
        if self._source_id is not None:
            self.reject_source.emit(int(self._source_id))
