"""Common builder for step parameter dialogs.

Usage (single-section flat dialog):
    run_param_dialog(self, "Title", MY_SPECS, on_save=lambda: self.persist_params())

Usage (multi-section with QGroupBox):
    form, widgets = build_param_form(params.P, MY_SPECS)
    layout.addLayout(form)
    ...
    read_param_form(widgets, params.P, MY_SPECS)
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Sequence

from PyQt5.QtWidgets import (
    QCheckBox,
    QDoubleSpinBox,
    QFormLayout,
    QLabel,
    QMessageBox,
    QSpinBox,
)

from apex.gui.workflow.ui_helpers import add_parameter_reset_button, build_scroll_param_dialog


@dataclass
class ParamSpec:
    """Descriptor for one form row in a parameter dialog.

    kind: "float" | "int" | "bool" | "sep"
          "sep" inserts a full-width QLabel (attr is ignored).
    write_also: additional attrs to receive the same value on save.
    """

    label: str
    attr: str = ""
    kind: str = "float"
    lo: float = 0.0
    hi: float = 100.0
    step: float = 1.0
    decimals: int = 2
    suffix: str = ""
    tooltip: str = ""
    default: Any = None
    write_also: tuple[str, ...] = ()


def build_param_form(
    params_P: Any,
    specs: Sequence[ParamSpec],
    overrides: dict[str, Any] | None = None,
) -> tuple[QFormLayout, dict[str, Any]]:
    """Build a QFormLayout populated from params_P.

    Returns (form_layout, widgets) mapping attr → widget.
    overrides maps attr → initial value, bypassing getattr on params_P.
    """
    form = QFormLayout()
    widgets: dict[str, Any] = {}

    for spec in specs:
        if spec.kind == "sep":
            form.addRow(QLabel(spec.label))
            continue

        if overrides and spec.attr in overrides:
            raw = overrides[spec.attr]
        else:
            raw = getattr(params_P, spec.attr, spec.default)

        if spec.kind == "float":
            w = QDoubleSpinBox()
            w.setRange(spec.lo, spec.hi)
            w.setSingleStep(spec.step)
            w.setDecimals(spec.decimals)
            if spec.suffix:
                w.setSuffix(spec.suffix)
            w.setValue(float(raw) if raw is not None else spec.lo)
            if spec.tooltip:
                w.setToolTip(spec.tooltip)
            form.addRow(spec.label + ":", w)

        elif spec.kind == "int":
            w = QSpinBox()
            w.setRange(int(spec.lo), int(spec.hi))
            w.setSingleStep(max(1, int(spec.step)))
            if spec.suffix:
                w.setSuffix(spec.suffix)
            w.setValue(int(raw) if raw is not None else int(spec.lo))
            if spec.tooltip:
                w.setToolTip(spec.tooltip)
            form.addRow(spec.label + ":", w)

        elif spec.kind == "bool":
            w = QCheckBox("Enable")
            w.setChecked(bool(raw) if raw is not None else False)
            if spec.tooltip:
                w.setToolTip(spec.tooltip)
            form.addRow(spec.label + ":", w)

        else:
            continue

        widgets[spec.attr] = w

    return form, widgets


def read_param_form(
    widgets: dict[str, Any],
    params_P: Any,
    specs: Sequence[ParamSpec],
) -> None:
    """Write current widget values back to params_P."""
    for spec in specs:
        if spec.kind == "sep" or spec.attr not in widgets:
            continue
        w = widgets[spec.attr]
        if spec.kind in ("float", "int"):
            value = w.value()
        elif spec.kind == "bool":
            value = w.isChecked()
        else:
            continue
        setattr(params_P, spec.attr, value)
        for extra in spec.write_also:
            setattr(params_P, extra, value)


def _spec_default_value(spec: ParamSpec) -> Any:
    if spec.default is not None:
        return spec.default
    if spec.kind == "bool":
        return False
    if spec.kind == "int":
        return int(spec.lo)
    if spec.kind == "float":
        return float(spec.lo)
    return None


def run_param_dialog(
    parent: Any,
    title: str,
    specs: Sequence[ParamSpec],
    *,
    on_save: Callable[[], bool | None],
    overrides: dict[str, Any] | None = None,
    resize: tuple[int, int] = (460, 400),
    info_text: str = "",
) -> None:
    """Build and exec a single-section modal parameter dialog.

    on_save() is called after widget values are written to params.P.
    Return True/False to show a save-result message; return None to skip it.
    """
    params_P = parent.params.P
    dialog, layout, buttons = build_scroll_param_dialog(
        parent, title, info_text=info_text, size=resize
    )

    form, widgets = build_param_form(params_P, specs, overrides)
    layout.addLayout(form)
    layout.addStretch(1)

    add_parameter_reset_button(
        buttons,
        [
            (widgets[spec.attr], _spec_default_value(spec))
            for spec in specs
            if spec.kind != "sep" and spec.attr in widgets
        ],
    )

    def _save() -> None:
        read_param_form(widgets, params_P, specs)
        result = on_save()
        if result is None:
            dialog.accept()
            return
        msg = "Parameters saved." if result else "Parameters updated, but TOML save failed."
        QMessageBox.information(dialog, "Saved", msg)
        dialog.accept()

    buttons.accepted.connect(_save)
    buttons.rejected.connect(dialog.reject)
    dialog.exec_()
