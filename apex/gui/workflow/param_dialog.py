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
from functools import lru_cache
from typing import Any, Callable, Sequence

from PyQt5.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QLabel,
    QMessageBox,
    QSpinBox,
)

from apex.config.parameter_map import toml_key_map_for_mode
from apex.gui.layout_rules import fit_combo
from apex.gui.workflow.ui_helpers import add_parameter_reset_button, build_scroll_param_dialog


@dataclass
class ParamSpec:
    """Descriptor for one form row in a parameter dialog.

    kind: "float" | "int" | "bool" | "choice" | "sep"
          "sep" inserts a full-width QLabel (attr is ignored).
          "choice" is a combo box; see `choices`.
    choices: for kind="choice", the allowed values in display order. Either
          plain strings, or (value, label) pairs when what is stored differs
          from what is shown. The stored value is what reaches the config, so
          it must be a string the loader and the calculation both understand —
          that is the whole reason these live on the map row and not in the
          window: a combo offering an option the loader does not accept is the
          same defect as a widget the config cannot save.
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
    choices: tuple = ()
    write_also: tuple[str, ...] = ()

    def choice_pairs(self) -> tuple[tuple[str, str], ...]:
        """`choices` as (stored value, shown label), however it was written."""
        pairs = []
        for entry in self.choices:
            if isinstance(entry, (tuple, list)) and len(entry) == 2:
                pairs.append((str(entry[0]), str(entry[1])))
            else:
                pairs.append((str(entry), str(entry)))
        return tuple(pairs)


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

        elif spec.kind == "choice":
            w = QComboBox()
            for value, shown in spec.choice_pairs():
                w.addItem(shown, value)
            # Select by stored value, not by position: the config holds the
            # value, and a list that gains an option must not silently change
            # what an existing workspace means.
            current = str(raw) if raw is not None else str(spec.default or "")
            index = w.findData(current)
            if index < 0 and current:
                # A workspace written before this option list existed. Keep the
                # value visible rather than snapping it to the first entry.
                w.addItem(f"{current} (설정 파일 값)", current)
                index = w.count() - 1
            w.setCurrentIndex(max(index, 0))
            fit_combo(w)
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
        elif spec.kind == "choice":
            # The stored value, not the shown label.
            value = w.currentData()
            if value is None:
                value = w.currentText()
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
    if spec.kind == "choice":
        pairs = spec.choice_pairs()
        return pairs[0][0] if pairs else None
    return None


@lru_cache(maxsize=4)
def _unread_for(mode: str) -> frozenset[str]:
    """Names no module reads, computed once per mode.

    Cached because it walks the source tree; a dialog opening must not pay for
    that twice. If the walk fails for any reason the guard opens rather than
    breaking the dialog — a missing check is better than an unusable window.
    """
    try:
        from apex.config.config_audit import unread_settings
        return frozenset(unread_settings(mode=mode)["dead"])
    except Exception:                                    # noqa: BLE001
        return frozenset()


def specs_from_map(
    attrs: Sequence[str],
    *,
    mode: str = "cmd",
    overrides: dict[str, dict[str, Any]] | None = None,
) -> tuple[ParamSpec, ...]:
    """Build the dialog rows from the key map instead of declaring them twice.

    A window used to repeat what the map already said — the attribute, its type
    and its default — and add the display bits on top. Two lists, no check that
    they agreed, and a third place (the loader) that had to agree with both.
    That is how a window came to offer settings it could not save: the widget
    existed, the map row did not, and `save_toml` skipped it in silence.

    Now the row is the declaration and the window only names which rows it
    shows, in order. A setting a window can display is one the loader builds and
    the map persists, by construction.

    `overrides` is for the rare row whose label or range is window-specific;
    anything passed here is a deliberate exception, not a second declaration.

    A row that no module reads is refused. Offering a knob that changes nothing
    is worse than not offering it: the user turns it, saves, sees the value in
    the file, and gets the same answer. Fifty-odd such rows exist (see
    `config_audit.unread_settings`); none of them should reach a dialog until
    someone either wires it up or deletes it.
    """
    rows = {}
    for row in toml_key_map_for_mode(mode):
        if len(row) >= 4 and row[1] not in rows:
            rows[row[1]] = row
        elif len(row) >= 5:
            rows[row[1]] = row

    specs: list[ParamSpec] = []
    for attr in attrs:
        if attr == "sep":
            specs.append(ParamSpec("", kind="sep"))
            continue
        row = rows.get(attr)
        if row is None:
            raise KeyError(
                f"{attr!r} 은 키 맵에 값이 없다 — 창에 띄우려면 먼저 행을 만들 것 "
                f"(그래야 저장이 파일까지 간다)"
            )
        if attr in _unread_for(mode):
            raise KeyError(
                f"{attr!r} 은 어떤 모듈도 읽지 않는다 — 창에 띄우면 사용자가 "
                f"돌려도 아무 일이 안 생긴다. 배선하거나 맵에서 지울 것 "
                f"(python -m apex.config.config_audit --unread)"
            )
        meta: dict[str, Any] = dict(row[4]) if len(row) >= 5 and row[4] else {}
        meta.update((overrides or {}).get(attr, {}))
        specs.append(ParamSpec(
            label=meta.pop("label", attr),
            attr=attr,
            kind=row[2],
            default=row[3],
            **meta,
        ))
    return tuple(specs)


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
