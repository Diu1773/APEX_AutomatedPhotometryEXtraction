"""Small shared UI helpers for workflow step controls."""
from __future__ import annotations

from PyQt5.QtCore import QEvent, QObject, Qt, QTimer
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (
    QAbstractSpinBox,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QGroupBox,
    QLabel,
    QLineEdit,
    QPushButton,
    QScrollArea,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)


PARAM_BUTTON_STYLE = """
QPushButton {
    background-color: #37474F;
    color: white;
    border: 1px solid #263238;
    border-radius: 4px;
    font-weight: bold;
    padding: 8px 14px;
}
QPushButton:hover {
    background-color: #455A64;
}
QPushButton:pressed {
    background-color: #263238;
}
QPushButton:disabled {
    background-color: #B0BEC5;
    color: #ECEFF1;
    border-color: #90A4AE;
}
"""


PARAM_DIALOG_STYLE = """
QDialog {
    background-color: #FAFAFA;
}
QGroupBox {
    font-weight: bold;
    border: 1px solid #D0D7DE;
    border-radius: 5px;
    margin-top: 10px;
    padding-top: 10px;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 8px;
    padding: 0 4px;
}
QLabel {
    color: #263238;
}
QDialogButtonBox QPushButton {
    min-width: 84px;
    padding: 6px 12px;
}
"""


_CACHE_CHECKBOX_STYLE = """
QCheckBox {
    color: #263238;
    font-weight: 600;
    padding: 2px 4px;
}
QCheckBox:disabled {
    color: #90A4AE;
}
"""


_CACHE_ACTION_STYLE = """
QPushButton {
    background-color: #455A64;
    color: white;
    border: 1px solid #37474F;
    border-radius: 4px;
    font-weight: bold;
    padding: 8px 12px;
}
QPushButton:hover {
    background-color: #546E7A;
}
QPushButton:pressed {
    background-color: #37474F;
}
QPushButton:disabled {
    background-color: #B0BEC5;
    color: #ECEFF1;
    border-color: #90A4AE;
}
"""


_COLLAPSE_BUTTON_STYLE = """
QPushButton {
    border: none;
    text-align: left;
    color: #455A64;
    font-size: 8pt;
    font-weight: 600;
    padding: 0px 0px 4px 0px;
}
QPushButton:checked {
    color: #1565C0;
}
"""


class CollapsibleSection(QGroupBox):
    """A lightweight collapsible section for parameter dialogs."""

    def __init__(
        self,
        title: str,
        *,
        initial_expanded: bool = False,
        parent: QWidget | None = None,
    ):
        super().__init__(title, parent)
        self.setCheckable(False)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(2)

        self.toggle_button = QPushButton(self)
        self.toggle_button.setCheckable(True)
        self.toggle_button.setStyleSheet(_COLLAPSE_BUTTON_STYLE)
        layout.addWidget(self.toggle_button)

        self.content_widget = QWidget(self)
        layout.addWidget(self.content_widget)

        self.toggle_button.toggled.connect(self._set_expanded)
        self.set_expanded(initial_expanded)

    def set_expanded(self, expanded: bool) -> None:
        self.toggle_button.setChecked(bool(expanded))
        self._set_expanded(bool(expanded))

    def _set_expanded(self, expanded: bool) -> None:
        self.content_widget.setVisible(bool(expanded))
        self.toggle_button.setText("Collapse" if expanded else "Expand")


def create_parameter_button(text: str = "Parameters") -> QPushButton:
    """Return a consistently styled workflow parameter button (neutral role)."""
    from apex.gui.theme import Tokens, style_button
    return style_button(QPushButton(text), height=Tokens.H_BUTTON)


def fit_parameter_dialog_width(dialog: QDialog) -> int:
    """Widen a parameter dialog until its scroll content fits, and report the gain.

    The `width` passed to `configure_parameter_dialog` is a hand-picked number,
    so a single longer label or wider combo silently pushes the content past the
    viewport and a horizontal scrollbar appears. That is a defect on its own —
    the reader has to drag sideways to see values (F-325 was the same family,
    fixed then by shrinking the content instead of by sizing the dialog).

    Measured on Step 8 (2026-09-15): content needed 611 px but the viewport was
    581 px, because the requested 620 px loses ~39 px to the frame and the
    vertical scrollbar. Sizing from the content removes that whole class.

    This only ever widens, never narrows, so dialogs that already fit are
    untouched. `clamp_to_screen` still caps it, and on a screen too small to fit
    the content the horizontal scrollbar correctly stays.
    """
    from apex.gui.layout_rules import clamp_to_screen

    # Measure only once the layout is current. Sections are added after the
    # dialog is constructed, so on the first pass the hints are still stale and
    # the fit silently does nothing — which is exactly what happened to Step 9
    # (no collapsible sections, so nothing later triggered a second pass).
    dialog.ensurePolished()
    if dialog.layout() is not None:
        dialog.layout().activate()

    needed = 0
    for area in dialog.findChildren(QScrollArea):
        content = area.widget()
        if content is None:
            continue
        content.ensurePolished()
        if content.layout() is not None:
            content.layout().activate()
        bar = area.verticalScrollBar()
        # Reserve the vertical scrollbar whenever the policy can produce one,
        # not only when it happens to be up right now. At the moment this runs
        # the dialog has just been laid out and the bar is often not visible
        # yet, so keying off `isVisible()` under-reserves by exactly its width —
        # measured on Step 8, that left 17 px of the original 30 px overflow.
        # Reserving a bar that never appears costs 17 px of width and nothing else.
        may_scroll = (bar is not None
                      and area.verticalScrollBarPolicy() != Qt.ScrollBarAlwaysOff)
        reserve = bar.sizeHint().width() if may_scroll else 0
        needed = max(needed,
                     content.minimumSizeHint().width() + reserve + 2 * area.frameWidth())
    if needed <= 0:
        return 0
    layout = dialog.layout()
    if layout is not None:
        m = layout.contentsMargins()
        needed += m.left() + m.right()
    if needed <= dialog.width():
        return 0
    before = dialog.width()
    dialog.resize(*clamp_to_screen(needed, dialog.height(), dialog))
    return dialog.width() - before


def configure_parameter_dialog(dialog: QDialog, title: str, width: int = 560, height: int = 620) -> None:
    """Apply common title, size, and visual style to a parameter dialog."""
    dialog.setWindowTitle(title)
    # Clamp to the monitor: the Save/Cancel row lives outside the scroll area,
    # so an over-tall dialog would clip those buttons on small screens.
    from apex.gui.layout_rules import clamp_to_screen
    dialog.resize(*clamp_to_screen(width, height, dialog))
    dialog.setStyleSheet(PARAM_DIALOG_STYLE)
    QTimer.singleShot(0, lambda: install_parameter_wheel_guard(dialog))
    # Deferred on purpose: sections are added after this call returns, so the
    # content width is only knowable once the event loop has laid the dialog out.
    QTimer.singleShot(0, lambda: install_parameter_dialog_autofit(dialog))


class _ShowRefit(QObject):
    """Re-fit the dialog width on its first show.

    Fitting before the dialog is shown reads stale content hints: measured on
    Step 9 (no collapsible sections), the content minimum was still small at
    that point, so the fit decided nothing was needed while the shown dialog
    then overflowed by 106 px. Step 8 only looked fixed because its sections
    fire the expand hook after the dialog is up.
    """

    def __init__(self, dialog: QDialog):
        super().__init__(dialog)
        self._dialog = dialog
        self._done = False

    def eventFilter(self, obj, event):  # noqa: N802 - Qt naming
        if not self._done and event.type() == QEvent.Show and obj is self._dialog:
            self._done = True
            # Next event-loop turn: on Show the layout has not run yet.
            QTimer.singleShot(0, lambda: fit_parameter_dialog_width(self._dialog))
        return False


def install_parameter_dialog_autofit(dialog: QDialog) -> None:
    """Fit the dialog to its content now, and again whenever a section expands.

    Fitting once is not enough. Sections open collapsed, so the first fit only
    sees the collapsed content, and the width a user actually needs appears the
    moment they expand one — which is exactly when they are trying to read the
    values. Re-fitting on expand keeps the horizontal scrollbar away for the
    state the reader is in.

    Widening only, never narrowing: collapsing a section leaves the dialog as
    wide as it was, which avoids the window jumping around while the reader
    opens and closes sections.
    """
    fit_parameter_dialog_width(dialog)
    # The dialog is usually not shown yet when this runs, and an unshown dialog
    # reports content hints that are too small — so fit again on first show.
    if getattr(dialog, "_apex_show_refit", None) is None:
        refit = _ShowRefit(dialog)
        setattr(dialog, "_apex_show_refit", refit)
        dialog.installEventFilter(refit)
    for section in dialog.findChildren(CollapsibleSection):
        button = getattr(section, "toggle_button", None)
        if button is None:
            continue
        # Queued so the re-fit runs after the layout has taken the newly shown
        # content into account; a direct call would measure the old width.
        button.toggled.connect(
            lambda _checked, d=dialog: QTimer.singleShot(
                0, lambda: fit_parameter_dialog_width(d)))


def create_collapsible_section(
    title: str, *, initial_expanded: bool = False
) -> tuple[CollapsibleSection, QWidget]:
    """Return a common collapsible parameter section and its content widget."""
    section = CollapsibleSection(title, initial_expanded=initial_expanded)
    return section, section.content_widget


def create_cache_checkbox(text: str, checked: bool, tooltip: str = "") -> QCheckBox:
    """Return a consistently styled cache/reuse checkbox."""
    checkbox = QCheckBox(text)
    checkbox.setChecked(bool(checked))
    checkbox.setStyleSheet(_CACHE_CHECKBOX_STYLE)
    if tooltip:
        checkbox.setToolTip(tooltip)
    return checkbox


def create_output_reuse_checkbox(checked: bool = True, tooltip: str = "") -> QCheckBox:
    """Return the standard step-output reuse checkbox."""
    return create_cache_checkbox(
        "Use existing output if complete",
        checked,
        tooltip or (
            "Load complete existing outputs when they match the current run. "
            "Disable to force recomputation."
        ),
    )


def create_detection_cache_checkbox(checked: bool = True, tooltip: str = "") -> QCheckBox:
    """Return the standard Step 4 detection-cache reuse checkbox."""
    return create_cache_checkbox(
        "Use detection cache",
        checked,
        tooltip or (
            "Skip frames with compatible Step 4 detection cache. "
            "Disable to force source detection for every selected frame."
        ),
    )


def create_cache_action_button(text: str) -> QPushButton:
    """Return a consistently styled cache maintenance button (neutral role)."""
    from apex.gui.theme import Tokens, style_button
    return style_button(QPushButton(text), height=Tokens.H_BUTTON)


_INFO_LABEL_STYLE = "QLabel { background-color: #E3F2FD; padding: 10px; margin-bottom: 10px; }"

STATUS_ROW_OK_BG = "#C8E6C9"
STATUS_ROW_WARN_BG = "#FFF3CD"
STATUS_ROW_FAIL_BG = "#FFCDD2"
STATUS_ROW_NEUTRAL_BG = "#ECEFF1"


def status_row_background(ok: bool | None, *, warning: bool = False) -> str:
    """Return the standard workflow table background for a row status."""
    if ok is None:
        return STATUS_ROW_NEUTRAL_BG
    if ok:
        return STATUS_ROW_WARN_BG if warning else STATUS_ROW_OK_BG
    return STATUS_ROW_FAIL_BG


def set_table_row_background(table: QTableWidget, row: int, color: str | QColor) -> None:
    """Apply one background color across an existing QTableWidget row."""
    bg = color if isinstance(color, QColor) else QColor(str(color))
    for col in range(table.columnCount()):
        item = table.item(row, col)
        if item is None:
            item = QTableWidgetItem("")
            table.setItem(row, col, item)
        item.setBackground(bg)


def set_parameter_widget_value(widget: QWidget, value) -> None:
    """Set a common parameter widget to *value* without knowing its concrete type."""
    if isinstance(widget, QCheckBox):
        widget.setChecked(bool(value))
        return
    if isinstance(widget, QComboBox):
        idx = widget.findData(value)
        if idx < 0:
            idx = widget.findText(str(value), Qt.MatchFixedString)
        if idx >= 0:
            widget.setCurrentIndex(idx)
        elif widget.isEditable():
            widget.setEditText(str(value))
        return
    if isinstance(widget, QLineEdit):
        widget.setText("" if value is None else str(value))
        return
    if hasattr(widget, "setValue"):
        widget.setValue(value)


def reset_parameter_widgets(defaults) -> None:
    """Reset parameter widgets from a mapping or sequence of ``(widget, value)`` pairs."""
    items = defaults.items() if hasattr(defaults, "items") else defaults
    for widget, value in list(items):
        if widget is None:
            continue
        previous = widget.blockSignals(True)
        try:
            set_parameter_widget_value(widget, value)
        finally:
            widget.blockSignals(previous)


def add_parameter_reset_button(
    buttons: QDialogButtonBox,
    defaults,
    *,
    on_reset=None,
    text: str = "Reset Defaults",
) -> QPushButton:
    """Attach a standard reset button to a parameter dialog button box.

    The reset only changes visible widget values. The caller's existing Save
    action remains responsible for persisting those values.
    """
    btn = buttons.addButton(text, QDialogButtonBox.ResetRole)

    def _reset() -> None:
        reset_parameter_widgets(defaults() if callable(defaults) else defaults)
        if callable(on_reset):
            on_reset()

    btn.clicked.connect(_reset)
    return btn


class _WheelGuard(QObject):
    """Prevent accidental wheel edits in parameter widgets while scrolling."""

    def eventFilter(self, obj, event):
        if event.type() == QEvent.Wheel and isinstance(obj, (QAbstractSpinBox, QComboBox)):
            if not obj.hasFocus():
                event.ignore()
                return True
        return super().eventFilter(obj, event)


def install_parameter_wheel_guard(root: QWidget) -> None:
    """Ignore mouse-wheel edits on spin boxes and combos unless focused.

    Parameter dialogs often live inside a scroll area. Without this guard, just
    scrolling the page over a spin box/combobox silently changes persisted
    runtime parameters.
    """
    guard = getattr(root, "_apex_wheel_guard", None)
    if guard is None:
        guard = _WheelGuard(root)
        setattr(root, "_apex_wheel_guard", guard)

    for widget in [root] + list(root.findChildren(QWidget)):
        if isinstance(widget, (QAbstractSpinBox, QComboBox)):
            widget.setFocusPolicy(Qt.StrongFocus)
            try:
                widget.installEventFilter(guard)
            except Exception:
                pass


def build_scroll_param_dialog(
    parent: QWidget,
    title: str,
    *,
    info_text: str = "",
    size: tuple[int, int] = (540, 620),
) -> tuple[QDialog, QVBoxLayout, QDialogButtonBox]:
    """Build a scrollable parameter dialog skeleton.

    Returns (dialog, content_layout, buttons).
    Add sections/groups to content_layout, then call dialog.exec_().
    Buttons (Save | Cancel) sit outside the scroll area.
    """
    dialog = QDialog(parent)
    configure_parameter_dialog(dialog, title, *size)

    outer_layout = QVBoxLayout(dialog)

    scroll = QScrollArea()
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QScrollArea.NoFrame)
    content = QWidget()
    content_layout = QVBoxLayout(content)
    content_layout.setContentsMargins(4, 4, 4, 4)
    content_layout.setSpacing(8)
    scroll.setWidget(content)
    outer_layout.addWidget(scroll, 1)

    if info_text:
        info = QLabel(info_text)
        info.setStyleSheet(_INFO_LABEL_STYLE)
        info.setWordWrap(True)
        content_layout.addWidget(info)

    buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
    outer_layout.addWidget(buttons)

    return dialog, content_layout, buttons
