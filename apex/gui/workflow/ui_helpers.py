"""Small shared UI helpers for workflow step controls."""
from __future__ import annotations

from PyQt5.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QGroupBox,
    QLabel,
    QPushButton,
    QScrollArea,
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
    """Return a consistently styled workflow parameter button."""
    btn = QPushButton(text)
    btn.setMinimumHeight(34)
    btn.setStyleSheet(PARAM_BUTTON_STYLE)
    return btn


def configure_parameter_dialog(dialog: QDialog, title: str, width: int = 560, height: int = 620) -> None:
    """Apply common title, size, and visual style to a parameter dialog."""
    dialog.setWindowTitle(title)
    dialog.resize(width, height)
    dialog.setStyleSheet(PARAM_DIALOG_STYLE)


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
    """Return a consistently styled cache maintenance button."""
    button = QPushButton(text)
    button.setMinimumHeight(34)
    button.setStyleSheet(_CACHE_ACTION_STYLE)
    return button


_INFO_LABEL_STYLE = "QLabel { background-color: #E3F2FD; padding: 10px; margin-bottom: 10px; }"


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
