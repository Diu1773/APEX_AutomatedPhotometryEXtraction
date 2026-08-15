"""A combo box must not paint a different value than the one it holds.

Step 12's colour combo held "B-V" and painted "B-\\": Qt's sizeHint for a
stylesheet-painted QComboBox under-reserves the native drop-down arrow, so the
widget asked for 52 px, the layout gave it exactly 52 px, and the text field
came out 18 px for a 20 px string. Nothing raised, nothing looked wrong in the
layout — the band name on screen was simply the wrong string.

That failure mode is worse than a clipped button, because the clipped thing
reads as a value. `fit_combo` measures the shortfall through the style that
will do the painting and raises the minimum width by exactly that much.
"""

from __future__ import annotations

import pytest

pytest.importorskip("PyQt5")

from PyQt5.QtCore import QRect                              # noqa: E402
from PyQt5.QtWidgets import (                               # noqa: E402
    QApplication, QComboBox, QStyle, QStyleOptionComboBox,
)

from apex.gui.layout_rules import fit_combo                 # noqa: E402
from apex.gui.theme import apply_theme                      # noqa: E402


@pytest.fixture(scope="module")
def app():
    application = QApplication.instance() or QApplication([])
    apply_theme(application)
    return application


def _text_field_width(box: QComboBox) -> int:
    """The width the style will actually paint the current text into."""
    hint = box.sizeHint()
    width = max(hint.width(), box.minimumWidth())
    option = QStyleOptionComboBox()
    option.initFrom(box)
    option.rect = QRect(0, 0, width, hint.height())
    return box.style().subControlRect(QStyle.CC_ComboBox, option,
                                      QStyle.SC_ComboBoxEditField, box).width()


@pytest.mark.parametrize("items", [
    ["B-V", "V-R"],                       # the case that failed
    ["B", "V", "R"],
    ["g-r", "r-i", "i-z"],
    ["Auto (PSF energy)", "Manual"],       # a long entry must still fit
])
def test_every_entry_fits_after_the_call(app, items):
    box = QComboBox()
    box.addItems(items)
    fit_combo(box)
    field = _text_field_width(box)
    widest = max(box.fontMetrics().horizontalAdvance(text) for text in items)
    assert field >= widest, f"{items}: 글자칸 {field}px < 필요 {widest}px"


def test_it_only_widens(app):
    box = QComboBox()
    box.addItems(["B-V", "V-R"])
    before = box.sizeHint().width()
    fit_combo(box)
    assert box.minimumWidth() >= 0
    assert max(box.sizeHint().width(), box.minimumWidth()) >= before


def test_a_combo_that_already_fits_is_left_alone(app):
    """No blanket widening — a roomy combo must keep its size."""
    box = QComboBox()
    box.addItems(["short"])
    box.setMinimumWidth(400)
    fit_combo(box)
    assert box.minimumWidth() == 400


def test_an_empty_combo_is_not_an_error(app):
    box = QComboBox()
    assert fit_combo(box) is box
