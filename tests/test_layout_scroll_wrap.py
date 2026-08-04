"""Layout guard: a tall tab page must not drive the window past the screen.

Root cause #1 in ``apex/gui/layout_rules.py`` — a window whose minimum is
larger than the monitor gets clamped by the OS, and the region below the cut
(the nav row, sliders, Run buttons) becomes permanently unreachable because Qt
refuses to lay the window out below its ``minimumSizeHint``.

``scroll_wrap`` is the remedy: a ``QScrollArea`` reports a small fixed minimum
whatever it holds, so wrapping the tallest page of a ``QTabWidget`` drops the
whole window's minimum to the next-tallest page. These tests pin that
mechanism, since the four windows that rely on it (Step 6, Step 7, Step 8,
variable-star tool) would silently regress if it changed.

Measured on the 1280x704 laptop box when the wrap was introduced:
Step 7 1156 -> 698 px, Step 6 727 -> 647, variable-star tool 1157 -> 605.
"""
from __future__ import annotations

import os

import pytest

pytest.importorskip("PyQt5")

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtWidgets import (  # noqa: E402
    QApplication, QLabel, QScrollArea, QTabWidget, QVBoxLayout, QWidget,
)

from apex.gui.layout_rules import scroll_wrap  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance() or QApplication([])
    yield app


def _tall_page(height: int) -> QWidget:
    page = QWidget()
    lay = QVBoxLayout(page)
    label = QLabel("content")
    label.setMinimumHeight(height)
    lay.addWidget(label)
    return page


def test_scroll_wrap_returns_a_resizable_scroll_area(qapp):
    page = _tall_page(800)
    wrapped = scroll_wrap(page)
    assert isinstance(wrapped, QScrollArea)
    assert wrapped.widgetResizable() is True
    assert wrapped.widget() is page


def test_scroll_wrap_hides_the_page_minimum_from_its_parent(qapp):
    """The whole point: a 800px page must stop demanding 800px."""
    page = _tall_page(800)
    assert page.minimumSizeHint().height() >= 800
    wrapped = scroll_wrap(page)
    assert wrapped.minimumSizeHint().height() < 200, (
        "scroll area should report a small minimum regardless of its content"
    )


def test_scroll_wrap_keeps_the_content_reachable(qapp):
    """Hiding the minimum must not hide the content — it has to scroll."""
    page = _tall_page(800)
    wrapped = scroll_wrap(page)
    wrapped.resize(400, 300)
    wrapped.show()
    qapp.processEvents()
    viewport = wrapped.viewport().height()
    reach = viewport + wrapped.verticalScrollBar().maximum()
    assert viewport > 0
    assert reach >= wrapped.widget().height() - 1, (
        "viewport plus scroll range must cover the wrapped page"
    )


def _wide_column(width: int) -> QWidget:
    page = QWidget()
    lay = QVBoxLayout(page)
    lay.setContentsMargins(0, 0, 0, 0)
    label = QLabel("column")
    lay.addWidget(label)
    page.setFixedWidth(width)
    return page


def test_plain_wrap_asks_for_less_than_a_wide_column_needs(qapp):
    """Baseline for ``preferred_width``: Qt caps a scroll area's width hint at a
    generic character-cell box, so a wrapped fixed-width column is handed less
    than its content and grows a horizontal scrollbar even on a wide screen."""
    wrapped = scroll_wrap(_wide_column(540), horizontal=True)
    assert wrapped.sizeHint().width() < 540


def test_preferred_width_asks_for_the_column_width(qapp):
    wrapped = scroll_wrap(_wide_column(540), horizontal=True, preferred_width=570)
    assert wrapped.sizeHint().width() == 570


def test_preferred_width_still_lets_the_column_shrink(qapp):
    """The width preference must not become a floor — that would hand the window
    back the hard minimum the wrap exists to remove."""
    wrapped = scroll_wrap(_wide_column(540), horizontal=True, preferred_width=570)
    assert wrapped.minimumSizeHint().width() < 300


def test_preferred_width_is_adjustable_after_wrapping(qapp):
    """Step 9 sets its column width at the end of the build, after the wrap."""
    wrapped = scroll_wrap(_wide_column(540), horizontal=True, preferred_width=1)
    wrapped.set_preferred_width(570)
    assert wrapped.sizeHint().width() == 570


def test_tab_widget_minimum_follows_the_tallest_page(qapp):
    """Baseline: a QTabWidget takes the MAX over its pages, so one fat tab
    drags the window up — this is the failure being guarded against."""
    tabs = QTabWidget()
    tabs.addTab(_tall_page(200), "small")
    tabs.addTab(_tall_page(900), "fat")
    assert tabs.minimumSizeHint().height() >= 900


def test_wrapping_the_fat_tab_drops_the_tab_widget_minimum(qapp):
    """After wrapping, the minimum falls back to the second-tallest page."""
    tabs = QTabWidget()
    tabs.addTab(_tall_page(200), "small")
    tabs.addTab(scroll_wrap(_tall_page(900)), "fat")
    height = tabs.minimumSizeHint().height()
    assert height < 400, f"expected the 200px page to set the floor, got {height}"


@pytest.mark.parametrize(
    "module_path, needle",
    [
        ("apex/gui/workflow/step6_ref_build.py", 'scroll_wrap(plot_tab)'),
        ("apex/gui/workflow/step7_forced_aperture_phot.py", 'scroll_wrap(tab0)'),
        # LC Step 9 wraps its fixed 540px control column instead of a tab: the
        # column plus the 841px plot splitter demanded 1451px, so every window
        # narrower than that clipped the column's buttons and let Phase Folding
        # overlap the light curve. Measured 1451 -> 1016 px.
        ("apex/gui/workflow/lc/step9_lightcurve_builder.py", "scroll_wrap(left_col"),
        # Step 8 (qc_tab) and the variable-star tool (mm_tab) get the same
        # treatment; add them here once their files land — both currently sit
        # inside a larger in-flight change.
    ],
)
def test_known_fat_tabs_stay_wrapped(module_path, needle):
    """These windows only fit a 1280x704 screen because their tallest tab is
    wrapped. Unwrapping one silently pushes its nav row off-screen again."""
    from pathlib import Path

    source = (Path(__file__).resolve().parents[1] / module_path).read_text(
        encoding="utf-8"
    )
    assert needle in source, (
        f"{module_path} no longer wraps its tallest tab ({needle}); the window "
        "will exceed the laptop screen and clip its bottom row"
    )
