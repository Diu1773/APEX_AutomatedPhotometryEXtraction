"""Central window-sizing and anti-clipping rules for every APEX window.

The recurring UI complaints — input fields and buttons getting clipped, plots
squished to a sliver, having to drag the window bigger by hand — almost always
trace back to four root causes that each window used to (mis)handle on its own:

1. A window asks for an initial size *larger than the screen*.  Windows then
   clamps it to the monitor and the bottom row (nav buttons / Save) is cut off,
   and the user cannot make it any bigger.
2. A window opens *smaller than the minimum its content needs*, so everything
   is cramped until the user enlarges it.
3. A matplotlib ``FigureCanvas`` reports ``minimumSizeHint() == 10x10`` so, when
   it shares a layout/splitter with a table or control column, it collapses to
   nothing.
4. ``QSplitter`` lets a pane be dragged (or initially laid out) to 0 px because
   ``childrenCollapsible`` defaults to True.

Instead of fixing those by hand in ~30 windows, the rules live here and are
applied centrally: the shared window bases call :func:`fit_window_to_content`
on first show, and plot-heavy windows call :func:`tame_canvas` /
:func:`prevent_collapse` on their canvases and splitters.

Every number here is a *rule*, not a per-window guess — change it once and all
windows follow.
"""

from __future__ import annotations

from PyQt5.QtCore import QSize, Qt
from PyQt5.QtWidgets import (
    QApplication, QDialog, QFrame, QScrollArea, QSizePolicy, QSplitter, QWidget,
)

# ── Tunables (the whole policy in one place) ─────────────────────────────────

# Margin kept between a window and the edges of the usable screen area so the
# title bar and a little breathing room are always reachable.
SCREEN_MARGIN_W = 48
SCREEN_MARGIN_H = 96

# Default minimum a matplotlib canvas may shrink to before scrollbars/clipping
# would be preferable to an unreadable plot.
CANVAS_MIN_W = 360
CANVAS_MIN_H = 260


def _available_geometry(widget: QWidget | None):
    """Usable screen rect (excludes taskbar) for *widget*'s monitor."""
    screen = None
    if widget is not None:
        handle = widget.windowHandle()
        if handle is not None:
            screen = handle.screen()
        if screen is None and hasattr(widget, "screen"):
            try:
                screen = widget.screen()
            except Exception:
                screen = None
    if screen is None:
        screen = QApplication.primaryScreen()
    if screen is None:
        return None
    return screen.availableGeometry()


def clamp_to_screen(width: int, height: int, widget: QWidget | None = None) -> tuple[int, int]:
    """Clamp a desired (width, height) to what fits on *widget*'s screen.

    Use when building an explicit ``resize()`` so a window never opens larger
    than the monitor (root cause #1).
    """
    avail = _available_geometry(widget)
    if avail is None:
        return int(width), int(height)
    max_w = max(360, avail.width() - SCREEN_MARGIN_W)
    max_h = max(300, avail.height() - SCREEN_MARGIN_H)
    return min(int(width), max_w), min(int(height), max_h)


def fit_window_to_content(window: QWidget, *, grow_only: bool = False,
                          recenter: bool = True) -> None:
    """Size a top-level *window* so its content fits, then keep it on-screen.

    The target size is ``clamp(max(current, sizeHint), screen)``:

    * if the content needs more than the current size, the window grows to it
      (fixes "I have to drag the window bigger", root cause #2);
    * if the requested/explicit size is larger than the screen, it shrinks to
      the screen so the bottom row is never clipped (root cause #1).

    Pass ``grow_only=True`` to never shrink below the current size (still
    clamped to the screen). Safe to call once on first show.
    """
    avail = _available_geometry(window)
    if avail is None:
        return
    max_w = max(360, avail.width() - SCREEN_MARGIN_W)
    max_h = max(300, avail.height() - SCREEN_MARGIN_H)

    hint = window.sizeHint()
    cur = window.size()
    want_w = max(cur.width(), hint.width()) if hint.width() > 0 else cur.width()
    want_h = max(cur.height(), hint.height()) if hint.height() > 0 else cur.height()
    if grow_only:
        want_w = max(want_w, cur.width())
        want_h = max(want_h, cur.height())

    target_w = min(want_w, max_w)
    target_h = min(want_h, max_h)
    if grow_only:
        target_w = max(target_w, min(cur.width(), max_w))
        target_h = max(target_h, min(cur.height(), max_h))

    if target_w != cur.width() or target_h != cur.height():
        window.resize(int(target_w), int(target_h))

    if recenter:
        frame = window.frameGeometry()
        frame.moveCenter(avail.center())
        x = min(max(avail.left(), frame.left()), avail.right() - window.width())
        y = min(max(avail.top(), frame.top()), avail.bottom() - window.height())
        window.move(int(x), int(y))


class AutoFitMixin:
    """Class-level first-show auto-fit for windows not based on WindowChromeMixin.

    Mix this in *before* ``QMainWindow`` (e.g. ``class W(AutoFitMixin,
    QMainWindow)``) so the launcher and the couple of legacy raw-``QMainWindow``
    tools get the same screen-clamped content fit the shared bases get.

    It must be a real method on the class — PyQt/sip routes the C++ ``showEvent``
    virtual to a Python reimplementation only when it is defined on a type in
    the MRO, never via an instance attribute, so monkeypatching ``self.showEvent``
    would silently do nothing here.
    """

    def showEvent(self, event):
        super().showEvent(event)
        if getattr(self, "_apex_autosize_done", False):
            return
        self._apex_autosize_done = True
        try:
            fit_window_to_content(self)
        except Exception:
            pass


class FittedDialog(AutoFitMixin, QDialog):
    """A QDialog that auto-fits its content and clamps to the monitor on show.

    Drop-in replacement for ``QDialog(parent)``. Modal dialogs (settings, plot
    popups, custom pickers) don't inherit the window bases, so without this they
    can open larger than the screen — and because Save/Cancel button rows live
    *outside* any internal scroll area, the bottom row gets clipped with no way
    to reach it. Using this everywhere makes every dialog obey the same rule.
    """

    pass


def tame_canvas(canvas: QWidget, *, min_w: int = CANVAS_MIN_W,
                min_h: int = CANVAS_MIN_H, expanding: bool = True) -> QWidget:
    """Stop a matplotlib canvas from collapsing to a sliver (root cause #3).

    Gives the canvas a readable minimum size and an Expanding size policy so it
    claims its share of the layout instead of yielding all space to neighbours.
    Returns the canvas for chaining.
    """
    canvas.setMinimumSize(QSize(int(min_w), int(min_h)))
    if expanding:
        canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
    return canvas


def scroll_wrap(widget: QWidget, *, horizontal: bool = False) -> QScrollArea:
    """Put *widget* in a scroll area so it stops driving the window minimum.

    A ``QScrollArea`` reports a small fixed ``minimumSizeHint`` (~69 px)
    whatever it contains, because scrolling is how it copes with overflow.
    That property is the cure for root cause #1: a ``QTabWidget`` takes the
    *maximum* minimum over its pages, so a single tall page drags the whole
    window past the screen and pushes the nav row out of reach. Wrapping only
    that page drops the window's minimum to the next-tallest page; the page
    scrolls instead of the window growing.

    Measured on the 1280x704 laptop box: Step 7 1156 -> 698 px, Step 6
    727 -> 647, the variable-star tool 1157 -> 605.

    Do **not** wrap a page that already contains a ``QScrollArea`` — nested
    scrolling leaves the user unsure which surface they are scrolling. Thin
    that page out (collapse groups, move rarely-used controls) instead.
    """
    scroll = QScrollArea()
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QFrame.NoFrame)
    if not horizontal:
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
    scroll.setWidget(widget)
    return scroll


def prevent_collapse(splitter: QSplitter, *, min_panes: int | None = None) -> QSplitter:
    """Keep splitter panes from being dragged/laid out to 0 px (root cause #4).

    Sets ``childrenCollapsible(False)`` and a small per-pane handle so neither
    side (typically a control column vs. a plot) can vanish. Returns the
    splitter for chaining.
    """
    splitter.setChildrenCollapsible(False)
    if splitter.handleWidth() < 4:
        splitter.setHandleWidth(6)
    if min_panes:
        for i in range(splitter.count()):
            w = splitter.widget(i)
            if w is not None and w.minimumWidth() < min_panes:
                w.setMinimumWidth(int(min_panes))
    return splitter
