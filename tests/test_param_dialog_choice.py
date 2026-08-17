"""A combo box that stores a value, not a position.

The parameter dialogs grew a `choice` kind so the remaining windows can be
built from map rows like the rest. Two failure modes are specific to combos and
neither is hypothetical in a config-driven app:

  * Storing the *shown* label. The config would then hold "ASTAP" where the
    loader expects "astap", and the setting silently stops working — the same
    shape as every defect this sweep chased.
  * Selecting by index. Add an option to the list and every existing workspace
    quietly means something else.

The third case is a workspace written before the list existed. Snapping it to
the first entry would rewrite the user's choice on merely opening the dialog,
which is exactly how Step 3's window overwrote Step 7's sky annulus.
"""

from __future__ import annotations

import os
from types import SimpleNamespace

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt5")

from apex.gui.workflow.param_dialog import (            # noqa: E402
    ParamSpec, build_param_form, read_param_form, _spec_default_value,
)

ENGINE = ParamSpec(
    "엔진", attr="wcs_engine", kind="choice", default="astap",
    choices=(("astap", "ASTAP"), ("astrometry", "Astrometry.net"), "none"),
)


@pytest.fixture(scope="module")
def qapp():
    from PyQt5.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_pairs_accept_both_spellings():
    assert ENGINE.choice_pairs() == (
        ("astap", "ASTAP"), ("astrometry", "Astrometry.net"), ("none", "none"))


def test_it_shows_the_label_and_keeps_the_value(qapp):
    P = SimpleNamespace(wcs_engine="astrometry")
    _, widgets = build_param_form(P, (ENGINE,))
    combo = widgets["wcs_engine"]
    assert combo.currentText() == "Astrometry.net"
    assert combo.currentData() == "astrometry"


def test_saving_writes_the_value_not_the_label(qapp):
    P = SimpleNamespace(wcs_engine="astrometry")
    _, widgets = build_param_form(P, (ENGINE,))
    widgets["wcs_engine"].setCurrentIndex(0)
    read_param_form(widgets, P, (ENGINE,))
    assert P.wcs_engine == "astap", "표시 라벨이 설정에 저장되면 로더가 못 읽는다"


def test_selection_follows_the_value_not_the_position(qapp):
    """Adding an option must not change what an existing workspace means."""
    grown = ParamSpec(
        "엔진", attr="wcs_engine", kind="choice", default="astap",
        choices=(("new_first", "새 항목"), ("astap", "ASTAP"),
                 ("astrometry", "Astrometry.net")),
    )
    P = SimpleNamespace(wcs_engine="astap")
    _, widgets = build_param_form(P, (grown,))
    assert widgets["wcs_engine"].currentData() == "astap"


def test_an_unknown_stored_value_survives_being_looked_at(qapp):
    """Opening a dialog must not rewrite a setting it does not recognise."""
    P = SimpleNamespace(wcs_engine="local_solver")
    _, widgets = build_param_form(P, (ENGINE,))
    combo = widgets["wcs_engine"]
    assert combo.currentData() == "local_solver"
    assert "설정 파일 값" in combo.currentText()

    read_param_form(widgets, P, (ENGINE,))
    assert P.wcs_engine == "local_solver", "창을 열었다 닫는 것만으로 값이 바뀌었다"


def test_a_missing_attribute_falls_back_to_the_declared_default(qapp):
    P = SimpleNamespace()
    _, widgets = build_param_form(P, (ENGINE,))
    assert widgets["wcs_engine"].currentData() == "astap"
    assert _spec_default_value(ENGINE) == "astap"


def test_reset_uses_the_first_option_when_no_default_is_declared():
    spec = ParamSpec("모드", attr="x", kind="choice", choices=("a", "b"))
    assert _spec_default_value(spec) == "a"
