"""파라미터 창의 폭은 손으로 정한 숫자였고, 내용이 넓어지면 조용히 넘쳤다.

Step 8 의 PSF 파라미터 창을 펼치면 **가로 스크롤바가 생겼다**(2026-09-15 실측:
내용 최소 611 px 인데 보이는 너비 581 px, 30 px 넘침). 가로 스크롤바는 값을 보려고
좌우로 밀어야 한다는 뜻이라 그 자체로 결함이다.

원인은 위젯이 아니라 **창 폭을 정하는 방식**이었다. `configure_parameter_dialog` 는
부르는 쪽이 건네는 숫자(Step 8 은 620)를 그대로 쓰는데, 그 숫자는

    내용의 실제 최소 너비        를 안 보고
    세로 스크롤바와 테두리가 먹는 폭   을 안 뺀 값

이라, 위젯이 하나 길어지거나 창이 세로로 길어져 스크롤바가 생기는 순간 넘친다.
F-325(Step 10, 12 px)도 같은 집안이었고 그때는 **내용을 줄여서** 피했다.

고침은 **폭을 내용에 맞추는 것**이다. `fit_parameter_dialog_width` 가 스크롤 내용의
최소 너비에 세로 스크롤바와 테두리·여백을 더해 창을 넓힌다(넓히기만 하고 줄이지
않으므로 이미 맞는 창은 그대로다).

**한 번 맞추는 것으로는 모자란다.** 섹션은 접힌 채로 열리므로 첫 측정은 접힌 내용만
보고, 사용자가 실제로 값을 읽으려고 **펼치는 순간** 필요한 폭이 드러난다. 그래서
섹션 토글에도 다시 맞춘다.

여기서 지키는 것은 그 두 가지다 — **내용에 맞춰 열리는가**, 그리고 **펼친 뒤에도
맞는가**. Step 8 창 전체를 띄우지 않고 같은 구조를 세워 잰다.
"""
from __future__ import annotations

import os

import pytest

pytest.importorskip("PyQt5")

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import QEventLoop, QTimer  # noqa: E402
from PyQt5.QtWidgets import (  # noqa: E402
    QApplication, QCheckBox, QScrollArea, QVBoxLayout,
)

from apex.gui.workflow.ui_helpers import (  # noqa: E402
    build_scroll_param_dialog, create_collapsible_section,
)

#: 창을 넘치게 만들 만큼 넓지만, offscreen 의 800×600 화면 안에는 들어가는 폭.
#: **화면보다 넓게 잡으면 `clamp_to_screen` 이 잘라서 시험이 고침과 무관하게
#: 실패한다** — 넘침은 화면 크기에 딸린 값이다(F-341).
WIDE_TEXT = "Share EPSF per filter " + "x" * 8
NARROW = (300, 400)


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance() or QApplication([])
    yield app


def _settle(app, seconds: float = 0.4) -> None:
    """지연 실행이 실제로 돌게 이벤트 루프를 돌린다.

    폭 맞추기는 `QTimer.singleShot` 으로 미뤄져 있다 — 섹션이 이 호출 뒤에 붙으므로
    그 전에는 내용 너비를 알 수 없기 때문이다. `processEvents` 한 번으로는 안 돈다.
    """
    loop = QEventLoop()
    QTimer.singleShot(int(seconds * 1000), loop.quit)
    loop.exec_()
    app.processEvents()


def _dialog_with_wide_section(expanded: bool):
    """넓은 위젯이 접힌 섹션 안에 든, 일부러 좁게 요청한 파라미터 창."""
    dialog, content_layout, _buttons = build_scroll_param_dialog(
        None, "autofit probe", size=NARROW)
    section, body = create_collapsible_section("probe", initial_expanded=expanded)
    QVBoxLayout(body).addWidget(QCheckBox(WIDE_TEXT))
    content_layout.addWidget(section)
    return dialog, section


def _overflow(dialog) -> int:
    """내용 최소 너비에서 보이는 너비를 뺀 값. 양수면 가로로 밀어야 한다는 뜻이다."""
    area = dialog.findChild(QScrollArea)
    assert area is not None and area.widget() is not None
    return area.widget().minimumSizeHint().width() - area.viewport().width()


def test_the_dialog_opens_wide_enough_for_its_content(qapp):
    """요청한 폭이 좁아도 내용에 맞춰 열린다."""
    dialog, _section = _dialog_with_wide_section(expanded=True)
    dialog.show()
    _settle(qapp)
    try:
        assert dialog.width() > NARROW[0], (
            f"요청한 {NARROW[0]} px 그대로다 — 내용에 맞추지 않았다")
        assert _overflow(dialog) <= 0, (
            f"열자마자 {_overflow(dialog)} px 넘친다")
    finally:
        dialog.deleteLater()


def test_it_fits_again_after_a_section_is_expanded(qapp):
    """**접힌 채 열고 펼치는 것이 실제 순서다.**

    첫 측정은 접힌 내용만 보므로, 펼칠 때 다시 안 맞추면 값을 읽으려는 바로 그
    순간에 가로 스크롤바가 생긴다.
    """
    dialog, section = _dialog_with_wide_section(expanded=False)
    dialog.show()
    _settle(qapp)
    collapsed_w = dialog.width()
    try:
        section.set_expanded(True)
        _settle(qapp)
        assert dialog.width() > collapsed_w, (
            "펼쳤는데 창이 그대로다 — 토글에 다시 맞추기가 안 걸렸다")
        assert _overflow(dialog) <= 0, (
            f"펼친 뒤 {_overflow(dialog)} px 넘친다")
    finally:
        dialog.deleteLater()


def test_it_fits_a_dialog_that_has_no_collapsible_sections(qapp):
    """섹션이 하나도 없어도 내용에 맞춰 열린다 — 맞추기가 섹션에 기대면 안 된다."""
    dialog, content_layout, _buttons = build_scroll_param_dialog(
        None, "no sections", size=NARROW)
    content_layout.addWidget(QCheckBox(WIDE_TEXT))
    dialog.show()
    _settle(qapp)
    try:
        assert dialog.width() > NARROW[0], (
            f"요청한 {NARROW[0]} px 그대로다 — 내용에 맞추지 않았다")
        assert _overflow(dialog) <= 0, (
            f"섹션 없는 창이 {_overflow(dialog)} px 넘친다")
    finally:
        dialog.deleteLater()


def test_the_dialog_is_refitted_when_it_is_first_shown(qapp):
    """창이 뜰 때 다시 맞추는 장치가 걸려 있어야 한다.

    **여기서는 장치가 붙어 있는지만 본다.** 창이 보이기 전에 재면 내용의 최소
    너비가 작게 보고되는 일이 있는데, 그 현상은 내용이 충분히 복잡할 때만 나오고
    이 파일이 세우는 간단한 구조로는 재현되지 않는다(2026-09-15 확인: 이 장치를
    꺼도 합성 창은 −14 px 로 안 넘쳤다). 그래서 이 시험은 **동작이 아니라 장치의
    존재**를 잠근다.

    실제 창에서의 증거는 따로 있다. Step 9(Editor Parameters)를 같은 조건에서
    끄고 켜 보니 **끈 쪽 460 px·넘침 +106 px, 켠 쪽 566 px·넘침 −26 px** 였고,
    그 측정이 이 장치가 필요하다는 근거다. 합성 시험이 그 자리를 대신하지 못한다.
    """
    from apex.gui.workflow.ui_helpers import _ShowRefit

    dialog, content_layout, _buttons = build_scroll_param_dialog(
        None, "refit hook", size=NARROW)
    content_layout.addWidget(QCheckBox(WIDE_TEXT))
    dialog.show()
    _settle(qapp)
    try:
        hook = getattr(dialog, "_apex_show_refit", None)
        assert isinstance(hook, _ShowRefit), (
            "창이 뜰 때 다시 맞추는 이벤트 필터가 안 걸렸다")
    finally:
        dialog.deleteLater()


def test_a_dialog_that_already_fits_is_left_alone(qapp):
    """넓히기만 하고 줄이지 않는다 — 멀쩡한 창을 건드리면 회귀다."""
    from apex.gui.workflow.ui_helpers import fit_parameter_dialog_width

    dialog, _content_layout, _buttons = build_scroll_param_dialog(
        None, "roomy", size=(700, 400))
    dialog.show()
    _settle(qapp)
    try:
        before = dialog.width()
        assert fit_parameter_dialog_width(dialog) == 0, "안 넘치는데 넓혔다"
        assert dialog.width() == before
    finally:
        dialog.deleteLater()
