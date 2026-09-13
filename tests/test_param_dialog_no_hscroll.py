"""폼 레이아웃의 이름 칸은 줄마다 공유된다 — 그래서 이름 없는 넓은 위젯이 넘친다.

Step 10 의 파라미터 창을 펼치면 **가로 스크롤바가 생겼다.** 가로 스크롤바는 값을
보려고 좌우로 밀어야 한다는 뜻이라 그 자체로 결함이다.

원인은 `QFormLayout` 의 성질이다. 이름 칸의 너비는 **그 폼의 모든 줄이 나눠 쓰는
하나의 값**이고, 값 칸도 마찬가지다. 그러니 폼의 최소 너비는

    가장 긴 이름 + 가장 넓은 값 + 사이 간격

이 된다. 「External Standard Anchor」 에서는 가장 긴 이름이
``Anchor mag_std to standards:``(161 px)이고 가장 넓은 값이
``Find catalogs for this field (VizieR)`` 버튼(252 px)이었다. **둘이 서로 다른
줄에 있는데도 더해진다.** 161 + 252 + 6 = 419 px 에 섹션 테두리가 붙어 465 px,
보이는 너비는 461 px 이라 12 px 넘쳤다.

고침은 텍스트를 줄이는 것이 아니라 **자리를 옮기는 것**이다. 그 버튼은 이름이
붙는 값이 아니라 동작이므로 줄을 가로지르게 두면 값 칸 경쟁에서 빠진다.
내용 최소 너비가 473 → 306 px 으로 줄고 스크롤바가 사라졌다(2026-09-14 실측).

여기서 지키는 것은 **성질**이다 — 「이름 없는 줄에 놓인 넓은 위젯은 값 칸의
너비를 혼자 정한다」. 그래서 Step 10 창 전체를 띄우지 않고 같은 구조를 세워 잰다.
"""
from __future__ import annotations

import os

import pytest

pytest.importorskip("PyQt5")

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtWidgets import (  # noqa: E402
    QApplication, QFormLayout, QPushButton, QWidget,
)


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance() or QApplication([])
    yield app


def _form_min_width(spanning: bool) -> int:
    """긴 이름 한 줄 + 넓은 버튼 한 줄짜리 폼의 최소 너비."""
    host = QWidget()
    form = QFormLayout(host)
    form.setContentsMargins(0, 0, 0, 0)
    form.addRow("Anchor mag_std to standards:", QPushButton("Enable"))
    wide = QPushButton("Find catalogs for this field (VizieR)")
    if spanning:
        form.addRow(wide)
    else:
        form.addRow("", wide)
    return host.minimumSizeHint().width()


def test_an_unlabelled_wide_widget_adds_to_the_longest_label(qapp):
    """이름 없는 줄에 놓으면 가장 긴 이름과 **더해진다** — 이것이 원인이었다."""
    labelled = _form_min_width(spanning=False)
    spanned = _form_min_width(spanning=True)
    assert labelled > spanned, (
        "이름 칸을 공유하는 성질이 사라졌다 — 이 시험의 전제가 바뀐 것이다")
    assert labelled - spanned > 80, (
        f"차이가 {labelled - spanned} px 밖에 안 난다 — 고침의 효과가 없다")


def test_the_step10_anchor_action_widgets_span_their_row(qapp):
    """Step 10 의 그 두 위젯이 실제로 줄을 가로지르게 놓여 있다.

    창을 통째로 띄우지 않고 **소스에 그렇게 적혀 있는지**를 본다. 창을 띄우면
    워크스페이스와 Gaia 조회가 필요해서 시험이 환경을 타기 때문이다.
    """
    from pathlib import Path

    src = (Path(__file__).absolute().parents[1]
           / "apex/gui/workflow/cmd/step10_zeropoint_calibration.py"
           ).read_text(encoding="utf-8")
    assert "anchor_form.addRow(self.param_anchor_find)" in src, (
        "Find 버튼이 다시 이름 칸을 쓰는 줄로 돌아갔다")
    assert "anchor_form.addRow(self.param_anchor_candidates)" in src, (
        "후보 목록이 다시 이름 칸을 쓰는 줄로 돌아갔다")
    assert 'anchor_form.addRow("", self.param_anchor_find)' not in src
    assert 'anchor_form.addRow("", self.param_anchor_candidates)' not in src


def test_the_candidate_combo_does_not_grow_with_its_contents(qapp):
    """찾은 카탈로그 설명이 길어도 창을 넓히지 않는다.

    VizieR 설명은 한 줄이 100 자를 넘기도 한다. 내용에 맞춰 늘어나는 정책이면
    후보를 채우는 순간 다시 넘친다.
    """
    from pathlib import Path

    src = (Path(__file__).absolute().parents[1]
           / "apex/gui/workflow/cmd/step10_zeropoint_calibration.py"
           ).read_text(encoding="utf-8")
    head = src[src.index("self.param_anchor_candidates = QComboBox()"):]
    head = head[:head.index("anchor_form.addRow(self.param_anchor_candidates)")]
    assert "AdjustToMinimumContentsLength" in head, (
        "후보 목록이 내용 길이에 맞춰 늘어나면 긴 설명 하나로 창이 넘친다")
