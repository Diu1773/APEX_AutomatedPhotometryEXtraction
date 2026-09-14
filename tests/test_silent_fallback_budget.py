"""근거 없이 조용한 갈래가 계산 층에 남아 있으면 안 된다.

사장님 지시(2026-09-14, `Main/OPERATOR.md` C-237): *"fallback들 있으면 항상 무슨
로직으로 들어가는지 디버깅관리도 해야돼; 알지 ?"*

**한 번 훑고 끝낼 일이 아니라 지켜야 할 일이다.** 2026-09-14 에 계산 층의 137 곳을
전부 처리해 **0 이 되었다.** 실제로 실패를 감추던 둘은 고쳤다 — 포화 별 수가 0 으로
읽혀 그 프레임이 오히려 기준 프레임 후보가 되던 것, 그리고 영점이 깨진 밴드의 색이
조용히 어긋나던 것. 나머지는 **왜 조용해도 되는지를 코드에 적었다.**

새로 만드는 갈래는 둘 중 하나여야 한다.

    실행 중에 남긴다    `note_fallback` · 사유 열(`drift_note`) · 실패 카운터
    코드에 적는다       `# fallback-ok: <왜 조용해도 되는지>`

**「없음을 표시하는 갈래」는 애초에 안 센다.** `None`·`NaN`·빈 표·사유를 담은 값은
아래로 흘러가도 없음이 그대로 드러나므로 조용해도 해롭지 않다. 해로운 것은 진짜
결과와 **겉보기가 같은 값**이다 — 색항 0 은 진짜 색항 0 과 숫자만 보고는 가를 수
없다.

세는 규칙은 `scripts/audit_silent_fallbacks.py`, 기준선은
`tests/silent_fallback_baseline.json` 이다.
"""
from __future__ import annotations

import ast
import json
import sys
from collections import Counter
from pathlib import Path

import pytest

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from audit_silent_fallbacks import MARKER, _justified, scan  # noqa: E402

_BASELINE_PATH = Path(__file__).absolute().parent / "silent_fallback_baseline.json"


@pytest.fixture(scope="module")
def measured() -> list[dict]:
    return [r for r in scan() if r["layer"] == "calc" and r["kind"] == "plausible"]


@pytest.fixture(scope="module")
def baseline() -> dict:
    return json.loads(_BASELINE_PATH.read_text(encoding="utf-8"))


def test_no_unjustified_silent_fallback_remains(measured, baseline):
    """계산 층에 근거 없는 조용한 갈래가 하나도 없어야 한다."""
    now, was = len(measured), int(baseline["total"])
    assert now <= was, (
        f"계산 층의 근거 없는 조용한 갈래가 {was} → {now} 로 늘었다. "
        f"실행 중에 사유를 남기거나(`apex.utils.fallback_log.note_fallback`), "
        f"없음을 뜻하는 값(None·NaN)을 쓰거나, "
        f"`# {MARKER} <왜 조용해도 되는지>` 를 적을 것. 어디인지는 "
        f"`python -X utf8 scripts/audit_silent_fallbacks.py "
        f"--layer calc --kind plausible --show` 로 본다")


def test_no_single_file_grows(measured, baseline):
    """전체가 같아도 한 파일이 늘고 다른 파일이 줄면 놓친다."""
    now = Counter(r["file"] for r in measured)
    was = {k: int(v) for k, v in baseline["per_file"].items()}
    grew = {f: (was.get(f, 0), n) for f, n in now.items() if n > was.get(f, 0)}
    assert not grew, (
        "이 파일들에서 근거 없는 조용한 갈래가 늘었다: "
        + ", ".join(f"{f} {a}→{b}" for f, (a, b) in sorted(grew.items())))


def test_the_baseline_is_not_stale(measured, baseline):
    """기준선이 실제보다 높으면 그만큼 다시 늘 여지를 열어 두는 셈이다."""
    now, was = len(measured), int(baseline["total"])
    assert was - now <= 0, (
        f"기준선이 {was} 인데 실제는 {now} 다. "
        f"`tests/silent_fallback_baseline.json` 을 지금 값으로 낮출 것")


def _handler(source: str) -> ast.ExceptHandler:
    return next(n for n in ast.walk(ast.parse(source))
                if isinstance(n, ast.ExceptHandler))


def test_a_marker_without_a_reason_does_not_count():
    """표시만 달고 이유를 안 적으면 안 쳐 준다.

    **이것은 검사를 끄는 장치가 아니라 근거를 남기게 하는 장치다.** 표시를 달려면
    한 줄을 써야 하고, 그 한 줄이 나중에 읽는 사람에게 답이 된다.
    """
    bare = "try:\n    x = f()\nexcept Exception:\n    # " + MARKER + "\n    x = 0\n"
    assert not _justified(bare.splitlines(), _handler(bare))


def test_a_marker_with_a_reason_counts():
    withy = ("try:\n    x = f()\nexcept Exception:\n    # " + MARKER
             + " 이 함수의 계약이 값이거나 기본값이다\n    x = 0\n")
    assert _justified(withy.splitlines(), _handler(withy))


def test_the_counter_itself_still_works():
    """세는 규칙이 망가지면 이 시험 전체가 무의미해진다."""
    rows = scan()
    assert rows, "한 곳도 못 셌다 — 세는 규칙이 깨졌다"
    assert {r["kind"] for r in rows} <= {"absence", "plausible", "justified"}
    assert {r["layer"] for r in rows} <= {"calc", "gui"}
    assert any(r["kind"] == "justified" for r in rows), (
        "근거를 적어 둔 갈래를 하나도 못 알아봤다 — 표시를 못 읽고 있다")
