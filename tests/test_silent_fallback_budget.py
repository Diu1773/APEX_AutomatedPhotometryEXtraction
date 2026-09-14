"""조용한 fallback 이 슬그머니 늘지 못하게 한다.

사장님 지시(2026-09-14, `Main/OPERATOR.md` C-237): *"fallback들 있으면 항상 무슨
로직으로 들어가는지 디버깅관리도 해야돼; 알지 ?"*

**한 번 훑고 끝낼 일이 아니라 지켜야 할 일이다.** 지금 계산 층에는 그럴듯한 값으로
갈아치우면서 아무 기록도 안 남기는 갈래가 144 곳 있다. 한 번에 다 고칠 수는 없지만
**늘어나게 둘 수는 없다.** 이 시험이 그 선을 지킨다.

    줄이는 것    언제든 좋다. 기준선을 같이 낮춰 주면 된다.
    늘리는 것    막는다. 새로 넣은 갈래에 기록을 붙이거나, 없음을 뜻하는 값
                (`None`·`NaN`)을 쓰면 통과한다.

**「없음을 표시하는 갈래」는 세지 않는다.** `None`·`NaN`·빈 표는 아래로 흘러가도
없음이 그대로 드러나므로 조용해도 해롭지 않다. 해로운 것은 진짜 결과와 **겉보기가
같은 값**이다 — 색항 0 은 진짜 색항 0 과 숫자만 보고는 가를 수 없다.

세는 규칙은 `scripts/audit_silent_fallbacks.py` 에 있고 기준선은
`tests/silent_fallback_baseline.json` 이다.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import pytest

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from audit_silent_fallbacks import scan  # noqa: E402

_BASELINE_PATH = Path(__file__).absolute().parent / "silent_fallback_baseline.json"


@pytest.fixture(scope="module")
def measured() -> list[dict]:
    return [r for r in scan() if r["layer"] == "calc" and r["kind"] == "plausible"]


@pytest.fixture(scope="module")
def baseline() -> dict:
    return json.loads(_BASELINE_PATH.read_text(encoding="utf-8"))


def test_the_total_does_not_grow(measured, baseline):
    """전체 수가 늘면 막는다."""
    now, was = len(measured), int(baseline["total"])
    assert now <= was, (
        f"계산 층의 조용한 fallback 이 {was} → {now} 로 늘었다. "
        f"새로 넣은 갈래에 기록을 붙이거나(`apex.utils.fallback_log.note_fallback`) "
        f"없음을 뜻하는 값(None·NaN)을 쓸 것. "
        f"어디인지는 `python -X utf8 scripts/audit_silent_fallbacks.py "
        f"--layer calc --kind plausible --show` 로 본다")


def test_no_single_file_grows(measured, baseline):
    """전체가 같아도 한 파일이 늘고 다른 파일이 줄면 놓친다."""
    now = Counter(r["file"] for r in measured)
    was = {k: int(v) for k, v in baseline["per_file"].items()}
    grew = {f: (was.get(f, 0), n) for f, n in now.items() if n > was.get(f, 0)}
    assert not grew, (
        "이 파일들에서 조용한 fallback 이 늘었다: "
        + ", ".join(f"{f} {a}→{b}" for f, (a, b) in sorted(grew.items())))


def test_the_baseline_is_not_stale(measured, baseline):
    """많이 줄었으면 기준선을 낮추라고 알린다.

    기준선이 실제보다 훨씬 높으면 그만큼 다시 늘어날 여지를 열어 두는 셈이다.
    """
    now, was = len(measured), int(baseline["total"])
    assert was - now <= 10, (
        f"기준선이 {was} 인데 실제는 {now} 다. "
        f"`tests/silent_fallback_baseline.json` 을 지금 값으로 낮출 것 — "
        f"안 낮추면 그만큼 다시 늘 수 있다")


def test_the_counter_itself_still_works():
    """세는 규칙이 망가지면 이 시험 전체가 무의미해진다."""
    rows = scan()
    assert rows, "한 곳도 못 셌다 — 세는 규칙이 깨졌다"
    kinds = {r["kind"] for r in rows}
    assert kinds <= {"absence", "plausible"}
    assert {r["layer"] for r in rows} <= {"calc", "gui"}
