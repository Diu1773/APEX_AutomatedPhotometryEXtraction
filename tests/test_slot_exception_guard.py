"""슬롯 안 예외가 앱을 조용히 죽이지 않아야 한다.

PyQt5 는 슬롯에서 처리되지 않은 파이썬 예외가 나면 `qFatal()` 로 프로세스를 끝낸다.
트레이스백도 로그도 안 남고 `faulthandler` 도 신호가 아니라 못 잡으므로, 화면에는
**앱이 그냥 사라진다.** 2026-09-15 에 그 뒤에 진짜 결함 둘이 숨어 있었다(빠진 임포트
하나, 남은 두 줄 — 둘 다 평범한 `NameError`).

`sys.excepthook` 을 걸면 PyQt 가 그것을 부르고 죽이지 않는다. 여기서 지키는 것은
**가드가 실제로 그 죽음을 막는가**이고, 자식 프로세스로 확인한다 — 죽는 쪽을 같은
프로세스에서 시험하면 시험 실행기까지 끝나기 때문이다.
"""
from __future__ import annotations

import subprocess
import sys
import textwrap

import pytest

pytest.importorskip("PyQt5")

PROBE = textwrap.dedent("""
    import os, sys
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    sys.path.insert(0, {repo!r})
    from PyQt5.QtWidgets import QApplication, QPushButton
    if {guard!r}:
        from apex.utils.app_setup import install_slot_exception_guard
        install_slot_exception_guard()
    app = QApplication(sys.argv[:1])
    b = QPushButton()
    b.clicked.connect(lambda: (_ for _ in ()).throw(NameError("boom")))
    b.click()
    print("SURVIVED")
""")


def _run(guard: bool) -> str:
    from pathlib import Path
    repo = str(Path(__file__).absolute().parents[1])
    p = subprocess.run([sys.executable, "-c", PROBE.format(repo=repo, guard=guard)],
                       capture_output=True, text=True, timeout=180)
    return p.stdout or ""


def test_without_the_guard_a_slot_exception_kills_the_process():
    """가드가 막아 주는 것이 실재하는 죽음인지부터 확인한다."""
    assert "SURVIVED" not in _run(False), (
        "가드 없이도 안 죽는다 — 이 시험이 지키는 것이 없다")


def test_with_the_guard_the_app_stays_alive():
    assert "SURVIVED" in _run(True), "가드를 걸었는데도 프로세스가 죽었다"
