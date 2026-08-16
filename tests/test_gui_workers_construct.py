"""The desktop workers must be constructible, not just importable.

Steps 8 and 10 moved their calculation into `apex.analysis` on 2026-08-17 and
the window kept a thin `QThread` subclass. Every headless test passed and a real
Qt-free run reproduced the stored tables exactly — but nothing constructed the
*GUI* worker, and it could not be constructed at all:

    TypeError: ZeropointCalibrationRunner.__init__() missing 4 required
               positional arguments

With `QThread` first in the bases, sip's initialiser walks the cooperative chain
and lands on the runner with no arguments. Putting the runner first fixes that
and its quieter twin — `QThread.run` shadowing the calculation, so `start()`
would return having done nothing at all.

The live checks run in their own interpreter. Building a `QApplication` and Qt
worker objects inside the shared pytest process ends in a Windows access
violation during fixture teardown; the identical code in its own process
finishes cleanly, so the crash is the harness, not the worker. Static base-order
checks need no Qt at all and stay in-process.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).absolute().parents[1]

pytest.importorskip("PyQt5")


def _probe(body: str, tmp_path: Path) -> str:
    source = f'''
import os
from pathlib import Path
from types import SimpleNamespace
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PyQt5.QtCore import QThread
from PyQt5.QtWidgets import QApplication

app = QApplication([])
tmp = Path(r"{tmp_path}")
params = SimpleNamespace(P=SimpleNamespace(
    result_dir=tmp, data_dir=tmp, cache_dir=tmp / "cache"))
{body}
'''
    out = subprocess.run([sys.executable, "-X", "utf8", "-c", source],
                         cwd=REPO, capture_output=True, text=True, timeout=600)
    assert out.returncode == 0, out.stderr[-2500:]
    return out.stdout


def test_the_zeropoint_worker_can_be_built_and_signals_cross(tmp_path):
    out = _probe('''
from apex.analysis.cmd.zeropoint_runner import ZeropointCalibrationRunner
from apex.gui.workflow.cmd.step10_zeropoint_calibration import ZeropointCalibrationWorker

worker = ZeropointCalibrationWorker(params, tmp, tmp, tmp / "cache")
assert isinstance(worker, QThread)
assert isinstance(worker, ZeropointCalibrationRunner)
# The window and a script run the same code, not similar code.
assert type(worker).run is ZeropointCalibrationRunner.run

seen = []
worker.log.connect(seen.append)
worker.on_log.send("안녕")
assert seen == ["안녕"], "러너 채널이 Qt 신호로 안 건너간다"
print("ZP OK")
''', tmp_path)
    assert "ZP OK" in out


def test_the_psf_worker_can_be_built_and_signals_cross(tmp_path):
    out = _probe('''
from apex.analysis.cmd.psf_photometry_runner import PsfPhotometryRunner
from apex.gui.workflow.cmd.step8_psf_photometry import Step6PSFWorker

worker = Step6PSFWorker([], params, tmp, tmp, tmp / "cache")
assert isinstance(worker, QThread)
assert isinstance(worker, PsfPhotometryRunner)
assert type(worker).run is PsfPhotometryRunner.run

seen = []
worker.log.connect(seen.append)
worker.on_log.send("안녕")
assert seen == ["안녕"]
print("PSF OK")
''', tmp_path)
    assert "PSF OK" in out


def test_a_thread_start_actually_runs_the_calculation(tmp_path):
    """The quieter failure: constructing fine and doing nothing.

    An empty workspace cannot succeed; what matters is that it *tried* — the
    calculation reported through its own channel instead of `QThread.run`
    silently returning.
    """
    out = _probe('''
from apex.gui.workflow.cmd.step10_zeropoint_calibration import ZeropointCalibrationWorker

worker = ZeropointCalibrationWorker(params, tmp, tmp, tmp / "cache")
calls = []
worker.on_error.subscribe(lambda *a: calls.append("error"))
worker.on_finished.subscribe(lambda *a: calls.append("finished"))
worker.start()
assert worker.wait(120_000), "스레드가 안 끝났다"
app.processEvents()
print("CALLS", len(calls))
''', tmp_path)
    count = int(out.split("CALLS")[1].split()[0])
    assert count > 0, "start() 가 계산을 안 돌렸다 (QThread.run 이 가렸다)"


@pytest.mark.parametrize("module_name, worker_name, runner_name", [
    ("apex.gui.workflow.cmd.step10_zeropoint_calibration",
     "ZeropointCalibrationWorker", "ZeropointCalibrationRunner"),
    ("apex.gui.workflow.cmd.step8_psf_photometry",
     "Step6PSFWorker", "PsfPhotometryRunner"),
])
def test_the_runner_precedes_qthread_in_the_bases(module_name, worker_name,
                                                  runner_name):
    """Base order is the fix; a future edit must not quietly restore it."""
    import os

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import importlib

    worker = getattr(importlib.import_module(module_name), worker_name)
    names = [c.__name__ for c in worker.__mro__]
    assert names.index(runner_name) < names.index("QThread"), (
        f"{worker_name} 의 기반 순서가 뒤집혔다 — 생성이 실패한다")
