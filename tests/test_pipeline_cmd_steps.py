"""CMD Steps 8-12 in the headless registry.

Two things are easy to get wrong here and neither shows up as an error.

First, PyQt5 is an *optional* dependency (the ``gui`` extra). Step 10's
calculation moved to ``apex.analysis.cmd.zeropoint_runner`` on 2026-08-16 and
needs no Qt at all; Step 8 still drives the GUI module's worker, so a base
install must be told what to install rather than meet an ImportError partway
through a pipeline.

Second, Step 8 writes its per-frame tables as it goes and its signature only
after every frame succeeded — so "the directory exists" is not the same
question as "this step finished". Step 10 switches to PSF magnitudes the
moment it sees a valid signature, which makes a signature left behind by a
half-finished run a silent change to the published numbers.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from apex.pipeline.base import DeferredStep, StepStatus
from apex.pipeline.context import RunContext
from apex.pipeline.registry import get_steps
from apex.pipeline.steps.psf import PsfPhotometryStep, _frames_from_step7
from apex.pipeline.steps.zeropoint import ZeropointStep, _photometry_source


def _ctx(tmp_path: Path) -> RunContext:
    return RunContext(
        mode="cmd", params=None, result_dir=tmp_path, data_dir=tmp_path,
        logger=logging.getLogger("test.pipeline.cmd"),
    )


def _step(mode: str, index: int):
    return next(s for s in get_steps(mode) if s.index == index)


def test_each_mode_gains_its_own_steps_after_the_shared_seven():
    """LC reached only 7 until 2026-08-19, and not because its science needed a
    window — the services were already Qt-free. Its target choice had no config
    surface, so a batch run could not say which star the curve was of. Step 8
    gave it one; step 9, the same day, walks through and builds the curve."""
    assert [s.index for s in get_steps("cmd")] == list(range(1, 13))
    assert [s.index for s in get_steps("lc")] == list(range(1, 10))
    assert [s.key for s in get_steps("lc")][-2:] == ["lctarget", "lclightcurve"]


def test_the_two_that_execute_are_psf_and_zeropoint():
    assert isinstance(_step("cmd", 8), PsfPhotometryStep)
    assert isinstance(_step("cmd", 10), ZeropointStep)


# The list has emptied down to one. 12 left on 2026-08-17, when the isochrone fit
# got a config surface for the settings that decide its answer. 11 left on
# 2026-08-19: it had been deferred as "a viewer with nothing to port", which was
# true of the viewer and false of the figure — a run that measures a cluster
# should leave the picture of it. 9 stays because it is interactive by nature and
# Step 10 does not need it.
@pytest.mark.parametrize("index", [9])
def test_the_rest_are_deferred_but_still_recognise_finished_work(index, tmp_path):
    """A deferred step must still see a GUI-produced result as complete."""
    step = _step("cmd", index)
    assert isinstance(step, DeferredStep)
    ctx = _ctx(tmp_path)
    assert not step.is_complete(ctx)
    for out in step.outputs(ctx):
        out.mkdir(parents=True, exist_ok=True)
        (out / "produced_by_the_gui.csv").write_text("x", encoding="utf-8")
    assert step.is_complete(ctx)


def test_the_master_id_editor_is_marked_interactive():
    assert _step("cmd", 9).interactive is True


def test_psf_needs_no_qt_at_all():
    """Step 8 used to report NOT_IMPLEMENTED without PyQt5 because it reached
    into the GUI module for its worker — so a script could not do PSF photometry
    at all. The calculation now lives in `apex.analysis.cmd.psf_photometry_runner`
    (2026-08-16). Probed in a fresh interpreter so an already-imported Qt from
    another test cannot mask the answer."""
    import subprocess
    import sys

    import apex.pipeline.steps.psf as step_module

    assert not hasattr(step_module, "_qt_available")

    probe = (
        "import sys;"
        "import apex.pipeline.registry;"
        "import apex.analysis.cmd.psf_photometry_runner;"
        "print([m for m in sys.modules if m.startswith(('PyQt5', 'apex.gui'))])"
    )
    out = subprocess.run([sys.executable, "-X", "utf8", "-c", probe],
                         cwd=Path(__file__).absolute().parents[1],
                         capture_output=True, text=True, timeout=180)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "[]", f"Qt/GUI 가 딸려 들어온다: {out.stdout}"


def test_the_window_and_the_script_run_the_same_psf_code():
    """The same function object, not a copy kept in agreement."""
    pytest.importorskip("PyQt5")
    from apex.analysis.cmd.psf_photometry_runner import PsfPhotometryRunner
    from apex.gui.workflow.cmd.step8_psf_photometry import Step6PSFWorker

    assert Step6PSFWorker.run is PsfPhotometryRunner.run
    differing = [name for name in dir(PsfPhotometryRunner)
                 if not name.startswith("__")
                 and getattr(Step6PSFWorker, name, None)
                 is not getattr(PsfPhotometryRunner, name, None)]
    assert not differing, f"GUI 가 코어와 다른 구현을 갖는다: {differing}"


def test_zeropoint_needs_no_qt_at_all():
    """Step 10 used to report NOT_IMPLEMENTED without PyQt5 because it reached
    into the GUI module for its worker. The calculation now lives in
    `apex.analysis.cmd.zeropoint_runner`, so a script can calibrate zeropoints
    on a Qt-free install (2026-08-16).

    Checked in a fresh interpreter: reloading the module in-process would hand
    the GUI subclass a stale function object and break the identity check below.
    """
    import subprocess
    import sys

    import apex.pipeline.steps.zeropoint as step_module

    source = Path(step_module.__file__).read_text(encoding="utf-8")
    assert "PyQt5" not in source
    assert not hasattr(step_module, "_qt_available")

    probe = (
        "import sys;"
        "import apex.pipeline.steps.zeropoint;"
        "import apex.analysis.cmd.zeropoint_runner;"
        "print([m for m in sys.modules if m.startswith('PyQt5')])"
    )
    out = subprocess.run([sys.executable, "-X", "utf8", "-c", probe],
                         cwd=Path(__file__).absolute().parents[1],
                         capture_output=True, text=True, timeout=180)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "[]", f"Qt 가 딸려 들어온다: {out.stdout}"


def test_the_window_and_the_script_run_the_same_zeropoint_code():
    """Not equivalent code — the same function object. The GUI subclass adds the
    thread and the Qt signals and binds the runner's `run` explicitly, because
    `QThread` sits first in the MRO and brings its own empty one."""
    pytest.importorskip("PyQt5")
    from apex.analysis.cmd.zeropoint_runner import ZeropointCalibrationRunner
    from apex.gui.workflow.cmd.step10_zeropoint_calibration import (
        ZeropointCalibrationWorker,
    )

    assert ZeropointCalibrationWorker.run is ZeropointCalibrationRunner.run
    assert issubclass(ZeropointCalibrationWorker, ZeropointCalibrationRunner)


def test_the_calculation_speaks_its_own_language_not_qt_s():
    """The core announces on channels it owns.

    The first cut had the core call `.emit()` so 1,767 lines could move without
    being touched; that left Qt's vocabulary in `apex.analysis`, where a reader
    would fairly assume Qt was involved. Renaming cost five call sites here and
    would have cost thirty-eight after Step 8 moves, so it was done first.
    """
    source = (Path(__file__).absolute().parents[1]
              / "apex/analysis/cmd/zeropoint_runner.py").read_text(encoding="utf-8")
    assert ".emit(" not in source
    assert "on_progress.send(" in source

    from apex.analysis.cmd.zeropoint_runner import ZeropointCalibrationRunner

    assert ZeropointCalibrationRunner._CHANNELS == ("progress", "log", "finished", "error")


def test_a_subscriber_that_raises_does_not_end_the_run():
    """An hour of photometry must not be lost to a progress bar throwing."""
    from apex.analysis.worker_signals import Channel

    seen = []
    channel = Channel("progress")
    channel.subscribe(lambda *args: seen.append(args))
    channel.subscribe(lambda *args: 1 / 0)
    channel.subscribe(lambda *args: seen.append(("after",) + args))
    channel.send(3, 10, "frame.fit")
    assert seen == [(3, 10, "frame.fit"), ("after", 3, 10, "frame.fit")]


def test_psf_completeness_is_the_signature_not_the_directory(tmp_path):
    """A part-run leaves tables behind; only the signature means finished."""
    step, ctx = PsfPhotometryStep(), _ctx(tmp_path)
    psf_dir = tmp_path / "cmd_psf"
    psf_dir.mkdir(parents=True)
    (psf_dir / "photometry_a.fit.tsv").write_text("partial", encoding="utf-8")
    assert not step.is_complete(ctx)
    (psf_dir / "psf_output_signature.json").write_text("{}", encoding="utf-8")
    assert step.is_complete(ctx)


def test_psf_is_blocked_without_step7_tables(tmp_path):
    (tmp_path / "step7_forced_phot").mkdir(parents=True)
    result = PsfPhotometryStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.BLOCKED


def test_frames_come_from_step7_and_keep_a_stable_order(tmp_path):
    out = tmp_path / "step7_forced_phot"
    out.mkdir(parents=True)
    for name in ("b.fit", "a.fit", "c.fit"):
        (out / f"photometry_{name}.tsv").write_text("x", encoding="utf-8")
    (out / "photometry_index.csv").write_text("x", encoding="utf-8")   # not a frame
    assert _frames_from_step7(tmp_path) == ["a.fit", "b.fit", "c.fit"]


def test_zeropoint_completeness_is_the_table_steps_11_and_12_read(tmp_path):
    step, ctx = ZeropointStep(), _ctx(tmp_path)
    zp = tmp_path / "cmd_zeropoint"
    zp.mkdir(parents=True)
    (zp / "frame_zeropoint.csv").write_text("x", encoding="utf-8")
    assert not step.is_complete(ctx)
    (zp / "median_by_ID_filter_wide_cmd.csv").write_text("x", encoding="utf-8")
    assert step.is_complete(ctx)


@pytest.mark.parametrize("body,expected", [
    ("filter,zp,photometry_source\nB,-3.9,psf\n", "psf"),
    ("filter,zp,photometry_source\nB,-3.9,aperture\n", "aperture"),
])
def test_the_photometry_source_is_read_back_not_assumed(tmp_path, body, expected):
    """Step 10 switches to PSF silently; the run report must say which it used."""
    zp = tmp_path / "cmd_zeropoint"
    zp.mkdir(parents=True)
    (zp / "zp_fit_coefficients.csv").write_text(body, encoding="utf-8")
    assert _photometry_source(zp) == expected


def test_an_absent_coefficients_file_is_unknown_not_a_guess(tmp_path):
    assert _photometry_source(tmp_path) == "unknown"


def test_nothing_still_reaches_for_a_qt_signal_on_the_runners():
    """The move renamed the announcements; every caller had to follow.

    Step 8 failed on its first Qt-free run with
    `'PsfPhotometryRunner' object has no attribute 'error'` — the pipeline was
    still doing `worker.error.connect(...)`, which only existed while the
    worker was a QThread subclass. The headless scripts had the same shape, and
    both also replaced `worker._log` wholesale instead of subscribing.
    """
    import re

    repo = Path(__file__).absolute().parents[1]
    watched = [
        "apex/pipeline/steps/psf.py",
        "apex/pipeline/steps/zeropoint.py",
        "scripts/run_step8_headless.py",
        "scripts/run_step10_headless.py",
    ]
    def code_only(text: str) -> str:
        """Drop docstrings and comments — the history is written there on
        purpose and naming Qt in it is not the same as importing it."""
        import io
        import tokenize

        kept = []
        for tok in tokenize.generate_tokens(io.StringIO(text).readline):
            if tok.type in (tokenize.COMMENT, tokenize.STRING):
                continue
            kept.append(tok.string)
        return " ".join(kept)

    offenders = []
    for name in watched:
        source = code_only((repo / name).read_text(encoding="utf-8"))
        for pattern, why in (
            (r"worker\.(progress|log|finished|error|worker_status|frame_done"
             r"|epsf_ready|residual_ready)\.connect\(", "Qt 신호에 connect"),
            (r"worker\._log\s*=", "메서드를 갈아끼움"),
            (r"\bPyQt5\b", "Qt 임포트"),
            (r"QCoreApplication", "Qt 앱 생성"),
        ):
            if re.search(pattern, source):
                offenders.append(f"{name}: {why}")
    assert not offenders, (
        "헤드리스 경로가 다시 Qt 를 붙들고 있다 — 채널을 subscribe 할 것: "
        + ", ".join(offenders)
    )


def test_the_ap_vs_psf_figure_is_drawn_by_shared_code():
    """`step8_ap_vs_psf_comparison.png` used to need a window.

    The drawing was a 180-line method on the Step 8 dialog, wired to three spin
    boxes and a canvas, so the figure existed only if somebody had opened that
    tab. A headless run wrote every other Step 8 product and not this one — the
    figure that says whether the two photometries agree. Measured 2026-08-18.
    """
    import io
    import tokenize
    from pathlib import Path

    from apex.analysis.cmd.psf_photometry_runner import (
        draw_ap_vs_psf, filter_ap_vs_psf,
    )

    assert callable(draw_ap_vs_psf) and callable(filter_ap_vs_psf)

    repo = Path(__file__).absolute().parents[1]
    window = (repo / "apex/gui/workflow/cmd/step8_psf_photometry.py").read_text(encoding="utf-8")
    code = " ".join(
        tok.string
        for tok in tokenize.generate_tokens(io.StringIO(window).readline)
        if tok.type not in (tokenize.COMMENT, tokenize.STRING)
    )
    assert "draw_ap_vs_psf" in code, "창이 공용 그리기를 안 쓴다"
    # The window must not have kept its own copy of the drawing.
    assert "add_subplot ( 121 )" not in code.replace("(", " ( ").replace(")", " ) ")

    runner = (repo / "apex/analysis/cmd/psf_photometry_runner.py").read_text(encoding="utf-8")
    assert "step8_ap_vs_psf_comparison.png" in runner, "배치가 그 그림을 안 쓴다"


def test_the_comparison_cuts_are_applied_once_and_reported(tmp_path):
    """The cut counts have to survive into the caption, or the figure lies."""
    import numpy as np
    import pandas as pd
    from matplotlib.figure import Figure

    from apex.analysis.cmd.psf_photometry_runner import (
        draw_ap_vs_psf, filter_ap_vs_psf,
    )

    n = 300
    df = pd.DataFrame({
        "mag_ap": np.linspace(15, 20, n),
        "mag_psf": np.linspace(15, 20, n) + 0.01,
        "FILTER": ["B"] * n,
        "snr_psf": np.linspace(1, 100, n),
        "qfit": np.full(n, 0.5),
    })
    kept, before = filter_ap_vs_psf(df, snr_min=50.0)
    assert before == n
    assert 0 < len(kept) < n, "SNR 절단이 안 걸렸다"

    text = draw_ap_vs_psf(Figure(figsize=(10, 4)), kept, before)
    assert f"N={len(kept)}/{before}" in text

    empty = draw_ap_vs_psf(Figure(figsize=(10, 4)), kept.iloc[0:0], 0)
    assert "No data" in empty, "빈 결과가 빈 축으로 보이면 안 된다"
