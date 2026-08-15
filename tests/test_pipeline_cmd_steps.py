"""CMD Steps 8-12 in the headless registry.

Two things are easy to get wrong here and neither shows up as an error.

First, PyQt5 is an *optional* dependency (the ``gui`` extra) while Steps 8 and
10 still drive the GUI module's worker. A base install must be told what to
install, not meet an ImportError partway through a pipeline.

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


def test_cmd_gains_steps_8_to_12_and_lc_does_not():
    cmd = [s.index for s in get_steps("cmd")]
    assert cmd == list(range(1, 13))
    assert [s.index for s in get_steps("lc")] == list(range(1, 8))


def test_the_two_that_execute_are_psf_and_zeropoint():
    assert isinstance(_step("cmd", 8), PsfPhotometryStep)
    assert isinstance(_step("cmd", 10), ZeropointStep)


@pytest.mark.parametrize("index", [9, 11, 12])
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


def test_missing_pyqt_says_what_to_install_rather_than_raising(tmp_path, monkeypatch):
    step = PsfPhotometryStep()
    monkeypatch.setattr("apex.pipeline.steps.psf._qt_available", lambda: False)
    result = step.run(_ctx(tmp_path))
    assert result.status == StepStatus.NOT_IMPLEMENTED
    assert "gui" in result.message and "PyQt5" in result.message


def test_zeropoint_says_the_same_without_pyqt(tmp_path, monkeypatch):
    monkeypatch.setattr("apex.pipeline.steps.zeropoint._qt_available", lambda: False)
    result = ZeropointStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.NOT_IMPLEMENTED
    assert "gui" in result.message


def test_psf_completeness_is_the_signature_not_the_directory(tmp_path):
    """A part-run leaves tables behind; only the signature means finished."""
    step, ctx = PsfPhotometryStep(), _ctx(tmp_path)
    psf_dir = tmp_path / "cmd_psf"
    psf_dir.mkdir(parents=True)
    (psf_dir / "photometry_a.fit.tsv").write_text("partial", encoding="utf-8")
    assert not step.is_complete(ctx)
    (psf_dir / "psf_output_signature.json").write_text("{}", encoding="utf-8")
    assert step.is_complete(ctx)


def test_psf_is_blocked_without_step7_tables(tmp_path, monkeypatch):
    monkeypatch.setattr("apex.pipeline.steps.psf._qt_available", lambda: True)
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
