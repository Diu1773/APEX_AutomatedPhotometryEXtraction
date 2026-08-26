"""Unit tests for the Qt-free pipeline orchestrator (no FITS data required)."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

from apex.pipeline.base import DeferredStep, PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.pipeline.registry import get_steps, parse_step_range
from apex.pipeline.runner import PipelineRunner


def _ctx(tmp_path: Path, **kw) -> RunContext:
    return RunContext(
        mode="cmd",
        params=kw.pop("params", None),
        result_dir=tmp_path,
        data_dir=tmp_path,
        logger=logging.getLogger("test.pipeline"),
        **kw,
    )


class _Stub(PipelineStep):
    def __init__(self, index, key="s", outs=None, ins=None, fail=False):
        self.index = index
        self.key = key
        self.title = "t"
        self._outs = list(outs or [])
        self._ins = list(ins or [])
        self._fail = fail
        self.ran = False

    def required_inputs(self, ctx):
        return list(self._ins)

    def outputs(self, ctx):
        return list(self._outs)

    def run(self, ctx):
        self.ran = True
        if self._fail:
            raise RuntimeError("boom")
        for p in self._outs:
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text("x", encoding="utf-8")
        return StepResult(self.index, self.key, StepStatus.OK, "done",
                          [str(x) for x in self._outs])


# ── parse_step_range ─────────────────────────────────────────────────────────

@pytest.mark.parametrize("spec,expected", [
    ("1-7", {1, 2, 3, 4, 5, 6, 7}),
    ("4", {4}),
    ("2,4,6", {2, 4, 6}),
    ("3-5", {3, 4, 5}),
    ("5-3", {3, 4, 5}),       # reversed range tolerated
    ("1, 2 , 4-5", {1, 2, 4, 5}),
])
def test_parse_step_range(spec, expected):
    assert parse_step_range(spec) == expected


# ── resolve_plan ─────────────────────────────────────────────────────────────

def test_resolve_plan_from_to_and_only():
    steps = [_Stub(i) for i in range(1, 8)]
    runner = PipelineRunner(steps)
    assert [s.index for s in runner.resolve_plan(from_index=2, to_index=4)] == [2, 3, 4]
    assert [s.index for s in runner.resolve_plan(only={4, 7})] == [4, 7]
    assert [s.index for s in runner.resolve_plan()] == [1, 2, 3, 4, 5, 6, 7]


# ── execution ────────────────────────────────────────────────────────────────

def test_run_all_ok_writes_manifest(tmp_path):
    steps = [_Stub(i, key=f"k{i}", outs=[tmp_path / f"out{i}.txt"]) for i in range(1, 4)]
    report = PipelineRunner(steps).run(_ctx(tmp_path))
    assert report.success
    assert [r.status for r in report.results] == [StepStatus.OK] * 3
    manifest = tmp_path / "pipeline_run.json"
    assert manifest.exists()
    data = json.loads(manifest.read_text(encoding="utf-8"))
    assert data["success"] is True
    assert len(data["steps"]) == 3


def test_skip_when_complete_and_force(tmp_path):
    out = tmp_path / "out1.txt"
    out.write_text("already", encoding="utf-8")
    step = _Stub(1, outs=[out])

    report = PipelineRunner([step]).run(_ctx(tmp_path))
    assert report.results[0].status == StepStatus.SKIPPED
    assert step.ran is False

    step2 = _Stub(1, outs=[out])
    report2 = PipelineRunner([step2]).run(_ctx(tmp_path, force=True))
    assert report2.results[0].status == StepStatus.OK
    assert step2.ran is True


def test_blocked_on_missing_input_stops_chain(tmp_path):
    missing = tmp_path / "upstream.json"
    s1 = _Stub(4, ins=[missing], outs=[tmp_path / "o4.txt"])
    s2 = _Stub(5, outs=[tmp_path / "o5.txt"])
    report = PipelineRunner([s1, s2]).run(_ctx(tmp_path))
    assert report.results[0].status == StepStatus.BLOCKED
    assert len(report.results) == 1  # chain stopped
    assert s1.ran is False and s2.ran is False
    assert not report.success


def test_failure_stops_chain(tmp_path):
    s1 = _Stub(1, fail=True)
    s2 = _Stub(2, outs=[tmp_path / "o2.txt"])
    report = PipelineRunner([s1, s2]).run(_ctx(tmp_path))
    assert report.results[0].status == StepStatus.FAILED
    assert len(report.results) == 1
    assert s2.ran is False


def test_dry_run_executes_nothing_no_manifest(tmp_path):
    # Step 2 needs step 1's output; in dry-run it must NOT block.
    up = tmp_path / "up.json"
    s1 = _Stub(1, outs=[up])
    s2 = _Stub(2, ins=[up], outs=[tmp_path / "o2.txt"])
    report = PipelineRunner([s1, s2]).run(_ctx(tmp_path, dry_run=True))
    assert [r.status for r in report.results] == [StepStatus.PENDING, StepStatus.PENDING]
    assert s1.ran is False and s2.ran is False
    assert not (tmp_path / "pipeline_run.json").exists()


def test_deferred_step_reports_not_implemented_and_stops(tmp_path):
    s1 = DeferredStep(1, "x", "deferred")
    s2 = _Stub(2, outs=[tmp_path / "o2.txt"])
    report = PipelineRunner([s1, s2]).run(_ctx(tmp_path))
    assert report.results[0].status == StepStatus.NOT_IMPLEMENTED
    assert len(report.results) == 1
    assert s2.ran is False


def test_project_state_marked_on_ok(tmp_path):
    marked = []

    class _PS:
        def mark_step_completed(self, idx):
            marked.append(idx)

    step = _Stub(1, outs=[tmp_path / "o1.txt"])
    PipelineRunner([step]).run(_ctx(tmp_path, project_state=_PS()))
    assert marked == [0]  # 1-based step index -> 0-based ProjectState


# ── registry wiring ──────────────────────────────────────────────────────────

def test_registry_shared_steps_shape():
    """Steps 1-7 are shared; each mode then continues into its own.

    This asserted `lc == [1..7]` until 2026-08-19 — and kept asserting it for a
    day after LC gained step 8, because the full suite was not run in between.
    Both mode lists now live here so a new step cannot land unnoticed.
    """
    for mode in ("cmd", "lc"):
        steps = get_steps(mode)
        assert [s.index for s in steps][:7] == [1, 2, 3, 4, 5, 6, 7]
        assert steps[0].key == "scan"
    # The one place that pins both full lists. Two other files used to pin the
    # LC one as well, so adding a step broke three tests instead of this one.
    assert [s.index for s in get_steps("lc")] == [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12]
    assert [s.key for s in get_steps("lc")][7:] == [
        "lctarget", "lclightcurve", "lcdetrend", "lcperiod"]
    assert [s.index for s in get_steps("cmd")] == list(range(1, 13))


# The window's step number and the runner's step number are the same number.
# `step_window_base._record_in_journal` writes `step_index + 1` and
# `runner.run()` maps `step.index - 1` onto ProjectState, so both records join
# on it — and when they disagree a window line lands on the step *after* the one
# it describes, silently.
#
# They did disagree, for LC steps 9-12, from 2026-07-15 (an optional PSF window
# was inserted at LC step 8, pushing the four after it down) until 2026-08-26.
# Nothing caught it because the GUI journal tests only used steps 1 and 4, where
# the two modes agree. Read the window list out of the source rather than
# importing it — the check must not need Qt or a workspace.
WINDOW_TO_STEP_KEY = {
    "File Selection": "scan",
    "Image Crop": "crop",
    "Sky Preview & QC": "sky",
    "Source Detection": "detect",
    "WCS Plate Solving": "wcs",
    "Master Catalog Build": "refbuild",
    "Forced Aperture Phot": "forcedphot",
    "PSF Photometry": "psf",
    "Master ID Editor": "masterid",
    "Zeropoint Calibration": "zeropoint",
    "CMD Plot": "cmdplot",
    "Isochrone Model": "isochrone",
    "Target/Comparison Selection": "lctarget",
    "Light Curve Builder": "lclightcurve",
    "Detrend & Night Merge": "lcdetrend",
    "Period Analysis": "lcperiod",
}


def _window_step_names() -> dict:
    """The two `self.step_names` literals in main_window.py, in order."""
    import ast

    source = Path(__file__).resolve().parents[1] / "apex" / "gui" / "main_window.py"
    found = []
    for node in ast.walk(ast.parse(source.read_text(encoding="utf-8"))):
        if (isinstance(node, ast.Assign) and isinstance(node.value, ast.List)
                and any(isinstance(t, ast.Attribute) and t.attr == "step_names"
                        for t in node.targets)):
            found.append([e.value for e in node.value.elts])
    assert len(found) == 2, f"expected one step list per mode, found {len(found)}"
    return {"cmd": found[0], "lc": found[1]}


def test_the_window_and_the_runner_number_the_same_step_the_same_way():
    names = _window_step_names()
    for mode in ("cmd", "lc"):
        window_names = names[mode]
        for step in get_steps(mode):
            position = step.index - 1
            assert 0 <= position < len(window_names), (
                f"{mode} step {step.index} ({step.key}) has no window at "
                f"position {position}")
            label = window_names[position]
            assert WINDOW_TO_STEP_KEY[label] == step.key, (
                f"{mode}: the window at step {step.index} is {label!r} but the "
                f"runner runs {step.key!r} there — the two records would "
                f"disagree about which step made a file")


def test_lc_step_8_is_a_window_only_step():
    """LC's optional PSF step has a window and no headless registration.

    That hole is why LC's headless indices are 9-12 with nothing at 8. Pinned so
    a future closing of the gap is a deliberate edit, not a renumber.
    """
    assert 8 not in {s.index for s in get_steps("lc")}
    assert _window_step_names()["lc"][7] == "PSF Photometry"


# ── the directory's own history ────────────────────────────────────────────
#
# The manifest is rewritten every run, so a later partial run erases the record
# of an earlier full one — that is how M13's steps 1-7 lost their record on
# 2026-08-23. These check that the journal the runner appends does not.

def test_the_runner_appends_a_journal_line_per_step(tmp_path):
    from apex.utils import run_journal as journal

    runner = PipelineRunner([_Stub(1, outs=[tmp_path / "o1.txt"]),
                             _Stub(2, outs=[tmp_path / "o2.txt"])])
    runner.run(_ctx(tmp_path))

    run = journal.history(tmp_path)[0]
    assert run["mode"] == "cmd" and run["success"] is True
    assert [s["index"] for s in run["steps"]] == [1, 2]
    assert run["environment"].get("apex")


def test_a_second_partial_run_adds_to_the_history_instead_of_replacing_it(tmp_path):
    """The M13 case: run 1-2, then run 2 alone, then ask what step 1 did."""
    from apex.utils import run_journal as journal

    o1, o2 = tmp_path / "o1.txt", tmp_path / "o2.txt"
    PipelineRunner([_Stub(1, outs=[o1]), _Stub(2, outs=[o2])]).run(_ctx(tmp_path))
    PipelineRunner([_Stub(2, key="s2", outs=[o2])]).run(_ctx(tmp_path, force=True),
                                                        only={2})

    manifest = json.loads((tmp_path / "pipeline_run.json").read_text(encoding="utf-8"))
    assert [s["index"] for s in manifest["steps"]] == [2], (
        "the manifest still only describes the last run — that is what it is for")

    assert len(journal.history(tmp_path)) == 2
    assert sorted(journal.latest_steps(tmp_path)) == [1, 2], (
        "step 1 ran, and the directory must still be able to say so")


def test_a_failed_step_is_on_the_record(tmp_path):
    """A run that died is the run most worth being able to read afterwards."""
    from apex.utils import run_journal as journal

    PipelineRunner([_Stub(1, fail=True)]).run(_ctx(tmp_path))

    steps = journal.latest_steps(tmp_path)
    assert steps[1]["status"] == "failed" and steps[1]["message"]
    assert journal.history(tmp_path)[0]["success"] is False


def test_the_journal_records_the_values_a_step_read_while_it_ran(tmp_path):
    """Not the names, and not the config as it stood afterwards — the values.

    A step reads `detect_thresh`; the config is edited; the run ends.
    `parameters_used.json` reports the edited value, because that is what it
    means. The journal has to report what the step was actually given.
    """
    from types import SimpleNamespace

    from apex.utils import run_journal as journal

    live = SimpleNamespace(detect_thresh=5.0, unread=1)
    params = SimpleNamespace(P=live, param_file=None)

    class _Reader(_Stub):
        def run(self, ctx):
            _ = ctx.params.P.detect_thresh     # the recorder sees this one
            return StepResult(index=self.index, key=self.key, status=StepStatus.OK)

    PipelineRunner([_Reader(1)]).run(_ctx(tmp_path, params=params))
    live.detect_thresh = 99.0                  # the config moves on

    recorded = journal.latest_steps(tmp_path)[1]["settings"]
    assert recorded == {"detect_thresh": 5.0}, (
        "'unread' was never asked for, and 99.0 was never used")


def test_a_dry_run_leaves_no_trace_in_the_history(tmp_path):
    """A dry run does nothing to the directory, so it is not part of its history."""
    from apex.utils import run_journal as journal

    PipelineRunner([_Stub(1, outs=[tmp_path / "o.txt"])]).run(_ctx(tmp_path, dry_run=True))
    assert not journal.journal_path(tmp_path).exists()
