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
        params=None,
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
    """Steps 1-7 are shared; CMD then continues into its own 8-12."""
    for mode in ("cmd", "lc"):
        steps = get_steps(mode)
        assert [s.index for s in steps][:7] == [1, 2, 3, 4, 5, 6, 7]
        assert steps[0].key == "scan"
    assert [s.index for s in get_steps("lc")] == [1, 2, 3, 4, 5, 6, 7]
    assert [s.index for s in get_steps("cmd")] == list(range(1, 13))
