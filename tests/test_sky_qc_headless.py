"""Headless sky-preview/QC (Step 3) tests.

Step 3 persists nothing in headless mode (parameters come from config). The
``SkyQCStep`` is a cheap optional no-op that returns OK and does not block the
chain. No Qt, no real data.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

pytest.importorskip("pandas")

from apex.config.parameters_cmd import read_params
from apex.pipeline.base import StepStatus
from apex.pipeline.context import RunContext
from apex.pipeline.steps.sky_qc import SkyQCStep
from apex.utils.step_paths import step1_dir

_EXAMPLE_TOML = Path(__file__).resolve().parents[1] / "parameters.example.toml"


def _ctx(tmp_path: Path) -> RunContext:
    data_dir = tmp_path / "data"
    result_dir = tmp_path / "result"
    cache_dir = result_dir / "cache"
    for d in (data_dir, result_dir, cache_dir):
        d.mkdir(parents=True, exist_ok=True)
    params = read_params(_EXAMPLE_TOML)
    params.P.data_dir = data_dir
    params.P.result_dir = result_dir
    params.P.cache_dir = cache_dir
    return RunContext(
        mode="cmd",
        params=params,
        result_dir=result_dir,
        data_dir=data_dir,
        logger=logging.getLogger("test.sky_qc"),
    )


def test_sky_qc_step_identity():
    step = SkyQCStep()
    assert step.index == 3
    assert step.key == "sky"
    assert step.interactive is False


def test_sky_qc_run_returns_ok_and_writes_nothing(tmp_path):
    ctx = _ctx(tmp_path)
    step = SkyQCStep()

    result = step.run(ctx)
    assert result.status == StepStatus.OK
    assert result.ok
    # No files written.
    assert result.outputs == []
    assert step.outputs(ctx) == []


def test_sky_qc_does_not_block_chain(tmp_path):
    ctx = _ctx(tmp_path)
    step = SkyQCStep()

    # A cheap optional no-op: never reports already-complete, so the chain
    # always proceeds and re-emits the OK no-op message.
    assert step.is_complete(ctx) is False

    # With Step 1 selection.json present, there are no missing required inputs.
    sel_dir = step1_dir(ctx.result_dir)
    sel_dir.mkdir(parents=True, exist_ok=True)
    (sel_dir / "selection.json").write_text(
        json.dumps({"filenames": ["a.fits"], "target": None}),
        encoding="utf-8",
    )
    assert step.missing_inputs(ctx) == []

    # Running it still returns OK (not blocked, not interactive-blocking).
    result = step.run(ctx)
    assert result.status == StepStatus.OK
