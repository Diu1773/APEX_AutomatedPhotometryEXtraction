"""A run has to record what it was made with.

Measured on 2026-08-17: the same code on byte-identical input gave a different
answer for five of 22,305 crowded measurements (<=5e-5 mag) because one machine
had scipy 1.18.0 and the other 1.17.1. Everything else — numpy, astropy,
photutils, pandas, sep — matched. The manifest recorded the steps and their
timings and nothing about the environment, so explaining the difference took a
full control run and a package-by-package diff.

Every measurement that moved was already flagged `CROWDING_UNRELIABLE` or
`NONCONVERGENCE`, so nothing published changed. The lesson is not about that
number, it is that a run which cannot say what produced it cannot be reproduced
— which is exactly the qualification the manuscript claim matrix puts on the
word "reproducible".
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

from apex.pipeline.runner import PipelineRunner


def test_the_environment_block_names_the_numerical_stack():
    env = PipelineRunner._environment()
    assert env["python"]
    assert env["platform"]
    packages = env["packages"]
    # The ones that can change an answer.
    for name in ("numpy", "scipy", "astropy", "photutils", "pandas"):
        assert name in packages, name
        assert packages[name], f"{name} 버전이 비어 있다"


def test_a_written_manifest_carries_it(tmp_path, monkeypatch):
    from apex.pipeline.base import StepStatus
    from apex.pipeline.context import RunContext

    class _Report:
        def to_dict(self):
            return {"mode": "cmd", "steps": [], "success": True}

    params = type("P", (), {"P": type("Inner", (), {})()})()
    ctx = RunContext(mode="cmd", params=params,
                     result_dir=tmp_path, data_dir=tmp_path,
                     logger=logging.getLogger("test"))
    path = PipelineRunner._write_manifest(ctx, _Report())
    assert path is not None

    written = json.loads(Path(path).read_text(encoding="utf-8"))
    assert "environment" in written, "실행이 자기 환경을 안 남긴다"
    assert written["environment"]["packages"]["scipy"]
    assert StepStatus  # imported to keep the pipeline contract in view
