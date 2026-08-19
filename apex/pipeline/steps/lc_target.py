"""LC Step 8 (headless): resolve the target and its comparison ensemble.

The LC branch stopped at Step 7 — the pipeline registry lists steps 1-7 for LC
and nothing after. Not because the science needs a window: fourteen Qt-free
services in `apex.analysis.light_curve` do the building, detrending and period
work already. What stopped it is that nothing in the config could say which star
the light curve is of.

This step is that gate, and it blocks rather than guesses, the same way the
isochrone step does. A batch that produces nothing beats one that produces a
clean light curve of the wrong object — because a light curve of the wrong
object does not look wrong.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import List

import pandas as pd

from apex.analysis.light_curve.target_config import (
    missing_target_settings, read_target, resolve_comparisons, write_selection,
)
from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths import step6_refbuild_dir
from apex.utils.step_paths_lc import step8_selection_dir

MASTER_TABLE = "master_sources.csv"


class LcTargetStep(PipelineStep):
    index = 8
    key = "lctarget"
    name = "Light-curve target"

    def inputs(self, ctx: RunContext) -> List[Path]:
        return [step6_refbuild_dir(ctx.result_dir) / MASTER_TABLE]

    def outputs(self, ctx: RunContext) -> List[Path]:
        return [step8_selection_dir(ctx.result_dir) / "lc_target_selection.json"]

    def is_complete(self, ctx: RunContext) -> bool:
        return self.outputs(ctx)[0].exists()

    def run(self, ctx: RunContext) -> StepResult:
        missing = missing_target_settings(ctx.params)
        if missing:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=("a light curve needs to know which star it is of, and "
                         "no default is defensible: " + "; ".join(missing)),
            )

        master = step6_refbuild_dir(ctx.result_dir) / MASTER_TABLE
        if not master.exists():
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=f"no master catalog at {master}",
            )

        started = time.perf_counter()
        target = read_target(ctx.params)
        try:
            catalog = pd.read_csv(master)
        except Exception as exc:                    # noqa: BLE001
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.FAILED,
                message=f"could not read the master catalog: {exc}",
            )

        known = set(pd.to_numeric(catalog.get("ID"), errors="coerce").dropna().astype(int))
        if target.target_id not in known:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=(f"lightcurve.target_id={target.target_id} is not in the "
                         f"master catalog ({len(known)} stars)"),
            )

        comparisons = resolve_comparisons(target, catalog)
        if not comparisons:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=("no comparison stars — set lightcurve.comparison_ids, "
                         "or check that the master catalog has more than the target"),
            )

        path = write_selection(ctx.result_dir, target, comparisons)
        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=(f"target ID {target.target_id}"
                     + (f" ({target.target_name})" if target.target_name else "")
                     + f"; {len(comparisons)} comparisons [{target.comparison_mode}]"),
            duration_s=time.perf_counter() - started,
            outputs=[str(path)],
        )
