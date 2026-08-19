"""CMD Step 11 (headless): the colour-magnitude figure.

This was deferred with the reason "nothing to port": the window loads the Step
10 table and opens an interactive viewer, so there was no calculation sitting
behind a `QThread` waiting to be freed.

That reason held for the viewer and not for the figure. The viewer is an
instrument — region selection, parallax and proper-motion sliders, quality
masks, a Teff colour bar — and translating it to a batch run would produce a
worse viewer, not a figure. But a run that measures a cluster should leave the
picture of it, and the picture needs none of that machinery: the same table, the
same axes, everything that passed.

So this is a new export rather than a move, which is why it took a decision
rather than a refactor.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import List

from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths_cmd import step10_zp_dir, step11_cmd_dir

CMD_TABLE = "median_by_ID_filter_wide_cmd.csv"


class CmdPlotStep(PipelineStep):
    index = 11
    key = "cmdplot"
    name = "CMD plot"

    def inputs(self, ctx: RunContext) -> List[Path]:
        return [step10_zp_dir(ctx.result_dir) / CMD_TABLE]

    def outputs(self, ctx: RunContext) -> List[Path]:
        return [step11_cmd_dir(ctx.result_dir) / "step11_cmd.png"]

    def is_complete(self, ctx: RunContext) -> bool:
        return self.outputs(ctx)[0].exists()

    def run(self, ctx: RunContext) -> StepResult:
        table = step10_zp_dir(ctx.result_dir) / CMD_TABLE
        if not table.exists():
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=f"no Step 10 table at {table}",
            )

        from apex.analysis.cmd.cmd_plot import export_cmd_plot

        started = time.perf_counter()
        written = export_cmd_plot(ctx.result_dir, ctx.params)
        elapsed = time.perf_counter() - started

        if not written:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=("could not draw a CMD — fewer than two calibrated "
                         "bands, or no finite colour/magnitude pairs"),
                duration_s=elapsed,
            )
        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=f"{len(written)} CMD figure", duration_s=elapsed,
            outputs=[str(p) for p in written],
        )
