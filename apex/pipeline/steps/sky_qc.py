"""Step 3 (headless): sky preview & QC.

In the GUI, Step 3 is an interactive imexamine-style QC tool: it lets the user
measure FWHM / SNR / magnitudes on the loaded frame, but it persists nothing to
disk — its parameters come from config. In headless mode there is nothing to
compute or write: the sky/QC parameters are taken from the TOML config. This
step is therefore a cheap, optional no-op that simply lets the chain proceed.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths import step1_dir


def _selection_path(result_dir: Path) -> Path:
    return step1_dir(result_dir) / "selection.json"


class SkyQCStep(PipelineStep):
    index = 3
    key = "sky"
    title = "Sky preview & QC"
    interactive = False

    def required_inputs(self, ctx: RunContext) -> List[Path]:
        return [_selection_path(ctx.result_dir)]

    def outputs(self, ctx: RunContext) -> List[Path]:
        return []  # optional QC step; persists nothing by default

    def is_complete(self, ctx: RunContext) -> bool:
        # Cheap optional no-op: never "already complete", so the chain proceeds
        # and the OK no-op message is emitted on each run.
        return False

    def run(self, ctx: RunContext) -> StepResult:
        ctx.logger.info(
            "[sky] headless mode: sky preview/QC parameters are taken from config; "
            "nothing to compute or write."
        )
        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message="sky/QC parameters taken from config (headless no-op; no files written)",
        )
