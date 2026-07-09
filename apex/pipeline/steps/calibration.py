"""Step 0 (headless stub): detector calibration (bias/dark/flat).

Off-chain OPTIONAL pre-stage. It is intentionally NOT part of the shared 1-7
chain returned by :func:`apex.pipeline.registry.get_steps` (so the runner's
1-based index math is untouched); it is exposed separately via
:func:`apex.pipeline.registry.get_calibration_step`.

The calibration algorithm lives in the Qt-free :mod:`apex.analysis.calibration`
and is driven interactively by the GUI Step-0 window. A headless implementation
of :meth:`run` (build masters + apply from config) can be added later; for now
it reports ``not_implemented`` so a GUI-produced ``calibration.json`` is still
recognised by :meth:`is_complete`.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths import step0_calibration_dir


class CalibrationStep(PipelineStep):
    index = 0
    key = "calibration"
    title = "Detector Calibration"
    interactive = True

    def outputs(self, ctx: RunContext) -> List[Path]:
        return [step0_calibration_dir(ctx.result_dir) / "calibration.json"]

    def run(self, ctx: RunContext) -> StepResult:
        return StepResult(
            index=self.index,
            key=self.key,
            status=StepStatus.NOT_IMPLEMENTED,
            message=("detector calibration runs in the GUI Step-0 window "
                     "(apex.analysis.calibration); headless port pending"),
        )
