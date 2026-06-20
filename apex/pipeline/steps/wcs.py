"""Step 5 (headless): WCS plate solving.

Delegates to the Qt-free :func:`apex.analysis.wcs_solve.run_wcs_solve`, the same
orchestration the GUI workers now use. Reads the file list from Step 1's
``selection.json`` and the detections from Step 4, runs the configured solver
(default: the internal Python engine), and writes ``step5_wcs/`` (per-frame FITS
WCS headers, sidecar JSONs, and ``wcs_solve_summary.csv`` / ``frame_wcs_qc.csv``).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List

from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths import (
    step1_dir,
    step2_cropped_dir,
    step4_dir,
    step5_wcs_dir,
    crop_is_active,
)


def _selection_path(result_dir: Path) -> Path:
    return step1_dir(result_dir) / "selection.json"


def _resolve_use_cropped(result_dir: Path) -> bool:
    """Mirror the GUI rule: crop active AND cropped dir holds FITS files."""
    cropped_dir = step2_cropped_dir(result_dir)
    if not crop_is_active(result_dir):
        return False
    if not cropped_dir.exists():
        return False
    return bool(list(cropped_dir.glob("*.fit*")))


class WcsStep(PipelineStep):
    index = 5
    key = "wcs"
    title = "WCS plate solving"
    interactive = False

    def required_inputs(self, ctx: RunContext) -> List[Path]:
        return [step4_dir(ctx.result_dir)]

    def outputs(self, ctx: RunContext) -> List[Path]:
        return [step5_wcs_dir(ctx.result_dir)]

    def is_complete(self, ctx: RunContext) -> bool:
        out = step5_wcs_dir(ctx.result_dir)
        if not out.exists():
            return False
        summary = out / "wcs_solve_summary.csv"
        return summary.exists() and summary.stat().st_size > 0

    def run(self, ctx: RunContext) -> StepResult:
        sel_path = _selection_path(ctx.result_dir)
        if not sel_path.exists():
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=f"missing required input: {sel_path}",
            )

        selection = json.loads(sel_path.read_text(encoding="utf-8"))
        file_list = list(selection.get("filenames", []))
        if not file_list:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message="selection.json contains no filenames",
            )

        use_cropped = _resolve_use_cropped(ctx.result_dir)

        from apex.analysis.wcs_solve import run_wcs_solve, resolve_wcs_engine

        engine = resolve_wcs_engine(ctx.params)

        summary = run_wcs_solve(
            file_list,
            ctx.params,
            ctx.data_dir,
            ctx.result_dir,
            ctx.params.P.cache_dir,
            engine=engine,
            use_cropped=use_cropped,
            logger=ctx.logger,
        )

        out_dir = step5_wcs_dir(ctx.result_dir)
        summary = summary if isinstance(summary, dict) else {}
        n_total = int(summary.get("total", 0) or 0)
        n_ok = int(summary.get("ok", 0) or 0)
        n_qc = int(summary.get("wcs_qc_pass", 0) or 0)
        msg = (
            f"{n_ok}/{len(file_list)} frames solved ({engine}); "
            f"WCS-QC pass={n_qc}"
        )
        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=msg, outputs=[str(out_dir)],
        )
