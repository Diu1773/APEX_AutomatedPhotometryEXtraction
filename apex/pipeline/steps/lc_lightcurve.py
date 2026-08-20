"""LC Step 9 (headless): the raw differential light curve.

Step 8 wrote down which star this is of. This builds the curve — target minus
comparison ensemble, frame by frame — and it does so by inheriting the window's
build rather than reimplementing it: `HeadlessLightCurveBuilder` supplies the
state a window would have set from widgets, and the calculation underneath is
`RawLightCurveBuilder._build_light_curve_core`, the same object the window runs.

That is why there is so little here. The port was moving 1,100 lines out of a
`QMainWindow` subclass; what remains at this layer is reading a selection off
disk and refusing to guess when it is absent.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import List

from apex.analysis.light_curve.target_config import read_selection
from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths_lc import step8_selection_dir, step9_lc_dir


def _diagnose_empty(ctx, builder, target_id: int, summary: dict) -> str:
    """Say why a curve came out empty, rather than guessing at it.

    An empty curve is the one failure here that does not look like one — the
    file is written, the row count is right, and every magnitude is blank. The
    common cause is not a missing star: it is a workspace from 2025, whose
    photometry tables identify sources by `det_uid` and carry neither `ID` nor
    `source_id`. Nothing downstream can match a star in those, and the old
    pipeline bridged the gap with a positional cross-match this one does not do.

    Reproduced on YZ Boo 2025-04-29: 77 rows, 0 usable points, no error (D-013).
    """
    n_total = int(summary.get("n_total", 0))
    columns: set[str] = set()
    try:
        index = builder._load_active_photometry_index(ctx.result_dir)
        frames = [str(v) for v in index.get("file", [])][:3]
        for frame in frames:
            table = builder._get_photometry_df(ctx.result_dir, frame)
            if table is not None and not table.empty:
                columns |= set(table.columns)
    except Exception:                                   # noqa: BLE001
        pass

    if columns and not ({"ID", "source_id"} & columns):
        return (f"built {n_total} rows and none carry a measurement, because the "
                f"photometry tables identify sources by det_uid and have neither "
                f"an ID nor a source_id column. That is a workspace from before "
                f"the current Step 7 — open it once in the GUI, which upgrades "
                f"it, or re-run Steps 1-7 headless.")
    return (f"built {n_total} rows but none had both the target and its "
            f"ensemble — check that ID {target_id} and the comparisons are in "
            f"the forced photometry")


class LcLightCurveStep(PipelineStep):
    index = 9
    key = "lclightcurve"
    name = "Light curve"

    def inputs(self, ctx: RunContext) -> List[Path]:
        return [step8_selection_dir(ctx.result_dir) / "lc_target_selection.json"]

    def outputs(self, ctx: RunContext) -> List[Path]:
        selection = read_selection(ctx.result_dir) or {}
        target_id = selection.get("target_id")
        if target_id is None:
            return []
        return [step9_lc_dir(ctx.result_dir) / f"lightcurve_ID{int(target_id)}_raw.csv"]

    def is_complete(self, ctx: RunContext) -> bool:
        out = self.outputs(ctx)
        return bool(out) and out[0].exists()

    def run(self, ctx: RunContext) -> StepResult:
        selection = read_selection(ctx.result_dir)
        if not selection:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=("no target selection — run LC Step 8, or set "
                         "lightcurve.target_id in the config"),
            )

        target_id = selection.get("target_id")
        comp_ids = [int(v) for v in (selection.get("comparison_ids") or [])]
        if target_id is None or int(target_id) <= 0:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=f"the selection names no target: {selection}",
            )
        if not comp_ids:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=("the selection names no comparison stars — a "
                         "differential light curve needs an ensemble"),
            )

        from apex.analysis.light_curve.raw_lightcurve import HeadlessLightCurveBuilder

        started = time.perf_counter()
        builder = HeadlessLightCurveBuilder(
            ctx.params, [ctx.result_dir],
            logger=getattr(ctx, "log", None),
            project_state=getattr(ctx, "project_state", None),
            comp_candidate_ids=comp_ids,
        )
        try:
            summary = builder._build_light_curve_core(int(target_id), comp_ids)
        except Exception as exc:                        # noqa: BLE001
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.FAILED,
                message=f"the light-curve build failed: {exc}",
                duration_s=time.perf_counter() - started,
            )
        elapsed = time.perf_counter() - started

        n_valid = int(summary.get("n_valid", 0))
        if n_valid == 0:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=_diagnose_empty(ctx, builder, target_id, summary),
                duration_s=elapsed,
            )

        written = [p for p in step9_lc_dir(ctx.result_dir).glob("lightcurve_*.csv")]
        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=(f"target ID {target_id}, {len(comp_ids)} comparisons; "
                     f"{n_valid}/{summary.get('n_total', 0)} points"),
            duration_s=elapsed,
            outputs=[str(p) for p in sorted(written)],
        )
