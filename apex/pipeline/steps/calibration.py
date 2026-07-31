"""Step 0 (headless): detector calibration (bias/dark/flat).

Off-chain OPTIONAL pre-stage. It is intentionally NOT part of the shared 1-7
chain returned by :func:`apex.pipeline.registry.get_steps` (so the runner's
1-based index math is untouched); it is exposed separately via
:func:`apex.pipeline.registry.get_calibration_step`.

The scan, the matching and the calibration all live in the Qt-free
:mod:`apex.analysis.calibration_scan` / :mod:`apex.analysis.calibration_run`,
which the GUI Step-0 window drives too — the same code path, so a headless run
and the window cannot produce different frames.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

from apex.analysis import calibration_scan as scan
from apex.analysis.calibration import CalibrationOptions
from apex.analysis.calibration_run import ALL_NIGHTS, run_calibration
from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths import step0_calibration_dir


class CalibrationStep(PipelineStep):
    index = 0
    key = "calibration"
    title = "Detector Calibration"
    interactive = False

    def outputs(self, ctx: RunContext) -> List[Path]:
        return [step0_calibration_dir(ctx.result_dir) / "calibration.json"]

    def run(self, ctx: RunContext) -> StepResult:
        P = getattr(ctx.params, "P", ctx.params)
        opts = CalibrationOptions.from_mapping(getattr(P, "calibration", None))
        log = ctx.logger.info

        roots = [str(ctx.data_dir)]
        # A shared bias/dark library usually lives outside data_dir.
        for extra in (getattr(P, "calibration_extra_dirs", None) or []):
            if extra:
                roots.append(str(extra))

        tz = getattr(P, "site_tz_offset_hours", None)
        frames: List[scan.FrameInfo] = []
        for root in roots:
            if root and Path(root).exists():
                frames.extend(scan.scan_folder(root, tz_offset_hours=tz, warn=log))

        if not any(f.ftype == "light" for f in frames):
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.SKIPPED,
                message=f"no light frames found under {roots}",
            )
        if not any(f.ftype in ("bias", "dark", "flat") for f in frames):
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.SKIPPED,
                message="no bias/dark/flat frames found — nothing to calibrate with",
            )

        out_dir = step0_calibration_dir(ctx.result_dir)
        summary = run_calibration(frames, ALL_NIGHTS, out_dir, opts, log=log)
        note = ""
        if summary.get("n_missing_flat"):
            note += f", {summary['n_missing_flat']} without a flat"
        if summary.get("n_temp_mismatch"):
            note += (f", {summary['n_temp_mismatch']} with a dark outside "
                     f"the {opts.temp_match_tol_c:g}°C tolerance")
        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=(f"calibrated {summary['n_calibrated']} frame(s) across "
                     f"{summary['n_nights']} night(s){note} → "
                     f"{summary['calibrated_root']}"),
            outputs=[str(p) for p in self.outputs(ctx)],
        )
