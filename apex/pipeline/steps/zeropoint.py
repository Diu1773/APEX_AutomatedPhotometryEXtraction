"""CMD Step 10 (headless): zero-point calibration and the wide CMD table.

Same shape as Step 8 — the compute is a ``QThread`` in the GUI module, driven
synchronously, guarded on PyQt5 being importable (it is an optional extra).

Step 10 does not need Steps 8 or 9. It reads star IDs straight from the Step 7
tables and the Step 6 master catalogue, and *switches* to PSF magnitudes the
moment a valid Step 8 signature exists. That switch is silent and it changes
which numbers the products carry, so this step records which source was used
and copies the previous solution aside before overwriting it.
"""

from __future__ import annotations

import csv
import shutil
import time
from pathlib import Path
from typing import List

from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths import step6_refbuild_dir, step7_forced_phot_dir
from apex.utils.step_paths_cmd import step10_zp_dir

# Files worth keeping a copy of before a re-run overwrites them. Learned the
# hard way on 2026-08-15: a re-run with a config still pointing at the source
# workspace overwrote M13's published zero-point, and only a backup made
# recovery checkable.
BACKUP_NAMES = (
    "frame_zeropoint.csv",
    "median_by_ID_filter.csv",
    "median_by_ID_filter_wide.csv",
    "median_by_ID_filter_wide_cmd.csv",
    "zp_fit_coefficients.csv",
)


def _qt_available() -> bool:
    try:
        import PyQt5.QtCore  # noqa: F401
    except Exception:  # noqa: BLE001
        return False
    return True


def _photometry_source(zp_dir: Path) -> str:
    """Whether the solution that just landed used PSF or aperture magnitudes.

    Read from the named column, not by looking for ``,psf,`` in the file: the
    column sits in the middle today but would go unmatched at the end of a
    line, and ``ref_source``/``mag_input_column`` carry values that contain
    the same word.
    """
    coefficients = zp_dir / "zp_fit_coefficients.csv"
    if not coefficients.exists():
        return "unknown"
    with coefficients.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        rows = csv.DictReader(handle)
        values = {
            str(row.get("photometry_source", "")).strip()
            for row in rows
            if str(row.get("photometry_source", "")).strip()
        }
    if not values:
        return "unknown"
    # One solution per filter; they should agree, and if they do not that is
    # the interesting thing to report rather than whichever came first.
    return values.pop() if len(values) == 1 else "mixed(" + "+".join(sorted(values)) + ")"


class ZeropointStep(PipelineStep):
    index = 10
    key = "zeropoint"
    title = "Zero-point calibration (CMD)"
    interactive = False

    def required_inputs(self, ctx: RunContext) -> List[Path]:
        return [step7_forced_phot_dir(ctx.result_dir), step6_refbuild_dir(ctx.result_dir)]

    def outputs(self, ctx: RunContext) -> List[Path]:
        return [step10_zp_dir(ctx.result_dir)]

    def is_complete(self, ctx: RunContext) -> bool:
        # The wide CMD table is what Steps 11 and 12 actually read.
        return (step10_zp_dir(ctx.result_dir)
                / "median_by_ID_filter_wide_cmd.csv").exists()

    def run(self, ctx: RunContext) -> StepResult:
        if not _qt_available():
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.NOT_IMPLEMENTED,
                message=("zero-point calibration still runs through the GUI "
                         "module's worker, which needs PyQt5. Install the extra "
                         "(pip install 'apex[gui]') or run this step in the app."),
            )

        zp_dir = step10_zp_dir(ctx.result_dir)
        if zp_dir.exists():
            backup = zp_dir / "_pre_run_backup"
            backup.mkdir(exist_ok=True)
            for name in BACKUP_NAMES:
                source = zp_dir / name
                if source.exists():
                    shutil.copyfile(source, backup / name)

        from PyQt5.QtCore import QCoreApplication

        QCoreApplication.instance() or QCoreApplication([])

        from apex.gui.workflow.cmd.step10_zeropoint_calibration import (
            ZeropointCalibrationWorker,
        )

        params = ctx.params
        worker = ZeropointCalibrationWorker(
            params, params.P.data_dir, params.P.result_dir, params.P.cache_dir)
        if ctx.logger is not None:
            worker._log = lambda message: ctx.logger.info("%s", message)

        started = time.perf_counter()
        worker.run()
        elapsed = time.perf_counter() - started

        summary = dict(getattr(worker, "last_summary", {}) or {})
        if not bool(summary.get("ok", False)):
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.FAILED,
                message=str(getattr(worker, "last_error", "")
                            or "zero-point calibration reported failure"),
                outputs=[str(zp_dir)], duration_s=elapsed,
            )

        source = _photometry_source(zp_dir)
        message = f"zero-points solved from {source} photometry"
        if source == "unknown":
            message = "zero-points solved; photometry source not recorded"
        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=message, outputs=[str(zp_dir)], duration_s=elapsed,
        )
