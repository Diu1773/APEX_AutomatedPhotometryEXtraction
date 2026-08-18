"""Step 7 (headless): forced aperture photometry.

Delegates to the Qt-free :func:`apex.analysis.forced_photometry.run_forced_photometry`,
the same compute used by the GUI ``ForcedPhotWorker``. Reads the file list from
Step 1's ``selection.json``, the master catalog from Step 6, and the WCS from
Step 5 (FITS headers / ASTAP sidecars), and writes
``step7_forced_phot/photometry_*.tsv`` plus the summary CSVs.
"""

from __future__ import annotations

import os

import json
from pathlib import Path
from typing import List

from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths import (
    step1_dir,
    step5_wcs_dir,
    step6_refbuild_dir,
    step7_forced_phot_dir,
)


def _selection_path(result_dir: Path) -> Path:
    return step1_dir(result_dir) / "selection.json"


class ForcedPhotStep(PipelineStep):
    index = 7
    key = "forcedphot"
    title = "Forced aperture photometry"
    interactive = False

    def required_inputs(self, ctx: RunContext) -> List[Path]:
        return [
            step6_refbuild_dir(ctx.result_dir),
            step5_wcs_dir(ctx.result_dir),
        ]

    def outputs(self, ctx: RunContext) -> List[Path]:
        return [step7_forced_phot_dir(ctx.result_dir)]

    def is_complete(self, ctx: RunContext) -> bool:
        out = step7_forced_phot_dir(ctx.result_dir)
        if not out.exists():
            return False
        return (out / "photometry_index.csv").exists()

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

        from apex.analysis.forced_photometry import run_forced_photometry

        # Processes by default *here* but not in the GUI. Forced photometry is
        # interpreter-lock bound (83 % CPU on one worker, 109 % on twelve), so
        # threads cost 4.2x wall time and buy nothing; processes cut Step 7 from
        # 244.2 s to 69.6 s with byte-identical per-frame tables. The GUI keeps
        # threads because its per-worker progress callbacks cannot cross a
        # process boundary — a headless run has no such display to lose.
        # Set APEX_FORCEDPHOT_PROCESSES=0 to force the thread path.
        use_processes = os.environ.get(
            "APEX_FORCEDPHOT_PROCESSES", "1").strip().lower() not in ("0", "false", "no", "off")

        # The growth curve is the one Step 7 product that never reaches
        # disk: the run hands it to `apcorr_cb` and the window plots it
        # live. Catch it here so a batch run can draw the same figure —
        # otherwise the curve the aperture correction is read off exists
        # only on somebody's screen.
        growth_curves: dict = {}

        def _catch_growth_curve(gc):
            try:
                name = str(gc.get('fname') or gc.get('file') or '')
                if name and gc.get('radii_px'):
                    growth_curves[name] = dict(gc)
            except Exception:  # noqa: BLE001 - QC must not disturb the run
                pass

        summary = run_forced_photometry(
            file_list,
            ctx.params,
            ctx.data_dir,
            ctx.params.P.cache_dir,
            result_dir=ctx.result_dir,
            logger=ctx.logger,
            use_processes=use_processes,
            apcorr_cb=_catch_growth_curve,
        )

        out_dir = step7_forced_phot_dir(ctx.result_dir)
        index_rows = summary.get("index_rows", []) if isinstance(summary, dict) else []
        n_ok = sum(1 for r in index_rows if r.get("status") == "ok")
        n_det = sum(int(r.get("n_detected", 0) or 0) for r in index_rows)
        n_forced = sum(int(r.get("n_forced", 0) or 0) for r in index_rows)
        msg = (
            f"{n_ok}/{len(file_list)} frames photometered; "
            f"detected={n_det} forced={n_forced}"
        )

        # Second-stage photometric QC (transparency offsets from matched
        # bright stars) — writes phot_quality.csv; never blocks the step.
        if n_ok > 0:
            try:
                from apex.analysis.photometric_qc import (
                    run_photometric_qc,
                    summarize_photometric_qc,
                )

                qc_df = run_photometric_qc(ctx.result_dir)
                if not qc_df.empty:
                    counts = summarize_photometric_qc(qc_df)
                    msg += (
                        f"; transparency QC PASS {counts['PASS']}"
                        f"/REVIEW {counts['REVIEW']}/FAIL {counts['FAIL']}"
                        f"/SKIP {counts['SKIP']}"
                    )
            except Exception as exc:  # noqa: BLE001 - QC must not fail the step
                if ctx.logger is not None:
                    ctx.logger.warning("photometric QC skipped: %s", exc)

        try:
            from apex.analysis.forcedphot_qc import (
                export_forcedphot_qc, one_per_filter,
            )
            qc_figures = export_forcedphot_qc(
                ctx.result_dir, ctx.params, one_per_filter(growth_curves))
            if qc_figures:
                msg += f'; {len(qc_figures)} QC figures'
        except Exception:  # noqa: BLE001 - QC must not fail finished photometry
            if ctx.logger is not None:
                ctx.logger.exception('Could not write Step 7 QC figures')

        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=msg, outputs=[str(out_dir)],
        )
