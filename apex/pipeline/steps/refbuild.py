"""Step 6 (headless): master catalog build.

Delegates to the Qt-free :func:`apex.analysis.refbuild.run_refbuild`, the same
compute used by the GUI ``RefBuildWorker``. Reads the file list from Step 1's
``selection.json`` (filtered by Step 4 frame QC, mirroring the GUI), the
detections from Step 4, and the WCS / Gaia catalog from Step 5, and writes the
master catalog TSVs plus ``ref_frame_stats.csv`` / ``ref_build_meta.json`` under
``step6_refbuild/``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List

from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths import (
    step1_dir,
    step4_dir,
    step5_wcs_dir,
    step6_refbuild_dir,
)
from apex.utils.common_helpers import normalize_filter_key
from apex.utils.qc_utils import filter_files_by_qc


def _selection_path(result_dir: Path) -> Path:
    return step1_dir(result_dir) / "selection.json"


class RefBuildStep(PipelineStep):
    index = 6
    key = "refbuild"
    title = "Master catalog build"
    interactive = False

    def required_inputs(self, ctx: RunContext) -> List[Path]:
        return [
            step5_wcs_dir(ctx.result_dir),
            step4_dir(ctx.result_dir),
        ]

    def outputs(self, ctx: RunContext) -> List[Path]:
        return [step6_refbuild_dir(ctx.result_dir)]

    def is_complete(self, ctx: RunContext) -> bool:
        out = step6_refbuild_dir(ctx.result_dir)
        if not out.exists():
            return False
        cat = out / "ref_catalog.tsv"
        return cat.exists() and cat.stat().st_size > 0

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

        P = ctx.params.P

        # Mirror RefBuildWindow.run_ref_build: optionally drop frames that fail
        # the Step 4 frame QC before building the master catalog.
        use_qc = bool(getattr(P, "wcs_require_qc_pass", True))
        file_list, qc_info = filter_files_by_qc(
            Path(ctx.result_dir), file_list, require_qc=use_qc
        )
        if not file_list:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message="no frames remaining after Step 4 frame QC filter",
            )

        ref_filter = normalize_filter_key(getattr(P, "global_ref_filter", ""))

        from apex.analysis.refbuild import run_refbuild

        summary = run_refbuild(
            params=ctx.params,
            data_dir=ctx.data_dir,
            result_dir=ctx.result_dir,
            cache_dir=P.cache_dir,
            file_list=file_list,
            ref_filter=ref_filter,
            sat_drop_pct=float(getattr(P, "ref_select_sat_pct", 20.0)),
            elong_drop_pct=float(getattr(P, "ref_select_elong_pct", 20.0)),
            ref_cat_max_sources=int(getattr(P, "ref_cat_max_sources", 0)),
            ref_cat_min_sources=int(getattr(P, "ref_cat_min_sources", 50)),
            ref_cat_max_elong=float(getattr(P, "ref_cat_max_elong", 1.5)),
            ref_cat_max_abs_round=float(getattr(P, "ref_cat_max_abs_round", 0.4)),
            ref_cat_sharp_min=float(getattr(P, "ref_cat_sharp_min", 0.2)),
            ref_cat_sharp_max=float(getattr(P, "ref_cat_sharp_max", 1.0)),
            ref_cat_min_peak_adu=float(getattr(P, "ref_cat_min_peak_adu", 0.0)),
            wcs_match_radius_arcsec=float(getattr(P, "ref_wcs_match_radius_arcsec", 2.0)),
            wcs_min_match_rate=float(getattr(P, "ref_wcs_min_match_rate", 0.2)),
            wcs_min_match_n=int(getattr(P, "ref_wcs_min_match_n", 50)),
            wcs_max_sep_med_arcsec=float(getattr(P, "ref_wcs_max_sep_med_arcsec", 1.5)),
            wcs_max_sep_p90_arcsec=float(getattr(P, "ref_wcs_max_sep_p90_arcsec", 2.5)),
            wcs_max_dup_rate=float(getattr(P, "ref_wcs_max_dup_rate", 0.1)),
            ref_per_date=bool(getattr(P, "ref_per_date", True)),
            ref_build_mode=str(getattr(P, "ref_build_mode", "hybrid")),
            gaia_mag_limit=float(
                getattr(P, "idmatch_gaia_g_limit", getattr(P, "gaia_mag_max", 18.0))
            ),
            ref_master_union=bool(getattr(P, "ref_master_union", True)),
            ref_union_min_frames=int(getattr(P, "ref_union_min_frames", 1)),
            logger=ctx.logger,
        )

        out_dir = step6_refbuild_dir(ctx.result_dir)
        summary = summary if isinstance(summary, dict) else {}
        n_sources = int(summary.get("n_sources", 0) or 0)
        ref_frame = summary.get("ref_frame", "")
        ref_filter_used = summary.get("ref_filter", "")
        msg = (
            f"master catalog built: {n_sources} sources "
            f"(ref={ref_frame}, filter={ref_filter_used})"
        )
        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=msg, outputs=[str(out_dir)],
        )
