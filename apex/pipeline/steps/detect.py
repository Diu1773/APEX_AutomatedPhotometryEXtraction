"""Step 4 (headless): source detection.

Delegates to the Qt-free :func:`apex.analysis.detection.run_detection`, the same
compute used by the GUI ``DetectionWorker``. Reads the file list from Step 1's
``selection.json`` and writes ``step4_detection/detect_*.{json,csv}``.
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
    crop_is_active,
)
from apex.utils.filter_parameters import normalize_filter_float_map


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


def _build_filter_sigma_map(P) -> dict:
    """Mirror DetectionWindow.load_filter_sigma_map (Qt-free subset).

    Reads the configured ``detect_sigma_by_filter`` map and any legacy
    ``detect_sigma_<filter>`` keys from the raw TOML, normalized to canonical
    filter keys. Per-filter expansion to active filters is unnecessary: missing
    keys fall back to ``detect_sigma`` inside run_detection anyway.
    """
    raw = getattr(P, "_raw", {}) or {}

    sigma_by_filter = getattr(P, "detect_sigma_by_filter", None)
    if not isinstance(sigma_by_filter, dict):
        sigma_by_filter = raw.get("detect_sigma_by_filter", {})
    filter_sigma_map = dict(normalize_filter_float_map(sigma_by_filter))

    for key, val in raw.items():
        if key.startswith("detect_sigma_") and key not in {
            "detect_sigma",
            "detect_sigma_by_filter",
            "detect_sigma_g",
            "detect_sigma_r",
            "detect_sigma_i",
        }:
            filt = key.replace("detect_sigma_", "")
            filter_sigma_map.update(normalize_filter_float_map({filt: val}))

    return filter_sigma_map


class DetectStep(PipelineStep):
    index = 4
    key = "detect"
    title = "Source detection"
    interactive = False

    def required_inputs(self, ctx: RunContext) -> List[Path]:
        return [_selection_path(ctx.result_dir)]

    def outputs(self, ctx: RunContext) -> List[Path]:
        return [step4_dir(ctx.result_dir)]

    def is_complete(self, ctx: RunContext) -> bool:
        out = step4_dir(ctx.result_dir)
        if not out.exists():
            return False
        return any(out.glob("detect_*.json"))

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
        filter_sigma_map = _build_filter_sigma_map(ctx.params.P)

        from apex.analysis.detection import run_detection

        summary = run_detection(
            file_list,
            ctx.params,
            ctx.data_dir,
            ctx.params.P.cache_dir,
            use_cropped,
            filter_sigma_map,
            logger=ctx.logger,
        )

        out_dir = step4_dir(ctx.result_dir)
        total_sources = int(summary.get("total_sources", 0))
        total_files = int(summary.get("total_files", 0))
        median_fwhm_px = summary.get("median_fwhm_px")
        msg = (
            f"{total_files}/{len(file_list)} frames detected; "
            f"{total_sources} sources total"
        )
        if isinstance(median_fwhm_px, (int, float)) and median_fwhm_px == median_fwhm_px:
            msg += f"; median FWHM {float(median_fwhm_px):.2f}px"
        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=msg, outputs=[str(out_dir)],
        )
