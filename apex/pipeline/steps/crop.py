"""Step 2 (headless): image crop (config-driven).

The GUI step gets the crop rectangle from a mouse selection; headless reads it
from config (``[crop]`` in parameters.toml) or from an existing
``step2_crop/crop_rect.json`` (the same schema the GUI writes). Delegates the
actual FITS cropping to the Qt-free :func:`apex.analysis.crop.run_crop`, the
same compute the GUI worker now uses.

If no crop is configured this step writes a skip marker
(``crop_rect.json`` with ``skipped: true``) and returns OK so the chain
proceeds on the original (uncropped) frames.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List, Optional, Tuple

from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths import step1_dir, step2_dir, crop_rect_path


def _selection_path(result_dir: Path) -> Path:
    return step1_dir(result_dir) / "selection.json"


def _resolve_crop_rect(ctx: RunContext) -> Optional[Tuple[int, int, int, int, str]]:
    """Resolve the crop rectangle for headless execution.

    Priority:
      1. An existing ``step2_crop/crop_rect.json`` (same schema the GUI writes);
         a file with ``skipped: true`` means "no crop".
      2. The ``[crop]`` config section (``enable`` + ``x0/y0/x1/y1``).

    Returns ``(x0, y0, x1, y1, ref_file)`` or ``None`` when crop is skipped /
    unconfigured.
    """
    # 1. Honor an existing crop_rect.json (GUI-produced or prior headless run).
    rect_path = crop_rect_path(ctx.result_dir)
    if rect_path.exists():
        try:
            data = json.loads(rect_path.read_text(encoding="utf-8"))
        except Exception:
            data = None
        if isinstance(data, dict):
            if data.get("skipped"):
                return None
            try:
                x0 = int(data["x0"])
                y0 = int(data["y0"])
                x1 = int(data["x1"])
                y1 = int(data["y1"])
                return x0, y0, x1, y1, str(data.get("ref_file", "") or "")
            except (KeyError, TypeError, ValueError):
                pass

    # 2. Config-driven crop rectangle.
    P = ctx.params.P
    if not bool(getattr(P, "crop_enable", False)):
        return None
    x0 = getattr(P, "crop_x0", None)
    y0 = getattr(P, "crop_y0", None)
    x1 = getattr(P, "crop_x1", None)
    y1 = getattr(P, "crop_y1", None)
    if None in (x0, y0, x1, y1):
        return None
    return int(x0), int(y0), int(x1), int(y1), ""


def _write_skip_marker(result_dir: Path, ref_file: str = "") -> Path:
    """Mirror the GUI skip state: a crop_rect.json with ``skipped: true``."""
    step2_out = step2_dir(result_dir)
    step2_out.mkdir(parents=True, exist_ok=True)
    rect_path = crop_rect_path(result_dir)
    rect_path.write_text(
        json.dumps({"skipped": True, "ref_file": ref_file}, indent=2),
        encoding="utf-8",
    )
    return rect_path


def _write_crop_marker(result_dir: Path, x0: int, y0: int, x1: int, y1: int,
                       ref_file: str = "") -> Path:
    """Mirror the GUI ``save_crop_info``: crop_rect.json with the rectangle."""
    step2_out = step2_dir(result_dir)
    step2_out.mkdir(parents=True, exist_ok=True)
    rect_path = crop_rect_path(result_dir)
    rect_path.write_text(
        json.dumps(
            {
                "x0": int(x0),
                "y0": int(y0),
                "x1": int(x1),
                "y1": int(y1),
                "ref_file": ref_file,
                "skipped": False,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return rect_path


class CropStep(PipelineStep):
    index = 2
    key = "crop"
    title = "Image crop"
    interactive = False  # config-driven in headless mode

    def required_inputs(self, ctx: RunContext) -> List[Path]:
        return [_selection_path(ctx.result_dir)]

    def outputs(self, ctx: RunContext) -> List[Path]:
        return [crop_rect_path(ctx.result_dir)]

    def is_complete(self, ctx: RunContext) -> bool:
        return crop_rect_path(ctx.result_dir).exists()

    def run(self, ctx: RunContext) -> StepResult:
        sel_path = _selection_path(ctx.result_dir)
        if not sel_path.exists():
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=f"missing required input: {sel_path}",
            )

        rect = _resolve_crop_rect(ctx)
        if rect is None:
            marker = _write_skip_marker(ctx.result_dir)
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.OK,
                message="crop skipped (no [crop] config); using original frames",
                outputs=[str(marker)],
            )

        x0, y0, x1, y1, ref_file = rect

        selection = json.loads(sel_path.read_text(encoding="utf-8"))
        file_list = list(selection.get("filenames", []))
        if not file_list:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message="selection.json contains no filenames",
            )

        from apex.analysis.crop import run_crop

        n_total, errors = run_crop(
            file_list, ctx.params, x0, y0, x1, y1, logger=ctx.logger,
        )

        if errors:
            first_err = errors[0][1]
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.FAILED,
                message=f"{len(errors)}/{n_total} frames failed to crop; "
                        f"first error: {first_err}",
            )

        marker = _write_crop_marker(ctx.result_dir, x0, y0, x1, y1, ref_file)
        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=f"cropped {n_total} frames to ({x0},{y0})-({x1},{y1})",
            outputs=[str(marker)],
        )
