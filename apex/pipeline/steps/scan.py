"""Step 1 (headless): scan FITS headers and resolve the target.

This is the first workflow step ported to headless execution. It reuses the
Qt-free :class:`~apex.core.file_manager.FileManager` (which writes
``step1_file_selection/headers.csv``) and the pure ``select_header_target``
helper. The target is taken from the config when supplied, otherwise inferred
from the FITS headers. SIMBAD name-only resolution is intentionally NOT done
here (it needs the network); supply RA/Dec in the config for headless runs.
"""

from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path
from typing import List

from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths import step1_dir


def _valid_radec(ra, dec) -> bool:
    try:
        ra = float(ra)
        dec = float(dec)
    except (TypeError, ValueError):
        return False
    if not (math.isfinite(ra) and math.isfinite(dec)):
        return False
    return 0.0 <= ra < 360.0 and -90.0 <= dec <= 90.0


class ScanStep(PipelineStep):
    index = 1
    key = "scan"
    title = "File scan & target resolution"
    interactive = False

    def outputs(self, ctx: RunContext) -> List[Path]:
        out = step1_dir(ctx.result_dir)
        return [out / "headers.csv", out / "selection.json"]

    def run(self, ctx: RunContext) -> StepResult:
        from apex.core.file_manager import FileManager
        from apex.gui.workflow.step1_header_target import select_header_target

        fm = FileManager(ctx.params)
        filenames = fm.scan_files()  # raises RuntimeError if none found
        df = fm.read_headers()       # writes step1_file_selection/headers.csv

        # Target: config first, then header inference.
        p = ctx.params.P
        ra = getattr(p, "target_ra_deg", None)
        dec = getattr(p, "target_dec_deg", None)
        if _valid_radec(ra, dec):
            target = {
                "filename": None,
                "object_name": getattr(p, "target_name", "") or "",
                "ra_deg": float(ra),
                "dec_deg": float(dec),
                "source": "config",
            }
        else:
            inferred = select_header_target(df, filenames)
            if inferred is not None:
                inferred["source"] = "header"
            target = inferred

        out_dir = step1_dir(ctx.result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        selection = {
            "generated": datetime.now().isoformat(),
            "data_dir": str(ctx.data_dir),
            "n_files": len(filenames),
            "filenames": list(filenames),
            "target": target,
        }
        sel_path = out_dir / "selection.json"
        sel_path.write_text(json.dumps(selection, indent=2, ensure_ascii=False), encoding="utf-8")

        if target is None:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.OK,
                message=(f"{len(filenames)} files scanned; no valid target in config "
                         f"or headers (set [target].ra_deg/dec_deg for downstream WCS)"),
                outputs=[str(out_dir / "headers.csv"), str(sel_path)],
            )
        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=(f"{len(filenames)} files scanned; target {target['object_name']} "
                     f"@ ({target['ra_deg']:.5f}, {target['dec_deg']:.5f}) [{target['source']}]"),
            outputs=[str(out_dir / "headers.csv"), str(sel_path)],
        )
