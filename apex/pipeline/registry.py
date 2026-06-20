"""Step registry: the ordered list of headless steps for each mode.

All shared steps 1-7 are now ported to real headless execution:
1 (scan), 2 (crop, config-driven), 3 (sky), 4 (detect), 5 (wcs), 6 (refbuild),
and 7 (forcedphot). No :class:`DeferredStep` placeholders remain for the shared
pipeline. Mode-specific Steps 8+ (CMD 8-12 / LC 8-11) will be appended as they
are ported.
"""

from __future__ import annotations

from typing import List, Optional, Set

from apex.pipeline.base import DeferredStep, PipelineStep
from apex.pipeline.steps.scan import ScanStep
from apex.pipeline.steps.crop import CropStep
from apex.pipeline.steps.sky_qc import SkyQCStep
from apex.pipeline.steps.detect import DetectStep
from apex.pipeline.steps.wcs import WcsStep
from apex.pipeline.steps.refbuild import RefBuildStep
from apex.pipeline.steps.forcedphot import ForcedPhotStep
from apex.utils import step_paths as sp


def _sel(rd):
    return sp.step1_dir(rd) / "selection.json"


def _shared_steps() -> List[PipelineStep]:
    return [
        ScanStep(),
        CropStep(),
        SkyQCStep(),
        DetectStep(),
        WcsStep(),
        RefBuildStep(),
        ForcedPhotStep(),
    ]


def get_steps(mode: str) -> List[PipelineStep]:
    """Ordered steps for a mode. Today both modes share Steps 1-7; mode-specific
    Steps 8+ (CMD 8-12 / LC 8-11) will be appended as they are ported."""
    if mode not in ("cmd", "lc"):
        raise ValueError(f"mode must be 'cmd' or 'lc', got {mode!r}")
    return _shared_steps()


def parse_step_range(spec: str) -> Set[int]:
    """Parse a step spec like '1-7', '4', '2,4,6', '3-5' into a set of indices."""
    indices: Set[int] = set()
    for part in str(spec).split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo_s, hi_s = part.split("-", 1)
            lo, hi = int(lo_s), int(hi_s)
            if lo > hi:
                lo, hi = hi, lo
            indices.update(range(lo, hi + 1))
        else:
            indices.add(int(part))
    return indices
