"""Step registry: the ordered list of headless steps for each mode.

Shared steps 1-7 all run headless: 1 (scan), 2 (crop, config-driven), 3 (sky),
4 (detect), 5 (wcs), 6 (refbuild), 7 (forcedphot).

CMD adds 8 (PSF photometry) and 10 (zero-point + wide CMD table), which take a
cluster from raw frames to a calibrated CMD table in one command. Their
calculation lives in ``apex.analysis.cmd`` and needs no Qt (2026-08-17); the
desktop windows subclass those runners to add a thread and Qt signals, so the
app and this pipeline drive the same objects.

Only Step 9 stays :class:`DeferredStep`, and not because it is unwritten: the
master-ID editor is interactive by nature and Step 10 does not need it.

Steps 11 and 12 were deferred and are not any more, for opposite reasons.
11 (CMD plot) was described as having "nothing to port" — true of the viewer,
which is an instrument, and false of the figure, which a run that measures a
cluster should leave behind (2026-08-19). 12 (isochrone MCMC) always had a
Qt-free service; what it lacked was a config surface for the settings that
decide the answer, and it now refuses to run rather than guess them.
"""

from __future__ import annotations

from typing import List, Optional, Set

from apex.pipeline.base import DeferredStep, PipelineStep
from apex.pipeline.steps.calibration import CalibrationStep
from apex.pipeline.steps.scan import ScanStep
from apex.pipeline.steps.crop import CropStep
from apex.pipeline.steps.sky_qc import SkyQCStep
from apex.pipeline.steps.detect import DetectStep
from apex.pipeline.steps.wcs import WcsStep
from apex.pipeline.steps.refbuild import RefBuildStep
from apex.pipeline.steps.forcedphot import ForcedPhotStep
from apex.pipeline.steps.psf import PsfPhotometryStep
from apex.pipeline.steps.isochrone import IsochroneStep
from apex.pipeline.steps.cmdplot import CmdPlotStep
from apex.pipeline.steps.lc_target import LcTargetStep
from apex.pipeline.steps.lc_lightcurve import LcLightCurveStep
from apex.pipeline.steps.zeropoint import ZeropointStep
from apex.utils import step_paths as sp
from apex.utils import step_paths_cmd as spc


def _sel(rd):
    return sp.step1_dir(rd) / "selection.json"


def _lc_steps() -> List[PipelineStep]:
    """LC's own steps. 8 and 9 so far — see the module docstring."""
    return [LcTargetStep(), LcLightCurveStep()]


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


def _cmd_steps() -> List[PipelineStep]:
    return [
        PsfPhotometryStep(),
        DeferredStep(
            9, "masterid", "Master ID editor",
            outputs_fn=lambda ctx: [spc.step9_selection_dir(ctx.result_dir)],
            interactive=True,
        ),
        ZeropointStep(),
        CmdPlotStep(),
        IsochroneStep(),
    ]


def get_steps(mode: str) -> List[PipelineStep]:
    """Ordered steps for a mode: shared 1-7, then the mode's own.

    CMD contributes 8-12; only 9 is deferred (interactive by nature — see the
    module docstring). LC now contributes 8 (target resolution); its 9-11 are
    not ported yet, though their calculation is already Qt-free.

    Detector calibration (index 0) is deliberately excluded here — it is an
    optional off-chain pre-stage; use :func:`get_calibration_step` for it. This
    keeps the runner's 1-based ``step.index - 1`` math from ever seeing 0."""
    if mode not in ("cmd", "lc"):
        raise ValueError(f"mode must be 'cmd' or 'lc', got {mode!r}")
    steps = _shared_steps()
    steps += _cmd_steps() if mode == "cmd" else _lc_steps()
    return steps


def get_calibration_step() -> PipelineStep:
    """The optional off-chain Step 0 (detector calibration), separate from the
    numbered 1-7 chain returned by :func:`get_steps`."""
    return CalibrationStep()


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
