"""Step registry: the ordered list of headless steps for each mode.

Shared steps 1-7 all run headless: 1 (scan), 2 (crop, config-driven), 3 (sky),
4 (detect), 5 (wcs), 6 (refbuild), 7 (forcedphot).

CMD adds 8 (PSF photometry) and 10 (zero-point + wide CMD table), which take a
cluster from raw frames to a calibrated CMD table in one command. Both drive
the GUI module's worker synchronously and say so when PyQt5 is absent.

Steps 9, 11 and 12 stay :class:`DeferredStep` on purpose, for different
reasons. 9 (master-ID editor) is interactive and Step 10 does not need it. 11
(CMD plot) has no headless path yet. 12 (isochrone MCMC) *does* have a Qt-free
service, but its answer is decided by settings the config does not carry —
colours, age bounds, priors — and running it on defaults produces a confident
wrong number: the default 0.2-6 Gyr window cannot reach a globular at all, and
without an E(B-V) prior an open cluster rails at the floor (measured, see
``validation/psf_crossinstrument/REPORT_UB_DEGENERACY.md``). Giving Step 12 a
config surface is a deliberate decision, not a default to be guessed here.
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
from apex.pipeline.steps.zeropoint import ZeropointStep
from apex.utils import step_paths as sp
from apex.utils import step_paths_cmd as spc


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


def _cmd_steps() -> List[PipelineStep]:
    return [
        PsfPhotometryStep(),
        DeferredStep(
            9, "masterid", "Master ID editor",
            outputs_fn=lambda ctx: [spc.step9_selection_dir(ctx.result_dir)],
            interactive=True,
        ),
        ZeropointStep(),
        DeferredStep(
            11, "cmdplot", "CMD plot",
            inputs_fn=lambda ctx: [spc.step10_zp_dir(ctx.result_dir)],
            outputs_fn=lambda ctx: [spc.step11_cmd_dir(ctx.result_dir)],
        ),
        DeferredStep(
            12, "isochrone", "Isochrone fit",
            inputs_fn=lambda ctx: [spc.step10_zp_dir(ctx.result_dir)],
            outputs_fn=lambda ctx: [spc.step12_iso_dir(ctx.result_dir)],
        ),
    ]


def get_steps(mode: str) -> List[PipelineStep]:
    """Ordered steps for a mode: shared 1-7, then the mode's own.

    CMD contributes 8-12 (8 and 10 execute; 9, 11, 12 are deferred — see the
    module docstring for why each). LC's 8-11 are not ported yet.

    Detector calibration (index 0) is deliberately excluded here — it is an
    optional off-chain pre-stage; use :func:`get_calibration_step` for it. This
    keeps the runner's 1-based ``step.index - 1`` math from ever seeing 0."""
    if mode not in ("cmd", "lc"):
        raise ValueError(f"mode must be 'cmd' or 'lc', got {mode!r}")
    steps = _shared_steps()
    if mode == "cmd":
        steps += _cmd_steps()
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
