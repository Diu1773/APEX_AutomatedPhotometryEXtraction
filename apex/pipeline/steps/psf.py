"""CMD Step 8 (headless): PSF photometry.

The compute lives in ``apex.analysis.cmd.psf_photometry_runner`` and needs no
Qt at all (2026-08-16). It used to sit inside the GUI module as a ``QThread``
subclass, so this step had to check for PyQt5 — an *optional* dependency — and
report NOT_IMPLEMENTED without it, which meant a script could not do PSF
photometry. The GUI now subclasses the same runner to add the thread and its
signals, so the window and this step drive the same object rather than two
copies that have to be kept in agreement.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import List

from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths import step7_forced_phot_dir
from apex.utils.step_paths_cmd import step8_psf_dir

SIGNATURE_NAME = "psf_output_signature.json"


def _frames_from_step7(result_dir: Path) -> List[str]:
    """Frame names Step 7 actually photometered, in a stable order.

    Step 8 measures what Step 7 produced, not what Step 1 selected — a frame
    that failed earlier has no forced catalogue to fit against.
    """
    prefix, suffix = "photometry_", ".tsv"
    return sorted(
        path.name[len(prefix):-len(suffix)]
        for path in step7_forced_phot_dir(result_dir).glob(f"{prefix}*{suffix}")
    )


class PsfPhotometryStep(PipelineStep):
    index = 8
    key = "psf"
    title = "PSF photometry (CMD)"
    interactive = False

    def required_inputs(self, ctx: RunContext) -> List[Path]:
        return [step7_forced_phot_dir(ctx.result_dir)]

    def outputs(self, ctx: RunContext) -> List[Path]:
        return [step8_psf_dir(ctx.result_dir)]

    def is_complete(self, ctx: RunContext) -> bool:
        """The signature, not the directory.

        Step 8 writes per-frame tables as it goes, so a half-finished run
        leaves a populated directory. The signature is written only after
        every frame succeeded — which is exactly the question being asked.
        """
        return (step8_psf_dir(ctx.result_dir) / SIGNATURE_NAME).exists()

    def run(self, ctx: RunContext) -> StepResult:
        frames = _frames_from_step7(ctx.result_dir)
        if not frames:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=(f"no Step 7 tables under "
                         f"{step7_forced_phot_dir(ctx.result_dir)}"),
            )

        from apex.analysis.cmd.psf_photometry_runner import (
            PsfPhotometryRunner,
            build_psf_output_signature,
            export_psf_qc_products,
            write_psf_output_signature,
        )

        params = ctx.params
        cache_dir = Path(params.P.cache_dir)
        if not cache_dir.is_absolute():
            cache_dir = ctx.result_dir / params.P.cache_dir

        worker = PsfPhotometryRunner(
            file_list=frames,
            params=params,
            data_dir=ctx.data_dir,
            result_dir=ctx.result_dir,
            cache_dir=cache_dir,
            use_cropped=False,
        )
        if ctx.logger is not None:
            worker._log = lambda message: ctx.logger.info("%s", message)
            worker.error.connect(lambda message: ctx.logger.error("%s", message))

        done: dict = {}
        worker.finished.connect(
            lambda payload: done.update(payload if isinstance(payload, dict) else {}))

        started = time.perf_counter()
        worker.run()
        elapsed = time.perf_counter() - started

        processed = int(done.get("processed", 0))
        stopped = int(done.get("stopped", 0))
        out_dir = step8_psf_dir(ctx.result_dir)
        signature_path = out_dir / SIGNATURE_NAME

        if processed == len(frames) and stopped == 0:
            signature_path = write_psf_output_signature(
                ctx.result_dir,
                build_psf_output_signature(params, frames, use_cropped=False,
                                           cache_dir=cache_dir),
            )
        else:
            # Never leave a signature that claims a partial run is complete —
            # Step 10 switches to PSF the moment it sees a valid one.
            signature_path.unlink(missing_ok=True)
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.FAILED,
                message=(f"{processed}/{len(frames)} frames measured "
                         f"(stopped={stopped}); signature withheld"),
                outputs=[str(out_dir)], duration_s=elapsed,
            )

        message = f"{processed}/{len(frames)} frames measured"
        try:
            qc = export_psf_qc_products(out_dir, params=params,
                                        result_dir=ctx.result_dir)
            if qc:
                message += f"; {len(qc)} QC products"
        except Exception as exc:  # noqa: BLE001 - QC must not fail the step
            if ctx.logger is not None:
                ctx.logger.warning("PSF QC export skipped: %s", exc)

        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=message, outputs=[str(out_dir), str(signature_path)],
            duration_s=elapsed,
        )
