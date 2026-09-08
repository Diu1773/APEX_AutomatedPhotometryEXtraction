"""The headless pipeline orchestrator.

Resolves a step plan, checks prerequisites, runs each step idempotently, marks
optional ProjectState, and writes a JSON run manifest. No Qt, no GUI.

Two records come out of a run, and they are not redundant. `pipeline_run.json`
is rewritten each time and describes *this* run — useful, and the reason a
`--steps 10` run once erased the record that steps 1-7 had ever happened.
`apex_journal.jsonl` is appended to as the run proceeds and describes
*everything the directory has been through*; see `apex/utils/run_journal.py`.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.pipeline.provenance import RecordingNamespace, write_parameter_record
from apex.utils import run_journal as journal


@dataclass
class RunReport:
    mode: str
    results: List[StepResult] = field(default_factory=list)
    started: str = ""
    ended: str = ""

    @property
    def success(self) -> bool:
        return all(r.ok for r in self.results)

    def to_dict(self) -> dict:
        from apex.utils.constants import get_worker_decisions

        return {
            "mode": self.mode,
            "started": self.started,
            "ended": self.ended,
            "success": self.success,
            # "auto" is not a reproducible record. Store the count each stage
            # actually got and the inputs that produced it (configured value,
            # stage ceiling, RAM headroom) so a run can be explained later.
            "worker_decisions": get_worker_decisions(),
            "steps": [r.to_dict() for r in self.results],
        }


class PipelineRunner:
    def __init__(self, steps: List[PipelineStep]):
        self.steps = sorted(steps, key=lambda s: s.index)
        self._by_index = {s.index: s for s in self.steps}

    def resolve_plan(
        self,
        from_index: Optional[int] = None,
        to_index: Optional[int] = None,
        only: Optional[set] = None,
    ) -> List[PipelineStep]:
        if only:
            return [s for s in self.steps if s.index in only]
        lo = from_index if from_index is not None else self.steps[0].index
        hi = to_index if to_index is not None else self.steps[-1].index
        return [s for s in self.steps if lo <= s.index <= hi]

    def run(
        self,
        ctx: RunContext,
        *,
        from_index: Optional[int] = None,
        to_index: Optional[int] = None,
        only: Optional[set] = None,
    ) -> RunReport:
        plan = self.resolve_plan(from_index, to_index, only)
        report = RunReport(mode=ctx.mode, started=datetime.now().isoformat())
        log = ctx.logger

        log.info("Pipeline plan (%s): steps %s",
                 ctx.mode, ", ".join(str(s.index) for s in plan) or "(none)")

        # The journal is appended to as the run goes, not assembled at the end,
        # so a run that is killed still leaves every step it finished. A dry run
        # changes nothing in the directory and so writes nothing to its history.
        keep_journal = not ctx.dry_run
        run_id = journal.new_run_id()
        titles = {s.index: s.title for s in plan}
        if keep_journal:
            journal.record_run_start(
                ctx.result_dir, run_id, mode=ctx.mode,
                plan=[s.index for s in plan], params=ctx.params,
                environment=self._environment(), force=ctx.force,
                dry_run=ctx.dry_run, logger=log)

        def emit(result: StepResult, settings: Optional[dict] = None) -> None:
            """Every outcome goes to both records — blocked and skipped too.

            "Step 4 was skipped because it was already complete" is the kind of
            line that explains a directory a year later, and it is exactly what
            a report holding only the steps that ran would leave out.
            """
            report.results.append(result)
            if keep_journal:
                journal.record_step(
                    ctx.result_dir, run_id, index=result.index, key=result.key,
                    title=titles.get(result.index, ""), status=result.status,
                    message=result.message, duration_s=result.duration_s,
                    outputs=result.outputs or (), settings=settings or {},
                    logger=log)

        settings_read: dict = {}
        for step in plan:
            label = f"Step {step.index} [{step.key}] {step.title}"

            # In dry-run nothing executes, so upstream outputs won't exist yet;
            # skip prerequisite enforcement and just preview the plan.
            missing = step.missing_inputs(ctx)
            if missing and not ctx.dry_run:
                log.error("%s -> BLOCKED (missing inputs: %s)",
                          label, ", ".join(str(p) for p in missing))
                emit(StepResult(
                    index=step.index, key=step.key, status=StepStatus.BLOCKED,
                    message="missing inputs: " + ", ".join(str(p) for p in missing),
                ))
                break  # downstream steps depend on this one

            if step.is_complete(ctx) and not ctx.force:
                log.info("%s -> skipped (already complete; use --force to rerun)", label)
                emit(StepResult(
                    index=step.index, key=step.key, status=StepStatus.SKIPPED,
                    message="already complete",
                    outputs=[str(p) for p in step.outputs(ctx)],
                ))
                continue

            if ctx.dry_run:
                note = " (interactive: needs config-supplied input)" if step.interactive else ""
                log.info("%s -> would run%s", label, note)
                emit(StepResult(
                    index=step.index, key=step.key, status=StepStatus.PENDING,
                    message="dry-run",
                ))
                continue

            log.info("%s -> running...", label)
            t0 = time.perf_counter()
            # Stand a recorder in front of the parameters for the duration of
            # the step, so the manifest can say which settings this step read
            # rather than which settings existed. Restored in `finally` — a step
            # that raises must not leave the proxy in place for the next one.
            real_P = getattr(ctx.params, "P", None)
            recorder = None
            if real_P is not None:
                recorder = RecordingNamespace(real_P)
                try:
                    ctx.params.P = recorder
                except Exception:                   # noqa: BLE001 - frozen params
                    recorder = None
            failure = None
            try:
                result = step.run(ctx)
            except Exception as exc:  # noqa: BLE001 - one bad step must not crash the run
                log.exception("%s -> FAILED", label)
                failure = exc
            finally:
                # A `finally`, which the comment above claimed since 2026-08-18
                # and the code never had: the restore lived in the two ordinary
                # branches, so anything that is not an `Exception` walked past
                # it. `KeyboardInterrupt` is the one that matters — stopping a
                # twelve-hour reprocess is the case this module is built for —
                # and it left the proxy standing in front of the parameters for
                # whatever ran next.
                if recorder is not None:
                    try:
                        ctx.params.P = real_P
                    except Exception:       # noqa: BLE001 - frozen params
                        log.debug("Could not restore the parameters after %s", label)

            used = {}
            if recorder is not None:
                settings_read[step.key] = recorder.seen
                # The values, not just the names. `parameters_used.json` holds
                # the settings as they stood when the run ended; this holds what
                # this step read while it ran, which is the same thing only if
                # nobody edited the config in between.
                used = journal.settings_snapshot(ctx.params, recorder.seen)

            if failure is not None:
                emit(StepResult(
                    index=step.index, key=step.key, status=StepStatus.FAILED,
                    message=f"{type(failure).__name__}: {failure}",
                    duration_s=time.perf_counter() - t0,
                ), used)
                break

            result.duration_s = time.perf_counter() - t0
            log.info("%s -> %s (%.2fs) %s",
                     label, result.status, result.duration_s, result.message)
            emit(result, used)

            # `step.index - 1` is the GUI's 0-based slot for the numbered chain
            # (Step 1 -> 0). Detector calibration is off-chain with index 0, so
            # the same arithmetic would append -1 to `completed_steps` — no
            # crash, but a slot the windows cannot read, written into the
            # workspace's saved progress. Off-chain steps keep their own state.
            if (result.status == StepStatus.OK and ctx.project_state is not None
                    and step.index >= 1):
                try:
                    ctx.project_state.mark_step_completed(step.index - 1)
                except Exception:  # noqa: BLE001 - state bookkeeping must not break the run
                    log.debug("Could not mark ProjectState for step %s", step.index)

            if result.status in (StepStatus.FAILED, StepStatus.NOT_IMPLEMENTED, StepStatus.BLOCKED):
                # Stop: downstream steps would only cascade-fail.
                break

        report.ended = datetime.now().isoformat()
        if keep_journal:
            from apex.utils.constants import get_worker_decisions

            journal.record_run_end(
                ctx.result_dir, run_id, success=report.success,
                worker_decisions=get_worker_decisions(), logger=log)
        if not ctx.dry_run:
            self._write_manifest(ctx, report)
            if getattr(ctx.params, "P", None) is None:
                log.debug("No parameters to record for this run")
                return report
            try:
                written = write_parameter_record(
                    ctx.result_dir, ctx.params, ctx.mode, settings_read,
                    environment=self._environment(),
                )
                log.info("Parameters recorded: %s",
                         ", ".join(p.name for p in written))
            except Exception:  # noqa: BLE001 - a record must not fail a finished run
                log.exception("Could not write the parameter record")
        return report

    @staticmethod
    def _environment() -> dict:
        """The versions that produced this run.

        Added 2026-08-17 after a measured surprise: the same code on the same
        input gave a different answer for five of 22,305 crowded measurements
        (<=5e-5 mag) purely because one machine had a newer scipy. The run had
        recorded nothing about its own environment, so the difference took a
        control run and a package-by-package diff to explain. Now every run says
        what it was made with.
        """
        import platform
        import sys

        versions = {}
        for name in ("numpy", "scipy", "astropy", "photutils", "pandas",
                     "sep", "bottleneck", "astroquery", "emcee", "matplotlib"):
            try:
                import importlib.metadata as md

                versions[name] = md.version(name)
            except Exception:                                    # noqa: BLE001
                versions[name] = None
        try:
            import apex

            apex_version = getattr(apex, "__version__", None)
        except Exception:                                        # noqa: BLE001
            apex_version = None
        return {
            "apex": apex_version,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "packages": versions,
        }

    @staticmethod
    def _write_manifest(ctx: RunContext, report: RunReport) -> Optional[Path]:
        try:
            ctx.result_dir.mkdir(parents=True, exist_ok=True)
            manifest = ctx.result_dir / "pipeline_run.json"
            payload = report.to_dict()
            payload["environment"] = PipelineRunner._environment()
            manifest.write_text(
                json.dumps(payload, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            ctx.logger.info("Run manifest written: %s", manifest)
            return manifest
        except OSError as exc:
            ctx.logger.warning("Could not write run manifest: %s", exc)
            return None
