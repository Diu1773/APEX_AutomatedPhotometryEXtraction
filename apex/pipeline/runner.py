"""The headless pipeline orchestrator.

Resolves a step plan, checks prerequisites, runs each step idempotently, marks
optional ProjectState, and writes a JSON run manifest. No Qt, no GUI.
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
        return {
            "mode": self.mode,
            "started": self.started,
            "ended": self.ended,
            "success": self.success,
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

        for step in plan:
            label = f"Step {step.index} [{step.key}] {step.title}"

            # In dry-run nothing executes, so upstream outputs won't exist yet;
            # skip prerequisite enforcement and just preview the plan.
            missing = step.missing_inputs(ctx)
            if missing and not ctx.dry_run:
                log.error("%s -> BLOCKED (missing inputs: %s)",
                          label, ", ".join(str(p) for p in missing))
                report.results.append(StepResult(
                    index=step.index, key=step.key, status=StepStatus.BLOCKED,
                    message="missing inputs: " + ", ".join(str(p) for p in missing),
                ))
                break  # downstream steps depend on this one

            if step.is_complete(ctx) and not ctx.force:
                log.info("%s -> skipped (already complete; use --force to rerun)", label)
                report.results.append(StepResult(
                    index=step.index, key=step.key, status=StepStatus.SKIPPED,
                    message="already complete",
                    outputs=[str(p) for p in step.outputs(ctx)],
                ))
                continue

            if ctx.dry_run:
                note = " (interactive: needs config-supplied input)" if step.interactive else ""
                log.info("%s -> would run%s", label, note)
                report.results.append(StepResult(
                    index=step.index, key=step.key, status=StepStatus.PENDING,
                    message="dry-run",
                ))
                continue

            log.info("%s -> running...", label)
            t0 = time.perf_counter()
            try:
                result = step.run(ctx)
            except Exception as exc:  # noqa: BLE001 - one bad step must not crash the run
                dt = time.perf_counter() - t0
                log.exception("%s -> FAILED", label)
                report.results.append(StepResult(
                    index=step.index, key=step.key, status=StepStatus.FAILED,
                    message=f"{type(exc).__name__}: {exc}", duration_s=dt,
                ))
                break
            result.duration_s = time.perf_counter() - t0
            log.info("%s -> %s (%.2fs) %s",
                     label, result.status, result.duration_s, result.message)
            report.results.append(result)

            if result.status == StepStatus.OK and ctx.project_state is not None:
                try:
                    ctx.project_state.mark_step_completed(step.index - 1)
                except Exception:  # noqa: BLE001 - state bookkeeping must not break the run
                    log.debug("Could not mark ProjectState for step %s", step.index)

            if result.status in (StepStatus.FAILED, StepStatus.NOT_IMPLEMENTED, StepStatus.BLOCKED):
                # Stop: downstream steps would only cascade-fail.
                break

        report.ended = datetime.now().isoformat()
        if not ctx.dry_run:
            self._write_manifest(ctx, report)
        return report

    @staticmethod
    def _write_manifest(ctx: RunContext, report: RunReport) -> Optional[Path]:
        try:
            ctx.result_dir.mkdir(parents=True, exist_ok=True)
            manifest = ctx.result_dir / "pipeline_run.json"
            manifest.write_text(
                json.dumps(report.to_dict(), indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            ctx.logger.info("Run manifest written: %s", manifest)
            return manifest
        except OSError as exc:
            ctx.logger.warning("Could not write run manifest: %s", exc)
            return None
