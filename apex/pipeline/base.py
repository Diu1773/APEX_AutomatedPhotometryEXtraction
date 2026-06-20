"""Step contract for the headless pipeline.

A :class:`PipelineStep` is a thin, Qt-free description of one workflow step:
what it needs, what it produces, whether it is already done, and how to run it.
The orchestrator (:mod:`apex.pipeline.runner`) drives a list of these.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, TYPE_CHECKING

if TYPE_CHECKING:  # avoid a runtime import cycle
    from apex.pipeline.context import RunContext


class StepStatus:
    """String status constants (kept as plain strings for trivial JSON dumps)."""

    OK = "ok"
    SKIPPED = "skipped"          # already complete, not forced
    BLOCKED = "blocked"          # required inputs missing
    FAILED = "failed"           # raised during run()
    PENDING = "pending"          # dry-run: would have executed
    NOT_IMPLEMENTED = "not_implemented"  # no headless implementation yet


@dataclass
class StepResult:
    index: int
    key: str
    status: str
    message: str = ""
    outputs: List[str] = field(default_factory=list)
    duration_s: float = 0.0

    @property
    def ok(self) -> bool:
        return self.status in (StepStatus.OK, StepStatus.SKIPPED, StepStatus.PENDING)

    def to_dict(self) -> dict:
        return {
            "index": self.index,
            "key": self.key,
            "status": self.status,
            "message": self.message,
            "outputs": list(self.outputs),
            "duration_s": round(self.duration_s, 3),
        }


class PipelineStep:
    """Base class for one headless workflow step.

    Subclasses set ``index`` (1-based), ``key``, ``title`` and override
    :meth:`run`. They may override :meth:`required_inputs`, :meth:`outputs`,
    and :meth:`is_complete` to teach the orchestrator about prerequisites and
    idempotent reuse.
    """

    index: int = 0
    key: str = ""
    title: str = ""
    interactive: bool = False  # True if the GUI step needs user input (config substitutes it)

    def required_inputs(self, ctx: "RunContext") -> List[Path]:
        """Paths from earlier steps that must exist before this step can run."""
        return []

    def outputs(self, ctx: "RunContext") -> List[Path]:
        """Paths/dirs this step is expected to produce."""
        return []

    def is_complete(self, ctx: "RunContext") -> bool:
        """Default: complete iff every declared output exists and is non-empty."""
        outs = self.outputs(ctx)
        if not outs:
            return False
        for path in outs:
            if not path.exists():
                return False
            if path.is_file() and path.stat().st_size == 0:
                return False
        return True

    def missing_inputs(self, ctx: "RunContext") -> List[Path]:
        return [p for p in self.required_inputs(ctx) if not p.exists()]

    def run(self, ctx: "RunContext") -> StepResult:  # pragma: no cover - abstract
        raise NotImplementedError


class DeferredStep(PipelineStep):
    """Placeholder for a step not yet ported to headless execution.

    It still participates in planning, prerequisite checks, and ``is_complete``
    detection (so a GUI-produced output is recognised), but running it reports
    ``not_implemented`` rather than failing silently.
    """

    def __init__(self, index: int, key: str, title: str,
                 inputs_fn=None, outputs_fn=None, interactive: bool = False):
        self.index = index
        self.key = key
        self.title = title
        self.interactive = interactive
        self._inputs_fn = inputs_fn
        self._outputs_fn = outputs_fn

    def required_inputs(self, ctx: "RunContext") -> List[Path]:
        return list(self._inputs_fn(ctx)) if self._inputs_fn else []

    def outputs(self, ctx: "RunContext") -> List[Path]:
        return list(self._outputs_fn(ctx)) if self._outputs_fn else []

    def run(self, ctx: "RunContext") -> StepResult:
        return StepResult(
            index=self.index,
            key=self.key,
            status=StepStatus.NOT_IMPLEMENTED,
            message=("headless execution not yet available for this step; "
                     "run it in the GUI, or it will be ported in a later iteration"),
        )
