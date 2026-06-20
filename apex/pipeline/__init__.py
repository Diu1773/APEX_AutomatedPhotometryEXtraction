"""APEX headless pipeline orchestration.

A Qt-free layer that runs the shared workflow steps without the PyQt5 GUI.
This is the backbone of ``apex run`` and of any scripted / batch / CI use.

Layering: this package may import ``apex.core``, ``apex.config``,
``apex.utils``, and ``apex.analysis`` only. It must never import ``apex.gui``
workers. Pure, Qt-free helpers that happen to live under ``apex.gui`` (e.g.
``select_header_target``) are tolerated but should migrate to ``apex.utils``.
"""

from apex.pipeline.base import PipelineStep, StepResult, StepStatus, DeferredStep
from apex.pipeline.context import RunContext
from apex.pipeline.runner import PipelineRunner, RunReport
from apex.pipeline.registry import get_steps, parse_step_range

__all__ = [
    "PipelineStep",
    "StepResult",
    "StepStatus",
    "DeferredStep",
    "RunContext",
    "PipelineRunner",
    "RunReport",
    "get_steps",
    "parse_step_range",
]
