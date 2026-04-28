"""
Project state management for workflow tracking and persistence.
"""

from __future__ import annotations
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional, List


class ProjectState:
    """
    Manages workflow state for step-by-step processing.
    Tracks completion status and saves/loads progress.

    The step list is injected at construction time so this class works
    unchanged for both cmd (13 steps) and lightcurve (12 steps) modes.
    """

    def __init__(self, project_dir: Path, steps: Optional[List[str]] = None):
        """
        Args:
            project_dir: Directory containing project data and state file.
            steps: Ordered list of step key names for this pipeline mode.
                   If None, an empty list is used (steps can be set later
                   via assign_steps()).
        """
        self.project_dir = Path(project_dir)
        self.state_file = self.project_dir / "project_state.json"
        self.steps: List[str] = list(steps) if steps is not None else []

        self.state: Dict[str, Any] = {
            "project_name": "Untitled Project",
            "created": datetime.now().isoformat(),
            "last_modified": datetime.now().isoformat(),
            "current_step": 0,
            "completed_steps": [],
            "step_data": {},
        }

        if self.state_file.exists():
            self.load()

    def assign_steps(self, steps: List[str]) -> None:
        """Set the step list after construction and re-normalize state."""
        self.steps = list(steps)
        self._normalize_state()

    # ── Step access helpers ───────────────────────────────────────────────────

    def is_step_completed(self, step_index: int) -> bool:
        return step_index in self.state["completed_steps"]

    def is_step_accessible(self, step_index: int) -> bool:
        if step_index == 0:
            return True
        return self.is_step_completed(step_index - 1)

    def mark_step_completed(self, step_index: int):
        if step_index not in self.state["completed_steps"]:
            self.state["completed_steps"].append(step_index)
            self.state["completed_steps"].sort()
        self.state["last_modified"] = datetime.now().isoformat()
        self.save()

    def mark_step_incomplete(self, step_index: int):
        if step_index in self.state["completed_steps"]:
            self.state["completed_steps"].remove(step_index)
        self.state["last_modified"] = datetime.now().isoformat()
        self.save()

    def set_current_step(self, step_index: int):
        if 0 <= step_index < len(self.steps):
            self.state["current_step"] = step_index
            self.state["last_modified"] = datetime.now().isoformat()
            self.save()

    def get_current_step(self) -> int:
        return self.state["current_step"]

    def get_next_incomplete_step(self) -> Optional[int]:
        for i in range(len(self.steps)):
            if not self.is_step_completed(i):
                return i
        return None

    def can_proceed_to_next(self) -> bool:
        current = self.state["current_step"]
        return self.is_step_completed(current) and current < len(self.steps) - 1

    def can_go_to_previous(self) -> bool:
        return self.state["current_step"] > 0

    # ── Step data storage ─────────────────────────────────────────────────────

    def store_step_data(self, step_name: str, data: Dict[str, Any]):
        if step_name not in self.state["step_data"]:
            self.state["step_data"][step_name] = {}
        self.state["step_data"][step_name].update(data)
        self.state["last_modified"] = datetime.now().isoformat()
        self.save()

    def get_step_data(self, step_name: str) -> Dict[str, Any]:
        return self.state["step_data"].get(step_name, {})

    # ── Persistence ───────────────────────────────────────────────────────────

    def save(self):
        self.project_dir.mkdir(parents=True, exist_ok=True)
        with open(self.state_file, "w", encoding="utf-8") as f:
            json.dump(self.state, f, indent=2, ensure_ascii=False)

    def load(self):
        try:
            with open(self.state_file, "r", encoding="utf-8") as f:
                self.state.update(json.load(f))
            self._normalize_state()
        except Exception as e:
            print(f"Warning: Could not load project state: {e}")

    def _normalize_state(self) -> None:
        max_index = len(self.steps) - 1
        completed = [
            i for i in self.state.get("completed_steps", [])
            if isinstance(i, int) and 0 <= i <= max_index
        ]
        self.state["completed_steps"] = sorted(set(completed))

        current = self.state.get("current_step", 0)
        if not isinstance(current, int) or current < 0:
            current = 0
        if current > max_index:
            current = max_index
        self.state["current_step"] = current

    def reset(self):
        project_name = self.state["project_name"]
        created = self.state["created"]
        self.state = {
            "project_name": project_name,
            "created": created,
            "last_modified": datetime.now().isoformat(),
            "current_step": 0,
            "completed_steps": [],
            "step_data": {},
        }
        self.save()

    def export_summary(self) -> str:
        lines = [
            f"Project: {self.state['project_name']}",
            f"Created: {self.state['created']}",
            f"Last Modified: {self.state['last_modified']}",
            f"\nProgress: {len(self.state['completed_steps'])}/{len(self.steps)} steps finished",
            "\nSteps:",
        ]
        for i, step_name in enumerate(self.steps):
            status = "✓" if i in self.state["completed_steps"] else "○"
            current = " (current)" if i == self.state["current_step"] else ""
            lines.append(f"  {status} {i+1}. {step_name.replace('_', ' ').title()}{current}")
        return "\n".join(lines)
