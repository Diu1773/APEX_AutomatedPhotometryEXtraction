"""Canonical filesystem layout for GUI tools.

New tool outputs should live under:

    <result_dir>/tools/<tool_id>/
        outputs/
        cache/
        logs/
        config.json

Existing tools may still read legacy paths during migration, but new writes
should prefer these helpers so outputs are discoverable and consistent.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


_SAFE_TOOL_ID_RE = re.compile(r"[^a-z0-9_]+")


def normalize_tool_id(tool_id: str) -> str:
    """Return a stable lowercase id safe for directory names."""
    value = str(tool_id or "").strip().lower().replace("-", "_").replace(" ", "_")
    value = _SAFE_TOOL_ID_RE.sub("_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value or "tool"


@dataclass(frozen=True)
class ToolWorkspace:
    """Canonical paths owned by one GUI tool."""

    tool_id: str
    root: Path
    outputs: Path
    cache: Path
    logs: Path
    config: Path

    def ensure(self) -> "ToolWorkspace":
        self.outputs.mkdir(parents=True, exist_ok=True)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.logs.mkdir(parents=True, exist_ok=True)
        return self


def tools_root(result_dir: Path | str) -> Path:
    return Path(result_dir) / "tools"


def tool_workspace(result_dir: Path | str, tool_id: str, *, create: bool = True) -> ToolWorkspace:
    """Build canonical paths for a tool under a result workspace."""
    safe_id = normalize_tool_id(tool_id)
    root = tools_root(result_dir) / safe_id
    workspace = ToolWorkspace(
        tool_id=safe_id,
        root=root,
        outputs=root / "outputs",
        cache=root / "cache",
        logs=root / "logs",
        config=root / "config.json",
    )
    return workspace.ensure() if create else workspace


def write_tool_error_log(
    result_dir: Path | str,
    tool_id: str,
    traceback_text: str,
    *,
    display_name: str | None = None,
) -> Path | None:
    """Append a launch/runtime traceback to the tool's canonical log folder."""
    try:
        workspace = tool_workspace(result_dir, tool_id, create=True)
        log_path = workspace.logs / "tool_errors.log"
        stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        name = display_name or workspace.tool_id
        with log_path.open("a", encoding="utf-8") as fh:
            fh.write(f"\n[{stamp}] {name}\n")
            fh.write(str(traceback_text).rstrip())
            fh.write("\n")
        return log_path
    except Exception:
        return None

