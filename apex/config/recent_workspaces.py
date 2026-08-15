"""The list of workspaces the user has opened, kept across restarts.

Switching target currently means editing a config file by hand or relaunching
with ``--params``. A Recent menu only helps if it survives the restart it is
meant to save, so the list lives in ``~/.apex/recent.json`` rather than in
memory or in a repo file — the repo is shared between checkouts and a user's
own history is not part of the project.

Entries are ordered most-recent-first, deduplicated by resolved path, and
pruned of files that no longer exist, so a menu built from this never offers a
dead path.
"""

from __future__ import annotations

import json
from pathlib import Path

RECENT_PATH = Path.home() / ".apex" / "recent.json"
MAX_ENTRIES = 10


def _resolve(path: str | Path) -> Path:
    """Absolute path without following junctions out of the repo.

    ``validation/`` and some workspaces are directory junctions; ``resolve()``
    would rewrite them to their target and make two names for one workspace
    look like two entries.
    """
    return Path(path).expanduser().absolute()


def load_recent(recent_path: Path | None = None) -> list[Path]:
    """Recent workspaces, newest first, with vanished ones dropped."""
    path = Path(recent_path) if recent_path is not None else RECENT_PATH
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    entries = raw.get("workspaces") if isinstance(raw, dict) else raw
    if not isinstance(entries, list):
        return []
    out: list[Path] = []
    seen: set[str] = set()
    for entry in entries:
        if not isinstance(entry, str) or not entry.strip():
            continue
        candidate = _resolve(entry)
        key = str(candidate).lower()
        if key in seen or not candidate.is_file():
            continue
        seen.add(key)
        out.append(candidate)
    return out[:MAX_ENTRIES]


def remember(config_path: str | Path, recent_path: Path | None = None) -> list[Path]:
    """Put one workspace at the top of the list and persist it."""
    path = Path(recent_path) if recent_path is not None else RECENT_PATH
    target = _resolve(config_path)
    entries = [p for p in load_recent(path)
               if str(p).lower() != str(target).lower()]
    entries.insert(0, target)
    entries = entries[:MAX_ENTRIES]
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps({"workspaces": [str(p) for p in entries]},
                       indent=2, ensure_ascii=False),
            encoding="utf-8")
    except OSError:
        # A missing home directory must not stop the user from opening a
        # workspace; the menu just will not remember it.
        pass
    return entries


def forget(config_path: str | Path, recent_path: Path | None = None) -> list[Path]:
    """Drop one workspace from the list, for when a path is gone for good."""
    path = Path(recent_path) if recent_path is not None else RECENT_PATH
    target = str(_resolve(config_path)).lower()
    entries = [p for p in load_recent(path) if str(p).lower() != target]
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps({"workspaces": [str(p) for p in entries]},
                       indent=2, ensure_ascii=False),
            encoding="utf-8")
    except OSError:
        pass
    return entries


def display_label(config_path: str | Path) -> str:
    """A menu label that distinguishes workspaces at a glance.

    Every workspace's config is called ``apex_config.json``, so the file name
    alone labels them all identically; the parent directory is the part that
    says which target this is.
    """
    path = _resolve(config_path)
    return f"{path.parent.name}  ({path.name})" if path.parent.name else path.name
