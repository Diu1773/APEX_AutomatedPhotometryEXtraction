"""What has ever been done in a result directory.

`pipeline_run.json` answers "what did the last run do". That turned out to be
the wrong question. It is written with `write_text`, so on 2026-08-23 a
`--steps 10` run over M13 replaced a seven-step record with a one-step one: the
directory still held the products of steps 1-7 and no longer held any statement
that they had run. Of twelve directories behind the paper's numbers, ten had no
parameter record at all, and the two that did had been written that morning by
accident.

So this file keeps the other record — the one that only grows:

    <result_dir>/apex_journal.jsonl

One JSON object per line, appended, never rewritten. A directory's journal is
the directory's history: every run that touched it, every step inside that run,
and **the parameter values that step actually read while it ran**. The existing
`parameters_used.json` records the settings as they stood at the end of the run;
that is a different claim, and it is silently wrong for any product made before
the last edit of `apex_config.json`.

Four kinds of line, distinguished by `event`:

`run_start`  the plan, the config file and its SHA-256, the environment, argv
`step`       one step's outcome, duration, outputs, and settings it read
`run_end`    success, wall time, the worker counts each stage was given
`note`       anything a human or the GUI wants on the record

Appending is the whole design. A run killed halfway still leaves every step it
finished, in order, which is exactly the case a rewritten manifest loses. The
cost is that a reader must fold the lines to get "current state" — `history()`
and `latest_steps()` do that.

Nothing here may raise into a run. A journal that breaks a twelve-hour
reprocess to report that it could not describe it is worse than no journal, so
every write is guarded and failures are logged and dropped.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Iterator, Optional

JOURNAL_NAME = "apex_journal.jsonl"

# A value that is not a number, a string or a bool gets `str()`, because the
# journal must stay readable by anything that can read JSON — including in five
# years, without this package installed.
_PLAIN = (int, float, bool, str, type(None))


def journal_path(result_dir: Path | str) -> Path:
    return Path(result_dir) / JOURNAL_NAME


_SEQUENCE = itertools.count(1)


def new_run_id() -> str:
    """Ties one run's lines together. Sortable, and unique.

    The timestamp alone is not: a script that runs two plans over one directory
    inside the same second gave both the same id, and the reader folded them
    into a single run — the exact collapse this file exists to prevent. The
    counter makes it unique within a process, the pid between processes.
    """
    return (f"{datetime.now().strftime('%Y%m%dT%H%M%S')}"
            f"-{os.getpid()}-{next(_SEQUENCE):03d}")


_SESSION_ID: Optional[str] = None


def session_id() -> str:
    """One id for the life of this process, for callers that have no "run".

    The desktop app does not run a plan; a person opens windows and finishes
    steps over an evening. Minting a fresh id at each step made every one look
    like a separate run that had been interrupted — the reader saw "no end" a
    dozen times and learned nothing. One id per session groups the sitting,
    which is the thing that actually happened.
    """
    global _SESSION_ID
    if _SESSION_ID is None:
        _SESSION_ID = new_run_id()
    return _SESSION_ID


def _plain(value: Any) -> Any:
    if isinstance(value, bool) or value is None or isinstance(value, (int, str)):
        return value
    if isinstance(value, float):
        # NaN and infinity are not JSON; they become nulls that a reader cannot
        # tell from "never set", so they are spelled out instead.
        return value if math.isfinite(value) else repr(value)
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_plain(v) for v in value]
    return str(value)


def config_fingerprint(params: Any) -> dict[str, Any]:
    """Which config file, and what was in it at the time.

    The path alone is not evidence: this session edited every workspace config
    three times. The digest is what lets a later reader say whether the file
    sitting there now is the file that was read.
    """
    path = (getattr(params, "param_file", None)
            or getattr(getattr(params, "P", None), "param_file", None))
    if not path:
        return {"path": None, "sha256": None}
    p = Path(path)
    if not p.exists():
        return {"path": str(p), "sha256": None}
    try:
        return {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
    except OSError:
        return {"path": str(p), "sha256": None}


def settings_snapshot(params: Any, names: Optional[Iterable[str]] = None) -> dict[str, Any]:
    """The values behind `names`, read off the live namespace.

    Called right after a step returns, with the names that step's recorder saw,
    so the journal holds `aperture_radius = 0.8` rather than the name alone. A
    name the namespace no longer answers to is recorded as unreadable — that is
    a real event (a step reached through `_raw`, or the namespace was swapped)
    and hiding it would make the record look complete when it is not.
    """
    P = getattr(params, "P", params)
    if names is None:
        try:
            names = [n for n in vars(P) if not n.startswith("_")]
        except TypeError:
            names = []
    out: dict[str, Any] = {}
    for name in sorted(names):
        try:
            out[name] = _plain(getattr(P, name))
        except Exception:                            # noqa: BLE001 - see docstring
            out[name] = "<unreadable>"
    return out


def append(result_dir: Path | str, event: str, payload: dict[str, Any],
           *, logger=None) -> bool:
    """Add one line. Returns whether it landed; never raises.

    **One writer per result directory.** Windows opens an append handle by
    seeking to the end and then writing, and those two steps are not atomic:
    measured with two processes appending 150 lines each, 300 attempts left
    261-281 lines on disk, every `append` having returned True. Nothing is
    corrupted — whole lines overwrite whole lines — so there is no signal at
    all. A single process is lossless (300/300).

    No lock is taken, deliberately. This module's first rule is that recording
    must never break a run, and a lock adds a way for it to: contention, a stale
    file after a kill, a new exception path in the one place that must not have
    one. The realistic collision is a desktop window and `apex run` on the same
    workspace at once, and nothing prevents that — so if you start doing it,
    this is what it costs.
    """
    try:
        target = Path(result_dir)
        target.mkdir(parents=True, exist_ok=True)
        line = {"at": datetime.now().isoformat(timespec="seconds"), "event": event}
        line.update(_plain(payload))
        with journal_path(target).open("a", encoding="utf-8", newline="\n") as handle:
            handle.write(json.dumps(line, ensure_ascii=False, default=str) + "\n")
        return True
    except Exception as exc:                         # noqa: BLE001 - see module docstring
        if logger is not None:
            logger.warning("Could not append to the journal: %s", exc)
        return False


def record_run_start(result_dir, run_id, *, mode, plan, params=None,
                     environment=None, force=False, dry_run=False,
                     source="headless", logger=None) -> bool:
    return append(result_dir, "run_start", {
        "run": run_id,
        "mode": mode,
        "source": source,
        "plan": list(plan),
        "force": bool(force),
        "dry_run": bool(dry_run),
        "config": config_fingerprint(params) if params is not None else {},
        "environment": environment or {},
        "argv": list(sys.argv),
        "cwd": str(Path.cwd()),
    }, logger=logger)


def record_step(result_dir, run_id, *, index, key, title="", status="",
                message="", duration_s=None, outputs=(), settings=None,
                settings_scope="read_by_step", source="headless",
                logger=None) -> bool:
    """One step, in whichever way it was run.

    `index` is 1-based — step 1 is the file scan — because that is what the
    pipeline, the CLI and the window titles all say. The GUI numbers its own
    windows from zero internally, so its caller has to add one; if it forgets,
    the two records disagree about which step made a file and the journal
    becomes worse than nothing.

    `settings_scope` says what the values mean, and the two callers mean
    different things. Headless stands a recorder in front of the namespace and
    records what the step *read* (`read_by_step`). The GUI has no such recorder,
    so it records the workspace's settings as they stood (`workspace`) — a
    weaker claim, and one a reader must be able to tell apart.
    """
    return append(result_dir, "step", {
        "run": run_id,
        "index": index,
        "key": key,
        "title": title,
        "status": str(status),
        "source": source,
        "message": message,
        "duration_s": round(duration_s, 3) if isinstance(duration_s, (int, float)) else None,
        "outputs": [str(p) for p in outputs],
        "settings_scope": settings_scope,
        "settings": settings or {},
    }, logger=logger)


def last_step(result_dir: Path | str, index: int) -> Optional[dict]:
    """The most recent line for one step, or None.

    Lets a caller that fires on something other than "the step just ran" — the
    GUI finalises a step whenever its window closes on a valid state — skip a
    line that would repeat the one before it. Twenty identical lines from
    opening and closing a window would bury the runs worth reading.
    """
    found = None
    for line in read(result_dir):
        if line.get("event") == "step" and line.get("index") == index:
            found = line
    return found


def record_run_end(result_dir, run_id, *, success, duration_s=None,
                   worker_decisions=None, logger=None) -> bool:
    return append(result_dir, "run_end", {
        "run": run_id,
        "success": bool(success),
        "duration_s": round(duration_s, 3) if isinstance(duration_s, (int, float)) else None,
        "worker_decisions": worker_decisions or {},
    }, logger=logger)


def record_note(result_dir, text, *, run_id=None, author="", logger=None) -> bool:
    """A line for something the pipeline cannot know — why this directory exists.

    `result_nocr`, `result_pre20260807`, `result_psf`: three directory names,
    three different meanings, none of them written down anywhere.
    """
    return append(result_dir, "note", {
        "run": run_id, "author": author or os.environ.get("USERNAME", ""),
        "text": str(text),
    }, logger=logger)


# ── reading ────────────────────────────────────────────────────────────────

def read(result_dir: Path | str) -> Iterator[dict]:
    """Every line, in order. A corrupt line is yielded as one, not skipped.

    A half-written last line is the normal shape of a killed run, and a reader
    that drops it silently would under-report exactly the runs worth looking at.
    """
    path = journal_path(result_dir)
    if not path.exists():
        return
    with path.open("r", encoding="utf-8") as handle:
        for number, raw in enumerate(handle, start=1):
            raw = raw.strip()
            if not raw:
                continue
            try:
                yield json.loads(raw)
            except json.JSONDecodeError:
                yield {"event": "unreadable", "line": number, "raw": raw[:200]}


def history(result_dir: Path | str) -> list[dict]:
    """One entry per run, oldest first, with its steps folded in."""
    runs: dict[str, dict] = {}
    order: list[str] = []
    for line in read(result_dir):
        run_id = line.get("run") or "(no run id)"
        if run_id not in runs:
            runs[run_id] = {"run": run_id, "steps": [], "notes": [],
                            "started": line.get("at"), "ended": None,
                            "mode": None, "source": None, "success": None,
                            "announced": False, "last": line.get("at"),
                            "config": {}, "environment": {}}
            order.append(run_id)
        entry = runs[run_id]
        entry["last"] = line.get("at")
        kind = line.get("event")
        if kind == "run_start":
            # Only a plan that announced itself can be said to have been cut
            # short. GUI lines never announce, and calling them interrupted
            # taught a reader the opposite of what happened.
            entry["announced"] = True
            entry.update({"started": line.get("at"), "mode": line.get("mode"),
                          "source": line.get("source"), "plan": line.get("plan"),
                          "config": line.get("config") or {},
                          "environment": line.get("environment") or {},
                          "force": line.get("force"), "dry_run": line.get("dry_run")})
        elif kind == "step":
            entry["steps"].append(line)
        elif kind == "run_end":
            entry.update({"ended": line.get("at"), "success": line.get("success"),
                          "duration_s": line.get("duration_s")})
        elif kind == "note":
            entry["notes"].append(line)
    return [runs[r] for r in order]


def latest_steps(result_dir: Path | str) -> dict[int, dict]:
    """The most recent run of each step index — what the products came from.

    This is the question `pipeline_run.json` was trying to answer and could not,
    because a later partial run erased the earlier full one. Here the earlier
    lines are still on disk, so the fold can prefer the newest per step instead
    of the newest per run.
    """
    newest: dict[int, dict] = {}
    for line in read(result_dir):
        if line.get("event") != "step":
            continue
        index = line.get("index")
        if not isinstance(index, int):
            continue
        if index not in newest or str(line.get("at")) >= str(newest[index].get("at")):
            newest[index] = line
    return dict(sorted(newest.items()))
