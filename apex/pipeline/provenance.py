"""What each step was actually run with.

Until 2026-08-18 a run recorded which steps ran, how long they took, and — since
that morning — the package versions. It did not record a single parameter
value. So a table in `cmd_zeropoint/` could not say what produced it: you had to
trust that the `apex_config.json` sitting next to it had not been edited since,
and this session edited all fifty of them three times.

That is the wrong footing for a paper. A number in a manuscript should be able
to name the settings that made it.

Two records are written, because they answer different questions:

`parameters_used.json`  every resolved setting, with the config file's path and
                        SHA-256. Answers "what was the state of the world".
`parameters_used.csv`   the same as a flat table, plus which step read which
                        setting. Answers "what does the methods section say" —
                        it is meant to be pasted into an appendix.

The per-step column is measured, not declared. A recording proxy stands in for
the parameter namespace while a step runs and notes every attribute the step
asks for, so the list cannot drift from the code the way a hand-maintained one
would. A setting read through `P._raw` is not seen; those are noted as such
rather than silently omitted.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any


class RecordingNamespace:
    """Delegates to the real namespace and remembers what was asked for.

    Deliberately not a subclass: `SimpleNamespace` attribute access is what the
    whole codebase uses, so wrapping it is enough, and staying out of the type
    hierarchy means an `isinstance` check somewhere cannot change behaviour.
    """

    def __init__(self, target: Any) -> None:
        object.__setattr__(self, "_target", target)
        object.__setattr__(self, "_seen", set())

    def __getattr__(self, name: str) -> Any:
        target = object.__getattribute__(self, "_target")
        if not name.startswith("_"):
            object.__getattribute__(self, "_seen").add(name)
        return getattr(target, name)

    def __setattr__(self, name: str, value: Any) -> None:
        setattr(object.__getattribute__(self, "_target"), name, value)

    def __dir__(self):
        return dir(object.__getattribute__(self, "_target"))

    @property
    def seen(self) -> set[str]:
        return set(object.__getattribute__(self, "_seen"))


def _settings_of(params: Any) -> dict[str, Any]:
    P = getattr(params, "P", params)
    out = {}
    for name, value in vars(P).items():
        if name.startswith("_"):
            continue
        out[name] = value if isinstance(value, (int, float, bool, str, type(None))) else str(value)
    return dict(sorted(out.items()))


def _config_fingerprint(params: Any) -> dict[str, Any]:
    path = getattr(params, "param_file", None) or getattr(getattr(params, "P", None), "param_file", None)
    if not path:
        return {"path": None, "sha256": None}
    p = Path(path)
    if not p.exists():
        return {"path": str(p), "sha256": None}
    return {"path": str(p),
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


def write_parameter_record(result_dir: Path, params: Any, mode: str,
                           per_step: dict[str, set[str]],
                           environment: dict | None = None) -> list[Path]:
    """Write both records. Returns the paths written."""
    result_dir = Path(result_dir)
    result_dir.mkdir(parents=True, exist_ok=True)
    settings = _settings_of(params)

    read_by: dict[str, list[str]] = {}
    for step_key, names in sorted(per_step.items()):
        for name in names:
            read_by.setdefault(name, []).append(step_key)

    payload = {
        "mode": mode,
        "config": _config_fingerprint(params),
        "environment": environment or {},
        "settings": settings,
        "read_by_step": {k: sorted(v) for k, v in sorted(read_by.items())},
        "steps_recorded": sorted(per_step),
    }
    json_path = result_dir / "parameters_used.json"
    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False,
                                    default=str), encoding="utf-8")

    csv_path = result_dir / "parameters_used.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["setting", "value", "read_by_steps"])
        for name, value in settings.items():
            writer.writerow([name, value, " ".join(read_by.get(name, []))])
    return [json_path, csv_path]
