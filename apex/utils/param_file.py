"""Write settings back into the user's ``parameters.toml``.

Extracted from the Step 1 window so any window can persist what the user
changed — a setting that resets when the window closes is a setting the user
has to re-enter every session.

Qt-free: the GUI passes a params object, this module only touches the file.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence

try:                                          # optional writer dependency
    import tomli_w
except ImportError:                           # pragma: no cover - env-dependent
    tomli_w = None


def _section(data: Dict[str, Any], path: Sequence[str]) -> Dict[str, Any]:
    """Fetch (creating as needed) a nested table inside the parsed TOML."""
    node = data
    for key in path:
        child = node.get(key)
        if not isinstance(child, dict):
            child = {}
            node[key] = child
        node = child
    return node


def update_param_file(param_file, sections: Mapping[Sequence[str], Mapping[str, Any]],
                      params: Optional[Any] = None) -> bool:
    """Merge ``sections`` into ``param_file`` in place.

    ``sections`` maps a TOML path (``("calibration",)``, ``("calibration",
    "overscan")``) to the keys to set there; a ``None`` value removes the key.
    Returns True when the file was rewritten.  Falls back to
    ``params.save_toml()`` when the TOML writer is unavailable.
    """
    if tomli_w is None:
        if params is not None and hasattr(params, "save_toml"):
            try:
                params.save_toml()
            except Exception:
                return False
        return False

    path = Path(param_file or "parameters.toml")
    if not path.exists():
        return False

    try:
        from apex.utils.io_utils import load_toml
        data = load_toml(path)                # BOM-tolerant
    except Exception:
        return False

    for keys, updates in sections.items():
        if not updates:
            continue
        block = _section(data, tuple(keys))
        for key, value in updates.items():
            if value is None:
                block.pop(key, None)
            else:
                block[key] = value

    try:
        with path.open("wb") as fh:
            tomli_w.dump(data, fh)
    except Exception:
        return False
    return True
