"""Write settings back into the user's workspace config (JSON).

Extracted from the Step 1 window so any window can persist what the user
changed — a setting that resets when the window closes is a setting the user
has to re-enter every session.

Qt-free: the GUI passes a params object, this module only touches the file.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence



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
    Returns True when the file was rewritten.
    """
    from apex.config.config_io import load_config_data, save_config_data

    try:
        data, path = load_config_data(param_file or "parameters.toml")
    except Exception:
        return False
    if not data and not path.exists():
        # Nothing to merge into — the workspace has no config yet.
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

    return save_config_data(path, data)
