"""Best-effort update check against GitHub Releases.

Qt-free and exception-safe: it never raises and never blocks the GUI on its own
(callers should run :func:`check_for_update` off the main thread). Used to show a
non-intrusive "a newer APEX is available" hint on startup.
"""

from __future__ import annotations

import re
from typing import Optional, Tuple

from apex import __version__

_RELEASES_API = (
    "https://api.github.com/repos/Diu1773/"
    "APEX_AutomatedPhotometryEXtraction/releases/latest"
)
_RELEASES_PAGE = (
    "https://github.com/Diu1773/"
    "APEX_AutomatedPhotometryEXtraction/releases/latest"
)


def parse_version(value: str) -> Tuple[int, ...]:
    """Parse a version string into a comparable tuple of ints.

    Tolerant of a leading 'v' and of pre-release/build suffixes: only the
    leading dotted-numeric core is used (e.g. 'v1.2.3-rc1' -> (1, 2, 3)).
    Returns (0,) on anything unparseable.
    """
    if not value:
        return (0,)
    s = str(value).strip().lstrip("vV")
    m = re.match(r"^(\d+(?:\.\d+)*)", s)
    if not m:
        return (0,)
    return tuple(int(p) for p in m.group(1).split("."))


def is_newer(remote: str, local: str) -> bool:
    """True if ``remote`` is a strictly newer version than ``local``."""
    a, b = parse_version(remote), parse_version(local)
    n = max(len(a), len(b))
    a = a + (0,) * (n - len(a))
    b = b + (0,) * (n - len(b))
    return a > b


def check_for_update(current: str = __version__, timeout: float = 6.0) -> Optional[dict]:
    """Query GitHub for the latest release. Best-effort; never raises.

    Returns a dict ``{"available", "latest", "current", "url"}`` on success, or
    ``None`` if the check could not be completed (offline, rate-limited, no
    releases yet, etc.).
    """
    try:
        import requests
    except Exception:
        return None
    try:
        resp = requests.get(
            _RELEASES_API,
            timeout=timeout,
            headers={"Accept": "application/vnd.github+json"},
        )
        if resp.status_code != 200:
            return None
        tag = (resp.json() or {}).get("tag_name")
        if not tag:
            return None
    except Exception:
        return None
    latest = str(tag).lstrip("vV")
    return {
        "available": is_newer(latest, current),
        "latest": latest,
        "current": current,
        "url": _RELEASES_PAGE,
    }
