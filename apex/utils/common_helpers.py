"""
Shared helper functions used across multiple step modules.

Consolidates duplicated utilities (_safe_float, _normalize_filter_key, _parse_jd)
that were previously defined independently in 6-8 step files.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.time import Time

from apex.utils.astro_utils import normalize_filter_name


def safe_float(value, default: float = np.nan) -> float:
    """Safely convert value to float, returning *default* on failure."""
    try:
        return float(value)
    except Exception:
        return default


def normalize_filter_key(value: str | None) -> str:
    """Normalize a filter name to its canonical APEX key.

    Delegates to :func:`apex.utils.astro_utils.normalize_filter_name`.
    Kept as an alias for the many call-sites that import this name.
    """
    return normalize_filter_name(value)


_JOHNSON_FILTERS = {"U", "B", "V", "R", "I"}
_SDSS_FILTERS = {"u", "g", "r", "i", "z"}


def photometric_system_label(*filters: str | None) -> str:
    """Infer a photometric-system label from filter names.

    Returns 'Johnson', 'SDSS', or 'Standard' (mixed/unknown). Used for plot
    titles where the standardized magnitude system varies with the data.
    """
    names = [str(f).strip() for f in filters if f]
    if not names:
        return "Standard"
    if all(n in _JOHNSON_FILTERS for n in names):
        return "Johnson"
    if all(n in _SDSS_FILTERS for n in names):
        return "SDSS"
    return "Standard"


def target_display_name(params_or_p=None, result_dir=None, default: str = "Target") -> str:
    """Return a compact target label for plot titles."""
    obj = getattr(params_or_p, "P", params_or_p)
    for attr in ("target_name", "object_name"):
        value = getattr(obj, attr, None)
        if value is not None and str(value).strip():
            return str(value).strip()

    raw = getattr(obj, "_raw", None)
    if isinstance(raw, dict):
        for key in ("target_name", "object_name"):
            value = raw.get(key)
            if value is not None and str(value).strip():
                return str(value).strip()
        target = raw.get("target")
        if isinstance(target, dict):
            value = target.get("name")
            if value is not None and str(value).strip():
                return str(value).strip()

    path_value = result_dir or getattr(obj, "result_dir", None)
    if path_value:
        try:
            skip = {"result", "cache", "light", "pp", "dark", "flat", "bias"}
            for part in reversed(Path(path_value).parts):
                label = str(part).strip()
                if label and label.lower() not in skip:
                    return label
        except Exception:
            pass
    return default


def format_cmd_title(
    params_or_p,
    mag_label: str,
    color_label: str,
    *,
    system_label: str | None = None,
    count_text: str | None = None,
    result_dir=None,
) -> str:
    """Build consistent CMD titles, e.g. 'M67 - Standard CMD B vs B-V (N=120)'."""
    target = target_display_name(params_or_p, result_dir=result_dir)
    system = f"{str(system_label).strip()} " if system_label else ""
    title = f"{target} - {system}CMD {mag_label} vs {color_label}"
    if count_text:
        title += f" ({count_text})"
    return title


def parse_jd(date_obs: str | None) -> float:
    """Parse a DATE-OBS string into Julian Date (float).

    Returns NaN on failure or empty input.
    """
    if not date_obs:
        return np.nan
    try:
        return float(Time(str(date_obs).strip()).jd)
    except Exception:
        return np.nan
