"""Helpers for numeric parameters keyed by canonical filter names."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping

from apex.utils.astro_utils import normalize_filter_name


def normalize_filter_float_map(values: Mapping | None) -> dict[str, float]:
    """Return finite positive values keyed by canonical APEX filter names."""
    if not isinstance(values, Mapping):
        return {}

    normalized: dict[str, float] = {}
    for raw_filter, raw_value in values.items():
        key = normalize_filter_name(raw_filter)
        if not key:
            continue
        try:
            value = float(raw_value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(value) and value > 0:
            normalized[key] = value
    return normalized


def build_active_filter_float_map(
    configured: Mapping | None,
    active_filters: Iterable[str] | None,
    default: float,
) -> dict[str, float]:
    """Materialize a filter map for the filters present in the active data.

    Existing canonical overrides are retained. Filters without an override use
    ``default``. If no active filters can be determined, the configured map is
    returned unchanged.
    """
    normalized = normalize_filter_float_map(configured)
    active = {
        key
        for raw_filter in (active_filters or ())
        if (key := normalize_filter_name(raw_filter))
    }
    if not active:
        return normalized

    default_value = float(default)
    if not math.isfinite(default_value) or default_value <= 0:
        raise ValueError("default filter value must be finite and positive")
    return {
        key: normalized.get(key, default_value)
        for key in sorted(active)
    }
