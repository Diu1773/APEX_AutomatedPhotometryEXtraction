"""Small statistical wrappers with optional Bottleneck acceleration."""

from __future__ import annotations

import numpy as np

try:
    import bottleneck as _bn
except Exception:  # pragma: no cover - optional runtime dependency
    _bn = None


HAS_BOTTLENECK = _bn is not None


def nanmedian(a, axis=None):
    if HAS_BOTTLENECK:
        return _bn.nanmedian(a, axis=axis)
    return np.nanmedian(a, axis=axis)


def median(a, axis=None):
    if HAS_BOTTLENECK:
        return _bn.median(a, axis=axis)
    return np.median(a, axis=axis)


def nanmean(a, axis=None):
    if HAS_BOTTLENECK:
        return _bn.nanmean(a, axis=axis)
    return np.nanmean(a, axis=axis)


def nanstd(a, axis=None, ddof: int = 0):
    if HAS_BOTTLENECK:
        return _bn.nanstd(a, axis=axis, ddof=ddof)
    return np.nanstd(a, axis=axis, ddof=ddof)


def nansum(a, axis=None):
    if HAS_BOTTLENECK:
        return _bn.nansum(a, axis=axis)
    return np.nansum(a, axis=axis)


def nanmax(a, axis=None):
    if HAS_BOTTLENECK:
        return _bn.nanmax(a, axis=axis)
    return np.nanmax(a, axis=axis)


def finite_values(a) -> np.ndarray:
    arr = np.asarray(a, dtype=float)
    return arr[np.isfinite(arr)]


def finite_nanmedian(a, default: float = np.nan) -> float:
    vals = finite_values(a)
    if vals.size == 0:
        return float(default)
    return float(nanmedian(vals))


def finite_nanstd(a, ddof: int = 0, default: float = np.nan) -> float:
    vals = finite_values(a)
    if vals.size <= ddof:
        return float(default)
    return float(nanstd(vals, ddof=ddof))


def robust_median_mad(a, default: float = np.nan) -> tuple[float, float]:
    vals = finite_values(a)
    if vals.size == 0:
        return float(default), float(default)
    med = float(nanmedian(vals))
    mad = float(nanmedian(np.abs(vals - med)))
    return med, mad
