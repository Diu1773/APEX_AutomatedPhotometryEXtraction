"""Small statistical wrappers with optional Bottleneck acceleration.

Bottleneck accumulates in the input dtype with a naive (non-pairwise) sum.
On a float32 full-frame reduction that is catastrophic: ``bn.nanstd`` came out
21 % high on a 3194x4788 sky-level frame (measured 2026-08-07) because the
running sum of squares exhausts float32 precision; NumPy escapes only because
it sums pairwise.  The accumulating wrappers (``nanmean``/``nanstd``/``nansum``)
therefore upcast sub-double floats to float64 for full (``axis=None``)
reductions.  Axis-wise reductions keep the input dtype — the pipeline only
reduces short axes there (frame stacks, N of order tens), where the
accumulation error is far below the noise, and keeping the dtype preserves
bit-identity with the 2026-08-07 baseline run.  Order statistics
(``median``/``nanmedian``/``nanmax``) are dtype-safe and untouched.
"""

from __future__ import annotations

import numpy as np

try:
    import bottleneck as _bn
except Exception:  # pragma: no cover - optional runtime dependency
    _bn = None


HAS_BOTTLENECK = _bn is not None


def _accum_safe(a, axis):
    """Upcast small floats for full reductions (see module docstring)."""
    arr = np.asarray(a)
    if axis is None and arr.dtype.kind == "f" and arr.dtype.itemsize < 8:
        return arr.astype(np.float64)
    return arr


def nanmedian(a, axis=None):
    if HAS_BOTTLENECK:
        return _bn.nanmedian(a, axis=axis)
    return np.nanmedian(a, axis=axis)


def median(a, axis=None):
    if HAS_BOTTLENECK:
        return _bn.median(a, axis=axis)
    return np.median(a, axis=axis)


def nanmean(a, axis=None):
    a = _accum_safe(a, axis)
    if HAS_BOTTLENECK:
        return _bn.nanmean(a, axis=axis)
    return np.nanmean(a, axis=axis)


def nanstd(a, axis=None, ddof: int = 0):
    a = _accum_safe(a, axis)
    if HAS_BOTTLENECK:
        return _bn.nanstd(a, axis=axis, ddof=ddof)
    return np.nanstd(a, axis=axis, ddof=ddof)


def nansum(a, axis=None):
    a = _accum_safe(a, axis)
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
