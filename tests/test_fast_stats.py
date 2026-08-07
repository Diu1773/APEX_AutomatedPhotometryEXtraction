"""Numerical safety of the Bottleneck-accelerated statistics wrappers.

Regression for the 2026-08-07 finding: ``bottleneck.nanstd`` on a float32
frame-sized array came out 21 % high (36.40 vs a true 30.00) because it
accumulates in the input dtype with a naive sum.  The wrappers must upcast
sub-double floats for full reductions so this can never reach a sky-noise or
QC estimate.
"""

from __future__ import annotations

import numpy as np
import pytest

from apex.utils import fast_stats


@pytest.fixture(scope="module")
def frame32():
    """A sky-level float32 frame large enough to break naive accumulation."""
    rng = np.random.default_rng(7)
    frame = rng.normal(1000.0, 30.0, (1200, 1600)).astype(np.float32)
    frame[::30, ::77] = np.nan
    return frame


def _ref(fn, frame, **kw):
    return fn(frame.astype(np.float64), **kw)


def test_full_nanstd_matches_float64_reference(frame32):
    got = float(fast_stats.nanstd(frame32))
    ref = float(_ref(np.nanstd, frame32))
    assert got == pytest.approx(ref, rel=1e-6)


def test_full_nanmean_matches_float64_reference(frame32):
    got = float(fast_stats.nanmean(frame32))
    ref = float(_ref(np.nanmean, frame32))
    assert got == pytest.approx(ref, rel=1e-9)


def test_full_nansum_matches_float64_reference(frame32):
    got = float(fast_stats.nansum(frame32))
    ref = float(_ref(np.nansum, frame32))
    assert got == pytest.approx(ref, rel=1e-9)


def test_axiswise_reduction_keeps_input_dtype(frame32):
    """Axis-wise stacks stay in their dtype — bit-compat with the baseline run.

    The pipeline only reduces short axes (frame stacks, tens of frames), where
    naive float32 accumulation is harmless; changing the dtype there would
    break bit-identity with the preserved 2026-08-07 baseline products.
    """
    stack = np.stack([frame32[:100], frame32[:100]], axis=0)
    out = fast_stats.nanmean(stack, axis=0)
    assert out.dtype == np.float32
    # and the short-axis result is still accurate
    ref = np.nanmean(stack.astype(np.float64), axis=0)
    assert np.allclose(out, ref, rtol=1e-5, equal_nan=True)


def test_nanmedian_is_dtype_safe(frame32):
    """Order statistics don't accumulate — must agree in any dtype."""
    got = float(fast_stats.nanmedian(frame32))
    ref = float(_ref(np.nanmedian, frame32))
    assert got == pytest.approx(ref, rel=1e-6)


def test_finite_nanstd_upcasts_via_finite_values(frame32):
    """The finite_* helpers were already safe (they upcast to float64)."""
    got = fast_stats.finite_nanstd(frame32, ddof=1)
    ref = float(np.nanstd(frame32.astype(np.float64), ddof=1))
    assert got == pytest.approx(ref, rel=1e-6)


def test_integer_and_float64_inputs_untouched():
    ints = np.arange(1000, dtype=np.int32)
    assert float(fast_stats.nansum(ints)) == float(ints.sum())
    doubles = np.linspace(0.0, 1.0, 1001)
    assert float(fast_stats.nanmean(doubles)) == pytest.approx(0.5)
