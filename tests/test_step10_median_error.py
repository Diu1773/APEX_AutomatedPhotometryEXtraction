"""A per-star magnitude error of exactly zero is a lie, not a measurement.

Step 10 reports each star's magnitude as the median across frames, with the
uncertainty taken from the scatter between them: `1.4826 * MAD / sqrt(n)`. With
one frame there is no scatter, so that is exactly 0 — and 0 says the magnitude
is known perfectly.

Found on 2026-08-15 while chasing an M5 star that moved 2.7 mag after the
Step 8 seed fix. It had not been re-measured: it dropped from five frames to
one, and the single survivor carried the whole answer at zero uncertainty.
Across all five Phase-3 clusters, every single-frame entry reported 0.000
(0.4-1.0 % of entries). The isochrone reader happens to guard with `err > 0`
and substitutes 0.03, so nothing divided by zero — but the exported column was
wrong, and the next consumer might not guard.
"""

from __future__ import annotations

import numpy as np
import pytest

from apex.gui.workflow.cmd.step10_zeropoint_calibration import (
    ZeropointCalibrationWorker as _Worker,
)

_median = _Worker._robust_median_and_err


def test_one_frame_uses_that_frame_s_own_error():
    med, err, n = _median([15.0], [0.04])
    assert (med, n) == (15.0, 1)
    assert err == pytest.approx(0.04)


def test_identical_values_across_frames_are_not_error_free():
    """MAD is zero whenever every frame agrees exactly, not just at n=1."""
    med, err, n = _median([15.0, 15.0, 15.0], [0.04, 0.05, 0.06])
    assert n == 3
    assert err > 0
    assert err == pytest.approx(0.05 / np.sqrt(3))


def test_without_a_fallback_it_is_unknown_not_perfect():
    med, err, n = _median([15.0])
    assert med == 15.0 and n == 1
    assert np.isnan(err)


def test_real_scatter_still_wins():
    """The measured spread is the better estimate whenever it exists."""
    values = [15.00, 15.10, 14.90, 15.05]
    med, err, n = _median(values, [0.001] * 4)
    assert n == 4
    mad = np.median(np.abs(np.asarray(values) - np.median(values)))
    assert err == pytest.approx(1.4826 * mad / np.sqrt(4))


def test_non_finite_measurements_are_dropped_from_both_arrays():
    """The error array must be filtered alongside the magnitudes, not zipped."""
    med, err, n = _median([15.0, np.nan, np.inf], [0.09, 0.01, 0.01])
    assert n == 1
    assert err == pytest.approx(0.09)


def test_no_finite_measurement_is_nan_all_round():
    med, err, n = _median([np.nan, np.nan], [0.02, 0.02])
    assert np.isnan(med) and np.isnan(err) and n == 0


def test_a_useless_error_column_does_not_resurrect_zero():
    """Zeros and NaNs in the fallback are not usable errors either."""
    med, err, n = _median([15.0], [0.0])
    assert np.isnan(err)
    med, err, n = _median([15.0], [np.nan])
    assert np.isnan(err)
