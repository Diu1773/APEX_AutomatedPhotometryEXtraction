"""Tests for the photon-transfer gain / read-noise measurement."""

from __future__ import annotations

import numpy as np
import pytest

from apex.analysis.detector_ptc import (
    PTCPoint, build_flat_pairs, fit_ptc_points, measure_ptc, FrameInfo,
    pairs_from_frames, robust_center, robust_diff_variance,
)

TRUE_GAIN = 0.68          # e-/ADU on the stored pixel
TRUE_RN_E = 2.35          # electrons
LEVELS = (2_000, 8_000, 16_000, 26_000, 38_000, 52_000)
SHAPE = (240, 240)


def _synthetic_pair(level_adu: float, rng: np.random.Generator):
    """Two flats at the same level with independent shot and read noise.

    Signal is Poisson in electrons, read noise Gaussian; both are converted
    back to ADU so the pair reproduces ``var = S/g + RN_adu**2``.
    """
    rn_adu = TRUE_RN_E / TRUE_GAIN
    electrons = level_adu * TRUE_GAIN
    # A fixed multiplicative pattern proves the pair difference cancels PRNU.
    prnu = 1.0 + 0.03 * rng.standard_normal(SHAPE)
    out = []
    for _ in range(2):
        signal = rng.poisson(electrons * prnu).astype(np.float64) / TRUE_GAIN
        out.append(signal + rng.normal(0.0, rn_adu, SHAPE))
    return out[0], out[1]


@pytest.fixture(scope="module")
def synthetic_points():
    rng = np.random.default_rng(20260806)
    points = []
    for level in LEVELS:
        a, b = _synthetic_pair(level, rng)
        points.append(PTCPoint(
            signal_adu=0.5 * (robust_center(a) + robust_center(b)),
            variance_adu2=robust_diff_variance(a, b) / 2.0,
        ))
    return points


def test_robust_diff_variance_ignores_outliers():
    rng = np.random.default_rng(7)
    a = rng.normal(100.0, 3.0, (200, 200))
    b = rng.normal(100.0, 3.0, (200, 200))
    clean = robust_diff_variance(a, b)

    spiked = b.copy()
    spiked[::37, ::41] += 5_000.0          # cosmic rays / hot pixels
    assert robust_diff_variance(a, spiked) == pytest.approx(clean, rel=0.05)


def test_free_fit_recovers_gain(synthetic_points):
    r = fit_ptc_points(synthetic_points, binning=1)
    assert r.ok
    assert r.gain_eff == pytest.approx(TRUE_GAIN, rel=0.03)


def test_read_noise_needs_bias_frames_not_the_intercept(synthetic_points):
    """Flats bright enough to flat-field cannot measure the read noise.

    At the faintest level here the shot variance is already a few thousand
    ADU^2 while RN^2 is about 12, so the intercept is a rounding error on the
    fit and comes back meaningless.  That is why the read noise is measured
    from a bias pair and, when the flats sit high, pinned into the fit.
    """
    from_intercept = fit_ptc_points(synthetic_points, binning=1)
    assert abs(from_intercept.read_noise_eff - TRUE_RN_E) > 0.5 * TRUE_RN_E

    from_bias = fit_ptc_points(synthetic_points, binning=1,
                               read_noise_adu=TRUE_RN_E / TRUE_GAIN)
    assert from_bias.read_noise_eff == pytest.approx(TRUE_RN_E, rel=0.05)


def test_pinned_intercept_matches_free_fit_when_range_is_wide(synthetic_points):
    rn_adu = TRUE_RN_E / TRUE_GAIN
    free = fit_ptc_points(synthetic_points, binning=1)
    pinned = fit_ptc_points(synthetic_points, read_noise_adu=rn_adu,
                            binning=1, fix_intercept=True)
    # With points down to 2 ke- the intercept is already well determined, so
    # pinning it must not move the slope appreciably.
    assert pinned.gain_eff == pytest.approx(free.gain_eff, rel=0.02)


def test_pinning_rescues_a_narrow_high_signal_range(synthetic_points):
    """High-signal-only flats leave the intercept unconstrained."""
    rn_adu = TRUE_RN_E / TRUE_GAIN
    narrow = [p for p in synthetic_points if p.signal_adu > 20_000]
    assert len(narrow) >= 3

    free = fit_ptc_points(narrow, binning=1)
    pinned = fit_ptc_points(narrow, read_noise_adu=rn_adu, binning=1,
                            fix_intercept=True)
    assert abs(pinned.gain_eff - TRUE_GAIN) <= abs(free.gain_eff - TRUE_GAIN)
    assert pinned.fit_intercept == pytest.approx(rn_adu ** 2)


def test_narrow_range_is_flagged():
    """Flats clustered at one exposure level must not look precise."""
    rn_adu = TRUE_RN_E / TRUE_GAIN
    clustered = [PTCPoint(signal_adu=s, variance_adu2=s / TRUE_GAIN + rn_adu ** 2)
                 for s in (44_000, 46_500, 49_000, 51_000)]
    r = fit_ptc_points(clustered, read_noise_adu=rn_adu, binning=1,
                       fix_intercept=True)
    assert r.lever_arm < 0.35
    assert "signal range" in r.message


def test_wide_range_is_not_flagged(synthetic_points):
    r = fit_ptc_points(synthetic_points, binning=1)
    assert r.lever_arm >= 0.35
    assert r.message == "ok"


def test_binning_conversion_is_reported():
    points = [PTCPoint(signal_adu=s, variance_adu2=s / 0.68 + 25.0)
              for s in (5_000, 20_000, 45_000)]
    r = fit_ptc_points(points, read_noise_adu=5.0, binning=2, fix_intercept=True)
    # An average-binned superpixel collects n^2 photosites, so its effective
    # gain is n^2 times the per-photosite value.
    assert r.gain_pixel == pytest.approx(r.gain_eff / 4.0)
    assert r.read_noise_pixel == pytest.approx(r.read_noise_eff / 2.0)


def test_header_ratio_exposes_a_wrong_egain(synthetic_points):
    r = fit_ptc_points(synthetic_points, binning=1, header_egain=0.05)
    assert r.header_ratio == pytest.approx(r.gain_eff / 0.05)
    assert r.header_ratio > 10          # the trap this tool exists to catch


def test_too_few_points_fails_cleanly():
    r = fit_ptc_points([PTCPoint(1000.0, 1500.0)], binning=1)
    assert not r.ok
    assert "at least 3" in r.message


def test_pairs_from_frames_pairs_consecutively():
    assert pairs_from_frames(["a", "b", "c", "d", "e"]) == [("a", "b"), ("c", "d")]


def test_build_flat_pairs_groups_by_filter_and_level():
    frames = [
        FrameInfo(path="v1", filter_name="V", level_adu=21_000),
        FrameInfo(path="v2", filter_name="V", level_adu=21_100),
        FrameInfo(path="b1", filter_name="B", level_adu=21_050),   # other filter
        FrameInfo(path="v3", filter_name="V", level_adu=40_000),   # unpaired
        FrameInfo(path="v4", filter_name="V", level_adu=5_000),    # below floor
    ]
    pairs = build_flat_pairs(frames, signal_floor=20_000, tolerance=0.02)
    assert [(a.path, b.path) for a, b in pairs] == [("v1", "v2")]


def test_build_flat_pairs_honours_the_ceiling():
    frames = [
        FrameInfo(path="a", filter_name="V", level_adu=60_000),
        FrameInfo(path="b", filter_name="V", level_adu=60_100),
    ]
    assert build_flat_pairs(frames, signal_floor=0, signal_ceiling=50_000) == []


def test_measure_ptc_accepts_arrays():
    rng = np.random.default_rng(11)
    pairs = [_synthetic_pair(level, rng) for level in LEVELS]
    bias_rn_adu = TRUE_RN_E / TRUE_GAIN
    bias = [(rng.normal(500.0, bias_rn_adu, SHAPE),
             rng.normal(500.0, bias_rn_adu, SHAPE)) for _ in range(4)]
    r = measure_ptc(pairs, bias_pairs=bias, bias_level=0.0, binning=1)
    assert r.ok
    assert r.gain_eff == pytest.approx(TRUE_GAIN, rel=0.05)
