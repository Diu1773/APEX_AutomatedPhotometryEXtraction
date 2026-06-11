"""Tests for SYSREM systematic-noise removal (Tamuz+2005)."""
from __future__ import annotations

import numpy as np
import pytest

from apex.analysis.light_curve.sysrem import (
    SysremComponent,
    SysremResult,
    apply_sysrem,
    sysrem,
)


def _make_synthetic(
    n_stars: int = 20,
    n_frames: int = 200,
    n_systematics: int = 2,
    noise: float = 0.003,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create a synthetic comparison star matrix with known systematics.

    Returns (mag_matrix, err_matrix, systematic_matrix) where
    systematic_matrix contains the injected systematics only.
    """
    rng = np.random.default_rng(seed)
    # Per-star systematic amplitudes
    a = rng.normal(0, 0.05, (n_stars, n_systematics))
    # Per-frame systematic trends
    c = rng.normal(0, 1.0, (n_systematics, n_frames))
    systematic = a @ c                      # (n_stars, n_frames)
    photon_noise = rng.normal(0, noise, (n_stars, n_frames))
    mag = systematic + photon_noise
    err = np.full_like(mag, noise)
    return mag, err, systematic


# ---------------------------------------------------------------------------
# Basic correctness
# ---------------------------------------------------------------------------

def test_sysrem_reduces_rms():
    """SYSREM residuals should have lower RMS than the raw input."""
    mag, err, _ = _make_synthetic(n_stars=20, n_frames=150, n_systematics=2, noise=0.003)
    result = sysrem(mag, err, n_iter=3)

    rms_raw = float(np.sqrt(np.nanmean(mag**2)))
    rms_final = float(np.sqrt(np.nanmean(result.residuals**2)))
    assert rms_final < rms_raw, f"RMS did not decrease: {rms_raw:.5f} → {rms_final:.5f}"


def test_sysrem_component_rms_decreasing():
    """Each successive component should lower RMS."""
    mag, err, _ = _make_synthetic()
    result = sysrem(mag, err, n_iter=4)
    rms_vals = [c.rms_after for c in result.components]
    for i in range(len(rms_vals) - 1):
        assert rms_vals[i] >= rms_vals[i + 1] - 1e-9, (
            f"RMS increased at iter {i+1}: {rms_vals[i]:.6f} → {rms_vals[i+1]:.6f}"
        )


def test_sysrem_returns_correct_shapes():
    n_s, n_f = 15, 100
    mag, err, _ = _make_synthetic(n_stars=n_s, n_frames=n_f, n_systematics=1)
    result = sysrem(mag, err, n_iter=2)

    assert result.residuals.shape == (n_s, n_f)
    for comp in result.components:
        assert comp.c.shape == (n_f,)
        assert comp.a.shape == (n_s,)


def test_sysrem_recovers_single_systematic():
    """With a single strong systematic, SYSREM iteration 1 should explain most variance."""
    rng = np.random.default_rng(0)
    n_s, n_f = 30, 300
    noise = 0.001
    c_true = np.sin(np.linspace(0, 4 * np.pi, n_f))
    a_true = rng.uniform(0.05, 0.15, n_s)
    mag = np.outer(a_true, c_true) + rng.normal(0, noise, (n_s, n_f))
    err = np.full_like(mag, noise)

    result = sysrem(mag, err, n_iter=2)
    comp = result.components[0]

    # Extracted c should correlate strongly with c_true
    corr = float(np.corrcoef(comp.c, c_true)[0, 1])
    assert abs(corr) > 0.95, f"Component correlation with injected systematic: {corr:.3f}"


# ---------------------------------------------------------------------------
# apply_sysrem
# ---------------------------------------------------------------------------

def test_apply_sysrem_recovers_target():
    """Target with a known systematic should be cleaned after apply_sysrem."""
    rng = np.random.default_rng(1)
    n_s, n_f = 25, 200
    noise = 0.002
    c_true = rng.normal(0, 1.0, n_f)

    # Comp stars: each has systematic + noise
    a_comp = rng.uniform(0.03, 0.12, n_s)
    mag = np.outer(a_comp, c_true) + rng.normal(0, noise, (n_s, n_f))
    err = np.full_like(mag, noise)

    # Target: has sinusoidal variability + the SAME systematic
    a_target = 0.08
    variability = 0.05 * np.sin(2 * np.pi * np.arange(n_f) / 30.0)
    target_raw = a_target * c_true + variability + rng.normal(0, noise, n_f)

    result = sysrem(mag, err, n_iter=2)
    corrected = apply_sysrem(target_raw, np.full(n_f, noise), result)

    # Systematic contribution to target should shrink
    systematic_only = a_target * c_true
    rms_before = float(np.std(systematic_only))
    residual_systematic = corrected - variability - np.mean(corrected - variability)
    rms_after = float(np.std(residual_systematic))
    assert rms_after < rms_before, (
        f"Systematic not reduced: {rms_before:.5f} → {rms_after:.5f}"
    )


def test_apply_sysrem_nan_handling():
    """NaN frames in target should pass through as NaN."""
    rng = np.random.default_rng(2)
    n_s, n_f = 15, 80
    mag = rng.normal(0, 0.01, (n_s, n_f))
    err = np.full_like(mag, 0.01)
    result = sysrem(mag, err, n_iter=1)

    target = rng.normal(0, 0.01, n_f)
    target[10] = np.nan
    target[50] = np.nan
    corrected = apply_sysrem(target, None, result)

    assert np.isnan(corrected[10])
    assert np.isnan(corrected[50])
    assert np.isfinite(corrected).sum() == n_f - 2


def test_apply_sysrem_n_components():
    """Applying fewer components should give intermediate correction."""
    mag, err, _ = _make_synthetic(n_stars=20, n_frames=150, n_systematics=3, noise=0.002)
    result = sysrem(mag, err, n_iter=3)
    target = mag[0].copy()
    target_err = err[0].copy()

    corr_1 = apply_sysrem(target, target_err, result, n_components=1)
    corr_3 = apply_sysrem(target, target_err, result, n_components=3)
    std_raw = float(np.std(target))
    std_1   = float(np.std(corr_1[np.isfinite(corr_1)]))
    std_3   = float(np.std(corr_3[np.isfinite(corr_3)]))

    assert std_1 <= std_raw + 1e-9
    assert std_3 <= std_1 + 1e-9


# ---------------------------------------------------------------------------
# Missing data
# ---------------------------------------------------------------------------

def test_sysrem_with_missing_data():
    """30% missing observations should not crash or produce all-NaN residuals."""
    rng = np.random.default_rng(3)
    n_s, n_f = 20, 100
    mag, err, _ = _make_synthetic(n_stars=n_s, n_frames=n_f, seed=3)
    missing = rng.random((n_s, n_f)) < 0.3
    mag[missing] = np.nan
    err[missing] = np.nan

    result = sysrem(mag, err, n_iter=2)
    assert result.residuals is not None
    finite_frac = np.isfinite(result.residuals).mean()
    assert finite_frac > 0.5, f"Too many NaN residuals: finite={finite_frac:.2f}"


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

def test_sysrem_single_star():
    """Single comparison star: should run without error."""
    rng = np.random.default_rng(4)
    mag = rng.normal(0, 0.01, (1, 50))
    err = np.full_like(mag, 0.01)
    result = sysrem(mag, err, n_iter=2)
    assert result.residuals.shape == (1, 50)


def test_sysrem_zero_iter():
    """n_iter=0: no components extracted, residuals equal whitened input."""
    mag, err, _ = _make_synthetic(n_stars=5, n_frames=30, n_systematics=1)
    result = sysrem(mag, err, n_iter=0)
    assert len(result.components) == 0
    assert result.residuals is not None
