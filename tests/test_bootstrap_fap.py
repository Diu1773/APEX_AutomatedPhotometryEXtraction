"""Tests for bootstrap_fap in period_analysis_service."""
from __future__ import annotations

import numpy as np
import pytest

from apex.analysis.light_curve.period_analysis_service import bootstrap_fap, compute_ls


def _sine_lc(period=3.14, n=300, noise=0.01, seed=0):
    rng = np.random.default_rng(seed)
    t = np.sort(rng.uniform(0, 30, n))
    m = 0.1 * np.sin(2 * np.pi * t / period) + rng.normal(0, noise, n)
    e = np.full(n, noise)
    return t, m, e


def _noise_lc(n=200, seed=1):
    rng = np.random.default_rng(seed)
    t = np.sort(rng.uniform(0, 20, n))
    m = rng.normal(0, 0.01, n)
    e = np.full(n, 0.01)
    return t, m, e


# ---------------------------------------------------------------------------

def test_bootstrap_fap_real_signal_low():
    """Strong periodic signal should yield a low FAP."""
    t, m, e = _sine_lc()
    ls = compute_ls(t, m, e, "test", 0.5, 20.0, 10)
    fap = bootstrap_fap(t, m, e, ls["best_power"], 0.5, 20.0,
                        n_bootstrap=200, seed=42)
    assert np.isfinite(fap)
    assert fap < 0.05, f"FAP={fap:.4f} expected < 0.05 for strong sinusoid"


def test_bootstrap_fap_pure_noise_high():
    """Pure noise should yield FAP close to 1 (or at least > 0.1)."""
    t, m, e = _noise_lc()
    ls = compute_ls(t, m, e, "test", 0.5, 15.0, 10)
    fap = bootstrap_fap(t, m, e, ls["best_power"], 0.5, 15.0,
                        n_bootstrap=300, seed=7)
    assert np.isfinite(fap)
    assert fap > 0.05, f"FAP={fap:.4f} expected high for pure noise"


def test_bootstrap_fap_returns_fraction_in_01():
    t, m, e = _sine_lc()
    ls = compute_ls(t, m, e, "test", 0.5, 20.0, 10)
    fap = bootstrap_fap(t, m, e, ls["best_power"], 0.5, 20.0,
                        n_bootstrap=100, seed=0)
    assert 0.0 <= fap <= 1.0


def test_bootstrap_fap_too_few_points():
    t = np.array([1.0, 2.0, 3.0])
    m = np.array([0.1, -0.1, 0.05])
    e = np.array([0.01, 0.01, 0.01])
    fap = bootstrap_fap(t, m, e, 0.5, 0.5, 5.0, n_bootstrap=10)
    assert not np.isfinite(fap)


def test_bootstrap_fap_progress_callback():
    t, m, e = _sine_lc(n=100)
    calls = []
    ls = compute_ls(t, m, e, "test", 0.5, 15.0, 5)
    bootstrap_fap(t, m, e, ls["best_power"], 0.5, 15.0,
                  n_bootstrap=20, seed=1,
                  progress_cb=lambda cur, tot: calls.append((cur, tot)))
    assert len(calls) == 20
    assert calls[-1] == (20, 20)


def test_bootstrap_fap_reproducible_with_seed():
    t, m, e = _sine_lc()
    ls = compute_ls(t, m, e, "test", 0.5, 20.0, 10)
    fap1 = bootstrap_fap(t, m, e, ls["best_power"], 0.5, 20.0,
                          n_bootstrap=100, seed=99)
    fap2 = bootstrap_fap(t, m, e, ls["best_power"], 0.5, 20.0,
                          n_bootstrap=100, seed=99)
    assert fap1 == fap2


def test_bootstrap_fap_higher_power_lower_fap():
    """Artificially inflating the threshold power should raise FAP."""
    t, m, e = _sine_lc()
    ls = compute_ls(t, m, e, "test", 0.5, 20.0, 10)
    bp = ls["best_power"]
    fap_low  = bootstrap_fap(t, m, e, bp,         0.5, 20.0, n_bootstrap=200, seed=5)
    fap_high = bootstrap_fap(t, m, e, bp * 10.0,  0.5, 20.0, n_bootstrap=200, seed=5)
    assert fap_low <= fap_high + 1e-9
