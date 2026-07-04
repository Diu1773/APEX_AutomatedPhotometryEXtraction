"""Algorithmic validation of iterative pre-whitening (Period04/Breger scheme).

These tests use clean synthetic signals (no spectral-window aliasing) so they
isolate the pre-whitening ALGORITHM from the separate aliasing problem: a
double-mode signal must yield both frequencies, a fundamental+harmonic signal
must yield both, and pure noise must yield nothing significant.
"""

import numpy as np

from apex.analysis.light_curve.period_analysis_service import prewhiten_frequencies


def _sample(freqs_amps, n=400, span=6.0, seed=1, noise=0.005):
    rng = np.random.default_rng(seed)
    t = np.sort(rng.uniform(0, span, n))
    y = np.zeros_like(t)
    for f, a, ph in freqs_amps:
        y += a * np.sin(2 * np.pi * f * t + ph)
    y += rng.normal(0, noise, n)
    dy = np.full(n, noise)
    return t, y, dy


def _closest(found, f_true, tol=0.05):
    return any(abs(g["freq_cd"] - f_true) < tol for g in found)


def test_double_mode_recovers_both_frequencies():
    # f0 and first overtone (P1/P0 ~ 0.77 -> f1/f0 ~ 1.30), like a double-mode HADS
    f0, f1 = 11.6256, 15.10
    t, y, dy = _sample([(f0, 0.20, 0.1), (f1, 0.12, 0.7)], noise=0.004)
    r = prewhiten_frequencies(t, y, dy, min_period=0.03, max_period=0.5, sn_stop=4.0)
    found = r["frequencies"]
    assert _closest(found, f0), f"f0 not recovered: {[round(g['freq_cd'],3) for g in found]}"
    assert _closest(found, f1), f"f1 not recovered: {[round(g['freq_cd'],3) for g in found]}"
    # amplitudes ordered so the fundamental (largest) is found first
    assert abs(found[0]["freq_cd"] - f0) < 0.05
    assert abs(found[0]["amplitude"] - 0.20) < 0.03


def test_fundamental_plus_harmonic():
    f0 = 9.6069  # YZ Boo fundamental
    t, y, dy = _sample([(f0, 0.19, 0.6), (2 * f0, 0.06, 0.5), (3 * f0, 0.02, 0.8)],
                       n=500, noise=0.004)
    r = prewhiten_frequencies(t, y, dy, min_period=0.02, max_period=0.5, sn_stop=4.0)
    found = r["frequencies"]
    assert _closest(found, f0)
    assert _closest(found, 2 * f0)


def test_pure_noise_yields_nothing_significant():
    rng = np.random.default_rng(3)
    t = np.sort(rng.uniform(0, 6, 300))
    y = rng.normal(0, 0.01, 300)
    dy = np.full(300, 0.01)
    r = prewhiten_frequencies(t, y, dy, min_period=0.03, max_period=0.5, sn_stop=4.0)
    # a couple of spurious low-S/N peaks may appear but none should be strong
    assert r["n_significant"] <= 2


def test_too_few_points_errors_gracefully():
    r = prewhiten_frequencies(np.arange(5.0), np.zeros(5), None, 0.03, 0.5)
    assert r["n_significant"] == 0 and "error" in r
