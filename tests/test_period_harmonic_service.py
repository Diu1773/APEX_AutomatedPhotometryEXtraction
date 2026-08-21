"""Choosing between a periodogram peak and its harmonics.

A Lomb-Scargle fits one sinusoid, so an eclipsing binary — two dips per orbit —
peaks at twice its orbital frequency and the search reports half the period.
Measured on ASAS-SN's variable star database, 35 eclipsing binaries of three
kinds: the periodogram peak matched the catalogue period for **none** of them.
Thirty-five pulsators of three kinds went in alongside, because the risk of a
rule like this is that it doubles periods that were already right — none were.

The rule added here doubles the peak when three things hold together, and each
one closes a different way of being wrong:

* the two sub-cycles of the doubled fold differ (there is a secondary eclipse),
* the difference repeats across the run (it is an eclipse, not a second
  unrelated frequency — a synthetic two-mode signal fails here at -0.23 after
  passing the first test at 249 sigma),
* the curve is flat between dips (it is an eclipsing star at all — this is what
  stops a pulsator being doubled).

These tests hold the shape of that rule. The archive measurements themselves
live in `validation/asassn_period_crosscheck.py`; nothing here reaches the
network.
"""

from __future__ import annotations

import numpy as np
import pytest

from apex.analysis.light_curve.period_harmonic_service import (
    choose_bins,
    flat_fraction,
    fold_profile,
    harmonic_periods,
    phase_dispersion,
    resolve_harmonic,
    subcycle_coherence,
    subcycle_difference,
)


def _sampled(period: float, n: int = 400, span: float = 60.0, seed: int = 5):
    """Irregular sampling, the way a survey actually observes."""
    rng = np.random.default_rng(seed)
    return np.sort(rng.uniform(0.0, span, n)), rng


def eclipsing(period=1.4, primary=0.40, secondary=0.18, width=0.035,
              n=400, span=60.0, seed=5, noise=0.006):
    """A detached binary: flat, with a deep dip and a shallower one."""
    t, rng = _sampled(period, n, span, seed)
    phase = np.mod(t / period, 1.0)

    def dip(centre, depth):
        d = np.minimum(np.abs(phase - centre), 1.0 - np.abs(phase - centre))
        return depth * np.exp(-0.5 * (d / width) ** 2)

    mag = dip(0.0, primary) + dip(0.5, secondary) + rng.normal(0, noise, n)
    return t, mag, np.full(n, noise)


def sinusoid(period=1.4, amplitude=0.30, n=400, span=60.0, seed=5, noise=0.006):
    """A pulsator: one brightening and fading per cycle, never flat."""
    t, rng = _sampled(period, n, span, seed)
    mag = amplitude * np.sin(2 * np.pi * t / period) + rng.normal(0, noise, n)
    return t, mag, np.full(n, noise)


# ── the pieces ─────────────────────────────────────────────────────────────

def test_harmonics_are_offered_within_the_search_window():
    got = harmonic_periods(1.0, (2, 3), min_period=0.2, max_period=4.0)
    assert sorted(round(p, 6) for p in got) == [0.333333, 0.5, 1.0, 2.0, 3.0]
    # Outside the window they are not offered at all — a candidate the search
    # could never return is not a candidate.
    assert harmonic_periods(1.0, (2,), min_period=0.9, max_period=1.1) == [1.0]
    assert harmonic_periods(float("nan"), (2,)) == []


def test_bins_follow_the_data_not_a_constant():
    """The bin count decides the answer, so it cannot be a fixed 20.

    A detached binary's eclipse is a few per cent of the cycle; at twenty bins
    the eclipse shares a bin with the flat part and folding at half the period
    barely changes anything. On ASAS-SN J234910 that flipped the conclusion.
    """
    assert choose_bins(40) < choose_bins(400) < choose_bins(4000)
    assert choose_bins(3) >= 10          # never below the coarsest rung
    for n in (50, 200, 800):
        assert choose_bins(n) * 5 <= n   # at least five points a bin


def test_folding_at_the_right_period_beats_folding_at_a_wrong_one():
    t, mag, _ = sinusoid(period=1.4)
    assert phase_dispersion(t, mag, 1.4) < 0.2
    assert phase_dispersion(t, mag, 1.4 * 1.3) > 0.6


def test_sub_cycle_significance_alone_is_not_a_discriminator():
    """It fires on a clean sinusoid too, which is why it is not used alone.

    A bin's standard error describes the scatter of points inside it, not the
    slope of the curve across it. Where the curve is steep, which half-cycle a
    point lands in shifts the bin mean by more than that error allows for, so a
    pure sinusoid folded at twice its period shows a sub-cycle difference of
    about 8 sigma with nothing wrong.

    The binary's signal is far larger, but "far larger" is not a threshold. The
    flatness condition is what actually separates them, and this test records
    why it has to exist.
    """
    t, mag, err = eclipsing(period=1.4)
    binary = subcycle_difference(t, mag, 1.4, 2, mag_err=err)
    assert binary["sigma"] > 3.0, binary

    t2, mag2, err2 = sinusoid(period=0.7)
    pulsator = subcycle_difference(t2, mag2, 1.4, 2, mag_err=err2)
    assert pulsator["sigma"] > 3.0, (
        "if a clean sinusoid ever stops tripping this, the flatness gate can "
        "be reconsidered — until then it is load-bearing")
    assert binary["sigma"] > pulsator["sigma"]


def test_flatness_tells_an_eclipsing_shape_from_a_pulsating_one():
    t, mag, err = eclipsing(period=1.4)
    t2, mag2, err2 = sinusoid(period=1.4)
    assert flat_fraction(t, mag, 1.4, mag_err=err) > 0.6
    assert flat_fraction(t2, mag2, 1.4, mag_err=err2) < 0.5


def test_the_eclipse_depth_diagnostic_is_gone():
    """It promised what it never delivered, over three rewrites.

    "How unequal are the two eclipses" was reported as 21 sigma for a star whose
    eclipses agree to 1.2, then as noise bumps, then smaller at half the period
    than at the period — the opposite of its premise. Nothing in the decision
    used it, so it was a number a reader would have trusted for nothing.
    """
    import apex.analysis.light_curve.period_harmonic_service as service

    assert not hasattr(service, "half_cycle_contrast")
    assert not hasattr(service, "alternating_minima")
    t, mag, err = eclipsing(period=1.4)
    verdict = resolve_harmonic(t, mag, 0.7, mag_err=err,
                               min_period=0.1, max_period=6.0)
    assert "alternating_sigma" not in verdict.as_dict()


def test_coherence_separates_an_eclipse_from_a_second_frequency():
    t, mag, err = eclipsing(period=1.4)
    assert subcycle_coherence(t, mag, 1.4, 2, err) > 0.5

    # Two unrelated frequencies. Their sub-cycles differ too, but the pattern
    # drifts, so the first half of the run and the second do not agree.
    rng = np.random.default_rng(11)
    t2 = np.sort(rng.uniform(0, 60, 400))
    mag2 = (0.20 * np.sin(2 * np.pi * t2 / 0.7)
            + 0.05 * np.sin(2 * np.pi * t2 / 0.543 + 0.8)
            + rng.normal(0, 0.006, t2.size))
    assert subcycle_coherence(t2, mag2, 1.4, 2, np.full(t2.size, 0.006)) < 0.5


# ── the rule ───────────────────────────────────────────────────────────────

def test_an_eclipsing_binary_is_doubled():
    t, mag, err = eclipsing(period=1.4)
    verdict = resolve_harmonic(t, mag, 0.7, mag_err=err,
                               min_period=0.1, max_period=6.0)
    assert verdict.factor == pytest.approx(2.0)
    assert verdict.adopted_period == pytest.approx(1.4)
    assert "flat between dips" in verdict.reason


def test_a_pulsator_is_left_alone():
    """The error that matters more: a correct period must not be doubled."""
    t, mag, err = sinusoid(period=0.7)
    verdict = resolve_harmonic(t, mag, 0.7, mag_err=err,
                               min_period=0.1, max_period=6.0)
    assert verdict.factor == pytest.approx(1.0)
    assert verdict.adopted_period == pytest.approx(0.7)
    assert not verdict.changed


def test_two_unrelated_frequencies_are_not_doubled():
    """249 sigma of sub-cycle difference, and still not a doubled period.

    This is the case that made significance alone unusable: a second pulsation
    frequency makes the sub-cycles differ more than any real binary measured.
    """
    rng = np.random.default_rng(3)
    t = np.sort(rng.uniform(0, 60, 400))
    mag = (0.20 * np.sin(2 * np.pi * t / 0.7)
           + 0.05 * np.sin(2 * np.pi * t / 0.543 + 0.8)
           + rng.normal(0, 0.006, t.size))
    verdict = resolve_harmonic(t, mag, 0.7, mag_err=np.full(t.size, 0.006),
                               min_period=0.1, max_period=6.0)
    assert verdict.factor == pytest.approx(1.0), verdict.reason


def test_only_doubling_is_adopted_however_deep_the_evidence():
    """x3 was tried, misfired, and is not adopted — only listed.

    On ASAS-SN J225021, whose two eclipses are indistinguishable and whose
    doubling test correctly failed, the x3 fold passed and tripled the period.
    """
    from apex.analysis.light_curve.period_harmonic_service import DEFAULT_MULTIPLES

    assert DEFAULT_MULTIPLES == (2,)
    t, mag, err = eclipsing(period=1.4)
    verdict = resolve_harmonic(t, mag, 0.7, mag_err=err,
                               min_period=0.1, max_period=6.0)
    assert verdict.factor in (pytest.approx(1.0), pytest.approx(2.0))


def test_the_verdict_carries_its_evidence():
    t, mag, err = eclipsing(period=1.4)
    verdict = resolve_harmonic(t, mag, 0.7, mag_err=err,
                               min_period=0.1, max_period=6.0)
    body = verdict.as_dict()
    assert set(body) >= {"adopted_period", "base_period", "factor", "reason",
                         "candidates"}
    doubled = next(r for r in verdict.candidates if r["factor"] == pytest.approx(2.0))
    assert doubled["subcycle_sigma"] > 3.0
    assert doubled["flat_fraction"] > 0.6
    assert "subcycle_coherence" in doubled


def test_too_few_points_returns_the_periodogram_period_unchanged():
    t, mag, err = eclipsing(period=1.4, n=12, span=6.0)
    verdict = resolve_harmonic(t, mag, 0.7, mag_err=err,
                               min_period=0.1, max_period=6.0)
    assert verdict.adopted_period == pytest.approx(0.7)
    assert not verdict.changed


def test_the_alias_candidate_list_now_contains_the_harmonics():
    """Adopting and listing are different acts; the table must hold the answer.

    On three of five ASAS-SN binaries the catalogue period was absent from the
    candidate table entirely, because the list was built from the periodogram
    peak and the sampling-window offsets alone.
    """
    from apex.analysis.light_curve.period_alias_service import (
        build_alias_candidates, compute_spectral_window,
    )
    from astropy.timeseries import LombScargle

    t, mag, err = eclipsing(period=1.4)
    freq, power = LombScargle(t, mag, err).autopower(
        minimum_frequency=1 / 6.0, maximum_frequency=1 / 0.1, samples_per_peak=10)
    candidates = build_alias_candidates(
        t, freq, power, min_period=0.1, max_period=6.0,
        window=compute_spectral_window(t), max_candidates=12)
    periods = np.array([c["period"] for c in candidates], dtype=float)
    assert periods.size
    assert (np.abs(periods - 1.4) / 1.4 < 0.05).any(), sorted(periods.round(4))
