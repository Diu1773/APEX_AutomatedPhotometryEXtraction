"""Regression tests for sampling-window aware period analysis."""

import numpy as np
from astropy.timeseries import LombScargle

from apex.analysis.light_curve.period_alias_service import (
    analyze_period_aliases,
    classify_frequency_relation,
    compute_spectral_window,
    diagnose_multimode_suspicion,
    fit_multimode_model,
    search_multimode_alias_solutions,
)


def _two_nights(gap_days, frequencies, spans=(0.14, 0.14), seed=7, noise=0.007):
    rng = np.random.default_rng(seed)
    t1 = np.sort(rng.uniform(0.0, spans[0], 100))
    t2 = np.sort(rng.uniform(gap_days, gap_days + spans[1], 100))
    time = np.concatenate([t1, t2])
    mag = np.zeros_like(time)
    for freq, amplitude, phase in frequencies:
        mag += amplitude * np.sin(2.0 * np.pi * freq * time + phase)
    mag += rng.normal(0.0, noise, len(time))
    err = np.full(len(time), noise)
    night_id = np.array(["n1"] * len(t1) + ["n2"] * len(t2))
    return time, mag, err, night_id


def _ls(time, mag, err):
    ls = LombScargle(time, mag, err)
    return ls.autopower(
        minimum_frequency=3.0,
        maximum_frequency=35.0,
        samples_per_peak=30,
    )


def test_two_day_gap_window_contains_half_cycle_per_day_peak():
    time, mag, err, _ = _two_nights(2.0, [(11.6256, 0.2, 0.3)])
    window = compute_spectral_window(time)
    offsets = [row["freq_cd"] for row in window["peaks"]]
    assert any(abs(offset - 0.5) < 0.08 for offset in offsets), offsets


def test_ae_like_alias_family_includes_half_day_alias():
    f0 = 11.6256
    time, mag, err, night_id = _two_nights(2.0, [(f0, 0.2, 0.3)])
    frequency, power = _ls(time, mag, err)
    result = analyze_period_aliases(
        time,
        mag,
        err,
        night_id,
        frequency,
        power,
        min_period=0.03,
        max_period=0.30,
        n_injections=4,
    )
    candidate_freqs = [row["freq_cd"] for row in result["candidates"]]
    assert any(abs(freq - f0) < 0.08 for freq in candidate_freqs), candidate_freqs
    assert any(abs(freq - (f0 + 0.5)) < 0.10 for freq in candidate_freqs), candidate_freqs
    assert any(row["relation_to_best"] == "window-alias" for row in result["candidates"][1:])


def test_single_mode_diagnostic_rejects_noise_residuals():
    time, mag, err, night_id = _two_nights(2.0, [(11.6256, 0.2, 0.3)], seed=3)
    frequency, power = _ls(time, mag, err)
    alias = analyze_period_aliases(
        time,
        mag,
        err,
        night_id,
        frequency,
        power,
        min_period=1.0 / 35.0,
        max_period=1.0 / 3.0,
        n_injections=2,
    )
    diagnostic = diagnose_multimode_suspicion(
        time,
        mag,
        err,
        night_id,
        alias,
        min_period=1.0 / 35.0,
        max_period=1.0 / 3.0,
    )
    assert diagnostic["status"] == "SINGLE-COMPATIBLE", diagnostic


def test_multimode_diagnostic_finds_independent_residual_mode():
    time, mag, err, night_id = _two_nights(
        2.0,
        [(11.6256, 0.20, 0.3), (15.03124, 0.05, 0.8)],
        seed=3,
    )
    frequency, power = _ls(time, mag, err)
    alias = analyze_period_aliases(
        time,
        mag,
        err,
        night_id,
        frequency,
        power,
        min_period=1.0 / 35.0,
        max_period=1.0 / 3.0,
        n_injections=2,
    )
    diagnostic = diagnose_multimode_suspicion(
        time,
        mag,
        err,
        night_id,
        alias,
        min_period=1.0 / 35.0,
        max_period=1.0 / 3.0,
    )
    assert diagnostic["status"] == "MULTIMODE-SUSPECT", diagnostic
    assert abs(float(diagnostic["candidate_frequency_cd"]) - 15.03124) < 0.08


def test_short_yz_like_nights_are_not_auto_resolved():
    f0 = 9.6069
    time, mag, err, night_id = _two_nights(1.0, [(f0, 0.19, 0.5)])
    frequency, power = _ls(time, mag, err)
    result = analyze_period_aliases(
        time,
        mag,
        err,
        night_id,
        frequency,
        power,
        min_period=0.03,
        max_period=0.30,
        n_injections=4,
    )
    assert result["status"] in {"AMBIGUOUS", "INSUFFICIENT"}, result
    assert result["longest_night_cycles"] < 1.5


def test_candidate_frequency_is_refined_beyond_the_ls_grid():
    rng = np.random.default_rng(23)
    true_frequency = 9.6069
    time = np.sort(rng.uniform(0.0, 0.22, 120))
    mag = 0.18 * np.sin(2.0 * np.pi * true_frequency * time + 0.2)
    mag += 0.04 * np.sin(4.0 * np.pi * true_frequency * time + 0.9)
    mag += rng.normal(0.0, 0.006, len(time))
    err = np.full(len(time), 0.006)
    ls = LombScargle(time, mag, err)
    frequency, power = ls.autopower(
        minimum_frequency=6.0,
        maximum_frequency=15.0,
        samples_per_peak=5,
    )
    raw_frequency = float(frequency[int(np.argmax(power))])
    result = analyze_period_aliases(
        time,
        mag,
        err,
        np.full(len(time), "n1"),
        frequency,
        power,
        min_period=1.0 / 15.0,
        max_period=1.0 / 6.0,
        n_injections=0,
    )
    refined_frequency = float(result["adopted_freq_cd"])
    assert abs(refined_frequency - true_frequency) < abs(raw_frequency - true_frequency)
    assert abs(refined_frequency - true_frequency) < 0.08


def test_many_short_nights_can_support_a_period_longer_than_one_day():
    rng = np.random.default_rng(11)
    period = 2.3
    chunks = [np.sort(day + rng.uniform(0.0, 0.16, 35)) for day in range(12)]
    time = np.concatenate(chunks)
    night_id = np.concatenate(
        [np.full(len(chunk), f"n{idx}") for idx, chunk in enumerate(chunks)]
    )
    err = np.full(len(time), 0.008)
    mag = 0.15 * np.sin(2.0 * np.pi * time / period + 0.4)
    mag += rng.normal(0.0, err[0], len(time))
    ls = LombScargle(time, mag, err)
    frequency, power = ls.autopower(
        minimum_frequency=0.2,
        maximum_frequency=1.5,
        samples_per_peak=30,
    )
    result = analyze_period_aliases(
        time,
        mag,
        err,
        night_id,
        frequency,
        power,
        min_period=0.5,
        max_period=5.0,
        n_injections=2,
    )
    assert result["global_cycles"] > 1.0
    assert result["occupied_fraction"] >= 0.35
    assert result["status"] != "INSUFFICIENT", result


def test_frequency_classification_uses_measured_window():
    relation, note = classify_frequency_relation(
        12.1256,
        [11.6256],
        window_peaks=[{"freq_cd": 0.5, "power": 0.95}],
        baseline_days=2.1,
    )
    assert relation == "alias"
    assert "0.5000" in note


def test_multimode_alias_search_is_joint_and_conservative():
    f0, f1 = 11.6256, 15.03124
    time, mag, err, night_id = _two_nights(
        2.0,
        [(f0, 0.20, 0.3), (f1, 0.05, 0.8)],
        seed=3,
    )
    window = compute_spectral_window(time)
    result = search_multimode_alias_solutions(
        time,
        mag,
        err,
        seed_periods=[1.0 / (f0 + 0.5), 1.0 / (f1 + 0.5)],
        harmonics=2,
        window_peaks=window["peaks"],
        night_id=night_id,
    )
    assert len(result["periods"]) == 2
    assert len(result["alias_solutions"]) > 1
    assert result["alias_status"] in {"AMBIGUOUS", "INSUFFICIENT"}

    fixed = fit_multimode_model(
        time,
        mag,
        err,
        periods=[1.0 / f0, 1.0 / f1],
        harmonics=2,
    )
    single = fit_multimode_model(time, mag, err, periods=[1.0 / f0], harmonics=2)
    assert fixed["rmse"] < single["rmse"]


def test_multimode_fit_can_include_nightly_zero_point_terms():
    f0, f1 = 11.6256, 15.03124
    time, mag, err, night_id = _two_nights(
        2.0,
        [(f0, 0.20, 0.3), (f1, 0.05, 0.8)],
        seed=19,
        noise=0.004,
    )
    mag = mag + np.where(night_id == "n2", 0.08, 0.0)
    without_offsets = fit_multimode_model(
        time,
        mag,
        err,
        periods=[1.0 / f0, 1.0 / f1],
        harmonics=2,
    )
    with_offsets = fit_multimode_model(
        time,
        mag,
        err,
        periods=[1.0 / f0, 1.0 / f1],
        harmonics=2,
        night_id=night_id,
        include_night_offsets=True,
    )
    assert with_offsets["rmse"] < 0.5 * without_offsets["rmse"]
    assert abs(with_offsets["night_offsets"]["n2"] - 0.08) < 0.02
