"""Regression tests for sampling-window aware period analysis."""

import numpy as np
import pandas as pd
import pytest
from astropy.timeseries import LombScargle

from apex.analysis.light_curve.period_alias_service import (
    analyze_period_aliases,
    classify_frequency_relation,
    compute_spectral_window,
    diagnose_multimode_suspicion,
    fit_multimode_model,
    search_multimode_alias_solutions,
)
from apex.analysis.light_curve.period_io_service import load_period_lightcurve_csv
from apex.analysis.light_curve.detrend_output_service import annotate_step10_output


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


def test_amplitude_modulated_single_mode_is_inconclusive():
    rng = np.random.default_rng(77)
    f0 = 11.6256
    chunks = [np.sort(day + rng.uniform(0.0, 0.15, 45)) for day in range(12)]
    time = np.concatenate(chunks)
    night_id = np.concatenate(
        [np.full(len(chunk), f"n{idx}") for idx, chunk in enumerate(chunks)]
    )
    err = np.full(len(time), 0.006)
    amplitude = 0.18 * (1.0 + 0.45 * np.sin(2.0 * np.pi * 0.18 * time + 0.3))
    mag = amplitude * np.sin(2.0 * np.pi * f0 * time + 0.2)
    mag += rng.normal(0.0, err[0], len(time))
    window = compute_spectral_window(time)
    diagnostic = diagnose_multimode_suspicion(
        time,
        mag,
        err,
        night_id,
        {
            "adopted_freq_cd": f0,
            "status": "RESOLVED",
            "window_peaks": window["peaks"],
        },
        min_period=0.04,
        max_period=0.15,
    )
    assert diagnostic["status"] == "INCONCLUSIVE", diagnostic
    assert "amplitude/phase modulation" in diagnostic["reason"]


def test_long_secondary_mode_reports_night_offset_sensitivity():
    rng = np.random.default_rng(42)
    f0 = 11.6256
    f1 = 1.0 / 2.3
    chunks = [np.sort(day + rng.uniform(0.0, 0.03, 30)) for day in range(12)]
    time = np.concatenate(chunks)
    night_id = np.concatenate(
        [np.full(len(chunk), f"n{idx}") for idx, chunk in enumerate(chunks)]
    )
    err = np.full(len(time), 0.007)
    mag = 0.18 * np.sin(2.0 * np.pi * f0 * time + 0.2)
    mag += 0.08 * np.sin(2.0 * np.pi * f1 * time + 0.7)
    mag += rng.normal(0.0, err[0], len(time))
    window = compute_spectral_window(time)
    diagnostic = diagnose_multimode_suspicion(
        time,
        mag,
        err,
        night_id,
        {
            "adopted_freq_cd": f0,
            "status": "RESOLVED",
            "window_peaks": window["peaks"],
        },
        min_period=0.03,
        max_period=3.0,
    )
    assert diagnostic["status"] == "INCONCLUSIVE", diagnostic
    assert diagnostic["offset_sensitivity"]["sensitive"] is True
    assert diagnostic["offset_sensitivity"]["status_without_offsets"] == "MULTIMODE-SUSPECT"
    assert diagnostic["offset_sensitivity"]["status_with_offsets"] == "SINGLE-COMPATIBLE"


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


def test_truncated_alias_search_cannot_report_resolved():
    time, mag, err, night_id = _two_nights(
        2.0,
        [(11.6256, 0.20, 0.3), (15.03124, 0.05, 0.8)],
        seed=5,
    )
    window = compute_spectral_window(time)
    result = search_multimode_alias_solutions(
        time,
        mag,
        err,
        seed_periods=[1.0 / 11.6256, 1.0 / 15.03124],
        harmonics=1,
        window_peaks=window["peaks"],
        night_id=night_id,
        max_alias_offsets=1,
        max_solutions=1,
    )
    assert result["alias_search_complete"] is False
    assert result["alias_status"] == "AMBIGUOUS"


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


def test_rank_deficient_multimode_fit_is_rejected():
    rng = np.random.default_rng(9)
    time = np.arange(30, dtype=float) + rng.uniform(0.0, 0.02, 30)
    night_id = np.array([f"n{idx}" for idx in range(len(time))])
    err = np.full(len(time), 0.01)
    mag = 0.2 * np.sin(2.0 * np.pi * 11.6256 * time)
    mag += rng.normal(0.0, err[0], len(time))
    with pytest.raises(ValueError, match="Underconstrained|Rank-deficient"):
        fit_multimode_model(
            time,
            mag,
            err,
            periods=[1.0 / 11.6256, 1.0 / 15.03124],
            harmonics=2,
            night_id=night_id,
            include_night_offsets=True,
        )


def test_relative_hour_time_is_converted_to_days(tmp_path):
    hours = np.linspace(0.0, 7.2, 120)
    path = tmp_path / "relative_hours.csv"
    pd.DataFrame(
        {
            "rel_time_hr": hours,
            "mag": np.sin(2.0 * np.pi * hours / 2.4),
            "mag_err": np.full(len(hours), 0.01),
        }
    ).to_csv(path, index=False)
    loaded = load_period_lightcurve_csv(path, "V", 1)
    assert np.ptp(loaded["time"]) == pytest.approx(7.2 / 24.0)
    assert loaded["time_unit"] == "day"
    assert loaded["source_time_unit"] == "hour"


def test_period_loader_combines_filters_after_median_alignment(tmp_path):
    path = tmp_path / "multi_filter.csv"
    pd.DataFrame(
        {
            "JD": [2450000.0, 2450000.1, 2450000.2, 2450000.3],
            "filter": ["g", "g", "r", "r"],
            "diff_mag_raw": [10.0, 10.2, 13.0, 13.4],
            "diff_mag": [9.9, 10.1, 12.8, 13.2],
            "diff_err": [0.01, 0.01, 0.02, 0.02],
        }
    ).to_csv(path, index=False)

    loaded = load_period_lightcurve_csv(path, "all", 1)

    assert loaded["filter"] == "all"
    assert loaded["filter_alignment"] == "per_filter_median"
    assert loaded["filter_values"].tolist() == ["g", "g", "r", "r"]
    assert loaded["mag_raw"] == pytest.approx([-0.1, 0.1, -0.2, 0.2])
    assert loaded["mag_corr"] == pytest.approx([-0.1, 0.1, -0.2, 0.2])


def test_period_loader_keeps_single_filter_zero_point(tmp_path):
    path = tmp_path / "multi_filter.csv"
    pd.DataFrame(
        {
            "JD": [2450000.0, 2450000.1, 2450000.2],
            "filter": ["g", "g", "r"],
            "diff_mag_raw": [10.0, 10.2, 13.0],
        }
    ).to_csv(path, index=False)

    loaded = load_period_lightcurve_csv(path, "g", 1)

    assert loaded["filter"] == "g"
    assert loaded["filter_alignment"] == "none"
    assert loaded["mag_raw"] == pytest.approx([10.0, 10.2])


def test_period_loader_preserves_photometry_provenance(tmp_path):
    path = tmp_path / "psf_lightcurve.csv"
    pd.DataFrame(
        {
            "JD": [2450000.0, 2450000.1],
            "filter": ["g", "g"],
            "diff_mag_raw": [0.01, -0.01],
            "photometry_source": ["psf", "psf"],
            "mag_input_column": ["mag_psf", "mag_psf"],
            "mag_error_input_column": ["mag_psf_err", "mag_psf_err"],
        }
    ).to_csv(path, index=False)

    loaded = load_period_lightcurve_csv(path, "g", 1)

    assert loaded["photometry_source"] == "psf"
    assert loaded["mag_input_column"] == "mag_psf"
    assert loaded["mag_error_input_column"] == "mag_psf_err"


@pytest.mark.parametrize(
    ("mode", "expected"),
    [("offset", False), ("color", False), ("global", True), ("sysrem", True)],
)
def test_step10_output_marks_nightly_baseline_preservation(mode, expected):
    output = annotate_step10_output(
        pd.DataFrame({"JD": [2450000.0], "diff_mag_raw": [0.1]}),
        mode,
        "test",
    )
    assert bool(output["correction_preserves_nightly_baseline"].iloc[0]) is expected


def test_period_loader_marks_offset_corrected_series_as_baseline_destructive(tmp_path):
    path = tmp_path / "lightcurve_ID1_offset.csv"
    pd.DataFrame(
        {
            "JD": [2450000.0, 2450000.1],
            "diff_mag_raw": [0.1, 0.2],
            "diff_mag_corr": [-0.05, 0.05],
            "correction_mode": ["offset", "offset"],
        }
    ).to_csv(path, index=False)
    loaded = load_period_lightcurve_csv(path, "V", 1)
    assert loaded["correction_preserves_nightly_baseline"] is False
