from __future__ import annotations

import numpy as np
import pandas as pd

from apex.analysis.light_curve.variable_analysis_contract import (
    ValidatedLightCurveBundle,
    VariableAnalysisRequest,
)
from apex.analysis.light_curve.variable_analysis_cache import (
    load_cached_result,
    request_cache_key,
    store_cached_result,
)
from apex.analysis.light_curve.variable_analysis_service import (
    refine_period_local,
    run_variable_analysis,
)


def _write_lightcurve(tmp_path, time, mag, mag_err, filters=None):
    if filters is None:
        filters = np.full(len(time), "g")
    path = tmp_path / "lightcurve_ID1_raw.csv"
    pd.DataFrame(
        {
            "ID": np.ones(len(time), dtype=int),
            "BJD_TDB": time,
            "filter": filters,
            "diff_mag_raw": mag,
            "diff_err": mag_err,
            "correction_mode": "raw",
        }
    ).to_csv(path, index=False)
    return path


def _bundle(path, time, *, alias_status="RESOLVED", multimode=None, release="APPROVED"):
    return ValidatedLightCurveBundle(
        workspace_dir=str(path.parent),
        source_file=str(path),
        target_id=1,
        analysis_filter="__all__",
        series_mode="raw",
        mag_col="diff_mag_raw",
        correction_mode="raw",
        correction_preserves_nightly_baseline=True,
        input_signature={
            "n_points": int(len(time)),
            "time_min": float(np.min(time)),
            "time_max": float(np.max(time)),
        },
        adopted_period=0.104,
        scan_results={},
        alias_analysis={
            "status": alias_status,
            "candidates": [
                {"rank": 1, "period": 0.104, "relative_power": 1.0},
                {"rank": 2, "period": 0.116, "relative_power": 0.8},
            ],
        },
        multimode_diagnostic=multimode or {"status": "SINGLE-COMPATIBLE"},
        search={"min_period": 0.04, "max_period": 0.18},
        release_status=release,
        main_qc={"status": "PASS"},
    )


def test_single_mode_service_refines_and_fits_each_filter(tmp_path):
    rng = np.random.default_rng(3)
    period = 0.10402
    time = np.sort(rng.uniform(0.0, 4.0, 240)) + 2460000.0
    filters = np.where(np.arange(len(time)) % 2 == 0, "g", "r")
    amplitudes = np.where(filters == "g", 0.18, 0.11)
    offsets = np.where(filters == "g", 0.0, 0.35)
    error = np.full(len(time), 0.008)
    mag = offsets + amplitudes * np.cos(2.0 * np.pi * time / period)
    mag += rng.normal(0.0, error)
    path = _write_lightcurve(tmp_path, time, mag, error, filters)

    result = run_variable_analysis(
        VariableAnalysisRequest(
            bundle=_bundle(path, time),
            bootstrap_resamples=12,
            single_harmonics=2,
            random_seed=7,
        )
    )

    assert result.status == "COMPLETE"
    assert result.branch == "single"
    assert abs(result.refined_period - period) < 2e-4
    assert set(result.per_filter_models) == {"g", "r"}
    assert result.per_filter_models["g"]["rmse"] < 0.02
    assert result.per_filter_models["r"]["rmse"] < 0.02
    assert result.diagnostics["main_qc_recomputed"] is False
    assert result.refinement["uncertainty_scope"] == "local_conditional_on_selected_alias"
    assert result.data_summary["observed_cycles"] > 30
    assert result.diagnostics["limited_cycle_coverage"] is False


def test_local_refinement_handles_short_asymmetric_multiband_waveform():
    rng = np.random.default_rng(99)
    period = 0.10402
    time = np.sort(rng.uniform(0.0, 0.23, 240)) + 2460000.0
    filters = np.where(np.arange(len(time)) % 2, "g", "r")
    amplitude_scale = np.where(filters == "g", 1.0, 0.7)
    phase = 2.0 * np.pi * time / period
    mag = amplitude_scale * (
        0.15 * np.cos(phase)
        + 0.06 * np.cos(2.0 * phase + 0.8)
        + 0.025 * np.sin(3.0 * phase - 0.3)
    )
    error = np.full(len(time), 0.005)

    result = refine_period_local(
        time,
        mag,
        error,
        0.1036,
        filter_values=filters,
        alias_frequency_offsets=[1.0],
        bootstrap_resamples=0,
        harmonics=4,
    )

    assert abs(result["refined_period"] - period) < 1e-5
    assert result["refinement_method"] == "per_filter_multi_harmonic_fourier"
    assert result["frequency_half_width_cd"] <= 0.45


def test_unresolved_alias_returns_review_before_reading_source(tmp_path):
    missing = tmp_path / "missing.csv"
    time = np.array([1.0, 2.0])
    bundle = _bundle(missing, time, alias_status="AMBIGUOUS")

    result = run_variable_analysis(
        VariableAnalysisRequest(bundle=bundle, bootstrap_resamples=0)
    )

    assert result.status == "REVIEW_REQUIRED"
    assert result.review.code == "ALIAS_SELECTION_REQUIRED"
    assert len(result.review.candidates) == 2


def test_multimode_service_uses_independent_second_period(tmp_path):
    rng = np.random.default_rng(11)
    primary = 0.10401
    secondary = 0.0713
    time = np.sort(rng.uniform(0.0, 5.0, 280)) + 2460100.0
    error = np.full(len(time), 0.006)
    mag = 0.16 * np.cos(2.0 * np.pi * time / primary)
    mag += 0.055 * np.sin(2.0 * np.pi * time / secondary)
    mag += rng.normal(0.0, error)
    path = _write_lightcurve(tmp_path, time, mag, error)
    multimode = {
        "status": "MULTIMODE-SUSPECT",
        "candidate_period": secondary,
        "candidate_delta_bic": 40.0,
        "candidate_fap": 1e-6,
    }

    result = run_variable_analysis(
        VariableAnalysisRequest(
            bundle=_bundle(path, time, multimode=multimode),
            bootstrap_resamples=8,
            multimode_harmonics=1,
            random_seed=5,
        )
    )

    assert result.status == "COMPLETE"
    assert result.branch == "multi"
    assert result.model["kind"] == "per_filter_multimode"
    assert abs(result.model["periods"][0] - primary) < 2e-4
    assert result.model["periods"][1] == secondary
    assert result.per_filter_models["g"]["rmse"] < 0.02


def test_approved_source_signature_change_is_rejected(tmp_path):
    time = np.linspace(1.0, 3.0, 60)
    error = np.full(len(time), 0.01)
    mag = 0.1 * np.sin(2.0 * np.pi * time / 0.104)
    path = _write_lightcurve(tmp_path, time, mag, error)
    bundle = _bundle(path, time)
    frame = pd.read_csv(path).iloc[:-1]
    frame.to_csv(path, index=False)

    result = run_variable_analysis(
        VariableAnalysisRequest(bundle=bundle, bootstrap_resamples=0)
    )

    assert result.status == "FAILED"
    assert "changed after Main-workflow validation" in result.error


def test_blocked_bundle_never_runs_advanced_analysis(tmp_path):
    missing = tmp_path / "missing.csv"
    bundle = _bundle(missing, np.array([1.0, 2.0]), release="BLOCKED")
    bundle.release_reasons = ["Check-star QC failed."]

    result = run_variable_analysis(VariableAnalysisRequest(bundle=bundle))

    assert result.status == "FAILED"
    assert "blocked" in result.error.lower()


def test_cache_round_trip_and_source_change_invalidate_key(tmp_path):
    time = np.linspace(1.0, 3.0, 80)
    error = np.full(len(time), 0.01)
    mag = 0.1 * np.sin(2.0 * np.pi * time / 0.104)
    path = _write_lightcurve(tmp_path, time, mag, error)
    request = VariableAnalysisRequest(
        bundle=_bundle(path, time),
        bootstrap_resamples=0,
        single_harmonics=1,
    )
    result = run_variable_analysis(request)
    initial_key = request_cache_key(request)

    stored = store_cached_result(request, result)
    restored = load_cached_result(request)

    assert stored.is_file()
    assert restored.status == "COMPLETE"
    assert isinstance(restored.per_filter_models["g"]["coeff"], np.ndarray)

    frame = pd.read_csv(path)
    frame.loc[0, "diff_mag_raw"] += 0.02
    frame.to_csv(path, index=False)

    assert request_cache_key(request) != initial_key
    assert load_cached_result(request) is None
