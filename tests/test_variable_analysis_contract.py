from __future__ import annotations

import json

import numpy as np

from apex.analysis.light_curve.variable_analysis_contract import (
    BUNDLE_SCHEMA,
    ReviewRequired,
    ValidatedLightCurveBundle,
    VariableAnalysisRequest,
    VariableAnalysisResult,
    coerce_validated_bundle,
)


def _bundle(**overrides) -> ValidatedLightCurveBundle:
    values = {
        "workspace_dir": r"C:\workspace",
        "source_file": r"C:\workspace\lightcurve_ID1_raw.csv",
        "target_id": 1,
        "analysis_filter": "__all__",
        "series_mode": "raw",
        "mag_col": "diff_mag_raw",
        "correction_mode": "raw",
        "correction_preserves_nightly_baseline": True,
        "input_signature": {"n_points": 120, "time_min": 1.0, "time_max": 2.0},
        "adopted_period": 0.104,
        "scan_results": {
            "raw_ls": {
                "frequency": np.array([8.0, 9.0, 10.0]),
                "power": np.array([0.1, 0.3, 0.9]),
                "best_period": 0.104,
                "best_power": 0.9,
            }
        },
        "alias_analysis": {"status": "RESOLVED", "adopted_period": 0.104},
        "multimode_diagnostic": {"status": "SINGLE-COMPATIBLE"},
        "search": {"min_period": 0.04, "max_period": 0.18, "methods": ["ls"]},
        "release_status": "APPROVED",
        "main_qc": {"check_star": {"status": "AVAILABLE", "check_id": 14}},
        "comparison_provenance": {"filters": {"g": {"comparison_ids": [5, 6]}}},
        "photometry_provenance": {"source": "aperture", "mag_column": "mag"},
    }
    values.update(overrides)
    return ValidatedLightCurveBundle(**values)


def test_bundle_json_round_trip_restores_runtime_period_arrays(tmp_path):
    path = tmp_path / "bundle.json"
    _bundle().write_json(path)

    serialized = json.loads(path.read_text(encoding="utf-8"))
    restored = ValidatedLightCurveBundle.read_json(path)
    legacy = restored.to_legacy_handoff()

    assert serialized["schema"] == BUNDLE_SCHEMA
    assert restored.release_status == "APPROVED"
    assert restored.main_qc["check_star"]["check_id"] == 14
    assert isinstance(legacy["scan_results"]["raw_ls"]["frequency"], np.ndarray)
    assert np.allclose(legacy["scan_results"]["raw_ls"]["power"], [0.1, 0.3, 0.9])


def test_legacy_handoff_is_accepted_as_unverified_bundle():
    legacy = {
        "workspace_dir": r"C:\workspace",
        "source_file": "lightcurve_ID1_raw.csv",
        "target_id": 1,
        "analysis_filter": "g",
        "series_mode": "raw",
        "mag_col": "diff_mag_raw",
        "correction_mode": "raw",
        "input_signature": {"n_points": 20},
        "adopted_period": 0.104,
        "scan_results": {},
        "alias_analysis": {},
        "multimode_diagnostic": {},
        "search": {},
    }

    bundle = coerce_validated_bundle(legacy)

    assert bundle.release_status == "UNVERIFIED"
    assert bundle.can_launch
    assert bundle.to_legacy_handoff()["analysis_filter"] == "g"


def test_blocked_bundle_cannot_launch():
    bundle = _bundle(
        release_status="BLOCKED",
        release_reasons=["Check-star stability result is missing."],
    )

    assert not bundle.can_launch
    assert "Check-star" in bundle.release_message


def test_request_and_result_are_json_safe():
    request = VariableAnalysisRequest(
        bundle=_bundle(),
        analysis_branch="single",
        bootstrap_resamples=20,
        random_seed=42,
    )
    result = VariableAnalysisResult(
        status="REVIEW_REQUIRED",
        adopted_period=0.104,
        refinement={"fine_power": np.array([0.1, 0.8])},
        review=ReviewRequired(
            code="ALIAS_SELECTION_REQUIRED",
            message="Select a candidate.",
            candidates=[{"period": np.float64(0.104)}],
            allowed_actions=["adopt_period_candidate"],
        ),
    )

    request_payload = request.to_dict(json_safe=True)
    result_payload = result.to_dict(json_safe=True)

    json.dumps(request_payload, allow_nan=False)
    json.dumps(result_payload, allow_nan=False)
    restored = VariableAnalysisResult.from_dict(result_payload)
    assert request_payload["analysis_branch"] == "single"
    assert result_payload["review"]["candidates"][0]["period"] == 0.104
    assert isinstance(restored.refinement["fine_power"], np.ndarray)
    assert restored.review.code == "ALIAS_SELECTION_REQUIRED"
