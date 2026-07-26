from __future__ import annotations

import numpy as np
import pandas as pd

from apex.analysis.light_curve.comparison_stability_service import (
    ComparisonSelectionConfig,
    build_target_difference,
    compute_leave_one_out_residuals,
    recommend_check_candidate,
    select_adaptive_ensemble,
    select_stable_comparisons,
)
from apex.gui.workflow.lc.step8_target_selection import (
    _build_comparison_preview_payload,
)


def _synthetic_measurements() -> pd.DataFrame:
    rng = np.random.default_rng(4)
    rows = []
    for frame_index in range(100):
        common = 0.10 * np.sin(2 * np.pi * frame_index / 37) + 0.03 * frame_index / 100
        night_id = frame_index // 25 + 1
        airmass = 1.0 + 0.8 * abs((frame_index % 25) - 12) / 12
        for star_id in range(1, 8):
            if star_id == 7 and frame_index % 3:
                continue
            magnitude = 12.0 + 0.25 * star_id + common + rng.normal(0.0, 0.004)
            if star_id == 6:
                magnitude += 0.045 * np.sin(2 * np.pi * frame_index / 22)
            rows.append(
                {
                    "frame": f"frame_{frame_index:03d}",
                    "time": frame_index / 24.0,
                    "night_id": night_id,
                    "star_id": star_id,
                    "mag": magnitude,
                    "mag_err": 0.004,
                    "airmass": airmass,
                    "fwhm": 2.0 + 0.2 * np.sin(frame_index),
                    "sky": 100.0,
                }
            )
    return pd.DataFrame(rows)


def test_leave_one_out_removes_common_transparency_but_retains_variable_star():
    measurements = _synthetic_measurements()
    residuals = compute_leave_one_out_residuals(measurements, range(1, 7))

    stable_sigma = residuals[residuals["star_id"] == 1]["residual"].std()
    variable_sigma = residuals[residuals["star_id"] == 6]["residual"].std()

    assert stable_sigma < 0.01
    assert variable_sigma > 0.025


def test_auto_selection_rejects_variable_and_low_coverage_candidates():
    result = select_stable_comparisons(
        _synthetic_measurements(),
        range(1, 8),
        target_mag=13.0,
        config=ComparisonSelectionConfig(target_count=5),
    )

    assert result["selected_ids"] == [1, 2, 3, 4, 5]
    assert set(result["removed_ids"]) == {6, 7}
    report = result["metrics"].set_index("star_id")
    assert "time_correlated" in report.loc[6, "reasons"]
    assert "low_coverage" in report.loc[7, "reasons"]


def test_target_difference_is_centered_and_keeps_time_metadata():
    measurements = _synthetic_measurements()
    output = build_target_difference(measurements, target_id=1, comparison_id=2)

    assert len(output) == 100
    assert {"frame", "time", "night_id", "value"} <= set(output.columns)
    assert abs(float(np.nanmedian(output["value"]))) < 1e-12


def test_comparison_preview_uses_only_requested_references_and_target():
    measurements = _synthetic_measurements()

    payload = _build_comparison_preview_payload(
        measurements,
        {"source": "aperture"},
        source_id=2,
        display_id=22,
        filter_name="g",
        target_id=1,
        reference_ids=[2, 3, 4, 3],
    )

    assert len(payload["candidate"]) == 100
    assert len(payload["loo"]) == 100
    assert len(payload["target_delta"]) == 100
    assert payload["default_mode"] == "loo"
    assert payload["metrics"]["n"] == 100
    assert payload["can_assign"] is True


def test_target_preview_defaults_to_raw_magnitude_and_disables_roles():
    payload = _build_comparison_preview_payload(
        _synthetic_measurements(),
        {"source": "psf"},
        source_id=1,
        display_id=11,
        filter_name="g",
        target_id=1,
        reference_ids=[2, 3, 4],
    )

    assert payload["default_mode"] == "mag"
    assert payload["loo"].empty
    assert payload["target_delta"].empty
    assert payload["metrics"]["status"] == "target"
    assert payload["can_assign"] is False


def test_rejected_final_candidates_are_not_used_to_fill_the_requested_count():
    measurements = _synthetic_measurements()
    result = select_stable_comparisons(
        measurements,
        [1, 2, 3],
        config=ComparisonSelectionConfig(
            min_coverage=1.01,
            min_comparisons=3,
            target_count=3,
        ),
    )

    assert result["selected_ids"] == []
    assert set(result["metrics"]["status"]) == {"reject"}


def test_large_gaia_source_ids_survive_the_stability_pipeline():
    measurements = _synthetic_measurements()
    base = 1_387_320_379_874_985_000
    source_ids = {small: base + small * 128 for small in range(1, 8)}
    measurements["star_id"] = measurements["star_id"].map(source_ids).astype("int64")

    result = select_stable_comparisons(
        measurements,
        source_ids.values(),
        target_mag=13.0,
        config=ComparisonSelectionConfig(target_count=5),
    )
    difference = build_target_difference(
        measurements,
        target_id=source_ids[1],
        comparison_id=source_ids[2],
    )

    assert result["selected_ids"] == [source_ids[index] for index in range(1, 6)]
    assert set(result["removed_ids"]) == {source_ids[6], source_ids[7]}
    assert result["metrics"]["star_id"].astype("int64").min() >= base
    assert len(difference) == 100


def test_check_candidate_is_stable_and_independent_of_selected_ensemble():
    metrics = pd.DataFrame(
        {
            "star_id": [101, 102, 103, 104, 105],
            "status": ["stable", "stable", "stable", "suspect", "stable"],
            "coverage": [1.0, 1.0, 0.85, 1.0, 1.0],
            "selection_score": [0.10, 0.20, 0.05, 0.01, 0.30],
            "robust_sigma": [0.01, 0.02, 0.01, 0.01, 0.03],
        }
    )

    candidate = recommend_check_candidate(metrics, excluded_ids={101})

    assert candidate is not None
    assert int(candidate["star_id"]) == 102


def test_check_candidate_returns_none_without_an_independent_stable_source():
    metrics = pd.DataFrame(
        {
            "star_id": [101, 102],
            "status": ["stable", "reject"],
            "coverage": [1.0, 1.0],
            "selection_score": [0.10, 0.20],
        }
    )

    assert recommend_check_candidate(metrics, excluded_ids={101}) is None


def test_adaptive_ensemble_reserves_one_stable_check_star():
    measurements = _synthetic_measurements()
    stability = select_stable_comparisons(
        measurements,
        range(1, 8),
        target_mag=13.0,
        config=ComparisonSelectionConfig(target_count=20),
    )

    result = select_adaptive_ensemble(
        measurements, stability["metrics"], max_comparisons=8
    )

    assert len(result["selected_ids"]) >= 3
    assert result["check_id"] is not None
    assert result["check_id"] not in result["selected_ids"]
    assert set(result["selected_ids"] + [result["check_id"]]) <= {1, 2, 3, 4, 5}
    assert not result["ensemble_trials"].empty


def test_adaptive_ensemble_requires_room_for_independent_check():
    metrics = pd.DataFrame(
        {
            "star_id": [1, 2, 3],
            "status": ["stable", "stable", "stable"],
            "selection_score": [0.1, 0.2, 0.3],
        }
    )

    result = select_adaptive_ensemble(_synthetic_measurements(), metrics)

    assert result["selected_ids"] == []
    assert result["check_id"] is None
    assert "plus one check star" in result["reason"]
