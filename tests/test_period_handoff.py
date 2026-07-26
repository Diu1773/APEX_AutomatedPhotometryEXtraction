from __future__ import annotations

from pathlib import Path

import pandas as pd

from apex.gui.tools.variable_star import (
    _detect_target_id_from_df,
    _handoff_mag_columns_match,
    _is_primary_lightcurve_path,
    _route_period_handoff,
    _select_period_handoff_series_key,
)


def _scan_result(period: float = 0.104) -> dict:
    return {
        "raw_ls": {
            "best_period": period,
            "best_power": 0.9,
            "top_periods": [period, 0.059],
            "top_powers": [0.9, 0.1],
        }
    }


def test_handoff_selects_exact_source_target_and_magnitude_series():
    options = {
        "wrong-target": {
            "source": "lightcurve_ID2_global.csv",
            "target_id": 2,
            "mag_col": "diff_mag_corr",
            "corr_tag": "Global ensemble",
        },
        "raw": {
            "source": "lightcurve_ID1_raw.csv",
            "target_id": 1,
            "mag_col": "diff_mag_raw",
            "corr_tag": "Raw",
        },
        "correct": {
            "source": "lightcurve_ID1_global.csv",
            "target_id": 1,
            "mag_col": "diff_mag_corr",
            "corr_tag": "Global ensemble",
        },
    }
    handoff = {
        "source_file": r"C:\workspace\lightcurve_ID1_global.csv",
        "target_id": 1,
        "mag_col": "diff_mag_corr",
        "correction_mode": "global",
    }

    assert _select_period_handoff_series_key(options, handoff) == "correct"


def test_raw_diff_mag_columns_are_handoff_compatible():
    assert _handoff_mag_columns_match("diff_mag", "diff_mag_raw", "raw")
    assert not _handoff_mag_columns_match("diff_mag_corr", "diff_mag_raw", "global")


def test_ambiguous_alias_handoff_stays_at_scan_review():
    mode, workflow, reason = _route_period_handoff(
        _scan_result(),
        {
            "status": "AMBIGUOUS",
            "reason": "Daily aliases remain.",
            "candidates": [{"period": 0.104}, {"period": 0.093}],
        },
        {"status": "SINGLE-COMPATIBLE", "reason": "No residual mode."},
    )

    assert mode == "single"
    assert workflow == "scan"
    assert "AMBIGUOUS" in reason
    assert "2 candidate(s)" in reason


def test_resolved_multimode_handoff_routes_to_multi_workflow():
    mode, workflow, reason = _route_period_handoff(
        _scan_result(),
        {"status": "RESOLVED", "candidates": [{"period": 0.104}]},
        {
            "status": "MULTIMODE-SUSPECT",
            "reason": "Independent residual frequency detected.",
        },
    )

    assert mode == "multi"
    assert workflow == "multi"
    assert reason == "Independent residual frequency detected."


def test_target_id_detection_does_not_treat_check_filename_as_target():
    assert _detect_target_id_from_df(
        pd.DataFrame({"check_id": [14]}),
        "lightcurve_check_g_ID14_raw.csv",
    ) is None
    assert not _is_primary_lightcurve_path(Path("lightcurve_check_g_ID14_raw.csv"))
    assert _is_primary_lightcurve_path(Path("lightcurve_ID1_raw.csv"))
