from pathlib import Path

import pandas as pd
import pytest

from apex.benchmark.cmd_combined_validation import (
    ValidationInput,
    build_dataset_summary,
    build_filter_summary,
    collect_combined_tables,
    parse_validation_input,
)


def _write_validation(root: Path, *, m50_loss: float, false_count: int, repeatability: float, zp: float):
    root.mkdir(parents=True)
    pd.DataFrame(
        [
            {
                "filter": "g",
                "m50_best": 15.0,
                "m50_worst": 15.0 - m50_loss,
                "m50_loss_best_minus_worst": m50_loss,
            }
        ]
    ).to_csv(root / "ast_m50_loss.csv", index=False)
    pd.DataFrame(
        [
            {
                "filter": "g",
                "condition": "best",
                "file": "a.fit",
                "fwhm_px": 5.0,
                "n_trials": 1,
                "n_injected": 100,
                "new_detections": 10,
                "new_false_detections": false_count,
                "false_per_trial": false_count,
                "false_per_1000_injected": false_count * 10.0,
                "false_fraction_of_new": false_count / 10.0,
                "run_dir": str(root / "run"),
            }
        ]
    ).to_csv(root / "false_positive_summary.csv", index=False)
    pd.DataFrame(
        [
            {
                "filter": "g",
                "mag_bin": 15.0,
                "n_stars": 3,
                "mag_median": 15.0,
                "repeatability_std_median": repeatability,
                "repeatability_mad_median": repeatability,
                "mag_err_median": 0.01,
                "snr_median": 50.0,
                "n_frames_median": 3,
            }
        ]
    ).to_csv(root / "repeatability_by_magnitude.csv", index=False)
    pd.DataFrame(
        [
            {
                "filter": "g",
                "n_raw": 10,
                "n_clipped": 9,
                "residual_median_raw": 0.0,
                "residual_rms_raw": zp,
                "residual_mad_sigma_raw": zp,
                "residual_p95_abs_raw": zp,
                "residual_median_clipped": 0.0,
                "residual_rms_clipped": zp,
                "residual_mad_sigma_clipped": zp,
                "residual_p95_abs_clipped": zp,
                "zp": 22.0,
                "ct": 0.0,
                "color_col": "none",
                "fit_scatter_rms": zp,
            }
        ]
    ).to_csv(root / "zeropoint_residual_summary.csv", index=False)
    pd.DataFrame(
        [
            {
                "filter": "g",
                "n_frames": 3,
                "zp_frame_median": 22.0,
                "zp_frame_std": 0.01,
                "zp_scatter_median": 0.02,
                "n_ref_median": 10,
                "snr_med_median": 50,
            }
        ]
    ).to_csv(root / "frame_zeropoint_summary.csv", index=False)
    pd.DataFrame(
        [
            {
                "filter": "g",
                "condition": "best",
                "file": "a.fit",
                "input_fits": "a.fit",
                "fwhm_px_step7": 5.0,
                "fwhm_px_benchmark": 5.0,
                "zp_frame": 22.0,
                "stage": "precision",
                "m90": 14.5,
                "m50": 15.0,
                "m50_ci95_low": 14.9,
                "m50_ci95_high": 15.1,
                "m10": 15.5,
                "completeness": 0.5,
                "forced_mag_bias_median": 0.0,
                "forced_mag_scatter_mad": 0.05,
                "run_dir": str(root / "run"),
                "condition_rank": 0,
            }
        ]
    ).to_csv(root / "precision.csv", index=False)


def test_parse_validation_input_supports_role():
    parsed = parse_validation_input("M5:crowded=benchmark/runs/m5")

    assert parsed.label == "M5"
    assert parsed.role == "crowded"
    assert parsed.path == "benchmark/runs/m5"


def test_combined_summaries(tmp_path):
    a = tmp_path / "a"
    b = tmp_path / "b"
    _write_validation(a, m50_loss=1.0, false_count=2, repeatability=0.01, zp=0.02)
    _write_validation(b, m50_loss=2.0, false_count=4, repeatability=0.03, zp=0.06)
    tables = collect_combined_tables(
        [
            ValidationInput(label="A", role="open", path=str(a)),
            ValidationInput(label="B", role="crowded", path=str(b)),
        ]
    )

    dataset_summary = build_dataset_summary(tables)
    filter_summary = build_filter_summary(tables)

    assert dataset_summary.set_index("dataset").loc["A", "false_per_1000_injected"] == pytest.approx(20.0)
    assert dataset_summary.set_index("dataset").loc["B", "m50_loss_median"] == pytest.approx(2.0)
    assert filter_summary.set_index(["dataset", "filter"]).loc[("A", "g"), "m50_loss"] == pytest.approx(1.0)
