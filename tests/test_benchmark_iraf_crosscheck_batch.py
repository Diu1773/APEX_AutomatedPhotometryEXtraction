from pathlib import Path

import pandas as pd
import pytest

from apex.benchmark.iraf_crosscheck_batch import (
    recommend_thresholds,
    select_iraf_batch_frames,
    summarize_apex_reference_repeatability,
    summarize_batch_filters,
    summarize_batch_thresholds,
    summarize_paired_repeatability,
)


def _write_project(root: Path) -> None:
    step7_dir = root / "result" / "step7_forced_phot"
    zp_dir = root / "result" / "cmd_zeropoint"
    step7_dir.mkdir(parents=True)
    zp_dir.mkdir(parents=True)
    index_rows = []
    zp_rows = []
    for filter_name in ("g", "r"):
        for number in range(2):
            filename = f"frame_{filter_name}_{number}.fit"
            (root / filename).touch()
            step7_path = step7_dir / f"photometry_{filename}.tsv"
            step7_path.write_text("x\ty\tmag_inst\n1\t1\t-10\n", encoding="utf-8")
            index_rows.append(
                {
                    "file": filename,
                    "filter": filter_name,
                    "status": "ok",
                    "wcs_ok": True,
                    "fwhm_px": 5.0 + number,
                    "path": str(step7_path),
                }
            )
            zp_rows.append({"file": filename, "filter": filter_name, "zp_frame": 22.0 + number})
    index_rows.append(
        {
            "file": "bad.fit",
            "filter": "g",
            "status": "failed",
            "wcs_ok": False,
            "fwhm_px": 5.0,
            "path": "",
        }
    )
    pd.DataFrame(index_rows).to_csv(step7_dir / "photometry_index.csv", index=False)
    pd.DataFrame(zp_rows).to_csv(zp_dir / "frame_zeropoint.csv", index=False)


def test_select_iraf_batch_frames_uses_ok_wcs_calibrated_frames(tmp_path):
    _write_project(tmp_path)

    selected = select_iraf_batch_frames(tmp_path, filters=["g", "r"])

    assert len(selected) == 4
    assert set(selected["filter_key"]) == {"g", "r"}
    assert all(Path(path).is_file() for path in selected["input_fits"])
    assert all(Path(path).is_file() for path in selected["step7_tsv"])


def test_summarize_batch_thresholds_and_recommendations():
    summary = pd.DataFrame(
        {
            "filter": ["g", "g", "g", "g"],
            "file": ["a", "b", "a", "b"],
            "threshold": [12.0, 12.0, 7.0, 7.0],
            "status": ["ok", "ok", "ok", "ok"],
            "n_detected": [100, 110, 200, 220],
            "n_matched": [90, 95, 96, 98],
            "n_apex_reference": [100, 100, 100, 100],
            "recall_vs_apex_reference": [0.90, 0.95, 0.96, 0.98],
            "matched_fraction_of_iraf": [0.90, 0.86, 0.48, 0.45],
            "mad_sigma_delta_mag_zp_aligned": [0.01, 0.02, 0.03, 0.04],
            "p95_abs_delta_mag_zp_aligned": [0.02, 0.03, 0.05, 0.06],
            "median_apex_mag_err": [0.03, 0.03, 0.03, 0.03],
            "median_iraf_merr": [0.02, 0.02, 0.02, 0.02],
            "apex_n_master": [120, 120, 120, 120],
            "apex_n_detected": [100, 105, 100, 105],
            "apex_n_forced": [20, 15, 20, 15],
            "apex_n_valid_phot": [120, 120, 120, 120],
            "apex_detected_rate": [0.83, 0.88, 0.83, 0.88],
            "apex_forced_rate": [0.17, 0.12, 0.17, 0.12],
        }
    )

    by_threshold = summarize_batch_thresholds(summary)
    by_filter = summarize_batch_filters(summary)
    recommendations = recommend_thresholds(by_threshold)

    assert len(by_threshold) == 2
    assert len(by_filter) == 1
    assert by_filter.loc[0, "apex_n_valid_phot_median"] == pytest.approx(120.0)
    assert by_filter.loc[0, "apex_n_forced_median"] == pytest.approx(17.5)
    assert by_filter.loc[0, "detection_qc"] == "likely over-detection; raise threshold"
    assert by_threshold.loc[by_threshold["threshold"] == 12.0, "n_frames"].iloc[0] == 2
    assert recommendations.loc[0, "threshold"] == pytest.approx(7.0)
    assert recommendations.loc[0, "reason"].startswith("highest purity")


def test_single_tested_threshold_is_not_claimed_as_optimized():
    by_threshold = pd.DataFrame(
        {
            "filter": ["i"],
            "threshold": [9.0],
            "recall_median": [0.99],
            "matched_fraction_median": [0.72],
        }
    )

    recommendations = recommend_thresholds(by_threshold)

    assert recommendations.loc[0, "threshold"] == pytest.approx(9.0)
    assert "not an independently optimized" in recommendations.loc[0, "reason"]


def test_repeatability_summaries_compare_apex_and_iraf_scatter():
    paired = pd.DataFrame(
        {
            "filter": ["g"] * 6,
            "threshold": [9.0] * 6,
            "master_id": [1, 1, 1, 2, 2, 2],
            "apex_mag_cal_frame": [12.00, 12.01, 11.99, 13.00, 13.02, 12.98],
            "iraf_mag": [14.00, 14.04, 13.96, 15.00, 15.08, 14.92],
            "iraf_mag_cal_apcorr_zp": [12.00, 12.015, 11.985, 13.00, 13.03, 12.97],
            "iraf_mag_aligned_frame": [12.00, 12.02, 11.98, 13.00, 13.04, 12.96],
            "mag_err": [0.01] * 6,
            "iraf_merr": [0.02] * 6,
        }
    )
    by_star, summary = summarize_paired_repeatability(paired, min_frames=3)

    assert len(by_star) == 2
    assert summary.loc[0, "n_paired_stars"] == 2
    assert summary.loc[0, "n_iraf_calibrated_frames_median"] == pytest.approx(3.0)
    assert summary.loc[0, "iraf_raw_repeat_mad_median"] > summary.loc[0, "iraf_repeat_mad_median"]
    assert summary.loc[0, "iraf_calibrated_repeat_mad_median"] > summary.loc[0, "apex_repeat_mad_median"]
    assert summary.loc[0, "iraf_calibrated_repeat_mad_median"] < summary.loc[0, "iraf_raw_repeat_mad_median"]
    assert summary.loc[0, "iraf_repeat_mad_median"] > summary.loc[0, "apex_repeat_mad_median"]
    assert summary.loc[0, "iraf_over_apex_repeat_mad_median"] == pytest.approx(2.0)


def test_apex_reference_repeatability_counts_high_snr_coverage():
    observations = pd.DataFrame(
        {
            "filter": ["g"] * 6,
            "master_id": [1, 1, 1, 2, 2, 2],
            "apex_mag_cal_frame": [12.00, 12.01, 11.99, 13.00, 13.02, 12.98],
        }
    )

    by_star, summary = summarize_apex_reference_repeatability(observations, min_frames=3)

    assert len(by_star) == 2
    assert summary.loc[0, "n_apex_reference_stars"] == 2
    assert summary.loc[0, "n_apex_frames_median"] == pytest.approx(3.0)
