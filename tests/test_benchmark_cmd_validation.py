import pandas as pd
import pytest

from apex.benchmark.cmd_validation import (
    compute_zeropoint_residuals,
    mark_zeropoint_inliers,
    summarize_repeatability,
    summarize_zeropoint_residuals,
)


def test_compute_zeropoint_residuals_uses_color_term():
    calibrators = pd.DataFrame(
        {
            "ID": [1, 2, 3],
            "delta_g": [22.1, 22.2, 22.3],
            "color_g_r": [0.0, 1.0, 2.0],
            "ref_g": [14.0, 15.0, 16.0],
            "snr_g": [100.0, 90.0, 80.0],
        }
    )
    coefficients = pd.DataFrame(
        [{"filter": "g", "zp": 22.1, "ct": 0.1, "color_col": "g_r"}]
    )

    residuals = compute_zeropoint_residuals(calibrators, coefficients)

    assert residuals["residual_mag"].tolist() == pytest.approx([0.0, 0.0, 0.0])


def test_summarize_zeropoint_residuals_reports_rms_and_mad():
    residuals = pd.DataFrame(
        {
            "filter": ["g", "g", "g"],
            "residual_mag": [-0.1, 0.0, 0.1],
            "zp": [22.0, 22.0, 22.0],
            "ct": [0.0, 0.0, 0.0],
            "color_col": ["none", "none", "none"],
            "fit_scatter_rms": [0.05, 0.05, 0.05],
        }
    )

    summary = summarize_zeropoint_residuals(mark_zeropoint_inliers(residuals))

    assert summary.loc[0, "filter"] == "g"
    assert summary.loc[0, "n_raw"] == 3
    assert summary.loc[0, "n_clipped"] == 3
    assert summary.loc[0, "residual_rms_clipped"] == pytest.approx((0.02 / 3) ** 0.5)
    assert summary.loc[0, "residual_mad_sigma_clipped"] == pytest.approx(0.14826)
    assert summary.loc[0, "fit_scatter_rms"] == pytest.approx(0.05)


def test_summarize_repeatability_bins_by_filter_and_magnitude():
    observations = pd.DataFrame(
        {
            "filter": ["g", "g", "g", "g", "g", "g"],
            "ID": [1, 1, 1, 2, 2, 2],
            "mag_cal_frame": [14.0, 14.1, 13.9, 15.0, 15.2, 14.8],
            "mag_err": [0.02, 0.02, 0.02, 0.04, 0.04, 0.04],
            "snr": [50, 50, 50, 25, 25, 25],
        }
    )

    stars, bins = summarize_repeatability(observations, min_frames=3, bin_width=0.5)

    assert len(stars) == 2
    assert set(bins["filter"]) == {"g"}
    assert bins["n_stars"].sum() == 2
