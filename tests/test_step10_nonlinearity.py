"""Unit tests for the Step 10 bridge-star detector-linearity diagnostic.

Builds a synthetic per-measurement ``obs`` table at two exposure levels and
checks that ``_write_nonlinearity_diag`` reports a flat delta for a linear
detector, a sloped delta when a brightness-dependent offset is injected, and
skips single-exposure data. Constructs a bare worker via ``__new__`` so no Qt
event loop is needed; skipped on the no-GUI CI where PyQt5 isn't installed.
"""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("PyQt5")

from apex.gui.workflow.cmd.step10_zeropoint_calibration import (
    ZeropointCalibrationWorker,
    build_zp_qc_summary,
    build_cmd_qc_summary,
    build_gaia_cmd_comparison,
    build_gaia_cmd_snr_sweep,
    export_cmd_qc_products,
    export_gaia_cmd_comparison_products,
    export_zp_qc_products,
    select_cmd_qc_axes,
)


def _worker():
    w = ZeropointCalibrationWorker.__new__(ZeropointCalibrationWorker)
    w._log = lambda *a, **k: None
    return w


def _obs(stars_mag, exptimes, delta_fn):
    rows = []
    for sid, mtrue in enumerate(stars_mag, start=1):
        for exp in exptimes:
            for _ in range(3):  # a few measurements per (star, level)
                rows.append({
                    "ID": sid, "FILTER": "V",
                    "mag_cal": mtrue + delta_fn(mtrue, exp),
                    "exptime": float(exp), "snr": 50.0,
                })
    return pd.DataFrame(rows)


STARS = [12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0]


def test_linear_detector_flat_delta(tmp_path):
    obs = _obs(STARS, [60, 240], lambda m, e: 0.0)
    _worker()._write_nonlinearity_diag(obs, tmp_path)

    summ = pd.read_csv(tmp_path / "nonlinearity_summary.csv")
    assert len(summ) == 1
    row = summ.iloc[0]
    assert abs(row["slope_mag_per_mag"]) < 0.02
    assert not bool(row["nonlinearity_flag"])
    assert int(row["n_bridge"]) == len(STARS)
    assert (tmp_path / "nonlinearity_by_exposure.csv").exists()
    assert (tmp_path / "nonlinearity_check.png").exists()


def test_nonlinear_detector_sloped_delta(tmp_path):
    # Long (240s) exposure carries a brightness-dependent offset: a star of mag m
    # reads 0.1*(m-15) fainter. The diagnostic plots delta vs mag_l (which itself
    # includes the offset), so the recovered slope is -0.1/1.1 = -1/11.
    obs = _obs(STARS, [60, 240], lambda m, e: (0.1 * (m - 15.0) if e == 240 else 0.0))
    _worker()._write_nonlinearity_diag(obs, tmp_path)

    summ = pd.read_csv(tmp_path / "nonlinearity_summary.csv")
    row = summ.iloc[0]
    assert row["slope_mag_per_mag"] == pytest.approx(-1.0 / 11.0, abs=1e-6)
    assert bool(row["nonlinearity_flag"])


def test_low_snr_faint_bias_does_not_trigger_nonlinearity(tmp_path):
    obs = _obs(
        STARS,
        [60, 240],
        lambda m, e: (0.25 * (m - 16.5) if e == 60 and m > 16.5 else 0.0),
    )
    faint_short = (obs["exptime"] == 60.0) & (obs["mag_cal"] > 16.5)
    obs.loc[faint_short, "snr"] = 4.0

    _worker()._write_nonlinearity_diag(obs, tmp_path)

    row = pd.read_csv(tmp_path / "nonlinearity_summary.csv").iloc[0]
    assert row["snr_min"] == pytest.approx(20.0)
    assert abs(row["slope_mag_per_mag"]) < 0.02
    assert not bool(row["nonlinearity_flag"])


def test_single_exposure_skips(tmp_path):
    obs = _obs(STARS, [60], lambda m, e: 0.0)
    _worker()._write_nonlinearity_diag(obs, tmp_path)
    # Returns before writing any output for a single exposure level.
    assert not (tmp_path / "nonlinearity_summary.csv").exists()
    assert not (tmp_path / "nonlinearity_by_exposure.csv").exists()


def test_build_zp_qc_summary_combines_global_and_frame_stats():
    coeff = pd.DataFrame(
        {
            "filter": ["g"],
            "zp": [25.1],
            "ct": [0.02],
            "N": [42],
            "scatter_rms": [0.031],
            "color_col": ["g_r"],
            "ref_source": ["Gaia synthetic"],
        }
    )
    frame = pd.DataFrame(
        {
            "file": ["a.fits", "b.fits", "c.fits"],
            "filter": ["g", "g", "g"],
            "zp_frame": [25.0, 25.1, 25.2],
            "zp_scatter": [0.04, 0.03, 0.05],
            "n_ref": [10, 12, 14],
            "outlier_fraction": [0.0, 0.1, 0.0],
            "snr_med": [80.0, 70.0, 90.0],
        }
    )
    cuts = pd.DataFrame(
        {
            "filter": ["g"],
            "n_total": [100],
            "n_kept": [80],
        }
    )
    rejects = pd.DataFrame({"file": ["bad.fits"], "filter": ["g"], "reason": ["n_ref_below_min"]})

    summary = build_zp_qc_summary(coeff, frame, cuts, rejects).set_index("filter")

    assert summary.loc["g", "global_zp"] == pytest.approx(25.1)
    assert summary.loc["g", "color_term"] == pytest.approx(0.02)
    assert summary.loc["g", "n_fit_calibrators"] == 42
    assert summary.loc["g", "n_frame_zp"] == 3
    assert summary.loc["g", "frame_zp_median"] == pytest.approx(25.1)
    assert summary.loc["g", "median_n_ref_per_frame"] == pytest.approx(12)
    assert summary.loc["g", "n_rejected_frames"] == 1
    assert summary.loc["g", "kept_measurement_fraction"] == pytest.approx(0.8)


def test_export_zp_qc_products_writes_summary_and_overview(tmp_path):
    pd.DataFrame(
        {
            "filter": ["g"],
            "zp": [25.1],
            "ct": [0.02],
            "N": [42],
            "scatter_rms": [0.031],
            "color_col": ["g_r"],
            "ref_source": ["Gaia synthetic"],
        }
    ).to_csv(tmp_path / "zp_fit_coefficients.csv", index=False)
    pd.DataFrame(
        {
            "file": ["a.fits", "b.fits", "c.fits"],
            "filter": ["g", "g", "g"],
            "zp_frame": [25.0, 25.1, 25.2],
            "zp_scatter": [0.04, 0.03, 0.05],
            "n_ref": [10, 12, 14],
            "outlier_fraction": [0.0, 0.1, 0.0],
            "snr_med": [80.0, 70.0, 90.0],
        }
    ).to_csv(tmp_path / "frame_zeropoint.csv", index=False)
    pd.DataFrame({"filter": ["g"], "n_total": [100], "n_kept": [80]}).to_csv(
        tmp_path / "frame_zeropoint_cut_summary.csv", index=False
    )

    saved = export_zp_qc_products(tmp_path)

    assert tmp_path / "zp_qc_summary.csv" in saved
    assert tmp_path / "step10_zp_qc_overview.png" in saved
    assert (tmp_path / "zp_qc_summary.csv").exists()
    assert (tmp_path / "step10_zp_qc_overview.png").exists()


def test_build_cmd_qc_summary_selects_best_calibrated_cmd_axes():
    df = pd.DataFrame(
        {
            "ID": range(1, 7),
            "mag_std_g": [15.0, 15.4, 16.0, 16.5, 17.2, np.nan],
            "mag_std_r": [14.5, 14.8, 15.2, 15.9, 16.6, 17.0],
            "mag_std_i": [14.2, 14.5, 14.9, np.nan, 16.1, 16.4],
            "mag_std_err_g": [0.01, 0.01, 0.02, 0.02, 0.03, np.nan],
            "mag_std_err_r": [0.01, 0.01, 0.02, 0.02, 0.03, 0.04],
            "mag_std_err_i": [0.01, 0.02, 0.02, np.nan, 0.04, 0.05],
            "snr_g": [100, 90, 80, 60, 40, np.nan],
            "snr_r": [110, 100, 90, 70, 45, 30],
            "snr_i": [120, 90, 70, np.nan, 35, 25],
        }
    )

    axes = select_cmd_qc_axes(df)
    summary = build_cmd_qc_summary(df).set_index("filter")

    assert axes["color_a"] == "g"
    assert axes["color_b"] == "r"
    assert axes["y_band"] == "r"
    assert axes["n"] == 5
    assert summary.loc["r", "n_finite_mag"] == 6
    assert summary.loc["r", "median_mag_err"] == pytest.approx(0.02)
    assert summary.loc["g", "cmd_n_points"] == 5


def test_export_cmd_qc_products_writes_summary_and_overview(tmp_path):
    pd.DataFrame(
        {
            "ID": range(1, 16),
            "mag_std_g": np.linspace(15.0, 20.0, 15),
            "mag_std_r": np.linspace(14.5, 19.2, 15),
            "mag_std_i": np.linspace(14.1, 18.7, 15),
            "mag_std_err_g": np.linspace(0.01, 0.08, 15),
            "mag_std_err_r": np.linspace(0.01, 0.06, 15),
            "mag_std_err_i": np.linspace(0.01, 0.05, 15),
            "snr_g": np.linspace(120.0, 25.0, 15),
            "snr_r": np.linspace(130.0, 30.0, 15),
            "snr_i": np.linspace(140.0, 35.0, 15),
        }
    ).to_csv(tmp_path / "median_by_ID_filter_wide_cmd.csv", index=False)

    saved = export_cmd_qc_products(tmp_path)

    assert tmp_path / "cmd_qc_summary.csv" in saved
    assert tmp_path / "step10_cmd_qc_overview.png" in saved
    assert (tmp_path / "cmd_qc_summary.csv").exists()
    assert (tmp_path / "step10_cmd_qc_overview.png").exists()


def test_build_gaia_cmd_comparison_reports_matched_residuals():
    df = pd.DataFrame(
        {
            "ID": [1, 2, 3, 4],
            "gaia_G": [15.0, 16.0, 17.0, 18.0],
            "gaia_BP": [15.6, 16.7, 17.8, np.nan],
            "gaia_RP": [14.7, 15.6, 16.7, 17.6],
            "gaia_G_syn": [15.1, 15.9, 17.2, 18.1],
            "gaia_BP_RP_syn": [0.95, 1.00, 1.05, 0.50],
        }
    )

    summary = build_gaia_cmd_comparison(df).iloc[0]

    assert int(summary["n_gaia_cmd"]) == 3
    assert int(summary["n_apex_synthetic_gaia_cmd"]) == 4
    assert int(summary["n_matched_cmd"]) == 3
    assert summary["median_delta_G"] == pytest.approx(0.1)
    assert summary["median_delta_BP_RP"] == pytest.approx(-0.05)


def test_build_gaia_cmd_comparison_falls_back_to_native_standard_axes():
    df = pd.DataFrame(
        {
            "ID": [1, 2, 3, 4],
            "gaia_G": [15.0, 16.0, 17.0, 18.0],
            "gaia_BP": [15.6, 16.7, 17.8, 18.8],
            "gaia_RP": [14.7, 15.6, 16.7, 17.7],
            "mag_std_B": [15.5, 16.5, 17.5, 18.5],
            "mag_std_V": [14.9, 15.9, 16.9, 17.9],
            "mag_std_R": [14.5, 15.5, 16.5, 17.5],
        }
    )

    summary = build_gaia_cmd_comparison(df).iloc[0]

    assert summary["mode"] == "native_standard"
    assert summary["cmd_color"] == "B - V"
    assert summary["cmd_y"] == "V"
    assert int(summary["n_matched_cmd"]) == 4
    assert np.isfinite(summary["median_delta_mag"])


def test_build_gaia_cmd_snr_sweep_tracks_cut_sensitivity():
    n = 8
    gaia_g = np.linspace(15.0, 18.0, n)
    gaia_color = np.linspace(0.5, 1.2, n)
    df = pd.DataFrame(
        {
            "ID": range(1, n + 1),
            "gaia_G": gaia_g,
            "gaia_BP": gaia_g + 0.5 * gaia_color,
            "gaia_RP": gaia_g - 0.5 * gaia_color,
            "mag_std_B": gaia_g + gaia_color + np.linspace(0.00, -0.12, n),
            "mag_std_V": gaia_g + np.linspace(0.00, 0.02, n),
            "snr_B": [4, 6, 12, 18, 25, 55, 80, 120],
            "snr_V": [8, 12, 15, 30, 45, 60, 90, 130],
        }
    )

    sweep = build_gaia_cmd_snr_sweep(df, snr_cuts=(5, 10, 20, 50, 100))

    assert list(sweep["snr_cut"]) == [5, 10, 20, 50, 100]
    assert list(sweep["n_used"]) == [7, 6, 4, 3, 1]
    assert sweep.loc[sweep["snr_cut"] == 20, "snr_columns"].iloc[0] == "snr_B,snr_V"
    assert np.isfinite(sweep.loc[sweep["snr_cut"] == 20, "median_delta_color"].iloc[0])


def test_export_gaia_cmd_comparison_products_writes_summary_and_overview(tmp_path):
    n = 20
    gaia_g = np.linspace(14.0, 19.0, n)
    gaia_color = np.linspace(0.4, 1.6, n)
    pd.DataFrame(
        {
            "ID": range(1, n + 1),
            "gaia_G": gaia_g,
            "gaia_BP": gaia_g + 0.5 * gaia_color,
            "gaia_RP": gaia_g - 0.5 * gaia_color,
            "gaia_G_syn": gaia_g + 0.03,
            "gaia_BP_RP_syn": gaia_color - 0.02,
            "mag_std_g": gaia_g + 0.5,
            "mag_std_i": gaia_g - 0.5,
            "snr_g": np.linspace(120.0, 15.0, n),
            "snr_i": np.linspace(110.0, 12.0, n),
        }
    ).to_csv(tmp_path / "median_by_ID_filter_wide_cmd.csv", index=False)

    saved = export_gaia_cmd_comparison_products(tmp_path)

    assert tmp_path / "gaia_cmd_comparison_summary.csv" in saved
    assert tmp_path / "step10_gaia_cmd_comparison.png" in saved
    assert tmp_path / "gaia_cmd_snr_sweep.csv" in saved
    assert tmp_path / "step10_gaia_cmd_snr_sweep.png" in saved
    assert (tmp_path / "gaia_cmd_comparison_summary.csv").exists()
    assert (tmp_path / "step10_gaia_cmd_comparison.png").exists()
    assert (tmp_path / "gaia_cmd_snr_sweep.csv").exists()
    assert (tmp_path / "step10_gaia_cmd_snr_sweep.png").exists()
