from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from apex.analysis.cmd.isochrone_fitter_v2 import FitResult, IsochroneFitterV2
from apex.gui.workflow.cmd.step10_zeropoint_calibration import CmdViewerWindow
from apex.gui.workflow.cmd.step12_isochrone_model import (
    _bands_from_df,
    _extinction_magnitude_shift,
    _legacy_b_calibration_reason,
    _preferred_err_col,
    _preferred_mag_col,
    _segmented_isochrone_line,
)


def _sample_iso_data() -> np.ndarray:
    rows = []
    for age in (8.0, 8.5):
        for mh in (-0.1, 0.0):
            for idx in range(12):
                mag = 20.0 - idx * 0.1
                row = np.zeros(6, dtype=float)
                row[1] = mh
                row[2] = age
                row[3] = mag
                row[4] = mag - 0.4
                row[5] = 0.2 + idx * 0.01
                rows.append(row)

    data = np.asarray(rows, dtype=float)
    return data[::-1]


def test_isochrone_fitter_accepts_preloaded_iso_data_without_reading_file(tmp_path):
    missing_file = tmp_path / "missing_isochrone.dat"

    fitter = IsochroneFitterV2(
        missing_file,
        col_g=3,
        col_r=4,
        col_mag=3,
        col_mass=5,
        iso_data=_sample_iso_data(),
    )

    assert len(fitter.ages) == 2
    assert len(fitter.metallicities) == 2
    assert len(fitter._iso_cache) == 4


def test_isochrone_cache_groups_once_and_sorts_by_magnitude(tmp_path):
    fitter = IsochroneFitterV2(
        tmp_path / "missing_isochrone.dat",
        col_g=3,
        col_r=4,
        col_mag=3,
        col_mass=5,
        iso_data=_sample_iso_data(),
    )

    for cached in fitter._iso_cache.values():
        assert len(cached) == 12
        assert np.all(np.diff(cached[:, 3]) >= 0)

    color, mag, mass = fitter._interpolate_isochrone(8.1, -0.08, dm=10.0, e_gr=0.0)
    assert len(color) == 12
    assert len(mag) == 12
    assert len(mass) == 12


def test_fit_result_summary_uses_active_cmd_axis_labels():
    result = FitResult(
        log_age=8.5,
        metallicity=0.0,
        distance_mod=10.0,
        extinction_gr=0.32,
        extinction_gr_err=0.01,
        chi2=12.0,
        reduced_chi2=1.2,
        n_stars=50,
        extinction_label="E(B-V)",
        color_label="B-V",
        mag_label="V",
    )

    summary = result.summary()

    assert "E(B-V)" in summary
    assert "E(g-r)" not in summary
    assert "CMD axes : color=B-V, mag=V" in summary


def test_cmd_magnitude_conventions_prefer_mag_cal_columns():
    df = pd.DataFrame(
        {
            "mag_inst_g": [-10.0],
            "mag_inst_r": [-10.2],
            "mag_std_g": [14.9],
            "mag_std_r": [14.7],
            "mag_cal_g": [15.0],
            "mag_cal_r": [14.8],
            "mag_cal_err_g": [0.01],
        }
    )

    pairs, bands = _bands_from_df(df)

    assert bands == ["g", "r"]
    assert pairs == [("g", "r")]
    assert _preferred_mag_col(df.columns, "g") == ("mag_cal_g", "Calibrated")
    assert _preferred_err_col(df.columns, "g") == "mag_cal_err_g"


def test_cmd_viewer_aliases_mag_cal_to_legacy_mag_std():
    df = pd.DataFrame({"mag_cal_g": [15.0], "mag_cal_err_g": [0.01]})

    out = CmdViewerWindow._with_calibrated_aliases(df)

    assert out["mag_std_g"].tolist() == [15.0]
    assert out["mag_std_err_g"].tolist() == [0.01]


def test_bv_color_excess_also_dims_b_magnitude():
    shift = _extinction_magnitude_shift(0.1, "B", "V", "B")

    assert shift == pytest.approx(0.4390642)


def test_segmented_isochrone_line_breaks_anomalous_model_jump():
    color = np.array([0.0, 0.1, 0.2, 0.3, 2.0, 2.1])
    magnitude = np.array([5.0, 4.9, 4.8, 4.7, 30.0, 29.9])

    line_color, line_magnitude = _segmented_isochrone_line(color, magnitude)

    break_indices = np.flatnonzero(np.isnan(line_color))
    assert break_indices.tolist() == [4]
    assert np.isnan(line_magnitude[4])
    assert np.allclose(line_color[~np.isnan(line_color)], color)


def test_segmented_isochrone_line_keeps_smooth_sparse_sequence_connected():
    color = np.array([1.5, 1.3, 1.1, 0.9])
    magnitude = np.array([18.0, 17.0, 16.0, 15.0])

    line_color, line_magnitude = _segmented_isochrone_line(color, magnitude)

    assert np.array_equal(line_color, color)
    assert np.array_equal(line_magnitude, magnitude)


def test_legacy_b_calibration_requires_step10_rerun(tmp_path):
    pd.DataFrame(
        [{"filter": "B", "zp": -3.5, "ct": -0.2}]
    ).to_csv(tmp_path / "zp_fit_coefficients.csv", index=False)

    reason = _legacy_b_calibration_reason(tmp_path, ("B", "V"))

    assert reason is not None
    assert "Rerun Step 10" in reason


def test_pancino_b_calibration_is_accepted(tmp_path):
    pd.DataFrame(
        [
            {
                "filter": "B",
                "zp": -3.5,
                "ct": -0.2,
                "ref_source": "Pancino+2022 dwarf",
            }
        ]
    ).to_csv(tmp_path / "zp_fit_coefficients.csv", index=False)

    assert _legacy_b_calibration_reason(tmp_path, ("B", "V")) is None
