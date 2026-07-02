from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("PyQt5")
from apex.gui.workflow.cmd.step8_psf_photometry import (
    _allstar_newton_group,
    _allstar_newton_one,
    _build_psf_frame_qc_table,
    _build_psf_qc_summary,
    _load_detect_positions,
)


def test_step8_detection_loader_prefers_step4_quality_flux(tmp_path):
    result_dir = tmp_path / "result"
    cache_dir = tmp_path / "cache"
    step4_dir = result_dir / "step4_detection"
    cache_dir.mkdir()
    step4_dir.mkdir(parents=True)

    pd.DataFrame(
        {
            "det_uid": [1, 2],
            "x": [float("nan"), 30.0],
            "y": [10.0, 40.0],
            "flux_for_quality": [111.0, 1234.0],
            "dao_flux": [222.0, 999.0],
            "peak_adu": [333.0, 888.0],
            "epsf_candidate": [True, True],
        }
    ).to_csv(step4_dir / "detect_frame.fits.csv", index=False)

    out = _load_detect_positions("frame.fits", cache_dir, result_dir)

    assert out is not None
    assert list(out.index) == [0]
    assert out.loc[0, "det_uid"] == 2
    assert out.loc[0, "flux_init"] == 1234.0


def test_allstar_newton_degenerate_paths_keep_five_value_contract():
    eval_psf = lambda x, y: np.ones_like(x, dtype=float)

    small = _allstar_newton_one(
        np.ones((2, 2)),
        1.0,
        1.0,
        0,
        0,
        10.0,
        eval_psf,
    )
    singular = _allstar_newton_one(
        np.ones((5, 5)),
        2.0,
        2.0,
        0,
        0,
        10.0,
        eval_psf,
    )
    grouped = _allstar_newton_group(
        np.ones((2, 2)),
        [(1.0, 1.0, 10.0), (2.0, 2.0, 8.0)],
        0,
        0,
        eval_psf,
    )

    assert len(small) == 5
    assert len(singular) == 5
    assert all(len(result) == 5 for result in grouped)


def test_build_psf_qc_summary_groups_photometry_and_comparison_stats():
    idx = pd.DataFrame(
        {
            "file": ["a.fits", "b.fits"],
            "filter": ["g", "r"],
            "n": [3, 2],
            "n_goodmag": [2, 1],
            "n_fail": [1, 1],
            "n_new_iter": [0, 1],
        }
    )
    phot = pd.DataFrame(
        {
            "FILTER": ["g", "g", "g", "r", "r"],
            "flags_psf": [0, 0, 1, 0, 2],
            "mag_psf": [15.0, 16.0, 17.0, 18.0, 19.0],
            "mag_psf_err": [0.01, 0.02, 0.5, 0.03, 0.4],
            "snr_psf": [100.0, 50.0, 3.0, 40.0, 2.0],
            "qfit": [1.0, 6.0, 8.0, 2.0, 9.0],
        }
    )
    meta = pd.DataFrame(
        {
            "filter": ["g", "g", "r"],
            "iter": [1, 2, 1],
            "residual_std": [4.0, 3.0, 8.0],
        }
    )
    cmp_df = pd.DataFrame(
        {
            "FILTER": ["g", "g", "r"],
            "mag_ap": [15.1, 15.8, 18.2],
            "mag_psf": [15.0, 16.0, 18.0],
        }
    )

    summary = _build_psf_qc_summary(idx, phot, meta, cmp_df).set_index("filter")

    assert summary.loc["ALL", "n_frames"] == 2
    assert summary.loc["ALL", "n_psf_rows"] == 5
    assert summary.loc["ALL", "n_clean"] == 3
    assert summary.loc["g", "n_clean"] == 2
    assert summary.loc["g", "clean_fraction"] == pytest.approx(2 / 3)
    assert summary.loc["g", "qfit_gt5_fraction"] == pytest.approx(0.5)
    assert summary.loc["g", "median_ap_minus_psf"] == pytest.approx(-0.05)
    assert summary.loc["g", "residual_std_iter1_mean"] == pytest.approx(4.0)
    assert summary.loc["g", "residual_std_iter2_mean"] == pytest.approx(3.0)


def test_build_psf_qc_summary_uses_canonical_filter_keys():
    idx = pd.DataFrame({"file": ["v_frame.fits"], "filter": ["v"], "n": [1]})
    phot = pd.DataFrame({"FILTER": ["V"], "flags_psf": [0], "mag_psf": [14.0]})

    summary = _build_psf_qc_summary(idx, phot).set_index("filter")

    assert "V" in summary.index
    assert "v" not in summary.index
    assert summary.loc["V", "n_frames"] == 1
    assert summary.loc["V", "n_psf_rows"] == 1


def test_build_psf_frame_qc_table_keeps_core_cut_and_residual_stats():
    idx = pd.DataFrame(
        {
            "file": ["a.fits"],
            "filter": ["g"],
            "n": [10],
            "n_goodmag": [8],
            "n_fail": [2],
            "n_new_iter": [3],
        }
    )
    meta = pd.DataFrame(
        {
            "file": ["a.fits", "a.fits"],
            "filter": ["g", "g"],
            "iter": [1, 2],
            "residual_std": [10.0, 7.5],
            "n_fit": [10, 13],
            "n_new_raw": [0, 4],
            "n_new_kept": [0, 3],
            "core_cut_enabled": [True, True],
            "core_cut_x_px": [500.0, 500.0],
            "core_cut_y_px": [510.0, 510.0],
            "core_cut_radius_px": [80.0, 80.0],
            "core_cut_method": ["auto", "auto"],
            "core_cut_reason": ["density_peak", "density_peak"],
            "n_core_excluded_init": [5, 5],
            "n_core_excluded_redetect": [2, 2],
            "n_core_excluded_result": [7, 7],
        }
    )

    table = _build_psf_frame_qc_table(idx, meta).set_index("file")

    assert table.loc["a.fits", "filter"] == "g"
    assert bool(table.loc["a.fits", "core_cut_enabled"]) is True
    assert table.loc["a.fits", "core_cut_radius_px"] == pytest.approx(80.0)
    assert table.loc["a.fits", "n_core_excluded_result"] == 7
    assert table.loc["a.fits", "residual_std_iter1"] == pytest.approx(10.0)
    assert table.loc["a.fits", "residual_std_final"] == pytest.approx(7.5)
    assert table.loc["a.fits", "residual_std_frac_change"] == pytest.approx(-0.25)
    assert table.loc["a.fits", "n_new_kept_iter2"] == 3


def test_build_psf_frame_qc_table_parses_string_false_core_cut():
    idx = pd.DataFrame(
        {
            "file": ["a.fits"],
            "filter": ["g"],
            "core_cut_enabled": ["False"],
            "core_cut_radius_px": [80.0],
        }
    )

    table = _build_psf_frame_qc_table(idx).set_index("file")

    assert bool(table.loc["a.fits", "core_cut_enabled"]) is False
