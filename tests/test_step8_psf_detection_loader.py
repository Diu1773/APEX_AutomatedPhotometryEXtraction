from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("PyQt5")
import apex.gui.workflow.cmd.step8_psf_photometry as psf_module
from apex.analysis.psf_iteration import PSFFitFlag
from apex.gui.workflow.cmd.step8_psf_photometry import (
    _allstar_newton_group,
    _allstar_newton_one,
    _allstar_apply_model_inplace,
    _allstar_build_model,
    _allstar_fit,
    _build_groups,
    _build_psf_frame_qc_table,
    _build_psf_qc_summary,
    _float32_difference,
    _load_detect_positions,
)


def _gaussian_psf(dx, dy, sigma=1.2):
    values = np.exp(-0.5 * (np.asarray(dx) ** 2 + np.asarray(dy) ** 2) / sigma**2)
    return values / (2.0 * np.pi * sigma**2)


def test_model_accumulation_and_difference_keep_float32_work_arrays():
    shape = (31, 31)
    xy = np.array([[14.2, 15.1], [18.0, 16.0]])
    flux = np.array([1200.0, 300.0])
    model = _allstar_build_model(
        shape, xy[:, 0], xy[:, 1], flux, _gaussian_psf, 17
    )
    work = np.zeros(shape, dtype=np.float32)

    _allstar_apply_model_inplace(
        work,
        xy[:, 0],
        xy[:, 1],
        flux,
        _gaussian_psf,
        17,
        subtract=False,
    )
    residual = _float32_difference(work, model)

    assert work.dtype == np.float32
    assert residual.dtype == np.float32
    assert np.allclose(work, model)
    assert np.allclose(residual, 0.0)


def test_apex_fit_emits_standard_quality_and_fixed_position_results():
    shape = (41, 41)
    truth_xy = np.array([[18.2, 20.4], [22.1, 20.0]])
    truth_flux = np.array([3000.0, 1200.0])
    image = _allstar_build_model(
        shape,
        truth_xy[:, 0],
        truth_xy[:, 1],
        truth_flux,
        _gaussian_psf,
        17,
    )
    seeds = truth_xy + np.array([[0.15, -0.12], [-0.1, 0.08]])
    result = _allstar_fit(
        image,
        seeds,
        truth_flux * np.array([0.9, 1.1]),
        _gaussian_psf,
        fit_shape=9,
        stamp_size=17,
        max_iter=4,
        flux_conv=0.001,
        max_shift=2.0,
        background_rms=1.0,
        gain=1.0,
        initial_positions=seeds,
        position_bound=2.0,
    )
    assert {"qfit", "cfit", "reduced_chi2", "flags", "n_pixels_fit"} <= set(
        result.colnames
    )
    assert np.nanmedian(np.asarray(result["qfit"], dtype=float)) < 0.01
    assert np.max(np.hypot(
        np.asarray(result["x_fit"], dtype=float) - truth_xy[:, 0],
        np.asarray(result["y_fit"], dtype=float) - truth_xy[:, 1],
    )) < 0.05

    fixed_xy = np.column_stack([result["x_fit"], result["y_fit"]]).astype(float)
    fixed = _allstar_fit(
        image,
        fixed_xy,
        np.asarray(result["flux_fit"], dtype=float),
        _gaussian_psf,
        fit_shape=9,
        stamp_size=17,
        max_iter=1,
        flux_conv=0.001,
        background_rms=1.0,
        gain=1.0,
        initial_positions=fixed_xy,
        position_fixed=True,
    )
    assert np.allclose(fixed["x_fit"], fixed_xy[:, 0])
    assert np.allclose(fixed["y_fit"], fixed_xy[:, 1])


def test_forced_fixed_fit_allows_signed_flux_without_moving_position():
    shape = (31, 31)
    forced_xy = np.array([[15.25, 14.75]])
    signed_flux = np.array([-180.0])
    image = _allstar_build_model(
        shape,
        forced_xy[:, 0],
        forced_xy[:, 1],
        signed_flux,
        _gaussian_psf,
        17,
    )

    result = _allstar_fit(
        image,
        forced_xy,
        np.array([75.0]),
        _gaussian_psf,
        fit_shape=9,
        stamp_size=17,
        max_iter=2,
        flux_conv=0.001,
        background_rms=1.0,
        gain=1.0,
        initial_positions=forced_xy,
        position_fixed_mask=np.array([True]),
        allow_negative_flux_mask=np.array([True]),
    )

    assert np.allclose(result["x_fit"], forced_xy[:, 0])
    assert np.allclose(result["y_fit"], forced_xy[:, 1])
    assert np.isclose(float(result["flux_fit"][0]), signed_flux[0], atol=1e-6)
    assert int(result["flags"][0]) & int(PSFFitFlag.NONPOSITIVE_FLUX)


def test_apex_refit_preserves_a_prior_valid_solution(monkeypatch):
    shape = (21, 21)
    xy = np.array([[10.0, 10.0]])
    flux = np.array([1000.0])
    image = _allstar_build_model(
        shape, xy[:, 0], xy[:, 1], flux, _gaussian_psf, 17
    )

    def failed_update(*args, **kwargs):
        return 10.0, 10.0, 1000.0, np.nan, False

    monkeypatch.setattr(psf_module, "_allstar_newton_one", failed_update)
    prior_valid = _allstar_fit(
        image,
        xy,
        flux,
        _gaussian_psf,
        fit_shape=9,
        stamp_size=17,
        max_iter=1,
        flux_conv=0.01,
        initial_positions=xy,
        initial_fit_valid=np.array([True]),
    )
    never_valid = _allstar_fit(
        image,
        xy,
        flux,
        _gaussian_psf,
        fit_shape=9,
        stamp_size=17,
        max_iter=1,
        flux_conv=0.01,
        initial_positions=xy,
    )

    assert int(prior_valid["flags"][0]) & int(PSFFitFlag.NONCONVERGENCE) == 0
    assert int(never_valid["flags"][0]) & int(PSFFitFlag.NONCONVERGENCE)


def test_apex_local_refit_updates_only_active_sources():
    shape = (41, 41)
    xy = np.array([[10.0, 20.0], [30.0, 20.0]])
    truth_flux = np.array([1000.0, 800.0])
    image = _allstar_build_model(
        shape, xy[:, 0], xy[:, 1], truth_flux, _gaussian_psf, 17
    )

    result = _allstar_fit(
        image,
        xy,
        np.array([700.0, 500.0]),
        _gaussian_psf,
        fit_shape=9,
        stamp_size=17,
        max_iter=3,
        flux_conv=0.001,
        background_rms=1.0,
        gain=1.0,
        initial_positions=xy,
        initial_fit_valid=np.array([True, False]),
        position_fixed=True,
        fit_active_mask=np.array([False, True]),
    )

    assert float(result["flux_fit"][0]) == 700.0
    assert np.isclose(float(result["flux_fit"][1]), 800.0, rtol=0.01)
    assert int(result["flags"][0]) & int(PSFFitFlag.NONCONVERGENCE) == 0


def test_local_psf_groups_respect_grouped_source_budget():
    groups = _build_groups(
        np.array([0.0, 1.0, 10.0, 11.0, 20.0]),
        np.zeros(5),
        np.array([100.0, 90.0, 80.0, 70.0, 60.0]),
        radius=2.0,
        max_size=3,
        max_grouped_sources=2,
    )

    flattened = [index for group in groups for index in group]
    assert sorted(flattened) == list(range(5))
    assert max(map(len, groups)) <= 3
    assert sum(len(group) for group in groups if len(group) > 1) <= 2


def test_large_newton_group_uses_sparse_solver_contract():
    group = [
        (8.0, 8.0, 1000.0),
        (11.0, 8.0, 800.0),
        (8.0, 11.0, 600.0),
        (11.0, 11.0, 500.0),
    ]
    yy, xx = np.mgrid[:20, :20]
    image = np.zeros((20, 20), dtype=float)
    for x, y, flux in group:
        image += _gaussian_psf(xx - x, yy - y) * flux

    result = _allstar_newton_group(
        image,
        group,
        0,
        0,
        _gaussian_psf,
        max_shift=2.0,
        weights=np.ones_like(image),
    )

    assert len(result) == 4
    assert all(len(source_result) == 5 for source_result in result)
    assert all(source_result[-1] for source_result in result)


def test_newton_group_supports_mixed_fixed_and_signed_sources():
    truth = [(14.25, 15.10, -180.0), (16.30, 15.20, 900.0)]
    yy, xx = np.mgrid[:31, :31]
    image = np.zeros((31, 31), dtype=float)
    for x, y, flux in truth:
        image += _gaussian_psf(xx - x, yy - y) * flux

    result = _allstar_newton_group(
        image,
        [(14.25, 15.10, 50.0), (16.15, 15.30, 800.0)],
        0,
        0,
        _gaussian_psf,
        max_shift=2.0,
        weights=np.ones_like(image),
        position_fixed=np.array([True, False]),
        allow_negative_flux=np.array([True, False]),
    )

    assert all(source_result[-1] for source_result in result)
    assert np.allclose(result[0][:2], truth[0][:2])
    assert np.isclose(result[0][2], truth[0][2], atol=10.0)
    assert np.hypot(result[1][0] - truth[1][0], result[1][1] - truth[1][1]) < 0.05
    assert np.isclose(result[1][2], truth[1][2], rtol=0.03)


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
            "flags_psf": [0, 0, 8193, 0, 2],
            "forced_psf": [False, False, True, False, False],
            "flux_psf_e": [1000.0, 500.0, -5.0, 300.0, 10.0],
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
    assert summary.loc["g", "n_forced_negative"] == 1
    assert summary.loc["g", "n_crowding_unreliable"] == 1


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
            "n_forced": [4],
            "n_forced_negative": [1],
            "n_crowding_unreliable": [2],
            "frame_fwhm_px": [6.5],
            "psf_qc_status": ["REVIEW"],
            "psf_qc_score": [73.0],
            "psf_qc_reasons": ["crowding_warning"],
            "psf_clean_fraction": [0.8],
            "frame_total_elapsed_s": [12.5],
            "fit_elapsed_s": [8.0],
            "epsf_elapsed_s": [2.0],
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
    assert table.loc["a.fits", "n_forced_negative"] == 1
    assert table.loc["a.fits", "n_crowding_unreliable"] == 2
    assert table.loc["a.fits", "frame_fwhm_px"] == pytest.approx(6.5)
    assert table.loc["a.fits", "psf_qc_status"] == "REVIEW"
    assert table.loc["a.fits", "psf_clean_fraction"] == pytest.approx(0.8)
    assert table.loc["a.fits", "frame_total_elapsed_s"] == pytest.approx(12.5)
    assert table.loc["a.fits", "fit_elapsed_s"] == pytest.approx(8.0)
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
