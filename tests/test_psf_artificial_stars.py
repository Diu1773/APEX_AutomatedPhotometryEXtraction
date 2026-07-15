from __future__ import annotations

import numpy as np
import pandas as pd

from apex.benchmark.psf_artificial_stars import (
    add_forced_truth_to_step7,
    aggregate_recovery_metrics,
    apply_recovery_quality_policy,
    inject_flux_catalog,
    match_injections_to_products,
    measure_preinjection_psf_residual,
    optimal_psf_flux_for_snr,
    oversampled_epsf_to_native_kernel,
    psf_noise_equivalent_area,
)


def test_oversampled_epsf_maps_to_normalized_native_kernel_with_expected_width():
    os = 2
    yy, xx = np.mgrid[-8:9, -8:9]
    data = np.exp(-0.5 * ((xx / (2.0 * os)) ** 2 + (yy / (2.0 * os)) ** 2))
    kernel = oversampled_epsf_to_native_kernel(
        data,
        header={"OVERSAMPL": os, "EPSFSIZE": 7},
    )
    assert kernel.shape == (7, 7)
    assert np.isclose(kernel.sum(), 1.0)
    y, x = np.indices(kernel.shape, dtype=float)
    sigma_x = np.sqrt(np.sum(kernel * (x - 3.0) ** 2))
    assert 1.5 < sigma_x < 2.5
    assert np.isclose(psf_noise_equivalent_area(kernel), 1.0 / np.sum(kernel * kernel))


def test_phase_aware_epsf_kernel_preserves_subpixel_centroid():
    os = 4
    yy, xx = np.mgrid[-24:25, -24:25]
    data = np.exp(-0.5 * ((xx / (0.9 * os)) ** 2 + (yy / (0.9 * os)) ** 2))
    phase_x, phase_y = 0.32, -0.27
    kernel = oversampled_epsf_to_native_kernel(
        data,
        header={"OVERSAMPL": os, "EPSFSIZE": 9},
        phase_x=phase_x,
        phase_y=phase_y,
    )

    y, x = np.indices(kernel.shape, dtype=float)
    center = kernel.shape[0] // 2
    assert np.isclose(kernel.sum(), 1.0)
    assert np.isclose(np.sum(kernel * (x - center)), phase_x, atol=0.03)
    assert np.isclose(np.sum(kernel * (y - center)), phase_y, atol=0.03)


def test_injection_uses_phase_aware_oversampled_kernel_without_native_reshift():
    os = 4
    yy, xx = np.mgrid[-24:25, -24:25]
    data = np.exp(-0.5 * ((xx / (0.8 * os)) ** 2 + (yy / (0.8 * os)) ** 2))
    reference = oversampled_epsf_to_native_kernel(
        data,
        header={"OVERSAMPL": os, "EPSFSIZE": 9},
    )

    def sample(phase_x, phase_y):
        return oversampled_epsf_to_native_kernel(
            data,
            header={"OVERSAMPL": os, "EPSFSIZE": 9},
            phase_x=phase_x,
            phase_y=phase_y,
        )

    catalog = pd.DataFrame({
        "x_true": [15.35],
        "y_true": [14.70],
        "true_flux_e": [100_000.0],
    })
    _, expected, _, _ = inject_flux_catalog(
        np.zeros((31, 31), dtype=float),
        reference,
        catalog,
        gain_e_per_adu=1.0,
        rng=np.random.default_rng(1),
        kernel_sampler=sample,
    )

    stamp = expected[11:20, 11:20] / 100_000.0
    assert np.allclose(stamp, sample(0.35, -0.30))


def test_injection_can_skip_full_frame_diagnostic_layers():
    kernel = np.ones((3, 3), dtype=float) / 9.0
    catalog = pd.DataFrame({
        "x_true": [8.0],
        "y_true": [8.0],
        "true_flux_e": [900.0],
    })

    injected, expected, realized, _ = inject_flux_catalog(
        np.zeros((17, 17), dtype=np.float32),
        kernel,
        catalog,
        gain_e_per_adu=1.0,
        rng=np.random.default_rng(2),
        return_layers=False,
    )

    assert injected.dtype == np.float32
    assert expected is None
    assert realized is None
    assert injected.sum() > 0


def test_preinjection_residual_projection_recovers_psf_amplitude_over_plane():
    yy, xx = np.mgrid[-4:5, -4:5]
    kernel = np.exp(-0.5 * (xx**2 + yy**2) / 1.1**2)
    kernel /= kernel.sum()
    image = np.zeros((31, 31), dtype=float)
    image += 8.0 + 0.02 * np.indices(image.shape)[1]
    image[11:20, 11:20] += 250.0 * kernel
    catalog = pd.DataFrame({
        "x_true": [15.0],
        "y_true": [15.0],
        "true_flux_adu": [1000.0],
        "target_snr": [50.0],
    })

    result = measure_preinjection_psf_residual(
        image,
        catalog,
        kernel_sampler=lambda _x, _y: kernel,
    )

    assert np.isclose(result.loc[0, "preinjection_psf_residual_adu"], 250.0)
    assert np.isclose(result.loc[0, "preinjection_psf_residual_frac"], 0.25)
    assert not bool(result.loc[0, "confusion_clean"])


def test_optimal_psf_snr_flux_inversion_includes_source_shot_noise():
    flux_e, flux_adu = optimal_psf_flux_for_snr(
        10.0,
        gain_e_per_adu=2.0,
        background_rms_adu=3.0,
        psf_nea_px=5.0,
    )
    expected_snr = flux_e / np.sqrt(flux_e + (2.0 * 3.0) ** 2 * 5.0)
    assert np.isclose(expected_snr, 10.0)
    assert np.isclose(flux_adu, flux_e / 2.0)
    assert optimal_psf_flux_for_snr(20.0, gain_e_per_adu=2.0, background_rms_adu=3.0, psf_nea_px=5.0)[0] > flux_e


def test_matching_is_one_to_one_and_does_not_use_product_ids():
    truth = pd.DataFrame({"injection_id": [1, 2], "x_true": [10.0, 10.4], "y_true": [10.0, 10.0], "true_flux_e": [100.0, 200.0], "gain_e_per_adu": [2.0, 2.0]})
    detections = pd.DataFrame({"det_uid": [900, 901], "x": [10.1, 10.6], "y": [10.0, 10.0]})
    psf = pd.DataFrame({"det_uid": [2, 1], "x_fit": [10.1, 10.6], "y_fit": [10.0, 10.0], "flux_psf_e": [110.0, 190.0], "flux_psf_err_e": [10.0, 10.0], "flux_psf_raw_e": [100.0, 200.0], "flags_psf": [0, 0]})
    recovery = match_injections_to_products(truth, detections, psf, radius_px=0.5)
    assert recovery["detection_recovered"].tolist() == [True, True]
    assert recovery["psf_recovered"].tolist() == [True, True]
    assert recovery["psf_row"].nunique() == 2
    assert np.allclose(recovery["flux_frac_error"], [0.1, -0.05])
    assert np.allclose(recovery["raw_flux_frac_error"], [0.0, 0.0])
    assert np.allclose(recovery["flux_pull"], [1.0, -1.0])


def test_aggregate_reports_robust_bias_rmse_and_completeness_by_stratum():
    recovery = pd.DataFrame({
        "target_snr": [5.0, 5.0, 20.0, 20.0],
        "radius_bin": ["0-0.25", "0-0.25", "0.5-0.75", "0.5-0.75"],
        "crowding_bin": ["0.75-1.5 FWHM"] * 4,
        "detection_recovered": [True, True, True, False],
        "psf_recovered": [True, False, True, False],
        "flux_frac_error": [0.1, np.nan, -0.2, np.nan],
        "delta_mag": [-0.1, np.nan, 0.2, np.nan],
    })
    summary = aggregate_recovery_metrics(recovery)
    overall = summary.loc[summary["scope"] == "overall"].iloc[0]
    snr5 = summary.loc[(summary["scope"] == "target_snr") & (summary["label"] == "5.0")].iloc[0]
    assert overall["n_truth"] == 4
    assert np.isclose(overall["psf_completeness"], 0.5)
    assert np.isclose(overall["flux_bias"], -0.05)
    assert overall["flux_rmse"] > 0
    assert snr5["n_valid_flux"] == 1
    assert np.isclose(snr5["psf_completeness"], 0.5)
    assert "snr_x_radius" in set(summary["scope"])
    assert "snr_x_crowding" in set(summary["scope"])


def test_aggregate_excludes_nonpositive_flux_failures_from_error_metrics():
    recovery = pd.DataFrame({
        "target_snr": [20.0, 20.0],
        "detection_recovered": [True, True],
        "psf_recovered": [True, False],
        "flux_frac_error": [0.1, -1.5],
        "delta_mag": [-0.1, np.nan],
    })

    overall = aggregate_recovery_metrics(recovery).iloc[0]

    assert overall["n_valid_flux"] == 1
    assert np.isclose(overall["flux_bias"], 0.1)
    assert np.isclose(overall["flux_rmse"], 0.1)


def test_artificial_truth_rows_become_forced_step7_seeds():
    step7 = pd.DataFrame({
        "file": ["frame.fit"],
        "filter": ["B"],
        "master_id": [12],
        "det_uid": [5],
        "x_fit": [1.0],
        "y_fit": [2.0],
        "flux_e": [300.0],
        "flux_net_adu": [200.0],
        "detected_flag": [True],
        "forced_flag": [False],
        "off_frame_flag": [False],
        "is_saturated": [False],
        "is_nonlinear": [False],
    })
    truth = pd.DataFrame({
        "x_true": [20.5],
        "y_true": [30.25],
        "true_flux_e": [900.0],
        "true_flux_adu": [600.0],
    })

    combined = add_forced_truth_to_step7(step7, truth, filename="frame.fit")
    injected = combined.iloc[-1]

    assert len(combined) == 2
    assert injected["master_id"] == 13
    assert bool(injected["forced_flag"])
    assert not bool(injected["detected_flag"])
    assert np.isclose(injected["x_fit"], 20.5)
    assert np.isclose(injected["flux_net_adu"], 600.0)

    generated = add_forced_truth_to_step7(
        pd.DataFrame(),
        truth,
        filename="frame.fit",
    )
    assert {"x_fit", "y_fit", "flux_net_adu", "forced_flag"} <= set(
        generated.columns
    )
    assert bool(generated.loc[0, "forced_flag"])


def test_recovery_quality_policy_uses_noise_normalized_qfit():
    recovery = pd.DataFrame({
        "detection_recovered": [True, True, True],
        "psf_positive_recovered": [True, True, True],
        "nearest_real_sep_fwhm": [1.49, 2.0, 2.0],
        "qfit": [0.1, 5.0, 1.01],
        "n_pixels_fit": [100, 100, 100],
        "snr_psf": [10.0, 10.0, 10.0],
        "psf_nea_px": [25.0, 25.0, 25.0],
        "cfit": [0.0, 0.0, 0.0],
        "reduced_chi2": [1.0, 1.0, 1.0],
        "flags_psf": [0, 0, 0],
    })

    quality = apply_recovery_quality_policy(recovery)

    assert quality["psf_recovered"].tolist() == [False, False, True]
    assert quality.loc[2, "qfit_noise_ratio"] < 1.0


def test_aggregate_reports_error_scale_and_systematic_floor():
    recovery = pd.DataFrame({
        "detection_recovered": [True, True, True],
        "psf_positive_recovered": [True, True, True],
        "psf_recovered": [True, True, True],
        "flux_frac_error": [-0.1, 0.0, 0.1],
        "delta_mag": [0.1, 0.0, -0.1],
        "flux_pull": [-5.0, 0.0, 5.0],
        "flux_frac_formal_err": [0.02, 0.02, 0.02],
    })

    overall = aggregate_recovery_metrics(recovery).iloc[0]

    assert overall["recommended_error_scale"] > 1.0
    assert overall["empirical_systematic_floor_frac"] > 0.0
