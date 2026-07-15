import numpy as np
import apex.analysis.psf_iteration as psf_iteration

from apex.analysis.psf_iteration import (
    PSFFitFlag,
    assess_psf_frame_quality,
    decide_residual_iteration,
    fit_parameters_changed,
    measure_psf_fit_quality,
    qfit_noise_diagnostics,
)


def _gaussian_psf(dx, dy, sigma=1.0):
    values = np.exp(-0.5 * (np.asarray(dx) ** 2 + np.asarray(dy) ** 2) / sigma**2)
    return values / (2.0 * np.pi * sigma**2)


def test_fit_parameter_convergence_includes_y_shift():
    assert fit_parameters_changed(10, 10, 100, 10, 10.1, 100, flux_fraction=0.01)
    assert not fit_parameters_changed(10, 10, 100, 10.001, 10.001, 100.1, flux_fraction=0.01)


def test_qfit_noise_normalization_accounts_for_snr_and_fit_pixels():
    expected, ratio = qfit_noise_diagnostics(
        np.array([2.0, 1.0, np.nan]),
        np.array([100.0, 50.0, 50.0]),
        np.array([10.0, 10.0, 0.0]),
        4.0,
    )

    np.testing.assert_allclose(expected[:2], [3.9894228, 1.9947114])
    np.testing.assert_allclose(ratio[:2], [0.50132565, 0.50132565])
    assert np.isnan(expected[2])
    assert np.isnan(ratio[2])


def test_psf_frame_assessment_passes_clean_well_sampled_fit():
    result = assess_psf_frame_quality(
        n_sources=1000,
        n_good=970,
        n_crowding_unreliable=5,
        median_qfit_noise_ratio=0.95,
        epsf_n_selected=30,
        epsf_median_contamination=0.08,
        frame_fwhm_px=7.0,
        frame_fwhm_max_px=10.0,
    )

    assert result.status == "PASS"
    assert result.clean_fraction == 0.97
    assert not result.reasons


def test_psf_frame_assessment_marks_difficult_but_usable_frame_for_review():
    result = assess_psf_frame_quality(
        n_sources=710,
        n_good=640,
        n_crowding_unreliable=27,
        median_qfit_noise_ratio=0.90,
        epsf_n_selected=13,
        epsf_median_contamination=0.061,
        frame_fwhm_px=10.32,
        frame_fwhm_max_px=10.18,
    )

    assert result.status == "REVIEW"
    assert "crowding_warning" in result.reasons
    assert "fwhm_above_config" in result.reasons


def test_psf_frame_assessment_fails_bad_epsf_and_fit_survival():
    result = assess_psf_frame_quality(
        n_sources=100,
        n_good=40,
        n_crowding_unreliable=20,
        median_qfit_noise_ratio=4.0,
        epsf_n_selected=2,
        epsf_median_contamination=0.7,
        frame_fwhm_px=14.0,
        frame_fwhm_max_px=10.0,
    )

    assert result.status == "FAIL"
    assert "low_clean_fraction" in result.reasons
    assert "high_qfit_noise_ratio" in result.reasons
    assert "too_few_epsf_stars" in result.reasons


def test_residual_convergence_uses_pre_cap_unique_candidates():
    decision = decide_residual_iteration(
        n_candidates_raw=100,
        n_candidates_unique=50,
        n_candidates_accepted=1,
        n_current=1000,
        convergence_fraction=0.02,
    )
    assert not decision.stop_now
    assert not decision.stop_after_refit
    assert decision.candidate_fraction == 0.05


def test_residual_iteration_stops_after_fitting_small_real_increment():
    decision = decide_residual_iteration(
        n_candidates_raw=10,
        n_candidates_unique=5,
        n_candidates_accepted=5,
        n_current=1000,
        convergence_fraction=0.01,
    )
    assert not decision.stop_now
    assert decision.stop_after_refit
    assert decision.reason == "candidate_fraction"


def test_quality_metrics_are_zero_for_perfect_model():
    yy, xx = np.mgrid[:21, :21]
    flux = np.array([1000.0])
    model = _gaussian_psf(xx - 10.0, yy - 10.0) * flux[0]
    quality = measure_psf_fit_quality(
        model,
        model,
        np.array([10.0]),
        np.array([10.0]),
        flux,
        _gaussian_psf,
        fit_shape=9,
        background_rms=2.0,
        gain=1.0,
    )
    assert quality.qfit[0] == 0.0
    assert quality.cfit[0] == 0.0
    assert quality.reduced_chi2[0] == 0.0
    assert np.isfinite(quality.flux_err[0])
    assert quality.flags[0] == 0


def test_quality_flags_failed_edge_fit_and_bound():
    yy, xx = np.mgrid[:11, :11]
    model = _gaussian_psf(xx - 1.0, yy - 1.0) * 100.0
    quality = measure_psf_fit_quality(
        model,
        model,
        np.array([1.0]),
        np.array([1.0]),
        np.array([100.0]),
        _gaussian_psf,
        fit_shape=9,
        background_rms=1.0,
        gain=1.0,
        fit_ok=np.array([False]),
        initial_xy=np.array([[0.0, 1.0]]),
        xy_bound=1.0,
    )
    flags = PSFFitFlag(int(quality.flags[0]))
    assert PSFFitFlag.INCOMPLETE_REGION in flags
    assert PSFFitFlag.NONCONVERGENCE in flags
    assert PSFFitFlag.NEAR_BOUND in flags


def test_quality_metrics_do_not_cast_full_frame_to_float64(monkeypatch):
    real_asarray = np.asarray

    def guarded_asarray(value, *args, **kwargs):
        dtype = kwargs.get("dtype", args[0] if args else None)
        array = real_asarray(value)
        if dtype is not None and array.ndim == 2 and array.size > 1000:
            raise AssertionError("full-frame dtype conversion")
        return real_asarray(value, *args, **kwargs)

    monkeypatch.setattr(psf_iteration.np, "asarray", guarded_asarray)
    yy, xx = np.mgrid[:64, :64]
    model = (_gaussian_psf(xx - 32.0, yy - 32.0) * 1000.0).astype(np.float32)

    quality = measure_psf_fit_quality(
        model,
        model,
        np.array([32.0]),
        np.array([32.0]),
        np.array([1000.0]),
        _gaussian_psf,
        fit_shape=17,
        background_rms=2.0,
        gain=1.0,
    )

    assert quality.qfit[0] == 0.0
    assert quality.reduced_chi2[0] == 0.0
