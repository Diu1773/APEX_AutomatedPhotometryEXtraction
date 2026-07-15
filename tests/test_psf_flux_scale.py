from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from apex.analysis.psf_flux_scale import (
    apply_psf_aperture_scale,
    estimate_psf_aperture_scale,
)


def _catalogues(count: int = 20, scale: float = 1.20):
    uid = np.arange(count)
    aperture_flux = np.linspace(1000.0, 5000.0, count)
    psf = pd.DataFrame({
        "det_uid": uid,
        "flux_psf_e": aperture_flux * scale,
        "flux_psf_err_e": aperture_flux * 0.02,
        "mag_psf": 25.0 - 2.5 * np.log10(aperture_flux * scale / 10.0),
        "flags_psf": 0,
        "neighbor_dist_fwhm": 5.0,
        "crowding_unreliable_psf": False,
    })
    aperture = pd.DataFrame({
        "det_uid": uid,
        "flux_e": aperture_flux,
        "snr": 100.0,
        "detected_flag": True,
        "bad_phot_flag": False,
        "off_frame_flag": False,
        "step4_apcorr_candidate": True,
    })
    return psf, aperture


def test_estimate_psf_aperture_scale_recovers_known_flux_ratio():
    psf, aperture = _catalogues(scale=1.20)
    psf.loc[0, "flux_psf_e"] *= 4.0

    result, references = estimate_psf_aperture_scale(psf, aperture, min_stars=8)

    assert result.applied
    assert result.scale == pytest.approx(1.0 / 1.20, rel=1e-6)
    assert result.n_used == 19
    assert references["eligible"].sum() == 20
    assert references["used"].sum() == 19


def test_estimate_psf_aperture_scale_rejects_crowded_and_low_snr_stars():
    psf, aperture = _catalogues(count=10)
    psf.loc[:4, "neighbor_dist_fwhm"] = 1.0
    aperture.loc[5:7, "snr"] = 10.0

    result, _ = estimate_psf_aperture_scale(psf, aperture, min_stars=3)

    assert not result.applied
    assert result.reason == "too_few_candidates"
    assert result.n_candidates == 2


def test_apply_psf_aperture_scale_preserves_raw_flux_and_rescales_error():
    psf, aperture = _catalogues(scale=1.25)
    result, _ = estimate_psf_aperture_scale(psf, aperture, min_stars=8)

    corrected = apply_psf_aperture_scale(psf, result, zeropoint=25.0, exptime=10.0)

    np.testing.assert_allclose(corrected["flux_psf_raw_e"], psf["flux_psf_e"])
    np.testing.assert_allclose(corrected["flux_psf_e"], aperture["flux_e"])
    np.testing.assert_allclose(
        corrected["flux_psf_err_e"], psf["flux_psf_err_e"] * result.scale
    )
    assert corrected["psf_aperture_scale_applied"].all()


def test_missing_legacy_apcorr_candidate_column_does_not_disable_scaling():
    psf, aperture = _catalogues(scale=1.10)
    aperture = aperture.drop(columns="step4_apcorr_candidate")

    result, _ = estimate_psf_aperture_scale(psf, aperture, min_stars=8)

    assert result.applied
    assert result.scale == pytest.approx(1.0 / 1.10, rel=1e-6)
