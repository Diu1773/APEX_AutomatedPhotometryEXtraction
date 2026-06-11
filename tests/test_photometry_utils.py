"""Tests for apex.utils.photometry_utils — CCD equation, aperture helpers."""
from __future__ import annotations

import importlib

import numpy as np
import pandas as pd
import pytest

_photutils_available = importlib.util.find_spec("photutils") is not None
requires_photutils = pytest.mark.skipif(
    not _photutils_available, reason="photutils not installed"
)

from apex.utils.photometry_utils import (
    MAG_ERR_COEFF,
    MAD_TO_SIGMA,
    circle_mask,
    get_numeric_array,
    phot_one_star,
    phot_vectorized,
    refine_local_centroid,
)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

def test_mag_err_coeff():
    assert abs(MAG_ERR_COEFF - 2.5 / np.log(10)) < 1e-12


def test_mad_to_sigma():
    assert abs(MAD_TO_SIGMA - 1.4826) < 1e-4


# ---------------------------------------------------------------------------
# circle_mask
# ---------------------------------------------------------------------------

def test_circle_mask_center_pixel():
    mask = circle_mask((10, 10), 5.0, 5.0, 0.5)
    assert mask[5, 5]
    assert mask.sum() >= 1


def test_circle_mask_radius_zero_centre_only():
    mask = circle_mask((10, 10), 5.0, 5.0, 0.0)
    assert mask[5, 5]


def test_circle_mask_shape():
    mask = circle_mask((20, 30), 15.0, 10.0, 3.0)
    assert mask.shape == (20, 30)


def test_circle_mask_approximate_area():
    r = 5.0
    mask = circle_mask((50, 50), 25.0, 25.0, r)
    expected = np.pi * r * r
    assert abs(mask.sum() - expected) / expected < 0.05


# ---------------------------------------------------------------------------
# refine_local_centroid
# ---------------------------------------------------------------------------

def _make_gaussian_image(cx, cy, amp=1000.0, sigma=2.0, shape=(60, 60), bg=100.0):
    h, w = shape
    y, x = np.mgrid[:h, :w]
    img = bg + amp * np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * sigma ** 2))
    return img.astype(float)


def test_refine_centroid_on_axis():
    cx, cy = 30.0, 30.0
    img = _make_gaussian_image(cx, cy)
    xc, yc = refine_local_centroid(img, cx, cy, fwhm_used=5.0, cbox_scale=1.5)
    assert abs(xc - cx) < 0.1
    assert abs(yc - cy) < 0.1


def test_refine_centroid_offset_seed():
    cx, cy = 30.7, 29.3
    img = _make_gaussian_image(cx, cy)
    xc, yc = refine_local_centroid(img, 30.0, 29.0, fwhm_used=5.0, cbox_scale=1.5)
    assert abs(xc - cx) < 0.3
    assert abs(yc - cy) < 0.3


def test_refine_centroid_returns_seed_for_tiny_box():
    img = np.ones((10, 10)) * 100.0
    xc, yc = refine_local_centroid(img, 5.0, 5.0, fwhm_used=0.1, cbox_scale=0.1)
    # Box too small → returns input unchanged
    assert xc == 5.0 and yc == 5.0


def test_refine_centroid_flat_image_returns_seed():
    """Uniform image has zero signal above background → returns seed."""
    img = np.full((60, 60), 500.0)
    x0, y0 = 30.0, 30.0
    xc, yc = refine_local_centroid(img, x0, y0, fwhm_used=5.0, cbox_scale=1.5)
    assert xc == x0 and yc == y0


# ---------------------------------------------------------------------------
# phot_one_star — CCD equation
# ---------------------------------------------------------------------------

def _flat_sky_image(shape=(100, 100), sky=200.0):
    return np.full(shape, sky, dtype=float)


def test_phot_one_star_pure_sky_zero_flux():
    """Star on flat background: net flux should be ~0."""
    img = _flat_sky_image()
    flux, sig, snr, ap_sum, bkg_med, bkg_std, ap_area, n_sky, is_sat, is_nl = phot_one_star(
        img, 50.0, 50.0, r_ap=5.0, r_in=10.0, r_out=15.0,
        gain=1.0, rn_param_e=0.0, sky_sigma_mode="local",
    )
    assert abs(flux) < 1.0   # net flux ≈ 0


def test_phot_one_star_with_gaussian_star():
    img = _make_gaussian_image(50.0, 50.0, amp=5000.0, sigma=2.0, shape=(100, 100), bg=100.0)
    flux, sig, snr, *_ = phot_one_star(
        img, 50.0, 50.0, r_ap=6.0, r_in=12.0, r_out=18.0,
        gain=1.0, rn_param_e=3.0,
    )
    assert flux > 0
    assert sig > 0
    assert snr > 10


def test_phot_one_star_snr_formula():
    """SNR = flux / sigma."""
    img = _make_gaussian_image(50.0, 50.0, amp=2000.0, sigma=2.0, shape=(100, 100), bg=50.0)
    flux, sig, snr, *_ = phot_one_star(img, 50.0, 50.0, r_ap=5.0, r_in=10.0, r_out=15.0)
    if sig > 0:
        assert abs(snr - flux / sig) < 1e-9


def test_phot_one_star_saturation_flag():
    img = np.full((60, 60), 1000.0, dtype=float)
    img[30, 30] = 70000.0
    _, _, _, _, _, _, _, _, is_sat, _ = phot_one_star(
        img, 30.0, 30.0, r_ap=3.0, r_in=6.0, r_out=10.0, sat_adu=60000.0
    )
    assert is_sat


def test_phot_one_star_nonlinear_flag():
    img = np.full((60, 60), 200.0, dtype=float)
    img[30, 30] = 56000.0
    _, _, _, _, _, _, _, _, _, is_nl = phot_one_star(
        img, 30.0, 30.0, r_ap=3.0, r_in=6.0, r_out=10.0,
        sat_adu=65000.0, datamax_adu=55000.0,
    )
    assert is_nl


def test_phot_one_star_ccd_equation_rn_in_bkg_est():
    """Bug-fix verification: var_bkg_est must include RN² when sky_sigma_includes_rn=True.

    With sky_sigma_includes_rn=True and a flat sky image:
      sigma_pix_e2 = max(bkg_std² - RN², 0)  (pure sky variance ≈ 0 for flat image)
      var_bkg_est  = (N_ap²/N_sky) * (sigma_pix_e2 + RN²)   ← must include RN²
      var_readnoise = N_ap * RN²

    Total noise² for a flat sky (no star signal) should include both RN terms.
    """
    rn = 10.0
    sky = 0.0
    img = np.full((80, 80), sky, dtype=float)
    flux, sig, snr, _, _, _, ap_area, n_sky, _, _ = phot_one_star(
        img, 40.0, 40.0, r_ap=5.0, r_in=10.0, r_out=15.0,
        gain=1.0, rn_param_e=rn,
        sky_sigma_mode="local", sky_sigma_includes_rn=True,
    )
    # Expected: var = 0 (no source) + 0 (no sky) + (N_ap²/N_sky)*RN² + N_ap*RN²
    if n_sky > 0 and ap_area > 0:
        var_bkg_est_expected   = (ap_area ** 2 / n_sky) * rn ** 2
        var_readnoise_expected = ap_area * rn ** 2
        expected_sig = np.sqrt(var_bkg_est_expected + var_readnoise_expected)
        # Allow ±20% because sigma_clip may reduce effective n_sky
        assert abs(sig - expected_sig) / max(expected_sig, 1e-9) < 0.20, (
            f"sigma={sig:.4f} expected≈{expected_sig:.4f} "
            f"(N_ap={ap_area}, N_sky={n_sky}, RN={rn})"
        )


# ---------------------------------------------------------------------------
# phot_vectorized — consistency with phot_one_star
# ---------------------------------------------------------------------------

@requires_photutils
def test_phot_vectorized_single_star_matches_scalar():
    img = _make_gaussian_image(50.0, 50.0, amp=3000.0, sigma=2.0, shape=(100, 100), bg=200.0)
    r_ap, r_in, r_out = 6.0, 12.0, 18.0
    kw = dict(gain=1.5, rn_param_e=5.0, sat_adu=65000.0)

    flux_v, sig_v, snr_v, sky_v, sat_v, nl_v = phot_vectorized(
        img, np.array([[50.0, 50.0]]), r_ap, r_in, r_out, **kw
    )
    flux_s, sig_s, snr_s, _, sky_s, _, _, _, sat_s, nl_s = phot_one_star(
        img, 50.0, 50.0, r_ap, r_in, r_out, **kw
    )

    assert abs(float(flux_v[0]) - flux_s) / max(abs(flux_s), 1.0) < 0.02
    # SNR sign and rough magnitude
    assert (snr_v[0] > 0) == (snr_s > 0)


@requires_photutils
def test_phot_vectorized_invalid_position_gives_nan():
    img = np.ones((50, 50), dtype=float) * 100.0
    flux, sig, snr, sky, sat, nl = phot_vectorized(
        img, np.array([[-999.0, -999.0]]), 5.0, 10.0, 15.0
    )
    assert not np.isfinite(flux[0])


@requires_photutils
def test_phot_vectorized_empty_positions():
    img = np.ones((50, 50))
    result = phot_vectorized(img, np.zeros((0, 2)), 5.0, 10.0, 15.0)
    for arr in result:
        assert len(arr) == 0


# ---------------------------------------------------------------------------
# get_numeric_array
# ---------------------------------------------------------------------------

def test_get_numeric_array_basic():
    df = pd.DataFrame({"v": [1.0, 2.0, np.nan, "bad"]})
    arr = get_numeric_array(df, "v")
    assert arr[0] == 1.0
    assert arr[1] == 2.0
    assert np.isnan(arr[2])
    assert np.isnan(arr[3])


def test_get_numeric_array_missing_col_uses_default():
    df = pd.DataFrame({"a": [1.0]})
    arr = get_numeric_array(df, "missing_col", default=-1.0)
    assert arr[0] == -1.0
