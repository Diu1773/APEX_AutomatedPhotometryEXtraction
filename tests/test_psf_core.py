import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from apex.utils.psf_core import (
    estimate_psf_core_cut,
    psf_core_keep_mask,
    target_pixel_from_wcs,
)


def test_estimate_psf_core_cut_auto_density_excludes_center():
    rng = np.random.default_rng(123)
    core = rng.normal(loc=(500.0, 500.0), scale=20.0, size=(250, 2))
    field = rng.uniform(0.0, 1000.0, size=(600, 2))
    xy = np.vstack([core, field])

    cut = estimate_psf_core_cut(
        xy,
        (1000, 1000),
        5.0,
        enabled=True,
        center_mode="auto",
        radius_px=80.0,
        auto_min_density_ratio=1.2,
    )

    assert cut.enabled
    assert abs(cut.center_x - 500.0) < 40.0
    assert abs(cut.center_y - 500.0) < 40.0
    assert cut.n_excluded > 100

    keep = psf_core_keep_mask(xy, cut)
    assert keep.dtype == bool
    assert np.sum(~keep) == cut.n_excluded


def test_estimate_psf_core_cut_stays_off_for_sparse_inputs():
    xy = np.array([[10.0, 10.0], [20.0, 20.0], [30.0, 30.0]])

    cut = estimate_psf_core_cut(xy, (100, 100), 4.0, enabled=True, auto_min_sources=10)

    assert not cut.enabled
    assert cut.reason == "too_few_sources"


def test_auto_radius_is_capped_in_fwhm_units():
    rng = np.random.default_rng(456)
    broad_cluster = rng.normal(loc=(500.0, 500.0), scale=180.0, size=(1200, 2))
    xy = np.clip(broad_cluster, 0.0, 999.0)

    cut = estimate_psf_core_cut(
        xy,
        (1000, 1000),
        5.0,
        enabled=True,
        center_mode="manual",
        manual_center=(500.0, 500.0),
        radius_fwhm_mult=20.0,
        max_exclude_frac=0.95,
    )

    assert cut.enabled
    assert cut.radius_px <= 100.0
    assert "fwhm_cap" in cut.method


def test_target_pixel_from_wcs_prefers_configured_sky_position():
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [51.0, 41.0]
    wcs.wcs.crval = [10.0, 20.0]
    wcs.wcs.cdelt = [-0.001, 0.001]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    pixel = target_pixel_from_wcs(wcs.to_header(), 10.0, 20.0, (100, 120))

    assert pixel is not None
    assert np.allclose(pixel, (50.0, 40.0))
    assert target_pixel_from_wcs(wcs.to_header(), 30.0, 20.0, (100, 120)) is None


def test_target_pixel_from_wcs_rejects_non_celestial_linear_header():
    assert target_pixel_from_wcs(fits.Header(), 229.638, 2.081, (3200, 4800)) is None


def test_auto_core_prefers_stable_surface_brightness_peak_when_available():
    rng = np.random.default_rng(789)
    detections = np.vstack([
        rng.normal(loc=(70.0, 60.0), scale=8.0, size=(180, 2)),
        rng.uniform(0.0, 200.0, size=(400, 2)),
    ])
    yy, xx = np.mgrid[:200, :200]
    image = rng.normal(0.0, 1.0, (200, 200))
    image += 100.0 * np.exp(-0.5 * ((xx - 100.0) ** 2 + (yy - 80.0) ** 2) / 15.0**2)

    cut = estimate_psf_core_cut(
        detections,
        image.shape,
        3.0,
        enabled=True,
        center_mode="auto",
        radius_px=30.0,
        auto_min_density_ratio=1.2,
        image=image,
    )

    assert cut.enabled
    assert cut.method.startswith("auto_surface_brightness")
    assert abs(cut.center_x - 100.0) < 8.0
    assert abs(cut.center_y - 80.0) < 8.0
