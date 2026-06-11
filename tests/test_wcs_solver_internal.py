from __future__ import annotations

import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

pytest.importorskip("scipy")

from apex.analysis.astrometry import solve as solve_wcs_internal


def _synthetic_wcs(
    *,
    ra_deg: float,
    dec_deg: float,
    scale_arcsec: float,
    shape: tuple[int, int],
    rotation_deg: float,
    crpix_offset: tuple[float, float] = (0.0, 0.0),
) -> WCS:
    ny, nx = shape
    scale_deg = scale_arcsec / 3600.0
    theta = np.deg2rad(rotation_deg)
    c = np.cos(theta)
    s = np.sin(theta)

    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crpix = [nx / 2.0 + crpix_offset[0], ny / 2.0 + crpix_offset[1]]
    w.wcs.crval = [ra_deg, dec_deg]
    w.wcs.cd = np.array(
        [
            [-scale_deg * c, scale_deg * s],
            [scale_deg * s, scale_deg * c],
        ]
    )
    return w


def test_internal_solver_recovers_shifted_center_hint():
    rng = np.random.default_rng(42)
    shape = (1024, 1024)
    approx_ra = 250.4218
    approx_dec = 36.4613
    scale = 0.42
    true_wcs = _synthetic_wcs(
        ra_deg=approx_ra,
        dec_deg=approx_dec,
        scale_arcsec=scale,
        shape=shape,
        rotation_deg=27.0,
        crpix_offset=(-210.0, 155.0),
    )

    n_true = 75
    true_xy = np.column_stack(
        [
            rng.uniform(80.0, shape[1] - 80.0, n_true),
            rng.uniform(80.0, shape[0] - 80.0, n_true),
        ]
    )
    sky = SkyCoord(true_wcs.pixel_to_world(true_xy[:, 0], true_xy[:, 1]))
    detected_xy = true_xy + rng.normal(0.0, 0.04, true_xy.shape)

    # Add unmatched detections to mimic a noisy/poor frame.  They are fainter,
    # so brightness ranking still has enough true stars to lock the pattern.
    false_xy = np.column_stack(
        [
            rng.uniform(0.0, shape[1], 30),
            rng.uniform(0.0, shape[0], 30),
        ]
    )
    sources_xy = np.vstack([detected_xy, false_xy])
    source_flux = np.concatenate([
        rng.uniform(1200.0, 5000.0, n_true),
        rng.uniform(20.0, 300.0, len(false_xy)),
    ])
    source_flux[::17] = np.nan

    result = solve_wcs_internal(
        sources_xy=sources_xy,
        source_flux=source_flux,
        gaia_ra=sky.ra.deg,
        gaia_dec=sky.dec.deg,
        gaia_flux=rng.uniform(12.0, 17.0, n_true),
        approx_ra=approx_ra,
        approx_dec=approx_dec,
        approx_scale_arcsec=scale,
        img_shape=shape,
        rotation_step_deg=6.0,
        scale_search_range=(0.90, 1.10),
        n_scale_steps=5,
        min_matches=12,
        translation_top_n=4,
        ransac_keep_candidates=4,
    )

    assert result.converged, "\n".join(result.log)
    assert result.n_matches >= 20
    assert result.rms_arcsec < 0.25
    assert any("RANSAC tried" in line for line in result.log)
    assert any("Selected candidate" in line for line in result.log)


def test_internal_solver_local_blind_recovers_large_center_offset():
    rng = np.random.default_rng(7)
    shape = (1024, 1024)
    approx_ra = 250.4218
    approx_dec = 36.4613
    scale = 0.42
    true_wcs = _synthetic_wcs(
        ra_deg=approx_ra,
        dec_deg=approx_dec,
        scale_arcsec=scale,
        shape=shape,
        rotation_deg=-38.0,
        crpix_offset=(780.0, -620.0),
    )

    n_true = 95
    true_xy = np.column_stack(
        [
            rng.uniform(90.0, shape[1] - 90.0, n_true),
            rng.uniform(90.0, shape[0] - 90.0, n_true),
        ]
    )
    sky = SkyCoord(true_wcs.pixel_to_world(true_xy[:, 0], true_xy[:, 1]))
    detected_xy = true_xy + rng.normal(0.0, 0.05, true_xy.shape)

    false_xy = np.column_stack(
        [
            rng.uniform(0.0, shape[1], 45),
            rng.uniform(0.0, shape[0], 45),
        ]
    )
    sources_xy = np.vstack([detected_xy, false_xy])
    source_flux = np.concatenate([
        rng.uniform(1500.0, 6000.0, n_true),
        rng.uniform(30.0, 400.0, len(false_xy)),
    ])

    result = solve_wcs_internal(
        sources_xy=sources_xy,
        source_flux=source_flux,
        gaia_ra=sky.ra.deg,
        gaia_dec=sky.dec.deg,
        gaia_flux=rng.uniform(12.0, 17.0, n_true),
        approx_ra=approx_ra,
        approx_dec=approx_dec,
        approx_scale_arcsec=scale,
        img_shape=shape,
        min_matches=12,
        n_brightest_src=120,
        n_brightest_cat=120,
        quad_k_neighbor=8,
        quad_neighbor_pool_factor=3,
        quad_max_per_side=2500,
        ransac_max_trials=3000,
        ransac_keep_candidates=5,
        local_blind=True,
        local_blind_radius_factor=2.8,
    )

    assert result.converged, "\n".join(result.log)
    assert result.n_matches >= 25
    assert result.rms_arcsec < 0.30
    assert any("local-blind catalog window" in line for line in result.log)


def test_internal_solver_adopts_sip_only_after_holdout_gain():
    rng = np.random.default_rng(9)
    shape = (1024, 1024)
    approx_ra = 250.4218
    approx_dec = 36.4613
    scale = 0.42
    true_wcs = _synthetic_wcs(
        ra_deg=approx_ra,
        dec_deg=approx_dec,
        scale_arcsec=scale,
        shape=shape,
        rotation_deg=0.0,
    )

    n_true = 160
    ideal_xy = np.column_stack(
        [
            rng.uniform(80.0, shape[1] - 80.0, n_true),
            rng.uniform(80.0, shape[0] - 80.0, n_true),
        ]
    )
    sky = SkyCoord(true_wcs.pixel_to_world(ideal_xy[:, 0], ideal_xy[:, 1]))

    center = np.array([shape[1] / 2.0, shape[0] / 2.0])
    delta = ideal_xy - center
    radius_norm = np.hypot(delta[:, 0] / center[0], delta[:, 1] / center[1])
    distorted_xy = (
        ideal_xy
        + delta * (0.010 * radius_norm[:, None] ** 2)
        + rng.normal(0.0, 0.03, ideal_xy.shape)
    )

    false_xy = np.column_stack(
        [
            rng.uniform(0.0, shape[1], 40),
            rng.uniform(0.0, shape[0], 40),
        ]
    )
    sources_xy = np.vstack([distorted_xy, false_xy])
    source_flux = np.concatenate([
        rng.uniform(2000.0, 8000.0, n_true),
        rng.uniform(20.0, 400.0, len(false_xy)),
    ])

    result = solve_wcs_internal(
        sources_xy=sources_xy,
        source_flux=source_flux,
        gaia_ra=sky.ra.deg,
        gaia_dec=sky.dec.deg,
        gaia_flux=rng.uniform(12.0, 17.0, n_true),
        approx_ra=approx_ra,
        approx_dec=approx_dec,
        approx_scale_arcsec=scale,
        img_shape=shape,
        min_matches=20,
        n_brightest_src=180,
        n_brightest_cat=160,
        quad_k_neighbor=8,
        quad_max_per_side=3000,
        ransac_inlier_radius_px=8.0,
        ransac_max_trials=4000,
        ransac_keep_candidates=5,
        sip_degree=3,
        sip_min_pairs=50,
        sip_holdout_fraction=0.25,
        sip_min_improvement=0.05,
    )

    assert result.converged, "\n".join(result.log)
    assert result.model == "SIP3"
    assert result.sip_order == 3
    assert result.rms_arcsec < 0.10
    assert any("SIP adopted" in line for line in result.log)
