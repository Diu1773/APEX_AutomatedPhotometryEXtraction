from __future__ import annotations

from types import SimpleNamespace

import pandas as pd

from apex.utils.source_quality import compute_source_quality


def _params():
    return SimpleNamespace(
        saturation_adu=60000.0,
        datamax_adu=58000.0,
        ref_cat_max_elong=1.5,
        ref_cat_max_abs_round=0.5,
        ref_cat_sharp_min=0.2,
        ref_cat_sharp_max=1.0,
        psf_epsf_sharp_lo=0.2,
        psf_epsf_sharp_hi=1.0,
        psf_epsf_round_abs_max=0.5,
        psf_epsf_elong_max=1.3,
        psf_flux_percentile_lo=50.0,
        psf_flux_percentile_hi=100.0,
        psf_isolation_fwhm_mult=3.0,
        apcorr_isolation_factor=2.5,
    )


def test_compute_source_quality_marks_role_candidates():
    df = pd.DataFrame(
        {
            "x": [20.0, 80.0, 100.0],
            "y": [20.0, 80.0, 100.0],
            "dao_flux": [1000.0, 4000.0, 9000.0],
            "peak_adu": [1000.0, 4000.0, 9000.0],
            "fwhm_px": [4.0, 4.0, 4.0],
            "fwhm_status": ["ok", "ok", "ok"],
            "sharpness": [0.5, 0.5, 0.5],
            "roundness": [0.1, 0.1, 0.1],
            "elongation": [1.1, 1.1, 1.1],
        }
    )

    out = compute_source_quality(df, frame_fwhm_px=4.0, image_shape=(128, 128), params=_params())

    assert out.loc[2, "quality_score"] == 100.0
    assert bool(out.loc[2, "anchor_candidate"])
    assert bool(out.loc[2, "apcorr_candidate"])
    assert bool(out.loc[2, "epsf_candidate"])
    assert out.loc[2, "nearest_neighbor_fwhm"] > 3.0


def test_compute_source_quality_rejects_saturated_and_close_neighbor():
    df = pd.DataFrame(
        {
            "x": [40.0, 43.0, 90.0],
            "y": [40.0, 40.0, 90.0],
            "dao_flux": [5000.0, 4500.0, 9000.0],
            "peak_adu": [5000.0, 4500.0, 61000.0],
            "fwhm_px": [4.0, 4.0, 4.0],
            "fwhm_status": ["ok", "ok", "saturated"],
            "sharpness": [0.5, 0.5, 0.5],
            "roundness": [0.1, 0.1, 0.1],
            "elongation": [1.1, 1.1, 1.1],
        }
    )

    out = compute_source_quality(df, frame_fwhm_px=4.0, image_shape=(128, 128), params=_params())

    assert "close_neighbor" in out.loc[0, "quality_flags"]
    assert not bool(out.loc[0, "anchor_candidate"])
    assert "saturated" in out.loc[2, "quality_flags"]
    assert out.loc[2, "quality_score"] == 0.0
    assert not bool(out.loc[2, "apcorr_candidate"])


def test_unmeasured_fwhm_uses_frame_median_for_candidates():
    df = pd.DataFrame(
        {
            "x": [20.0, 80.0, 120.0],
            "y": [20.0, 80.0, 120.0],
            "dao_flux": [1000.0, 5000.0, 9000.0],
            "peak_adu": [1000.0, 5000.0, 9000.0],
            "fwhm_px": [float("nan"), float("nan"), float("nan")],
            "fwhm_status": ["not_measured", "not_measured", "not_measured"],
            "sharpness": [0.5, 0.5, 0.5],
            "roundness": [0.1, 0.1, 0.1],
            "elongation": [1.1, 1.1, 1.1],
        }
    )

    out = compute_source_quality(df, frame_fwhm_px=4.0, image_shape=(160, 160), params=_params())

    assert "fwhm_missing" in out.loc[2, "quality_flags"]
    assert out.loc[2, "fwhm_ratio_to_frame"] == 1.0
    assert bool(out.loc[2, "anchor_candidate"])
    assert bool(out.loc[2, "apcorr_candidate"])
    assert bool(out.loc[2, "epsf_candidate"])
