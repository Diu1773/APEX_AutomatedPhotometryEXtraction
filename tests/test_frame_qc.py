import numpy as np
import pandas as pd

from apex.analysis.frame_qc import (
    FAIL,
    PASS,
    REVIEW,
    FrameQCThresholds,
    evaluate_frame_qc,
    summarize_frame_qc,
)


def _base_frames(n=12):
    return pd.DataFrame(
        {
            "file": [f"frame_{i:03d}.fits" for i in range(n)],
            "filter": ["V"] * n,
            "time_index": np.arange(n),
            "airmass": np.linspace(1.0, 1.5, n),
            "sky_med": np.full(n, 1000.0),
            "sky_sigma": np.full(n, 32.0),
            "gain_e_per_adu": np.ones(n),
            "rdnoise_e": np.full(n, 5.0),
            "fwhm_med": np.full(n, 5.0),
            "n_sources": np.full(n, 120),
            "elong_med": np.full(n, 1.05),
            "n_anchor_candidates": np.full(n, 20),
            "n_epsf_candidates": np.full(n, 8),
            "n_psf_seed_candidates": np.full(n, 50),
            "quality_score_median": np.full(n, 95.0),
            "sat_star_count": np.zeros(n),
        }
    )


def test_evaluate_frame_qc_passes_stable_frames():
    out = evaluate_frame_qc(_base_frames())

    assert set(out["qc_status"]) == {PASS}
    assert summarize_frame_qc(out) == {PASS: 12, REVIEW: 0, FAIL: 0}
    assert np.isfinite(out["fwhm_model_px"]).all()
    assert np.isfinite(out["sky_sigma_expected_e"]).all()


def test_evaluate_frame_qc_fails_missing_fwhm_and_high_elongation():
    df = _base_frames()
    df.loc[2, "fwhm_med"] = np.nan
    df.loc[5, "elong_med"] = 1.6

    out = evaluate_frame_qc(df)

    assert out.loc[2, "qc_status"] == FAIL
    assert "fwhm_missing" in out.loc[2, "qc_reasons"]
    assert out.loc[5, "qc_status"] == FAIL
    assert "high_elongation" in out.loc[5, "qc_reasons"]


def test_evaluate_frame_qc_keeps_sky_outlier_as_review_when_sources_survive():
    df = _base_frames()
    df.loc[4, "sky_med"] = 9000.0
    df.loc[4, "sky_sigma"] = 280.0

    out = evaluate_frame_qc(df)

    assert out.loc[4, "qc_status"] == REVIEW
    assert "sky_warning" in out.loc[4, "qc_reasons"]


def test_evaluate_frame_qc_fails_source_count_collapse():
    df = _base_frames()
    df.loc[7, "n_sources"] = 2

    out = evaluate_frame_qc(df)

    assert out.loc[7, "qc_status"] == FAIL
    assert "low_n_sources" in out.loc[7, "qc_reasons"]


def test_evaluate_frame_qc_missing_candidate_columns_do_not_force_review():
    df = _base_frames().drop(
        columns=["n_anchor_candidates", "n_epsf_candidates", "n_psf_seed_candidates"]
    )

    out = evaluate_frame_qc(df)

    assert set(out["qc_status"]) == {PASS}


def test_evaluate_frame_qc_tight_config_elongation_keeps_review_below_fail():
    class _Params:
        fwhm_elong_max = 1.15  # tighter than the default review cut (1.22)

    df = _base_frames()
    df.loc[3, "elong_med"] = 1.10  # between (fail - 0.08) and fail

    out = evaluate_frame_qc(df, params=_Params())

    # Without the clamp the review threshold (1.22) would sit ABOVE the fail
    # threshold (1.15) and 1.10 could never be flagged for review.
    assert out.loc[3, "qc_status"] == REVIEW
    assert "elongation_warning" in out.loc[3, "qc_reasons"]

    df.loc[3, "elong_med"] = 1.18  # above the configured fail cut
    out = evaluate_frame_qc(df, params=_Params())
    assert out.loc[3, "qc_status"] == FAIL
    assert "high_elongation" in out.loc[3, "qc_reasons"]


def test_evaluate_frame_qc_uses_custom_z_thresholds():
    df = _base_frames()
    df.loc[6, "fwhm_med"] = 7.0

    loose = evaluate_frame_qc(
        df,
        thresholds=FrameQCThresholds(fwhm_z_review=8.0, fwhm_z_fail=9.0),
    )
    strict = evaluate_frame_qc(
        df,
        thresholds=FrameQCThresholds(fwhm_z_review=1.0, fwhm_z_fail=1.5),
    )

    assert loose.loc[6, "qc_status"] in {PASS, REVIEW}
    assert strict.loc[6, "qc_status"] == FAIL
