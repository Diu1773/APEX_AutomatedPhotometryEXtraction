import numpy as np
import pandas as pd

from apex.analysis.photometric_qc import (
    FAIL,
    PASS,
    REVIEW,
    SKIP,
    PhotometricQCThresholds,
    evaluate_photometric_qc,
    summarize_photometric_qc,
)


def _synthetic_night(
    n_frames=20,
    n_stars=30,
    seed=42,
    cloud={},           # fname index -> mag offset added (dimming)
    noisy=set(),        # fname indices with inflated scatter
):
    rng = np.random.default_rng(seed)
    base = rng.uniform(12.0, 15.0, size=n_stars)
    frames = {}
    for j in range(n_frames):
        sigma = 0.15 if j in noisy else 0.005
        mags = base + rng.normal(0.0, sigma, size=n_stars) + cloud.get(j, 0.0)
        frames[f"frame_{j:03d}.fits"] = pd.DataFrame(
            {
                "source_id": np.arange(n_stars, dtype=np.int64),
                "mag_inst": mags,
                "snr": np.full(n_stars, 80.0),
                "filter": ["V"] * n_stars,
            }
        )
    return frames


def test_clean_night_all_pass_with_zero_offsets():
    qc = evaluate_photometric_qc(_synthetic_night())

    assert set(qc["phot_qc_status"]) == {PASS}
    assert np.nanmax(np.abs(qc["transparency_offset_mag"].to_numpy())) < 0.02
    assert summarize_photometric_qc(qc)[PASS] == 20


def test_cloud_frames_flagged_by_offset():
    qc = evaluate_photometric_qc(
        _synthetic_night(cloud={7: 0.5, 12: 1.0})
    ).set_index("file")

    row_review = qc.loc["frame_007.fits"]
    row_fail = qc.loc["frame_012.fits"]
    assert row_review["phot_qc_status"] == REVIEW
    assert "transparency_warning" in row_review["phot_qc_reasons"]
    assert abs(row_review["transparency_offset_mag"] - 0.5) < 0.05
    assert row_fail["phot_qc_status"] == FAIL
    assert "transparency_loss" in row_fail["phot_qc_reasons"]
    assert abs(row_fail["transparency_offset_mag"] - 1.0) < 0.05
    # Clear frames must not be dragged into REVIEW by the cloudy ones.
    clear = qc.drop(index=["frame_007.fits", "frame_012.fits"])
    assert set(clear["phot_qc_status"]) == {PASS}


def test_scatter_blowup_flagged():
    qc = evaluate_photometric_qc(_synthetic_night(noisy={4})).set_index("file")

    assert qc.loc["frame_004.fits", "phot_qc_status"] == REVIEW
    assert "frame_scatter" in qc.loc["frame_004.fits", "phot_qc_reasons"]


def test_too_few_stars_is_skip_not_fail():
    frames = _synthetic_night(n_stars=3)
    qc = evaluate_photometric_qc(frames)

    assert set(qc["phot_qc_status"]) == {SKIP}
    assert (qc["phot_qc_reasons"] == "too_few_qc_stars").all()


def test_faint_stars_excluded_by_snr_cut():
    frames = _synthetic_night()
    for df in frames.values():
        df.loc[df.index[:20], "snr"] = 5.0  # 20 of 30 stars too faint for QC
    qc = evaluate_photometric_qc(
        frames, PhotometricQCThresholds(min_stars=5)
    )

    assert set(qc["phot_qc_status"]) == {PASS}
    assert (qc["n_qc_stars"] <= 10).all()
