"""Frame QC: the excess-detection gate.

Step 4's count test was one-sided — it failed frames with too FEW sources
(the Dragonfly NOBJ gate's family, danieli2020 S3.3.4) but had no counterpart
for too many. Over-detection is the failure that matters most for APEX
specifically: the master catalogue is the union of single-frame detections, so
a noise peak detected once is catalogued for good, where ALLFRAME's median
stack would have removed it before detection ever ran (stetson1994).

The thresholds are a measured envelope, not a guess. Over 194 real frames in
24 filter groups (M5, M13, M37, M67, NGC 6811; both the working and the
reprocessed reductions) the ratio n_sources / group-median tops out at 1.53x,
99th percentile 1.49x — and that 1.53 is itself a mixed-exposure group. Review
fires at 1.8x, fail at 2.5x. These tests pin the envelope so a future change
to the gate has to argue with the data.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from apex.analysis.frame_qc import FrameQCThresholds, evaluate_frame_qc

# Widest normal ratio seen across the 194-frame survey described above.
OBSERVED_NORMAL_MAX_RATIO = 1.53


def _group(counts, *, raw=None, filt="r", fwhm=5.0, sky=1000.0, sky_sigma=30.0):
    """A single filter group of otherwise healthy frames."""
    n = len(counts)
    frame = pd.DataFrame({
        "file": [f"f{i}.fit" for i in range(n)],
        "filter": [filt] * n,
        "n_sources": list(counts),
        "fwhm_med": [fwhm] * n,
        "sky_med": [sky] * n,
        "sky_sigma": [sky_sigma] * n,
        "elong_med": [1.05] * n,
        "gain_e_per_adu": [0.689] * n,
        "rdnoise_e": [2.1] * n,
        "sat_star_count": [0] * n,
        "quality_score_median": [85.0] * n,
    })
    if raw is not None:
        frame["n_raw_detections"] = list(raw)
    return frame


def _reasons(row) -> str:
    return str(row.get("qc_reasons", "") or "")


def test_thresholds_sit_above_the_observed_normal_envelope():
    """A gate that fires inside the measured normal range is a false-alarm
    generator; this is the arithmetic that keeps it outside."""
    thr = FrameQCThresholds()
    assert thr.nsrc_excess_review > OBSERVED_NORMAL_MAX_RATIO
    assert thr.nsrc_excess_fail > thr.nsrc_excess_review


def test_uniform_group_is_clean():
    out = evaluate_frame_qc(_group([1000, 1010, 990, 1005, 995]))
    assert not out["qc_reasons"].fillna("").str.contains("excess_detections").any()


def test_normal_envelope_does_not_trip_the_gate():
    """1.53x — the widest ratio in the real 194-frame survey — must pass."""
    counts = [1000, 1000, 1000, 1000, int(1000 * OBSERVED_NORMAL_MAX_RATIO)]
    out = evaluate_frame_qc(_group(counts))
    assert not out["qc_reasons"].fillna("").str.contains("excess_detections").any()
    ratio = out.loc[4, "n_sources_excess_ratio"]
    assert ratio == pytest.approx(OBSERVED_NORMAL_MAX_RATIO, abs=0.01)


def test_blowup_fails_and_names_itself():
    """The user's recalled case: ~4000 detections where ~1200 was normal."""
    out = evaluate_frame_qc(_group([1200, 1200, 1200, 1200, 4000]))
    row = out.loc[4]
    assert row["n_sources_excess_ratio"] == pytest.approx(4000 / 1200, abs=0.01)
    assert row["qc_status"] == "FAIL"
    assert "excess_detections" in _reasons(row)


def test_between_thresholds_is_review_not_fail():
    """Ambiguous excess is kept, per this module's conservative policy."""
    out = evaluate_frame_qc(_group([1000, 1000, 1000, 1000, 2000]))
    row = out.loc[4]
    assert row["qc_status"] == "REVIEW"
    assert "excess_detections_warning" in _reasons(row)


def test_raw_count_drives_the_gate_when_present():
    """n_sources is post-filter, so a blow-up can hide behind a normal-looking
    stored count — the raw yield is what has to be judged."""
    counts = [1000, 1000, 1000, 1000, 1000]          # all look normal
    raw = [1000, 1000, 1000, 1000, 4000]             # one extracted 4000
    out = evaluate_frame_qc(_group(counts, raw=raw))
    assert out.loc[4, "n_sources_excess_ratio"] == pytest.approx(4.0, abs=0.01)
    assert "excess_detections" in _reasons(out.loc[4])
    # and the frames whose raw count is normal stay clean
    assert not _reasons(out.loc[0])


def test_missing_raw_column_falls_back_to_n_sources():
    """Products written before n_raw_detections existed must still be judged."""
    out = evaluate_frame_qc(_group([1200, 1200, 1200, 1200, 4000]))
    assert "n_raw_detections" not in _group([1]).columns
    assert "excess_detections" in _reasons(out.loc[4])


def test_tiny_groups_are_not_judged():
    """With a group median below the guard the ratio is meaningless."""
    out = evaluate_frame_qc(_group([4, 4, 4, 40]))
    assert out["n_sources_excess_ratio"].isna().all()
    assert not out["qc_reasons"].fillna("").str.contains("excess_detections").any()


def test_low_side_gate_still_works():
    """The pre-existing low-count failure must survive the addition."""
    out = evaluate_frame_qc(_group([1000, 1000, 1000, 1000, 100]))
    assert "low_n_sources" in _reasons(out.loc[4])
    assert not np.isnan(out.loc[4, "n_sources_excess_ratio"])
