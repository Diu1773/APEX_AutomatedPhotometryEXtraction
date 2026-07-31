"""LC Step 1 night classification (apex.gui.workflow.lc.step1_night_setup).

Covers the local-noon split, its JD-gap fallback, and the handling of frames
whose observing time could not be read.
"""

from __future__ import annotations

from datetime import datetime

from apex.gui.workflow.lc.step1_night_setup import (
    NIGHT_ID_UNKNOWN,
    _build_night_summary,
    _classify_nights_by_jd_gap,
    classify_nights,
)
from apex.utils.night_utils import parse_dateobs

KST = 9.0                      # site tz offset used by the fixtures


def _jd(stamp: str) -> float:
    dt = parse_dateobs(stamp)
    return 2440587.5 + (dt - datetime(1970, 1, 1)).total_seconds() / 86400.0


def _rec(name: str, stamp: str | None, filt: str = "V") -> dict:
    return {
        "filename": name,
        "jd": _jd(stamp) if stamp else None,
        "date_obs": stamp or "",
        "filter": filt,
    }


def test_noon_split_keeps_one_session_together():
    """Evening / post-midnight / dawn frames of one night share a night_id."""
    records = [
        _rec("evening.fit", "2026-05-15T10:47:00"),   # 19:47 KST
        _rec("midnight.fit", "2026-05-15T15:10:00"),  # 00:10 KST (next day UTC-wise)
        _rec("dawn.fit", "2026-05-15T19:30:00"),      # 04:30 KST
    ]
    out = classify_nights(records, gap_hours=8.0, tz_offset_hours=KST)
    assert {r["night_id"] for r in out} == {1}
    assert {r["night_key"] for r in out} == {"20260515"}


def test_noon_split_separates_consecutive_nights():
    records = [
        _rec("n1_a.fit", "2026-05-15T15:10:00"),
        _rec("n2_a.fit", "2026-05-16T15:10:00"),
        _rec("n1_b.fit", "2026-05-15T19:00:00"),
    ]
    out = classify_nights(records, tz_offset_hours=KST)
    by_name = {r["filename"]: r["night_id"] for r in out}
    assert by_name["n1_a.fit"] == by_name["n1_b.fit"] == 1
    assert by_name["n2_a.fit"] == 2


def test_noon_split_survives_a_short_gap_across_midnight():
    """A JD-gap classifier with a small gap setting would split this in two;
    the noon split cannot, because the cut is fixed at local noon."""
    records = [
        _rec("a.fit", "2026-05-15T14:50:00"),   # 23:50 KST
        _rec("b.fit", "2026-05-15T15:10:00"),   # 00:10 KST
    ]
    out = classify_nights(records, gap_hours=8.0, tz_offset_hours=KST)
    assert {r["night_id"] for r in out} == {1}


def test_without_local_reference_jd_gap_is_kept():
    """tz unset (0.0) and no longitude: the Greenwich noon rule would mis-bucket,
    so the JD-gap classifier stays in charge."""
    records = [
        _rec("a.fit", "2026-05-15T10:47:00"),
        _rec("b.fit", "2026-05-15T15:10:00"),
    ]
    out = classify_nights(records, gap_hours=8.0, tz_offset_hours=0.0)
    assert {r["night_id"] for r in out} == {1}     # 4.4 h apart, same night
    assert "night_key" not in out[0]               # took the JD-gap path


def test_unreadable_time_becomes_night_zero_not_night_one():
    records = [
        _rec("good.fit", "2026-05-15T15:10:00"),
        _rec("broken.fit", None),
    ]
    noon = classify_nights(records, tz_offset_hours=KST)
    assert {r["filename"]: r["night_id"] for r in noon} == {
        "good.fit": 1, "broken.fit": NIGHT_ID_UNKNOWN}

    gap = _classify_nights_by_jd_gap(records, 8.0)
    assert {r["filename"]: r["night_id"] for r in gap} == {
        "good.fit": 1, "broken.fit": NIGHT_ID_UNKNOWN}


def test_summary_flags_the_unknown_night():
    records = [
        _rec("good.fit", "2026-05-15T15:10:00"),
        _rec("broken.fit", None),
    ]
    out = classify_nights(records, tz_offset_hours=KST)
    summary = {row["night_id"]: row for row in _build_night_summary(out, KST)}
    assert summary[NIGHT_ID_UNKNOWN]["unknown"] is True
    assert summary[NIGHT_ID_UNKNOWN]["label"] == "N?"
    assert summary[1]["unknown"] is False
    assert summary[1]["label"] == "N1"
