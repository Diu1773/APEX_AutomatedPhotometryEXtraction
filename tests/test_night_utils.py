"""Tests for the shared observing-night definition (apex.utils.night_utils)."""

from __future__ import annotations

from datetime import datetime

import pytest

from apex.utils import night_utils as nu


# --- longitude parsing -----------------------------------------------------

def test_parse_lon_east_forms():
    assert 127.3 < nu.parse_lon_east("127 21 37") < 127.4
    assert 127.3 < nu.parse_lon_east("127:21:37 E") < 127.4
    assert nu.parse_lon_east("70 W") == pytest.approx(-70.0)
    assert nu.parse_lon_east(-70.5) == -70.5
    assert nu.parse_lon_east(None) is None
    assert nu.parse_lon_east("") is None
    assert nu.parse_lon_east("not a longitude") is None
    assert nu.parse_lon_east(999.0) is None          # out of range


# --- DATE-OBS parsing ------------------------------------------------------

def test_parse_dateobs_formats():
    assert nu.parse_dateobs("2026-05-15T20:00:00") == datetime(2026, 5, 15, 20, 0, 0)
    assert nu.parse_dateobs("2026-05-15T20:00:00.512") == datetime(2026, 5, 15, 20, 0, 0, 512000)
    assert nu.parse_dateobs("2026-05-15T20:00:00Z") == datetime(2026, 5, 15, 20, 0, 0)
    assert nu.parse_dateobs("2026-05-15 20:00:00") == datetime(2026, 5, 15, 20, 0, 0)
    assert nu.parse_dateobs("2026-05-15") == datetime(2026, 5, 15)
    assert nu.parse_dateobs("garbage") is None
    assert nu.parse_dateobs(None) is None


def test_parse_dateobs_converts_zone_offset_to_utc():
    # A '+09:00' stamp is 09:00 ahead of UTC; dropping it instead of converting
    # would move the frame by 9 h and could flip its observing night.
    assert nu.parse_dateobs("2026-05-16T05:00:00+09:00") == datetime(2026, 5, 15, 20, 0, 0)
    assert nu.parse_dateobs("2026-05-15T15:00:00-05:00") == datetime(2026, 5, 15, 20, 0, 0)


# --- the noon split --------------------------------------------------------

def test_night_never_splits_at_midnight():
    """Evening, post-midnight and dawn frames of one session share a night."""
    lon = 127.36                                     # East-Asian site
    evening = nu.observing_night("2026-05-15T10:47:00", lon)   # ~19:47 local
    midnight = nu.observing_night("2026-05-15T15:10:00", lon)  # ~00:10 local
    dawn = nu.observing_night("2026-05-15T19:30:00", lon)      # ~04:30 local
    assert evening == midnight == dawn == "20260515"


def test_greenwich_rule_tears_east_asian_evening():
    """Why the longitude/tz reference matters: the bare UTC rule mis-buckets."""
    stamp = "2026-05-15T10:47:00"
    assert nu.observing_night(stamp, 127.36) == "20260515"
    assert nu.observing_night(stamp) == "20260514"    # no local reference


def test_tz_offset_is_the_civil_fallback():
    stamp = "2026-05-15T10:47:00"
    assert nu.observing_night(stamp, None, tz_offset_hours=9.0) == "20260515"
    night, method = nu.observing_night_detail(stamp, None, tz_offset_hours=9.0)
    assert (night, method) == ("20260515", nu.METHOD_CIVIL)


def test_longitude_outranks_tz_offset():
    _night, method = nu.observing_night_detail(
        "2026-05-15T10:47:00", 127.36, tz_offset_hours=9.0)
    assert method == nu.METHOD_SOLAR


def test_zero_tz_offset_counts_as_unset():
    # 0.0 is the config default and cannot be told apart from a real Greenwich
    # site — but for Greenwich the UTC fallback gives the same answer anyway.
    assert not nu.has_local_reference(tz_offset_hours=0.0)
    assert nu.has_local_reference(tz_offset_hours=9.0)
    assert nu.has_local_reference(lon_east_deg=127.36)
    assert nu.observing_night("2026-05-15T10:47:00", None, 0.0) == \
        nu.observing_night("2026-05-15T10:47:00")


def test_unparsable_dateobs_returns_empty():
    assert nu.observing_night("garbage", 127.36) == ""
    assert nu.observing_night(None) == ""


# --- JD entry point --------------------------------------------------------

def test_jd_to_datetime_roundtrip():
    # JD 2440587.5 is the Unix epoch by definition.
    assert nu.jd_to_datetime(2440587.5) == datetime(1970, 1, 1)
    assert nu.jd_to_datetime(None) is None
    assert nu.jd_to_datetime(float("nan")) is None
    assert nu.jd_to_datetime(-1) is None


def test_observing_night_from_jd_matches_dateobs():
    stamp = "2026-05-15T15:10:00"
    dt = nu.parse_dateobs(stamp)
    jd = 2440587.5 + (dt - datetime(1970, 1, 1)).total_seconds() / 86400.0
    assert nu.observing_night_from_jd(jd, 127.36) == nu.observing_night(stamp, 127.36)


def test_observing_night_from_jd_bad_input():
    assert nu.observing_night_from_jd(None) == ""
    assert nu.observing_night_from_jd("nonsense") == ""


# --- headless night fallback ------------------------------------------------

def test_fallback_night_key_uses_the_noon_split_with_a_reference():
    # Post-midnight KST frame stays on the previous observing night.
    assert nu.fallback_night_key("2025-04-30T15:10:00", 9.0) == "20250430"
    assert nu.fallback_night_key("2025-04-29T17:47:00", 9.0) == "20250429"


def test_fallback_night_key_without_reference_is_the_calendar_date():
    """No tz/longitude: the Greenwich noon rule would tear an evening apart,
    so the plain DATE-OBS date is used (photometry_source_service rule)."""
    assert nu.fallback_night_key("2025-04-30T15:10:00") == "2025-04-30"
    assert nu.fallback_night_key("2025-04-30 15:10:00", 0.0) == "2025-04-30"
    assert nu.fallback_night_key("", 9.0) == ""
    assert nu.fallback_night_key(None) == ""


def test_fill_missing_night_ids_numbers_chronologically():
    ids, assigned = nu.fill_missing_night_ids(
        [0, 0, 0, 0],
        ["20250430", "20250429", "20250430", "20250429"],
    )
    assert assigned == {"20250429": 1, "20250430": 2}
    assert ids == [2, 1, 2, 1]


def test_fill_missing_night_ids_keeps_existing_and_starts_above_them():
    ids, assigned = nu.fill_missing_night_ids(
        [3, 0, 0], ["", "20250430", "20250429"], start_after=3)
    assert ids == [3, 5, 4]
    assert assigned == {"20250429": 4, "20250430": 5}


def test_fill_missing_night_ids_leaves_keyless_frames_at_zero():
    """A frame with no DATE-OBS must not be folded into an arbitrary night."""
    ids, assigned = nu.fill_missing_night_ids([0, 0], ["", "20250430"])
    assert ids == [0, 1]
    assert assigned == {"20250430": 1}


def test_fill_missing_night_ids_noop_when_all_assigned():
    ids, assigned = nu.fill_missing_night_ids([1, 2], ["x", "y"])
    assert ids == [1, 2] and assigned == {}


# --- epoch span -------------------------------------------------------------

def test_night_span_days():
    assert nu.night_span_days(["20250429", "20250430"]) == 1
    assert nu.night_span_days(["20241106", "20250428"]) == 173
    assert nu.night_span_days(["20250429"]) is None          # one epoch
    assert nu.night_span_days([]) is None
    assert nu.night_span_days(["garbage", "20250429"]) is None
    # order must not matter
    assert nu.night_span_days(["20250428", "20241106"]) == 173
