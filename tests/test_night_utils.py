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


# --- epoch span -------------------------------------------------------------

def test_night_span_days():
    assert nu.night_span_days(["20250429", "20250430"]) == 1
    assert nu.night_span_days(["20241106", "20250428"]) == 173
    assert nu.night_span_days(["20250429"]) is None          # one epoch
    assert nu.night_span_days([]) is None
    assert nu.night_span_days(["garbage", "20250429"]) is None
    # order must not matter
    assert nu.night_span_days(["20250428", "20241106"]) == 173
