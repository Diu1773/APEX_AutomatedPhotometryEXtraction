"""Observing-night determination — one definition for the whole app.

An observing night must not be split at midnight: the evening flats, the
lights taken after midnight and the dawn flats all belong to the same session.
The cut therefore has to fall while the Sun is up.  We cut at **local noon**,
which outside the polar regions always lies between sunrise and sunset, so no
real sunrise/sunset computation — and no latitude — is needed: any cut made
during daylight yields the same grouping.

The local reference is resolved in this order:

1. **site longitude** (FITS ``SITELONG`` & friends) → local *solar* time
   (``UTC + lon/15 h``).  Most faithful: noon is the Sun's transit.
2. **configured tz offset** (``[site] tz_offset_hours`` →
   ``params.P.site_tz_offset_hours``) → local *civil* time.
3. neither → ``UTC - 12 h``.  This is only correct near Greenwich; for an
   East-Asian site it tears an evening off onto the previous night, so callers
   that have another signal (a date in the path, a JD-gap classifier) should
   prefer that signal instead of this fallback.  :func:`has_local_reference`
   tells them whether step 1 or 2 succeeded.

``tz_offset_hours == 0.0`` is treated as "unset" because it is also the config
default and cannot be told apart from a genuine Greenwich site — and for a
genuine Greenwich site the step-3 fallback returns exactly the same answer, so
nothing is lost by that choice.
"""

from __future__ import annotations

import re
from datetime import datetime, timedelta
from typing import Optional, Tuple

__all__ = [
    "LON_HEADER_KEYS",
    "DATE_HEADER_KEYS",
    "METHOD_SOLAR",
    "METHOD_CIVIL",
    "METHOD_UTC",
    "parse_lon_east",
    "parse_dateobs",
    "night_offset_hours",
    "has_local_reference",
    "observing_night",
    "observing_night_detail",
    "observing_night_from_jd",
    "jd_to_datetime",
]

LON_HEADER_KEYS = ("SITELONG", "SITELON", "LONG-OBS", "LONGITUD", "OBSGEO-L",
                   "OBSLONG", "TELLONG")
DATE_HEADER_KEYS = ("DATE-OBS", "DATE", "DATEOBS")

METHOD_SOLAR = "solar"    # local solar time from the site longitude
METHOD_CIVIL = "civil"    # local civil time from the configured tz offset
METHOD_UTC = "utc"        # no local reference — UTC - 12 h

# JD 2440587.5 is 1970-01-01T00:00:00 UTC.
_UNIX_EPOCH_JD = 2440587.5
_UNIX_EPOCH = datetime(1970, 1, 1)

_TZ_SUFFIX_RE = re.compile(r"([+-])(\d{2}):?(\d{2})?$")


def parse_lon_east(value) -> Optional[float]:
    """Parse a site-longitude header to signed degrees east (west negative).

    Accepts a float, ``'127 21 37'``, ``'127:21:37 E'``, ``'70 W'``.  Returns
    ``None`` when the value is missing or cannot be read as a longitude.
    """
    if value is None:
        return None
    if isinstance(value, (int, float)):
        lon = float(value)
        return lon if abs(lon) <= 180.0 else None
    s = str(value).strip()
    if not s:
        return None
    upper = s.upper()
    sign = -1.0 if (s[0] == "-" or "W" in upper) else 1.0
    nums = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", s)]
    if not nums:
        return None
    deg = nums[0]
    if len(nums) > 1:
        deg += nums[1] / 60.0
    if len(nums) > 2:
        deg += nums[2] / 3600.0
    return sign * deg if abs(deg) <= 180.0 else None


def parse_dateobs(value) -> Optional[datetime]:
    """Parse a FITS ``DATE-OBS`` into a naive UTC :class:`datetime`.

    FITS dates are UTC by convention; a trailing ``+HH:MM`` / ``-HH:MM`` offset
    written by some capture software is honoured and converted to UTC rather
    than dropped (dropping it could shift a frame by up to 14 h — exactly the
    kind of silent mis-bucketing this module exists to prevent).
    """
    if value is None:
        return None
    s = str(value).strip()
    if not s:
        return None
    s = s.replace("Z", "").replace("z", "").strip()
    if "T" not in s and " " in s:
        s = s.replace(" ", "T", 1)

    offset = timedelta(0)
    # Only a zone suffix on a full timestamp; a bare '2026-05-15' has no time.
    if "T" in s:
        match = _TZ_SUFFIX_RE.search(s)
        if match:
            hours = int(match.group(2))
            minutes = int(match.group(3) or 0)
            offset = timedelta(hours=hours, minutes=minutes)
            if match.group(1) == "-":
                offset = -offset
            s = s[: match.start()].strip()

    for fmt in ("%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S",
                "%Y-%m-%dT%H:%M", "%Y-%m-%d"):
        try:
            dt = datetime.strptime(s, fmt)
        except ValueError:
            continue
        return dt - offset
    return None


def jd_to_datetime(jd) -> Optional[datetime]:
    """Julian Date (UTC) to a naive UTC :class:`datetime`."""
    try:
        jd_val = float(jd)
    except (TypeError, ValueError):
        return None
    if jd_val <= 0 or jd_val != jd_val:            # non-positive or NaN
        return None
    try:
        return _UNIX_EPOCH + timedelta(seconds=(jd_val - _UNIX_EPOCH_JD) * 86400.0)
    except (OverflowError, OSError, ValueError):
        return None


def night_offset_hours(lon_east_deg=None,
                       tz_offset_hours=None) -> Tuple[float, str]:
    """Resolve the local-time offset used for the noon split.

    Returns ``(hours_east_of_utc, method)`` where method is one of
    :data:`METHOD_SOLAR`, :data:`METHOD_CIVIL`, :data:`METHOD_UTC`.
    """
    lon = parse_lon_east(lon_east_deg)
    if lon is not None:
        return lon / 15.0, METHOD_SOLAR
    if tz_offset_hours is not None:
        try:
            tz = float(tz_offset_hours)
        except (TypeError, ValueError):
            tz = 0.0
        if tz != 0.0 and abs(tz) <= 14.0:
            return tz, METHOD_CIVIL
    return 0.0, METHOD_UTC


def has_local_reference(lon_east_deg=None, tz_offset_hours=None) -> bool:
    """True when a site longitude or a non-zero tz offset is available.

    When this is False the noon split degenerates to ``UTC - 12 h``; callers
    with another signal (path date, JD gaps) should use that instead.
    """
    return night_offset_hours(lon_east_deg, tz_offset_hours)[1] != METHOD_UTC


def _night_key(dt_utc: datetime, offset_hours: float) -> str:
    return (dt_utc + timedelta(hours=offset_hours)
            - timedelta(hours=12)).strftime("%Y%m%d")


def observing_night_detail(date_obs,
                           lon_east_deg=None,
                           tz_offset_hours=None) -> Tuple[str, str]:
    """``(night_key, method)`` for one ``DATE-OBS``; ``("", method)`` if unparsable."""
    offset, method = night_offset_hours(lon_east_deg, tz_offset_hours)
    dt = parse_dateobs(date_obs)
    if dt is None:
        return "", method
    return _night_key(dt, offset), method


def observing_night(date_obs, lon_east_deg=None, tz_offset_hours=None) -> str:
    """Observing night as ``YYYYMMDD`` from ``DATE-OBS``; ``""`` if unparsable."""
    return observing_night_detail(date_obs, lon_east_deg, tz_offset_hours)[0]


def observing_night_from_jd(jd, lon_east_deg=None, tz_offset_hours=None) -> str:
    """Observing night as ``YYYYMMDD`` from a Julian Date; ``""`` if unusable."""
    dt = jd_to_datetime(jd)
    if dt is None:
        return ""
    offset, _method = night_offset_hours(lon_east_deg, tz_offset_hours)
    return _night_key(dt, offset)
