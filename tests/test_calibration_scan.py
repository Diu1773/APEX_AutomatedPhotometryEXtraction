"""Tests for Step-0 auto-scan / classification (apex.analysis.calibration_scan)."""

from __future__ import annotations

import numpy as np
import pytest
from astropy.io import fits

from apex.analysis import calibration_scan as cs


def _write(path, imagetyp, exptime=0.0, filt="", temp=-10.0, dateobs="2026-05-15T20:00:00"):
    hdu = fits.PrimaryHDU(np.zeros((4, 4), dtype=np.float32))
    hdu.header["IMAGETYP"] = imagetyp
    hdu.header["EXPTIME"] = exptime
    if filt:
        hdu.header["FILTER"] = filt
    hdu.header["SET-TEMP"] = temp
    hdu.header["DATE-OBS"] = dateobs
    hdu.writeto(path, overwrite=True)
    return str(path)


def test_classify_type_by_imagetyp():
    assert cs.classify_type("Bias Frame", "x.fit") == "bias"
    assert cs.classify_type("Dark Frame", "x.fit") == "dark"
    assert cs.classify_type("Flat Field", "x.fit") == "flat"
    assert cs.classify_type("Light Frame", "x.fit") == "light"
    assert cs.classify_type("MASTER BIAS", "x.fit") == "bias"
    assert cs.classify_type("Object", "x.fit") == "light"


def test_classify_type_filename_fallback():
    assert cs.classify_type("", "bias-001.fit") == "bias"
    assert cs.classify_type(None, "my_dark_10s.fits") == "dark"
    assert cs.classify_type("WEIRD", "flatfield.fit") == "flat"
    assert cs.classify_type("", "random.fit") is None


def test_night_from_path_and_dateobs():
    assert cs.night_from_path(r"E:\obs\M13_20260515\dark\d1.fit") == "20260515"
    assert cs.night_from_path(r"E:\obs\M13\d1.fit") == ""
    # noon split: 03:00 belongs to the previous evening's night
    assert cs.night_from_dateobs("2026-05-16T03:00:00") == "20260515"
    assert cs.night_from_dateobs("2026-05-15T20:00:00") == "20260515"


def test_night_rollover_uses_longitude():
    # East-Asian site (SITELONG 127.36 E): the UTC night sits ~10-20 h UTC,
    # straddling UTC noon. A plain -12 h split tears the evening flat (10:47 UTC)
    # onto the previous night; the longitude-aware split keeps the whole session
    # on one night, matching the lights (13:22 UTC).
    lon = cs._parse_lon_east("127 21 37")
    assert 127.3 < lon < 127.4
    flat = cs.night_from_dateobs("2026-05-15T10:47:52", lon)
    light = cs.night_from_dateobs("2026-05-15T13:22:30", lon)
    assert flat == light == "20260515"
    # without the longitude, the same evening frame mis-buckets a day earlier
    assert cs.night_from_dateobs("2026-05-15T10:47:52") == "20260514"
    # W longitude parses negative; bare float passes through
    assert cs._parse_lon_east("70 W") < 0
    assert cs._parse_lon_east(-70.5) == -70.5


def _header(dateobs=None, sitelong=None):
    hdr = fits.Header()
    if dateobs:
        hdr["DATE-OBS"] = dateobs
    if sitelong is not None:
        hdr["SITELONG"] = sitelong
    return hdr


def test_resolve_night_dateobs_outranks_path_date():
    """A post-midnight frame stamped with the NEXT day in its path must not
    start a second night — DATE-OBS with a local reference wins."""
    path = r"E:\obs\M13_20260516\light\L_0001.fit"
    hdr = _header("2026-05-15T15:10:00", "127 21 37")     # ~00:10 local
    night, method, conflict = cs.resolve_night(path, hdr)
    assert night == "20260515"
    assert method == "solar"
    assert conflict == "20260516"                          # reported, not used


def test_resolve_night_tz_offset_fallback():
    """No SITELONG in the header: the configured tz offset carries the split."""
    path = r"E:\obs\M13_20260516\light\L_0001.fit"
    hdr = _header("2026-05-15T15:10:00")
    night, method, _conflict = cs.resolve_night(path, hdr, tz_offset_hours=9.0)
    assert (night, method) == ("20260515", "civil")


def test_resolve_night_falls_back_to_path_without_local_reference():
    """Without longitude or tz the noon split degenerates to the Greenwich rule,
    which tears East-Asian evenings — so the path date is preferred over it."""
    path = r"E:\obs\M13_20260515\flat\F_0001.fit"
    hdr = _header("2026-05-15T10:47:00")                   # ~19:47 local
    night, method, _conflict = cs.resolve_night(path, hdr)
    assert (night, method) == ("20260515", "path")


def test_resolve_night_last_resort_is_utc_rule():
    hdr = _header("2026-05-15T20:00:00")
    night, method, _conflict = cs.resolve_night(r"E:\obs\M13\x.fit", hdr)
    assert (night, method) == ("20260515", "utc")


def test_resolve_night_evening_and_dawn_share_one_night():
    """The whole session buckets together regardless of the path dates."""
    hdr_kwargs = {"sitelong": "127 21 37"}
    flat = cs.resolve_night(r"E:\obs\20260515\f.fit",
                            _header("2026-05-15T10:47:00", **hdr_kwargs))[0]
    light = cs.resolve_night(r"E:\obs\20260516\l.fit",
                             _header("2026-05-15T16:00:00", **hdr_kwargs))[0]
    dawn = cs.resolve_night(r"E:\obs\20260516\d.fit",
                            _header("2026-05-15T19:30:00", **hdr_kwargs))[0]
    assert flat == light == dawn == "20260515"


def test_temp_bucket():
    assert cs.temp_bucket(4.8) == 5
    assert cs.temp_bucket(-9.9) == -10
    assert cs.temp_bucket(None) is None


def test_scan_and_group(tmp_path):
    (tmp_path / "n1").mkdir()
    for i in range(3):
        _write(tmp_path / "n1" / f"bias{i}.fit", "Bias Frame", 0.0)
        _write(tmp_path / "n1" / f"dark{i}.fit", "Dark Frame", 60.0, temp=5.0)
        _write(tmp_path / "n1" / f"flatV{i}.fit", "Flat Field", 3.0, filt="V")
        _write(tmp_path / "n1" / f"lightV{i}.fit", "Light Frame", 60.0, filt="V", temp=5.0)
    _write(tmp_path / "n1" / "lightB.fit", "Light Frame", 60.0, filt="B", temp=5.0)

    frames = cs.scan_folder(str(tmp_path))
    assert len(frames) == 13
    assert cs.nights(frames) == ["20260515"]

    g = cs.group_for_night(frames, "20260515")
    assert len(g["bias"]) == 3
    assert g["dark"][(60.0, 5)]
    assert len(g["flat"]["V"]) == 3
    assert len(g["light"]) == 4          # 3 V + 1 B

    # matching
    assert cs.match_dark(g["dark"], 60.0, 5.1) == (60.0, 5)
    assert cs.match_flat(g["flat"], "V") == "V"
    assert cs.match_flat(g["flat"], "B") is None    # no B flat available


def test_match_dark_nearest_exposure():
    darks = {(60.0, 5): [1], (120.0, 5): [1]}
    assert cs.match_dark(darks, 60.0, 5.0) == (60.0, 5)
    assert cs.match_dark(darks, 110.0, 5.0) == (120.0, 5)   # nearest exp


# --- dark temperature tolerance (P2) ---------------------------------------

def _dark_frame(exp, temp):
    return cs.FrameInfo(path=f"d_{exp}_{temp}.fit", ftype="dark", exp=exp, temp=temp)


def test_match_dark_uses_actual_temperature_not_the_bucket():
    """-10.4 and -10.6 land in different 1 °C buckets but are 0.2 °C apart;
    matching must rank on the real temperatures."""
    darks = {
        (60.0, -10): [_dark_frame(60.0, -10.4)],
        (60.0, -11): [_dark_frame(60.0, -10.6)],
    }
    match = cs.match_dark_detail(darks, 60.0, -10.55)
    assert match.key == (60.0, -11)                     # -10.6 is nearer
    assert match.delta_temp_c == pytest.approx(0.05)
    assert match.within_temp_tol


def test_match_dark_reports_temperature_mismatch():
    darks = {(60.0, 5): [_dark_frame(60.0, 5.0)]}
    match = cs.match_dark_detail(darks, 60.0, -10.0, tol_c=1.0)
    assert match.key == (60.0, 5)                       # still the only option
    assert match.delta_temp_c == pytest.approx(15.0)
    assert not match.within_temp_tol                    # caller warns / refuses


def test_match_dark_tolerance_is_configurable():
    darks = {(60.0, -10): [_dark_frame(60.0, -10.5)]}
    tight = cs.match_dark_detail(darks, 60.0, -10.0, tol_c=0.1)
    loose = cs.match_dark_detail(darks, 60.0, -10.0, tol_c=1.0)
    assert tight.delta_temp_c == loose.delta_temp_c == pytest.approx(0.5)
    assert not tight.within_temp_tol                    # 0.1 °C observer
    assert loose.within_temp_tol


def test_match_dark_reports_exposure_mismatch():
    darks = {(300.0, 5): [_dark_frame(300.0, 5.0)]}
    match = cs.match_dark_detail(darks, 10.0, 5.0)
    assert match.delta_exp_s == pytest.approx(290.0)


def test_match_dark_unknown_temperature_is_not_a_violation():
    darks = {(60.0, None): [_dark_frame(60.0, None)]}
    match = cs.match_dark_detail(darks, 60.0, None)
    assert match.delta_temp_c is None
    assert match.within_temp_tol                        # nothing to compare


def test_group_temperature_falls_back_to_bucket():
    assert cs.group_temperature([_dark_frame(60.0, -10.4)], (60.0, -10)) == pytest.approx(-10.4)
    assert cs.group_temperature([], (60.0, -10)) == pytest.approx(-10.0)
    assert cs.group_temperature([], (60.0, None)) is None
