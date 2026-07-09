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
