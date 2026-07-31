"""Masters stacked across epochs must say so (D-007).

When a night has no calibration frames of its own, group_for_night falls back to
the shared pool, which groups on (exposure, temperature) only — so darks taken
months apart end up in one master and the hot-pixel pattern blurs. That used to
happen with no trace at all.
"""

from __future__ import annotations

import numpy as np
import pytest
from astropy.io import fits

from apex.analysis import calibration_scan as cs
from apex.analysis.calibration import CalibrationOptions
from apex.analysis.calibration_run import (
    EPOCH_WARN_DAYS,
    ALL_NIGHTS,
    epoch_info,
    run_calibration,
)


def _frame(night, exp=30.0, temp=-10.0):
    return cs.FrameInfo(path=f"d_{night}.fit", ftype="dark", exp=exp,
                        temp=temp, night=night)


def test_epoch_info_single_night():
    info = epoch_info([_frame("20250428"), _frame("20250428")])
    assert info["nights"] == ["20250428"]
    assert info["span_days"] is None
    assert info["mixed"] is False


def test_epoch_info_adjacent_nights_are_not_flagged():
    """Two nights of one observing run are normal."""
    info = epoch_info([_frame("20250429"), _frame("20250430")])
    assert info["span_days"] == 1
    assert info["mixed"] is False


def test_epoch_info_flags_a_months_apart_mix():
    """The real YZ Boo case: 2024-11-06 darks stacked with 2025-04-28 ones."""
    info = epoch_info([_frame("20241106"), _frame("20250428")])
    assert info["span_days"] == 173 > EPOCH_WARN_DAYS
    assert info["mixed"] is True
    assert info["counts"] == {"20241106": 1, "20250428": 1}


def test_epoch_info_ignores_frames_without_a_night():
    info = epoch_info([_frame(""), _frame("20250428")])
    assert info["nights"] == ["20250428"]
    assert info["mixed"] is False
    assert epoch_info([])["nights"] == []


# --- end to end -------------------------------------------------------------

def _write(path, imagetyp, exptime=0.0, filt="", temp=-10.0, dateobs=None):
    hdu = fits.PrimaryHDU(np.full((8, 8), 100.0, dtype=np.float32))
    hdu.header["IMAGETYP"] = imagetyp
    hdu.header["EXPTIME"] = exptime
    if filt:
        hdu.header["FILTER"] = filt
    hdu.header["SET-TEMP"] = temp
    hdu.header["DATE-OBS"] = dateobs
    hdu.header["SITELONG"] = "127 21 37"
    hdu.writeto(path, overwrite=True)


def test_run_warns_when_the_master_dark_mixes_epochs(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    # Lights on one night, with no darks of their own.
    for i in range(2):
        _write(raw / f"light{i}.fit", "Light Frame", 30.0, filt="V",
               dateobs="2025-04-30T15:10:00")
        _write(raw / f"flat{i}.fit", "Flat Field", 3.0, filt="V",
               dateobs="2025-04-30T10:00:00")
    # Darks from two epochs, 6 months apart, same exposure and temperature.
    _write(raw / "dark_a.fit", "Dark Frame", 30.0, dateobs="2025-04-28T18:39:00")
    _write(raw / "dark_b.fit", "Dark Frame", 30.0, dateobs="2024-11-06T18:39:00")

    frames = cs.scan_folder(str(raw), tz_offset_hours=9.0)
    lines: list[str] = []
    run_calibration(frames, ALL_NIGHTS, tmp_path / "out",
                    CalibrationOptions(cosmetic_enable=False),
                    log=lines.append)

    warned = [ln for ln in lines if "days" in ln and "apart" in ln]
    assert warned, f"no epoch warning in: {lines}"
    assert "master dark" in warned[0]
    assert "blend" in warned[0]


def test_run_is_quiet_when_all_frames_share_an_epoch(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    for i in range(2):
        _write(raw / f"light{i}.fit", "Light Frame", 30.0, filt="V",
               dateobs="2025-04-30T15:10:00")
        _write(raw / f"flat{i}.fit", "Flat Field", 3.0, filt="V",
               dateobs="2025-04-30T10:00:00")
        _write(raw / f"dark{i}.fit", "Dark Frame", 30.0,
               dateobs="2025-04-30T18:39:00")

    frames = cs.scan_folder(str(raw), tz_offset_hours=9.0)
    lines: list[str] = []
    run_calibration(frames, ALL_NIGHTS, tmp_path / "out",
                    CalibrationOptions(cosmetic_enable=False),
                    log=lines.append)
    assert not [ln for ln in lines if "apart" in ln], lines


def test_epoch_provenance_lands_in_the_summary(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    _write(raw / "light.fit", "Light Frame", 30.0, filt="V",
           dateobs="2025-04-30T15:10:00")
    _write(raw / "flat.fit", "Flat Field", 3.0, filt="V",
           dateobs="2025-04-30T10:00:00")
    _write(raw / "dark_a.fit", "Dark Frame", 30.0, dateobs="2025-04-28T18:39:00")
    _write(raw / "dark_b.fit", "Dark Frame", 30.0, dateobs="2024-11-06T18:39:00")

    frames = cs.scan_folder(str(raw), tz_offset_hours=9.0)
    summary = run_calibration(frames, ALL_NIGHTS, tmp_path / "out",
                              CalibrationOptions(cosmetic_enable=False))
    dark = next(m for m in summary["masters"] if m.get("type") == "dark")
    assert dark["epoch_nights"] == ["20241106", "20250428"]
    assert dark["epoch_span_days"] == 173
