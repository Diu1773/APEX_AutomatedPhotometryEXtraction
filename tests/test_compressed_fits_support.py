"""A frame from a public archive must read the same as one off our own camera.

Archives serve frames as ``.fits.fz`` (fpack), which puts an **empty**
``PrimaryHDU`` at index 0 and the image in a ``CompImageHDU`` at index 1.
APEX listed ``.fits.fz`` among the file types it accepts, but every reader
took ``hdul[0]``, so every keyword lookup silently returned its default and
Step 0 reported "no light frames found".

The second half is naming: LCO writes ``OBSTYPE = EXPOSE`` for a science
frame, which matched none of LIGHT / SCIENCE / OBJECT.

Both were found by actually running APEX on LCO MuSCAT3 data
(``Main/FAILURES.md`` F-302).
"""

from __future__ import annotations

import numpy as np
import pytest
from astropy.io import fits

from apex.analysis.calibration_scan import classify_type, read_frame_info
from apex.utils.io_utils import (
    read_fits_header,
    read_fits_image,
    science_hdu_index,
)


def _write_compressed(path, obstype="EXPOSE", exptime=10.0, filt="gp"):
    """A file shaped like the archive's: empty primary, image in extension 1."""
    data = np.arange(64, dtype=np.int16).reshape(8, 8)
    hdr = fits.Header()
    hdr["OBSTYPE"] = obstype
    hdr["EXPTIME"] = exptime
    hdr["FILTER"] = filt
    hdr["DATE-OBS"] = "2021-03-18T08:58:44.051"
    hdr["GAIN"] = 1.9
    hdul = fits.HDUList([fits.PrimaryHDU(), fits.CompImageHDU(data, hdr, name="SCI")])
    hdul.writeto(path, overwrite=True)
    return path


def _write_plain(path, obstype="LIGHT"):
    data = np.arange(64, dtype=np.int16).reshape(8, 8)
    hdr = fits.Header()
    hdr["OBSTYPE"] = obstype
    hdr["EXPTIME"] = 10.0
    hdr["FILTER"] = "V"
    fits.PrimaryHDU(data, hdr).writeto(path, overwrite=True)
    return path


# --- the empty primary HDU ------------------------------------------------

def test_the_hazard_is_real_so_the_test_means_something(tmp_path):
    """Guard the premise: HDU 0 of a compressed file really is empty."""
    p = _write_compressed(tmp_path / "c.fits.fz")
    assert "OBSTYPE" not in fits.getheader(p)          # what APEX used to read
    assert fits.getheader(p).get("NAXIS") == 0


def test_science_hdu_is_found_past_an_empty_primary(tmp_path):
    p = _write_compressed(tmp_path / "c.fits.fz")
    with fits.open(p) as hdul:
        assert science_hdu_index(hdul) == 1


def test_a_plain_file_still_reads_from_hdu_zero(tmp_path):
    p = _write_plain(tmp_path / "p.fits")
    with fits.open(p) as hdul:
        assert science_hdu_index(hdul) == 0


@pytest.mark.parametrize("writer,name", [(_write_compressed, "c.fits.fz"),
                                         (_write_plain, "p.fits")])
def test_header_reads_the_same_either_way(tmp_path, writer, name):
    p = writer(tmp_path / name)
    h = read_fits_header(p)
    assert h.get("OBSTYPE")
    assert float(h.get("EXPTIME")) == 10.0
    assert h.get("FILTER")


@pytest.mark.parametrize("writer,name", [(_write_compressed, "c.fits.fz"),
                                         (_write_plain, "p.fits")])
def test_image_reads_the_same_either_way(tmp_path, writer, name):
    p = writer(tmp_path / name)
    data, header = read_fits_image(p, dtype=np.float32)
    assert data.shape == (8, 8)
    assert data.dtype == np.float32
    assert float(data[0, 0]) == 0.0
    assert float(data[-1, -1]) == 63.0
    assert header.get("OBSTYPE")


def test_read_fits_image_refuses_a_file_with_no_image(tmp_path):
    p = tmp_path / "empty.fits"
    fits.HDUList([fits.PrimaryHDU()]).writeto(p, overwrite=True)
    with pytest.raises(ValueError):
        read_fits_image(p)


# --- the OBSTYPE vocabulary -----------------------------------------------

@pytest.mark.parametrize("obstype", ["EXPOSE", "LIGHT", "SCIENCE", "OBJECT",
                                     "SCI", "TARGET"])
def test_science_frames_are_recognised_whatever_the_archive_calls_them(obstype):
    assert classify_type(obstype, "frame.fits") == "light"


@pytest.mark.parametrize("obstype,want", [("BIAS", "bias"), ("ZERO", "bias"),
                                          ("DARK", "dark"), ("SKYFLAT", "flat"),
                                          ("DOMEFLAT", "flat")])
def test_calibration_frames_keep_their_meaning(obstype, want):
    assert classify_type(obstype, "frame.fits") == want


def test_a_calibration_name_still_wins_over_the_science_words():
    """SKYFLAT contains no science word, but check the order holds anyway."""
    assert classify_type("SKYFLAT EXPOSE", "frame.fits") == "flat"


# --- the two together, through the real scanner ---------------------------

def test_the_scanner_classifies_an_archive_frame_as_light(tmp_path):
    """Both defects at once: compressed file whose OBSTYPE is EXPOSE."""
    p = _write_compressed(tmp_path / "ogg2m001-ep04-20210317-0043-e00.fits.fz")
    info = read_frame_info(str(p))
    assert info is not None
    assert info.ftype == "light", f"scanner said {info.ftype!r}"
    assert info.exp == pytest.approx(10.0)
    assert info.filt
    assert info.shape == (8, 8)          # HDU 0 에서 읽었으면 None 이다


def test_the_scanner_reads_a_compressed_bias_too(tmp_path):
    p = _write_compressed(tmp_path / "ogg2m001-ep04-20210317-0001-b00.fits.fz",
                          obstype="BIAS", exptime=0.0, filt="gp")
    info = read_frame_info(str(p))
    assert info is not None
    assert info.ftype == "bias"
