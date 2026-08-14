"""IRAF and APEX disagree about which pixel is the first one.

IRAF numbers the first pixel 1. APEX, numpy and photutils number it 0. Handing
APEX's coordinates to IRAF unchanged therefore asks it to measure one pixel down
and to the left of every star — 1.41 px away — and reading its answers back
without the reverse shift leaves every returned position 1.41 px from the truth.

Both halves were missing until 2026-08-14, in the step-7 cross-check harness and
in the DAOPHOT/ALLSTAR driver, and neither produced an error. In the artificial-
star comparison the damage was severe: positions were matched to truth within
1.5 px, and the convention alone accounts for 1.414 of that, so ALLSTAR appeared
to lose a third of its stars (236 of 400 rather than 357) and the survivors were
whichever ones happened to fall inside the radius anyway.

Nothing about a one-pixel shift looks wrong in a log. These tests pin both
directions in both harnesses so the silence cannot come back.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from apex.benchmark.iraf_crosscheck import (
    IRAF_ORIGIN_OFFSET, parse_txdump, write_iraf_coords,
)


def test_offset_is_one_whole_pixel():
    assert IRAF_ORIGIN_OFFSET == 1.0


def test_written_coordinates_are_shifted_into_irafs_frame(tmp_path):
    reference = pd.DataFrame({"iraf_id": [1, 2], "x": [100.25, 7.0],
                              "y": [200.75, 9.5]})
    path = write_iraf_coords(reference, tmp_path / "stars.coo")
    written = pd.read_csv(path, sep=r"\s+", header=None, names=["x", "y"])
    assert written["x"].tolist() == pytest.approx([101.25, 8.0])
    assert written["y"].tolist() == pytest.approx([201.75, 10.5])


def test_returned_coordinates_are_shifted_back(tmp_path):
    dump = tmp_path / "txdump.txt"
    dump.write_text("1 101.25 201.75 12.0 0.01 100.0 1.0 50\n", encoding="ascii")
    table = parse_txdump(dump)
    assert table["iraf_x"].iloc[0] == pytest.approx(100.25)
    assert table["iraf_y"].iloc[0] == pytest.approx(200.75)


def test_round_trip_leaves_a_star_where_apex_put_it(tmp_path):
    """Write then read must be the identity, or `iraf_x - x` is meaningless."""
    reference = pd.DataFrame({"iraf_id": [1, 2, 3],
                              "x": [10.0, 512.5, 3999.125],
                              "y": [20.0, 256.25, 2048.875]})
    path = write_iraf_coords(reference, tmp_path / "stars.coo")
    written = pd.read_csv(path, sep=r"\s+", header=None, names=["x", "y"])
    dump = tmp_path / "txdump.txt"
    dump.write_text(
        "".join(f"{i + 1} {x} {y} 12.0 0.01 100.0 1.0 50\n"
                for i, (x, y) in enumerate(zip(written["x"], written["y"]))),
        encoding="ascii")
    back = parse_txdump(dump)
    assert back["iraf_x"].to_numpy() == pytest.approx(reference["x"].to_numpy())
    assert back["iraf_y"].to_numpy() == pytest.approx(reference["y"].to_numpy())


def test_daophot_driver_converts_both_ends():
    """The ALLSTAR driver carries its own copy of the same conversion."""
    from validation.psf_engines import daophot_allstar as driver

    assert driver.IRAF_ORIGIN_OFFSET == IRAF_ORIGIN_OFFSET
    source = __import__("inspect").getsource(driver.run)
    assert "+ IRAF_ORIGIN_OFFSET" in source, "입력 좌표 변환이 사라졌다"
    assert "- IRAF_ORIGIN_OFFSET" in source, "출력 좌표 변환이 사라졌다"


def test_an_unconverted_match_would_fall_at_the_radius_used():
    """Why this mattered: the error is the size of the matching radius.

    The comparison scripts pair an engine's star to an implanted one within
    1.5 px. A whole-pixel shift in each axis puts every star 1.414 px away — so
    the convention alone nearly exhausts the tolerance, and which stars survive
    becomes a matter of which way their fit happened to wander.
    """
    assert np.hypot(IRAF_ORIGIN_OFFSET, IRAF_ORIGIN_OFFSET) == pytest.approx(1.4142, abs=1e-4)
    assert np.hypot(IRAF_ORIGIN_OFFSET, IRAF_ORIGIN_OFFSET) < 1.5
