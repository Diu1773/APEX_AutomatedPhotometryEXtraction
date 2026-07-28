"""Neighbour distances must ignore single-frame detections.

Why this exists (2026-07-28). ``neighbor_dist_px`` / ``crowding_flag`` are
computed from the master catalogue against itself. Sources detected in exactly
one frame are mostly cosmic rays and hot pixels, so including them mislabels the
*real* star beside them as crowded -- and unlike a bad photometric row, that
label is not filtered downstream: it is baked into the master and travels with
the real star. Measured on ten frames each at ``ref_fwhm * 2.5``:
6/732 real stars affected in M67 g', 12/830 in M3 B.

These pin the two halves that matter: the tree excludes single-frame rows, and
single-frame rows still get a distance instead of being dropped.
"""

import numpy as np
import pandas as pd
import pytest

from apex.analysis.refbuild import neighbor_distances

# One real star at the origin, a cosmic ray 5 px away, and two more real stars
# far off. Without the filter the real star's nearest neighbour is the cosmic
# ray; with it, the nearest real neighbour.
_XY = np.array([[0.0, 0.0], [100.0, 0.0], [5.0, 0.0], [200.0, 0.0]])
_ND = [10, 10, 1, 10]


def test_single_frame_source_does_not_crowd_a_real_star():
    dist, n_trusted = neighbor_distances(_XY, _ND)
    assert n_trusted == 3
    assert dist[0] == pytest.approx(100.0), "cosmic ray still counted as neighbour"


def test_single_frame_source_still_gets_a_distance():
    """Excluded from the tree, not from the answer."""
    dist, _ = neighbor_distances(_XY, _ND)
    assert dist[2] == pytest.approx(5.0)
    assert np.isfinite(dist).all()
    assert len(dist) == len(_XY)


def test_without_n_det_frames_behaviour_is_unchanged():
    """The single-anchor path (ref_master_union=False) has no such column."""
    dist, n_trusted = neighbor_distances(_XY, None)
    assert n_trusted == len(_XY)
    assert dist[0] == pytest.approx(5.0)


def test_filter_never_empties_the_tree():
    """A run where nothing is seen twice falls back to unfiltered behaviour."""
    dist, n_trusted = neighbor_distances(_XY, [1, 1, 1, 1])
    assert n_trusted == len(_XY)
    assert dist[0] == pytest.approx(5.0)


def test_accepts_a_pandas_column():
    """The caller passes a DataFrame column, not a list."""
    dist, n_trusted = neighbor_distances(_XY, pd.Series(_ND, name="n_det_frames"))
    assert n_trusted == 3
    assert dist[0] == pytest.approx(100.0)


def test_non_numeric_n_det_frames_does_not_crash():
    dist, n_trusted = neighbor_distances(_XY, ["10", "10", None, "10"])
    assert n_trusted == 3
    assert np.isfinite(dist).all()
