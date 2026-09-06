"""Unit tests for the Step 6 union master-catalog build (B1).

Exercises ``build_union_master`` orchestration (per-frame collection, sky-position
dedup via ``merge_ref_catalogs``, ``n_det_frames`` counting, the min-frames
filter, and dense re-numbering) without Qt or FITS fixtures, by passing stubs
for the per-frame catalog and WCS lookups.

Until 2026-09-06 these tests reached into ``RefBuildWorker`` and drove a copy of
this logic that no longer ran: the engine had moved to ``apex.analysis.refbuild``
and the Qt worker only forwards to it. The functions now live at module level so
the code under test is the code that runs.
"""

import numpy as np
import pandas as pd

from apex.analysis.refbuild import build_union_master


def _cat(rows):
    """rows: list of (ra, dec). Returns a per-frame catalog like the single-frame build."""
    df = pd.DataFrame(rows, columns=["ra_deg", "dec_deg"])
    df["x_ref"] = df["ra_deg"] * 1000.0
    df["y_ref"] = df["dec_deg"] * 1000.0
    df["source_id"] = np.arange(1, len(df) + 1, dtype=int)
    df["ID"] = df["source_id"]
    return df[["ID", "source_id", "ra_deg", "dec_deg", "x_ref", "y_ref"]]


def test_union_dedups_and_counts_detections():
    # Frame A (anchor): two stars. Frame B: one matches A's first, one new.
    cats = {
        "A": _cat([(10.0, 20.0), (10.0100, 20.0)]),   # ~34" apart -> distinct
        "B": _cat([(10.0, 20.0), (10.0200, 20.0)]),   # first matches A, second new
    }
    group = pd.DataFrame({"file": ["A", "B"]})
    master, stats = build_union_master(
        group, anchor_fname="A", match_radius_arcsec=2.0,
        build_master_catalog=lambda f: (cats[f].copy(), {"n_ref_used": len(cats[f])}),
    )

    # 3 unique stars: shared(10.0), A-only(10.01), B-only(10.02)
    assert len(master) == 3
    assert stats["n_union_frames"] == 2
    assert stats["n_master_union"] == 3

    # n_det_frames: shared star seen in 2 frames, the others in 1 each.
    counts = sorted(master["n_det_frames"].tolist())
    assert counts == [1, 1, 2]
    # The star at ra=10.0 is the one detected twice.
    shared = master.loc[np.isclose(master["ra_deg"], 10.0)]
    assert int(shared["n_det_frames"].iloc[0]) == 2

    # Dense, unique source_ids matching ID.
    assert sorted(master["source_id"].tolist()) == [1, 2, 3]
    assert master["source_id"].tolist() == master["ID"].tolist()


def test_union_min_frames_filter_drops_singletons():
    cats = {
        "A": _cat([(10.0, 20.0), (10.0100, 20.0)]),
        "B": _cat([(10.0, 20.0), (10.0200, 20.0)]),
    }
    group = pd.DataFrame({"file": ["A", "B"]})
    master, _ = build_union_master(
        group, anchor_fname="A", match_radius_arcsec=2.0,
        build_master_catalog=lambda f: (cats[f].copy(), {}),
        min_frames=2,
    )

    # Only the star detected in >= 2 frames survives.
    assert len(master) == 1
    assert np.isclose(master["ra_deg"].iloc[0], 10.0)
    assert int(master["n_det_frames"].iloc[0]) == 2
    assert master["source_id"].tolist() == [1]


def test_union_falls_back_when_all_frames_fail():
    # union build raises for every frame -> falls back to single-frame build,
    # which we stub to succeed for the anchor.
    anchor_cat = _cat([(10.0, 20.0)])
    calls = {"n": 0}

    def _build(f):
        calls["n"] += 1
        if calls["n"] == 1:
            # first call inside the loop (anchor) raises to force the union to
            # find nothing; subsequent fallback call returns a catalog.
            raise RuntimeError("missing detections")
        return anchor_cat.copy(), {"n_ref_used": 1}

    group = pd.DataFrame({"file": ["A"]})
    master, _ = build_union_master(
        group, anchor_fname="A", match_radius_arcsec=2.0,
        build_master_catalog=_build,
    )
    assert len(master) == 1


def test_union_stops_when_asked():
    """A stop request between frames leaves the frames already merged in place."""
    cats = {
        "A": _cat([(10.0, 20.0)]),
        "B": _cat([(10.0200, 20.0)]),
    }
    seen = []

    def _build(f):
        seen.append(f)
        return cats[f].copy(), {}

    group = pd.DataFrame({"file": ["A", "B"]})
    master, stats = build_union_master(
        group, anchor_fname="A", match_radius_arcsec=2.0,
        build_master_catalog=_build,
        should_stop=lambda: len(seen) >= 1,       # stop after the anchor
    )
    assert seen == ["A"]
    assert len(master) == 1
    assert stats["n_union_frames"] == 1
