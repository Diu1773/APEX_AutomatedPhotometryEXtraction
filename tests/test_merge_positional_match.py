"""Positional ID reconciliation in the multi-night merger (P5).

The old assignment was greedy in row order, so in a crowded field a distant row
could claim a canonical source that a much closer row needed, and that closer
row was written out as a spurious new star.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from apex.analysis.merge.id_match import (
    reconcile_workspace_catalogs,
    resolve_positional_pairs,
)


def test_closest_pair_wins_a_contested_source():
    """Row 0 is 1.8" away, row 1 is 0.3" away, both nearest the same source."""
    nn_sid = {0: 100, 1: 100}
    nn_sep = {0: 1.8, 1: 0.3}
    assigned = resolve_positional_pairs(nn_sid, nn_sep, tol_arcsec=2.0)
    assert assigned == {1: 100}          # the far row does not get to claim it


def test_assignment_is_independent_of_row_order():
    forward = resolve_positional_pairs({0: 7, 1: 7}, {0: 1.8, 1: 0.3}, 2.0)
    reverse = resolve_positional_pairs({1: 7, 0: 7}, {1: 0.3, 0: 1.8}, 2.0)
    assert forward == reverse == {1: 7}


def test_pairs_beyond_the_tolerance_are_dropped():
    assert resolve_positional_pairs({0: 5}, {0: 3.5}, tol_arcsec=2.0) == {}
    assert resolve_positional_pairs({0: 5}, {0: float("nan")}, 2.0) == {}
    assert resolve_positional_pairs({0: None}, {0: 0.1}, 2.0) == {}


def test_already_taken_sources_are_not_reassigned():
    assert resolve_positional_pairs({0: 9}, {0: 0.1}, 2.0, taken={9}) == {}


def test_one_source_per_row_and_one_row_per_source():
    nn_sid = {0: 1, 1: 1, 2: 2}
    nn_sep = {0: 0.5, 1: 0.2, 2: 0.4}
    assigned = resolve_positional_pairs(nn_sid, nn_sep, 2.0)
    assert assigned == {1: 1, 2: 2}
    assert len(set(assigned.values())) == len(assigned)


# --- end-to-end through the reconciler -------------------------------------

def _catalog(rows):
    return pd.DataFrame(
        [{"ID": i + 1, "source_id": sid, "ra_deg": ra, "dec_deg": dec}
         for i, (sid, ra, dec) in enumerate(rows)]
    )


def _offset_deg(arcsec):
    """Arcsec -> degrees, applied to DEC so the separation is exact (an RA
    offset would be compressed by cos(dec) and no longer equal the input)."""
    return arcsec / 3600.0


def test_crowded_field_does_not_invent_a_duplicate(tmp_path):
    """Two catalogue rows near one canonical star: the nearer one must match it
    and only the other becomes new — the old order-dependent pass could match
    the far one and create a duplicate for the near one."""
    base = tmp_path / "base"
    other = tmp_path / "other"
    canon_ra, canon_dec = 250.0, 36.0

    base_cat = _catalog([(None, canon_ra, canon_dec)])
    # Row order deliberately puts the far star first.
    other_cat = _catalog([
        (None, canon_ra, canon_dec + _offset_deg(1.8)),
        (None, canon_ra, canon_dec + _offset_deg(0.3)),
    ])

    result = reconcile_workspace_catalogs(
        [base, other],
        {str(base): {"V": base_cat}, str(other): {"V": other_cat}},
        {str(base): "F01", str(other): "F02"},
        pos_tol_arcsec=2.0,
    )

    records = [r for r in result["match_records"] if r["folder"] == other.name]
    by_local = {r["local_id"]: r for r in records}
    assert by_local[2]["method"] == "position"    # the 0.3" row matched
    assert by_local[2]["sep_arcsec"] == pytest.approx(0.3, abs=0.02)
    assert by_local[1]["method"] == "new"         # the 1.8" row is a new star
    # One canonical star existed; exactly one new one was added.
    assert len(result["canonical_by_filter"]["V"]) == 2


def test_source_id_match_is_never_stolen_by_a_positional_candidate(tmp_path):
    """An exact Gaia source_id match outranks any positional claim on it."""
    base = tmp_path / "base"
    other = tmp_path / "other"
    ra, dec = 250.0, 36.0

    base_cat = _catalog([(4242, ra, dec)])
    other_cat = _catalog([
        (None, ra, dec + _offset_deg(0.2)),       # very close, but no source_id
        (4242, ra, dec + _offset_deg(1.0)),       # same star by identity
    ])

    result = reconcile_workspace_catalogs(
        [base, other],
        {str(base): {"V": base_cat}, str(other): {"V": other_cat}},
        {str(base): "F01", str(other): "F02"},
        pos_tol_arcsec=2.0,
    )
    by_local = {r["local_id"]: r for r in result["match_records"]
                if r["folder"] == other.name}
    assert by_local[2]["method"] == "source_id"
    assert by_local[2]["merged_source_id"] == 4242
    assert by_local[1]["method"] == "new"         # the 0.2" row could not take it


def test_identical_catalogs_still_match_one_to_one(tmp_path):
    """Regression guard: the common case (same field, same stars) is unchanged."""
    base = tmp_path / "base"
    other = tmp_path / "other"
    rows = [(None, 250.0 + i * 0.01, 36.0) for i in range(5)]
    cat = _catalog(rows)

    result = reconcile_workspace_catalogs(
        [base, other],
        {str(base): {"V": cat.copy()}, str(other): {"V": cat.copy()}},
        {str(base): "F01", str(other): "F02"},
        pos_tol_arcsec=2.0,
    )
    summary = {(r["folder"], r["filter"]): r for r in result["match_summary_rows"]}
    assert summary[(other.name, "V")]["pos"] == 5
    assert summary[(other.name, "V")]["new"] == 0
    assert len(result["canonical_by_filter"]["V"]) == 5
