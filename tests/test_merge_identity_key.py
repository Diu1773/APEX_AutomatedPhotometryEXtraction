"""Cross-folder identity must run on the Gaia id, never on ``source_id``.

Step 6 puts the Gaia DR3 identifier in ``source_id`` only when the source
matched Gaia. Without a Gaia catalogue — an offline observatory, a failed
query, a field Gaia does not cover — it leaves a counter that means nothing
outside its own workspace, and the merger used to accept that counter as an
identity. Two folders pointing at opposite halves of the sky then bound row 5
to row 5 without ever measuring a separation (``Main/FAILURES.md`` F-285).

Every test here fails on the pre-fix reconciler.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from apex.analysis.merge.id_match import (
    global_identity_series,
    reconcile_workspace_catalogs,
)

GAIA_BASE = 4295806720000000000        # a real DR3 id is ~19 digits
ARCSEC = 1.0 / 3600.0


def _catalog(rows):
    """rows: (source_id, gaia_source_id, ra_deg, dec_deg)

    ``pd.NA``, never ``None``: a dict column mixing ints with ``None`` becomes
    float64 at construction and a 19-digit id loses its last three digits
    before anything under test has run.
    """
    return pd.DataFrame(
        [{"ID": i + 1, "source_id": sid,
          "gaia_source_id": pd.NA if gsid is None else gsid,
          "ra_deg": ra, "dec_deg": dec}
         for i, (sid, gsid, ra, dec) in enumerate(rows)]
    )


def _merge(tmp_path, cat_a, cat_b, tol=2.0):
    a, b = tmp_path / "night1", tmp_path / "night2"
    return reconcile_workspace_catalogs(
        [a, b],
        {str(a): {"V": cat_a}, str(b): {"V": cat_b}},
        {str(a): "N1", str(b): "N2"},
        pos_tol_arcsec=tol,
    ), b


# --- source_id must not be an identity ------------------------------------

@pytest.mark.parametrize("sid_sign", [1, -1])
def test_two_fields_sharing_local_numbers_stay_separate(tmp_path, sid_sign):
    """Opposite halves of the sky, same workspace counters, no Gaia ids."""
    a = _catalog([(sid_sign * (i + 1), None, 250.0 + i * 0.01, 36.0) for i in range(5)])
    b = _catalog([(sid_sign * (i + 1), None, 100.0 + i * 0.01, -20.0) for i in range(5)])

    res, folder_b = _merge(tmp_path, a, b)

    canon = res["canonical_by_filter"]["V"]
    assert len(canon) == 10, "ten distinct stars must stay ten rows"
    methods = {r["method"] for r in res["match_records"]
               if r["folder"] == folder_b.name}
    assert methods == {"new"}
    # Nothing may be reported as a Gaia match.
    assert canon["gaia_id"].notna().sum() == 0
    assert set(canon["match_status"]) == {"no_gaia_match"}


def test_local_counters_do_not_leak_into_the_merged_catalogue(tmp_path):
    """A merged row's positive source_id is read back as a Gaia id, so a
    workspace counter must be replaced by a fresh negative one."""
    a = _catalog([(7, None, 250.0, 36.0)])
    b = _catalog([(9, None, 100.0, -20.0)])

    res, _ = _merge(tmp_path, a, b)

    sids = res["canonical_by_filter"]["V"]["source_id"].astype("int64").tolist()
    assert all(s < 0 for s in sids), sids


# --- the Gaia id is the identity ------------------------------------------

def test_same_gaia_id_binds_across_folders(tmp_path):
    """Local numbering differs; the Gaia id decides. The non-Gaia star still
    matches positionally."""
    a = _catalog([
        (1, GAIA_BASE + 1, 250.00, 36.0),
        (2, GAIA_BASE + 2, 250.01, 36.0),
        (3, None,          250.02, 36.0),
    ])
    b = _catalog([
        (7, GAIA_BASE + 2, 250.01, 36.0),
        (8, GAIA_BASE + 1, 250.00, 36.0),
        (9, None,          250.02, 36.0),
    ])

    res, folder_b = _merge(tmp_path, a, b)

    canon = res["canonical_by_filter"]["V"]
    assert len(canon) == 3
    by_local = {r["local_id"]: r for r in res["match_records"]
                if r["folder"] == folder_b.name}
    assert by_local[1]["method"] == "gaia_id"
    assert by_local[1]["merged_source_id"] == GAIA_BASE + 2
    assert by_local[2]["method"] == "gaia_id"
    assert by_local[2]["merged_source_id"] == GAIA_BASE + 1
    assert by_local[3]["method"] == "position"


def test_gaia_id_outranks_the_positional_tolerance(tmp_path):
    """Same star, 5" apart, tolerance 2": identity still wins."""
    a = _catalog([(1, GAIA_BASE + 1, 250.0, 36.0)])
    b = _catalog([(5, GAIA_BASE + 1, 250.0, 36.0 + 5 * ARCSEC)])

    res, folder_b = _merge(tmp_path, a, b)

    assert len(res["canonical_by_filter"]["V"]) == 1
    rec = [r for r in res["match_records"] if r["folder"] == folder_b.name][0]
    assert rec["method"] == "gaia_id"


# --- 19-digit ids must survive the round trip ------------------------------

def test_gaia_ids_keep_every_digit_through_the_merge(tmp_path):
    """float64 runs out of mantissa at ~9e15. Two ids three apart must not
    collapse into one value on the way through the canonical DataFrame."""
    ids = [GAIA_BASE + 1, GAIA_BASE + 2, GAIA_BASE + 3]
    a = _catalog([(i + 1, gid, 250.0 + i * 0.01, 36.0) for i, gid in enumerate(ids)])
    b = _catalog([(i + 1, gid, 250.0 + i * 0.01, 36.0) for i, gid in enumerate(ids)])

    res, _ = _merge(tmp_path, a, b)

    canon = res["canonical_by_filter"]["V"]
    assert len(canon) == 3, "three ids must not collapse to one"
    assert sorted(canon["source_id"].astype("int64").tolist()) == ids
    assert sorted(canon["gaia_id"].astype("int64").tolist()) == ids
    assert sorted(canon["gaia_source_id"].astype("int64").tolist()) == ids


# --- the helper itself -----------------------------------------------------

def test_global_identity_prefers_gaia_columns_and_rejects_counters():
    df = pd.DataFrame([
        {"source_id": 1, "gaia_source_id": GAIA_BASE + 1},
        {"source_id": 2, "gaia_source_id": pd.NA},
        {"source_id": 3, "gaia_source_id": -4},        # a negative is not a Gaia id
    ])
    got = global_identity_series(df).tolist()
    assert got[0] == GAIA_BASE + 1
    assert pd.isna(got[1]) and pd.isna(got[2])


def test_global_identity_falls_back_to_the_gaia_id_column():
    df = pd.DataFrame([{"source_id": 5, "gaia_id": GAIA_BASE + 9}])
    assert global_identity_series(df).tolist() == [GAIA_BASE + 9]


def test_global_identity_without_any_gaia_column_is_all_missing():
    df = pd.DataFrame([{"source_id": 1}, {"source_id": 2}])
    assert global_identity_series(df).isna().all()
