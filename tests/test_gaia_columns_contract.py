"""Four lists of Gaia column names had to agree, and nobody checked.

The two ADQL SELECTs, the name normaliser and refbuild's carry-through list
each named the same columns separately, and they drifted:

* `ruwe` was in the VizieR SELECT and not the ESA one, so step10's RUWE cut ran
  or not depending on which server answered — M67's calibrator count came out
  564 one run and 503 the next, from identical code on identical data;
* `phot_bp_rp_excess_factor` was in neither SELECT, so the Riello+2021 C* cut
  that `gaia_quality.py` implements had never executed on any target;
* `visibility_periods_used` was in the carry-through list and neither SELECT,
  so that entry did nothing.

They are now all derived from `apex.utils.gaia_columns`. These tests are what
stops the derivation from being quietly bypassed again — they assert the shape
of what each consumer gets, not just that the module imports.
"""

from __future__ import annotations

import pandas as pd

from apex.analysis.refbuild import _GAIA_REF_EXTRA_COLS
from apex.utils.gaia_catalog_service import _normalise_gaia_columns
from apex.utils.gaia_columns import (
    GAIA_COLUMNS,
    canonical_rename_map,
    carry_through_columns,
    esa_select_clause,
    missing_columns,
    vizier_select_clause,
)


def test_both_services_request_the_same_logical_columns():
    """The drift that started it all: one SELECT knowing a column the other did not."""
    esa = {c.name for c in GAIA_COLUMNS if c.esa}
    vizier = {c.name for c in GAIA_COLUMNS if c.vizier}
    esa_only = esa - vizier
    assert esa_only == {"phot_variable_flag"}, (
        "the only column one service may have and the other may not is "
        f"phot_variable_flag (VizieR has no equivalent); found {esa_only}"
    )
    assert not vizier - esa, "VizieR must not request anything ESA does not"


def test_the_columns_the_quality_cuts_need_are_requested_from_both():
    """RUWE and the BP/RP excess factor drive step10's two calibrator cuts."""
    for name in ("ruwe", "phot_bp_rp_excess_factor"):
        col = next(c for c in GAIA_COLUMNS if c.name == name)
        assert col.esa, f"{name} missing from the ESA query"
        assert col.vizier, f"{name} missing from the VizieR query"
        assert col.carry, f"{name} must survive into the master catalogue"


def test_select_clauses_actually_contain_those_columns():
    """Guards the generation, not just the table."""
    esa, viz = esa_select_clause(), vizier_select_clause()
    for name in ("ruwe", "phot_bp_rp_excess_factor", "source_id", "parallax"):
        assert name in esa, f"{name} absent from the generated ESA SELECT"
        assert f"AS {name}" in viz, f"{name} absent from the generated VizieR SELECT"
    # The one-character trap: "E(BP/RP)" is the flux excess factor,
    # "E(BP-RP)" is GSP-Phot reddening. Selecting the latter would silently
    # feed the C* cut a completely different quantity.
    assert '"E(BP/RP)"' in viz
    assert "E(BP-RP)" not in viz


def test_refbuild_carries_exactly_what_the_contract_declares():
    assert tuple(_GAIA_REF_EXTRA_COLS) == carry_through_columns()
    assert "phot_bp_rp_excess_factor" in _GAIA_REF_EXTRA_COLS


def test_carry_list_never_names_a_column_no_query_fetches():
    """`visibility_periods_used` used to be carried and never fetched."""
    fetched = {c.name for c in GAIA_COLUMNS if c.esa or c.vizier}
    assert set(carry_through_columns()) <= fetched


def test_normaliser_maps_both_vendors_onto_canonical_names():
    vizier_style = pd.DataFrame({
        "Source": [1], "RA_ICRS": [10.0], "DE_ICRS": [20.0], "Gmag": [15.0],
        "BPmag": [15.5], "RPmag": [14.5], "RUWE": [1.0], "Plx": [1.0],
        "E(BP/RP)": [1.2],
    })
    out = _normalise_gaia_columns(vizier_style)
    for name in ("source_id", "ra", "dec", "phot_g_mean_mag", "ruwe",
                 "parallax", "phot_bp_rp_excess_factor"):
        assert name in out.columns, f"{name} not recovered from VizieR naming"

    esa_style = pd.DataFrame({
        "source_id": [1], "ra": [10.0], "dec": [20.0],
        "phot_g_mean_mag": [15.0], "ruwe": [1.0],
        "phot_bp_rp_excess_factor": [1.2],
    })
    out2 = _normalise_gaia_columns(esa_style)
    assert {"source_id", "ra", "dec", "phot_g_mean_mag", "ruwe",
            "phot_bp_rp_excess_factor"} <= set(out2.columns)


def test_every_alias_resolves_to_a_declared_column():
    names = {c.name for c in GAIA_COLUMNS}
    for alias, canonical in canonical_rename_map().items():
        assert canonical in names, f"alias {alias!r} points at unknown {canonical!r}"


def test_missing_columns_ignores_esa_only_absences():
    """A VizieR catalogue lacking phot_variable_flag is not stale.

    Treating it as stale is what made every VizieR-written cache miss, so each
    run re-fetched the catalogue and paid the ESA timeout again first.
    """
    vizier_cache = pd.DataFrame(columns=[
        c.name for c in GAIA_COLUMNS if c.vizier])
    assert missing_columns(vizier_cache, only_used=True) == []
    assert [c.name for c in missing_columns(vizier_cache, only_used=False)] == \
        ["phot_variable_flag"]


def test_missing_columns_reports_a_genuinely_stale_cache():
    old = pd.DataFrame(columns=["source_id", "ra", "dec", "phot_g_mean_mag"])
    missing = {c.name for c in missing_columns(old, only_used=True)}
    assert "ruwe" in missing
    assert "phot_bp_rp_excess_factor" in missing
