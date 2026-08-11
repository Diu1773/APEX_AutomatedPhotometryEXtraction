"""The one declaration of which Gaia columns APEX needs, and what they feed.

Two ADQL queries (ESA `gaiadr3.gaia_source`, VizieR `I/355/gaiadr3`) are
supposed to return the same catalogue, a normaliser renames the two vendors'
column names to one set, and `refbuild` decides which of those survive into the
master catalogue. Four hand-maintained lists, no check that they agreed — and
they did not:

* `ruwe` was in the VizieR SELECT and not the ESA one, so step10's RUWE <= 1.4
  cut ran or not depending on **which server answered**. Measured on M67: the
  calibrator count moved 564 vs 503 between two runs of identical code on
  identical data, with the zero point shifting only 7 mmag — small enough that
  nothing looked wrong.
* `phot_bp_rp_excess_factor` was in **neither** query, so the Riello et al.
  (2021) C* cut — which `gaia_quality.py` documents as "the published detector
  for exactly that failure mode" of BP/RP contamination in crowded fields — has
  never run on any target. APEX's targets include three globular clusters.
* `visibility_periods_used` was in refbuild's carry-through list and in neither
  query, so that entry was inert.

`gaia_quality_mask` is deliberately permissive about absent columns, which is
right — an old catalogue should still work — but it means a missing column
degrades the science silently instead of failing. The defence is not to be less
permissive; it is to stop the lists from drifting apart. Everything below is
derived from this table, and `tests/test_gaia_columns_contract.py` asserts that
the derived forms still agree.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

__all__ = [
    "GaiaColumn",
    "GAIA_COLUMNS",
    "esa_select_clause",
    "vizier_select_clause",
    "canonical_rename_map",
    "carry_through_columns",
    "missing_columns",
]

_VIZIER_TABLE = '"I/355/gaiadr3"'


@dataclass(frozen=True)
class GaiaColumn:
    """One logical column, and how each service spells it.

    `name` is the canonical form and deliberately equals ESA's own naming, so
    the ESA path needs no renaming and the VizieR path is the one that adapts.
    `aliases` are additional lower-case spellings the normaliser should accept,
    for catalogues written by earlier versions.
    """

    name: str
    esa: str | None
    vizier: str | None
    used_for: str
    carry: bool = False
    aliases: tuple[str, ...] = ()


# `vizier` entries are ADQL expressions, quoted where the column name contains
# characters ADQL will not take bare. Note the trap next to the excess factor:
# VizieR exposes BOTH "E(BP/RP)" (the flux excess factor, what C* needs) and
# "E(BP-RP)" (GSP-Phot reddening, a completely different quantity). The names
# differ by one character.
GAIA_COLUMNS: tuple[GaiaColumn, ...] = (
    GaiaColumn("source_id", "source_id", "Source",
               "cross-match key", aliases=("source",)),
    GaiaColumn("ra", "ra", "RA_ICRS", "astrometry", aliases=("ra_deg", "ra_icrs")),
    GaiaColumn("dec", "dec", "DE_ICRS", "astrometry", aliases=("dec_deg", "de_icrs")),
    GaiaColumn("phot_g_mean_mag", "phot_g_mean_mag", "Gmag",
               "magnitude limit, C* sigma(G)", aliases=("gmag",)),
    GaiaColumn("phot_bp_mean_mag", "phot_bp_mean_mag", "BPmag",
               "BP-RP colour, reference transformations", aliases=("bpmag",)),
    GaiaColumn("phot_rp_mean_mag", "phot_rp_mean_mag", "RPmag",
               "BP-RP colour, reference transformations", aliases=("rpmag",)),
    GaiaColumn("phot_variable_flag", "phot_variable_flag", None,
               "variable-star rejection (ESA only; VizieR has no equivalent)"),
    GaiaColumn("ruwe", "ruwe", "RUWE",
               "step10 calibrator quality cut RUWE <= 1.4", carry=True),
    GaiaColumn("phot_bp_rp_excess_factor", "phot_bp_rp_excess_factor",
               f'{_VIZIER_TABLE}."E(BP/RP)"',
               "step10 calibrator quality cut, Riello+2021 |C*| <= 3 sigma",
               carry=True, aliases=("e_bp_rp_", "e(bp/rp)")),
    GaiaColumn("visibility_periods_used", "visibility_periods_used", "Nper",
               "astrometric reliability (carried for diagnostics)", carry=True,
               aliases=("nper",)),
    GaiaColumn("parallax", "parallax", "Plx",
               "cluster membership, distance prior", carry=True, aliases=("plx",)),
    GaiaColumn("parallax_error", "parallax_error", "e_Plx",
               "cluster membership, distance prior", carry=True, aliases=("e_plx",)),
    GaiaColumn("pmra", "pmra", "pmRA", "cluster membership", carry=True),
    GaiaColumn("pmra_error", "pmra_error", "e_pmRA", "cluster membership",
               carry=True, aliases=("e_pmra",)),
    GaiaColumn("pmdec", "pmdec", "pmDE", "cluster membership", carry=True,
               aliases=("pmde",)),
    GaiaColumn("pmdec_error", "pmdec_error", "e_pmDE", "cluster membership",
               carry=True, aliases=("e_pmde",)),
)


def esa_select_clause() -> str:
    """Comma-separated SELECT list for `gaiadr3.gaia_source`."""
    return ", ".join(c.esa for c in GAIA_COLUMNS if c.esa)


def vizier_select_clause() -> str:
    """SELECT list for `I/355/gaiadr3`, aliased to the canonical names."""
    parts = []
    for col in GAIA_COLUMNS:
        if not col.vizier:
            continue
        expr = col.vizier if col.vizier.startswith('"') else f"{_VIZIER_TABLE}.{col.vizier}"
        parts.append(f"{expr} AS {col.name}")
    return ",\n  ".join(parts)


def canonical_rename_map() -> dict[str, str]:
    """Lower-case spelling -> canonical name, for every service and vintage."""
    mapping: dict[str, str] = {}
    for col in GAIA_COLUMNS:
        mapping[col.name.lower()] = col.name
        if col.vizier:
            bare = col.vizier.rsplit(".", 1)[-1].strip('"')
            mapping[bare.lower()] = col.name
        for alias in col.aliases:
            mapping[alias.lower()] = col.name
    return mapping


def carry_through_columns() -> tuple[str, ...]:
    """Columns `refbuild` must copy from the Gaia catalogue into the master."""
    return tuple(c.name for c in GAIA_COLUMNS if c.carry)


def missing_columns(df: "pd.DataFrame", *, only_used: bool = True) -> list[GaiaColumn]:
    """Contract columns absent from a catalogue, so a caller can say so.

    `only_used` skips ESA-only columns when reporting on a VizieR catalogue
    would just be noise — a caller that wants the full picture passes False.
    """
    present = {str(c).strip().lower() for c in df.columns}
    out = []
    for col in GAIA_COLUMNS:
        if col.name.lower() in present:
            continue
        if only_used and col.vizier is None and col.esa is not None:
            # ESA-only by nature; absent from VizieR is expected, not drift.
            continue
        out.append(col)
    return out
