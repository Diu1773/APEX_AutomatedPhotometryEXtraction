"""Pure helpers for merged ID reconciliation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

from apex.utils.io_utils import coerce_int64_source_id, normalize_id_columns


def extract_row_float(row: pd.Series, *cols: str) -> float:
    for col in cols:
        if col in row.index:
            val = pd.to_numeric(pd.Series([row[col]]), errors="coerce").iloc[0]
            if np.isfinite(val):
                return float(val)
    return float("nan")


def row_radec(row: pd.Series) -> tuple[float, float]:
    return (
        extract_row_float(row, "ra_deg", "ra", "RA"),
        extract_row_float(row, "dec_deg", "dec", "DEC"),
    )


def append_folder_tag(existing: str, tag: str) -> str:
    parts = [p for p in str(existing or "").split(",") if p]
    if tag not in parts:
        parts.append(tag)
    return ",".join(parts)


def best_positional_match(row: pd.Series, canonical_df: pd.DataFrame, tol_arcsec: float) -> tuple[int | None, float]:
    ra, dec = row_radec(row)
    if not (np.isfinite(ra) and np.isfinite(dec)):
        return None, float("nan")
    if canonical_df is None or canonical_df.empty or "ra_deg" not in canonical_df.columns or "dec_deg" not in canonical_df.columns:
        return None, float("nan")

    cand = canonical_df.copy()
    cand_ra = pd.to_numeric(cand["ra_deg"], errors="coerce")
    cand_dec = pd.to_numeric(cand["dec_deg"], errors="coerce")
    mask = cand_ra.notna() & cand_dec.notna()
    if not mask.any():
        return None, float("nan")

    sc = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    csc = SkyCoord(cand_ra[mask].to_numpy(float) * u.deg, cand_dec[mask].to_numpy(float) * u.deg, frame="icrs")
    sep = sc.separation(csc).arcsec
    if len(sep) == 0:
        return None, float("nan")
    best_i = int(np.argmin(sep))
    best_sep = float(sep[best_i])
    if not np.isfinite(best_sep) or best_sep > tol_arcsec:
        return None, best_sep
    best_rows = cand.loc[mask].reset_index(drop=True)
    sid_val = pd.to_numeric(best_rows.loc[best_i, "source_id"], errors="coerce")
    if pd.isna(sid_val):
        return None, best_sep
    return int(sid_val), best_sep


def _pick_radec(frame: pd.DataFrame, names: tuple[str, ...]) -> pd.Series:
    for n in names:
        if n in frame.columns:
            return pd.to_numeric(frame[n], errors="coerce")
    return pd.Series(np.nan, index=frame.index, dtype="float64")


def resolve_positional_pairs(
    nn_sid: dict[int, int | None],
    nn_sep: dict[int, float],
    tol_arcsec: float,
    taken: set[int] | None = None,
) -> dict[int, int]:
    """Assign each incoming row its canonical partner, closest pair first.

    Row order used to decide who wins a contested canonical source: whichever
    row happened to be iterated first claimed it, so in a crowded field a row
    1.8" away could take the source that another row sat 0.3" from, and that
    second row — never offered a runner-up — was written out as a brand new
    star. Sorting the candidate pairs by separation makes the assignment
    order-independent and gives every source to its nearest claimant.

    Returns ``{row position: canonical source_id}`` for the pairs that are
    within ``tol_arcsec``; rows left out are new sources.
    """
    used: set[int] = set(taken or ())
    pairs = sorted(
        (
            (sep, pos, nn_sid[pos])
            for pos, sep in nn_sep.items()
            if nn_sid.get(pos) is not None and np.isfinite(sep) and sep <= tol_arcsec
        ),
        key=lambda item: (item[0], item[1]),
    )
    assigned: dict[int, int] = {}
    for _sep, pos, sid in pairs:
        if sid in used:
            continue
        used.add(int(sid))
        assigned[int(pos)] = int(sid)
    return assigned


def _nearest_canon_matches(
    df: pd.DataFrame, canon: pd.DataFrame
) -> tuple[dict[int, int | None], dict[int, float]]:
    """Nearest canonical source per incoming row, by sky position.

    Returns ``(pos -> source_id|None, pos -> separation_arcsec)`` keyed by the
    positional index of rows in ``df`` (iteration order). Vectorized equivalent
    of calling :func:`best_positional_match` per row — one
    ``match_to_catalog_sky`` over all rows instead of a fresh ``SkyCoord`` build
    per row. ``source_id`` is ``None`` when the nearest canonical row has no
    usable source_id (mirroring ``best_positional_match``).
    """
    nn_sid: dict[int, int | None] = {}
    nn_sep: dict[int, float] = {}
    if (
        canon is None
        or canon.empty
        or "ra_deg" not in canon.columns
        or "dec_deg" not in canon.columns
        or "source_id" not in canon.columns
    ):
        return nn_sid, nn_sep

    canon_ra = pd.to_numeric(canon["ra_deg"], errors="coerce")
    canon_dec = pd.to_numeric(canon["dec_deg"], errors="coerce")
    cmask = (canon_ra.notna() & canon_dec.notna()).to_numpy()
    if not cmask.any():
        return nn_sid, nn_sep
    canon_sid = coerce_int64_source_id(canon["source_id"]).to_numpy()[cmask]
    ccat = SkyCoord(
        canon_ra.to_numpy(dtype=float)[cmask] * u.deg,
        canon_dec.to_numpy(dtype=float)[cmask] * u.deg,
        frame="icrs",
    )

    src_ra = _pick_radec(df, ("ra_deg", "ra", "RA"))
    src_dec = _pick_radec(df, ("dec_deg", "dec", "DEC"))
    smask = (src_ra.notna() & src_dec.notna()).to_numpy()
    if not smask.any():
        return nn_sid, nn_sep
    scat = SkyCoord(
        src_ra.to_numpy(dtype=float)[smask] * u.deg,
        src_dec.to_numpy(dtype=float)[smask] * u.deg,
        frame="icrs",
    )

    match_idx, sep2d, _ = scat.match_to_catalog_sky(ccat)
    seps = sep2d.arcsec
    positions = np.flatnonzero(smask)
    for k, pos in enumerate(positions):
        sid_val = canon_sid[int(match_idx[k])]
        nn_sid[int(pos)] = None if pd.isna(sid_val) else int(sid_val)
        nn_sep[int(pos)] = float(seps[k])
    return nn_sid, nn_sep


#: Columns that hold a Gaia DR3 identifier, in the order they are trusted.
#: ``gaia_source_id`` is written by Step 6; ``gaia_id`` is what a previously
#: merged workspace carries.
GAIA_ID_COLUMNS = ("gaia_source_id", "gaia_id")


def global_identity_series(df: pd.DataFrame) -> pd.Series:
    """Return the per-row identifier that means the same thing in every folder.

    ``source_id`` does not qualify. Step 6 puts the Gaia DR3 identifier there
    only when the source matched Gaia; otherwise it is a counter that means
    nothing outside its own workspace. Keying the merge on it made two folders
    pointing at opposite halves of the sky bind row 5 to row 5 without ever
    measuring a separation (``Main/FAILURES.md`` F-285). The Gaia columns are
    the only globally valid key, and Step 6 writes them in both build modes.
    """
    empty = pd.Series(pd.array([pd.NA] * len(df), dtype="Int64"), index=df.index)
    if df is None or df.empty:
        return empty
    out = empty
    for col in GAIA_ID_COLUMNS:
        if col not in df.columns:
            continue
        vals = coerce_int64_source_id(df[col]).astype("Int64")
        vals = vals.where(vals.notna() & (vals > 0))   # a Gaia DR3 id is positive
        out = out.where(out.notna(), vals)
    return out


def canonicalize_catalog_row(
    row: pd.Series,
    merged_id: int,
    merged_source_id: int,
    folder_tag: str,
) -> dict:
    data = row.to_dict()
    data["ID"] = int(merged_id)
    data["source_id"] = int(merged_source_id)
    # pd.NA, not np.nan: a dict column mixing ints with np.nan becomes float64
    # when the rows are turned into a DataFrame, and a 19-digit Gaia id does not
    # survive that. pd.NA keeps the column object-typed and the digits intact.
    data["gaia_id"] = int(merged_source_id) if int(merged_source_id) > 0 else pd.NA
    # Rewrite gaia_source_id from the merged id rather than trusting what the
    # row carried. `df.iterrows()` casts an all-numeric row to one dtype, so the
    # value in `row` may already be a rounded float; `merged_source_id` came
    # from the DataFrame column and is exact.
    data["gaia_source_id"] = data["gaia_id"]
    data["match_status"] = "matched" if int(merged_source_id) > 0 else "no_gaia_match"
    data["folder_count"] = 1
    data["folder_tags"] = folder_tag
    return data


def next_generated_negative_source_id(current_catalogs: dict[str, pd.DataFrame]) -> int:
    min_sid = 0
    for df in current_catalogs.values():
        if df is None or df.empty or "source_id" not in df.columns:
            continue
        sid_vals = coerce_int64_source_id(df["source_id"]).dropna().astype("int64")
        if not sid_vals.empty:
            min_sid = min(min_sid, int(sid_vals.min()))
    return min_sid - 1 if min_sid <= 0 else -1


def reconcile_workspace_catalogs(
    folders: list[Path],
    catalogs_by_folder: dict[str, dict[str, pd.DataFrame]],
    folder_tags: dict[str, str],
    pos_tol_arcsec: float,
    logger=None,
) -> dict:
    if not folders:
        return {
            "canonical_by_filter": {},
            "local_id_maps": {},
            "match_summary_rows": [],
            "match_records": [],
        }

    all_filters = sorted({
        flt for folder in folders
        for flt in catalogs_by_folder.get(str(folder), {}).keys()
    })
    if not all_filters:
        return {
            "canonical_by_filter": {},
            "local_id_maps": {},
            "match_summary_rows": [],
            "match_records": [],
        }

    base_folder = folders[0]
    next_negative_sid = next_generated_negative_source_id({})
    canonical_by_filter: dict[str, pd.DataFrame] = {}
    next_id_by_filter: dict[str, int] = {}
    local_id_maps: dict[str, dict[str, dict[int, dict[str, int]]]] = {}
    match_summary_rows: list[dict] = []
    match_records: list[dict] = []

    for folder in folders:
        folder_key = str(folder)
        folder_tag = folder_tags[folder_key]
        local_id_maps.setdefault(folder_key, {})
        filter_catalogs = catalogs_by_folder.get(folder_key, {})

        for flt in all_filters:
            df = filter_catalogs.get(flt)
            if df is None or df.empty:
                continue
            # Work on a copy whose identifier columns are Int64. A caller that
            # built the table itself may have string ids; this makes those
            # exact. It cannot undo a float64 column — those digits are already
            # gone before the reconciler sees them.
            df = normalize_id_columns(df.copy())

            if flt not in canonical_by_filter:
                canonical_by_filter[flt] = pd.DataFrame()
                next_id_by_filter[flt] = 1

            canon = canonical_by_filter[flt].copy()
            local_map: dict[int, dict[str, int]] = {}
            n_exact = 0
            n_pos = 0
            n_new = 0

            if folder == base_folder and canon.empty:
                seeded_rows = []
                max_id = 0
                base_gaia = global_identity_series(df)
                for row_pos, (_, row) in enumerate(df.iterrows()):
                    local_id = pd.to_numeric(pd.Series([row.get("ID")]), errors="coerce").iloc[0]
                    if not np.isfinite(local_id):
                        continue
                    # The canonical id is the Gaia one when the source has it;
                    # a workspace-local counter must not survive into the merged
                    # catalogue, where a positive id is read back as a Gaia id.
                    sid_val = base_gaia.iloc[row_pos]
                    if pd.isna(sid_val):
                        sid = next_negative_sid
                        next_negative_sid -= 1
                    else:
                        sid = int(sid_val)
                    merged_id = int(local_id)
                    max_id = max(max_id, merged_id)
                    seeded_rows.append(canonicalize_catalog_row(row, merged_id, sid, folder_tag))
                    local_map[int(local_id)] = {
                        "merged_id": merged_id,
                        "merged_source_id": sid,
                    }
                    match_records.append({
                        "folder": folder.name,
                        "folder_tag": folder_tag,
                        "filter": flt,
                        "local_id": int(local_id),
                        "local_source_id": None if pd.isna(sid_val) else int(sid_val),
                        "merged_id": merged_id,
                        "merged_source_id": sid,
                        "method": "base",
                        "sep_arcsec": np.nan,
                        "status": "base",
                    })
                canon = pd.DataFrame(seeded_rows)
                next_id_by_filter[flt] = max_id + 1 if max_id > 0 else 1
                canonical_by_filter[flt] = canon
                local_id_maps[folder_key][flt] = local_map
                match_summary_rows.append({
                    "folder": folder.name,
                    "filter": flt,
                    "exact": len(seeded_rows),
                    "pos": 0,
                    "new": 0,
                    "total": len(seeded_rows),
                    "status": "base",
                })
                continue

            canon_sid_map = {}
            if not canon.empty and "source_id" in canon.columns:
                sid_vals = coerce_int64_source_id(canon["source_id"]).astype("Int64")
                for idx_row, sid_val in enumerate(sid_vals):
                    if pd.notna(sid_val) and int(sid_val) not in canon_sid_map:
                        canon_sid_map[int(sid_val)] = idx_row

            # Identity across folders runs on the Gaia id, never on source_id.
            canon_gaia_map: dict[int, int] = {}
            for idx_row, gid in enumerate(global_identity_series(canon)):
                if pd.notna(gid) and int(gid) not in canon_gaia_map:
                    canon_gaia_map[int(gid)] = idx_row

            used_canonical_sids: set[int] = set()
            # Collect new rows in a list; concat once at the end to avoid O(N²) copies
            new_canon_rows: list[dict] = []

            # `canon` is fixed during this loop (only metadata columns mutate),
            # so positionally match every incoming row against it in a single
            # vectorized pass instead of one SkyCoord build per row.
            pos_nn_sid, pos_nn_sep = _nearest_canon_matches(df, canon)
            # Rows that already match a canonical source by Gaia id take it
            # before the positional pass runs, so a positional candidate can
            # never steal a source out from under an exact identity match.
            df_gaia = global_identity_series(df)
            canon_sids = (coerce_int64_source_id(canon["source_id"]).astype("Int64")
                          if not canon.empty and "source_id" in canon.columns
                          else pd.Series(dtype="Int64"))
            gaia_to_canon_sid: dict[int, int] = {}
            for gid in df_gaia.dropna():
                idx_row = canon_gaia_map.get(int(gid))
                if idx_row is None or idx_row >= len(canon_sids):
                    continue
                csid = canon_sids.iloc[idx_row]
                if pd.notna(csid):
                    gaia_to_canon_sid[int(gid)] = int(csid)
            exact_sids = set(gaia_to_canon_sid.values())
            pos_assigned = resolve_positional_pairs(
                pos_nn_sid, pos_nn_sep, pos_tol_arcsec, taken=exact_sids)

            for pos, (_, row) in enumerate(df.iterrows()):
                local_id = pd.to_numeric(pd.Series([row.get("ID")]), errors="coerce").iloc[0]
                if not np.isfinite(local_id):
                    continue
                local_id = int(local_id)
                sid_val = coerce_int64_source_id(pd.Series([row.get("source_id")])).iloc[0]
                sid_int = None if pd.isna(sid_val) else int(sid_val)
                gid_val = df_gaia.iloc[pos]
                gid_int = None if pd.isna(gid_val) else int(gid_val)

                matched_sid = None
                match_method = ""
                sep_arcsec = float("nan")

                if gid_int is not None and gid_int in gaia_to_canon_sid:
                    matched_sid = gaia_to_canon_sid[gid_int]
                    match_method = "gaia_id"
                else:
                    sep_arcsec = pos_nn_sep.get(pos, float("nan"))
                    cand_sid = pos_assigned.get(pos)
                    if cand_sid is not None and cand_sid not in used_canonical_sids:
                        matched_sid = cand_sid
                        match_method = "position"
                    else:
                        matched_sid = None

                if matched_sid is not None and matched_sid in canon_sid_map:
                    canon_idx = canon_sid_map[matched_sid]
                    merged_id = int(pd.to_numeric(pd.Series([canon.iloc[canon_idx]["ID"]]), errors="coerce").iloc[0])
                    local_map[local_id] = {
                        "merged_id": merged_id,
                        "merged_source_id": int(matched_sid),
                    }
                    used_canonical_sids.add(int(matched_sid))
                    if match_method == "gaia_id":
                        n_exact += 1
                    else:
                        n_pos += 1
                    canon.at[canon_idx, "folder_count"] = int(pd.to_numeric(pd.Series([canon.iloc[canon_idx].get("folder_count", 1)]), errors="coerce").iloc[0] or 1) + 1
                    canon.at[canon_idx, "folder_tags"] = append_folder_tag(canon.iloc[canon_idx].get("folder_tags", ""), folder_tag)
                    match_records.append({
                        "folder": folder.name,
                        "folder_tag": folder_tag,
                        "filter": flt,
                        "local_id": local_id,
                        "local_source_id": sid_int,
                        "merged_id": merged_id,
                        "merged_source_id": int(matched_sid),
                        "method": match_method,
                        "sep_arcsec": sep_arcsec,
                        "status": "matched",
                    })
                    continue

                merged_id = next_id_by_filter.get(flt, 1)
                next_id_by_filter[flt] = merged_id + 1

                # A new canonical row keeps the Gaia id when the source has
                # one. Without it the row gets a fresh negative id, because a
                # positive id in the merged catalogue is read back as Gaia's.
                if gid_int is not None and gid_int not in canon_sid_map:
                    merged_source_id = gid_int
                else:
                    merged_source_id = next_negative_sid
                    next_negative_sid -= 1

                new_row = canonicalize_catalog_row(row, merged_id, merged_source_id, folder_tag)
                # Defer concat: track position in final merged df (len(canon) + pending rows)
                canon_sid_map[int(merged_source_id)] = len(canon) + len(new_canon_rows)
                new_canon_rows.append(new_row)
                local_map[local_id] = {
                    "merged_id": merged_id,
                    "merged_source_id": int(merged_source_id),
                }
                n_new += 1
                match_records.append({
                    "folder": folder.name,
                    "folder_tag": folder_tag,
                    "filter": flt,
                    "local_id": local_id,
                    "local_source_id": sid_int,
                    "merged_id": merged_id,
                    "merged_source_id": int(merged_source_id),
                    "method": "new",
                    "sep_arcsec": np.nan,
                    "status": "new",
                })

            # Single concat for all new rows in this folder/filter
            if new_canon_rows:
                canon = pd.concat(
                    [canon, pd.DataFrame(new_canon_rows)],
                    ignore_index=True, sort=False,
                )

            canon = canon.sort_values("ID").reset_index(drop=True)
            canonical_by_filter[flt] = canon
            local_id_maps[folder_key][flt] = local_map
            match_summary_rows.append({
                "folder": folder.name,
                "filter": flt,
                "exact": n_exact,
                "pos": n_pos,
                "new": n_new,
                "total": len(local_map),
                "status": "OK" if local_map else "empty",
            })
            if logger is not None:
                logger(
                    f"[MATCH] {folder.name} / {flt}: exact={n_exact} positional={n_pos} new={n_new} total={len(local_map)}"
                )

    return {
        "canonical_by_filter": canonical_by_filter,
        "local_id_maps": local_id_maps,
        "match_summary_rows": match_summary_rows,
        "match_records": match_records,
    }
