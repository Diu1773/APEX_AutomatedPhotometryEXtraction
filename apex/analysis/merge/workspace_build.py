"""Helpers for materializing full merged workspaces."""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from apex.utils.io_utils import (
    coerce_int64_source_id,
    load_file_path_map,
    load_headers_table,
    load_night_assignments,
)
from apex.utils.photometry_loader import _load_source_to_id_map, load_frame_photometry
from apex.utils.run_workspace import infer_result_workspace_label, write_run_manifest
from apex.utils.step_paths_lc import (
    step1_dir,
    step7_forced_phot_dir,
    step8_selection_dir,
    step9_lc_dir,
    step10_detrend_dir,
    step11_period_dir,
)
from .workspace_scan import normalize_filter_key, read_step7_index


def selection_to_id_map(
    merged_catalogs: dict[str, pd.DataFrame],
    flt: str,
    source_ids: set[int],
) -> dict[int, int]:
    df = merged_catalogs.get(flt)
    if df is None or df.empty or "source_id" not in df.columns or "ID" not in df.columns:
        return {}
    sid_vals = coerce_int64_source_id(df["source_id"])
    id_vals = pd.to_numeric(df["ID"], errors="coerce").astype("Int64")
    out: dict[int, int] = {}
    for sid_val, id_val in zip(sid_vals, id_vals):
        if pd.isna(sid_val) or pd.isna(id_val):
            continue
        sid_int = int(sid_val)
        if sid_int in source_ids:
            out[sid_int] = int(id_val)
    return out


def write_selection_outputs(
    out_dir: Path,
    merged_catalogs: dict[str, pd.DataFrame],
    selection_target_by_filter: dict[str, int | None],
    selection_comp_by_filter: dict[str, set[int]],
    selection_check_by_filter: dict[str, int | None],
    stamp: str | None = None,
) -> None:
    s9 = step8_selection_dir(out_dir)
    s9.mkdir(parents=True, exist_ok=True)
    stamp = stamp or time.strftime("%Y-%m-%d %H:%M:%S")

    for flt, df in merged_catalogs.items():
        if df is None or df.empty:
            continue
        df_out = df.copy()
        target_sid = selection_target_by_filter.get(flt)
        comp_sids = selection_comp_by_filter.get(flt, set())
        check_sid = selection_check_by_filter.get(flt)

        def _role(sid: int) -> str:
            sid = int(sid)
            if target_sid is not None and sid == int(target_sid):
                return "T"
            if sid in comp_sids:
                return "C"
            if check_sid is not None and sid == int(check_sid):
                return "K"
            return ""

        sid_series = coerce_int64_source_id(df_out["source_id"]).fillna(-999999).astype("int64")
        df_out["role"] = [_role(int(sid)) for sid in sid_series]
        if "gaia_id" not in df_out.columns:
            df_out["gaia_id"] = [int(sid) if int(sid) > 0 else np.nan for sid in sid_series]
        if "match_status" not in df_out.columns:
            df_out["match_status"] = [
                "matched" if int(sid) > 0 else "no_gaia_match" for sid in sid_series
            ]

        output_cols = ["ID", "x_ref", "y_ref", "ra_deg", "dec_deg", "role", "gaia_id", "match_status"]
        for col in [
            "gaia_G", "gaia_BP", "gaia_RP",
            "gaia_g", "gaia_bp", "gaia_rp", "color_gr",
            "folder_count", "folder_tags", "source_id",
        ]:
            if col in df_out.columns and col not in output_cols:
                output_cols.append(col)
        for col in df_out.columns:
            if col not in output_cols:
                output_cols.append(col)
        df_out = df_out[[c for c in output_cols if c in df_out.columns]]
        df_out = df_out.sort_values("ID")

        cat_path = s9 / f"master_catalog_{flt}.tsv"
        df_out.to_csv(cat_path, sep="\t", index=False, na_rep="NaN", encoding="utf-8-sig")

        id_map_path = s9 / f"id_mapping_{flt}.csv"
        id_map_df = df_out[
            [c for c in ["ID", "source_id", "gaia_id", "role", "x_ref", "y_ref"] if c in df_out.columns]
        ].copy()
        id_map_df.to_csv(id_map_path, index=False, na_rep="NaN")

        sel_sids = set(int(s) for s in comp_sids if s is not None)
        if target_sid is not None:
            sel_sids.add(int(target_sid))
        if check_sid is not None:
            sel_sids.add(int(check_sid))
        sid_to_id = selection_to_id_map(merged_catalogs, flt, sel_sids)

        data = {
            "filter": flt,
            "target_id": sid_to_id.get(int(target_sid)) if target_sid is not None else None,
            "target_source_id": int(target_sid) if target_sid is not None else None,
            "comparison_ids": sorted(int(sid_to_id[sid]) for sid in comp_sids if sid in sid_to_id),
            "comparison_source_ids": sorted(int(sid) for sid in comp_sids),
            "check_id": sid_to_id.get(int(check_sid)) if check_sid is not None else None,
            "check_source_id": int(check_sid) if check_sid is not None else None,
            "timestamp": stamp,
        }
        (s9 / f"selection_{flt}.json").write_text(
            json.dumps(data, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )


def materialize_merged_workspace(
    out_dir: Path,
    folders: list[Path],
    folder_tags: dict[str, str],
    local_id_maps: dict[str, dict[str, dict[int, dict[str, int]]]],
    merged_catalogs: dict[str, pd.DataFrame],
    selection_target_by_filter: dict[str, int | None],
    selection_comp_by_filter: dict[str, set[int]],
    selection_check_by_filter: dict[str, int | None],
    match_records: list[dict],
) -> dict:
    s1 = step1_dir(out_dir)
    forced_phot_dir = step7_forced_phot_dir(out_dir)
    selection_dir = step8_selection_dir(out_dir)
    s1.mkdir(parents=True, exist_ok=True)
    forced_phot_dir.mkdir(parents=True, exist_ok=True)
    selection_dir.mkdir(parents=True, exist_ok=True)
    step9_lc_dir(out_dir).mkdir(parents=True, exist_ok=True)
    step10_detrend_dir(out_dir).mkdir(parents=True, exist_ok=True)
    step11_period_dir(out_dir).mkdir(parents=True, exist_ok=True)

    merged_headers_rows: list[dict] = []
    merged_index_rows: list[dict] = []
    merged_night_assignments: dict[str, int] = {}
    merged_path_map: dict[str, str] = {}
    merged_filter_frames: dict[str, list[str]] = {}

    next_merged_night = 1
    for folder in folders:
        folder_key = str(folder)
        folder_tag = folder_tags.get(folder_key, folder.name)
        source_path_map = load_file_path_map(folder)
        idx = read_step7_index(folder)
        if idx.empty or "file" not in idx.columns:
            continue

        headers_df = load_headers_table(folder)
        header_lookup: dict[str, dict] = {}
        if not headers_df.empty and "Filename" in headers_df.columns:
            # Warn on duplicates: set_index silently keeps last; log which files collide.
            dup_mask = headers_df["Filename"].duplicated(keep=False)
            if dup_mask.any():
                dup_names = sorted(headers_df.loc[dup_mask, "Filename"].astype(str).unique().tolist())
                import warnings
                warnings.warn(
                    f"[workspace_build] {folder.name}: {len(dup_names)} duplicate Filename(s) "
                    f"in headers.csv — keeping last occurrence. Duplicates: {dup_names[:5]}"
                    + (" …" if len(dup_names) > 5 else ""),
                    UserWarning,
                    stacklevel=4,
                )
            header_lookup = {str(fn): row.to_dict() for fn, row in headers_df.set_index("Filename").iterrows()}

        night_map_raw = load_night_assignments(folder)
        # to_dict("records") is much faster than iterrows() and yields plain
        # dicts reused for both the night-id pass and the main pass below.
        idx_records = idx.to_dict("records")
        local_night_ids: dict[str, int] = {}
        for row in idx_records:
            fname = str(row.get("file", "")).strip()
            if not fname:
                continue
            night_id = night_map_raw.get(fname)
            if night_id is None:
                nid_val = pd.to_numeric(pd.Series([row.get("night_id")]), errors="coerce").iloc[0]
                night_id = int(nid_val) if np.isfinite(nid_val) and int(nid_val) > 0 else 1
            local_night_ids[fname] = int(night_id)

        local_to_merged: dict[int, int] = {}
        for local_night in sorted(set(local_night_ids.values())):
            local_to_merged[local_night] = next_merged_night
            next_merged_night += 1

        # The Step 8 source_id→ID catalog is identical for every frame of a
        # filter within a folder; read it once per filter instead of letting
        # load_frame_photometry re-glob and re-parse it on each frame.
        sid_map_by_filter: dict[str, dict[int, int]] = {}
        # Per-filter {local_id: merged_*} lookups are constant across all frames
        # of a filter; build once and reuse instead of per-cell lambdas.
        id_lookup_by_filter: dict[str, dict[int, int]] = {}
        sourceid_lookup_by_filter: dict[str, dict[int, int]] = {}

        for row in idx_records:
            fname = str(row.get("file", "")).strip()
            if not fname:
                continue
            flt = normalize_filter_key(row.get("filter", row.get("FILTER", "")))
            local_map = local_id_maps.get(folder_key, {}).get(flt, {})
            if not local_map:
                continue
            if flt not in id_lookup_by_filter:
                id_lookup_by_filter[flt] = {
                    int(lid): int(m["merged_id"]) for lid, m in local_map.items()
                }
                sourceid_lookup_by_filter[flt] = {
                    int(lid): int(m["merged_source_id"]) for lid, m in local_map.items()
                }

            if flt not in sid_map_by_filter:
                sid_map_by_filter[flt] = _load_source_to_id_map(folder, flt)
            phot_df = load_frame_photometry(
                folder, fname, flt, sid_map=sid_map_by_filter[flt]
            )
            if phot_df is None or phot_df.empty or "ID" not in phot_df.columns:
                continue

            local_ids = pd.to_numeric(phot_df["ID"], errors="coerce").astype("Int64")
            phot_df = phot_df.loc[local_ids.notna()].copy()
            phot_df["ID_local"] = local_ids[local_ids.notna()].astype(int)
            phot_df = phot_df[phot_df["ID_local"].isin(local_map.keys())].copy()
            if phot_df.empty:
                continue

            phot_df["source_id"] = phot_df["ID_local"].map(sourceid_lookup_by_filter[flt])
            phot_df["ID"] = phot_df["ID_local"].map(id_lookup_by_filter[flt])

            merged_fname = f"{folder_tag}__{fname}"
            phot_df["file"] = merged_fname
            phot_df["source_folder"] = folder_tag
            phot_df["original_file"] = fname

            out_phot_path = forced_phot_dir / f"photometry_{merged_fname}.tsv"
            phot_df.to_csv(out_phot_path, sep="\t", index=False, na_rep="NaN")

            original_path = source_path_map.get(fname)
            if original_path:
                merged_path_map[merged_fname] = str(original_path)
            merged_night_id = local_to_merged.get(local_night_ids.get(fname, 1), 1)
            merged_night_assignments[merged_fname] = merged_night_id
            merged_filter_frames.setdefault(flt, []).append(merged_fname)

            row_dict = dict(row)
            row_dict["file"] = merged_fname
            row_dict["filter"] = flt
            row_dict["night_id"] = merged_night_id
            row_dict["source_folder"] = folder_tag
            row_dict["original_file"] = fname
            row_dict["path"] = str(out_phot_path)
            merged_index_rows.append(row_dict)

            header_row = header_lookup.get(fname, {}).copy()
            if not header_row:
                header_row = {"Filename": merged_fname, "FILTER": flt}
            else:
                header_row["Filename"] = merged_fname
            header_row["SourceFolder"] = folder_tag
            header_row["OriginalFilename"] = fname
            merged_headers_rows.append(header_row)

    if not merged_index_rows:
        raise RuntimeError("병합 가능한 forced photometry rows를 만들지 못했습니다.")

    pd.DataFrame(merged_index_rows).to_csv(forced_phot_dir / "photometry_index.csv", index=False)
    if merged_headers_rows:
        pd.DataFrame(merged_headers_rows).drop_duplicates(
            subset=["Filename"], keep="first"
        ).to_csv(s1 / "headers.csv", index=False)
    (s1 / "night_assignments.json").write_text(
        json.dumps({"night_assignments": merged_night_assignments}, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    (s1 / "file_path_map.json").write_text(
        json.dumps(merged_path_map, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    (forced_phot_dir / "filter_frames.json").write_text(
        json.dumps(merged_filter_frames, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    master_source_rows: list[dict] = []
    seen_master_sids: set[int] = set()
    for flt, df in merged_catalogs.items():
        if df is None or df.empty or "source_id" not in df.columns:
            continue
        sid_vals = coerce_int64_source_id(df["source_id"])
        for idx, sid_val in sid_vals.items():
            if pd.isna(sid_val):
                continue
            sid_int = int(sid_val)
            if sid_int in seen_master_sids:
                continue
            seen_master_sids.add(sid_int)
            row = {"source_id": sid_int}
            for col in (
                "n_frames",
                "ra_deg",
                "dec_deg",
                "x_ref",
                "y_ref",
                "x_mean",
                "y_mean",
                "x_std",
                "y_std",
                "gaia_G",
                "gaia_BP",
                "gaia_RP",
                "gaia_g",
                "gaia_bp",
                "gaia_rp",
                "match_status",
            ):
                if col in df.columns:
                    row[col] = df.at[idx, col]
            master_source_rows.append(row)
    if master_source_rows:
        pd.DataFrame(master_source_rows).to_csv(
            forced_phot_dir / "master_sources.csv", index=False, na_rep="NaN"
        )

    write_selection_outputs(
        out_dir=out_dir,
        merged_catalogs=merged_catalogs,
        selection_target_by_filter=selection_target_by_filter,
        selection_comp_by_filter=selection_comp_by_filter,
        selection_check_by_filter=selection_check_by_filter,
    )

    manifest = {
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "input_folders": [str(p) for p in folders],
        "folder_tags": folder_tags,
        "filters": sorted(merged_catalogs.keys()),
        "merged_result_dir": str(out_dir),
        "records": len(match_records),
    }
    (out_dir / "merge_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    pd.DataFrame(match_records).to_csv(out_dir / "merge_id_map.csv", index=False)
    write_run_manifest(
        out_dir,
        run_type="merged",
        root_dir=out_dir.parent,
        input_result_dirs=folders,
        target_name=infer_result_workspace_label(folders),
        storage_mode="full",
    )

    return {
        "night_assignments": merged_night_assignments,
        "path_map": merged_path_map,
        "index_rows": merged_index_rows,
    }
