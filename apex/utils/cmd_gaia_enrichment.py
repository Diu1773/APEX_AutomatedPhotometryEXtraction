"""Shared Gaia enrichment helpers for CMD workflow tables."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from apex.utils.io_utils import parse_int64_series, read_ecsv_int64_source_id
from apex.utils.step_paths import step5_wcs_dir, step6_refbuild_dir


GAIA_ENRICH_COLS = (
    "phot_g_mean_mag",
    "phot_bp_mean_mag",
    "phot_rp_mean_mag",
    "pmra",
    "pmdec",
    "pmra_error",
    "pmdec_error",
    "parallax",
    "parallax_error",
    "ruwe",
    "visibility_periods_used",
    "gaia_pmem",
    "pmem_gaia",
    "membership_prob_gaia",
)


def first_existing_col(cols, candidates):
    colset = set(cols)
    for col in candidates:
        if col in colset:
            return col
    return None


def normalize_master_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    ren = {}
    id_col = first_existing_col(out.columns, ["ID", "id", "Id", "star_id", "master_id"])
    if id_col and id_col != "ID":
        ren[id_col] = "ID"
    sid_col = first_existing_col(
        out.columns, ["source_id", "gaia_source_id", "SOURCE_ID", "Source_ID"]
    )
    if sid_col and sid_col != "source_id":
        ren[sid_col] = "source_id"
    g_col = first_existing_col(
        out.columns, ["gaia_G", "gmag", "phot_g_mean_mag", "Gmag", "G_MAG"]
    )
    if g_col and g_col != "gaia_G":
        ren[g_col] = "gaia_G"
    bp_col = first_existing_col(
        out.columns, ["gaia_BP", "bpmag", "phot_bp_mean_mag", "BPmag", "BP_MAG"]
    )
    if bp_col and bp_col != "gaia_BP":
        ren[bp_col] = "gaia_BP"
    rp_col = first_existing_col(
        out.columns, ["gaia_RP", "rpmag", "phot_rp_mean_mag", "RPmag", "RP_MAG"]
    )
    if rp_col and rp_col != "gaia_RP":
        ren[rp_col] = "gaia_RP"
    if ren:
        out = out.rename(columns=ren)
    return out


def load_gaia_enrichment_table(result_dir: Path, needed_cols=None) -> pd.DataFrame | None:
    result_dir = Path(result_dir)
    needed = [c for c in (needed_cols or GAIA_ENRICH_COLS) if c != "source_id"]
    candidates = [
        step5_wcs_dir(result_dir) / "gaia_derived.csv",
        step5_wcs_dir(result_dir) / "gaia_fov.ecsv",
        result_dir / "gaia_derived.csv",
        result_dir / "gaia_fov.ecsv",
    ]
    best_df = None
    best_score = -1
    for path in candidates:
        if not path.exists():
            continue
        try:
            if path.suffix.lower() == ".ecsv":
                gdf = read_ecsv_int64_source_id(path)
            else:
                gdf = pd.read_csv(path, dtype={"source_id": str})
        except Exception:
            continue
        if gdf is None or gdf.empty or "source_id" not in gdf.columns:
            continue
        gdf = gdf.copy()
        gdf["source_id"] = parse_int64_series(gdf["source_id"]).astype("Int64")
        gdf = gdf[gdf["source_id"].notna()].copy()
        if gdf.empty:
            continue
        keep_cols = ["source_id"] + [c for c in GAIA_ENRICH_COLS if c in gdf.columns]
        score = sum(1 for c in needed if c in gdf.columns)
        if score > best_score:
            best_score = score
            best_df = gdf[keep_cols].drop_duplicates(subset=["source_id"], keep="first")
    return best_df


def merge_gaia_columns_from_catalog(
    df: pd.DataFrame,
    result_dir: Path,
    needed_cols=None,
) -> pd.DataFrame:
    if df is None or df.empty:
        return df
    requested_cols = [c for c in (needed_cols or GAIA_ENRICH_COLS) if c != "source_id"]
    gdf = load_gaia_enrichment_table(result_dir, requested_cols)
    if gdf is None or gdf.empty:
        return df

    add_cols = [c for c in (needed_cols or GAIA_ENRICH_COLS) if c in gdf.columns and c != "source_id"]
    if not add_cols:
        return df

    out = df.copy()
    join_keys = []
    if "gaia_source_id" in out.columns:
        join_keys.append("gaia_source_id")
    if "source_id" in out.columns:
        join_keys.append("source_id")
    if not join_keys:
        return out

    for left_key in join_keys:
        join_col = f"__gaia_join_{left_key}"
        out[join_col] = parse_int64_series(out[left_key]).astype("Int64")
        use = gdf[["source_id"] + add_cols].rename(columns={"source_id": join_col})
        before_cols = set(out.columns)
        try:
            merged = out.merge(use, on=join_col, how="left", suffixes=("", "__gaia"))
        except Exception:
            out = out.drop(columns=[join_col], errors="ignore")
            continue

        for col in add_cols:
            if col in before_cols:
                aux = f"{col}__gaia"
                if aux not in merged.columns:
                    continue
                missing = merged[col].isna()
                merged.loc[missing, col] = merged.loc[missing, aux]
                merged = merged.drop(columns=[aux])
            elif col in merged.columns:
                pass
        out = merged.drop(columns=[join_col], errors="ignore")

    return out


def load_master_table(result_dir: Path) -> tuple[pd.DataFrame, str, Path]:
    refbuild = step6_refbuild_dir(result_dir)
    per_filter = sorted(refbuild.glob("ref_catalog_*.tsv")) if refbuild.exists() else []
    candidates = [
        ("ref_catalog", refbuild / "ref_catalog.tsv", "\t"),
    ] + [
        ("ref_catalog", p, "\t") for p in per_filter
    ] + [
        ("master_catalog", refbuild / "master_catalog.tsv", "\t"),
        ("master_catalog", result_dir / "master_catalog.tsv", "\t"),
        ("master_gaia_map", result_dir / "master_gaia_map.csv", ","),
    ]
    tried = []
    for source_name, path, sep in candidates:
        if not path.exists():
            continue
        tried.append(str(path))
        try:
            df = pd.read_csv(path, sep=sep)
        except Exception:
            continue
        if df.empty:
            continue
        return normalize_master_columns(df), source_name, path
    msg = "master_catalog.tsv or master_gaia_map.csv missing/invalid"
    if tried:
        msg += f" (tried: {', '.join(tried)})"
    raise FileNotFoundError(msg)
