"""Workspace scanning helpers for RESULT_/MERGED_ inputs."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from apex.utils.io_utils import coerce_int64_source_id, read_csv_int64_source_id
from apex.utils.run_workspace import (
    build_merged_workspace_dir,
    infer_workspace_date_range,
    infer_workspace_label,
    load_run_manifest,
)
from apex.utils.step_paths_lc import step7_forced_phot_dir, step8_selection_dir, step9_lc_dir


def normalize_filter_key(value) -> str:
    return str(value or "").strip().lower() or "unknown"


def folder_tag(index: int, folder: Path) -> str:
    safe = "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in folder.name)
    return f"F{index + 1:02d}_{safe}"


def default_merged_output_dir(folders: list[Path]) -> Path:
    return build_merged_workspace_dir(folders)


def run_meta(result_dir: Path) -> dict:
    meta = load_run_manifest(result_dir)
    if meta:
        return meta
    start_date, end_date = infer_workspace_date_range(result_dir)
    label = infer_workspace_label(result_dir)
    return {
        "run_type": "result",
        "label": label,
        "date_start": start_date,
        "date_end": end_date,
        "result_dir": str(result_dir),
    }


def read_step7_index(result_dir: Path) -> pd.DataFrame:
    idx_path = step7_forced_phot_dir(result_dir) / "photometry_index.csv"
    if not idx_path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(idx_path)
    except Exception:
        return pd.DataFrame()


def load_selection_payloads(result_dir: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    s9 = step8_selection_dir(result_dir)
    if not s9.exists():
        return out
    for path in sorted(s9.glob("selection_*.json")):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        flt = normalize_filter_key(data.get("filter") or path.stem.replace("selection_", ""))
        out[flt] = data
    return out


def load_master_catalogs_by_filter(result_dir: Path) -> dict[str, pd.DataFrame]:
    catalogs: dict[str, pd.DataFrame] = {}
    s9 = step8_selection_dir(result_dir)
    if not s9.exists():
        return catalogs
    for path in sorted(s9.glob("master_catalog_*.tsv")):
        flt = normalize_filter_key(path.stem.replace("master_catalog_", ""))
        try:
            df = read_csv_int64_source_id(path, sep="\t")
        except Exception:
            continue
        if df is None or df.empty or "ID" not in df.columns:
            continue
        df = df.copy()
        if "source_id" in df.columns:
            df["source_id"] = coerce_int64_source_id(df["source_id"]).astype("Int64")
        df["ID"] = pd.to_numeric(df["ID"], errors="coerce").astype("Int64")
        for col in ("ra_deg", "dec_deg", "gaia_G", "gaia_id"):
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        catalogs[flt] = df
    return catalogs


def workspace_scan_signature(result_dir: Path) -> tuple:
    """Cheap signature for cache invalidation during folder scan."""
    s7_idx = step7_forced_phot_dir(result_dir) / "photometry_index.csv"
    s9 = step8_selection_dir(result_dir)
    s10 = step9_lc_dir(result_dir)
    manifest = result_dir / "run_manifest.json"
    return (
        str(result_dir),
        manifest.stat().st_mtime if manifest.exists() else None,
        s7_idx.stat().st_mtime if s7_idx.exists() else None,
        max((p.stat().st_mtime for p in s9.glob("master_catalog_*.tsv")), default=None),
        max((p.stat().st_mtime for p in s9.glob("selection_*.json")), default=None),
        max((p.stat().st_mtime for p in s10.glob("lightcurve_*.csv")), default=None),
    )


def _scan_filter_names(result_dir: Path) -> list[str]:
    filters: set[str] = set()
    s9 = step8_selection_dir(result_dir)
    if s9.exists():
        for path in s9.glob("master_catalog_*.tsv"):
            filters.add(normalize_filter_key(path.stem.replace("master_catalog_", "")))
        for path in s9.glob("selection_*.json"):
            filters.add(normalize_filter_key(path.stem.replace("selection_", "")))
    return sorted(f for f in filters if f)


def scan_merge_input_workspace(result_dir: Path) -> dict:
    """Fast scan for Step 1 merger input validation without loading large TSVs."""
    meta = run_meta(result_dir)
    s7_idx = step7_forced_phot_dir(result_dir) / "photometry_index.csv"
    filters = _scan_filter_names(result_dir)
    s9_lc = step9_lc_dir(result_dir)
    has_step9_lc = any(s9_lc.glob("lightcurve_*.csv"))
    has_step7_forced = s7_idx.exists()
    has_master = any(step8_selection_dir(result_dir).glob("master_catalog_*.tsv"))
    has_sel = any(step8_selection_dir(result_dir).glob("selection_*.json"))
    has_step8_selection = bool(has_master and has_sel)
    merge_ready = bool(has_step7_forced and has_step8_selection and has_step9_lc)
    return {
        "folder": Path(result_dir),
        "label": str(meta.get("label") or Path(result_dir).name),
        "run_type": str(meta.get("run_type", "result")),
        "date_start": meta.get("date_start") or "—",
        "date_end": meta.get("date_end") or "—",
        "has_step7": has_step7_forced,
        "has_step8": has_step8_selection,
        "has_step9": has_step9_lc,
        "filters": filters,
        "merge_ready": merge_ready,
    }
