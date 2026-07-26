"""Shared I/O helpers for Step 10 check-star light-curve outputs."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from apex.utils.common_helpers import normalize_filter_key as _normalize_filter_key
from apex.utils.step_paths_lc import step8_selection_dir, step9_lc_dir


def load_check_star_meta_by_filter(result_dir: Path) -> dict[str, dict]:
    """Return ``{filter_key: {"check_id": int, "check_source_id": int}}`` from
    ``lc_selection/selection_*.json`` files.

    Only entries where at least one of check_id / check_source_id is present
    are included.
    """
    s9 = step8_selection_dir(result_dir)
    out: dict[str, dict] = {}
    if not s9.exists():
        return out
    for sel_path in sorted(s9.glob("selection_*.json")):
        raw_flt = sel_path.stem.replace("selection_", "")
        flt = _normalize_filter_key(raw_flt) or raw_flt
        try:
            data = json.loads(sel_path.read_text(encoding="utf-8"))
        except Exception:
            continue
        entry: dict = {}
        check_id = data.get("check_id")
        check_source_id = data.get("check_source_id")
        if check_id is not None:
            entry["check_id"] = int(check_id)
        if check_source_id is not None:
            entry["check_source_id"] = int(check_source_id)
        if entry:
            out[flt] = entry
    return out


def load_check_star_id(
    result_dir: Path, filt: str | None = None
) -> int | None:
    """Return the check-star integer ID, optionally for a specific filter."""
    meta = load_check_star_meta_by_filter(result_dir)
    if filt:
        entry = meta.get(_normalize_filter_key(filt), {})
        cid = entry.get("check_id")
        return int(cid) if cid is not None else None
    for entry in meta.values():
        cid = entry.get("check_id")
        if cid is not None:
            return int(cid)
    return None


def load_check_star_csv(
    result_dir: Path, filt: str | None = None
) -> tuple[int | None, pd.DataFrame]:
    """Load the check-star light curve CSV from ``lc_lightcurve/``.

    Parameters
    ----------
    result_dir : Path
        Pipeline result directory.
    filt : str, optional
        If given, filter the combined multi-dataset curve to that band before
        falling back to legacy per-filter files.

    Returns
    -------
    (check_id, DataFrame). The DataFrame is empty when no file is found. The
    ID is ``None`` when the combined curve contains multiple local check IDs.
    """
    out_dir = step9_lc_dir(result_dir)
    if not out_dir.exists():
        return None, pd.DataFrame()

    if filt:
        filt_key = _normalize_filter_key(filt)
        check_id = load_check_star_id(result_dir, filt_key)
        candidates: list[tuple[Path, bool]] = [
            (out_dir / "lightcurve_check_combined_raw.csv", False)
        ]
        if check_id is not None and filt_key:
            candidates.append(
                (out_dir / f"lightcurve_check_{filt_key}_ID{check_id}_raw.csv", True)
            )
            candidates.append((out_dir / f"lightcurve_check_ID{check_id}_raw.csv", True))
        for path, require_check_id in candidates:
            if not path.exists():
                continue
            try:
                df = pd.read_csv(path)
            except Exception:
                continue
            if "filter" in df.columns and filt_key:
                df = df[df["filter"].astype(str).map(_normalize_filter_key) == filt_key].copy()
            if require_check_id and "check_id" in df.columns and check_id is not None:
                df = df[pd.to_numeric(df["check_id"], errors="coerce") == int(check_id)].copy()
            if not df.empty:
                ids = []
                if "check_id" in df.columns:
                    ids = sorted({
                        int(value)
                        for value in pd.to_numeric(df["check_id"], errors="coerce")
                        .dropna()
                        .astype(int)
                        .tolist()
                    })
                loaded_id = ids[0] if len(ids) == 1 else (check_id if not ids else None)
                return loaded_id, df
        return check_id, pd.DataFrame()

    # No filter specified: prefer combined CSV, then any single-check CSV
    combined_path = out_dir / "lightcurve_check_combined_raw.csv"
    if combined_path.exists():
        try:
            df = pd.read_csv(combined_path)
            cid: int | None = None
            if "check_id" in df.columns:
                ids = sorted({
                    int(x)
                    for x in pd.to_numeric(df["check_id"], errors="coerce").dropna().astype(int).tolist()
                })
                if len(ids) == 1:
                    cid = ids[0]
            return cid, df
        except Exception:
            pass

    check_id = load_check_star_id(result_dir)
    if check_id is not None:
        p = out_dir / f"lightcurve_check_ID{check_id}_raw.csv"
        if p.exists():
            try:
                return check_id, pd.read_csv(p)
            except Exception:
                pass

    for p in sorted(out_dir.glob("lightcurve_check_ID*_raw.csv")):
        try:
            cid = int(p.stem.replace("lightcurve_check_ID", "").replace("_raw", ""))
            return cid, pd.read_csv(p)
        except Exception:
            continue

    return None, pd.DataFrame()
