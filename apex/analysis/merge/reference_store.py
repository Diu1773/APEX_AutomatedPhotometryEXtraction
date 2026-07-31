"""Reference-mode storage for merged workspaces.

A merged workspace normally copies every frame's forced-photometry TSV into
itself ("full" mode), which duplicates the inputs on disk — a 3000-frame
multi-night merge costs as much again as the nights it merges.

Reference mode stores a pointer instead: which input workspace each merged
frame came from, plus the ID remap to apply when it is read. The loader
(:mod:`apex.utils.photometry_loader`) resolves that on the fly, so every
consumer keeps calling ``load_frame_photometry`` and sees the same columns.

The trade-off is deliberate and must be surfaced in the UI: a reference
workspace breaks if the input folders are moved or deleted, while a full one
stands alone.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Optional, Tuple

import pandas as pd

from apex.utils.io_utils import coerce_int64_source_id

STORAGE_FULL = "full"
STORAGE_REFERENCE = "reference"

REFERENCE_INDEX = "merge_reference.json"
ID_MAP_FILE = "merge_id_map.csv"

# {result_dir: (mtime, payload)} — the loader hits this once per frame, so the
# index and the remap tables are parsed once per workspace, not per frame.
_INDEX_CACHE: Dict[str, Tuple[float, dict]] = {}
_REMAP_CACHE: Dict[str, Tuple[float, dict]] = {}


def write_reference_index(out_dir: Path, frames: Dict[str, dict]) -> Path:
    """Record where each merged frame's photometry actually lives."""
    path = Path(out_dir) / REFERENCE_INDEX
    path.write_text(
        json.dumps({"storage_mode": STORAGE_REFERENCE, "frames": frames},
                   indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return path


def read_reference_index(result_dir: Path) -> dict:
    """``{merged filename: {dir, file, tag, filter}}``; empty when not reference."""
    path = Path(result_dir) / REFERENCE_INDEX
    if not path.exists():
        return {}
    key = str(path)
    try:
        mtime = path.stat().st_mtime
    except OSError:
        return {}
    cached = _INDEX_CACHE.get(key)
    if cached and cached[0] == mtime:
        return cached[1]
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    frames = data.get("frames", {}) if isinstance(data, dict) else {}
    frames = frames if isinstance(frames, dict) else {}
    _INDEX_CACHE[key] = (mtime, frames)
    return frames


def is_reference_workspace(result_dir: Path) -> bool:
    return bool(read_reference_index(result_dir))


def read_source_id_remap(result_dir: Path) -> Dict[Tuple[str, str], Dict[int, int]]:
    """``{(folder_tag, filter): {local source_id: merged source_id}}``.

    Built from ``merge_id_map.csv``, which the merge already writes.
    """
    path = Path(result_dir) / ID_MAP_FILE
    if not path.exists():
        return {}
    key = str(path)
    try:
        mtime = path.stat().st_mtime
    except OSError:
        return {}
    cached = _REMAP_CACHE.get(key)
    if cached and cached[0] == mtime:
        return cached[1]

    try:
        df = pd.read_csv(path)
    except Exception:
        return {}
    remap: Dict[Tuple[str, str], Dict[int, int]] = {}
    required = {"folder_tag", "filter", "local_source_id", "merged_source_id"}
    if required <= set(df.columns):
        local = coerce_int64_source_id(df["local_source_id"])
        merged = coerce_int64_source_id(df["merged_source_id"])
        for tag, flt, local_sid, merged_sid in zip(
            df["folder_tag"].astype(str), df["filter"].astype(str), local, merged
        ):
            if pd.isna(local_sid) or pd.isna(merged_sid):
                continue
            remap.setdefault((tag, str(flt)), {})[int(local_sid)] = int(merged_sid)
    _REMAP_CACHE[key] = (mtime, remap)
    return remap


def resolve_frame(result_dir: Path, merged_filename: str) -> Optional[dict]:
    """The input-workspace entry for one merged frame, or None."""
    entry = read_reference_index(result_dir).get(str(merged_filename))
    return entry if isinstance(entry, dict) else None


def source_photometry_path(entry: dict) -> Optional[Path]:
    """Path of the original per-frame TSV described by ``entry``."""
    source_dir = entry.get("dir")
    filename = entry.get("file")
    if not source_dir or not filename:
        return None
    from apex.utils.step_paths import forced_phot_input_dir

    forced_dir = forced_phot_input_dir(Path(source_dir))
    for name in (f"photometry_{filename}.tsv", f"{filename}_photometry.tsv"):
        candidate = forced_dir / name
        if candidate.exists():
            return candidate
    return None


def apply_reference_remap(df: pd.DataFrame, entry: dict,
                          remap: Dict[Tuple[str, str], Dict[int, int]],
                          merged_filename: str) -> pd.DataFrame:
    """Rewrite an input workspace's frame into merged identity.

    Mirrors exactly what full mode bakes into the copied TSV: merged
    ``source_id``, the merged filename in ``file``, and the provenance columns.
    ``ID`` is left to the caller's Step 8 ``sid_map`` pass, which maps merged
    source_id to the merged display ID.
    """
    tag = str(entry.get("tag", ""))
    flt = str(entry.get("filter", ""))
    table = remap.get((tag, flt), {})
    out = df.copy()
    if table and "source_id" in out.columns:
        local = coerce_int64_source_id(out["source_id"]).astype("Int64")
        out["source_id"] = local.map(table).astype("Int64")
        # A row whose source_id is not in the map did not survive the merge.
        out = out[out["source_id"].notna()].copy()
        if "ID" in out.columns:
            # The local ID is meaningless in the merged workspace; drop it so
            # the Step 8 catalog assigns the merged one.
            out = out.drop(columns=["ID"])
    out["file"] = merged_filename
    out["source_folder"] = tag
    out["original_file"] = str(entry.get("file", ""))
    return out
