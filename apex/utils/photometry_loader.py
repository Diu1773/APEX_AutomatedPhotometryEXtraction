from __future__ import annotations

import re
from collections import OrderedDict
from pathlib import Path

import pandas as pd

from .common_helpers import normalize_filter_key
from .io_utils import read_csv_int64_source_id, coerce_int64_source_id
from .step_paths import forced_phot_input_dir
from .step_paths_lc import selection_input_dirs


_DATE_RE = re.compile(r"(20\d{6})")


def _extract_date_key(filename: str) -> str:
    match = _DATE_RE.search(str(filename))
    return match.group(1) if match else ""


def _read_table(path: Path) -> pd.DataFrame | None:
    suffix = path.suffix.lower()
    preferred_seps = ["\t", ","] if suffix == ".tsv" else [",", "\t"]
    seen: set[str] = set()

    for sep in preferred_seps:
        if sep in seen:
            continue
        seen.add(sep)
        try:
            df = read_csv_int64_source_id(path, sep=sep)
        except Exception:
            continue
        if df is None:
            continue

        # Retry with the alternate separator when a CSV/TSV was parsed into a
        # single unsplit header column (common after the Step 5/8 refactor).
        if len(df.columns) == 1:
            only_col = str(df.columns[0])
            alt_sep = "," if sep == "\t" else "\t"
            if alt_sep in only_col:
                continue
        return df

    return None


def _resolve_photometry_path(result_dir: Path, fname: str) -> Path | None:
    forced_dir = forced_phot_input_dir(result_dir)
    for name in (f"photometry_{fname}.tsv", f"{fname}_photometry.tsv"):
        p = forced_dir / name
        if p.exists():
            return p
    return None


def _resolve_reference_frame(result_dir: Path, fname: str):
    """``(path, entry)`` for a frame stored by reference, else ``(None, None)``.

    A merged workspace built in reference mode keeps no per-frame TSV of its
    own; it points back at the input workspace it came from
    (:mod:`apex.analysis.merge.reference_store`).
    """
    from apex.analysis.merge import reference_store

    entry = reference_store.resolve_frame(Path(result_dir), fname)
    if not entry:
        return None, None
    return reference_store.source_photometry_path(entry), entry


def _resolve_idmatch_path(result_dir: Path, fname: str) -> Path | None:
    # New pipeline embeds source_id directly in the forced-phot TSV, so no
    # separate idmatch file is needed (callers use source_id straight from the
    # TSV). Legacy workspaces kept a per-frame det->source_id map under
    # cache/idmatch/; fall back to it so those frames can still be identified.
    root = Path(result_dir)
    date_key = _extract_date_key(fname)
    candidates = [root / "cache" / "idmatch" / f"idmatch_{fname}.csv"]
    if date_key:
        candidates.append(root / "step8_idmatch" / date_key / f"idmatch_{fname}.csv")
    candidates.append(root / "step8_idmatch" / f"idmatch_{fname}.csv")
    for path in candidates:
        if path.exists():
            return path
    matches = list((root / "step8_idmatch").glob(f"*/idmatch_{fname}.csv"))
    return matches[0] if matches else None


# Load-path counters for the LC benchmark (docs/audit/APEX_PERF_DEV_PLAN.md
# T0.4). Cheap module-level ints under the GIL; an experiment resets them, runs
# a build, and reads how many files/rows were actually parsed — the number the
# shared-cache optimisation must reduce.
LOAD_COUNTERS = {
    "frames_loaded": 0,     # Step 7 tables actually read from disk
    "rows_loaded": 0,       # rows parsed across those tables
    "sidmap_built": 0,      # source_id→ID map built by re-reading catalogs
    "sidmap_reused": 0,     # map supplied by the caller (the cheap path)
}


def reset_load_counters() -> None:
    for key in LOAD_COUNTERS:
        LOAD_COUNTERS[key] = 0


def get_load_counters() -> dict[str, int]:
    return dict(LOAD_COUNTERS)


def _load_source_to_id_map(result_dir: Path, filt_hint: str | None = None) -> dict[int, int]:
    selection_dirs = selection_input_dirs(result_dir)
    if not any(path.exists() for path in selection_dirs):
        return {}

    candidates: list[tuple[Path, str]] = []
    filt_key = normalize_filter_key(filt_hint)
    if filt_key:
        for step9_out in selection_dirs:
            candidates.extend(
                [
                    (step9_out / f"master_catalog_{filt_key}.tsv", "\t"),
                    (step9_out / f"id_mapping_{filt_key}.csv", ","),
                ]
            )
    for step9_out in selection_dirs:
        candidates.extend((p, "\t") for p in sorted(step9_out.glob("master_catalog_*.tsv")))
        candidates.extend((p, ",") for p in sorted(step9_out.glob("id_mapping_*.csv")))

    mapping: dict[int, int] = {}
    for path, sep in candidates:
        if not path.exists():
            continue
        try:
            df = read_csv_int64_source_id(path, sep=sep)
        except Exception:
            continue
        if not {"source_id", "ID"} <= set(df.columns):
            continue
        sid_vals = coerce_int64_source_id(df["source_id"])
        id_vals = pd.to_numeric(df["ID"], errors="coerce").astype("Int64")
        for sid_val, id_val in zip(sid_vals, id_vals):
            if pd.isna(sid_val) or pd.isna(id_val):
                continue
            sid_int = int(sid_val)
            if sid_int not in mapping:
                mapping[sid_int] = int(id_val)
    return mapping


def load_frame_photometry(
    result_dir: Path,
    fname: str,
    filt_hint: str | None = None,
    sid_map: dict[int, int] | None = None,
) -> pd.DataFrame | None:
    """Load Step 7 forced photometry and enrich it with Step 8 final IDs.

    Forced photometry already carries source identity. Downstream steps still
    need the final stable display `ID`, so this loader applies Step 8 selection
    catalogs when available.

    ``sid_map`` lets callers that load many frames from the same workspace pass
    a pre-built ``{source_id: ID}`` map so the Step 8 selection catalog is read
    once rather than re-globbed and re-parsed on every frame. Pass an empty dict
    to mean "no enrichment" while still skipping the per-frame catalog read.
    """

    reference_entry = None
    phot_path = _resolve_photometry_path(result_dir, fname)
    if phot_path is None:
        phot_path, reference_entry = _resolve_reference_frame(result_dir, fname)
    if phot_path is None:
        return None

    df = _read_table(phot_path)
    if df is not None:
        LOAD_COUNTERS["frames_loaded"] += 1
        LOAD_COUNTERS["rows_loaded"] += len(df)
    if df is None or df.empty:
        return df

    if reference_entry is not None:
        from apex.analysis.merge import reference_store

        df = reference_store.apply_reference_remap(
            df, reference_entry,
            reference_store.read_source_id_remap(Path(result_dir)), fname)
        if df.empty:
            return df

    df = df.copy()
    if "id" in df.columns and "ID" not in df.columns:
        df = df.rename(columns={"id": "ID"})
    if "det_idx" in df.columns and "det_uid" not in df.columns:
        df = df.rename(columns={"det_idx": "det_uid"})
    if "file" not in df.columns:
        df["file"] = fname

    if "mag" not in df.columns:
        for column in ("mag_inst", "mag_ap"):
            if column in df.columns:
                df["mag"] = pd.to_numeric(df[column], errors="coerce")
                break
    if "mag_err" not in df.columns and "mag_ap_err" in df.columns:
        df["mag_err"] = pd.to_numeric(df["mag_ap_err"], errors="coerce")

    if "det_uid" in df.columns:
        df["det_uid"] = pd.to_numeric(df["det_uid"], errors="coerce").astype("Int64")

    need_source_id = "source_id" not in df.columns or coerce_int64_source_id(df["source_id"]).notna().sum() == 0
    if need_source_id and "det_uid" in df.columns:
        idmatch_path = _resolve_idmatch_path(result_dir, fname)
        if idmatch_path is not None:
            idm = _read_table(idmatch_path)
            if idm is not None and not idm.empty:
                if "det_idx" in idm.columns and "det_uid" not in idm.columns:
                    idm = idm.rename(columns={"det_idx": "det_uid"})
                if "det_uid" in idm.columns and "source_id" in idm.columns:
                    idm = idm.copy()
                    idm["det_uid"] = pd.to_numeric(idm["det_uid"], errors="coerce").astype("Int64")
                    idm["source_id"] = coerce_int64_source_id(idm["source_id"]).astype("Int64")
                    merge_cols = ["det_uid", "source_id"]
                    for extra in (
                        "x", "y", "ra_deg", "dec_deg", "sep_arcsec",
                        "match_confidence",
                    ):
                        if extra in idm.columns and extra not in df.columns:
                            merge_cols.append(extra)
                    df = df.merge(idm[merge_cols], on="det_uid", how="left")

    if not {"x", "y"} <= set(df.columns):
        for x_col, y_col in (
            ("xcenter", "ycenter"),
            ("x_fit", "y_fit"),
            ("x_det", "y_det"),
        ):
            if {x_col, y_col} <= set(df.columns):
                df["x"] = pd.to_numeric(df[x_col], errors="coerce")
                df["y"] = pd.to_numeric(df[y_col], errors="coerce")
                break

    if "source_id" in df.columns:
        df["source_id"] = coerce_int64_source_id(df["source_id"]).astype("Int64")

    if "source_id" in df.columns:
        if sid_map is None:
            filt_key = filt_hint
            if not filt_key:
                for col in ("FILTER", "filter"):
                    if col in df.columns and not df.empty:
                        filt_key = normalize_filter_key(df[col].iloc[0])
                        break
            else:
                filt_key = normalize_filter_key(filt_key)
            sid_map = _load_source_to_id_map(result_dir, filt_key)
            LOAD_COUNTERS["sidmap_built"] += 1
        else:
            LOAD_COUNTERS["sidmap_reused"] += 1
        if sid_map:
            mapped_ids = df["source_id"].map(sid_map).astype("Int64")
            if "ID" in df.columns:
                # Legacy forced-phot tables may still carry a stale per-frame/local ID.
                existing_ids = pd.to_numeric(df["ID"], errors="coerce").astype("Int64")
                df["ID"] = mapped_ids.where(mapped_ids.notna(), existing_ids).astype("Int64")
            else:
                df["ID"] = mapped_ids

    return df


class FramePhotometryCache:
    """LRU cache of per-frame photometry tables, keyed by (result_dir, frame).

    The light-curve builder used to keep **one directory at a time** and throw
    the whole cache away whenever it changed::

        if self._photometry_cache_dir != result_dir:
            self._photometry_cache.clear()

    A multi-night workspace interleaves frames from several result directories,
    so every switch discarded everything and re-read it. What that costs is
    measured: three stars over 124 frames took 372 frame reads and 504,804 rows
    instead of 124 and 168,268 — the amplification is exactly linear in stars
    (benchmark/perf/20260809/b3_lc_load.json).

    Bounded by **bytes, not entries**: a frame table scales with the master
    catalog, so an entry count that is safe at 1,357 stars (0.7 MB each) is not
    safe at 5,000. The budget defaults to a quarter of the RAM that is free
    when the cache is first used, which is the same reasoning as the worker
    admission control in ``apex.utils.constants`` — the rest of memory belongs
    to the OS page cache holding the frames themselves.
    """

    _BUDGET_FRACTION = 0.25
    _BUDGET_MIN_MB = 256.0
    _BUDGET_MAX_MB = 4096.0

    def __init__(self, budget_mb: float | None = None):
        self._entries: OrderedDict[tuple[str, str], object] = OrderedDict()
        self._sizes: dict[tuple[str, str], int] = {}
        self._bytes = 0
        self._budget_mb = budget_mb
        self.hits = 0
        self.misses = 0
        self.evictions = 0

    @property
    def budget_bytes(self) -> int:
        if self._budget_mb is None:
            from apex.utils.constants import available_ram_mb

            available = available_ram_mb()
            share = (available * self._BUDGET_FRACTION) if available else self._BUDGET_MIN_MB
            self._budget_mb = min(max(share, self._BUDGET_MIN_MB), self._BUDGET_MAX_MB)
        return int(self._budget_mb * 1e6)

    @staticmethod
    def _key(result_dir, fname: str) -> tuple[str, str]:
        return (str(Path(result_dir)), str(fname))

    @staticmethod
    def _sizeof(df) -> int:
        try:
            return int(df.memory_usage(deep=True).sum())
        except Exception:
            return 0

    def get(self, result_dir, fname: str):
        """Cached table, or None. A hit moves the entry to the MRU end."""
        key = self._key(result_dir, fname)
        if key not in self._entries:
            self.misses += 1
            return None
        self._entries.move_to_end(key)
        self.hits += 1
        return self._entries[key]

    def put(self, result_dir, fname: str, df) -> None:
        key = self._key(result_dir, fname)
        if key in self._entries:
            self._bytes -= self._sizes.pop(key, 0)
            del self._entries[key]
        size = self._sizeof(df)
        self._entries[key] = df
        self._sizes[key] = size
        self._bytes += size
        budget = self.budget_bytes
        while self._bytes > budget and len(self._entries) > 1:
            old_key, _ = self._entries.popitem(last=False)
            self._bytes -= self._sizes.pop(old_key, 0)
            self.evictions += 1

    def __contains__(self, key) -> bool:
        result_dir, fname = key
        return self._key(result_dir, fname) in self._entries

    def __len__(self) -> int:
        return len(self._entries)

    def clear(self) -> None:
        self._entries.clear()
        self._sizes.clear()
        self._bytes = 0

    def stats(self) -> dict:
        return {"entries": len(self._entries),
                "bytes": self._bytes,
                "budget_bytes": self.budget_bytes,
                "hits": self.hits,
                "misses": self.misses,
                "evictions": self.evictions}
