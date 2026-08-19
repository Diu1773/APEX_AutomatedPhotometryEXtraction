"""The differential light curve itself — target minus comparison ensemble.

This is what LC actually computes, and until now it lived as a 271-line method
on a `QMainWindow` subclass. Not because it needed a window: it reads a
photometry index, looks up one star and its ensemble in each frame, and returns
a table. Every dependency it had on the window was data — the parameter model,
two caches, the night assignments — never a widget.

Moving it here is the same move Step 8 and Step 10 got, and it buys the same
thing: the window now *inherits* the calculation instead of owning a copy, so
`LightCurveBuilderWindow._build_ensemble_series is
RawLightCurveBuilder._build_ensemble_series` is true. Identity stops being
something a test has to re-establish after every edit.

One piece genuinely belonged to the window and stayed behind the seam:
preloading reports progress to a worker panel. That is two hook methods here,
no-ops by default, and the window overrides them with its Qt signals — so the
batch path skips the reporting and runs the same loads.
"""

from __future__ import annotations

import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time

from apex.analysis.light_curve.photometry_source_service import (
    load_lightcurve_frame_photometry,
)
from apex.utils.astro_utils import compute_airmass_from_header, compute_bjd_tdb_array
from apex.utils.common_helpers import (
    normalize_filter_key as _normalize_filter_key,
    parse_jd as _parse_jd,
    safe_float as _safe_float,
)
from apex.utils.constants import get_parallel_workers
from apex.utils.io_utils import (
    coerce_int64_source_id,
    read_csv_int64_source_id,
    load_headers_table as _load_headers_table_util,
    load_night_assignments as _load_night_assignments_util,
)
from apex.utils.night_utils import (
    fallback_night_key as _nu_fallback_night_key,
    fill_missing_night_ids as _nu_fill_missing_night_ids,
)
from apex.utils.photometry_provenance import build_photometry_provenance
from apex.utils.qc_utils import load_frame_excludes
from apex.utils.step_paths import forced_phot_input_dir
from apex.utils.step_paths_lc import (
    step2_cropped_dir,
    step6_refbuild_dir,
    step8_selection_dir,
    step9_lc_dir,
)
from apex.analysis.light_curve.lightcurve_output_service import annotate_raw_lightcurve, save_combined_raw_outputs, save_dataset_raw_outputs
from apex.analysis.light_curve.photometry_source_service import resolve_lightcurve_photometry_source
from apex.utils.constants import MAD_TO_SIGMA
from apex.utils.photometry_provenance import format_photometry_provenance


def _select_rows_by_source_id(df: pd.DataFrame, source_id: int | None) -> pd.DataFrame:
    if source_id is None or "source_id" not in df.columns:
        return pd.DataFrame()
    sid_series = coerce_int64_source_id(df["source_id"])
    return df.loc[sid_series == int(source_id)]

def _date_from_dateobs(date_obs: str | None) -> str:
    if not date_obs:
        return "unknown"
    try:
        t = Time(str(date_obs).strip())
        return t.to_datetime().strftime("%Y-%m-%d")
    except Exception:
        return "unknown"

def _load_headers_table(result_dir: Path) -> pd.DataFrame:
    return _load_headers_table_util(result_dir)

def _load_headers_map(result_dir: Path) -> dict:
    df = _load_headers_table(result_dir)
    if df.empty:
        return {}
    if "Filename" in df.columns and "DATE-OBS" in df.columns:
        return dict(zip(df["Filename"].astype(str), df["DATE-OBS"].astype(str)))
    return {}

def _load_night_assignments(result_dir: Path) -> dict[str, int]:
    """Load filename -> night_id mapping from step1 night_assignments.json."""
    return _load_night_assignments_util(result_dir)

def _extract_date_from_path(path: Path | str | None = None, fname: str = "") -> str:
    """폴더 경로 또는 파일명에서 날짜 추출 (YYYY-MM-DD 또는 YYYYMMDD)

    우선순위:
    1. 폴더 경로에서 날짜 추출 (result_dir 이름)
    2. 파일명에서 날짜 추출 (날짜__파일명 형식)
    """
    import re

    # 1. 폴더 경로에서 날짜 추출
    if path is not None:
        path_str = str(path)
        # 경로의 각 부분에서 날짜 패턴 찾기
        for part in reversed(Path(path_str).parts):
            # YYYY-MM-DD 패턴
            m = re.match(r"(\d{4}-\d{2}-\d{2})", part)
            if m:
                return m.group(1)
            # YYYYMMDD 패턴
            m = re.match(r"(\d{8})", part)
            if m:
                d = m.group(1)
                return f"{d[:4]}-{d[4:6]}-{d[6:8]}"

    # 2. 파일명에서 날짜 추출 (날짜__파일명 형식)
    if fname and "__" in fname:
        folder_part = fname.split("__")[0]
        m = re.match(r"(\d{4}-\d{2}-\d{2})", folder_part)
        if m:
            return m.group(1)
        m = re.match(r"(\d{8})", folder_part)
        if m:
            d = m.group(1)
            return f"{d[:4]}-{d[4:6]}-{d[6:8]}"

    return "unknown"

def _resolve_fits_path(
    data_dir: Path,
    result_dir: Path,
    fname: str,
    file_path_map: dict | None = None,
) -> Path | None:
    if isinstance(file_path_map, dict):
        mapped = file_path_map.get(fname)
        if mapped:
            return Path(mapped)
    candidates = [
        data_dir / fname,
        result_dir / fname,
        step2_cropped_dir(result_dir) / fname,
    ]
    for cand in candidates:
        if cand.exists():
            return cand
    return None

def _active_comparison_ids_for_filter(
    selection: dict,
    active_comp_ids: list[int],
) -> list[int]:
    """Limit the global Step 10 candidate list to one filter's selection."""
    active_set = {int(value) for value in active_comp_ids if value is not None}
    selected_ids = [
        int(value)
        for value in selection.get("comparison_ids", [])
        if value is not None
    ]
    if not selected_ids:
        return [int(value) for value in active_comp_ids if value is not None]
    return [value for value in selected_ids if value in active_set]

def _load_selection_ids_by_filter(result_dir: Path) -> dict:
    """필터별 selection 로드 (Step 9에서 저장한 selection_{filter}.json)"""
    step9_out = step8_selection_dir(result_dir)
    filter_selections = {}

    if not step9_out.exists():
        return {}

    id_map_cache: dict[str, dict[int, int]] = {}

    def _load_step8_id_map(flt: str) -> dict[int, int]:
        key = _normalize_filter_key(flt)
        if key in id_map_cache:
            return id_map_cache[key]
        mapping: dict[int, int] = {}
        candidates = [
            (step9_out / f"master_catalog_{key}.tsv", "\t"),
            (step9_out / f"id_mapping_{key}.csv", ","),
        ]
        for path, sep in candidates:
            if not path.exists():
                continue
            try:
                df = read_csv_int64_source_id(path, sep=sep)
            except Exception:
                continue
            id_map = _build_source_to_id_map(df)
            for sid_int, id_val in id_map.items():
                if sid_int not in mapping:
                    mapping[sid_int] = id_val
        id_map_cache[key] = mapping
        return mapping

    def _load_step6_refbuild_id_map() -> dict[int, int]:
        key = "__step6_refbuild__"
        if key in id_map_cache:
            return id_map_cache[key]
        candidates = [step6_refbuild_dir(result_dir) / "sourceid_to_ID.csv"]
        mapping: dict[int, int] = {}
        for path in candidates:
            if not path.exists():
                continue
            try:
                df = read_csv_int64_source_id(path)
            except Exception:
                continue
            id_map = _build_source_to_id_map(df)
            for sid_int, id_val in id_map.items():
                mapping[sid_int] = id_val
            if mapping:
                break
        id_map_cache[key] = mapping
        return mapping

    for sel_path in sorted(step9_out.glob("selection_*.json")):
        flt = sel_path.stem.replace("selection_", "")
        try:
            data = json.loads(sel_path.read_text(encoding="utf-8"))
            target_id = data.get("target_id")
            comp_ids = data.get("comparison_ids", [])
            target_source_id = data.get("target_source_id")
            comp_source_ids = data.get("comparison_source_ids", [])

            target_id_val = int(target_id) if target_id is not None else None
            comp_id_vals = [int(x) for x in comp_ids if x is not None]
            comp_source_id_vals = [int(x) for x in comp_source_ids if x is not None]

            # Recover IDs from source IDs when selection JSON has null/empty final IDs.
            if target_id_val is None and target_source_id is not None:
                sid_map = _load_step8_id_map(flt)
                if int(target_source_id) in sid_map:
                    target_id_val = int(sid_map[int(target_source_id)])
                else:
                    refbuild_map = _load_step6_refbuild_id_map()
                    if int(target_source_id) in refbuild_map:
                        target_id_val = int(refbuild_map[int(target_source_id)])
            if comp_source_id_vals:
                sid_map = _load_step8_id_map(flt)
                if not sid_map:
                    refbuild_map = _load_step6_refbuild_id_map()
                    sid_map = refbuild_map
                resolved_pairs: list[tuple[int, int]] = []
                seen_ids: set[int] = set()
                for source_id in comp_source_id_vals:
                    if source_id not in sid_map:
                        continue
                    resolved_id = int(sid_map[source_id])
                    if resolved_id in seen_ids:
                        continue
                    seen_ids.add(resolved_id)
                    resolved_pairs.append((resolved_id, source_id))
                if resolved_pairs:
                    comp_id_vals = [pair[0] for pair in resolved_pairs]
                    comp_source_id_vals = [pair[1] for pair in resolved_pairs]

            # Check star
            check_source_id = data.get("check_source_id")
            check_id = data.get("check_id")
            check_id_val = int(check_id) if check_id is not None else None
            if check_id_val is None and check_source_id is not None:
                sid_map = _load_step8_id_map(flt)
                if int(check_source_id) in sid_map:
                    check_id_val = int(sid_map[int(check_source_id)])
                else:
                    refbuild_map = _load_step6_refbuild_id_map()
                    if int(check_source_id) in refbuild_map:
                        check_id_val = int(refbuild_map[int(check_source_id)])

            filter_selections[flt] = {
                "target_id": target_id_val,
                "comparison_ids": comp_id_vals,
                "target_source_id": int(target_source_id) if target_source_id is not None else None,
                "comparison_source_ids": comp_source_id_vals,
                "check_id": check_id_val,
                "check_source_id": int(check_source_id) if check_source_id is not None else None,
            }
        except Exception:
            continue

    return filter_selections

def _load_target_radec(result_dir: Path, target_id: int) -> tuple[float, float]:
    """Look up target RA/Dec from master_catalog.tsv.

    Returns (ra_deg, dec_deg), or (nan, nan) if not found.
    """
    step9_out = step8_selection_dir(result_dir)
    candidates = list(step9_out.glob("master_catalog_*.tsv")) if step9_out.exists() else []
    candidates += [step9_out / "master_catalog.tsv"] if step9_out.exists() else []
    for path in candidates:
        if not path.exists():
            continue
        try:
            df = read_csv_int64_source_id(path, sep="\t")
            if "ID" not in df.columns:
                continue
            row = df[pd.to_numeric(df["ID"], errors="coerce") == int(target_id)]
            if row.empty:
                continue
            ra = float(pd.to_numeric(row["ra_deg"].values[0], errors="coerce")) if "ra_deg" in df.columns else np.nan
            dec = float(pd.to_numeric(row["dec_deg"].values[0], errors="coerce")) if "dec_deg" in df.columns else np.nan
            if np.isfinite(ra) and np.isfinite(dec):
                return ra, dec
        except Exception:
            continue
    return np.nan, np.nan


def _build_source_to_id_map(df: pd.DataFrame) -> dict[int, int]:
    if not {"source_id", "ID"} <= set(df.columns):
        return {}
    sid_vals = coerce_int64_source_id(df["source_id"])
    id_vals = pd.to_numeric(df["ID"], errors="coerce").astype("Int64")
    mapping: dict[int, int] = {}
    for sid_val, id_val in zip(sid_vals, id_vals):
        if pd.isna(sid_val) or pd.isna(id_val):
            continue
        sid_int = int(sid_val)
        if sid_int not in mapping:
            mapping[sid_int] = int(id_val)
    return mapping

def _load_check_star_meta_by_filter(result_dir: Path) -> dict[str, dict[str, int]]:
    """Load per-filter check star metadata from Step 9 selection JSONs."""
    out: dict[str, dict[str, int]] = {}
    filter_sel = _load_selection_ids_by_filter(result_dir)
    for flt, sel in sorted(filter_sel.items()):
        key = _normalize_filter_key(flt)
        if not key:
            continue
        entry: dict[str, int] = {}
        check_id = sel.get("check_id")
        check_source_id = sel.get("check_source_id")
        if check_id is not None:
            entry["check_id"] = int(check_id)
        if check_source_id is not None:
            entry["check_source_id"] = int(check_source_id)
        if entry:
            out[key] = entry
    return out

class RawLightCurveBuilder:
    """Builds the raw differential light curve of one target.

    A mixin rather than a standalone service, because the window already holds
    the state this needs — `params`, the photometry and header caches, the file
    manager's night assignments — and re-plumbing all of it through arguments
    would have meant rewriting the body rather than moving it. Rewriting is how
    the two paths drift apart.

    A user of this class supplies:
      `params`              the parameter model (reads `params.P`)
      `_photometry_cache`   a `FramePhotometryCache`
      `_header_cache`       a dict for FITS headers
      `file_manager`        optional; contributes night assignments
      `log`                 optional; a line-oriented logger
    """


    # Instance state a user supplies. Named here so a batch runner can be built
    # without reading the window to find out what it happens to set.
    datasets: list = []                      # [(label, result_dir), ...]
    comp_candidate_ids: list = []
    comp_ids_list: list = []
    runtime_mode: bool = False
    qc_sigma: float = 3.0
    project_state = None
    file_manager = None
    _force_aperture_for_datasets = False

    def log(self, message: str) -> None:
        """Where progress goes when nobody is watching a window."""
        print(message)

    # -- progress reporting: the one thing that stayed with the window --------

    def _preload_progress(self, index: int, message: str) -> None:
        """Called once per worker as preloading starts. No-op off the GUI."""

    def _preload_finished(self) -> None:
        """Called when preloading ends. No-op off the GUI."""


    def _photometry_source_for_dir(self, result_dir: Path) -> dict:
        source = self._resolved_photometry_source_for_dir(result_dir)
        if not self._force_aperture_for_datasets or source.get("source") != "psf":
            return source
        aperture_dir = forced_phot_input_dir(result_dir)
        return {
            **source,
            **build_photometry_provenance("aperture", "mag", "mag_err"),
            "directory": aperture_dir,
            "index_path": aperture_dir / "photometry_index.csv",
            "reason": "PSF unavailable for every dataset; using aperture for all",
        }

    def _load_active_photometry_index(self, result_dir: Path) -> pd.DataFrame:
        source = self._photometry_source_for_dir(result_dir)
        idx_path = Path(source["index_path"])
        try:
            index = pd.read_csv(idx_path)
        except Exception:
            return pd.DataFrame()
        if source.get("source") == "psf" and "file" in index.columns:
            allowed = {Path(str(frame)).name for frame in source.get("frames", [])}
            basenames = index["file"].astype(str).map(lambda value: Path(value).name)
            index = index[basenames.isin(allowed)].reset_index(drop=True)
        return index

    def _get_photometry_df(self, result_dir: Path, fname: str) -> pd.DataFrame | None:
        """Load photometry TSV with caching."""
        cached = self._photometry_cache.get(result_dir, fname)
        if cached is not None:
            return cached

        source = self._photometry_source_for_dir(result_dir)
        df = load_lightcurve_frame_photometry(result_dir, fname, source)

        self._photometry_cache.put(result_dir, fname, df)
        return df

    def _preload_photometry_cache(self, result_dir: Path, filenames: list[str]):
        """Bulk-preload photometry TSVs using ThreadPoolExecutor."""
        to_load = [fn for fn in filenames
                   if (result_dir, fn) not in self._photometry_cache]
        if not to_load:
            return

        n_workers = min(get_parallel_workers(self.params), len(to_load))
        for i in range(n_workers):
            self._preload_progress(i, f"Preloading {len(to_load)} files...")

        source = self._photometry_source_for_dir(result_dir)

        def _load_one(fname):
            df = load_lightcurve_frame_photometry(result_dir, fname, source)
            return fname, df

        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            for fname, df in pool.map(_load_one, to_load):
                self._photometry_cache.put(result_dir, fname, df)

        self._preload_finished()

    def _get_header(self, result_dir: Path, fname: str, cache: dict) -> fits.Header | None:
        if fname in cache:
            return cache[fname]
        fpath = _resolve_fits_path(
            Path(self.params.P.data_dir),
            result_dir,
            fname,
            getattr(self.params.P, "file_path_map", None),
        )
        if fpath is None:
            cache[fname] = None
            return None
        try:
            hdr = fits.getheader(fpath)
        except Exception:
            hdr = None
        cache[fname] = hdr
        return hdr

    def _map_comp_source_id(self, sel: dict, comp_id: int) -> int | None:
        comp_ids_all = sel.get("comparison_ids", [])
        comp_source_ids_all = sel.get("comparison_source_ids", [])
        try:
            idx = comp_ids_all.index(comp_id)
        except ValueError:
            return None
        if 0 <= idx < len(comp_source_ids_all):
            return comp_source_ids_all[idx]
        return None

    def _get_frame_exclude_map(self, result_dir: Path) -> dict[str, set[str]]:
        key = str(result_dir)
        if key not in self._frame_exclude_cache:
            # LC excludes live in lc_lightcurve/ (mode-specific); falls back to
            # result_dir/ for backward compat with projects saved before this change.
            self._frame_exclude_cache[key] = load_frame_excludes(
                result_dir, exclude_dir=step9_lc_dir(result_dir)
            )
        return self._frame_exclude_cache[key]

    def _compute_airmass(self, header: fits.Header | None) -> float:
        if header is None:
            return np.nan
        lat = float(getattr(self.params.P, "site_lat_deg", 0.0))
        lon = float(getattr(self.params.P, "site_lon_deg", 0.0))
        alt = float(getattr(self.params.P, "site_alt_m", 0.0))
        tz = float(getattr(self.params.P, "site_tz_offset_hours", 0.0))
        formula = getattr(self.params.P, "airmass_formula", None)
        try:
            hdr = header.copy()
            if "AIRMASS" in hdr:
                del hdr["AIRMASS"]
        except Exception:
            hdr = header
        info = compute_airmass_from_header(hdr, lat, lon, alt, tz, formula=formula)
        am = _safe_float(info.get("airmass", np.nan))
        if np.isfinite(am):
            return float(am)
        return _safe_float(header.get("AIRMASS", np.nan))

    def _compute_bjd_array(self, tarr: np.ndarray, result_dir: Path, target_id: int) -> np.ndarray:
        """Compute BJD_TDB array; returns NaN array if site coords or target RA/Dec are missing."""
        site_lat = float(getattr(self.params.P, "site_lat_deg", np.nan))
        site_lon = float(getattr(self.params.P, "site_lon_deg", np.nan))
        site_alt = float(getattr(self.params.P, "site_alt_m", 0.0))
        tgt_ra, tgt_dec = _load_target_radec(result_dir, target_id)
        if not np.isfinite(tgt_ra):
            # 설정의 [target] 좌표로 폴백한다. 예전 코드는 `P.target.ra_deg` 를
            # 봤는데 파라미터 모델은 그런 중첩 속성을 만들지 않는다
            # (`P.target_ra_deg`) — hasattr 이 늘 False 라 **폴백이 죽어 있었다**.
            # 그 탓에 Step 8 을 거치지 않은 워크스페이스에서는 설정에 좌표가
            # 멀쩡히 있는데도 BJD_TDB 가 전부 NaN 으로 나왔다.
            def _cfg(*names):
                for n in names:
                    try:
                        v = float(getattr(self.params.P, n, np.nan))
                    except (TypeError, ValueError):
                        continue
                    if np.isfinite(v):
                        return v
                return np.nan

            tgt_ra = _cfg("target_ra_deg", "ra_deg")
            tgt_dec = _cfg("target_dec_deg", "dec_deg")
        if np.isfinite(site_lat) and np.isfinite(site_lon) and np.isfinite(tgt_ra) and np.isfinite(tgt_dec):
            return compute_bjd_tdb_array(tarr, tgt_ra, tgt_dec, site_lat, site_lon, site_alt)
        return np.full(len(tarr), np.nan)

    def _build_ensemble_series(
        self,
        result_dir: Path,
        target_id: int,
        comp_ids: list[int],
        verbose: bool = True,
        target_source_id_by_filter: dict[str, int] | None = None,
    ) -> pd.DataFrame:
        source_info = self._photometry_source_for_dir(result_dir)
        photometry_source = source_info["source"]
        idx_path = Path(source_info["index_path"])
        if not idx_path.exists():
            if verbose:
                self.log(f"[WARN] photometry_index.csv not found in {result_dir}")
            return pd.DataFrame()
        idx = self._load_active_photometry_index(result_dir)
        if "file" not in idx.columns:
            if verbose:
                self.log(f"[WARN] photometry_index.csv missing 'file' column")
            return pd.DataFrame()
        # Apply frame exclusions (manual QC from step9 D-key) before building series
        exclude_map = self._get_frame_exclude_map(result_dir)
        if exclude_map:
            before = len(idx)
            idx = idx[~idx["file"].astype(str).isin(exclude_map.keys())].reset_index(drop=True)
            if verbose and before != len(idx):
                self.log(f"[Frame QC] Excluded {before - len(idx)} frame(s) from ensemble series")
        files = idx["file"].astype(str).tolist()
        self._preload_photometry_cache(result_dir, files)
        headers_map = _load_headers_map(result_dir)
        headers_df = _load_headers_table(result_dir)

        filter_selections = _load_selection_ids_by_filter(result_dir)

        filter_map = {}
        if "filter" in idx.columns:
            filter_map = dict(zip(idx["file"].astype(str), idx["filter"].astype(str)))
        elif "FILTER" in idx.columns:
            filter_map = dict(zip(idx["file"].astype(str), idx["FILTER"].astype(str)))

        header_filter_map = {}
        if not headers_df.empty and "Filename" in headers_df.columns:
            for col in ("FILTER", "filter"):
                if col in headers_df.columns:
                    header_filter_map = dict(zip(headers_df["Filename"].astype(str), headers_df[col].astype(str)))
                    break

        # Night assignment map: basename→night_id, merged from all sources
        night_id_map: dict[str, int] = {}
        raw_na = _load_night_assignments(result_dir)
        night_id_map.update({Path(str(k)).name: int(v) for k, v in raw_na.items() if int(v) > 0})
        fm = getattr(self, "file_manager", None)
        if fm is not None:
            for k, v in getattr(fm, "night_assignments", {}).items():
                if int(v) > 0:
                    night_id_map[Path(str(k)).name] = int(v)
        if "night_id" in idx.columns:
            for fn, nid in zip(idx["file"].astype(str), pd.to_numeric(idx["night_id"], errors="coerce")):
                bn = Path(fn).name
                if not pd.isna(nid) and int(nid) > 0 and bn not in night_id_map:
                    night_id_map[bn] = int(nid)

        # Use instance-level header cache (self._header_cache)
        times = []
        dates = []
        night_ids = []
        night_fallback_keys = []          # DATE-OBS-derived night, headless case
        _tz = float(getattr(self.params.P, "site_tz_offset_hours", 0.0) or 0.0)
        filters = []
        airmasses = []
        mags = []
        mag_errs = []
        comp_avgs = []
        comp_errs = []
        diffs = []
        diff_errs = []

        # headers.csv에서 airmass 컬럼 확인
        header_airmass_map = {}
        if not headers_df.empty and "Filename" in headers_df.columns:
            for col in ("AIRMASS", "airmass", "AM"):
                if col in headers_df.columns:
                    for _, row in headers_df.iterrows():
                        fn = str(row["Filename"])
                        am_val = pd.to_numeric(row[col], errors="coerce")
                        if np.isfinite(am_val):
                            header_airmass_map[fn] = float(am_val)
                    break

        n_target_found = 0
        n_comp_found = 0

        for fname in files:
            # 1) DATE-OBS, FILTER, AIRMASS: headers.csv에서 먼저 시도
            date_obs = headers_map.get(fname) if headers_map else None
            jd = _parse_jd(date_obs)

            filt_val = filter_map.get(fname, "")
            if not filt_val:
                filt_val = header_filter_map.get(fname, "")

            am = header_airmass_map.get(fname, np.nan)

            # 2) 정보가 부족한 경우에만 FITS 헤더 읽기 (lazy load)
            need_fits = (not date_obs) or (not filt_val) or (not np.isfinite(am))
            if need_fits:
                hdr = self._get_header(result_dir, fname, self._header_cache)
                if hdr is not None:
                    if not date_obs:
                        date_obs = hdr.get("DATE-OBS")
                        jd = _parse_jd(date_obs)
                    if not filt_val:
                        filt_val = hdr.get("FILTER", hdr.get("FILTER1", hdr.get("FILTER2", "")))
                    if not np.isfinite(am):
                        am = self._compute_airmass(hdr)

            times.append(jd)
            dates.append(_date_from_dateobs(date_obs) if date_obs else _extract_date_from_path(result_dir, fname))
            night_ids.append(night_id_map.get(fname, 0))
            night_fallback_keys.append(
                "" if night_id_map.get(fname, 0) > 0
                else _nu_fallback_night_key(date_obs, _tz))
            filt_key = _normalize_filter_key(filt_val)
            filters.append(filt_key)
            airmasses.append(am if np.isfinite(am) else np.nan)

            df = self._get_photometry_df(result_dir, fname)
            if df is None or df.empty:
                mags.append(np.nan)
                mag_errs.append(np.nan)
                comp_avgs.append(np.nan)
                comp_errs.append(np.nan)
                diffs.append(np.nan)
                diff_errs.append(np.nan)
                continue

            use_source_id = False
            target_source_id = None
            comp_source_map: dict[int, int] = {}
            frame_comp_ids = [int(cid) for cid in comp_ids]
            if filt_key in filter_selections:
                sel = filter_selections[filt_key]
                frame_comp_ids = _active_comparison_ids_for_filter(sel, comp_ids)
                override_target_sid = None
                if target_source_id_by_filter is not None:
                    override_target_sid = target_source_id_by_filter.get(filt_key)
                if override_target_sid is not None:
                    target_source_id = int(override_target_sid)
                else:
                    sel_target_id = sel.get("target_id")
                    if sel_target_id is not None and int(sel_target_id) == int(target_id):
                        target_source_id = sel.get("target_source_id")
                    else:
                        target_source_id = self._map_comp_source_id(sel, target_id)
                for cid in frame_comp_ids:
                    sid = self._map_comp_source_id(sel, cid)
                    if sid is not None:
                        comp_source_map[int(cid)] = int(sid)
                use_source_id = True

            # Target
            row_t = pd.DataFrame()
            if use_source_id and target_source_id is not None and "source_id" in df.columns:
                row_t = _select_rows_by_source_id(df, int(target_source_id))
            if row_t.empty and "ID" in df.columns:
                row_t = df[df["ID"] == int(target_id)]

            if not row_t.empty:
                n_target_found += 1
                tmag = _safe_float(row_t["mag"].values[0])
                terr = _safe_float(row_t["mag_err"].values[0]) if "mag_err" in row_t.columns else np.nan
            else:
                tmag = np.nan
                terr = np.nan

            # Comparison ensemble
            cmags = []
            cerrs = []
            for cid in frame_comp_ids:
                row_c = pd.DataFrame()
                if use_source_id and cid in comp_source_map and "source_id" in df.columns:
                    row_c = _select_rows_by_source_id(df, int(comp_source_map[cid]))
                if row_c.empty and "ID" in df.columns:
                    row_c = df[df["ID"] == int(cid)]
                if not row_c.empty and np.isfinite(_safe_float(row_c["mag"].values[0])):
                    cmags.append(_safe_float(row_c["mag"].values[0]))
                    cerrs.append(_safe_float(row_c["mag_err"].values[0]) if "mag_err" in row_c.columns else np.nan)
            if cmags:
                n_comp_found += 1
                cmags_arr = np.array(cmags, dtype=float)
                cerrs_arr = np.array(cerrs, dtype=float)
                valid_w = np.isfinite(cerrs_arr) & (cerrs_arr > 0)
                if np.any(valid_w):
                    w = 1.0 / (cerrs_arr[valid_w] ** 2)
                    comp_mean = float(np.sum(cmags_arr[valid_w] * w) / np.sum(w))
                    comp_err = float(1.0 / np.sqrt(np.sum(w)))
                else:
                    comp_mean = float(np.nanmean(cmags_arr))
                    comp_err = float(np.nanmean(cerrs_arr)) if cerrs else np.nan
            else:
                comp_mean = np.nan
                comp_err = np.nan

            mags.append(tmag)
            mag_errs.append(terr)
            comp_avgs.append(comp_mean)
            comp_errs.append(comp_err)

            if np.isfinite(tmag) and np.isfinite(comp_mean):
                diff = tmag - comp_mean
                diffs.append(diff)
                if np.isfinite(terr) and np.isfinite(comp_err):
                    diff_errs.append(float(np.sqrt(terr * terr + comp_err * comp_err)))
                else:
                    diff_errs.append(terr if np.isfinite(terr) else np.nan)
            else:
                diffs.append(np.nan)
                diff_errs.append(np.nan)

        # Headless workspaces carry no night assignments (the classifier is a
        # GUI Step 1 mixin) — every frame then landed on night_id 0 and the
        # nightly-offset correction treated all nights as one. Fill the gap
        # from DATE-OBS with the app's one night definition (P1 noon split).
        filled_ids, inferred = _nu_fill_missing_night_ids(
            night_ids, night_fallback_keys,
            start_after=max(night_id_map.values(), default=0))
        if inferred:
            night_ids = filled_ids
            if verbose:
                self.log(f"[NIGHT] no stored night assignments — inferred "
                         f"{len(inferred)} night(s) from DATE-OBS: "
                         + ", ".join(f"N{v}={k}" for k, v in sorted(
                             inferred.items(), key=lambda kv: kv[1])))

        tarr = np.array(times, float)
        if np.all(~np.isfinite(tarr)):
            tarr = np.arange(len(files), dtype=float)
        t0 = np.nanmedian(tarr)
        rel_time_hr = (tarr - t0) * 24.0

        # BJD_TDB 계산
        bjd_arr = self._compute_bjd_array(tarr, result_dir, target_id)

        if verbose:
            total = len(files)
            self.log(f"Ensemble series (Target={target_id}): frames={total}, target={n_target_found}, comp={n_comp_found}")
            if np.any(np.isfinite(bjd_arr)):
                delta = np.nanmedian(bjd_arr - tarr) * 86400
                self.log(f"[BJD] BJD_TDB computed, median correction {delta:+.1f}s")
            else:
                self.log("[BJD] BJD_TDB not computed (missing site coords or target RA/Dec)")

        return pd.DataFrame({
            "file": files,
            "photometry_source": [photometry_source] * len(files),
            "mag_input_column": [source_info["mag_column"]] * len(files),
            "mag_error_input_column": [source_info["mag_error_column"]] * len(files),
            "filter": filters,
            "date": dates,
            "night_id": night_ids,
            "JD": tarr,
            "BJD_TDB": bjd_arr,
            "rel_time_hr": rel_time_hr,
            "mag": np.array(mags, float),
            "mag_err": np.array(mag_errs, float),
            "comp_avg": np.array(comp_avgs, float),
            "comp_err": np.array(comp_errs, float),
            "diff_mag_raw": np.array(diffs, float),
            "diff_err": np.array(diff_errs, float),
            "airmass": np.array(airmasses, float),
        })

    def _refresh_photometry_source_policy(self, *, log: bool = False) -> None:
        sources = [
            self._resolved_photometry_source_for_dir(Path(result_dir)).get("source")
            for _, result_dir in self.datasets
        ]
        force_aperture = len(sources) > 1 and any(source != "psf" for source in sources)
        changed = force_aperture != self._force_aperture_for_datasets
        self._force_aperture_for_datasets = force_aperture
        if changed:
            self._photometry_cache.clear()
            self._diff_series_cache.clear()
            self._check_series_cache.clear()
        source_label = self.__dict__.get("photometry_source_label")
        if sources and source_label is not None:
            active = "aperture" if force_aperture else str(sources[0])
            reason = (
                "PSF unavailable for every dataset; aperture enforced for all"
                if force_aperture
                else self._photometry_source_for_dir(Path(self.datasets[0][1])).get(
                    "reason", ""
                )
            )
            active_info = self._photometry_source_for_dir(Path(self.datasets[0][1]))
            source_label.setText(format_photometry_provenance(active_info))
            source_label.setToolTip(reason)
        if log and sources:
            active = "aperture" if force_aperture else str(sources[0])
            detail = ", ".join(str(source) for source in sources)
            self.log(f"[Photometry] Active source: {active} (datasets: {detail})")

    def _resolved_photometry_source_for_dir(self, result_dir: Path) -> dict:
        key = str(Path(result_dir).resolve())
        if key not in self._photometry_source_cache:
            current = Path(self.params.P.result_dir).resolve()
            state = self.project_state if Path(result_dir).resolve() == current else None
            self._photometry_source_cache[key] = resolve_lightcurve_photometry_source(
                result_dir, state
            )
        return self._photometry_source_cache[key]

    def _build_star_mag_series(
        self,
        result_dir: Path,
        star_id: int,
        verbose: bool = True,
        include_excluded: bool = False,
        files_override: list[str] | None = None,
        preload: bool = True,
    ) -> pd.DataFrame:
        source_info = self._photometry_source_for_dir(result_dir)
        photometry_source = source_info["source"]
        idx_path = Path(source_info["index_path"])
        if not idx_path.exists():
            if verbose:
                self.log(f"[WARN] photometry_index.csv not found in {result_dir}")
            return pd.DataFrame()
        idx = self._load_active_photometry_index(result_dir)
        if "file" not in idx.columns:
            if verbose:
                self.log(f"[WARN] photometry_index.csv missing 'file' column")
            return pd.DataFrame()
        exclude_map = self._get_frame_exclude_map(result_dir)
        if exclude_map and not include_excluded:
            before = len(idx)
            idx = idx[~idx["file"].astype(str).isin(exclude_map.keys())]
            if verbose and before != len(idx):
                self.log(f"[Frame QC] Excluded frames: {before} → {len(idx)}")
        available_files = idx["file"].astype(str).tolist()
        if files_override is not None:
            allowed = set(available_files)
            files = [str(fname) for fname in files_override if str(fname) in allowed]
        else:
            files = available_files
        if preload:
            self._preload_photometry_cache(result_dir, files)
        headers_map = _load_headers_map(result_dir)
        headers_df = _load_headers_table(result_dir)

        filter_selections = _load_selection_ids_by_filter(result_dir)

        filter_map = {}
        if "filter" in idx.columns:
            filter_map = dict(zip(idx["file"].astype(str), idx["filter"].astype(str)))
        elif "FILTER" in idx.columns:
            filter_map = dict(zip(idx["file"].astype(str), idx["FILTER"].astype(str)))

        header_filter_map = {}
        if not headers_df.empty and "Filename" in headers_df.columns:
            for col in ("FILTER", "filter"):
                if col in headers_df.columns:
                    header_filter_map = dict(zip(headers_df["Filename"].astype(str), headers_df[col].astype(str)))
                    break

        # Use instance-level header cache (self._header_cache)
        times = []
        dates = []
        filters = []
        mags = []
        mag_errs = []

        n_found = 0
        n_fits_read = 0

        for fname in files:
            # 1) DATE-OBS, FILTER: headers.csv에서 먼저 시도
            date_obs = headers_map.get(fname) if headers_map else None
            jd = _parse_jd(date_obs)

            filt_val = filter_map.get(fname, "")
            if not filt_val:
                filt_val = header_filter_map.get(fname, "")

            # 2) 정보가 부족한 경우에만 FITS 헤더 읽기 (lazy load)
            need_fits = (not np.isfinite(jd)) or (not filt_val) or (not date_obs)
            hdr = None
            if need_fits:
                hdr = self._get_header(result_dir, fname, self._header_cache)
                if hdr is not None:
                    n_fits_read += 1
                    if not np.isfinite(jd):
                        jd = _parse_jd(hdr.get("DATE-OBS"))
                    if not date_obs:
                        date_obs = hdr.get("DATE-OBS")
                    if not filt_val:
                        filt_val = hdr.get("FILTER", hdr.get("FILTER1", hdr.get("FILTER2", "")))

            times.append(jd)
            if date_obs:
                dates.append(_date_from_dateobs(date_obs))
            else:
                dates.append(_extract_date_from_path(result_dir, fname))
            filt_key = _normalize_filter_key(filt_val)
            filters.append(filt_key)

            # Use photometry cache
            df = self._get_photometry_df(result_dir, fname)
            if df is None or df.empty:
                mags.append(np.nan)
                mag_errs.append(np.nan)
                continue

            use_source_id = False
            comp_source_id = None
            if filt_key in filter_selections:
                sel = filter_selections[filt_key]
                comp_source_id = self._map_comp_source_id(sel, star_id)
                use_source_id = True

            row = pd.DataFrame()
            if use_source_id and comp_source_id is not None and "source_id" in df.columns:
                row = _select_rows_by_source_id(df, int(comp_source_id))
            if row.empty and "ID" in df.columns:
                row = df[df["ID"] == int(star_id)]

            if row.empty:
                mags.append(np.nan)
                mag_errs.append(np.nan)
            else:
                n_found += 1
                mags.append(_safe_float(row["mag"].values[0]))
                mag_errs.append(_safe_float(row["mag_err"].values[0]) if "mag_err" in row.columns else np.nan)

        tarr = np.array(times, float)
        if np.all(~np.isfinite(tarr)):
            tarr = np.arange(len(files), dtype=float)
        t0 = np.nanmedian(tarr)
        rel_time_hr = (tarr - t0) * 24.0

        if verbose:
            total = len(files)
            self.log(f"Star series ID={star_id}: {n_found}/{total} frames")

        return pd.DataFrame({
            "file": files,
            "photometry_source": [photometry_source] * len(files),
            "mag_input_column": [source_info["mag_column"]] * len(files),
            "mag_error_input_column": [source_info["mag_error_column"]] * len(files),
            "filter": filters,
            "date": dates,
            "JD": tarr,
            "rel_time_hr": rel_time_hr,
            "mag": np.array(mags, float),
            "mag_err": np.array(mag_errs, float),
        })

    def _build_check_star_series(
        self,
        result_dir: Path,
        comp_ids: list[int],
        verbose: bool = False,
    ) -> tuple[dict[str, int], pd.DataFrame]:
        # 캐시 확인
        _ck_key = (str(result_dir), tuple(sorted(comp_ids)))
        if _ck_key in self._check_series_cache:
            return self._check_series_cache[_ck_key]

        meta_by_filter = _load_check_star_meta_by_filter(result_dir)
        check_ids_by_filter = {
            flt: int(meta["check_id"])
            for flt, meta in meta_by_filter.items()
            if meta.get("check_id") is not None
        }
        check_source_ids_by_filter = {
            flt: int(meta["check_source_id"])
            for flt, meta in meta_by_filter.items()
            if meta.get("check_source_id") is not None
        }
        if not check_ids_by_filter and not check_source_ids_by_filter:
            return {}, pd.DataFrame()

        fallback_id = next(iter(check_ids_by_filter.values()), -1)
        check_df = self._build_ensemble_series(
            result_dir,
            int(fallback_id),
            comp_ids,
            verbose=verbose,
            target_source_id_by_filter=check_source_ids_by_filter or None,
        )
        if check_df.empty:
            return check_ids_by_filter, check_df

        valid_filters = set(check_ids_by_filter) | set(check_source_ids_by_filter)
        if "filter" in check_df.columns and valid_filters:
            filter_keys = check_df["filter"].astype(str).map(_normalize_filter_key)
            check_df = check_df[filter_keys.isin(valid_filters)].copy()
            if not check_df.empty:
                row_filter_keys = filter_keys.loc[check_df.index]
                check_df["check_id"] = row_filter_keys.map(check_ids_by_filter).astype("Int64")
                check_df["check_source_id"] = row_filter_keys.map(
                    check_source_ids_by_filter
                ).astype("Int64")
        identity_mask = pd.Series(False, index=check_df.index)
        for column in ("check_id", "check_source_id"):
            if column in check_df.columns:
                identity_mask |= pd.to_numeric(check_df[column], errors="coerce").notna()
        if not check_df.empty:
            check_df = check_df[identity_mask].copy()
        if "JD" in check_df.columns and not check_df.empty:
            check_df = check_df.sort_values("JD").reset_index(drop=True)
        result = (check_ids_by_filter, check_df)
        self._check_series_cache[_ck_key] = result
        return result

    def _compute_comp_qc(
        self,
        result_dir: Path,
        target_id: int,
        comp_ids: list[int],
        files_override: list[str] | None = None,
        verbose: bool = True,
    ) -> list[dict]:
        """Check-star QC: each comp is treated as target against remaining comps as ensemble."""
        rows = []
        active_set = set(self.comp_ids_list)
        if not hasattr(self, "_qc_checkstar_cache"):
            self._qc_checkstar_cache: dict[int, pd.DataFrame] = {}
        self._frame_exclude_cache: dict[str, dict[str, set[str]]] = {}

        for comp_id in comp_ids:
            other_comps = [c for c in comp_ids if c != comp_id]
            if not other_comps:
                # Only 1 comp — fall back to absolute mag scatter
                df = self._build_star_mag_series(result_dir, comp_id, verbose=False,
                                                  files_override=files_override)
                if df.empty:
                    rows.append({"comp_id": int(comp_id), "n": 0, "rms": np.nan, "mad": np.nan,
                                 "sigma_nights": np.nan, "outliers": 0, "outlier_frac": np.nan,
                                 "use": comp_id in active_set})
                    continue
                y = df["mag"].to_numpy(float)
                m = np.isfinite(y)
                yv = y[m]
                med = np.nanmedian(yv)
                yv = yv - med
                rms = float(np.nanstd(yv))
                mad = float(np.nanmedian(np.abs(yv - np.nanmedian(yv))))
                rows.append({"comp_id": int(comp_id), "n": int(m.sum()), "rms": rms, "mad": mad,
                             "sigma_nights": np.nan, "outliers": 0, "outlier_frac": np.nan,
                             "use": comp_id in active_set})
                continue

            # Check star: build diff series with this comp as "target"
            df = self._build_ensemble_series(result_dir, comp_id, other_comps, verbose=False)
            if files_override:
                df = df[df["file"].isin(set(files_override))].copy()
            self._qc_checkstar_cache[comp_id] = df  # cache for preview plot

            if df.empty:
                rows.append({"comp_id": int(comp_id), "n": 0, "rms": np.nan, "mad": np.nan,
                             "sigma_nights": np.nan, "outliers": 0, "outlier_frac": np.nan,
                             "use": comp_id in active_set})
                continue

            y = pd.to_numeric(df["diff_mag_raw"], errors="coerce").to_numpy(float)
            night_ids = df["night_id"].to_numpy(int) if "night_id" in df.columns else np.zeros(len(df), int)
            m = np.isfinite(y)
            n = int(m.sum())
            if n <= 1:
                rows.append({"comp_id": int(comp_id), "n": n, "rms": np.nan, "mad": np.nan,
                             "sigma_nights": np.nan, "outliers": 0, "outlier_frac": np.nan,
                             "use": comp_id in active_set})
                continue

            yv = y[m]
            rms = float(np.nanstd(yv))
            med = float(np.nanmedian(yv))
            mad = float(np.nanmedian(np.abs(yv - med)))

            # Per-night medians → σ(nights) = night-to-night stability
            night_medians = []
            for nid in sorted(set(night_ids[m])):
                if nid <= 0:
                    continue
                nm = m & (night_ids == nid)
                if nm.sum() > 0:
                    night_medians.append(float(np.nanmedian(y[nm])))
            sigma_nights = float(np.std(night_medians, ddof=1)) if len(night_medians) >= 2 else np.nan

            # Outliers (3σ from median)
            robust_sigma = MAD_TO_SIGMA * mad
            outlier_count = int(
                np.sum(np.abs(yv - med) > self.qc_sigma * robust_sigma)
            ) if np.isfinite(robust_sigma) and robust_sigma > 0 else 0
            outlier_frac = outlier_count / max(n, 1)

            if verbose:
                sn_str = f"{sigma_nights:.4f}" if np.isfinite(sigma_nights) else "nan"
                self.log(f"[QC] Comp {comp_id}: RMS={rms:.4f} σ_nights={sn_str}")

            rows.append({
                "comp_id": int(comp_id),
                "n": n,
                "rms": rms,
                "mad": mad,
                "sigma_nights": sigma_nights,
                "night_medians": night_medians,
                "outliers": outlier_count,
                "outlier_frac": float(outlier_frac),
                "use": comp_id in active_set,
            })

        max_points = max((int(row.get("n", 0)) for row in rows), default=0)
        for row in rows:
            row["coverage_fraction"] = (
                float(int(row.get("n", 0)) / max_points) if max_points > 0 else 0.0
            )
        return rows

    def _save_comp_qc_summary(self, result_dir: Path, rows: list[dict]) -> None:
        if not rows:
            return
        out_dir = step9_lc_dir(result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        df = pd.DataFrame(rows)
        path = out_dir / "comp_qc_summary.csv"
        try:
            df.to_csv(path, index=False)
            self.log(f"[QC] Saved {path.name}")
        except Exception as e:
            self.log(f"[QC] Failed to save summary: {e}")

    def _build_light_curve_core(self, target_id: int, active_comp_ids: list[int]) -> dict:
        self.log("=" * 60)
        self.log("[BUILD] Starting Light Curve Build (RAW)")
        self.log(f"[BUILD] Target ID: {target_id}")
        self.log(f"[BUILD] Active Comp IDs: {active_comp_ids}")
        self.log(f"[BUILD] Datasets: {len(self.datasets)}")
        self._refresh_photometry_source_policy(log=True)

        if active_comp_ids and not self.runtime_mode:
            qc_rows = self._compute_comp_qc(self.datasets[0][1], target_id, active_comp_ids, verbose=False)
            self._save_comp_qc_summary(Path(self.datasets[0][1]), qc_rows)
        elif active_comp_ids and self.runtime_mode:
            self.log("[BUILD] Runtime mode: skip precomputing QC summary")

        combined_raw: list[pd.DataFrame] = []
        combined_check: list[pd.DataFrame] = []
        combined_check_ids_by_filter: dict[str, int] = {}
        single_dataset_mode = len(self.datasets) == 1
        for label, result_dir in self.datasets:
            result_dir = Path(result_dir)
            raw_df = self._build_ensemble_series(result_dir, target_id, active_comp_ids, verbose=True)
            if raw_df.empty:
                self.log(f"[{label}] Raw light curve empty")
                continue
            raw_df = annotate_raw_lightcurve(raw_df, label, logger=self.log)

            check_ids_by_filter = {}
            check_df = pd.DataFrame()
            try:
                check_ids_by_filter, check_df = self._build_check_star_series(
                    result_dir,
                    active_comp_ids,
                    verbose=False,
                )
                for filt_key, check_id in check_ids_by_filter.items():
                    combined_check_ids_by_filter.setdefault(filt_key, int(check_id))
                if not check_df.empty:
                    check_df = check_df.copy()
                    check_df["dataset"] = str(label)
                    combined_check.append(check_df)
                if check_ids_by_filter and check_df.empty:
                    self.log("  Check star configured but no usable check-star light curve was built")
            except Exception as e:
                self.log(f"  Check star export failed: {e}")

            save_dataset_raw_outputs(
                result_dir=result_dir,
                target_id=target_id,
                raw_df=raw_df,
                check_ids_by_filter=check_ids_by_filter,
                check_df=check_df,
                logger=self.log,
            )
            combined_raw.append(raw_df)

        save_combined_raw_outputs(
            base_result_dir=Path(self.params.P.result_dir),
            target_id=target_id,
            combined_raw=combined_raw,
            single_dataset_mode=single_dataset_mode,
            comp_candidate_ids=self.comp_candidate_ids,
            active_comp_ids=active_comp_ids,
            combined_check=combined_check,
            check_ids_by_filter=combined_check_ids_by_filter,
            logger=self.log,
        )

        self.log("=" * 60)
        self.log("[BUILD] Light Curve Build Complete (RAW)")
        summary = {
            "target_id": int(target_id),
            "n_datasets": len(self.datasets),
            "n_outputs": len(combined_raw),
            "n_total": 0,
            "n_valid": 0,
        }
        if combined_raw:
            all_data = pd.concat(combined_raw, ignore_index=True)
            valid_y = all_data["diff_mag_raw"].dropna()
            n_total = len(all_data)
            n_valid = len(valid_y)
            summary["n_total"] = int(n_total)
            summary["n_valid"] = int(n_valid)
            if n_valid > 0:
                y_mean = valid_y.mean()
                y_std = valid_y.std()
                y_range = valid_y.max() - valid_y.min()
                self.log(f"[RESULT] RAW: {n_valid}/{n_total} valid points")
                self.log(f"[RESULT] RAW: mean={y_mean:.4f}, std={y_std:.4f}, range={y_range:.4f} mag")
            else:
                self.log(f"[RESULT] RAW: 0/{n_total} valid points - CHECK DETECTION!")
        self.log("=" * 60)
        return summary


class HeadlessLightCurveBuilder(RawLightCurveBuilder):
    """A batch runner for the same build the window performs.

    The window sets its state up across a constructor and a dozen widget
    callbacks. A batch run has none of those, so this gathers the same fields in
    one place — and gathering them is the whole of it. There is no second
    implementation of the build here: `_build_light_curve_core` is inherited,
    which is the point of having moved it.
    """

    def __init__(self, params, result_dirs, *, logger=None, project_state=None,
                 file_manager=None, qc_sigma: float = 3.0,
                 comp_candidate_ids=None, runtime_mode: bool = False):
        from apex.utils.photometry_loader import FramePhotometryCache

        self.params = params
        self.project_state = project_state
        self.file_manager = file_manager
        self.runtime_mode = bool(runtime_mode)
        self.qc_sigma = float(qc_sigma)

        dirs = [Path(d) for d in result_dirs]
        self.datasets = [(d.name, d) for d in dirs]
        self.comp_candidate_ids = list(comp_candidate_ids or [])
        self.comp_ids_list = list(self.comp_candidate_ids)

        self._photometry_cache = FramePhotometryCache()
        self._photometry_source_cache: dict[str, dict] = {}
        self._header_cache: dict[str, fits.Header | None] = {}
        self._diff_series_cache: dict[tuple, pd.DataFrame] = {}
        self._check_series_cache: dict[tuple, tuple[dict, pd.DataFrame]] = {}
        self._qc_checkstar_cache: dict[int, pd.DataFrame] = {}
        self._frame_exclude_cache: dict[str, dict[str, set[str]]] = {}
        self._force_aperture_for_datasets = False

        self._logger = logger

    def log(self, message: str) -> None:
        if self._logger is not None:
            self._logger(message)
