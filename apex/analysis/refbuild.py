"""Qt-free master catalog build (headless port of RefBuildWorker.run()).

This module is a VERBATIM RELOCATION of the compute body of
``apex.gui.workflow.step6_ref_build.RefBuildWorker`` (its ``run`` / ``_run_impl``
and all the reference-selection / master-catalog helpers it calls) with only the
Qt coupling substituted by optional callbacks. The numerical algorithm, parameter
lookups, rounding, and file writes are unchanged.

LAYER RULE: this module must NOT import apex.gui or PyQt5.
"""

from __future__ import annotations

import json
import hashlib  # noqa: F401  (kept for parity with the GUI module imports)
import re
import time
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

from apex.utils.step_paths import (
    step5_wcs_dir,
    step6_refbuild_dir,
    step2_cropped_dir,
    step4_dir,
    crop_is_active,
    crop_rect_path,
)
from apex.utils.common_helpers import normalize_filter_key, safe_float as _safe_float
from apex.utils.io_utils import coerce_int64_source_id
from apex.utils.cache_utils import (
    norm_path_key,
    build_file_signature,
    cache_schema_value,
    detection_cache_signature_matches,
    file_signature_matches,
    file_signature_matches_relaxed,
    astap_wcs_candidates,
    parse_astap_wcs_file,
)


# ── Qt-free module helpers (relocated byte-for-byte from step6_ref_build.py) ─────

_FILTER_RE = re.compile(r"[-_]([ugrizbvUGRIZBV])[-_.]", re.IGNORECASE)
_DATE_RE = re.compile(r"(20\d{6})")
_GAIA_REF_EXTRA_COLS = (
    "pmra",
    "pmdec",
    "pmra_error",
    "pmdec_error",
    "parallax",
    "parallax_error",
    "ruwe",
    "visibility_periods_used",
)


def _get_filter_from_filename(filename: str) -> Optional[str]:
    match = _FILTER_RE.search(str(filename))
    return normalize_filter_key(match.group(1)) if match else None


def _parse_date_key(value: str, params) -> Optional[str]:
    mode = str(getattr(params.P, "night_parse_mode", "regex") or "regex").strip().lower()
    if mode == "split":
        delim = str(getattr(params.P, "night_parse_split_delim", "_"))
        parts = value.split(delim) if delim else [value]
        idx = int(getattr(params.P, "night_parse_split_index", -1))
        if idx < 0:
            idx = len(parts) + idx
        if idx < 0 or idx >= len(parts):
            return None
        return parts[idx]
    if mode == "last_digits":
        n_digits = max(1, int(getattr(params.P, "night_parse_last_digits", 8)))
        m = re.search(rf"(\\d{{{n_digits}}})$", value)
        return m.group(1) if m else None
    try:
        pattern = str(getattr(params.P, "night_parse_regex", r".*_(\d{8})"))
        m = re.search(pattern, value)
    except re.error:
        return None
    if not m:
        return None
    if m.groupdict().get("date"):
        return m.group("date")
    if m.groups():
        return m.group(1)
    return m.group(0)


def _extract_date_key(filename: str, params=None) -> str:
    if params is None or not hasattr(params, "P"):
        match = _DATE_RE.search(str(filename))
        return match.group(1) if match else "unknown_date"
    date_key = None
    try:
        data_dir = Path(getattr(params.P, "data_dir", "."))
        file_path = Path(params.get_file_path(filename))
        if file_path.parent != data_dir:
            date_key = _parse_date_key(file_path.parent.name, params)
        if not date_key:
            date_key = _parse_date_key(file_path.name, params)
    except Exception:
        date_key = None
    if not date_key:
        date_key = _parse_date_key(str(filename), params)
    return date_key or "unknown_date"


def run_refbuild(
    params,
    data_dir,
    result_dir,
    cache_dir,
    file_list,
    ref_filter,
    sat_drop_pct,
    elong_drop_pct,
    ref_cat_max_sources,
    ref_cat_min_sources,
    ref_cat_max_elong,
    ref_cat_max_abs_round,
    ref_cat_sharp_min,
    ref_cat_sharp_max,
    ref_cat_min_peak_adu,
    wcs_match_radius_arcsec,
    wcs_min_match_rate,
    wcs_min_match_n,
    wcs_max_sep_med_arcsec,
    wcs_max_sep_p90_arcsec,
    wcs_max_dup_rate,
    ref_per_date,
    ref_build_mode="hybrid",
    gaia_mag_limit=18.0,
    ref_master_union=True,
    ref_union_min_frames=1,
    *,
    progress_cb=None,
    log_cb=None,
    error_cb=None,
    should_stop=None,
    logger=None,
) -> dict:
    """Headless master catalog build over ``file_list``.

    Verbatim relocation of ``RefBuildWorker.run`` / ``_run_impl``: same algorithm,
    parameter lookups, rounding, and on-disk outputs. Qt signals are replaced by
    optional callbacks; the summary dict is returned instead of emitted.
    """
    # Bind call/state shims so the relocated body keeps using the same names.
    data_dir = Path(data_dir)
    result_dir = Path(result_dir)
    cache_dir = Path(cache_dir)
    file_list = list(file_list)
    sat_drop_pct = float(sat_drop_pct)
    elong_drop_pct = float(elong_drop_pct)
    ref_cat_max_sources = int(ref_cat_max_sources)
    ref_cat_min_sources = int(ref_cat_min_sources)
    ref_cat_max_elong = float(ref_cat_max_elong)
    ref_cat_max_abs_round = float(ref_cat_max_abs_round)
    ref_cat_sharp_min = float(ref_cat_sharp_min)
    ref_cat_sharp_max = float(ref_cat_sharp_max)
    ref_cat_min_peak_adu = float(ref_cat_min_peak_adu)
    wcs_match_radius_arcsec = float(wcs_match_radius_arcsec)
    wcs_min_match_rate = float(wcs_min_match_rate)
    wcs_min_match_n = int(wcs_min_match_n)
    wcs_max_sep_med_arcsec = float(wcs_max_sep_med_arcsec)
    wcs_max_sep_p90_arcsec = float(wcs_max_sep_p90_arcsec)
    wcs_max_dup_rate = float(wcs_max_dup_rate)
    ref_per_date = bool(ref_per_date)
    ref_build_mode = str(ref_build_mode).lower()
    gaia_mag_limit = float(gaia_mag_limit)
    ref_master_union = bool(ref_master_union)
    ref_union_min_frames = max(1, int(ref_union_min_frames))
    # Cache for WCS headers (path -> fits.Header)
    _wcs_header_cache: Dict[str, fits.Header] = {}
    _wcs_summary_by_file: Optional[Dict[str, dict]] = None
    _meta_missing_files: set = set()
    _meta_incompatible_files: set = set()
    _detcsv_missing_files: set = set()
    _detcsv_incompatible_files: set = set()

    def _stop_requested() -> bool:
        return should_stop() if should_stop else False

    def _log(msg: str):
        (log_cb(msg) if log_cb else None)

    def _resolve_fits_path(fname: str) -> Optional[Path]:
        cropped_dir = step2_cropped_dir(result_dir)
        if crop_is_active(result_dir):
            cand = cropped_dir / fname
            if cand.exists():
                return cand
        # Prefer original/source path first. Step5 summary signatures are based on
        # the true source frame, and stale copies can remain in older cache folders.
        try:
            orig = Path(params.get_file_path(fname))
            if orig.exists():
                return orig
        except Exception:
            pass
        step6_out = step5_wcs_dir(result_dir)
        cand = step6_out / fname
        if cand.exists():
            return cand
        try:
            return Path(params.get_file_path(fname))
        except Exception:
            return None

    def _current_file_signature(fname: str) -> Optional[dict]:
        path = _resolve_fits_path(fname)
        if path is None or not path.exists():
            return None
        return build_file_signature(path, use_cropped=bool(crop_is_active(result_dir)))

    def _detect_meta_compatible(fname: str, payload: dict, meta_path: Path) -> bool:
        if not isinstance(payload, dict):
            return False
        schema = cache_schema_value(payload)
        if schema < 2:
            return _schema1_detect_cache_allowed(meta_path)
        sig_now = _current_file_signature(fname)
        if sig_now is None:
            return False
        if detection_cache_signature_matches(
            payload,
            sig_now,
            min_schema=2,
            allow_mtime_drift=False,
        ):
            payload["__compat_relaxed_mtime"] = False
            payload["__compat_relaxed_size"] = False
            return True
        if detection_cache_signature_matches(
            payload,
            sig_now,
            min_schema=2,
            allow_mtime_drift=True,
        ):
            payload["__compat_relaxed_mtime"] = True
            payload["__compat_relaxed_size"] = True
            return True
        return False

    def _schema1_detect_cache_allowed(marker_path: Path) -> bool:
        try:
            marker_mtime = int(marker_path.stat().st_mtime_ns)
        except Exception:
            return False
        if crop_is_active(result_dir):
            rect_path = crop_rect_path(result_dir)
            if rect_path.exists():
                try:
                    rect_mtime = int(rect_path.stat().st_mtime_ns)
                    if marker_mtime < rect_mtime:
                        return False
                except Exception:
                    return False
        return True

    def _wcs_row_compatible(fname: str, row: dict) -> bool:
        sig_now = _current_file_signature(fname)
        if sig_now is None:
            return False
        if not isinstance(row, dict):
            return False
        src_path_val = row.get("source_path")
        src_path_txt = str(src_path_val).strip().lower() if src_path_val is not None else ""
        if src_path_txt not in ("", "nan", "none", "null"):
            if file_signature_matches(row, sig_now):
                row["__compat_relaxed_mtime"] = False
                row["__compat_relaxed_size"] = False
                return True
            if not file_signature_matches_relaxed(row, sig_now):
                return False
            row["__compat_relaxed_mtime"] = True
            row["__compat_relaxed_size"] = True
            return True

        # Schema-1 rows did not record full signatures; accept by path match only.
        fits_path = norm_path_key(row.get("fits_path", ""))
        if not fits_path:
            return False
        if fits_path != norm_path_key(sig_now.get("source_path", "")):
            return False
        row["__compat_relaxed_mtime"] = True
        return True

    def _load_meta(fname: str) -> Optional[dict]:
        candidates = [
            cache_dir / f"detect_{fname}.json",
            step4_dir(result_dir) / f"detect_{fname}.json",
        ]
        candidates = [p for p in candidates if p.exists()]
        if not candidates:
            _meta_missing_files.add(str(fname))
            return None
        candidates.sort(key=lambda p: p.stat().st_mtime_ns, reverse=True)
        for meta_path in candidates:
            try:
                payload = json.loads(meta_path.read_text(encoding="utf-8"))
            except Exception:
                continue
            if _detect_meta_compatible(fname, payload, meta_path):
                return payload
        _meta_incompatible_files.add(str(fname))
        return None

    def _resolve_detect_csv(fname: str) -> Optional[Path]:
        candidates = [
            cache_dir / f"detect_{fname}.csv",
            step4_dir(result_dir) / f"detect_{fname}.csv",
        ]
        candidates = [p for p in candidates if p.exists()]
        if not candidates:
            _detcsv_missing_files.add(str(fname))
            return None
        candidates.sort(key=lambda p: p.stat().st_mtime_ns, reverse=True)
        for cand in candidates:
            meta_path = cand.with_suffix(".json")
            if meta_path.exists():
                try:
                    payload = json.loads(meta_path.read_text(encoding="utf-8"))
                except Exception:
                    continue
                if _detect_meta_compatible(fname, payload, meta_path):
                    return cand
        _detcsv_incompatible_files.add(str(fname))
        return None

    def _load_wcs_summary_map() -> Dict[str, dict]:
        nonlocal _wcs_summary_by_file
        if _wcs_summary_by_file is not None:
            return _wcs_summary_by_file

        mapping: Dict[str, dict] = {}
        candidates = [
            step5_wcs_dir(result_dir) / "wcs_solve_summary.csv",
        ]
        existing = [p for p in candidates if p.exists()]
        existing.sort(key=lambda p: p.stat().st_mtime_ns, reverse=True)
        for path in existing:
            try:
                st = path.stat()
                rec_mtime_ns = int(st.st_mtime_ns)
            except Exception:
                rec_mtime_ns = 0
            try:
                df = pd.read_csv(path)
            except Exception:
                continue
            key_col = None
            if "file" in df.columns:
                key_col = "file"
            elif "fname" in df.columns:
                key_col = "fname"
            if key_col is None:
                continue
            for rec in df.to_dict(orient="records"):
                key = str(rec.get(key_col, "")).strip()
                if key:
                    rec["__record_mtime_ns"] = rec_mtime_ns
                    mapping[key] = rec
            break

        _wcs_summary_by_file = mapping
        return mapping

    def _load_wcs_summary_row(fname: str) -> dict:
        row = _load_wcs_summary_map().get(str(fname), {})
        if not row:
            return {}
        if not _wcs_row_compatible(fname, row):
            return {}
        return row

    def _load_wcs_for_frame(fname: str) -> Optional[WCS]:
        fits_path = _resolve_fits_path(fname)
        candidates = []
        if fits_path is not None and fits_path.exists():
            candidates.append(fits_path)
        # fallback: original file if cropped path lacks WCS
        try:
            orig = Path(params.get_file_path(fname))
            if orig.exists() and orig not in candidates:
                candidates.append(orig)
        except Exception:
            pass
        for path in candidates:
            try:
                # Use cached header if available
                path_key = str(path)
                if path_key in _wcs_header_cache:
                    hdr = _wcs_header_cache[path_key]
                else:
                    hdr = fits.getheader(path)
                    _wcs_header_cache[path_key] = hdr
                w = WCS(hdr, relax=True)
                if w.has_celestial:
                    return w
                # fallback: ASTAP .wcs sidecar
                for wcs_path in astap_wcs_candidates(path):
                    if not wcs_path.exists():
                        continue
                    wcs_dict = parse_astap_wcs_file(wcs_path)
                    if not wcs_dict:
                        continue
                    hdr2 = fits.Header()
                    for k, v in wcs_dict.items():
                        try:
                            hdr2[k] = v
                        except Exception:
                            pass
                    w2 = WCS(hdr2, relax=True)
                    if w2.has_celestial:
                        return w2
            except Exception:
                continue
        return None

    def _load_gaia_catalog() -> Optional[SkyCoord]:
        try:
            df = _load_gaia_table()
            if df is None or df.empty:
                return None
            ra = pd.to_numeric(df["ra"], errors="coerce")
            dec = pd.to_numeric(df["dec"], errors="coerce")
            m = ra.notna() & dec.notna()
            if not m.any():
                return None
            return SkyCoord(ra[m].to_numpy(float) * u.deg, dec[m].to_numpy(float) * u.deg, frame="icrs")
        except Exception:
            return None

    def _load_gaia_table() -> Optional[pd.DataFrame]:
        gaia_path = step5_wcs_dir(result_dir) / "gaia_fov.ecsv"
        if not gaia_path.exists():
            return None
        try:
            tab = Table.read(str(gaia_path), format="ascii.ecsv")
            cols = [c.lower() for c in tab.colnames]
            if cols != list(tab.colnames):
                tab.rename_columns(tab.colnames, cols)
            df = tab.to_pandas()
            if "ra" not in df.columns or "dec" not in df.columns:
                return None
            return df
        except Exception:
            return None

    def _attach_gaia_photometry(master_df: pd.DataFrame, gaia_df: Optional[pd.DataFrame]) -> pd.DataFrame:
        if gaia_df is None or gaia_df.empty:
            return master_df
        if "ra_deg" not in master_df.columns or "dec_deg" not in master_df.columns:
            return master_df

        gaia_ra = pd.to_numeric(gaia_df.get("ra"), errors="coerce")
        gaia_dec = pd.to_numeric(gaia_df.get("dec"), errors="coerce")
        mask = gaia_ra.notna() & gaia_dec.notna()
        if not mask.any():
            return master_df

        gaia_df = gaia_df.loc[mask].copy()
        gaia_ra = gaia_ra[mask].to_numpy(float)
        gaia_dec = gaia_dec[mask].to_numpy(float)
        gaia_sky = SkyCoord(gaia_ra * u.deg, gaia_dec * u.deg, frame="icrs")

        src_ra = pd.to_numeric(master_df["ra_deg"], errors="coerce")
        src_dec = pd.to_numeric(master_df["dec_deg"], errors="coerce")
        src_mask = src_ra.notna() & src_dec.notna()
        if not src_mask.any():
            return master_df

        src_sky = SkyCoord(src_ra[src_mask].to_numpy(float) * u.deg,
                           src_dec[src_mask].to_numpy(float) * u.deg,
                           frame="icrs")
        idx, sep2d, _ = src_sky.match_to_catalog_sky(gaia_sky)
        match_r = max(0.5, float(wcs_match_radius_arcsec))
        ok = sep2d.arcsec <= match_r
        if not np.any(ok):
            return master_df

        out = master_df.copy()
        # Keep Gaia source_id as nullable Int64 to avoid float64 precision loss
        # on 64-bit IDs (which can introduce +/-256 rounding drift).
        out["gaia_source_id"] = pd.Series(pd.array([pd.NA] * len(out), dtype="Int64"), index=out.index)
        out["gaia_ra_deg"] = np.nan
        out["gaia_dec_deg"] = np.nan
        out["phot_g_mean_mag"] = np.nan
        out["phot_bp_mean_mag"] = np.nan
        out["phot_rp_mean_mag"] = np.nan
        for col in _GAIA_REF_EXTRA_COLS:
            if col in gaia_df.columns:
                out[col] = np.nan

        src_idx = np.where(src_mask)[0]
        match_idx = src_idx[ok]
        gaia_idx = idx[ok]

        if "source_id" in gaia_df.columns:
            gaia_sid = coerce_int64_source_id(gaia_df["source_id"])
            out.loc[match_idx, "gaia_source_id"] = pd.array(
                gaia_sid.iloc[np.asarray(gaia_idx, dtype=int)].tolist(),
                dtype="Int64",
            )

        out.loc[match_idx, "gaia_ra_deg"] = gaia_ra[gaia_idx]
        out.loc[match_idx, "gaia_dec_deg"] = gaia_dec[gaia_idx]

        if "phot_g_mean_mag" in gaia_df.columns:
            g = pd.to_numeric(gaia_df["phot_g_mean_mag"], errors="coerce").to_numpy()
            out.loc[match_idx, "phot_g_mean_mag"] = g[gaia_idx]
        if "phot_bp_mean_mag" in gaia_df.columns:
            bp = pd.to_numeric(gaia_df["phot_bp_mean_mag"], errors="coerce").to_numpy()
            out.loc[match_idx, "phot_bp_mean_mag"] = bp[gaia_idx]
        if "phot_rp_mean_mag" in gaia_df.columns:
            rp = pd.to_numeric(gaia_df["phot_rp_mean_mag"], errors="coerce").to_numpy()
            out.loc[match_idx, "phot_rp_mean_mag"] = rp[gaia_idx]
        for col in _GAIA_REF_EXTRA_COLS:
            if col in gaia_df.columns:
                vals = pd.to_numeric(gaia_df[col], errors="coerce").to_numpy()
                out.loc[match_idx, col] = vals[gaia_idx]

        out["gaia_G"] = out["phot_g_mean_mag"]
        out["gaia_BP"] = out["phot_bp_mean_mag"]
        out["gaia_RP"] = out["phot_rp_mean_mag"]
        try:
            out["color_gr"] = pd.to_numeric(out["gaia_BP"], errors="coerce") - pd.to_numeric(out["gaia_RP"], errors="coerce")
        except Exception:
            pass
        return out

    def _apply_hybrid_source_ids(df: pd.DataFrame, gaia_mag_limit: float = 18.0) -> tuple[pd.DataFrame, dict[int, int], dict[int, int]]:
        """Apply hybrid source_id assignment: Gaia ID for matched sources, negative ID for non-Gaia.

        In hybrid mode:
        - Sources matched to Gaia: use gaia_source_id (positive, from Gaia DR3)
        - Sources not matched: assign negative local IDs (-1, -2, ...)

        This ensures consistent source_id across all frames for Gaia-matched sources.
        """
        out = df.copy()
        old_ids = coerce_int64_source_id(out["source_id"]) if "source_id" in out.columns else None

        # Check if gaia_source_id column exists
        if "gaia_source_id" not in out.columns:
            _log("[REF] No gaia_source_id column; hybrid mode not applied.")
            return out, {}, {}

        # Filter by magnitude limit if gaia_G is available
        n_trimmed = 0
        gaia_sid = coerce_int64_source_id(out["gaia_source_id"])
        if "gaia_G" in out.columns and gaia_mag_limit > 0:
            gaia_g = pd.to_numeric(out["gaia_G"], errors="coerce")
            too_faint = gaia_g > gaia_mag_limit
            n_trimmed = int((too_faint & gaia_sid.notna() & (gaia_sid > 0)).sum())
            # Clear Gaia ID for sources fainter than limit
            out.loc[too_faint, "gaia_source_id"] = pd.NA
            gaia_sid = coerce_int64_source_id(out["gaia_source_id"])

        # Identify sources with valid Gaia source_id
        has_gaia = gaia_sid.notna() & (gaia_sid > 0)

        # Create new source_id column
        new_source_id = pd.Series(pd.array([pd.NA] * len(out), dtype="Int64"), index=out.index)
        next_local_id = -1

        for i in out.index:
            if bool(has_gaia.loc[i]):
                # Use Gaia source_id (positive)
                new_source_id.loc[i] = int(gaia_sid.loc[i])
            else:
                # Assign negative local ID
                new_source_id.loc[i] = next_local_id
                next_local_id -= 1

        out["source_id"] = new_source_id.astype("Int64")
        # Keep ID as sequential for display purposes
        out["ID"] = range(1, len(out) + 1)

        # Build mapping from old source_id to new source_id and display ID
        sid_map = {}
        id_map = {}
        if old_ids is not None:
            new_vals = coerce_int64_source_id(out["source_id"])
            id_vals = pd.to_numeric(out["ID"], errors="coerce").astype("Int64")
            for o, n, i in zip(old_ids, new_vals, id_vals):
                if pd.isna(o):
                    continue
                old_i = int(o)
                sid_map[old_i] = int(n) if pd.notna(n) else old_i
                id_map[old_i] = int(i) if pd.notna(i) else old_i

        n_gaia = int(has_gaia.sum())
        n_local = len(out) - n_gaia
        if n_trimmed > 0:
            _log(f"[REF] Gaia IDs excluded by mag limit (G>{gaia_mag_limit:.2f}): {n_trimmed}")
        _log(f"[REF] Hybrid IDs assigned: {n_gaia} Gaia, {n_local} local (negative)")

        return out, sid_map, id_map

    def _merge_ref_catalogs(
        base_df: Optional[pd.DataFrame],
        new_df: pd.DataFrame,
        match_radius_arcsec: float,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        if base_df is None or base_df.empty:
            out = new_df.copy()
            out = out.reset_index(drop=True)
            out["source_id"] = np.arange(1, len(out) + 1, dtype=int)
            out["ID"] = out["source_id"]
            return out, out

        base = base_df.copy().reset_index(drop=True)
        new = new_df.copy().reset_index(drop=True)

        base_ra = pd.to_numeric(base.get("ra_deg"), errors="coerce")
        base_dec = pd.to_numeric(base.get("dec_deg"), errors="coerce")
        new_ra = pd.to_numeric(new.get("ra_deg"), errors="coerce")
        new_dec = pd.to_numeric(new.get("dec_deg"), errors="coerce")
        base_mask = base_ra.notna() & base_dec.notna()
        new_mask = new_ra.notna() & new_dec.notna()
        if not base_mask.any() or not new_mask.any():
            return base, new

        base_sky = SkyCoord(base_ra[base_mask].to_numpy(float) * u.deg,
                            base_dec[base_mask].to_numpy(float) * u.deg,
                            frame="icrs")
        new_sky = SkyCoord(new_ra[new_mask].to_numpy(float) * u.deg,
                           new_dec[new_mask].to_numpy(float) * u.deg,
                           frame="icrs")
        idx, sep2d, _ = new_sky.match_to_catalog_sky(base_sky)
        ok = sep2d.arcsec <= match_radius_arcsec

        new["source_id"] = np.nan
        new["ID"] = np.nan

        base_ids = coerce_int64_source_id(base.loc[base_mask, "source_id"]).to_numpy(dtype=np.int64, na_value=0)
        match_idx = np.where(new_mask)[0]
        ok_idx = match_idx[ok]
        if len(ok_idx):
            new.loc[ok_idx, "source_id"] = base_ids[idx[ok]]
            new.loc[ok_idx, "ID"] = base_ids[idx[ok]]

        base_sid = coerce_int64_source_id(base["source_id"]).dropna()
        next_id = int(base_sid.max() if not base_sid.empty else 0) + 1
        new_rows = []
        for i in match_idx[~ok]:
            sid = next_id
            next_id += 1
            new.loc[i, "source_id"] = sid
            new.loc[i, "ID"] = sid
            new_rows.append(new.loc[i:i])

        if new_rows:
            base = pd.concat([base] + new_rows, ignore_index=True)

        new["source_id"] = coerce_int64_source_id(new["source_id"])
        new["ID"] = pd.to_numeric(new["ID"], errors="coerce").astype("Int64")
        return base, new

    def _load_wcs_meta(fname: str) -> dict:
        meta_path = Path(cache_dir) / "wcs_solve" / f"wcs_{fname}.json"
        if not meta_path.exists():
            return {}
        try:
            return json.loads(meta_path.read_text(encoding="utf-8"))
        except Exception:
            return {}

    def _combined_wcs_meta(fname: str) -> dict:
        summary_row = _load_wcs_summary_row(fname)
        cache_meta = _load_wcs_meta(fname)
        merged = {}
        if summary_row:
            merged.update(summary_row)
        if cache_meta:
            merged.update(cache_meta)
        return merged

    def _wcs_meta_match_stats(meta: dict) -> dict:
        def _num(key: str, default=np.nan):
            return _safe_float(meta.get(key), default)

        n_match = meta.get("n_match", meta.get("match_n", 0))
        try:
            n_match = int(float(n_match)) if pd.notna(n_match) else 0
        except Exception:
            n_match = 0

        n_cat = meta.get("n_catalog_in_fov", 0)
        try:
            n_cat = int(float(n_cat)) if pd.notna(n_cat) else 0
        except Exception:
            n_cat = 0

        return {
            "n_match": n_match,
            "n_catalog_in_fov": n_cat,
            "match_rate": _num("match_rate"),
            "match_rate_cat": _num("match_rate_cat"),
            "match_rate_eff": _num("match_rate_eff"),
            "sep_med_arcsec": np.nan,
            "sep_p90_arcsec": np.nan,
            "dup_rate": np.nan,
        }

    def _compute_match_stats(det_xy: np.ndarray, wcs: WCS, gaia_sky: SkyCoord) -> dict:
        n_det = len(det_xy)
        n_cat = int(len(gaia_sky)) if gaia_sky is not None else 0
        if n_det == 0:
            return dict(
                n_det=0,
                n_catalog_in_fov=n_cat,
                n_match=0,
                match_rate=0.0,
                match_rate_cat=0.0,
                match_rate_eff=0.0,
                sep_med_arcsec=np.nan,
                sep_p90_arcsec=np.nan,
                dup_rate=np.nan,
            )

        try:
            ra, dec = wcs.all_pix2world(det_xy[:, 0], det_xy[:, 1], 0)
            det_sky = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
        except Exception:
            return dict(
                n_det=n_det,
                n_catalog_in_fov=n_cat,
                n_match=0,
                match_rate=0.0,
                match_rate_cat=0.0,
                match_rate_eff=0.0,
                sep_med_arcsec=np.nan,
                sep_p90_arcsec=np.nan,
                dup_rate=np.nan,
            )

        idx, sep2d, _ = det_sky.match_to_catalog_sky(gaia_sky)
        sep_arcsec = sep2d.arcsec
        r_match = max(0.5, float(wcs_match_radius_arcsec))
        ok = sep_arcsec <= r_match
        n_match = int(np.sum(ok))
        if n_match == 0:
            return dict(
                n_det=n_det,
                n_catalog_in_fov=n_cat,
                n_match=0,
                match_rate=0.0,
                match_rate_cat=0.0,
                match_rate_eff=0.0,
                sep_med_arcsec=np.nan,
                sep_p90_arcsec=np.nan,
                dup_rate=np.nan,
            )

        sep_ok = sep_arcsec[ok]
        sep_med = float(np.nanmedian(sep_ok)) if len(sep_ok) else np.nan
        sep_p90 = float(np.nanpercentile(sep_ok, 90)) if len(sep_ok) else np.nan

        # duplicate matches: multiple detections -> same Gaia index
        dup_rate = np.nan
        try:
            counts = pd.Series(idx[ok]).value_counts()
            dup = counts[counts > 1].sum()
            dup_rate = float(dup / max(n_match, 1))
        except Exception:
            dup_rate = np.nan

        return dict(
            n_det=n_det,
            n_catalog_in_fov=n_cat,
            n_match=n_match,
            match_rate=float(n_match / max(n_det, 1)),
            match_rate_cat=float(n_match / max(n_cat, 1)),
            match_rate_eff=float(max(n_match / max(n_det, 1), n_match / max(n_cat, 1))),
            sep_med_arcsec=sep_med,
            sep_p90_arcsec=sep_p90,
            dup_rate=dup_rate,
        )

    def _frame_metrics(fname: str) -> Optional[dict]:
        meta = _load_meta(fname)
        if not meta:
            return None
        filt = normalize_filter_key(meta.get("filter", ""))
        if not filt:
            filt = _get_filter_from_filename(fname) or "unknown"
        fwhm_px = _safe_float(meta.get("fwhm_px"), np.nan)
        n_sources = int(meta.get("n_sources", 0) or 0)
        sat_count = int(meta.get("sat_star_count", 0) or 0)
        med_elong = _safe_float(meta.get("median_elongation"), np.nan)
        med_round = _safe_float(meta.get("median_roundness"), np.nan)
        shape_metric = med_elong
        if not np.isfinite(shape_metric):
            shape_metric = abs(med_round) if np.isfinite(med_round) else np.nan

        # Sky background statistics (from detection metadata)
        sky_med = _safe_float(meta.get("sky_med"), np.nan)
        sky_sigma = _safe_float(meta.get("sky_sigma"), np.nan)
        # Alternative names that might be used
        if not np.isfinite(sky_med):
            sky_med = _safe_float(meta.get("bkg_median"), np.nan)
        if not np.isfinite(sky_sigma):
            sky_sigma = _safe_float(meta.get("bkg_rms"), np.nan)

        return {
            "file": fname,
            "filter": filt,
            "date_key": _extract_date_key(fname, params),
            "fwhm_px": fwhm_px,
            "n_sources": n_sources,
            "sat_star_count": sat_count,
            "shape_metric": shape_metric,
            "sky_med": sky_med,
            "sky_sigma": sky_sigma,
        }

    def _select_reference(metrics: pd.DataFrame, ref_filter: str) -> str:
        if metrics.empty:
            raise RuntimeError("No detection metrics available")

        df = metrics.copy()
        filt = normalize_filter_key(ref_filter)
        if filt:
            cand = df[df["filter"] == filt].copy()
            if cand.empty:
                _log(f"[REF][QC] Filter '{filt}' not found; using all filters.")
                cand = df.copy()
        else:
            cand = df.copy()

        def _drop_top_percent(sub: pd.DataFrame, col: str, pct: float) -> pd.DataFrame:
            if sub.empty or pct <= 0:
                return sub
            vals = pd.to_numeric(sub[col], errors="coerce")
            if vals.notna().sum() == 0:
                return sub
            n = len(sub)
            drop_n = int(np.ceil(n * pct / 100.0))
            drop_n = max(0, min(n - 1, drop_n))
            if drop_n == 0:
                return sub
            sub = sub.copy()
            sub["_metric"] = vals
            sub = sub.sort_values("_metric", ascending=False)
            kept = sub.iloc[drop_n:]
            return kept.drop(columns=["_metric"])

        n_start = len(cand)
        cand = _drop_top_percent(cand, "sat_star_count", sat_drop_pct)
        n_sat = len(cand)
        cand = _drop_top_percent(cand, "shape_metric", elong_drop_pct)
        n_shape = len(cand)
        if cand.empty:
            _log("[REF][QC] All candidates dropped by sat/shape filters; using full set.")
            cand = df.copy()

        cand = cand.copy()
        cand["fwhm_px"] = pd.to_numeric(cand["fwhm_px"], errors="coerce")
        cand["shape_metric"] = pd.to_numeric(cand["shape_metric"], errors="coerce")
        cand["sat_star_count"] = pd.to_numeric(cand["sat_star_count"], errors="coerce")
        cand["n_sources"] = pd.to_numeric(cand["n_sources"], errors="coerce")

        # Apply WCS match quality filters when available
        if "match_rate" in cand.columns:
            cand["match_rate"] = pd.to_numeric(cand["match_rate"], errors="coerce")
        if "match_rate_cat" in cand.columns:
            cand["match_rate_cat"] = pd.to_numeric(cand["match_rate_cat"], errors="coerce")
        if "match_rate_eff" in cand.columns:
            cand["match_rate_eff"] = pd.to_numeric(cand["match_rate_eff"], errors="coerce")
        rate_cols = [c for c in ("match_rate", "match_rate_cat", "match_rate_eff") if c in cand.columns]
        if rate_cols:
            cand["match_rate_eff"] = pd.concat(
                [cand[c] for c in rate_cols], axis=1
            ).max(axis=1)
        if "n_match" in cand.columns:
            cand["n_match"] = pd.to_numeric(cand["n_match"], errors="coerce")
        if "sep_med_arcsec" in cand.columns:
            cand["sep_med_arcsec"] = pd.to_numeric(cand["sep_med_arcsec"], errors="coerce")
        if "sep_p90_arcsec" in cand.columns:
            cand["sep_p90_arcsec"] = pd.to_numeric(cand["sep_p90_arcsec"], errors="coerce")
        if "dup_rate" in cand.columns:
            cand["dup_rate"] = pd.to_numeric(cand["dup_rate"], errors="coerce")

        cand_wcs = cand.copy()
        if "wcs_ok" in cand_wcs.columns:
            cand_wcs = cand_wcs[cand_wcs["wcs_ok"] == True]
        if "match_rate_eff" in cand_wcs.columns and wcs_min_match_rate > 0:
            if cand_wcs["match_rate_eff"].notna().any():
                cand_wcs = cand_wcs[cand_wcs["match_rate_eff"] >= wcs_min_match_rate]
        if "n_match" in cand_wcs.columns and wcs_min_match_n > 0:
            if cand_wcs["n_match"].notna().any():
                cand_wcs = cand_wcs[cand_wcs["n_match"] >= wcs_min_match_n]
        if "sep_med_arcsec" in cand_wcs.columns and np.isfinite(wcs_max_sep_med_arcsec) and wcs_max_sep_med_arcsec > 0:
            if cand_wcs["sep_med_arcsec"].notna().any():
                cand_wcs = cand_wcs[cand_wcs["sep_med_arcsec"] <= wcs_max_sep_med_arcsec]
        if "sep_p90_arcsec" in cand_wcs.columns and np.isfinite(wcs_max_sep_p90_arcsec) and wcs_max_sep_p90_arcsec > 0:
            if cand_wcs["sep_p90_arcsec"].notna().any():
                cand_wcs = cand_wcs[cand_wcs["sep_p90_arcsec"] <= wcs_max_sep_p90_arcsec]
        if "dup_rate" in cand_wcs.columns and np.isfinite(wcs_max_dup_rate) and wcs_max_dup_rate > 0:
            if cand_wcs["dup_rate"].notna().any():
                cand_wcs = cand_wcs[cand_wcs["dup_rate"] <= wcs_max_dup_rate]

        n_wcs = len(cand_wcs)
        if not cand_wcs.empty:
            cand = cand_wcs
        _log(
            "[REF][QC] candidates: start={s} sat_drop={sat:.1f}% -> {n1} "
            "shape_drop={elong:.1f}% -> {n2} wcs_pass -> {n3}".format(
                s=n_start,
                sat=sat_drop_pct,
                n1=n_sat,
                elong=elong_drop_pct,
                n2=n_shape,
                n3=n_wcs,
            )
        )
        _log(
            "[REF][QC] wcs thresholds: match_r={r:.2f}\" min_rate={mr:.2f} "
            "min_match={mn} max_sep_med={smed:.2f}\" max_sep_p90={sp90:.2f}\" max_dup={dup:.2f}".format(
                r=wcs_match_radius_arcsec,
                mr=wcs_min_match_rate,
                mn=wcs_min_match_n,
                smed=wcs_max_sep_med_arcsec,
                sp90=wcs_max_sep_p90_arcsec,
                dup=wcs_max_dup_rate,
            )
        )

        sort_cols = []
        sort_asc = []
        if "match_rate_eff" in cand.columns:
            sort_cols.append("match_rate_eff")
            sort_asc.append(False)
        if "sep_med_arcsec" in cand.columns:
            sort_cols.append("sep_med_arcsec")
            sort_asc.append(True)
        sort_cols.extend(["fwhm_px", "shape_metric", "sat_star_count", "n_sources"])
        sort_asc.extend([True, True, True, False])

        cand = cand.sort_values(
            sort_cols,
            ascending=sort_asc,
            na_position="last",
        )
        try:
            top = cand.head(5)
            for _, row in top.iterrows():
                _log(
                    "[REF][QC] cand {file} f={flt} mr_eff={mre} mr_det={mr} mr_cat={mrc} sep={sep} fwhm={fwhm} "
                    "shape={shape} sat={sat} n={n} wcs_ok={wcs} match_n={mn}".format(
                        file=row.get("file", ""),
                        flt=row.get("filter", ""),
                        mre=_safe_float(row.get("match_rate_eff"), np.nan),
                        mr=_safe_float(row.get("match_rate"), np.nan),
                        mrc=_safe_float(row.get("match_rate_cat"), np.nan),
                        sep=_safe_float(row.get("sep_med_arcsec"), np.nan),
                        fwhm=_safe_float(row.get("fwhm_px"), np.nan),
                        shape=_safe_float(row.get("shape_metric"), np.nan),
                        sat=_safe_float(row.get("sat_star_count"), np.nan),
                        n=_safe_float(row.get("n_sources"), np.nan),
                        wcs=bool(row.get("wcs_ok", False)),
                        mn=int(row.get("n_match", 0) or 0),
                    )
                )
        except Exception:
            pass
        return str(cand.iloc[0]["file"])

    def _build_master_catalog(ref_fname: str) -> tuple[pd.DataFrame, dict]:
        """Build master catalog from reference frame detections.

        Returns:
            Tuple of (catalog DataFrame, stats dict with n_ref_total/after_cuts/used)
        """
        det_path = _resolve_detect_csv(ref_fname)
        if det_path is None:
            raise RuntimeError(f"Missing detection file: detect_{ref_fname}.csv")
        df = pd.read_csv(det_path)
        if not {"x", "y"} <= set(df.columns):
            raise RuntimeError(f"Detection file missing x/y: {det_path}")

        df = df.copy()
        for col in (
            "x",
            "y",
            "elongation",
            "roundness",
            "sharpness",
            "dao_flux",
            "dao_peak",
            "peak_adu",
            "fwhm_px",
        ):
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        base_all = df[["x", "y"]].rename(columns={"x": "x_ref", "y": "y_ref"})
        base_all = base_all.dropna(subset=["x_ref", "y_ref"]).copy()
        if base_all.empty:
            raise RuntimeError("Reference detection list is empty")

        # Track n_ref_total
        n_ref_total = len(base_all)

        cand = df.dropna(subset=["x", "y"]).copy()
        n_before = len(cand)

        if np.isfinite(ref_cat_max_elong) and ref_cat_max_elong > 0 and "elongation" in cand.columns:
            cand = cand[cand["elongation"] <= ref_cat_max_elong]

        if np.isfinite(ref_cat_max_abs_round) and ref_cat_max_abs_round > 0 and "roundness" in cand.columns:
            cand = cand[cand["roundness"].abs() <= ref_cat_max_abs_round]

        if "sharpness" in cand.columns:
            if np.isfinite(ref_cat_sharp_min):
                cand = cand[cand["sharpness"] >= ref_cat_sharp_min]
            if np.isfinite(ref_cat_sharp_max) and ref_cat_sharp_max > 0:
                cand = cand[cand["sharpness"] <= ref_cat_sharp_max]

        # Track n_ref_after_qualitycuts (after shape/quality filters, before brightness)
        n_ref_after_qualitycuts = len(cand)

        brightness = None
        for col in ("dao_flux", "dao_peak", "peak_adu"):
            if col in cand.columns:
                brightness = cand[col]
                break
        if brightness is not None:
            cand = cand.copy()
            cand["brightness"] = brightness
            if np.isfinite(ref_cat_min_peak_adu) and ref_cat_min_peak_adu > 0:
                cand = cand[cand["brightness"] >= ref_cat_min_peak_adu]
            if ref_cat_max_sources > 0:
                cand = cand.sort_values("brightness", ascending=False)
                cand = cand.head(ref_cat_max_sources)

        base = cand[["x", "y"]].rename(columns={"x": "x_ref", "y": "y_ref"})
        base = base.dropna(subset=["x_ref", "y_ref"]).copy()

        used_full_detections = False
        if base.empty or len(base) < max(10, ref_cat_min_sources):
            _log(
                f"[REF] Ref catalog filter too strict "
                f"({len(base)}/{n_before}); using full detections."
            )
            base = base_all.copy()
            used_full_detections = True
        else:
            _log(
                f"[REF] Ref catalog filtered: {len(base)}/{n_before} sources kept."
            )

        # Track n_ref_used (final count)
        n_ref_used = len(base)

        wcs = _load_wcs_for_frame(ref_fname)
        if wcs is None:
            raise RuntimeError(f"No WCS for reference frame: {ref_fname}")
        try:
            ra, dec = wcs.all_pix2world(base["x_ref"].to_numpy(float), base["y_ref"].to_numpy(float), 0)
        except Exception as e:
            raise RuntimeError(f"WCS conversion failed for {ref_fname}: {e}")

        base = base.sort_values(["y_ref", "x_ref"]).reset_index(drop=True)
        base["source_id"] = range(1, len(base) + 1)
        base["ID"] = range(1, len(base) + 1)
        base["ra_deg"] = ra
        base["dec_deg"] = dec

        # Log ref catalog stats
        _log(
            f"[REF][STATS] n_ref_total={n_ref_total} n_ref_after_qualitycuts={n_ref_after_qualitycuts} "
            f"n_ref_used={n_ref_used} used_full={used_full_detections}"
        )

        ref_stats = {
            "n_ref_total": n_ref_total,
            "n_ref_after_qualitycuts": n_ref_after_qualitycuts,
            "n_ref_used": n_ref_used,
            "used_full_detections": used_full_detections,
        }

        return base[["ID", "source_id", "ra_deg", "dec_deg", "x_ref", "y_ref"]].copy(), ref_stats

    def _build_master_for_group(
        group_metrics: pd.DataFrame, anchor_fname: str, match_radius_arcsec: float
    ) -> tuple[pd.DataFrame, dict]:
        """Build the master catalog for a group of frames.

        With ``ref_master_union`` enabled, detections from *all* frames in the
        group are projected to sky coordinates and merged (deduped by position)
        so that stars only reachable in some frames — faint stars in deep
        exposures, and bright stars that saturate in long exposures but are
        clean in short ones — all enter the master. Otherwise the catalog is
        built from the single anchor (reference) frame only.
        """
        if not ref_master_union:
            return _build_master_catalog(anchor_fname)
        return _build_union_master(group_metrics, anchor_fname, match_radius_arcsec)

    def _build_union_master(
        group_metrics: pd.DataFrame, anchor_fname: str, match_radius_arcsec: float
    ) -> tuple[pd.DataFrame, dict]:
        """Union per-frame detections into one position-deduped master catalog.

        Each frame's detections (after the same per-frame quality cuts as the
        single-frame build) are matched to the accumulating master by sky
        position; matches reuse the existing source_id, new sources are added.
        ``n_det_frames`` records how many frames each source was detected in.
        """
        files = [str(f) for f in group_metrics["file"].tolist()]
        # Anchor first: its detections seed the source_ids and the union x_ref/
        # y_ref share the anchor's pixel system.
        ordered = [anchor_fname] + [f for f in files if f != anchor_fname]

        master: Optional[pd.DataFrame] = None
        det_counts: Counter = Counter()
        last_stats: dict = {}
        n_used = 0
        for fname in ordered:
            if _stop_requested():
                break
            try:
                cat, stats = _build_master_catalog(fname)
            except RuntimeError as e:
                _log(f"[REF][UNION] skip {fname}: {e}")
                continue
            last_stats = stats
            n_used += 1
            master, merged = _merge_ref_catalogs(master, cat, match_radius_arcsec)
            sids = coerce_int64_source_id(merged.get("source_id")).dropna()
            if len(sids):
                det_counts.update(int(s) for s in sids.to_numpy(dtype=np.int64))

        if master is None or master.empty:
            _log("[REF][UNION] union produced no sources; falling back to anchor frame.")
            return _build_master_catalog(anchor_fname)

        master = master.copy().reset_index(drop=True)
        sid_int = coerce_int64_source_id(master["source_id"])
        master["n_det_frames"] = [
            int(det_counts.get(int(s), 0)) if pd.notna(s) else 0 for s in sid_int
        ]

        if ref_union_min_frames > 1:
            before = len(master)
            master = master[master["n_det_frames"] >= ref_union_min_frames].copy()
            _log(
                f"[REF][UNION] min_frames>={ref_union_min_frames}: "
                f"{before}->{len(master)} sources"
            )

        # Express x_ref/y_ref in a single (anchor) pixel system so neighbor /
        # crowding geometry is consistent across stars added from different frames.
        anchor_wcs = _load_wcs_for_frame(anchor_fname)
        if anchor_wcs is not None and {"ra_deg", "dec_deg"} <= set(master.columns):
            try:
                xr, yr = anchor_wcs.all_world2pix(
                    master["ra_deg"].to_numpy(float),
                    master["dec_deg"].to_numpy(float),
                    0,
                )
                master["x_ref"] = xr
                master["y_ref"] = yr
            except Exception as e:
                _log(f"[REF][UNION] anchor reprojection failed: {e}")

        # Dense, position-ordered IDs (matches the single-frame build contract).
        sort_cols = [c for c in ("y_ref", "x_ref") if c in master.columns]
        if sort_cols:
            master = master.sort_values(sort_cols).reset_index(drop=True)
        master["source_id"] = np.arange(1, len(master) + 1, dtype=int)
        master["ID"] = master["source_id"]

        stats = dict(last_stats)
        stats["n_master_union"] = int(len(master))
        stats["n_union_frames"] = int(n_used)
        _log(
            f"[REF][UNION] {n_used} frames -> {len(master)} master sources "
            f"(anchor={anchor_fname}, match_r={match_radius_arcsec:.2f}\")"
        )
        return master, stats

    # ── Orchestration (verbatim relocation of RefBuildWorker._run_impl) ──────────

    if not file_list:
        raise RuntimeError("No frames available")

    _log(f"[REF] Gaia hybrid ID mag limit: G<={gaia_mag_limit:.2f}")

    metrics_rows = []
    total = len(file_list)
    for i, fname in enumerate(file_list, 1):
        if _stop_requested():
            return
        row = _frame_metrics(fname)
        if row:
            metrics_rows.append(row)
        (progress_cb(i, total, fname) if progress_cb else None)

    if not metrics_rows:
        raise RuntimeError("No detection metrics found. Run Source Detection first.")

    metrics = pd.DataFrame(metrics_rows)

    gaia_df = _load_gaia_table()
    gaia_sky = _load_gaia_catalog()
    if gaia_sky is None:
        _log("[REF] Gaia catalog not available; WCS match stats will be skipped.")
    else:
        _log(f"[REF] Gaia catalog loaded: {len(gaia_sky)} sources")

    stats_rows = []
    for row in metrics.to_dict(orient="records"):
        fname = row["file"]
        det_path = _resolve_detect_csv(fname)
        det_xy = np.zeros((0, 2), float)
        if det_path is not None and det_path.exists():
            try:
                df_det = pd.read_csv(det_path)
                if {"x", "y"} <= set(df_det.columns):
                    det_xy = df_det[["x", "y"]].to_numpy(float)
                    det_xy = det_xy[np.isfinite(det_xy).all(axis=1)]
            except Exception:
                det_xy = np.zeros((0, 2), float)

        wcs = _load_wcs_for_frame(fname)
        row["wcs_ok"] = bool(wcs is not None)
        wcs_meta = _combined_wcs_meta(fname)
        row["wcs_resid_med"] = _safe_float(wcs_meta.get("resid_med"), np.nan)
        row["wcs_resid_max"] = _safe_float(wcs_meta.get("resid_max"), np.nan)
        row["wcs_resid_med_px"] = _safe_float(wcs_meta.get("resid_med_px"), np.nan)
        row["wcs_rms_px"] = _safe_float(wcs_meta.get("rms_px"), np.nan)
        row["wcs_center_offset_arcsec"] = _safe_float(wcs_meta.get("center_offset_arcsec"), np.nan)
        qc_pass_raw = wcs_meta.get("wcs_qc_pass", False)
        if isinstance(qc_pass_raw, str):
            row["wcs_qc_pass"] = qc_pass_raw.strip().lower() in {"1", "true", "t", "yes", "y", "pass"}
        else:
            row["wcs_qc_pass"] = bool(qc_pass_raw) if not pd.isna(qc_pass_raw) else False
        row["wcs_qc_reason"] = str(wcs_meta.get("wcs_qc_reason", "") or "")
        row["gaia_source"] = str(wcs_meta.get("gaia_source", "") or "")
        row["wcs_match_n"] = int(wcs_meta.get("match_n", wcs_meta.get("n_match", 0)) or 0)

        if wcs is None or gaia_sky is None:
            row.update(_wcs_meta_match_stats(wcs_meta))
        else:
            row.update(_compute_match_stats(det_xy, wcs, gaia_sky))

        mr_candidates = [
            _safe_float(row.get("match_rate"), np.nan),
            _safe_float(row.get("match_rate_cat"), np.nan),
            _safe_float(row.get("match_rate_eff"), np.nan),
        ]
        mr_candidates = [v for v in mr_candidates if np.isfinite(v)]
        if mr_candidates:
            row["match_rate_eff"] = float(max(mr_candidates))

        stats_rows.append(row)

    metrics = pd.DataFrame(stats_rows)

    try:
        _log(f"[REF][QC] total frames={len(metrics)}")
        if "filter" in metrics.columns:
            counts = metrics["filter"].value_counts(dropna=False)
            parts = [f"{k}:{v}" for k, v in counts.items()]
            _log(f"[REF][QC] filter counts: {', '.join(parts)}")
        if "wcs_ok" in metrics.columns:
            wcs_ok = int(metrics["wcs_ok"].fillna(False).astype(bool).sum())
            _log(f"[REF][QC] wcs_ok={wcs_ok}/{len(metrics)}")
        if "match_rate" in metrics.columns and metrics["match_rate"].notna().any():
            mr = metrics["match_rate"].median()
            _log(f"[REF][QC] match_rate median={mr:.3f}")
        if "match_rate_cat" in metrics.columns and metrics["match_rate_cat"].notna().any():
            mrc = metrics["match_rate_cat"].median()
            _log(f"[REF][QC] match_rate_cat median={mrc:.3f}")
        if "match_rate_eff" in metrics.columns and metrics["match_rate_eff"].notna().any():
            mre = metrics["match_rate_eff"].median()
            _log(f"[REF][QC] match_rate_eff median={mre:.3f}")
        if "sep_med_arcsec" in metrics.columns and metrics["sep_med_arcsec"].notna().any():
            sm = metrics["sep_med_arcsec"].median()
            _log(f"[REF][QC] sep_med_arcsec median={sm:.3f}")
        if "fwhm_px" in metrics.columns and metrics["fwhm_px"].notna().any():
            fmed = metrics["fwhm_px"].median()
            _log(f"[REF][QC] fwhm_px median={fmed:.3f}")
    except Exception:
        pass

    ref_frames_by_date: Dict[str, str] = {}
    ref_filters_by_date: Dict[str, str] = {}
    ref_catalogs_by_date: Dict[str, pd.DataFrame] = {}
    master_df: Optional[pd.DataFrame] = None
    ref_catalog_stats: dict = {}  # Track ref catalog build stats

    match_r = max(0.5, float(wcs_match_radius_arcsec))
    if ref_per_date:
        metrics["date_key"] = metrics.get(
            "date_key",
            metrics["file"].apply(lambda x: _extract_date_key(x, params)),
        )
        for date_key, group in metrics.groupby("date_key", dropna=False):
            group = group.copy()
            ref_fname_date = _select_reference(group, ref_filter)
            ref_filter_date = str(group.loc[group["file"] == ref_fname_date, "filter"].iloc[0])
            ref_frames_by_date[str(date_key)] = ref_fname_date
            ref_filters_by_date[str(date_key)] = ref_filter_date
            _log(f"[REF][QC] date={date_key} ref={ref_fname_date} (filter={ref_filter_date})")

            date_df, date_ref_stats = _build_master_for_group(group, ref_fname_date, match_r)
            date_df = _attach_gaia_photometry(date_df, gaia_df)
            if ref_catalogs_by_date:
                master_df, date_df = _merge_ref_catalogs(
                    master_df, date_df, match_r
                )
            else:
                master_df, date_df = _merge_ref_catalogs(
                    None, date_df, match_r
                )
            ref_catalogs_by_date[str(date_key)] = date_df

    ref_fname = _select_reference(metrics, ref_filter)
    ref_filter_val = str(metrics.loc[metrics["file"] == ref_fname, "filter"].iloc[0])

    _log("=" * 60)
    _log(f"[REF] Selected reference frame: {ref_fname} (filter={ref_filter_val})")

    if not ref_per_date:
        master_df, ref_catalog_stats = _build_master_for_group(metrics, ref_fname, match_r)
        master_df = _attach_gaia_photometry(master_df, gaia_df)

    # Apply hybrid source_id assignment if mode is "hybrid"
    if ref_build_mode == "hybrid":
        master_df, sid_map, id_map = _apply_hybrid_source_ids(master_df, gaia_mag_limit)
        # Also apply to date catalogs if ref_per_date (map to master IDs for consistency)
        if ref_per_date and sid_map:
            for date_key in ref_catalogs_by_date:
                df_date = ref_catalogs_by_date[date_key].copy()
                old_sid = coerce_int64_source_id(df_date["source_id"]) if "source_id" in df_date.columns else None
                if old_sid is not None:
                    mapped_sid = old_sid.map(sid_map).astype("Int64")
                    mapped_id = old_sid.map(id_map).astype("Int64")
                    # Fallback to original IDs if mapping missing
                    df_date["source_id"] = mapped_sid.where(mapped_sid.notna(), old_sid).astype("Int64")
                    fallback_id = coerce_int64_source_id(df_date["ID"]) if "ID" in df_date.columns else old_sid
                    df_date["ID"] = mapped_id.where(mapped_id.notna(), fallback_id).astype("Int64")
                ref_catalogs_by_date[date_key] = df_date

    if "phot_g_mean_mag" in master_df.columns:
        try:
            n_g = int(pd.to_numeric(master_df["phot_g_mean_mag"], errors="coerce").notna().sum())
            _log(f"[REF] Gaia photometry attached: {n_g}/{len(master_df)}")
        except Exception:
            pass

    # ── Neighbor distance + crowding_flag ────────────────────────────────
    if {"x_ref", "y_ref"} <= set(master_df.columns) and len(master_df) > 1:
        try:
            from scipy.spatial import KDTree
            xy = master_df[["x_ref", "y_ref"]].to_numpy(float)
            tree = KDTree(xy)
            dists, _ = tree.query(xy, k=2)  # k=2: self + nearest neighbour
            master_df["neighbor_dist_px"] = dists[:, 1]
            ref_fwhm_row = (
                metrics.loc[metrics["file"] == ref_fname, "fwhm_px"]
                if "fwhm_px" in metrics.columns else None
            )
            ref_fwhm = (
                float(ref_fwhm_row.iloc[0])
                if ref_fwhm_row is not None and ref_fwhm_row.notna().any()
                else float(getattr(params, "fwhm_guess_arcsec", 6.0))
            )
            crowding_mult = float(getattr(params, "crowding_fwhm_mult", 2.5))
            crowding_thresh_px = ref_fwhm * crowding_mult
            master_df["crowding_flag"] = master_df["neighbor_dist_px"] < crowding_thresh_px
            n_crowded = int(master_df["crowding_flag"].sum())
            _log(
                f"[REF] crowding_flag: {n_crowded}/{len(master_df)} crowded "
                f"(fwhm={ref_fwhm:.2f}px, thresh={crowding_thresh_px:.2f}px)"
            )
        except Exception as exc:
            _log(f"[REF] crowding_flag skipped: {exc}")
            master_df["neighbor_dist_px"] = np.nan
            master_df["crowding_flag"] = False

    out_dir = step6_refbuild_dir(result_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    filters = sorted(metrics["filter"].dropna().astype(str).unique().tolist())
    if not filters:
        filters = [ref_filter_val]

    for flt in filters:
        out_path = out_dir / f"ref_catalog_{flt}.tsv"
        master_df.to_csv(out_path, sep="\t", index=False, na_rep="NaN", encoding="utf-8-sig")
        map_path = out_dir / f"sourceid_to_ID_{flt}.csv"
        master_df[["source_id", "ID"]].to_csv(map_path, index=False, encoding="utf-8-sig")
        if ref_per_date:
            for date_key, date_df in ref_catalogs_by_date.items():
                date_path = out_dir / f"ref_catalog_{flt}_{date_key}.tsv"
                date_df.to_csv(date_path, sep="\t", index=False, na_rep="NaN", encoding="utf-8-sig")
                date_map = out_dir / f"sourceid_to_ID_{flt}_{date_key}.csv"
                date_df[["source_id", "ID"]].to_csv(date_map, index=False, encoding="utf-8-sig")

    # Default (no-filter) copies for downstream compatibility
    master_df.to_csv(
        out_dir / "ref_catalog.tsv",
        sep="\t",
        index=False,
        na_rep="NaN",
        encoding="utf-8-sig",
    )
    master_df[["source_id", "ID"]].to_csv(
        out_dir / "sourceid_to_ID.csv", index=False, encoding="utf-8-sig"
    )

    metrics["selected"] = metrics["file"] == ref_fname
    if ref_per_date:
        metrics["date_key"] = metrics.get(
            "date_key",
            metrics["file"].apply(lambda x: _extract_date_key(x, params)),
        )
        metrics["selected_date"] = metrics.apply(
            lambda r: r["file"] == ref_frames_by_date.get(str(r.get("date_key")), ""),
            axis=1
        )
    metrics_path = out_dir / "ref_frame_stats.csv"
    metrics.to_csv(metrics_path, index=False, encoding="utf-8-sig")

    # Compute sky statistics summary from frame metrics
    sky_med_median = np.nan
    sky_sigma_median = np.nan
    if "sky_med" in metrics.columns:
        sky_med_vals = pd.to_numeric(metrics["sky_med"], errors="coerce")
        if sky_med_vals.notna().any():
            sky_med_median = float(sky_med_vals.median())
    if "sky_sigma" in metrics.columns:
        sky_sigma_vals = pd.to_numeric(metrics["sky_sigma"], errors="coerce")
        if sky_sigma_vals.notna().any():
            sky_sigma_median = float(sky_sigma_vals.median())

    meta = {
        "ref_frame": ref_fname,
        "ref_filter": ref_filter_val,
        "ref_per_date": bool(ref_per_date),
        "ref_frames_by_date": ref_frames_by_date,
        "ref_filters_by_date": ref_filters_by_date,
        "sat_drop_pct": float(sat_drop_pct),
        "elong_drop_pct": float(elong_drop_pct),
        "ref_cat_max_sources": int(ref_cat_max_sources),
        "ref_cat_min_sources": int(ref_cat_min_sources),
        "ref_cat_max_elong": float(ref_cat_max_elong),
        "ref_cat_max_abs_round": float(ref_cat_max_abs_round),
        "ref_cat_sharp_min": float(ref_cat_sharp_min),
        "ref_cat_sharp_max": float(ref_cat_sharp_max),
        "ref_cat_min_peak_adu": float(ref_cat_min_peak_adu),
        "wcs_match_radius_arcsec": float(wcs_match_radius_arcsec),
        "wcs_min_match_rate": float(wcs_min_match_rate),
        "wcs_min_match_n": int(wcs_min_match_n),
        "wcs_max_sep_med_arcsec": float(wcs_max_sep_med_arcsec),
        "wcs_max_sep_p90_arcsec": float(wcs_max_sep_p90_arcsec),
        "wcs_max_dup_rate": float(wcs_max_dup_rate),
        "filters": filters,
        # Reference catalog statistics
        "n_ref_total": ref_catalog_stats.get("n_ref_total"),
        "n_ref_after_qualitycuts": ref_catalog_stats.get("n_ref_after_qualitycuts"),
        "n_ref_used": ref_catalog_stats.get("n_ref_used"),
        "used_full_detections": ref_catalog_stats.get("used_full_detections"),
        # Sky statistics (median across frames)
        "sky_med_median": sky_med_median if np.isfinite(sky_med_median) else None,
        "sky_sigma_median": sky_sigma_median if np.isfinite(sky_sigma_median) else None,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    (out_dir / "ref_build_meta.json").write_text(
        json.dumps(meta, indent=2), encoding="utf-8"
    )

    return {
        "ref_frame": ref_fname,
        "ref_filter": ref_filter_val,
        "n_sources": len(master_df),
        "filters": filters,
        "ref_per_date": bool(ref_per_date),
        "ref_frames_by_date": ref_frames_by_date,
    }
