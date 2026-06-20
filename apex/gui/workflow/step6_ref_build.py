"""
Step 6: Master Catalog Build (WCS-based reference catalog)

- Select reference frame using detection-based quality metrics
- Build a fixed master star list from the reference frame detections
- Write master catalogs for downstream steps
"""

from __future__ import annotations

import json
import hashlib
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

from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QGroupBox, QMessageBox,
    QTextEdit, QFormLayout, QProgressBar, QDoubleSpinBox, QSpinBox,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView, QWidget,
    QDialog, QDialogButtonBox, QTabWidget, QCheckBox, QComboBox
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from .step_window_base import StepWindowBase
from .run_control import RunControlBar, format_duration, progress_status_text
from .param_dialog import ParamSpec, run_param_dialog
from .log_panel import WorkflowLogWindow, WorkerStatusPanel, append_timestamped_log, show_raised
from .ui_helpers import (
    create_output_reuse_checkbox,
    create_parameter_button,
    set_table_row_background,
    status_row_background,
)
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
from apex.utils.qc_utils import filter_files_by_qc
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
_REF_SIGNATURE_FILE = "ref_build_signature.json"
_REF_SIGNATURE_VERSION = 2
_REF_SIGNATURE_PARAMS = (
    "wcs_require_qc_pass",
    "global_ref_filter",
    "ref_select_sat_pct",
    "ref_select_elong_pct",
    "ref_cat_max_sources",
    "ref_cat_min_sources",
    "ref_cat_max_elong",
    "ref_cat_max_abs_round",
    "ref_cat_sharp_min",
    "ref_cat_sharp_max",
    "ref_cat_min_peak_adu",
    "ref_wcs_match_radius_arcsec",
    "ref_wcs_min_match_rate",
    "ref_wcs_min_match_n",
    "ref_wcs_max_sep_med_arcsec",
    "ref_wcs_max_sep_p90_arcsec",
    "ref_wcs_max_dup_rate",
    "ref_per_date",
    "ref_build_mode",
    "ref_master_union",
    "ref_union_min_frames",
    "idmatch_gaia_g_limit",
    "gaia_mag_max",
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


# _safe_float imported from utils.common_helpers


class RefBuildWorker(QThread):
    progress = pyqtSignal(int, int, str)
    log = pyqtSignal(str)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str, str)

    def __init__(
        self,
        params,
        data_dir: Path,
        result_dir: Path,
        cache_dir: Path,
        file_list: List[str],
        ref_filter: str,
        sat_drop_pct: float,
        elong_drop_pct: float,
        ref_cat_max_sources: int,
        ref_cat_min_sources: int,
        ref_cat_max_elong: float,
        ref_cat_max_abs_round: float,
        ref_cat_sharp_min: float,
        ref_cat_sharp_max: float,
        ref_cat_min_peak_adu: float,
        wcs_match_radius_arcsec: float,
        wcs_min_match_rate: float,
        wcs_min_match_n: int,
        wcs_max_sep_med_arcsec: float,
        wcs_max_sep_p90_arcsec: float,
        wcs_max_dup_rate: float,
        ref_per_date: bool,
        ref_build_mode: str = "hybrid",
        gaia_mag_limit: float = 18.0,
        ref_master_union: bool = True,
        ref_union_min_frames: int = 1,
    ):
        super().__init__()
        self.params = params
        self.data_dir = Path(data_dir)
        self.result_dir = Path(result_dir)
        self.cache_dir = Path(cache_dir)
        self.file_list = list(file_list)
        self.ref_filter = ref_filter
        self.sat_drop_pct = float(sat_drop_pct)
        self.elong_drop_pct = float(elong_drop_pct)
        self.ref_cat_max_sources = int(ref_cat_max_sources)
        self.ref_cat_min_sources = int(ref_cat_min_sources)
        self.ref_cat_max_elong = float(ref_cat_max_elong)
        self.ref_cat_max_abs_round = float(ref_cat_max_abs_round)
        self.ref_cat_sharp_min = float(ref_cat_sharp_min)
        self.ref_cat_sharp_max = float(ref_cat_sharp_max)
        self.ref_cat_min_peak_adu = float(ref_cat_min_peak_adu)
        self.wcs_match_radius_arcsec = float(wcs_match_radius_arcsec)
        self.wcs_min_match_rate = float(wcs_min_match_rate)
        self.wcs_min_match_n = int(wcs_min_match_n)
        self.wcs_max_sep_med_arcsec = float(wcs_max_sep_med_arcsec)
        self.wcs_max_sep_p90_arcsec = float(wcs_max_sep_p90_arcsec)
        self.wcs_max_dup_rate = float(wcs_max_dup_rate)
        self.ref_per_date = bool(ref_per_date)
        self.ref_build_mode = str(ref_build_mode).lower()
        self.gaia_mag_limit = float(gaia_mag_limit)
        self.ref_master_union = bool(ref_master_union)
        self.ref_union_min_frames = max(1, int(ref_union_min_frames))
        self._stop_requested = False
        # Cache for WCS headers (path -> fits.Header)
        self._wcs_header_cache: Dict[str, fits.Header] = {}
        self._wcs_summary_by_file: Optional[Dict[str, dict]] = None
        self._meta_missing_files: set = set()
        self._meta_incompatible_files: set = set()
        self._detcsv_missing_files: set = set()
        self._detcsv_incompatible_files: set = set()

    def stop(self):
        self._stop_requested = True

    def _log(self, msg: str):
        self.log.emit(msg)

    def _resolve_fits_path(self, fname: str) -> Optional[Path]:
        cropped_dir = step2_cropped_dir(self.result_dir)
        if crop_is_active(self.result_dir):
            cand = cropped_dir / fname
            if cand.exists():
                return cand
        # Prefer original/source path first. Step5 summary signatures are based on
        # the true source frame, and stale copies can remain in older cache folders.
        try:
            orig = Path(self.params.get_file_path(fname))
            if orig.exists():
                return orig
        except Exception:
            pass
        step6_out = step5_wcs_dir(self.result_dir)
        cand = step6_out / fname
        if cand.exists():
            return cand
        try:
            return Path(self.params.get_file_path(fname))
        except Exception:
            return None

    def _current_file_signature(self, fname: str) -> Optional[dict]:
        path = self._resolve_fits_path(fname)
        if path is None or not path.exists():
            return None
        return build_file_signature(path, use_cropped=bool(crop_is_active(self.result_dir)))

    def _detect_meta_compatible(self, fname: str, payload: dict, meta_path: Path) -> bool:
        if not isinstance(payload, dict):
            return False
        schema = cache_schema_value(payload)
        if schema < 2:
            return self._schema1_detect_cache_allowed(meta_path)
        sig_now = self._current_file_signature(fname)
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

    def _schema1_detect_cache_allowed(self, marker_path: Path) -> bool:
        try:
            marker_mtime = int(marker_path.stat().st_mtime_ns)
        except Exception:
            return False
        if crop_is_active(self.result_dir):
            rect_path = crop_rect_path(self.result_dir)
            if rect_path.exists():
                try:
                    rect_mtime = int(rect_path.stat().st_mtime_ns)
                    if marker_mtime < rect_mtime:
                        return False
                except Exception:
                    return False
        return True

    def _wcs_row_compatible(self, fname: str, row: dict) -> bool:
        sig_now = self._current_file_signature(fname)
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

    def _load_meta(self, fname: str) -> Optional[dict]:
        candidates = [
            self.cache_dir / f"detect_{fname}.json",
            step4_dir(self.result_dir) / f"detect_{fname}.json",
        ]
        candidates = [p for p in candidates if p.exists()]
        if not candidates:
            self._meta_missing_files.add(str(fname))
            return None
        candidates.sort(key=lambda p: p.stat().st_mtime_ns, reverse=True)
        for meta_path in candidates:
            try:
                payload = json.loads(meta_path.read_text(encoding="utf-8"))
            except Exception:
                continue
            if self._detect_meta_compatible(fname, payload, meta_path):
                return payload
        self._meta_incompatible_files.add(str(fname))
        return None

    def _resolve_detect_csv(self, fname: str) -> Optional[Path]:
        candidates = [
            self.cache_dir / f"detect_{fname}.csv",
            step4_dir(self.result_dir) / f"detect_{fname}.csv",
        ]
        candidates = [p for p in candidates if p.exists()]
        if not candidates:
            self._detcsv_missing_files.add(str(fname))
            return None
        candidates.sort(key=lambda p: p.stat().st_mtime_ns, reverse=True)
        for cand in candidates:
            meta_path = cand.with_suffix(".json")
            if meta_path.exists():
                try:
                    payload = json.loads(meta_path.read_text(encoding="utf-8"))
                except Exception:
                    continue
                if self._detect_meta_compatible(fname, payload, meta_path):
                    return cand
        self._detcsv_incompatible_files.add(str(fname))
        return None

    def _load_wcs_summary_map(self) -> Dict[str, dict]:
        if self._wcs_summary_by_file is not None:
            return self._wcs_summary_by_file

        mapping: Dict[str, dict] = {}
        candidates = [
            step5_wcs_dir(self.result_dir) / "wcs_solve_summary.csv",
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

        self._wcs_summary_by_file = mapping
        return mapping

    def _load_wcs_summary_row(self, fname: str) -> dict:
        row = self._load_wcs_summary_map().get(str(fname), {})
        if not row:
            return {}
        if not self._wcs_row_compatible(fname, row):
            return {}
        return row

    def _load_wcs_for_frame(self, fname: str) -> Optional[WCS]:
        fits_path = self._resolve_fits_path(fname)
        candidates = []
        if fits_path is not None and fits_path.exists():
            candidates.append(fits_path)
        # fallback: original file if cropped path lacks WCS
        try:
            orig = Path(self.params.get_file_path(fname))
            if orig.exists() and orig not in candidates:
                candidates.append(orig)
        except Exception:
            pass
        for path in candidates:
            try:
                # Use cached header if available
                path_key = str(path)
                if path_key in self._wcs_header_cache:
                    hdr = self._wcs_header_cache[path_key]
                else:
                    hdr = fits.getheader(path)
                    self._wcs_header_cache[path_key] = hdr
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

    def _load_gaia_catalog(self) -> Optional[SkyCoord]:
        try:
            df = self._load_gaia_table()
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

    def _load_gaia_table(self) -> Optional[pd.DataFrame]:
        gaia_path = step5_wcs_dir(self.result_dir) / "gaia_fov.ecsv"
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

    def _attach_gaia_photometry(self, master_df: pd.DataFrame, gaia_df: Optional[pd.DataFrame]) -> pd.DataFrame:
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
        match_r = max(0.5, float(self.wcs_match_radius_arcsec))
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

    def _apply_hybrid_source_ids(self, df: pd.DataFrame, gaia_mag_limit: float = 18.0) -> tuple[pd.DataFrame, dict[int, int], dict[int, int]]:
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
            self._log("[REF] No gaia_source_id column; hybrid mode not applied.")
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
            self._log(f"[REF] Gaia IDs excluded by mag limit (G>{gaia_mag_limit:.2f}): {n_trimmed}")
        self._log(f"[REF] Hybrid IDs assigned: {n_gaia} Gaia, {n_local} local (negative)")

        return out, sid_map, id_map

    def _merge_ref_catalogs(
        self,
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

    def _load_wcs_meta(self, fname: str) -> dict:
        meta_path = Path(self.cache_dir) / "wcs_solve" / f"wcs_{fname}.json"
        if not meta_path.exists():
            return {}
        try:
            return json.loads(meta_path.read_text(encoding="utf-8"))
        except Exception:
            return {}

    def _combined_wcs_meta(self, fname: str) -> dict:
        summary_row = self._load_wcs_summary_row(fname)
        cache_meta = self._load_wcs_meta(fname)
        merged = {}
        if summary_row:
            merged.update(summary_row)
        if cache_meta:
            merged.update(cache_meta)
        return merged

    @staticmethod
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

    def _compute_match_stats(self, det_xy: np.ndarray, wcs: WCS, gaia_sky: SkyCoord) -> dict:
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
        r_match = max(0.5, float(self.wcs_match_radius_arcsec))
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

    def _frame_metrics(self, fname: str) -> Optional[dict]:
        meta = self._load_meta(fname)
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
            "date_key": _extract_date_key(fname, self.params),
            "fwhm_px": fwhm_px,
            "n_sources": n_sources,
            "sat_star_count": sat_count,
            "shape_metric": shape_metric,
            "sky_med": sky_med,
            "sky_sigma": sky_sigma,
        }

    def _select_reference(self, metrics: pd.DataFrame, ref_filter: str) -> str:
        if metrics.empty:
            raise RuntimeError("No detection metrics available")

        df = metrics.copy()
        filt = normalize_filter_key(ref_filter)
        if filt:
            cand = df[df["filter"] == filt].copy()
            if cand.empty:
                self._log(f"[REF][QC] Filter '{filt}' not found; using all filters.")
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
        cand = _drop_top_percent(cand, "sat_star_count", self.sat_drop_pct)
        n_sat = len(cand)
        cand = _drop_top_percent(cand, "shape_metric", self.elong_drop_pct)
        n_shape = len(cand)
        if cand.empty:
            self._log("[REF][QC] All candidates dropped by sat/shape filters; using full set.")
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
        if "match_rate_eff" in cand_wcs.columns and self.wcs_min_match_rate > 0:
            if cand_wcs["match_rate_eff"].notna().any():
                cand_wcs = cand_wcs[cand_wcs["match_rate_eff"] >= self.wcs_min_match_rate]
        if "n_match" in cand_wcs.columns and self.wcs_min_match_n > 0:
            if cand_wcs["n_match"].notna().any():
                cand_wcs = cand_wcs[cand_wcs["n_match"] >= self.wcs_min_match_n]
        if "sep_med_arcsec" in cand_wcs.columns and np.isfinite(self.wcs_max_sep_med_arcsec) and self.wcs_max_sep_med_arcsec > 0:
            if cand_wcs["sep_med_arcsec"].notna().any():
                cand_wcs = cand_wcs[cand_wcs["sep_med_arcsec"] <= self.wcs_max_sep_med_arcsec]
        if "sep_p90_arcsec" in cand_wcs.columns and np.isfinite(self.wcs_max_sep_p90_arcsec) and self.wcs_max_sep_p90_arcsec > 0:
            if cand_wcs["sep_p90_arcsec"].notna().any():
                cand_wcs = cand_wcs[cand_wcs["sep_p90_arcsec"] <= self.wcs_max_sep_p90_arcsec]
        if "dup_rate" in cand_wcs.columns and np.isfinite(self.wcs_max_dup_rate) and self.wcs_max_dup_rate > 0:
            if cand_wcs["dup_rate"].notna().any():
                cand_wcs = cand_wcs[cand_wcs["dup_rate"] <= self.wcs_max_dup_rate]

        n_wcs = len(cand_wcs)
        if not cand_wcs.empty:
            cand = cand_wcs
        self._log(
            "[REF][QC] candidates: start={s} sat_drop={sat:.1f}% -> {n1} "
            "shape_drop={elong:.1f}% -> {n2} wcs_pass -> {n3}".format(
                s=n_start,
                sat=self.sat_drop_pct,
                n1=n_sat,
                elong=self.elong_drop_pct,
                n2=n_shape,
                n3=n_wcs,
            )
        )
        self._log(
            "[REF][QC] wcs thresholds: match_r={r:.2f}\" min_rate={mr:.2f} "
            "min_match={mn} max_sep_med={smed:.2f}\" max_sep_p90={sp90:.2f}\" max_dup={dup:.2f}".format(
                r=self.wcs_match_radius_arcsec,
                mr=self.wcs_min_match_rate,
                mn=self.wcs_min_match_n,
                smed=self.wcs_max_sep_med_arcsec,
                sp90=self.wcs_max_sep_p90_arcsec,
                dup=self.wcs_max_dup_rate,
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
                self._log(
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

    def _build_master_catalog(self, ref_fname: str) -> tuple[pd.DataFrame, dict]:
        """Build master catalog from reference frame detections.

        Returns:
            Tuple of (catalog DataFrame, stats dict with n_ref_total/after_cuts/used)
        """
        det_path = self._resolve_detect_csv(ref_fname)
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

        if np.isfinite(self.ref_cat_max_elong) and self.ref_cat_max_elong > 0 and "elongation" in cand.columns:
            cand = cand[cand["elongation"] <= self.ref_cat_max_elong]

        if np.isfinite(self.ref_cat_max_abs_round) and self.ref_cat_max_abs_round > 0 and "roundness" in cand.columns:
            cand = cand[cand["roundness"].abs() <= self.ref_cat_max_abs_round]

        if "sharpness" in cand.columns:
            if np.isfinite(self.ref_cat_sharp_min):
                cand = cand[cand["sharpness"] >= self.ref_cat_sharp_min]
            if np.isfinite(self.ref_cat_sharp_max) and self.ref_cat_sharp_max > 0:
                cand = cand[cand["sharpness"] <= self.ref_cat_sharp_max]

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
            if np.isfinite(self.ref_cat_min_peak_adu) and self.ref_cat_min_peak_adu > 0:
                cand = cand[cand["brightness"] >= self.ref_cat_min_peak_adu]
            if self.ref_cat_max_sources > 0:
                cand = cand.sort_values("brightness", ascending=False)
                cand = cand.head(self.ref_cat_max_sources)

        base = cand[["x", "y"]].rename(columns={"x": "x_ref", "y": "y_ref"})
        base = base.dropna(subset=["x_ref", "y_ref"]).copy()

        used_full_detections = False
        if base.empty or len(base) < max(10, self.ref_cat_min_sources):
            self._log(
                f"[REF] Ref catalog filter too strict "
                f"({len(base)}/{n_before}); using full detections."
            )
            base = base_all.copy()
            used_full_detections = True
        else:
            self._log(
                f"[REF] Ref catalog filtered: {len(base)}/{n_before} sources kept."
            )

        # Track n_ref_used (final count)
        n_ref_used = len(base)

        wcs = self._load_wcs_for_frame(ref_fname)
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
        self._log(
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
        self, group_metrics: pd.DataFrame, anchor_fname: str, match_radius_arcsec: float
    ) -> tuple[pd.DataFrame, dict]:
        """Build the master catalog for a group of frames.

        With ``ref_master_union`` enabled, detections from *all* frames in the
        group are projected to sky coordinates and merged (deduped by position)
        so that stars only reachable in some frames — faint stars in deep
        exposures, and bright stars that saturate in long exposures but are
        clean in short ones — all enter the master. Otherwise the catalog is
        built from the single anchor (reference) frame only.
        """
        if not self.ref_master_union:
            return self._build_master_catalog(anchor_fname)
        return self._build_union_master(group_metrics, anchor_fname, match_radius_arcsec)

    def _build_union_master(
        self, group_metrics: pd.DataFrame, anchor_fname: str, match_radius_arcsec: float
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
            if self._stop_requested:
                break
            try:
                cat, stats = self._build_master_catalog(fname)
            except RuntimeError as e:
                self._log(f"[REF][UNION] skip {fname}: {e}")
                continue
            last_stats = stats
            n_used += 1
            master, merged = self._merge_ref_catalogs(master, cat, match_radius_arcsec)
            sids = coerce_int64_source_id(merged.get("source_id")).dropna()
            if len(sids):
                det_counts.update(int(s) for s in sids.to_numpy(dtype=np.int64))

        if master is None or master.empty:
            self._log("[REF][UNION] union produced no sources; falling back to anchor frame.")
            return self._build_master_catalog(anchor_fname)

        master = master.copy().reset_index(drop=True)
        sid_int = coerce_int64_source_id(master["source_id"])
        master["n_det_frames"] = [
            int(det_counts.get(int(s), 0)) if pd.notna(s) else 0 for s in sid_int
        ]

        if self.ref_union_min_frames > 1:
            before = len(master)
            master = master[master["n_det_frames"] >= self.ref_union_min_frames].copy()
            self._log(
                f"[REF][UNION] min_frames>={self.ref_union_min_frames}: "
                f"{before}->{len(master)} sources"
            )

        # Express x_ref/y_ref in a single (anchor) pixel system so neighbor /
        # crowding geometry is consistent across stars added from different frames.
        anchor_wcs = self._load_wcs_for_frame(anchor_fname)
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
                self._log(f"[REF][UNION] anchor reprojection failed: {e}")

        # Dense, position-ordered IDs (matches the single-frame build contract).
        sort_cols = [c for c in ("y_ref", "x_ref") if c in master.columns]
        if sort_cols:
            master = master.sort_values(sort_cols).reset_index(drop=True)
        master["source_id"] = np.arange(1, len(master) + 1, dtype=int)
        master["ID"] = master["source_id"]

        stats = dict(last_stats)
        stats["n_master_union"] = int(len(master))
        stats["n_union_frames"] = int(n_used)
        self._log(
            f"[REF][UNION] {n_used} frames -> {len(master)} master sources "
            f"(anchor={anchor_fname}, match_r={match_radius_arcsec:.2f}\")"
        )
        return master, stats

    def run(self):
        """Run the master catalog build.

        Thin wrapper that delegates the compute to the Qt-free
        ``apex.analysis.refbuild.run_refbuild``, re-emitting its callbacks
        through the existing Qt signals. Behavior is identical to the prior
        inline ``_run_impl`` implementation.
        """
        from apex.analysis.refbuild import run_refbuild

        try:
            summary = run_refbuild(
                params=self.params,
                data_dir=self.data_dir,
                result_dir=self.result_dir,
                cache_dir=self.cache_dir,
                file_list=self.file_list,
                ref_filter=self.ref_filter,
                sat_drop_pct=self.sat_drop_pct,
                elong_drop_pct=self.elong_drop_pct,
                ref_cat_max_sources=self.ref_cat_max_sources,
                ref_cat_min_sources=self.ref_cat_min_sources,
                ref_cat_max_elong=self.ref_cat_max_elong,
                ref_cat_max_abs_round=self.ref_cat_max_abs_round,
                ref_cat_sharp_min=self.ref_cat_sharp_min,
                ref_cat_sharp_max=self.ref_cat_sharp_max,
                ref_cat_min_peak_adu=self.ref_cat_min_peak_adu,
                wcs_match_radius_arcsec=self.wcs_match_radius_arcsec,
                wcs_min_match_rate=self.wcs_min_match_rate,
                wcs_min_match_n=self.wcs_min_match_n,
                wcs_max_sep_med_arcsec=self.wcs_max_sep_med_arcsec,
                wcs_max_sep_p90_arcsec=self.wcs_max_sep_p90_arcsec,
                wcs_max_dup_rate=self.wcs_max_dup_rate,
                ref_per_date=self.ref_per_date,
                ref_build_mode=self.ref_build_mode,
                gaia_mag_limit=self.gaia_mag_limit,
                ref_master_union=self.ref_master_union,
                ref_union_min_frames=self.ref_union_min_frames,
                progress_cb=self.progress.emit,
                log_cb=self._log,
                error_cb=self.error.emit,
                should_stop=lambda: self._stop_requested,
            )
        except Exception as e:
            import traceback
            self._log(f"[ERROR] {e}\n{traceback.format_exc()}")
            self.error.emit("WORKER", str(e))
            self.finished.emit({})
            return
        self.finished.emit(summary if summary else {})


_STEP6_SPECS: tuple[ParamSpec, ...] = (
    ParamSpec("Drop top saturation frames", "ref_select_sat_pct", "float", 0.0, 100.0, 1.0, 1, "%", default=20.0),
    ParamSpec("Drop top elongation frames", "ref_select_elong_pct", "float", 0.0, 100.0, 1.0, 1, "%", default=20.0),
    ParamSpec("Per-date reference", "ref_per_date", "bool", default=True),
    ParamSpec("Union master (all frames)", "ref_master_union", "bool", default=True),
    ParamSpec("Union min detections/star", "ref_union_min_frames", "int", 1, 1000, default=1),
    ParamSpec("Ref catalog max sources (0=all)", "ref_cat_max_sources", "int", 0, 50000, default=0),
    ParamSpec("Ref catalog min sources", "ref_cat_min_sources", "int", 0, 50000, default=50),
    ParamSpec("Ref max elongation", "ref_cat_max_elong", "float", 0.0, 10.0, 0.1, 2, default=1.5),
    ParamSpec("Ref max |roundness|", "ref_cat_max_abs_round", "float", 0.0, 5.0, 0.05, 2, default=0.4),
    ParamSpec("Ref sharpness min", "ref_cat_sharp_min", "float", -5.0, 5.0, 0.1, 2, default=0.2),
    ParamSpec("Ref sharpness max", "ref_cat_sharp_max", "float", -5.0, 5.0, 0.1, 2, default=1.0),
    ParamSpec("Ref min peak/flux (0=off)", "ref_cat_min_peak_adu", "float", 0.0, 1e9, 1.0, 1, default=0.0),
    ParamSpec("WCS match radius (arcsec)", "ref_wcs_match_radius_arcsec", "float", 0.1, 30.0, 0.1, 2, default=2.0),
    ParamSpec("WCS min match rate", "ref_wcs_min_match_rate", "float", 0.0, 1.0, 0.05, 2, default=0.2),
    ParamSpec("WCS min match count", "ref_wcs_min_match_n", "int", 0, 100000, default=50),
    ParamSpec("WCS max sep median (arcsec)", "ref_wcs_max_sep_med_arcsec", "float", 0.0, 30.0, 0.1, 2, default=1.5),
    ParamSpec("WCS max sep p90 (arcsec)", "ref_wcs_max_sep_p90_arcsec", "float", 0.0, 60.0, 0.1, 2, default=2.5),
    ParamSpec("WCS max duplicate rate", "ref_wcs_max_dup_rate", "float", 0.0, 1.0, 0.05, 2, default=0.1),
    ParamSpec("Gaia G limit (hybrid ID)", "idmatch_gaia_g_limit", "float", 10.0, 25.0, 0.5, 2, default=18.0),
)


class RefBuildWindow(StepWindowBase):
    """Step 6: Master Catalog Build (WCS-based)."""

    def __init__(self, params, file_manager, project_state, main_window):
        self.file_manager = file_manager
        self.worker = None
        self.log_window = None
        self.results = {}
        self._current_ref_signature: dict | None = None

        super().__init__(
            step_index=5,
            step_name="Master Catalog Build",
            params=params,
            project_state=project_state,
            main_window=main_window,
        )

        self.setup_step_ui()
        self.restore_state()

    def setup_step_ui(self):
        info = QLabel(
            "Build a fixed master catalog using WCS-solved frames.\n"
            "Selection prefers good WCS match stats, then saturation/elongation/FWHM."
        )
        info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 10px; border-radius: 5px; }")
        info.setWordWrap(True)
        self.content_layout.addWidget(info)

        status_group = QGroupBox("WCS/Detection Status")
        status_layout = QVBoxLayout(status_group)
        self.status_label = QLabel("Checking...")
        status_layout.addWidget(self.status_label)
        self.content_layout.addWidget(status_group)

        control_layout = QHBoxLayout()

        btn_params = create_parameter_button("Master Catalog Parameters")
        btn_params.clicked.connect(self.open_parameters_dialog)
        control_layout.addWidget(btn_params)

        self.chk_use_existing_output = create_output_reuse_checkbox(
            not bool(getattr(self.params.P, "force_master_build", False)),
            "When enabled, Step 6 loads the existing reference-build summary and master catalogs "
            "if they are complete. Disable to rebuild the master catalog.",
        )
        control_layout.addWidget(self.chk_use_existing_output)

        self.run_bar = RunControlBar(
            "Run Master Catalog Build", "Show Log",
            run_cb=self.run_ref_build,
            stop_cb=self.stop_ref_build,
            log_cb=self.show_log_window,
        )
        control_layout.addWidget(self.run_bar)
        self.btn_run = self.run_bar.btn_run
        self.btn_stop = self.run_bar.btn_stop

        self.content_layout.addLayout(control_layout)

        progress_layout = QHBoxLayout()
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        progress_layout.addWidget(self.progress_bar)

        self.progress_label = QLabel("Ready")
        self.progress_label.setMinimumWidth(300)
        progress_layout.addWidget(self.progress_label)
        self.content_layout.addLayout(progress_layout)

        self.tabs = QTabWidget()

        summary_tab = QWidget()
        summary_layout = QVBoxLayout(summary_tab)
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(6)
        self.results_table.setHorizontalHeaderLabels([
            "Date", "Filter", "Ref Frame", "Sources", "FWHM (px)", "Sat Count"
        ])
        self.results_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.results_table.horizontalHeader().setStretchLastSection(True)
        self.results_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.results_table.setMinimumHeight(120)
        summary_layout.addWidget(self.results_table)
        self.tabs.addTab(summary_tab, "Summary")

        stats_tab = QWidget()
        stats_layout = QVBoxLayout(stats_tab)
        self.stats_table = QTableWidget()
        self.stats_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.stats_table.setMinimumHeight(220)
        stats_layout.addWidget(self.stats_table)
        self.tabs.addTab(stats_tab, "Stats")

        plot_tab = QWidget()
        plot_layout = QVBoxLayout(plot_tab)
        plot_controls = QHBoxLayout()
        plot_controls.addWidget(QLabel("Date:"))
        self.plot_date_combo = QComboBox()
        self.plot_date_combo.addItem("All")
        self.plot_date_combo.currentIndexChanged.connect(self._on_plot_date_changed)
        plot_controls.addWidget(self.plot_date_combo)
        plot_controls.addStretch()
        plot_layout.addLayout(plot_controls)
        self.plot_canvas = FigureCanvas(Figure(figsize=(8, 4)))
        self.plot_canvas.setMinimumHeight(260)
        plot_layout.addWidget(self.plot_canvas)
        self.tabs.addTab(plot_tab, "Plot")

        self.content_layout.addWidget(self.tabs)

        _worker_group = QGroupBox("Workers")
        _worker_group.setMinimumWidth(300)
        _wg_layout = QVBoxLayout(_worker_group)
        _wg_layout.setContentsMargins(5, 5, 5, 5)
        self.worker_panel = WorkerStatusPanel(_worker_group)
        _wg_layout.addWidget(self.worker_panel)

        self.log_window = WorkflowLogWindow(
            self, "Master Catalog Build Log", width=900, height=500,
            side_widget=_worker_group,
        )
        self.log_text = self.log_window.log_text

        self.check_detection_status()

    def check_detection_status(self):
        cache_dir = Path(self.params.P.cache_dir)
        metas = list(cache_dir.glob("detect_*.json"))
        if not metas:
            step4_out = step4_dir(self.params.P.result_dir)
            metas = list(step4_out.glob("detect_*.json"))
        if not metas:
            self.status_label.setText("No detection cache found. Run Source Detection first.")
            self.status_label.setStyleSheet("color: red;")
            return
        self.status_label.setText(f"Detection cache found: {len(metas)} frames")
        self.status_label.setStyleSheet("color: green;")

    def open_parameters_dialog(self):
        run_param_dialog(
            self, "Master Catalog Build Parameters", _STEP6_SPECS,
            on_save=self.persist_params,
            overrides={"idmatch_gaia_g_limit": getattr(self.params.P, "idmatch_gaia_g_limit",
                        getattr(self.params.P, "gaia_mag_max", 18.0))},
            resize=(460, 560),
            info_text="Adjust reference catalog and WCS quality filter parameters. Changes apply to the next build run.",
        )

    @staticmethod
    def _signature_value(value):
        if isinstance(value, np.generic):
            value = value.item()
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, float):
            return float(value) if np.isfinite(value) else None
        if isinstance(value, (bool, int, str)) or value is None:
            return value
        if isinstance(value, (list, tuple, set)):
            return [RefBuildWindow._signature_value(v) for v in value]
        if isinstance(value, dict):
            return {
                str(k): RefBuildWindow._signature_value(v)
                for k, v in sorted(value.items(), key=lambda item: str(item[0]))
            }
        try:
            if pd.isna(value):
                return None
        except Exception:
            pass
        return str(value)

    @staticmethod
    def _file_signature(path: Path | None) -> dict | None:
        if path is None:
            return None
        try:
            p = Path(path)
            if not p.exists():
                return None
            st = p.stat()
            try:
                path_text = str(p.resolve())
            except Exception:
                path_text = str(p)
            return {
                "path": path_text,
                "size": int(st.st_size),
                "mtime_ns": int(st.st_mtime_ns),
            }
        except Exception:
            return None

    @staticmethod
    def _first_existing(candidates: list[Path]) -> Path | None:
        for path in candidates:
            try:
                if path.exists() and path.stat().st_size > 0:
                    return path
            except Exception:
                continue
        return None

    def _build_ref_output_signature(self, files: list[str]) -> dict:
        result_dir = Path(self.params.P.result_dir)
        cache_dir = Path(self.params.P.cache_dir)
        s4_dir = step4_dir(result_dir)
        s5_dir = step5_wcs_dir(result_dir)

        frame_inputs = []
        for fname in files:
            detect_json = self._first_existing([
                cache_dir / f"detect_{fname}.json",
                s4_dir / f"detect_{fname}.json",
            ])
            detect_csv = self._first_existing([
                cache_dir / f"detect_{fname}.csv",
                s4_dir / f"detect_{fname}.csv",
            ])
            frame_inputs.append({
                "file": str(fname),
                "detect_json": self._file_signature(detect_json),
                "detect_csv": self._file_signature(detect_csv),
            })

        payload = {
            "signature_version": _REF_SIGNATURE_VERSION,
            "step": "step6_ref_build",
            "frames": [str(f) for f in files],
            "params": {
                k: self._signature_value(getattr(self.params.P, k, None))
                for k in _REF_SIGNATURE_PARAMS
            },
            "inputs": {
                "frame_quality": self._file_signature(s4_dir / "frame_quality.csv"),
                "wcs_summary": self._file_signature(s5_dir / "wcs_solve_summary.csv"),
                "gaia_fov": self._file_signature(s5_dir / "gaia_fov.ecsv"),
                "frames": frame_inputs,
            },
        }
        encoded = json.dumps(payload, sort_keys=True, ensure_ascii=False, allow_nan=False)
        payload["signature_hash"] = hashlib.sha1(encoded.encode("utf-8")).hexdigest()
        return payload

    def _stored_ref_signature(self) -> dict | None:
        sig_path = step6_refbuild_dir(self.params.P.result_dir) / _REF_SIGNATURE_FILE
        if not sig_path.exists():
            return None
        try:
            data = json.loads(sig_path.read_text(encoding="utf-8"))
            return data if isinstance(data, dict) else None
        except Exception:
            return None

    def _write_ref_signature(self, signature: dict) -> None:
        out_dir = step6_refbuild_dir(self.params.P.result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / _REF_SIGNATURE_FILE).write_text(
            json.dumps(signature, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False),
            encoding="utf-8",
        )

    def _ref_signature_matches(self, signature: dict) -> tuple[bool, str]:
        stored = self._stored_ref_signature()
        if not stored:
            return False, "missing signature"
        if stored.get("signature_version") != _REF_SIGNATURE_VERSION:
            return False, "signature version mismatch"
        if stored.get("signature_hash") != signature.get("signature_hash"):
            return False, "signature hash mismatch"
        return True, "ok"

    def _current_ref_cache_status(self) -> tuple[bool, str, dict | None]:
        try:
            files = list(self.file_manager.get_file_list())
        except Exception as exc:
            return False, f"file list unavailable: {exc}", None
        if not files:
            return False, "no current frames", None
        use_qc = bool(getattr(self.params.P, "wcs_require_qc_pass", True))
        files, _ = filter_files_by_qc(
            Path(self.params.P.result_dir),
            files,
            require_qc=use_qc,
        )
        if not files:
            return False, "no current frames after QC", None
        signature = self._build_ref_output_signature(list(files))
        matches, reason = self._ref_signature_matches(signature)
        if not matches:
            return False, reason, None
        summary = self._load_existing_summary()
        if not summary or int(summary.get("n_sources", 0) or 0) <= 0:
            return False, "reference output missing or empty", None
        return True, "ok", summary

    def _sync_cache_params(self) -> None:
        use_existing = bool(
            getattr(self, "chk_use_existing_output", None)
            and self.chk_use_existing_output.isChecked()
        )
        self.params.P.force_master_build = not use_existing
        if hasattr(self, "persist_params"):
            self.persist_params()

    def run_ref_build(self):
        if self.worker and self.worker.isRunning():
            return
        self._sync_cache_params()
        files = self.file_manager.get_file_list() if self.file_manager else []
        if not files:
            QMessageBox.warning(self, "Warning", "No frames found")
            return

        use_qc = bool(getattr(self.params.P, "wcs_require_qc_pass", True))
        files, qc_info = filter_files_by_qc(Path(self.params.P.result_dir), files, require_qc=use_qc)
        if use_qc:
            if qc_info.get("applied"):
                self.log(f"[REF][QC] Frame QC filter: {qc_info['kept']}/{qc_info['total']} kept.")
            elif qc_info.get("path") is None:
                self.log("[REF][QC] frame_quality.csv not found; using all frames.")
            else:
                self.log(f"[REF][QC] frame_quality.csv ignored ({qc_info['reason']}); using all frames.")
        if not files:
            QMessageBox.warning(self, "Warning", "No frames after QC filter.")
            return

        ref_filter = normalize_filter_key(getattr(self.params.P, "global_ref_filter", ""))
        self.params.P.ref_select_sat_pct = float(getattr(self.params.P, "ref_select_sat_pct", 20.0))
        self.params.P.ref_select_elong_pct = float(getattr(self.params.P, "ref_select_elong_pct", 20.0))
        self.params.P.ref_cat_max_sources = int(getattr(self.params.P, "ref_cat_max_sources", 0))
        self.params.P.ref_cat_min_sources = int(getattr(self.params.P, "ref_cat_min_sources", 50))
        self.params.P.ref_cat_max_elong = float(getattr(self.params.P, "ref_cat_max_elong", 1.5))
        self.params.P.ref_cat_max_abs_round = float(getattr(self.params.P, "ref_cat_max_abs_round", 0.4))
        self.params.P.ref_cat_sharp_min = float(getattr(self.params.P, "ref_cat_sharp_min", 0.2))
        self.params.P.ref_cat_sharp_max = float(getattr(self.params.P, "ref_cat_sharp_max", 1.0))
        self.params.P.ref_cat_min_peak_adu = float(getattr(self.params.P, "ref_cat_min_peak_adu", 0.0))
        self.params.P.ref_wcs_match_radius_arcsec = float(getattr(self.params.P, "ref_wcs_match_radius_arcsec", 2.0))
        self.params.P.ref_wcs_min_match_rate = float(getattr(self.params.P, "ref_wcs_min_match_rate", 0.2))
        self.params.P.ref_wcs_min_match_n = int(getattr(self.params.P, "ref_wcs_min_match_n", 50))
        self.params.P.ref_wcs_max_sep_med_arcsec = float(getattr(self.params.P, "ref_wcs_max_sep_med_arcsec", 1.5))
        self.params.P.ref_wcs_max_sep_p90_arcsec = float(getattr(self.params.P, "ref_wcs_max_sep_p90_arcsec", 2.5))
        self.params.P.ref_wcs_max_dup_rate = float(getattr(self.params.P, "ref_wcs_max_dup_rate", 0.1))
        self.params.P.ref_per_date = bool(getattr(self.params.P, "ref_per_date", True))
        self.params.P.ref_master_union = bool(getattr(self.params.P, "ref_master_union", True))
        self.params.P.ref_union_min_frames = int(getattr(self.params.P, "ref_union_min_frames", 1))

        signature = self._build_ref_output_signature(files)
        self._current_ref_signature = signature
        if not bool(getattr(self.params.P, "force_master_build", False)):
            sig_ok, sig_reason = self._ref_signature_matches(signature)
            if sig_ok:
                cached_summary = self._load_existing_summary()
                if cached_summary:
                    self.log("[REF][CACHE] Existing Step 6 reference build matches current inputs; using cached output.")
                    self.results = cached_summary
                    self._update_results_table(cached_summary)
                    self._update_stats_table()
                    self._update_plot_tab(cached_summary)
                    self.save_state(cached_summary)
                    self.progress_label.setText("Cached Step 6 output loaded")
                    self.update_navigation_buttons()
                    self._current_ref_signature = None
                    return
            self.log(f"[REF][CACHE] Existing Step 6 output not reusable ({sig_reason}); rebuilding.")

        self.worker = RefBuildWorker(
            params=self.params,
            data_dir=self.params.P.data_dir,
            result_dir=self.params.P.result_dir,
            cache_dir=self.params.P.cache_dir,
            file_list=files,
            ref_filter=ref_filter,
            sat_drop_pct=self.params.P.ref_select_sat_pct,
            elong_drop_pct=self.params.P.ref_select_elong_pct,
            ref_cat_max_sources=self.params.P.ref_cat_max_sources,
            ref_cat_min_sources=self.params.P.ref_cat_min_sources,
            ref_cat_max_elong=self.params.P.ref_cat_max_elong,
            ref_cat_max_abs_round=self.params.P.ref_cat_max_abs_round,
            ref_cat_sharp_min=self.params.P.ref_cat_sharp_min,
            ref_cat_sharp_max=self.params.P.ref_cat_sharp_max,
            ref_cat_min_peak_adu=self.params.P.ref_cat_min_peak_adu,
            wcs_match_radius_arcsec=self.params.P.ref_wcs_match_radius_arcsec,
            wcs_min_match_rate=self.params.P.ref_wcs_min_match_rate,
            wcs_min_match_n=self.params.P.ref_wcs_min_match_n,
            wcs_max_sep_med_arcsec=self.params.P.ref_wcs_max_sep_med_arcsec,
            wcs_max_sep_p90_arcsec=self.params.P.ref_wcs_max_sep_p90_arcsec,
            wcs_max_dup_rate=self.params.P.ref_wcs_max_dup_rate,
            ref_per_date=self.params.P.ref_per_date,
            ref_master_union=self.params.P.ref_master_union,
            ref_union_min_frames=self.params.P.ref_union_min_frames,
            ref_build_mode=getattr(self.params.P, "ref_build_mode", "hybrid"),
            gaia_mag_limit=float(
                getattr(
                    self.params.P,
                    "idmatch_gaia_g_limit",
                    getattr(self.params.P, "gaia_mag_max", 18.0),
                )
            ),
        )
        self.worker.progress.connect(self.on_progress)
        self.worker.log.connect(self.log)
        self.worker.finished.connect(self.on_finished)
        self.worker.error.connect(self.on_error)

        self.run_bar.set_running(True)
        self.progress_bar.setValue(0)
        self.progress_bar.setMaximum(len(files))
        self._ref_start_time = time.monotonic()
        self.progress_label.setText(
            progress_status_text(0, len(files), self._ref_start_time, message="Starting...")
        )
        self.worker.start()
        self.show_log_window()

    def stop_ref_build(self):
        if self.worker and self.worker.isRunning():
            self.worker.stop()

    def on_progress(self, current, total, filename):
        self.progress_bar.setValue(current)
        if hasattr(self, "worker_panel"):
            pct = int(100 * current / max(1, total))
            self.worker_panel.update_worker(0, filename, f"{current}/{total}", pct)
        self.progress_label.setText(
            progress_status_text(current, total, getattr(self, "_ref_start_time", None), message=filename)
        )

    def on_finished(self, summary: dict):
        self.run_bar.set_running(False)
        elapsed_txt = ""
        if hasattr(self, "_ref_start_time"):
            elapsed_txt = f" | elapsed {format_duration(time.monotonic() - self._ref_start_time)}"
        self.progress_label.setText(f"Done{elapsed_txt}")
        if summary:
            self.results = summary
            self._update_results_table(summary)
            self._update_stats_table()
            self._update_plot_tab(summary)
            self.save_state(summary)
            if self._current_ref_signature:
                try:
                    self._write_ref_signature(self._current_ref_signature)
                    self.log("[REF][CACHE] Output signature saved for future reuse.")
                except Exception as exc:
                    self.log(f"[REF][CACHE] Signature write failed: {exc}")
        self._current_ref_signature = None
        self.update_navigation_buttons()

    def on_error(self, filename, error):
        self.log(f"ERROR {filename}: {error}")

    def _update_results_table(self, summary: dict):
        self.results_table.setRowCount(0)
        ref_frame = summary.get("ref_frame", "")
        ref_filter = summary.get("ref_filter", "")
        n_sources = summary.get("n_sources", 0)
        ref_frames_by_date = summary.get("ref_frames_by_date", {}) or {}

        metrics_path = step6_refbuild_dir(self.params.P.result_dir) / "ref_frame_stats.csv"
        metrics_df = None
        if metrics_path.exists():
            try:
                metrics_df = pd.read_csv(metrics_path)
            except Exception:
                metrics_df = None

        def _row_stats(fname: str):
            fwhm = "-"
            sat = "-"
            if metrics_df is not None and fname:
                row = metrics_df[metrics_df["file"] == fname]
                if not row.empty:
                    fwhm = f"{float(row.iloc[0].get('fwhm_px', np.nan)):.2f}" if np.isfinite(row.iloc[0].get('fwhm_px', np.nan)) else "-"
                    sat = str(int(row.iloc[0].get('sat_star_count', 0) or 0))
            return fwhm, sat

        def _has_sources(value) -> bool:
            try:
                return int(float(value)) > 0
            except (TypeError, ValueError):
                return False

        if ref_frames_by_date:
            for date_key, fname in sorted(ref_frames_by_date.items()):
                fwhm, sat = _row_stats(fname)
                row = self.results_table.rowCount()
                self.results_table.insertRow(row)
                self.results_table.setItem(row, 0, QTableWidgetItem(str(date_key)))
                self.results_table.setItem(row, 1, QTableWidgetItem(str(ref_filter)))
                self.results_table.setItem(row, 2, QTableWidgetItem(str(fname)))
                self.results_table.setItem(row, 3, QTableWidgetItem(str(n_sources)))
                self.results_table.setItem(row, 4, QTableWidgetItem(str(fwhm)))
                self.results_table.setItem(row, 5, QTableWidgetItem(str(sat)))
                set_table_row_background(
                    self.results_table,
                    row,
                    status_row_background(bool(fname) and _has_sources(n_sources)),
                )
        else:
            fwhm, sat = _row_stats(ref_frame)
            row = self.results_table.rowCount()
            self.results_table.insertRow(row)
            self.results_table.setItem(row, 0, QTableWidgetItem(_extract_date_key(ref_frame, self.params)))
            self.results_table.setItem(row, 1, QTableWidgetItem(str(ref_filter)))
            self.results_table.setItem(row, 2, QTableWidgetItem(str(ref_frame)))
            self.results_table.setItem(row, 3, QTableWidgetItem(str(n_sources)))
            self.results_table.setItem(row, 4, QTableWidgetItem(str(fwhm)))
            self.results_table.setItem(row, 5, QTableWidgetItem(str(sat)))
            set_table_row_background(
                self.results_table,
                row,
                status_row_background(bool(ref_frame) and _has_sources(n_sources)),
            )

    def _update_stats_table(self):
        stats_path = step6_refbuild_dir(self.params.P.result_dir) / "ref_frame_stats.csv"
        if not stats_path.exists():
            self.stats_table.setRowCount(0)
            self.stats_table.setColumnCount(0)
            return
        try:
            df = pd.read_csv(stats_path)
        except Exception:
            self.stats_table.setRowCount(0)
            self.stats_table.setColumnCount(0)
            return

        preferred_cols = [
            "date_key",
            "file",
            "filter",
            "wcs_ok",
            "match_rate_eff",
            "match_rate",
            "match_rate_cat",
            "n_match",
            "sep_med_arcsec",
            "sep_p90_arcsec",
            "dup_rate",
            "wcs_resid_med",
            "wcs_resid_med_px",
            "wcs_rms_px",
            "wcs_center_offset_arcsec",
            "wcs_match_n",
            "wcs_qc_pass",
            "wcs_qc_reason",
            "gaia_source",
            "fwhm_px",
            "sat_star_count",
            "n_sources",
            "selected",
            "selected_date",
        ]
        cols = [c for c in preferred_cols if c in df.columns]
        if not cols:
            cols = list(df.columns)
        df = df[cols].copy()

        self.stats_table.setColumnCount(len(cols))
        self.stats_table.setRowCount(len(df))
        self.stats_table.setHorizontalHeaderLabels(cols)
        self.stats_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.stats_table.horizontalHeader().setStretchLastSection(True)

        for r, (_, row) in enumerate(df.iterrows()):
            for c, col in enumerate(cols):
                val = row.get(col, "")
                if isinstance(val, float):
                    if "rate" in col or "sep_" in col:
                        text = f"{val:.3f}" if np.isfinite(val) else ""
                    else:
                        text = f"{val:.2f}" if np.isfinite(val) else ""
                else:
                    text = str(val)
                self.stats_table.setItem(r, c, QTableWidgetItem(text))

    def _update_plot_tab(self, summary: Optional[dict] = None) -> None:
        if not hasattr(self, "plot_canvas") or self.plot_canvas is None:
            return
        fig = self.plot_canvas.figure
        fig.clear()
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)

        stats_path = step6_refbuild_dir(self.params.P.result_dir) / "ref_frame_stats.csv"
        if not stats_path.exists():
            ax1.text(0.5, 0.5, "No ref stats available", ha="center", va="center")
            ax2.axis("off")
            self.plot_canvas.draw_idle()
            return
        try:
            df = pd.read_csv(stats_path)
        except Exception:
            ax1.text(0.5, 0.5, "Failed to read ref stats", ha="center", va="center")
            ax2.axis("off")
            self.plot_canvas.draw_idle()
            return
        if df.empty:
            ax1.text(0.5, 0.5, "No ref stats available", ha="center", va="center")
            ax2.axis("off")
            self.plot_canvas.draw_idle()
            return

        if "date_key" in df.columns and hasattr(self, "plot_date_combo"):
            dates = sorted(df["date_key"].fillna("unknown_date").astype(str).unique().tolist())
            items = ["All"] + dates
            prev = self.plot_date_combo.currentText() if self.plot_date_combo.count() else "All"
            self.plot_date_combo.blockSignals(True)
            self.plot_date_combo.clear()
            self.plot_date_combo.addItems(items)
            if prev in items:
                self.plot_date_combo.setCurrentText(prev)
            else:
                self.plot_date_combo.setCurrentIndex(0)
            self.plot_date_combo.blockSignals(False)

        selected_date = None
        if hasattr(self, "plot_date_combo"):
            text = self.plot_date_combo.currentText()
            if text and text != "All":
                selected_date = text

        if "match_rate" in df.columns:
            df["match_rate"] = pd.to_numeric(df["match_rate"], errors="coerce")
        else:
            df["match_rate"] = np.nan
        if "match_rate_eff" in df.columns:
            df["match_rate_eff"] = pd.to_numeric(df["match_rate_eff"], errors="coerce")
        else:
            df["match_rate_eff"] = np.nan
        if "n_match" in df.columns:
            df["n_match"] = pd.to_numeric(df["n_match"], errors="coerce")
        else:
            df["n_match"] = np.nan
        if "sep_med_arcsec" in df.columns:
            df["sep_med_arcsec"] = pd.to_numeric(df["sep_med_arcsec"], errors="coerce")
        else:
            df["sep_med_arcsec"] = np.nan
        if "wcs_resid_med_px" in df.columns:
            df["wcs_resid_med_px"] = pd.to_numeric(df["wcs_resid_med_px"], errors="coerce")
        else:
            df["wcs_resid_med_px"] = np.nan
        if "wcs_rms_px" in df.columns:
            df["wcs_rms_px"] = pd.to_numeric(df["wcs_rms_px"], errors="coerce")
        else:
            df["wcs_rms_px"] = np.nan
        if "fwhm_px" in df.columns:
            df["fwhm_px"] = pd.to_numeric(df["fwhm_px"], errors="coerce")
        else:
            df["fwhm_px"] = np.nan
        if "n_sources" in df.columns:
            df["n_sources"] = pd.to_numeric(df["n_sources"], errors="coerce")
        else:
            df["n_sources"] = np.nan
        if "sat_star_count" in df.columns:
            df["sat_star_count"] = pd.to_numeric(df["sat_star_count"], errors="coerce")
        else:
            df["sat_star_count"] = np.nan

        if "filter" in df.columns:
            filters = sorted(df["filter"].fillna("").astype(str).unique().tolist())
        else:
            filters = [""]
        color_cycle = ["#1E88E5", "#43A047", "#F4511E", "#8E24AA", "#00897B", "#6D4C41"]
        color_map = {f: color_cycle[i % len(color_cycle)] for i, f in enumerate(filters)}

        df_plot = df
        if selected_date and "date_key" in df.columns:
            df_plot = df[df["date_key"].astype(str) == str(selected_date)]
            if df_plot.empty:
                df_plot = df

        def _finite_any(col: str) -> bool:
            return col in df_plot.columns and pd.to_numeric(df_plot[col], errors="coerce").notna().any()

        x_col = "sep_med_arcsec"
        y_col = "match_rate"
        title1 = "Match Rate vs Sep (arcsec)"
        xlabel1 = "Sep med (arcsec)"
        ylabel1 = "Match rate"
        if not _finite_any(x_col):
            if _finite_any("wcs_resid_med_px"):
                x_col = "wcs_resid_med_px"
                xlabel1 = "WCS resid med (px)"
                title1 = "WCS QC Residual vs Match"
            elif _finite_any("wcs_rms_px"):
                x_col = "wcs_rms_px"
                xlabel1 = "WCS RMS (px)"
                title1 = "WCS QC RMS vs Match"
            else:
                x_col = "fwhm_px"
                xlabel1 = "FWHM (px)"
                title1 = "Detection Stats (Gaia match unavailable)"
        if not _finite_any(y_col):
            if _finite_any("match_rate_eff"):
                y_col = "match_rate_eff"
                ylabel1 = "Effective match rate"
            elif _finite_any("n_match"):
                y_col = "n_match"
                ylabel1 = "Matched stars"
            else:
                y_col = "n_sources"
                ylabel1 = "Detected sources"

        for flt in filters:
            sub = df_plot[df_plot["filter"] == flt] if "filter" in df_plot.columns else df_plot
            ax1.scatter(
                sub[x_col],
                sub[y_col],
                s=28,
                alpha=0.75,
                color=color_map.get(flt, "#90A4AE"),
                label=flt or "unknown",
                edgecolors="none",
            )
            ax2.scatter(
                sub["fwhm_px"],
                sub["sat_star_count"],
                s=28,
                alpha=0.75,
                color=color_map.get(flt, "#90A4AE"),
                edgecolors="none",
            )

        selected_frames: Dict[str, str] = {}
        if summary and isinstance(summary, dict):
            ref_by_date = summary.get("ref_frames_by_date", {}) or {}
            if ref_by_date:
                selected_frames.update({str(k): str(v) for k, v in ref_by_date.items()})
            elif summary.get("ref_frame"):
                selected_frames["ref"] = str(summary.get("ref_frame"))
        if not selected_frames:
            meta_path = step6_refbuild_dir(self.params.P.result_dir) / "ref_build_meta.json"
            if meta_path.exists():
                try:
                    meta = json.loads(meta_path.read_text(encoding="utf-8"))
                    ref_by_date = meta.get("ref_frames_by_date", {}) or {}
                    if ref_by_date:
                        selected_frames.update({str(k): str(v) for k, v in ref_by_date.items()})
                    elif meta.get("ref_frame"):
                        selected_frames["ref"] = str(meta.get("ref_frame"))
                except Exception:
                    pass

        if selected_date and selected_date in selected_frames:
            selected_frames = {str(selected_date): selected_frames[str(selected_date)]}

        for _, fname in selected_frames.items():
            if "file" in df_plot.columns:
                row = df_plot[df_plot["file"] == fname]
            else:
                row = df_plot.iloc[0:0]
            if row.empty:
                continue
            r = row.iloc[0]
            ax1.scatter(
                r[x_col], r[y_col],
                s=140, marker="*", color="#FF5252", edgecolors="#212121", linewidths=0.8, zorder=5
            )
            ax2.scatter(
                r["fwhm_px"], r["sat_star_count"],
                s=140, marker="*", color="#FF5252", edgecolors="#212121", linewidths=0.8, zorder=5
            )

        ax1.set_title(title1)
        ax1.set_xlabel(xlabel1)
        ax1.set_ylabel(ylabel1)
        ax1.grid(True, alpha=0.2)
        ax1.legend(fontsize=7, loc="best")

        ax2.set_title("FWHM vs Saturation")
        ax2.set_xlabel("FWHM (px)")
        ax2.set_ylabel("Sat star count")
        ax2.grid(True, alpha=0.2)

        fig.tight_layout()
        self.plot_canvas.draw_idle()

    def _on_plot_date_changed(self, index: int) -> None:
        if index < 0:
            return
        summary = self.results if isinstance(self.results, dict) else None
        self._update_plot_tab(summary)

    def _load_existing_summary(self) -> Optional[dict]:
        out_dir = step6_refbuild_dir(self.params.P.result_dir)
        meta_path = out_dir / "ref_build_meta.json"
        if not meta_path.exists():
            return None
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except Exception as exc:
            self.log(f"[REF][CACHE] Existing ref_build_meta.json ignored: {exc}")
            return None

        ref_filter = str(meta.get("ref_filter", "") or "")
        catalog_candidates = [out_dir / "master_catalog.tsv"]
        if ref_filter:
            catalog_candidates.append(out_dir / f"ref_catalog_{ref_filter}.tsv")
        catalog_candidates.append(out_dir / "ref_catalog.tsv")
        catalog_path = next((p for p in catalog_candidates if p.exists()), None)
        if catalog_path is None:
            self.log("[REF][CACHE] ref_build_meta.json exists but no reference catalog was found; rebuilding.")
            return None

        n_sources = int(meta.get("n_ref_used") or 0)
        if n_sources <= 0:
            try:
                n_sources = len(pd.read_csv(catalog_path, sep="\t"))
            except Exception:
                n_sources = 0

        return {
            "ref_frame": meta.get("ref_frame", ""),
            "ref_filter": ref_filter,
            "n_sources": n_sources,
            "filters": meta.get("filters", []),
            "ref_per_date": bool(meta.get("ref_per_date", False)),
            "ref_frames_by_date": meta.get("ref_frames_by_date", {}) or {},
        }

    def validate_step(self) -> bool:
        valid, _, _ = self._current_ref_cache_status()
        return valid

    def save_state(self, summary: Optional[dict] = None):
        if summary is None:
            summary = self.results if isinstance(self.results, dict) else {}
        summary = summary or {}

        ref_frame = summary.get("ref_frame")
        ref_filter = summary.get("ref_filter")
        n_sources = summary.get("n_sources", 0)

        if not ref_frame or not ref_filter:
            meta_path = step6_refbuild_dir(self.params.P.result_dir) / "ref_build_meta.json"
            if meta_path.exists():
                try:
                    meta = json.loads(meta_path.read_text(encoding="utf-8"))
                    ref_frame = ref_frame or meta.get("ref_frame")
                    ref_filter = ref_filter or meta.get("ref_filter")
                except Exception:
                    pass

        if (not n_sources) and ref_filter:
            ref_path = step6_refbuild_dir(self.params.P.result_dir) / f"ref_catalog_{ref_filter}.tsv"
            if ref_path.exists():
                try:
                    n_sources = len(pd.read_csv(ref_path, sep="\t"))
                except Exception:
                    n_sources = 0

        state_data = {
            "ref_frame": ref_frame,
            "ref_filter": ref_filter,
            "n_sources": n_sources,
        }
        self.project_state.store_step_data("ref_build", state_data)
        if ref_frame:
            self.file_manager.ref_filename = ref_frame

    def restore_state(self):
        state = self.project_state.get_step_data("ref_build")
        valid, reason, summary = self._current_ref_cache_status()
        if valid and summary:
            self.results = summary
            ref_frame = summary.get("ref_frame")
            if ref_frame:
                self.file_manager.ref_filename = ref_frame
            self._update_results_table(summary)
            self._update_stats_table()
            self._update_plot_tab(summary)
            n_sources = summary.get("n_sources", 0)
            self.progress_label.setText(f"Loaded previous reference build ({n_sources} sources)")
            self.update_navigation_buttons()
            try:
                self.log("[REF][CACHE] Loaded previous Step 6 reference build from disk.")
            except Exception:
                pass
            return

        if state or self._load_existing_summary():
            try:
                self.log(f"[REF][CACHE] Previous Step 6 output not restored ({reason}).")
            except Exception:
                pass

    def log(self, msg: str):
        append_timestamped_log(self.log_text, msg)

    def show_log_window(self):
        show_raised(self.log_window)
