"""
Step 7 (UI index 6): Forced Aperture Photometry

Master-catalog-driven forced aperture photometry.
For every frame: project master RA/Dec → pixel coords via WCS,
re-center on detected stars (±max_recenter_shift px), measure flux
at a fixed small aperture, apply aperture correction from bright
isolated detected stars.

Outputs (step_forced_phot/):
  photometry_{fname}.tsv   — per-frame, per-source measurements
  photometry_index.csv     — summary row per frame
  apcorr_summary.csv       — per-frame aperture correction values
"""

from __future__ import annotations

import json
import time
import traceback
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.spatial import cKDTree as KDTree

from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QGroupBox,
    QTextEdit, QProgressBar, QTableWidget, QTableWidgetItem,
    QHeaderView, QAbstractItemView, QWidget, QMessageBox,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal

from .step_window_base import StepWindowBase
from apex.utils.step_paths import (
    step2_cropped_dir,
    step4_dir,
    step6_wcs_dir,
    step7_refbuild_dir,
    step_forced_phot_dir,
    crop_is_active,
)
from apex.utils.photometry_utils import (
    phot_one_star,
    refine_local_centroid,
)
from apex.utils.cache_utils import astap_wcs_candidates, parse_astap_wcs_file

MAG_ERR_COEFF = 2.5 / np.log(10)


# ── Scalar helpers ─────────────────────────────────────────────────────────────

def _to_float(val, default: float) -> float:
    try:
        v = float(val) if val is not None else float(default)
        return v if np.isfinite(v) else float(default)
    except Exception:
        return float(default)


def _to_int(val, default: int) -> int:
    try:
        return int(float(val)) if val is not None else int(default)
    except Exception:
        return int(default)


def _mag_from_flux(flux: float, zp: float = 0.0) -> float:
    if not np.isfinite(flux) or flux <= 0:
        return np.nan
    return zp - 2.5 * np.log10(flux)


def _catalog_series(df: pd.DataFrame, col: str, fallback) -> pd.Series:
    if col in df.columns:
        return df[col].reset_index(drop=True)
    return pd.Series(fallback)


def _normalize_filter_value(value) -> str:
    raw = str(value or "").strip().lower()
    if not raw:
        return ""
    import re
    tokens = [t for t in re.split(r"[^a-z0-9]+", raw) if t]
    for token in reversed(tokens or [raw]):
        if token in {"ha", "halpha", "h-alpha"}:
            return "ha"
        if token in {"u", "g", "r", "i", "z", "b", "v"}:
            return token
    if raw in {"u", "g", "r", "i", "z", "b", "v"}:
        return raw
    return raw


# ── ForcedPhotWorker ───────────────────────────────────────────────────────────

class ForcedPhotWorker(QThread):
    """Per-frame forced aperture photometry worker."""

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
        output_dir: Optional[Path] = None,
    ):
        super().__init__()
        self.params = params
        self.data_dir = Path(data_dir)
        self.result_dir = Path(result_dir)
        self.cache_dir = Path(cache_dir)
        self.file_list = list(file_list)
        self.output_dir = Path(output_dir) if output_dir is not None else step_forced_phot_dir(result_dir)
        self._stop_requested = False
        self._wcs_header_cache: Dict[str, fits.Header] = {}

    def stop(self):
        self._stop_requested = True

    def _log(self, msg: str):
        self.log.emit(msg)

    # ── Path resolution ────────────────────────────────────────────────────────

    def _resolve_fits_path(self, fname: str) -> Optional[Path]:
        if crop_is_active(self.result_dir):
            p = step2_cropped_dir(self.result_dir) / fname
            if p.exists():
                return p
        try:
            orig = Path(self.params.get_file_path(fname))
            if orig.exists():
                return orig
        except Exception:
            pass
        cand = self.data_dir / fname
        return cand if cand.exists() else None

    def _load_wcs_for_frame(self, fname: str) -> Optional[WCS]:
        fits_path = self._resolve_fits_path(fname)
        candidates: List[Path] = []
        if fits_path is not None and fits_path.exists():
            candidates.append(fits_path)
        try:
            orig = Path(self.params.get_file_path(fname))
            if orig not in candidates and orig.exists():
                candidates.append(orig)
        except Exception:
            pass

        for path in candidates:
            try:
                key = str(path)
                if key in self._wcs_header_cache:
                    hdr = self._wcs_header_cache[key]
                else:
                    hdr = fits.getheader(path)
                    self._wcs_header_cache[key] = hdr
                w = WCS(hdr, relax=True)
                if w.has_celestial:
                    return w
                # Try ASTAP sidecar .wcs file
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

    def _load_image(self, fname: str) -> Optional[np.ndarray]:
        path = self._resolve_fits_path(fname)
        if path is None:
            return None
        try:
            with fits.open(path) as hdul:
                for hdu in hdul:
                    if hdu.data is not None and hdu.data.ndim == 2:
                        return hdu.data.astype(float)
        except Exception:
            pass
        return None

    def _load_fwhm(self, fname: str, fallback: float = 6.0) -> float:
        for cand in [
            self.cache_dir / f"detect_{fname}.json",
            step4_dir(self.result_dir) / f"detect_{fname}.json",
        ]:
            if not cand.exists():
                continue
            try:
                meta = json.loads(cand.read_text(encoding="utf-8"))
                for key in ("fwhm_med_rad_px", "fwhm_med_px", "fwhm_px", "fwhm_med"):
                    v = meta.get(key)
                    if v is not None:
                        f = float(v)
                        if np.isfinite(f) and f > 0:
                            return f
            except Exception:
                continue
        return fallback

    def _load_sky_stats(self, fname: str) -> Tuple[float, float]:
        """Load frame-level sky median and sigma (electrons) from detect meta."""
        for cand in [
            self.cache_dir / f"detect_{fname}.json",
            step4_dir(self.result_dir) / f"detect_{fname}.json",
        ]:
            if not cand.exists():
                continue
            try:
                meta = json.loads(cand.read_text(encoding="utf-8"))
                sky_med = float(meta.get("sky_med") or meta.get("bkg_median") or np.nan)
                sky_sig = float(meta.get("sky_sigma") or meta.get("bkg_rms") or np.nan)
                return sky_med, sky_sig
            except Exception:
                continue
        return np.nan, np.nan

    def _load_detect_positions(self, fname: str) -> Optional[pd.DataFrame]:
        """Load detected source positions (det_uid, x, y) for a frame."""
        for cand in [
            self.cache_dir / f"detect_{fname}.csv",
            step4_dir(self.result_dir) / f"detect_{fname}.csv",
        ]:
            if not cand.exists() or cand.stat().st_size == 0:
                continue
            try:
                df = pd.read_csv(cand)
                xc = "x" if "x" in df.columns else ("xcenter" if "xcenter" in df.columns else None)
                yc = "y" if "y" in df.columns else ("ycenter" if "ycenter" in df.columns else None)
                if xc is None or yc is None:
                    continue
                out = pd.DataFrame({
                    "det_uid": (
                        pd.to_numeric(df["det_uid"], errors="coerce").fillna(
                            pd.Series(range(len(df)))
                        ).astype(int)
                        if "det_uid" in df.columns
                        else pd.Series(range(len(df)))
                    ),
                    "x": pd.to_numeric(df[xc], errors="coerce"),
                    "y": pd.to_numeric(df[yc], errors="coerce"),
                })
                return out.dropna(subset=["x", "y"]).reset_index(drop=True)
            except Exception:
                continue
        return None

    def _load_master_catalog(self, filt: str) -> Optional[pd.DataFrame]:
        """Load master catalog for a given filter from step7_refbuild."""
        refbuild_dir = step7_refbuild_dir(self.result_dir)
        candidates = [
            refbuild_dir / f"ref_catalog_{filt}.tsv",
            refbuild_dir / f"ref_catalog_{filt.lower()}.tsv",
            refbuild_dir / f"ref_catalog_{filt.upper()}.tsv",
            refbuild_dir / "ref_catalog.tsv",
        ]
        for path in candidates:
            if not path.exists():
                continue
            try:
                df = pd.read_csv(path, sep="\t")
                # Must have at least ra_deg/dec_deg or x_ref/y_ref
                has_radec = {"ra_deg", "dec_deg"} <= set(df.columns)
                has_xy = {"x_ref", "y_ref"} <= set(df.columns)
                if not has_radec and not has_xy:
                    continue
                if "master_id" not in df.columns:
                    if "ID" in df.columns:
                        df.insert(0, "master_id", df["ID"])
                    elif "source_id" in df.columns:
                        df.insert(0, "master_id", df["source_id"])
                    else:
                        df.insert(0, "master_id", range(1, len(df) + 1))
                if "source_id" not in df.columns:
                    df["source_id"] = df["master_id"]
                if "ID" not in df.columns:
                    df["ID"] = df["master_id"]
                for col in ("crowding_flag", "neighbor_dist_px"):
                    if col not in df.columns:
                        df[col] = False if col == "crowding_flag" else np.nan
                return df.reset_index(drop=True)
            except Exception as exc:
                self._log(f"[FORCED] catalog load failed ({path.name}): {exc}")
        return None

    # ── Per-frame photometry ───────────────────────────────────────────────────

    def _phot_frame(
        self,
        fname: str,
        master_df: pd.DataFrame,
        img: np.ndarray,
        wcs: Optional[WCS],
        fwhm_px: float,
        P,
    ) -> pd.DataFrame:
        """Measure forced aperture photometry for all master positions in one frame."""
        n = len(master_df)
        h, w = img.shape

        # Photometry parameters
        r_ap_scale    = _to_float(getattr(P, "forced_r_ap_scale",     0.8), 0.8)
        ref_ap_scale  = _to_float(getattr(P, "forced_ref_ap_scale",   2.4), 2.4)
        min_r_ap      = _to_float(getattr(P, "min_r_ap_px",           4.0), 4.0)
        ann_scale     = _to_float(getattr(P, "fitsky_annulus_scale",   4.0), 4.0)
        dann_scale    = _to_float(getattr(P, "fitsky_dannulus_scale",  2.0), 2.0)
        ann_gap       = _to_float(getattr(P, "annulus_min_gap_px",     6.0), 6.0)
        max_shift     = _to_float(getattr(P, "max_recenter_shift",     2.0), 2.0)
        gain          = _to_float(getattr(P, "gain_e_per_adu",         1.0), 1.0)
        rn_e          = _to_float(getattr(P, "rdnoise_e",              7.5), 7.5)
        sat_adu       = _to_float(getattr(P, "saturation_adu",     65000.0), 65000.0)
        datamax_adu   = _to_float(getattr(P, "datamax_adu",        55000.0), 55000.0)
        sigma_clip    = _to_float(getattr(P, "phot_sigma_clip",        3.0), 3.0)
        max_iter      = _to_int(  getattr(P, "phot_max_iter",            5), 5)
        sky_mode      = str(getattr(P, "sky_sigma_mode", "local") or "local").strip().lower()
        sky_incl_rn   = bool(getattr(P, "sky_sigma_includes_rn", True))
        min_n_sky     = _to_int(getattr(P, "sky_sigma_min_n_sky", 50), 50)
        cbox_scale    = _to_float(getattr(P, "center_cbox_scale", 1.5), 1.5)
        apcorr_min_snr = _to_float(getattr(P, "apcorr_min_snr", 40.0), 40.0)
        apcorr_min_n   = _to_int(getattr(P, "apcorr_use_min_n", 12), 12)

        r_ap  = max(min_r_ap, r_ap_scale * fwhm_px)
        r_ref = max(r_ap + 2.0, ref_ap_scale * fwhm_px)
        r_in  = max(r_ref + ann_gap, ann_scale * fwhm_px)
        r_out = r_in + max(ann_gap, dann_scale * fwhm_px)

        sky_med_frame, sky_sigma_frame = self._load_sky_stats(fname)
        # Convert sky sigma to electrons if given in ADU
        sky_e = sky_sigma_frame * gain if np.isfinite(sky_sigma_frame) else np.nan

        # --- Project master positions to pixels via WCS ---
        has_radec = {"ra_deg", "dec_deg"} <= set(master_df.columns)
        has_xy    = {"x_ref", "y_ref"} <= set(master_df.columns)

        x_pred = np.full(n, np.nan)
        y_pred = np.full(n, np.nan)

        if has_radec and wcs is not None:
            ra  = pd.to_numeric(master_df["ra_deg"],  errors="coerce").to_numpy(float)
            dec = pd.to_numeric(master_df["dec_deg"], errors="coerce").to_numpy(float)
            valid = np.isfinite(ra) & np.isfinite(dec)
            if valid.any():
                try:
                    xy = wcs.all_world2pix(
                        np.column_stack([ra[valid], dec[valid]]), 0
                    )
                    x_pred[valid] = xy[:, 0]
                    y_pred[valid] = xy[:, 1]
                except Exception as exc:
                    self._log(f"[FORCED][{fname}] WCS project failed: {exc}")

        # Fallback: use reference frame pixel coords
        if has_xy:
            x_ref = pd.to_numeric(master_df["x_ref"], errors="coerce").to_numpy(float)
            y_ref = pd.to_numeric(master_df["y_ref"], errors="coerce").to_numpy(float)
            missing = ~np.isfinite(x_pred)
            x_pred[missing] = x_ref[missing]
            y_pred[missing] = y_ref[missing]

        # --- Match to detections for recentering ---
        det_df = self._load_detect_positions(fname)
        detected_flag  = np.zeros(n, dtype=bool)
        det_uid_arr    = np.full(n, -1, dtype=np.int64)

        if det_df is not None and len(det_df) > 0:
            det_xy = det_df[["x", "y"]].to_numpy(float)
            master_xy_pred = np.column_stack([x_pred, y_pred])
            valid_pred = np.isfinite(x_pred) & np.isfinite(y_pred)
            if valid_pred.any():
                tree = KDTree(det_xy)
                dists, idxs = tree.query(master_xy_pred[valid_pred], k=1)
                match_mask = dists <= max_shift
                valid_indices = np.where(valid_pred)[0]
                for k, (vi, matched, idx) in enumerate(zip(valid_indices, match_mask, idxs)):
                    if matched:
                        detected_flag[vi] = True
                        det_uid_arr[vi] = int(det_df["det_uid"].iloc[idx])

        # --- Per-source photometry ---
        x_fit = x_pred.copy()
        y_fit = y_pred.copy()
        centroid_shift = np.full(n, np.nan)

        # Recenter detected sources
        for i in range(n):
            if not detected_flag[i]:
                continue
            xi, yi = x_pred[i], y_pred[i]
            if not (np.isfinite(xi) and np.isfinite(yi)):
                continue
            try:
                xc, yc = refine_local_centroid(img, xi, yi, fwhm_px, cbox_scale)
                shift = np.hypot(xc - xi, yc - yi)
                if shift <= max_shift:
                    x_fit[i] = xc
                    y_fit[i] = yc
                    centroid_shift[i] = shift
            except Exception:
                pass

        # Run aperture photometry for every master position
        flux_arr     = np.full(n, np.nan)
        flux_err_arr = np.full(n, np.nan)
        snr_arr      = np.full(n, np.nan)
        sky_arr      = np.full(n, np.nan)
        is_sat_arr   = np.zeros(n, dtype=bool)
        is_nl_arr    = np.zeros(n, dtype=bool)

        # Also measure at the reference aperture for apcorr
        flux_ref_arr = np.full(n, np.nan)

        for i in range(n):
            xi, yi = x_fit[i], y_fit[i]
            if not (np.isfinite(xi) and np.isfinite(yi)):
                continue
            if not (0 <= xi < w and 0 <= yi < h):
                continue
            try:
                (flux_e, sigma_e, snr, _, bkg_med, _, _, _, is_sat, is_nl) = phot_one_star(
                    img, xi, yi, r_ap, r_in, r_out,
                    sigma_clip_val=sigma_clip,
                    maxiters=max_iter,
                    gain=gain,
                    rn_param_e=rn_e,
                    sky_frame_e=sky_e,
                    sky_sigma_mode=sky_mode,
                    sky_sigma_includes_rn=sky_incl_rn,
                    min_n_sky_for_local=min_n_sky,
                    sat_adu=sat_adu,
                    datamax_adu=datamax_adu,
                )
                flux_arr[i]     = flux_e
                flux_err_arr[i] = sigma_e
                snr_arr[i]      = snr
                sky_arr[i]      = bkg_med
                is_sat_arr[i]   = is_sat
                is_nl_arr[i]    = is_nl

                # Reference aperture measurement (for apcorr)
                (flux_ref_e, *_) = phot_one_star(
                    img, xi, yi, r_ref, r_in, r_out,
                    sigma_clip_val=sigma_clip,
                    maxiters=max_iter,
                    gain=gain,
                    rn_param_e=rn_e,
                    sky_frame_e=sky_e,
                    sky_sigma_mode=sky_mode,
                    sky_sigma_includes_rn=sky_incl_rn,
                    min_n_sky_for_local=min_n_sky,
                    sat_adu=sat_adu,
                    datamax_adu=datamax_adu,
                )
                flux_ref_arr[i] = flux_ref_e
            except Exception:
                pass

        # --- Aperture correction ---
        crowding_flag = master_df["crowding_flag"].to_numpy(bool) if "crowding_flag" in master_df.columns else np.zeros(n, dtype=bool)
        apcorr_candidates = (
            detected_flag &
            ~crowding_flag &
            np.isfinite(flux_arr) & (flux_arr > 0) &
            np.isfinite(flux_ref_arr) & (flux_ref_arr > 0) &
            np.isfinite(snr_arr) & (snr_arr >= apcorr_min_snr) &
            ~is_sat_arr
        )
        apcorr_values = flux_ref_arr[apcorr_candidates] / flux_arr[apcorr_candidates]
        apcorr_values = apcorr_values[np.isfinite(apcorr_values) & (apcorr_values > 0)]

        if len(apcorr_values) >= apcorr_min_n:
            # Sigma-clip the apcorr values
            from astropy.stats import sigma_clipped_stats
            _, apcorr_med, _ = sigma_clipped_stats(apcorr_values, sigma=3.0)
            apcorr = float(apcorr_med) if np.isfinite(apcorr_med) else 1.0
        else:
            apcorr = 1.0

        # Apply apcorr
        flux_corr = flux_arr * apcorr
        flux_err_corr = flux_err_arr * apcorr

        # Compute instrumental magnitudes
        mag_inst = np.where(
            np.isfinite(flux_corr) & (flux_corr > 0),
            -2.5 * np.log10(np.where(flux_corr > 0, flux_corr, np.nan)),
            np.nan,
        )
        mag_err = np.where(
            np.isfinite(flux_corr) & (flux_corr > 0) & np.isfinite(flux_err_corr) & (flux_err_corr >= 0),
            MAG_ERR_COEFF * flux_err_corr / np.where(flux_corr > 0, flux_corr, np.nan),
            np.nan,
        )

        # Bad phot flag: saturated, nonlinear, or no valid flux
        bad_phot = is_sat_arr | is_nl_arr | ~np.isfinite(flux_corr)

        # source_id is the stable master/Gaia identifier used downstream.
        # det_uid is only the frame-local detection id used for recentering.
        master_id = _catalog_series(master_df, "master_id", np.arange(1, n + 1))
        source_id = _catalog_series(master_df, "source_id", master_id)
        display_id = _catalog_series(master_df, "ID", master_id)

        out = pd.DataFrame({
            "master_id":       master_id,
            "source_id":       source_id,
            "ID":              display_id,
            "det_uid":         det_uid_arr,
            "x":               x_fit,
            "y":               y_fit,
            "x_pred":          x_pred,
            "y_pred":          y_pred,
            "x_fit":           x_fit,
            "y_fit":           y_fit,
            "detected_flag":   detected_flag,
            "forced_flag":     ~detected_flag,
            "centroid_shift_px": centroid_shift,
            "flux":            flux_corr,
            "flux_e":          flux_corr,
            "flux_net_adu":    flux_corr / max(gain, 1e-12),
            "flux_err":        flux_err_corr,
            "mag_inst":        mag_inst,
            "mag_err":         mag_err,
            "snr":             snr_arr,
            "sky":             sky_arr,
            "apcorr":          apcorr,
            "is_saturated":    is_sat_arr,
            "is_nonlinear":    is_nl_arr,
            "bad_phot_flag":   bad_phot,
        })
        for col in ("ra_deg", "dec_deg", "gaia_source_id"):
            if col in master_df.columns and col not in out.columns:
                out[col] = master_df[col].reset_index(drop=True)
        return out, apcorr

    # ── Main run ──────────────────────────────────────────────────────────────

    def run(self):
        try:
            self._run_impl()
        except Exception as exc:
            self._log(f"[FORCED] Unhandled error: {exc}\n{traceback.format_exc()}")
            self.error.emit("ForcedPhot", str(exc))
            self.finished.emit({})

    def _run_impl(self):  # noqa: C901
        P = self.params.P
        out_dir = self.output_dir
        out_dir.mkdir(parents=True, exist_ok=True)

        fwhm_guess = _to_float(getattr(P, "fwhm_pix_guess", 6.0), 6.0)

        def _get_filter(fname: str) -> str:
            fits_path = self._resolve_fits_path(fname)
            if fits_path is not None and fits_path.exists():
                try:
                    hdr = fits.getheader(fits_path)
                    for key in ("FILTER", "FILTER1", "FILTER2", "FILTNAM", "FILTERID"):
                        filt = _normalize_filter_value(hdr.get(key))
                        if filt:
                            return filt
                except Exception:
                    pass
            import re
            m = re.search(r"[-_]([ugrizbvUGRIZBV])[-_.]", str(fname), re.IGNORECASE)
            return m.group(1).lower() if m else "unknown"

        filter_map: Dict[str, List[str]] = {}
        for fname in self.file_list:
            filt = _get_filter(fname)
            filter_map.setdefault(filt, []).append(fname)
        try:
            # Compatibility artifact: LC/CMD readers still use the old Step 8
            # frame map name even though forced phot now replaces IDMatch.
            (out_dir / "step8_filter_frames.json").write_text(
                json.dumps(filter_map, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
        except Exception as exc:
            self._log(f"[FORCED] step8_filter_frames write failed: {exc}")

        # Load master catalogs
        master_catalogs: Dict[str, Optional[pd.DataFrame]] = {}
        for filt in filter_map:
            cat = self._load_master_catalog(filt)
            if cat is None:
                # Try loading a generic catalog for any filter
                cat = self._load_master_catalog("r") or self._load_master_catalog("unknown")
            master_catalogs[filt] = cat
            if cat is not None:
                self._log(f"[FORCED] Loaded master catalog for filter '{filt}': {len(cat)} sources")
            else:
                self._log(f"[FORCED] WARNING: No master catalog found for filter '{filt}'")
        try:
            master_rows = []
            for filt, cat in master_catalogs.items():
                if cat is None or cat.empty:
                    continue
                tmp = cat.copy()
                tmp.insert(0, "filter", filt)
                master_rows.append(tmp)
            if master_rows:
                pd.concat(master_rows, ignore_index=True).drop_duplicates(
                    subset=[c for c in ("filter", "source_id") if c in master_rows[0].columns]
                ).to_csv(out_dir / "step8_master_sources.csv", index=False)
        except Exception as exc:
            self._log(f"[FORCED] step8_master_sources write failed: {exc}")

        total = len(self.file_list)
        index_rows: List[dict] = []
        apcorr_rows: List[dict] = []
        n_done = 0

        for filt, files in filter_map.items():
            master_df = master_catalogs.get(filt)
            if master_df is None or master_df.empty:
                self._log(f"[FORCED] Skipping filter '{filt}' — no master catalog")
                n_done += len(files)
                continue

            for fname in files:
                if self._stop_requested:
                    self._log("[FORCED] Stopped by user.")
                    self.finished.emit({"index_rows": index_rows, "apcorr_rows": apcorr_rows})
                    return

                n_done += 1
                self.progress.emit(n_done, total, fname)

                # Load image
                img = self._load_image(fname)
                if img is None:
                    self._log(f"[FORCED] [{fname}] Image not found, skipping")
                    index_rows.append({"file": fname, "filter": filt, "status": "no_image",
                                       "n_master": 0, "n_detected": 0, "n_forced": 0, "apcorr": np.nan})
                    continue

                # Load WCS
                wcs = self._load_wcs_for_frame(fname)
                if wcs is None:
                    self._log(f"[FORCED] [{fname}] No WCS — using reference pixel coords")

                # Load FWHM
                fwhm_px = self._load_fwhm(fname, fallback=fwhm_guess)

                # Run photometry
                try:
                    phot_df, apcorr_val = self._phot_frame(fname, master_df, img, wcs, fwhm_px, P)
                except Exception as exc:
                    self._log(f"[FORCED] [{fname}] phot_frame failed: {exc}\n{traceback.format_exc()}")
                    index_rows.append({"file": fname, "filter": filt, "status": "error",
                                       "n_master": len(master_df), "n_detected": 0, "n_forced": 0,
                                       "wcs_ok": bool(wcs is not None), "apcorr": np.nan})
                    continue

                # Write per-frame TSV
                out_path = out_dir / f"photometry_{fname}.tsv"
                try:
                    phot_df.insert(0, "file", fname)
                    phot_df.insert(1, "filter", filt)
                    phot_df.insert(2, "FILTER", filt)
                    phot_df.to_csv(out_path, sep="\t", index=False, float_format="%.6f")
                except Exception as exc:
                    self._log(f"[FORCED] [{fname}] write failed: {exc}")

                n_det    = int(phot_df["detected_flag"].sum())
                n_forced = int(phot_df["forced_flag"].sum())
                n_valid  = int(np.isfinite(phot_df["flux"]).sum())

                self._log(
                    f"[FORCED] [{fname}] fwhm={fwhm_px:.1f}px  "
                    f"det={n_det}  forced={n_forced}  "
                    f"valid_phot={n_valid}  apcorr={apcorr_val:.4f}"
                )

                index_rows.append({
                    "file":       fname,
                    "filter":     filt,
                    "status":     "ok",
                    "n_master":   len(master_df),
                    "n_detected": n_det,
                    "n_forced":   n_forced,
                    "n_valid_phot": n_valid,
                    "fwhm_px":    fwhm_px,
                    "wcs_ok":     bool(wcs is not None),
                    "apcorr":     apcorr_val,
                    "path":       str(out_path),
                })
                apcorr_rows.append({"file": fname, "filter": filt, "apcorr": apcorr_val,
                                    "n_apcorr_stars": n_det})

        # Write index
        if index_rows:
            idx_df = pd.DataFrame(index_rows)
            try:
                idx_df.to_csv(out_dir / "photometry_index.csv", index=False, float_format="%.6f")
                idx_df.to_csv(out_dir / "step8_frame_stats.csv", index=False, float_format="%.6f")
            except Exception as exc:
                self._log(f"[FORCED] photometry_index write failed: {exc}")

        if apcorr_rows:
            apc_df = pd.DataFrame(apcorr_rows)
            try:
                apc_df.to_csv(out_dir / "apcorr_summary.csv", index=False, float_format="%.6f")
            except Exception as exc:
                self._log(f"[FORCED] apcorr_summary write failed: {exc}")

        n_ok     = sum(1 for r in index_rows if r.get("status") == "ok")
        n_det_t  = sum(r.get("n_detected", 0) for r in index_rows)
        n_frc_t  = sum(r.get("n_forced", 0) for r in index_rows)
        self._log(
            f"[FORCED] Done. {n_ok}/{total} frames OK. "
            f"Total detected={n_det_t} forced={n_frc_t}"
        )

        self.finished.emit({"index_rows": index_rows, "apcorr_rows": apcorr_rows})


# ── ForcedPhotWindow ───────────────────────────────────────────────────────────

class ForcedPhotWindow(StepWindowBase):
    """Step 7: Forced Aperture Photometry (master-driven)."""

    def __init__(self, params, file_manager, project_state, main_window):
        self.file_manager = file_manager
        self.worker = None
        self.results: dict = {}

        super().__init__(
            step_index=6,
            step_name="Forced Aperture Phot",
            params=params,
            project_state=project_state,
            main_window=main_window,
        )

        self.setup_step_ui()
        self.restore_state()

    def setup_step_ui(self):
        info = QLabel(
            "Project master catalog positions onto each frame via WCS and measure forced "
            "aperture photometry.\nForces measurement even for stars not individually detected "
            "in a frame — recovers crowded-field detections lost to IDMatch."
        )
        info.setStyleSheet("QLabel { background-color: #E8F5E9; padding: 10px; border-radius: 5px; }")
        info.setWordWrap(True)
        self.content_layout.addWidget(info)

        # Status
        status_group = QGroupBox("Prerequisites")
        status_layout = QVBoxLayout(status_group)
        self.status_label = QLabel("Checking...")
        status_layout.addWidget(self.status_label)
        self.content_layout.addWidget(status_group)
        self._check_prerequisites()

        # Controls
        ctrl_layout = QHBoxLayout()
        self.btn_run = QPushButton("Run Forced Photometry")
        self.btn_run.setStyleSheet(
            "QPushButton { background-color: #4CAF50; color: white; font-weight: bold; padding: 8px 20px; }"
        )
        self.btn_run.clicked.connect(self.run_forced_phot)
        ctrl_layout.addWidget(self.btn_run)

        self.btn_stop = QPushButton("Stop")
        self.btn_stop.setStyleSheet(
            "QPushButton { background-color: #f44336; color: white; font-weight: bold; padding: 8px 15px; }"
        )
        self.btn_stop.clicked.connect(self.stop_forced_phot)
        self.btn_stop.setEnabled(False)
        ctrl_layout.addWidget(self.btn_stop)

        ctrl_layout.addStretch()

        btn_log = QPushButton("Show Log")
        btn_log.setStyleSheet(
            "QPushButton { background-color: #607D8B; color: white; font-weight: bold; padding: 8px 15px; }"
        )
        btn_log.clicked.connect(self._show_log)
        ctrl_layout.addWidget(btn_log)
        self.content_layout.addLayout(ctrl_layout)

        # Progress
        prog_layout = QHBoxLayout()
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        prog_layout.addWidget(self.progress_bar)
        self.progress_label = QLabel("Ready")
        self.progress_label.setMinimumWidth(300)
        prog_layout.addWidget(self.progress_label)
        self.content_layout.addLayout(prog_layout)

        # Results table
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(7)
        self.results_table.setHorizontalHeaderLabels(
            ["File", "Filter", "Status", "N master", "N detected", "N forced", "Apcorr"]
        )
        self.results_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.results_table.horizontalHeader().setStretchLastSection(True)
        self.results_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.results_table.setMinimumHeight(200)
        self.content_layout.addWidget(self.results_table)

        # Log window (floating)
        self._log_win = QWidget(self, Qt.Window)
        self._log_win.setWindowTitle("Forced Phot Log")
        self._log_win.resize(700, 400)
        log_layout = QVBoxLayout(self._log_win)
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setFont(self.log_text.font())
        log_layout.addWidget(self.log_text)

        # Load existing results if step already done
        self._try_load_existing_results()

    def _check_prerequisites(self):
        refbuild_dir = step7_refbuild_dir(self.params.P.result_dir)
        wcs_dir = step6_wcs_dir(self.params.P.result_dir)
        has_wcs = wcs_dir.exists() and any(wcs_dir.glob("wcs_solve_summary.csv"))
        has_cat = refbuild_dir.exists() and any(refbuild_dir.glob("ref_catalog*.tsv"))
        parts = []
        if not has_wcs:
            parts.append("WCS (step6_wcs/) not found")
        if not has_cat:
            parts.append("Master catalog (step7_refbuild/) not found")
        if parts:
            self.status_label.setText("Missing: " + "; ".join(parts))
            self.status_label.setStyleSheet("QLabel { color: #f44336; }")
        else:
            self.status_label.setText("Prerequisites OK — ready to run")
            self.status_label.setStyleSheet("QLabel { color: #4CAF50; }")

    def _try_load_existing_results(self):
        idx_path = step_forced_phot_dir(self.params.P.result_dir) / "photometry_index.csv"
        if not idx_path.exists():
            return
        try:
            df = pd.read_csv(idx_path)
            self._populate_table(df.to_dict(orient="records"))
            self.progress_label.setText(f"Loaded {len(df)} frames from previous run")
            self.update_navigation_buttons()
        except Exception:
            pass

    def _populate_table(self, rows: list):
        self.results_table.setRowCount(len(rows))
        for ri, row in enumerate(rows):
            vals = [
                str(row.get("file", "")),
                str(row.get("filter", "")),
                str(row.get("status", "")),
                str(row.get("n_master", "")),
                str(row.get("n_detected", "")),
                str(row.get("n_forced", "")),
                f"{row.get('apcorr', float('nan')):.4f}" if pd.notna(row.get("apcorr")) else "—",
            ]
            for ci, val in enumerate(vals):
                item = QTableWidgetItem(val)
                item.setTextAlignment(Qt.AlignCenter)
                self.results_table.setItem(ri, ci, item)

    def _show_log(self):
        self._log_win.show()
        self._log_win.raise_()
        self._log_win.activateWindow()

    # ── Worker control ─────────────────────────────────────────────────────────

    def run_forced_phot(self):
        P = self.params.P
        try:
            result_dir = Path(P.result_dir)
            data_dir   = Path(P.data_dir)
            cache_dir  = Path(P.cache_dir) if not Path(P.cache_dir).is_absolute() else Path(P.cache_dir)
        except Exception as exc:
            QMessageBox.critical(self, "Config Error", f"Cannot read paths from params: {exc}")
            return

        file_list = []
        try:
            file_list = list(self.file_manager.get_file_list())
        except Exception:
            pass
        if not file_list:
            QMessageBox.warning(self, "No Files", "File list is empty. Run File Selection first.")
            return

        # Check for master catalog
        refbuild_dir = step7_refbuild_dir(result_dir)
        if not any(refbuild_dir.glob("ref_catalog*.tsv")):
            QMessageBox.warning(
                self, "No Master Catalog",
                "Master catalog not found in step7_refbuild/. Run Master Catalog Build first."
            )
            return

        self.log_text.clear()
        self.worker = ForcedPhotWorker(
            params=self.params,
            data_dir=data_dir,
            result_dir=result_dir,
            cache_dir=cache_dir if cache_dir.is_absolute() else result_dir / cache_dir,
            file_list=file_list,
        )
        self.worker.progress.connect(self._on_progress)
        self.worker.log.connect(self._on_log)
        self.worker.finished.connect(self._on_finished)
        self.worker.error.connect(self._on_error)
        self.worker.start()

        self.btn_run.setEnabled(False)
        self.btn_stop.setEnabled(True)
        self.progress_bar.setValue(0)
        self.progress_label.setText("Running…")

    def stop_forced_phot(self):
        if self.worker and self.worker.isRunning():
            self.worker.stop()
            self.worker.wait(5000)
        self.btn_run.setEnabled(True)
        self.btn_stop.setEnabled(False)
        self.progress_label.setText("Stopped")

    def _on_progress(self, current: int, total: int, fname: str):
        self.progress_bar.setMaximum(total)
        self.progress_bar.setValue(current)
        self.progress_label.setText(f"[{current}/{total}] {fname}")

    def _on_log(self, msg: str):
        self.log_text.append(msg)

    def _on_finished(self, results: dict):
        self.results = results
        rows = results.get("index_rows", [])
        self._populate_table(rows)
        self.btn_run.setEnabled(True)
        self.btn_stop.setEnabled(False)
        n_ok = sum(1 for r in rows if r.get("status") == "ok")
        self.progress_label.setText(f"Done — {n_ok}/{len(rows)} frames OK")
        self.update_navigation_buttons()

    def _on_error(self, error_type: str, msg: str):
        QMessageBox.critical(self, error_type, msg)
        self.btn_run.setEnabled(True)
        self.btn_stop.setEnabled(False)
        self.progress_label.setText("Error — see log")

    # ── StepWindowBase overrides ───────────────────────────────────────────────

    def validate_step(self) -> bool:
        idx_path = step_forced_phot_dir(self.params.P.result_dir) / "photometry_index.csv"
        return idx_path.exists()

    def restore_state(self):
        self._try_load_existing_results()
        self._check_prerequisites()
