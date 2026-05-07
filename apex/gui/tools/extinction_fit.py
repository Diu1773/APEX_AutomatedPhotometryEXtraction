"""
Extinction (Airmass Fit) Tool Window
Fits per-filter extinction coefficients using instrumental magnitudes.
"""

from __future__ import annotations

import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QTextEdit,
    QGroupBox, QMessageBox, QComboBox, QCheckBox, QSplitter, QTableWidget,
    QTableWidgetItem, QHeaderView, QDoubleSpinBox, QProgressBar, QDialog,
    QDialogButtonBox, QFormLayout, QSpinBox, QTabWidget, QFileDialog,
    QLineEdit, QAbstractItemView,
)

from apex.utils.astro_utils import compute_airmass_from_header
from apex.utils.photometry_utils import (
    get_numeric_array,
    MAG_ERR_COEFF,
    MAD_TO_SIGMA,
)
from apex.utils.io_utils import coerce_int64_source_id, load_file_path_map, read_csv_int64_source_id
from apex.utils.photometry_loader import load_frame_photometry
from apex.utils.qc_utils import load_frame_excludes as _load_frame_excludes
from apex.utils.step_paths_lc import (
    step2_cropped_dir,
    step5_wcs_dir,
    step6_refbuild_dir,
    step7_forced_phot_dir,
    step8_selection_dir,
    tool_extinction_dir,
)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

K2_SIGNIFICANCE = 0.005        # quadratic extinction term significance threshold
KCOLOR_SIGNIFICANCE = 0.01    # color-extinction term significance threshold

# Gaia DR3 -> SDSS polynomial transforms (Evans+ 2018, A&A 616, A4)
GAIA_TO_SDSS = {
    "g": {"coeffs": [0.2199, -0.6365, -0.1548, 0.0064], "color_range": (0.3, 3.0)},
    "r": {"coeffs": [-0.09837, 0.08592, 0.1907, -0.1701, 0.02263], "color_range": (0.0, 3.0)},
    "i": {"coeffs": [-0.293, 0.6404, -0.09609, -0.002104], "color_range": (0.5, 2.0)},
}


def _ext_group_key(date_val: str, filt: str) -> str:
    return f"{str(date_val)}::{str(filt).strip().lower()}"


def _finite_range(vals) -> float:
    arr = np.asarray(vals, float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.nan
    return float(np.nanmax(arr) - np.nanmin(arr))


class ExtinctionFitWorker(QThread):
    progress = pyqtSignal(int, int, str)
    log = pyqtSignal(str)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)

    def __init__(self, params, data_dir: Path, result_dir: Path, *,
                 mode: str = "ensemble", task: str = "all", phot_df: pd.DataFrame | None = None,
                 source_result_dir: Path | None = None,
                 selected_star_map: dict[str, list[int]] | None = None,
                 rejected_star_map: dict[str, list[int]] | None = None,
                 ap_scale: float = 1.0,
                 ann_in_scale: float = 4.0, ann_out_scale: float = 2.0):
        super().__init__()
        self.params = params
        self.data_dir = Path(data_dir)
        self.result_dir = Path(result_dir)
        self.source_result_dir = Path(source_result_dir) if source_result_dir is not None else Path(result_dir)
        self.mode = mode  # "ensemble" or "gaia"
        self.task = task  # "photometry", "fit", or "all"
        self.phot_df = phot_df
        self.selected_star_map = {
            str(k): sorted({int(v) for v in vals if v is not None})
            for k, vals in (selected_star_map or {}).items()
            if vals
        }
        self.rejected_star_map = {
            str(k): sorted({int(v) for v in vals if v is not None})
            for k, vals in (rejected_star_map or {}).items()
            if vals
        }
        self.ap_scale = ap_scale
        self.ann_in_scale = ann_in_scale
        self.ann_out_scale = ann_out_scale
        self._stop_requested = False
        self._source_file_path_map: dict[str, str] | None = None

    def stop(self):
        self._stop_requested = True

    def _log(self, msg: str):
        self.log.emit(msg)

    def _filter_variable_stars(
        self,
        sub: pd.DataFrame,
        scatter_col: str,
        method: str,
        sigma: float,
        min_frames: int,
        label: str = "",
    ) -> pd.DataFrame:
        """Remove variable stars based on per-star scatter in *scatter_col*.

        Returns the filtered DataFrame (or the original if nothing removed).
        """
        if sigma <= 0:
            return sub
        star_counts = sub.groupby("source_id")["file"].nunique()
        good_ids = star_counts[star_counts >= min_frames].index
        sub_var = sub[sub["source_id"].isin(good_ids)]
        if sub_var.empty:
            return sub
        if method == "std":
            star_scatter = sub_var.groupby("source_id")[scatter_col].std(ddof=1)
        else:
            star_scatter = sub_var.groupby("source_id")[scatter_col].apply(
                lambda s: float(np.nanmedian(np.abs(s - np.nanmedian(s)))))
        vals = star_scatter.to_numpy(float)
        med_s = float(np.nanmedian(vals)) if vals.size else np.nan
        mad_s = float(np.nanmedian(np.abs(vals - med_s))) if vals.size else np.nan
        sig_s = MAD_TO_SIGMA * mad_s if mad_s > 0 else float(np.nanstd(vals))
        if not (np.isfinite(sig_s) and sig_s > 0):
            return sub
        thresh = med_s + sigma * sig_s
        bad_ids = set(star_scatter[star_scatter > thresh].index)
        if bad_ids:
            before = sub["source_id"].nunique()
            sub = sub[~sub["source_id"].isin(bad_ids)].copy()
            after = sub["source_id"].nunique()
            self._log(f"{label} var-star filter ({method}, {sigma:.1f}σ): "
                      f"{before - after} stars removed (thr={thresh:.4f})")
        return sub

    @staticmethod
    def _as_float(value, default: float) -> float:
        try:
            if value is None:
                return float(default)
            v = float(value)
            return v if np.isfinite(v) else float(default)
        except Exception:
            return float(default)

    @staticmethod
    def _as_int(value, default: int) -> int:
        try:
            if value is None:
                return int(default)
            return int(value)
        except Exception:
            return int(default)

    @staticmethod
    def _pick_first_column(columns, candidates):
        colset = set(columns)
        for cand in candidates:
            if cand in colset:
                return cand
        return None

    @staticmethod
    def _normalize_path(value) -> str:
        try:
            return str(Path(value).expanduser().resolve())
        except Exception:
            return str(Path(value).expanduser())

    def _fit_clip_sigma(self) -> float:
        P = self.params.P
        value = getattr(P, "extfit_clip_sigma", None)
        if value is None:
            value = getattr(P, "zp_clip_sigma", 3.0)
        return self._as_float(value, 3.0)

    def _fit_iters(self) -> int:
        P = self.params.P
        value = getattr(P, "extfit_fit_iters", None)
        if value is None:
            value = getattr(P, "zp_fit_iters", 5)
        return max(1, self._as_int(value, 5))

    def _fit_min_points(self) -> int:
        P = self.params.P
        value = getattr(P, "extfit_min_points", None)
        if value is None:
            value = getattr(P, "min_master_gaia_matches", 10)
        return max(3, self._as_int(value, 10))

    def _extinction_input_paths(self) -> tuple[Path, Path]:
        out_dir = tool_extinction_dir(self.result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        return (
            out_dir / "step7_extinction_input.csv",
            out_dir / "ensemble_allstar_phot.csv",
        )

    def _source_path_map_data(self) -> dict[str, str]:
        if self._source_file_path_map is None:
            raw = None
            try:
                if self._normalize_path(self.source_result_dir) == self._normalize_path(self.result_dir):
                    raw = getattr(self.params.P, "file_path_map", None)
            except Exception:
                raw = None
            if isinstance(raw, dict) and raw:
                self._source_file_path_map = {str(k): str(v) for k, v in raw.items() if k and v}
            else:
                self._source_file_path_map = load_file_path_map(self.source_result_dir)
        return self._source_file_path_map

    def _selected_ids_for_group(self, date_val: str, filt: str) -> set[int]:
        key = _ext_group_key(date_val, filt)
        vals = self.selected_star_map.get(key, [])
        return {int(v) for v in vals if v is not None}

    def _rejected_ids_for_group(self, date_val: str, filt: str) -> set[int]:
        key = _ext_group_key(date_val, filt)
        vals = self.rejected_star_map.get(key, [])
        return {int(v) for v in vals if v is not None}

    def _resolve_source_fits_path(self, fname: str) -> Path | None:
        fname = str(fname).strip()
        if not fname:
            return None
        candidates: list[Path] = []
        mapped = self._source_path_map_data().get(fname)
        if mapped:
            candidates.append(Path(mapped))
        try:
            if self._normalize_path(self.source_result_dir) == self._normalize_path(self.result_dir):
                candidates.append(Path(self.params.get_file_path(fname)))
        except Exception:
            pass
        candidates.extend([
            self.data_dir / fname,
            self.source_result_dir / fname,
            step2_cropped_dir(self.source_result_dir) / fname,
        ])
        seen: set[str] = set()
        for cand in candidates:
            key = str(cand)
            if key in seen:
                continue
            seen.add(key)
            if cand.exists():
                return cand
        return None

    @staticmethod
    def _robust_median_and_err(arr):
        x = np.asarray(arr, float)
        x = x[np.isfinite(x)]
        if len(x) == 0:
            return (np.nan, np.nan, 0)
        med = float(np.median(x))
        mad = float(np.median(np.abs(x - med)))
        err = float(MAD_TO_SIGMA * mad / np.sqrt(max(len(x), 1)))
        return (med, err, int(len(x)))

    def _poly_eval(self, x, coeffs):
        x = np.asarray(x, float)
        y = np.zeros_like(x, dtype=float)
        p = np.ones_like(x, dtype=float)
        for a in coeffs:
            y += a * p
            p *= x
        return y

    def _robust_linfit(self, x, y, clip_sigma=3.0, iters=5, min_n=10):
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        m0 = np.isfinite(x) & np.isfinite(y)
        x = x[m0]
        y = y[m0]
        if len(x) < min_n:
            return (np.nan, np.nan, np.nan, 0, np.nan)
        k, zp = np.polyfit(x, y, 1)
        base_n = len(x)
        for _ in range(int(iters)):
            yhat = zp + k * x
            r = y - yhat
            med = np.nanmedian(r)
            mad = np.nanmedian(np.abs(r - med)) + 1e-12
            sig = MAD_TO_SIGMA * mad
            keep = np.abs(r - med) <= float(clip_sigma) * sig
            if keep.sum() < min_n:
                break
            if keep.sum() == len(x):
                break
            x, y = x[keep], y[keep]
            k, zp = np.polyfit(x, y, 1)
        yhat = zp + k * x
        scatter = float(np.nanstd(y - yhat)) if len(x) else np.nan
        outlier_frac = float(1.0 - (len(x) / max(base_n, 1)))
        return (float(k), float(zp), scatter, int(len(x)), outlier_frac)

    def _robust_weighted_linfit(self, x, y, w=None, clip_sigma=3.0, iters=5, min_n=10):
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        if w is not None:
            w = np.asarray(w, float)
            m0 = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
        else:
            m0 = np.isfinite(x) & np.isfinite(y)
        x = x[m0]
        y = y[m0]
        w = w[m0] if w is not None else None
        if len(x) < min_n:
            return (np.nan, np.nan, np.nan, 0, np.nan)

        def _fit(xx, yy, ww):
            if ww is None:
                k, zp = np.polyfit(xx, yy, 1)
            else:
                k, zp = np.polyfit(xx, yy, 1, w=ww)
            return k, zp

        k, zp = _fit(x, y, w)
        base_n = len(x)
        for _ in range(int(iters)):
            yhat = zp + k * x
            r = y - yhat
            med = np.nanmedian(r)
            mad = np.nanmedian(np.abs(r - med)) + 1e-12
            sig = MAD_TO_SIGMA * mad
            keep = np.abs(r - med) <= float(clip_sigma) * sig
            if keep.sum() < min_n:
                break
            if keep.sum() == len(x):
                break
            x = x[keep]
            y = y[keep]
            if w is not None:
                w = w[keep]
            k, zp = _fit(x, y, w)

        yhat = zp + k * x
        scatter = float(np.nanstd(y - yhat)) if len(x) else np.nan
        outlier_frac = float(1.0 - (len(x) / max(base_n, 1)))
        return (float(k), float(zp), scatter, int(len(x)), outlier_frac)

    def _fit_ensemble_once(
        self,
        sub: pd.DataFrame,
        clip_sigma: float,
        iters: int,
        min_n: int,
    ):
        """Fit m_ij = s_i + z_j + k1*X_j for one date+filter group."""
        if sub.empty:
            return None

        sub = sub.reset_index(drop=True)
        x = get_numeric_array(sub, "airmass")
        y = get_numeric_array(sub, "mag_inst")
        star_ids = coerce_int64_source_id(sub["source_id"]).astype("int64").to_numpy()
        frame_ids = sub["file"].astype(str).to_numpy()

        base_mask = np.isfinite(x) & np.isfinite(y)
        if base_mask.sum() < min_n:
            return None

        unique_stars = np.unique(star_ids[base_mask])
        unique_frames = np.unique(frame_ids[base_mask])
        n_star = len(unique_stars)
        n_frame = len(unique_frames)
        if n_star == 0 or n_frame == 0:
            return None

        star_index = {sid: i for i, sid in enumerate(unique_stars)}
        frame_index = {fid: i for i, fid in enumerate(unique_frames)}
        ref_frame = unique_frames[0]
        ref_idx = frame_index[ref_frame]

        star_idx = np.array([star_index[s] for s in star_ids], dtype=int)
        frame_idx = np.array([frame_index[f] for f in frame_ids], dtype=int)

        # SNR-based weights (if available)
        has_snr = "snr" in sub.columns
        if has_snr:
            snr_arr = get_numeric_array(sub, "snr", default=1.0)
            sig_mag = MAG_ERR_COEFF / np.clip(snr_arr, 1e-6, None)
            w_all = 1.0 / np.clip(sig_mag ** 2, 1e-12, None)
        else:
            w_all = np.ones(len(sub), dtype=float)

        def build_matrix(mask):
            idx = np.where(mask)[0]
            n_obs = len(idx)
            n_params = n_star + (n_frame - 1) + 1
            A = np.zeros((n_obs, n_params), dtype=float)
            rows = np.arange(n_obs)
            A[rows, star_idx[idx]] = 1.0
            if n_frame > 1:
                fidx = frame_idx[idx]
                m = fidx != ref_idx
                if m.any():
                    fcol = fidx[m].copy()
                    fcol = np.where(fcol > ref_idx, fcol - 1, fcol)
                    A[rows[m], n_star + fcol] = 1.0
            A[:, -1] = x[idx]
            return A, y[idx], idx

        mask_fit = base_mask.copy()
        base_n = int(mask_fit.sum())
        if base_n < min_n:
            return None

        for _ in range(int(iters)):
            A, yv, idx = build_matrix(mask_fit)
            if len(yv) < min_n:
                break
            W = np.sqrt(w_all[idx])
            coef, _, _, _ = np.linalg.lstsq(W[:, None] * A, W * yv, rcond=None)
            r = yv - (A @ coef)
            med = np.nanmedian(r)
            mad = np.nanmedian(np.abs(r - med)) + 1e-12
            sig = MAD_TO_SIGMA * mad
            keep = np.abs(r - med) <= float(clip_sigma) * sig
            if keep.sum() < min_n:
                break
            if keep.sum() == len(yv):
                break
            new_mask = mask_fit.copy()
            new_mask[idx] = keep
            mask_fit = new_mask

        A, yv, idx = build_matrix(mask_fit)
        if len(yv) < min_n:
            return None
        W = np.sqrt(w_all[idx])
        coef, _, _, _ = np.linalg.lstsq(W[:, None] * A, W * yv, rcond=None)

        s_vals = coef[:n_star]
        z_vals = np.zeros(n_frame, dtype=float)
        if n_frame > 1:
            frame_cols = coef[n_star:-1]
            for fi in range(n_frame):
                if fi == ref_idx:
                    z_vals[fi] = 0.0
                else:
                    col = fi if fi < ref_idx else fi - 1
                    z_vals[fi] = frame_cols[col]
        k1 = float(coef[-1])

        # residuals for all points in sub
        yhat = s_vals[star_idx] + z_vals[frame_idx] + k1 * x
        resid = y - yhat
        delta_m = y - (s_vals[star_idx] + z_vals[frame_idx])
        scatter = float(np.nanstd(resid[base_mask])) if base_n else np.nan
        outlier_frac = float(1.0 - (mask_fit.sum() / max(base_n, 1)))

        return {
            "k1": k1,
            "s_vals": s_vals,
            "z_vals": z_vals,
            "star_index": star_index,
            "frame_index": frame_index,
            "ref_frame": ref_frame,
            "resid": resid,
            "delta_m": delta_m,
            "scatter": scatter,
            "n_used": int(mask_fit.sum()),
            "outlier_frac": outlier_frac,
            "n_star": n_star,
            "n_frame": n_frame,
        }

    def _robust_quadfit(self, x, y, clip_sigma=3.0, iters=5, min_n=15):
        """2차 다항식 피팅: y = k2*x^2 + k1*x + zp"""
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        m0 = np.isfinite(x) & np.isfinite(y)
        x = x[m0]
        y = y[m0]
        if len(x) < min_n:
            return (np.nan, np.nan, np.nan, np.nan, 0, np.nan, False)
        try:
            coeffs = np.polyfit(x, y, 2)  # [k2, k1, zp]
            k2, k1, zp = coeffs
        except Exception:
            return (np.nan, np.nan, np.nan, np.nan, 0, np.nan, False)
        base_n = len(x)
        for _ in range(int(iters)):
            yhat = zp + k1 * x + k2 * x * x
            r = y - yhat
            med = np.nanmedian(r)
            mad = np.nanmedian(np.abs(r - med)) + 1e-12
            sig = MAD_TO_SIGMA * mad
            keep = np.abs(r - med) <= float(clip_sigma) * sig
            if keep.sum() < min_n:
                break
            if keep.sum() == len(x):
                break
            x, y = x[keep], y[keep]
            try:
                coeffs = np.polyfit(x, y, 2)
                k2, k1, zp = coeffs
            except Exception:
                break
        yhat = zp + k1 * x + k2 * x * x
        scatter = float(np.nanstd(y - yhat)) if len(x) else np.nan
        outlier_frac = float(1.0 - (len(x) / max(base_n, 1)))
        # k2 significance check
        k2_significant = abs(k2) > K2_SIGNIFICANCE if np.isfinite(k2) else False
        return (float(k1), float(k2), float(zp), scatter, int(len(x)), outlier_frac, k2_significant)

    def _robust_color_extinction_fit(self, X, C, y, clip_sigma=3.0, iters=5, min_n=15):
        """
        색 의존 소광 피팅: delta = k' * X + k'' * C * X + zp

        Parameters:
        - X: airmass 배열
        - C: 색지수 배열 (예: g-r, BP-RP)
        - y: delta (m_ref - m_inst) 배열

        Returns:
        - k_prime: 평균 소광 계수
        - k_double_prime: 색 의존 소광 계수
        - zp: 제로포인트
        - scatter, n_used, outlier_frac, significant
        """
        X = np.asarray(X, float)
        C = np.asarray(C, float)
        y = np.asarray(y, float)

        m0 = np.isfinite(X) & np.isfinite(C) & np.isfinite(y)
        X = X[m0]
        C = C[m0]
        y = y[m0]

        if len(X) < min_n:
            return (np.nan, np.nan, np.nan, np.nan, 0, np.nan, False)

        # 다중 선형 회귀: y = a0 + a1*X + a2*(C*X)
        # Design matrix: [1, X, C*X]
        try:
            A = np.column_stack([np.ones_like(X), X, C * X])
            coeffs, residuals, rank, s = np.linalg.lstsq(A, y, rcond=None)
            zp, k_prime, k_double_prime = coeffs
        except Exception:
            return (np.nan, np.nan, np.nan, np.nan, 0, np.nan, False)

        base_n = len(X)
        for _ in range(int(iters)):
            yhat = zp + k_prime * X + k_double_prime * C * X
            r = y - yhat
            med = np.nanmedian(r)
            mad = np.nanmedian(np.abs(r - med)) + 1e-12
            sig = MAD_TO_SIGMA * mad
            keep = np.abs(r - med) <= float(clip_sigma) * sig
            if keep.sum() < min_n:
                break
            if keep.sum() == len(X):
                break
            X, C, y = X[keep], C[keep], y[keep]
            try:
                A = np.column_stack([np.ones_like(X), X, C * X])
                coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
                zp, k_prime, k_double_prime = coeffs
            except Exception:
                break

        yhat = zp + k_prime * X + k_double_prime * C * X
        scatter = float(np.nanstd(y - yhat)) if len(X) else np.nan
        outlier_frac = float(1.0 - (len(X) / max(base_n, 1)))

        # k''가 유의미한지 체크
        k_double_significant = abs(k_double_prime) > KCOLOR_SIGNIFICANCE if np.isfinite(k_double_prime) else False

        return (float(k_prime), float(k_double_prime), float(zp), scatter, int(len(X)), outlier_frac, k_double_significant)

    # ------------------------------------------------------------------
    # Shared input build: Step 7 forced photometry + frame airmass
    # ------------------------------------------------------------------

    def _run_step7_source_load(self) -> pd.DataFrame:
        """Load Step 7 forced photometry tables and enrich them for extinction fitting."""
        result_dir = self.source_result_dir
        idx_path = step7_forced_phot_dir(result_dir) / "photometry_index.csv"
        if not idx_path.exists():
            raise FileNotFoundError("Step 7 forced photometry_index.csv not found. Run Step 7 forced photometry first.")

        idx = pd.read_csv(idx_path)
        if "file" not in idx.columns:
            for cand in ("fname", "frame", "image", "fits", "name"):
                if cand in idx.columns:
                    idx = idx.rename(columns={cand: "file"})
                    break
        if "file" not in idx.columns:
            raise RuntimeError("photometry_index.csv is missing a file column")

        idx["file"] = idx["file"].astype(str).str.strip()
        idx = idx[idx["file"] != ""].copy()
        if "filter" in idx.columns:
            idx["filter"] = idx["filter"].astype(str).str.strip().str.lower()
        elif "FILTER" in idx.columns:
            idx["filter"] = idx["FILTER"].astype(str).str.strip().str.lower()
        else:
            idx["filter"] = ""

        qc_exclude: set[str] = set()
        for qp in [step7_forced_phot_dir(result_dir) / "frame_quality.csv", result_dir / "frame_quality.csv"]:
            if not qp.exists():
                continue
            try:
                dfq = pd.read_csv(qp)
            except Exception:
                continue
            if {"file", "passed"} <= set(dfq.columns):
                bad = dfq.loc[dfq["passed"] == False, "file"].astype(str)
                qc_exclude |= set(bad.tolist())
            break

        if qc_exclude:
            before = len(idx)
            idx = idx[~idx["file"].isin(qc_exclude)].reset_index(drop=True)
            self._log(f"Frame QC: excluded {before - len(idx)} frame(s) from Step 7 quality flags")

        step10_excludes = _load_frame_excludes(result_dir)
        if step10_excludes:
            before = len(idx)
            idx = idx[~idx["file"].isin(set(step10_excludes.keys()))].reset_index(drop=True)
            self._log(f"Frame QC: excluded {before - len(idx)} frame(s) from manual frame exclusion")

        if idx.empty:
            raise RuntimeError("No Step 7 frames remain after exclusions")

        frame_airmass = self._build_frame_airmass(idx)
        if frame_airmass is None or frame_airmass.empty:
            raise RuntimeError("No usable airmass could be resolved from Step 7 frames")

        frame_airmass = frame_airmass.drop_duplicates(subset=["file"], keep="last").reset_index(drop=True)
        airmass_map = dict(zip(frame_airmass["file"].astype(str), pd.to_numeric(frame_airmass["airmass"], errors="coerce")))
        date_map = dict(
            zip(
                frame_airmass["file"].astype(str),
                frame_airmass.get("date", pd.Series(["unknown"] * len(frame_airmass))).astype(str),
            )
        )
        filter_map = dict(
            zip(
                frame_airmass["file"].astype(str),
                frame_airmass.get("filter", pd.Series([""] * len(frame_airmass))).astype(str),
            )
        )

        rows: list[pd.DataFrame] = []
        total = len(idx)
        missing_tables = 0
        missing_mag_column = 0
        dropped_bad_mag = 0
        dropped_saturated = 0
        missing_source_frames = 0

        for i, row in idx.iterrows():
            if self._stop_requested:
                break

            fname = str(row.get("file", "")).strip()
            filt_hint = str(row.get("filter", "")).strip().lower()
            dfp = load_frame_photometry(result_dir, fname, filt_hint)
            if dfp is None or dfp.empty:
                missing_tables += 1
                self.progress.emit(i + 1, total, fname)
                continue

            dfp = dfp.copy()
            mag_col = self._pick_first_column(dfp.columns, ("mag_inst", "mag", "mag_ap", "mag_apcorr"))
            if mag_col is None:
                missing_mag_column += 1
                self._log(f"[WARN] {fname}: no usable magnitude column in Step 7 forced photometry")
                self.progress.emit(i + 1, total, fname)
                continue

            err_col = self._pick_first_column(dfp.columns, ("mag_err", "emag", "emag_inst", "magerr"))
            snr_col = self._pick_first_column(dfp.columns, ("snr",))
            filter_col = self._pick_first_column(dfp.columns, ("FILTER", "filter"))

            filter_value = str(filter_map.get(fname, filt_hint)).strip().lower()
            if filter_col is not None:
                valid_filter = dfp[filter_col].dropna()
                if len(valid_filter):
                    filter_value = str(valid_filter.iloc[0]).strip().lower()

            frame = pd.DataFrame(index=dfp.index.copy())
            frame["file"] = fname
            frame["filter"] = filter_value
            frame["date"] = str(date_map.get(fname, self._extract_date_from_file(fname)))
            frame["airmass"] = float(airmass_map.get(fname, np.nan))
            frame["mag_inst"] = pd.to_numeric(dfp[mag_col], errors="coerce")
            frame["mag_err"] = pd.to_numeric(dfp[err_col], errors="coerce") if err_col is not None else np.nan
            frame["snr"] = pd.to_numeric(dfp[snr_col], errors="coerce") if snr_col is not None else np.nan

            if "source_id" in dfp.columns:
                frame["source_id"] = coerce_int64_source_id(dfp["source_id"]).astype("Int64")
            else:
                frame["source_id"] = pd.Series(pd.array([pd.NA] * len(dfp), dtype="Int64"), index=dfp.index)
            if frame["source_id"].isna().all():
                missing_source_frames += 1

            if "ID" in dfp.columns:
                frame["ID"] = pd.to_numeric(dfp["ID"], errors="coerce").astype("Int64")
            else:
                frame["ID"] = pd.Series(pd.array([pd.NA] * len(dfp), dtype="Int64"), index=dfp.index)

            bad_signal = np.zeros(len(frame), dtype=bool)
            for bad_col in ("is_saturated", "is_nonlinear"):
                if bad_col not in dfp.columns:
                    continue
                series = dfp[bad_col]
                if pd.api.types.is_bool_dtype(series):
                    bad_signal |= series.fillna(False).to_numpy(dtype=bool)
                else:
                    bad_signal |= pd.to_numeric(series, errors="coerce").fillna(0).to_numpy(dtype=float) > 0
            if bad_signal.any():
                dropped_saturated += int(bad_signal.sum())
                frame = frame.loc[~bad_signal].copy()

            bad_mag = ~np.isfinite(frame["mag_inst"].to_numpy(float))
            if bad_mag.any():
                dropped_bad_mag += int(bad_mag.sum())
                frame = frame.loc[~bad_mag].copy()

            if not frame.empty:
                rows.append(frame.reset_index(drop=True))
            self.progress.emit(i + 1, total, fname)

        if self._stop_requested:
            self.finished.emit({"stopped": True})
            return pd.DataFrame()

        if not rows:
            raise RuntimeError("No usable Step 7 forced photometry rows were loaded")

        phot_df = pd.concat(rows, ignore_index=True)
        phot_df["source_workspace"] = self._normalize_path(result_dir)
        finite_airmass = int(np.isfinite(pd.to_numeric(phot_df["airmass"], errors="coerce")).sum())
        source_rows = int(phot_df["source_id"].notna().sum()) if "source_id" in phot_df.columns else 0

        self._log(
            "Step 7 source load complete: "
            f"{len(phot_df)} rows, {phot_df['file'].nunique()} frames, "
            f"{source_rows} matched-source rows"
        )
        if missing_tables:
            self._log(f"Skipped {missing_tables} frame(s) with missing/empty Step 7 forced photometry tables")
        if missing_mag_column:
            self._log(f"Skipped {missing_mag_column} frame(s) without a usable magnitude column")
        if dropped_bad_mag:
            self._log(f"Dropped {dropped_bad_mag} row(s) with NaN instrumental magnitude")
        if dropped_saturated:
            self._log(f"Dropped {dropped_saturated} saturated/nonlinear row(s) from Step 7")
        if missing_source_frames:
            self._log(f"Warning: {missing_source_frames} frame(s) had no resolved source_id values")
        if finite_airmass == 0:
            raise RuntimeError("All loaded Step 7 rows have NaN airmass")

        return phot_df

    def _run_ensemble_fit(self, phot_df: pd.DataFrame):
        """Ensemble extinction fit with frame zeropoint separation.

        Model per date+filter:
            m_ij = s_i + z_j + k1 * X_j + eps
        where s_i is star constant, z_j is frame zeropoint, and k1 is extinction.
        """
        P = self.params.P
        clip_sigma = self._fit_clip_sigma()
        fit_iters = self._fit_iters()
        min_match = self._fit_min_points()
        snr_cut = self._as_float(getattr(P, "extinction_snr_min", 10.0), 10.0)

        df = phot_df.copy()
        if "source_id" in df.columns:
            before_sid = len(df)
            df = df[df["source_id"].notna()].copy()
            if before_sid != len(df):
                self._log(f"Dropped {before_sid - len(df)} rows without source_id")

        # SNR filter
        if "snr" in df.columns:
            n_before = len(df)
            df = df[df["snr"] >= snr_cut]
            self._log(f"SNR filter (>={snr_cut}): {n_before} → {len(df)}")

        # Remove stars with too few frames (< 3)
        star_counts = df.groupby("source_id")["file"].nunique()
        good_stars = set(star_counts[star_counts >= 3].index)
        df = df[df["source_id"].isin(good_stars)]
        self._log(f"Stars with >=3 frames: {len(good_stars)}")

        if df.empty:
            raise RuntimeError("No data after quality filters")

        # Also reject frames with no valid airmass
        n_before_am = len(df)
        nan_am_frames = df.loc[~np.isfinite(df["airmass"]), "file"].nunique()
        df = df[np.isfinite(df["airmass"])].copy()
        if nan_am_frames > 0:
            self._log(f"Airmass NaN: dropped {n_before_am - len(df)} points from {nan_am_frames} frames")
        if df.empty:
            raise RuntimeError("No data after airmass filtering (all frames had NaN airmass)")
        self._log(f"After SNR/airmass filters: {len(df)} points, {df['file'].nunique()} frames, "
                  f"{df['source_id'].nunique()} stars")

        # Step 3: fit per date+filter
        fit_rows = []
        point_rows = []

        dates_seen = sorted(df["date"].unique())
        filters_seen = sorted(df["filter"].unique())
        self._log(f"Dates: {dates_seen}")
        self._log(f"Filters: {filters_seen}")

        var_method_raw = getattr(P, "extinction_varstar_method", "mad")
        var_method = str(var_method_raw if var_method_raw else "mad").strip().lower()
        var_sigma = self._as_float(getattr(P, "extinction_varstar_sigma", 3.0), 3.0)
        var_min_frames = self._as_int(getattr(P, "extinction_varstar_min_frames", 5), 5)

        qc_method_raw = getattr(P, "extinction_frame_qc_method", "mad")
        qc_method = str(qc_method_raw if qc_method_raw else "mad").strip().lower()
        qc_sigma = self._as_float(getattr(P, "extinction_frame_qc_sigma", 3.0), 3.0)

        groups = list(df.groupby(["date", "filter"]))
        n_groups = len(groups)
        for i_grp, ((date_val, filt), sub) in enumerate(groups):
            self.progress.emit(i_grp, n_groups, f"Ensemble {date_val}/{filt}")
            sub = sub.copy()
            n_pts = len(sub)

            use_dx = bool(getattr(P, "extinction_delta_x_enable", True))
            dx_min = self._as_float(getattr(P, "extinction_delta_x_min", 0.3), 0.3)
            if use_dx:
                x = get_numeric_array(sub, "airmass")
                dx = float(np.nanmax(x) - np.nanmin(x)) if np.isfinite(x).any() else np.nan
                if dx < dx_min:
                    self._log(f"[{date_val}][{filt}] skipped: ΔX={dx:.3f} < {dx_min:.3f}")
                    continue

            if n_pts < min_match:
                self._log(f"[{date_val}][{filt}] skipped: only {n_pts} points (min={min_match})")
                continue

            fit = self._fit_ensemble_once(sub, clip_sigma, fit_iters, min_match)
            if fit is None or not np.isfinite(fit.get("k1", np.nan)):
                self._log(f"[{date_val}][{filt}] FAILED: fit did not converge")
                continue

            sub["resid"] = fit["resid"]
            sub["delta_m"] = fit["delta_m"]

            # Variable star filter (residual-based)
            if var_sigma > 0:
                n_before = sub["source_id"].nunique()
                sub = self._filter_variable_stars(
                    sub, "resid", var_method, var_sigma, var_min_frames,
                    label=f"[{date_val}][{filt}]")
                if sub["source_id"].nunique() < n_before:
                    fit = self._fit_ensemble_once(sub, clip_sigma, fit_iters, min_match)
                    if fit is None:
                        continue
                    sub["resid"] = fit["resid"]
                    sub["delta_m"] = fit["delta_m"]

            # Frame QC (residual median)
            if qc_sigma > 0:
                frame_stats = sub.groupby("file").agg(
                    resid_med=("resid", "median"),
                    resid_mad=("resid", lambda s: float(np.nanmedian(np.abs(s - np.nanmedian(s))))),
                    n_stars=("resid", "count"),
                )
                if qc_method == "std":
                    global_med = float(np.nanmean(frame_stats["resid_med"]))
                    global_sig = float(np.nanstd(frame_stats["resid_med"]))
                else:
                    global_med = float(np.nanmedian(frame_stats["resid_med"]))
                    global_mad = float(np.nanmedian(np.abs(frame_stats["resid_med"] - global_med)))
                    global_sig = MAD_TO_SIGMA * global_mad if global_mad > 0 else float(np.nanstd(frame_stats["resid_med"]))
                if global_sig > 0:
                    bad_mask = np.abs(frame_stats["resid_med"] - global_med) > qc_sigma * global_sig
                    bad_frames = set(frame_stats.index[bad_mask])
                    if bad_frames:
                        self._log(f"[{date_val}][{filt}] Frame QC ({qc_method}, {qc_sigma:.1f}σ): "
                                  f"{len(bad_frames)} frames removed")
                        sub = sub[~sub["file"].isin(bad_frames)].copy()
                        fit = self._fit_ensemble_once(sub, clip_sigma, fit_iters, min_match)
                        if fit is None:
                            continue
                        sub["resid"] = fit["resid"]
                        sub["delta_m"] = fit["delta_m"]

            if sub.empty:
                continue

            k1 = float(fit["k1"])
            scatter = float(fit.get("scatter", np.nan))
            n_used = int(fit.get("n_used", len(sub)))
            out_frac = float(fit.get("outlier_frac", np.nan))
            n_star = int(sub["source_id"].nunique())
            n_frame = int(sub["file"].nunique())
            m0 = float(np.nanmean(fit.get("z_vals"))) if fit.get("z_vals") is not None else np.nan

            fit_rows.append({
                "date": date_val,
                "filter": filt,
                "k1": k1,
                "k2": 0.0,
                "zp": m0,
                "m0": m0,
                "scatter": scatter,
                "n_total": n_pts,
                "n_used": n_used,
                "outlier_fraction": out_frac,
                "fit_order": 1,
                "n_stars": n_star,
                "n_frames": n_frame,
                "method": "ensemble",
            })

            for _, r in sub.iterrows():
                point_rows.append({
                    "date": date_val, "filter": filt,
                    "airmass": float(r["airmass"]),
                    "delta_m": float(r["delta_m"]),
                    "resid": float(r["resid"]),
                })

            self._log(f"[{date_val}][{filt}] k1={k1:.5f}, scatter={scatter:.4f}, "
                      f"n={n_used}/{n_pts} ({n_star} stars, {n_frame} frames)")

        if not fit_rows:
            raise RuntimeError("No valid ensemble extinction fits produced")

        fit_df = pd.DataFrame(fit_rows)
        points_df = pd.DataFrame(point_rows)

        # Save results
        out_dir = tool_extinction_dir(self.result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        fit_path = out_dir / "ensemble_extinction_by_filter.csv"
        pts_path = out_dir / "ensemble_phot_points.csv"
        fit_df.to_csv(fit_path, index=False)
        points_df.to_csv(pts_path, index=False)
        self._log(f"Saved {fit_path.name} ({len(fit_df)} rows)")
        self._log(f"Saved {pts_path.name} ({len(points_df)} rows)")

        # Shared filename for downstream loaders
        compat_fit = out_dir / "extinction_fit_by_filter.csv"
        try:
            fit_df.to_csv(compat_fit, index=False)
            self._log(f"Saved {compat_fit.name} ({len(fit_df)} rows)")
        except Exception:
            self._log("Warning: failed to write extinction_fit_by_filter.csv")

        # Also refresh the cached Step 7 extinction input table
        self._save_photometry(phot_df)

        self.finished.emit({
            "fit": fit_df, "points": points_df,
            "mode": "ensemble",
        })

    def _run_per_star_fit(self, phot_df: pd.DataFrame):
        """Per-star Bouguer line fitting (k_i), then aggregate to k1."""
        P = self.params.P
        clip_sigma = self._fit_clip_sigma()
        fit_iters = self._fit_iters()
        min_points = self._fit_min_points()
        snr_cut = self._as_float(getattr(P, "extinction_snr_min", 10.0), 10.0)

        star_min_frames = self._as_int(getattr(P, "extinction_star_min_frames", 8), 8)
        star_rms_max = self._as_float(getattr(P, "extinction_star_rms_max", 0.10), 0.10)
        star_snr_med_min = self._as_float(getattr(P, "extinction_star_snr_med_min", 10.0), 10.0)
        min_good_stars = self._as_int(getattr(P, "extinction_min_good_stars", 3), 3)
        use_weights = bool(getattr(P, "extinction_star_use_weights", True))

        use_dx = bool(getattr(P, "extinction_delta_x_enable", True))
        dx_min_global = self._as_float(getattr(P, "extinction_delta_x_min", 0.3), 0.3)

        var_method_raw = getattr(P, "extinction_varstar_method", "mad")
        var_method = str(var_method_raw if var_method_raw else "mad").strip().lower()
        var_sigma = self._as_float(getattr(P, "extinction_varstar_sigma", 3.0), 3.0)
        var_min_frames = self._as_int(getattr(P, "extinction_varstar_min_frames", 5), 5)

        df = phot_df.copy()
        if "source_id" in df.columns:
            before_sid = len(df)
            df = df[df["source_id"].notna()].copy()
            if before_sid != len(df):
                self._log(f"Dropped {before_sid - len(df)} rows without source_id")
        if "snr" in df.columns:
            n_before = len(df)
            df = df[df["snr"] >= snr_cut]
            self._log(f"SNR filter (>={snr_cut}): {n_before} → {len(df)}")

        nan_am = ~np.isfinite(df["airmass"])
        nan_mag = ~np.isfinite(df["mag_inst"])
        n_drop = int((nan_am | nan_mag).sum())
        if n_drop > 0:
            self._log(f"Dropped {n_drop} points with NaN airmass/mag_inst "
                      f"({int(nan_am.sum())} airmass, {int(nan_mag.sum())} mag)")
        df = df[~nan_am & ~nan_mag].copy()
        if df.empty:
            raise RuntimeError("No data after SNR/airmass filters")

        fit_rows = []
        point_rows = []
        star_rows = []
        frame_rows = []

        groups = list(df.groupby(["date", "filter"]))
        n_groups = len(groups)
        for i_grp, ((date_val, filt), sub) in enumerate(groups):
            self.progress.emit(i_grp, n_groups, f"Per-star {date_val}/{filt}")
            sub = sub.copy()
            n_pts = len(sub)

            selected_ids = self._selected_ids_for_group(date_val, filt)
            rejected_ids = self._rejected_ids_for_group(date_val, filt)
            if rejected_ids:
                before_rej = int(sub["source_id"].nunique())
                sub = sub[~sub["source_id"].isin(rejected_ids)].copy()
                after_rej = int(sub["source_id"].nunique())
                self._log(
                    f"[{date_val}][{filt}] manual reject: "
                    f"{before_rej - after_rej} stars removed"
                )
                if sub.empty:
                    self._log(f"[{date_val}][{filt}] skipped: no rows after manual reject")
                    continue
            if selected_ids:
                before_sel = int(sub["source_id"].nunique())
                sub = sub[sub["source_id"].isin(selected_ids)].copy()
                after_sel = int(sub["source_id"].nunique())
                self._log(
                    f"[{date_val}][{filt}] manual selection: "
                    f"{after_sel}/{before_sel} stars kept"
                )
                if sub.empty:
                    self._log(f"[{date_val}][{filt}] skipped: no rows after manual selection")
                    continue

            if use_dx:
                x_all = get_numeric_array(sub, "airmass")
                if np.isfinite(x_all).any():
                    dx_all = float(np.nanmax(x_all) - np.nanmin(x_all))
                    if dx_all < dx_min_global:
                        self._log(f"[{date_val}][{filt}] skipped: ΔX={dx_all:.3f} < {dx_min_global:.3f}")
                        continue

            if n_pts < min_points:
                self._log(f"[{date_val}][{filt}] skipped: only {n_pts} points (min={min_points})")
                continue

            # Variable star pre-filter (mag_inst scatter)
            if var_sigma > 0:
                sub = self._filter_variable_stars(
                    sub, "mag_inst", var_method, var_sigma, var_min_frames,
                    label=f"[{date_val}][{filt}]")

            k_list = []
            m0_map: dict[int, float] = {}
            star_rms: dict[int, float] = {}
            stats = {"total": 0, "frames": 0, "snr": 0, "fit": 0, "rms": 0}

            for sid, ssub in sub.groupby("source_id"):
                sid_int = int(sid)
                stats["total"] += 1
                n_i = len(ssub)
                med_snr = float(np.nanmedian(get_numeric_array(ssub, "snr"))) if "snr" in ssub.columns else np.nan
                x = get_numeric_array(ssub, "airmass")
                y = get_numeric_array(ssub, "mag_inst")
                dx_i = _finite_range(x)
                reject_reason = ""
                if n_i < star_min_frames:
                    reject_reason = "min_frames"
                    star_rows.append({
                        "date": date_val, "filter": filt, "source_id": sid_int,
                        "n_frames": int(n_i), "delta_x": dx_i, "snr_med": med_snr,
                        "k_i": np.nan, "m0_i": np.nan, "rms_i": np.nan,
                        "used": False, "reject_reason": reject_reason,
                    })
                    continue
                stats["frames"] += 1

                if "snr" in ssub.columns:
                    if np.isfinite(star_snr_med_min) and star_snr_med_min > 0 and med_snr < star_snr_med_min:
                        reject_reason = "median_snr"
                        star_rows.append({
                            "date": date_val, "filter": filt, "source_id": sid_int,
                            "n_frames": int(n_i), "delta_x": dx_i, "snr_med": med_snr,
                            "k_i": np.nan, "m0_i": np.nan, "rms_i": np.nan,
                            "used": False, "reject_reason": reject_reason,
                        })
                        continue
                stats["snr"] += 1

                if not np.isfinite(x).any() or not np.isfinite(y).any():
                    reject_reason = "nan_xy"
                    star_rows.append({
                        "date": date_val, "filter": filt, "source_id": sid_int,
                        "n_frames": int(n_i), "delta_x": dx_i, "snr_med": med_snr,
                        "k_i": np.nan, "m0_i": np.nan, "rms_i": np.nan,
                        "used": False, "reject_reason": reject_reason,
                    })
                    continue

                w = None
                if use_weights and "snr" in ssub.columns:
                    snr = get_numeric_array(ssub, "snr")
                    sig_mag = MAG_ERR_COEFF / np.clip(snr, 1e-6, None)
                    w = 1.0 / np.clip(sig_mag ** 2, 1e-12, None)

                k_i, m0_i, rms_i, n_used, _ = self._robust_weighted_linfit(
                    x, y, w=w, clip_sigma=clip_sigma, iters=fit_iters, min_n=star_min_frames
                )
                if not np.isfinite(k_i):
                    reject_reason = "fit_failed"
                    star_rows.append({
                        "date": date_val, "filter": filt, "source_id": sid_int,
                        "n_frames": int(n_i), "delta_x": dx_i, "snr_med": med_snr,
                        "k_i": np.nan, "m0_i": np.nan, "rms_i": np.nan,
                        "used": False, "reject_reason": reject_reason,
                    })
                    continue
                stats["fit"] += 1
                if np.isfinite(star_rms_max) and star_rms_max > 0 and np.isfinite(rms_i) and rms_i > star_rms_max:
                    reject_reason = "rms_max"
                    star_rows.append({
                        "date": date_val, "filter": filt, "source_id": sid_int,
                        "n_frames": int(n_i), "delta_x": dx_i, "snr_med": med_snr,
                        "k_i": float(k_i), "m0_i": float(m0_i), "rms_i": float(rms_i),
                        "used": False, "reject_reason": reject_reason,
                    })
                    continue
                stats["rms"] += 1

                k_list.append(float(k_i))
                m0_map[sid_int] = float(m0_i)
                star_rms[sid_int] = float(rms_i) if np.isfinite(rms_i) else np.nan
                star_rows.append({
                    "date": date_val, "filter": filt, "source_id": sid_int,
                    "n_frames": int(n_i), "delta_x": dx_i, "snr_med": med_snr,
                    "k_i": float(k_i), "m0_i": float(m0_i), "rms_i": float(rms_i) if np.isfinite(rms_i) else np.nan,
                    "used": True, "reject_reason": "",
                })

            if len(k_list) < min_good_stars:
                dx_note = f"ΔX≥{dx_min_global:.2f}" if use_dx else "ΔX off"
                self._log(
                    f"[{date_val}][{filt}] skipped: good stars={len(k_list)} < {min_good_stars} "
                    f"(min_frames={star_min_frames}, {dx_note}, RMS≤{star_rms_max:.2f}, "
                    f"medianSNR≥{star_snr_med_min:.1f}) | "
                    f"stats: total={stats['total']}, frames_ok={stats['frames']}, "
                    f"snr_ok={stats['snr']}, fit_ok={stats['fit']}, rms_ok={stats['rms']}"
                )
                continue

            k_arr = np.asarray(k_list, float)
            k1 = float(np.nanmedian(k_arr))
            k_mad = float(np.nanmedian(np.abs(k_arr - k1)))
            k1_err = float(MAD_TO_SIGMA * k_mad / np.sqrt(max(len(k_arr), 1)))

            # Build points for plotting using global k1
            good_ids = set(m0_map.keys())
            sub_good = sub[sub["source_id"].isin(good_ids)].copy()
            sub_good["m0_i"] = sub_good["source_id"].map(m0_map)
            sub_good["delta_m"] = sub_good["mag_inst"] - sub_good["m0_i"]
            sub_good["resid"] = sub_good["mag_inst"] - (sub_good["m0_i"] + k1 * sub_good["airmass"])

            frame_summary = sub_good.groupby("file").agg(
                date=("date", "first"),
                filter=("filter", "first"),
                airmass=("airmass", "median"),
                n_star_used=("source_id", "nunique"),
                resid_med=("resid", "median"),
                resid_mad=("resid", lambda s: float(np.nanmedian(np.abs(s - np.nanmedian(s))))),
            ).reset_index()
            for _, fr in frame_summary.iterrows():
                frame_rows.append({
                    "date": fr["date"],
                    "filter": fr["filter"],
                    "file": fr["file"],
                    "airmass": float(fr["airmass"]) if np.isfinite(fr["airmass"]) else np.nan,
                    "n_star_used": int(fr["n_star_used"]),
                    "resid_med": float(fr["resid_med"]) if np.isfinite(fr["resid_med"]) else np.nan,
                    "resid_mad": float(fr["resid_mad"]) if np.isfinite(fr["resid_mad"]) else np.nan,
                })

            for _, r in sub_good.iterrows():
                point_rows.append({
                    "date": date_val, "filter": filt,
                    "file": str(r["file"]),
                    "source_id": int(r["source_id"]),
                    "airmass": float(r["airmass"]),
                    "mag_inst": float(r["mag_inst"]),
                    "model_mag": float(r["m0_i"] + k1 * r["airmass"]),
                    "delta_m": float(r["delta_m"]),
                    "resid": float(r["resid"]),
                })

            fit_rows.append({
                "date": date_val,
                "filter": filt,
                "k1": k1,
                "k1_err": k1_err,
                "k2": 0.0,
                "zp": 0.0,
                "m0": np.nan,
                "scatter": float(np.nanmedian(list(star_rms.values()))) if star_rms else np.nan,
                "n_total": n_pts,
                "n_used": int(len(sub_good)),
                "outlier_fraction": float(1.0 - (len(sub_good) / max(n_pts, 1))),
                "fit_order": 1,
                "n_stars": int(len(k_list)),
                "n_frames": int(sub_good["file"].nunique()),
                "method": "per_star",
            })

            self._log(f"[{date_val}][{filt}] k1={k1:.5f} ±{k1_err:.5f}, "
                      f"stars={len(k_list)}, frames={int(sub_good['file'].nunique())}")

        if not fit_rows:
            raise RuntimeError("No valid per-star extinction fits produced")

        fit_df = pd.DataFrame(fit_rows)
        points_df = pd.DataFrame(point_rows)
        star_df = pd.DataFrame(star_rows)
        frame_df = pd.DataFrame(frame_rows)

        # Save results
        out_dir = tool_extinction_dir(self.result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        fit_path = out_dir / "per_star_extinction_by_filter.csv"
        pts_path = out_dir / "per_star_extinction_points.csv"
        star_path = out_dir / "extinction_fit_star_stats.csv"
        frame_path = out_dir / "extinction_fit_frame_stats.csv"
        fit_df.to_csv(fit_path, index=False)
        points_df.to_csv(pts_path, index=False)
        if not star_df.empty:
            star_df.to_csv(star_path, index=False)
        if not frame_df.empty:
            frame_df.to_csv(frame_path, index=False)
        self._log(f"Saved {fit_path.name} ({len(fit_df)} rows)")
        self._log(f"Saved {pts_path.name} ({len(points_df)} rows)")
        if not star_df.empty:
            self._log(f"Saved {star_path.name} ({len(star_df)} rows)")
        if not frame_df.empty:
            self._log(f"Saved {frame_path.name} ({len(frame_df)} rows)")

        compat_fit = out_dir / "extinction_fit_by_filter.csv"
        try:
            fit_df.to_csv(compat_fit, index=False)
            self._log(f"Saved {compat_fit.name} ({len(fit_df)} rows)")
        except Exception:
            self._log("Warning: failed to write extinction_fit_by_filter.csv")

        self.finished.emit({
            "fit": fit_df, "points": points_df,
            "mode": "per_star",
        })

    def _run_median_fit(self, phot_df: pd.DataFrame):
        """Median-subtract Δm vs X fit."""
        P = self.params.P
        clip_sigma = self._fit_clip_sigma()
        fit_iters = self._fit_iters()
        min_points = self._fit_min_points()
        snr_cut = self._as_float(getattr(P, "extinction_snr_min", 10.0), 10.0)
        use_weights = bool(getattr(P, "extinction_star_use_weights", True))

        use_dx = bool(getattr(P, "extinction_delta_x_enable", True))
        dx_min_global = self._as_float(getattr(P, "extinction_delta_x_min", 0.3), 0.3)

        var_method_raw = getattr(P, "extinction_varstar_method", "mad")
        var_method = str(var_method_raw if var_method_raw else "mad").strip().lower()
        var_sigma = self._as_float(getattr(P, "extinction_varstar_sigma", 3.0), 3.0)
        var_min_frames = self._as_int(getattr(P, "extinction_varstar_min_frames", 5), 5)

        qc_method_raw = getattr(P, "extinction_frame_qc_method", "mad")
        qc_method = str(qc_method_raw if qc_method_raw else "mad").strip().lower()
        qc_sigma = self._as_float(getattr(P, "extinction_frame_qc_sigma", 3.0), 3.0)

        df = phot_df.copy()
        if "source_id" in df.columns:
            before_sid = len(df)
            df = df[df["source_id"].notna()].copy()
            if before_sid != len(df):
                self._log(f"Dropped {before_sid - len(df)} rows without source_id")
        if "snr" in df.columns:
            n_before = len(df)
            df = df[df["snr"] >= snr_cut]
            self._log(f"SNR filter (>={snr_cut}): {n_before} → {len(df)}")

        nan_am = ~np.isfinite(df["airmass"])
        nan_mag = ~np.isfinite(df["mag_inst"])
        n_drop = int((nan_am | nan_mag).sum())
        if n_drop > 0:
            self._log(f"Dropped {n_drop} points with NaN airmass/mag_inst "
                      f"({int(nan_am.sum())} airmass, {int(nan_mag.sum())} mag)")
        df = df[~nan_am & ~nan_mag].copy()
        if df.empty:
            raise RuntimeError("No data after SNR/airmass filters")

        fit_rows = []
        point_rows = []

        groups = list(df.groupby(["date", "filter"]))
        n_groups = len(groups)
        for i_grp, ((date_val, filt), sub) in enumerate(groups):
            self.progress.emit(i_grp, n_groups, f"Median {date_val}/{filt}")
            sub = sub.copy()
            n_pts = len(sub)

            if use_dx:
                x_all = get_numeric_array(sub, "airmass")
                if np.isfinite(x_all).any():
                    dx_all = float(np.nanmax(x_all) - np.nanmin(x_all))
                    if dx_all < dx_min_global:
                        self._log(f"[{date_val}][{filt}] skipped: ΔX={dx_all:.3f} < {dx_min_global:.3f}")
                        continue

            if n_pts < min_points:
                self._log(f"[{date_val}][{filt}] skipped: only {n_pts} points (min={min_points})")
                continue

            # per-star median within this date+filter
            med_map = sub.groupby("source_id")["mag_inst"].median()
            sub["delta_m"] = sub["mag_inst"] - sub["source_id"].map(med_map)

            # Variable star prefilter (delta_m scatter)
            if var_sigma > 0:
                sub = self._filter_variable_stars(
                    sub, "delta_m", var_method, var_sigma, var_min_frames,
                    label=f"[{date_val}][{filt}]")

            # Frame QC (delta_m median)
            if qc_sigma > 0:
                frame_stats = sub.groupby("file").agg(
                    delta_med=("delta_m", "median"),
                    delta_mad=("delta_m", lambda s: float(np.nanmedian(np.abs(s - np.nanmedian(s))))),
                    n_stars=("delta_m", "count"),
                )
                if qc_method == "std":
                    global_med = float(np.nanmean(frame_stats["delta_med"]))
                    global_sig = float(np.nanstd(frame_stats["delta_med"]))
                else:
                    global_med = float(np.nanmedian(frame_stats["delta_med"]))
                    global_mad = float(np.nanmedian(np.abs(frame_stats["delta_med"] - global_med)))
                    global_sig = MAD_TO_SIGMA * global_mad if global_mad > 0 else float(np.nanstd(frame_stats["delta_med"]))
                if global_sig > 0:
                    bad_mask = np.abs(frame_stats["delta_med"] - global_med) > qc_sigma * global_sig
                    bad_frames = set(frame_stats.index[bad_mask])
                    if bad_frames:
                        self._log(f"[{date_val}][{filt}] Frame QC ({qc_method}, {qc_sigma:.1f}σ): "
                                  f"{len(bad_frames)} frames removed")
                        sub = sub[~sub["file"].isin(bad_frames)]

            if sub.empty:
                continue

            x = get_numeric_array(sub, "airmass")
            y = get_numeric_array(sub, "delta_m")
            w = None
            if use_weights and "snr" in sub.columns:
                snr = get_numeric_array(sub, "snr")
                sig_mag = MAG_ERR_COEFF / np.clip(snr, 1e-6, None)
                w = 1.0 / np.clip(sig_mag ** 2, 1e-12, None)

            k1, zp, scatter, n_used, out_frac = self._robust_weighted_linfit(
                x, y, w=w, clip_sigma=clip_sigma, iters=fit_iters, min_n=min_points
            )
            if not np.isfinite(k1):
                self._log(f"[{date_val}][{filt}] FAILED: fit did not converge")
                continue

            resid = y - (zp + k1 * x)
            for xi, yi, ri in zip(x, y, resid):
                point_rows.append({
                    "date": date_val, "filter": filt,
                    "airmass": float(xi),
                    "delta_m": float(yi),
                    "resid": float(ri),
                })

            fit_rows.append({
                "date": date_val,
                "filter": filt,
                "k1": float(k1),
                "k2": 0.0,
                "zp": float(zp),
                "m0": np.nan,
                "scatter": float(scatter),
                "n_total": n_pts,
                "n_used": int(n_used),
                "outlier_fraction": float(out_frac),
                "fit_order": 1,
                "n_stars": int(sub["source_id"].nunique()),
                "n_frames": int(sub["file"].nunique()),
                "method": "median",
            })

            self._log(f"[{date_val}][{filt}] k1={k1:.5f}, scatter={scatter:.4f}, "
                      f"n={n_used}/{n_pts} ({sub['source_id'].nunique()} stars, {sub['file'].nunique()} frames)")

        if not fit_rows:
            raise RuntimeError("No valid median-subtract extinction fits produced")

        fit_df = pd.DataFrame(fit_rows)
        points_df = pd.DataFrame(point_rows)

        out_dir = tool_extinction_dir(self.result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        fit_path = out_dir / "ensemble_extinction_by_filter.csv"
        pts_path = out_dir / "ensemble_phot_points.csv"
        fit_df.to_csv(fit_path, index=False)
        points_df.to_csv(pts_path, index=False)
        self._log(f"Saved {fit_path.name} ({len(fit_df)} rows)")
        self._log(f"Saved {pts_path.name} ({len(points_df)} rows)")

        compat_fit = out_dir / "extinction_fit_by_filter.csv"
        try:
            fit_df.to_csv(compat_fit, index=False)
            self._log(f"Saved {compat_fit.name} ({len(fit_df)} rows)")
        except Exception:
            self._log("Warning: failed to write extinction_fit_by_filter.csv")

        self.finished.emit({
            "fit": fit_df, "points": points_df,
            "mode": "median",
        })

    def _save_photometry(self, phot_df: pd.DataFrame):
        primary_path, alternate_path = self._extinction_input_paths()
        phot_df.to_csv(primary_path, index=False)
        self._log(f"Saved {primary_path.name} ({len(phot_df)} rows)")
        if alternate_path != primary_path:
            try:
                phot_df.to_csv(alternate_path, index=False)
                self._log(f"Saved {alternate_path.name} ({len(phot_df)} rows)")
            except Exception:
                self._log("Warning: failed to write ensemble_allstar_phot.csv")

    def _load_photometry(self, expected_source_dir: Path | None = None) -> pd.DataFrame:
        primary_path, alternate_path = self._extinction_input_paths()
        for path in (primary_path, alternate_path):
            if path.exists():
                phot_df = pd.read_csv(path)
                if expected_source_dir is not None:
                    expected = self._normalize_path(expected_source_dir)
                    if "source_workspace" in phot_df.columns:
                        cached = {
                            self._normalize_path(v)
                            for v in phot_df["source_workspace"].dropna().astype(str)
                            if str(v).strip()
                        }
                        if cached and expected not in cached:
                            cached_label = sorted(cached)[0]
                            raise FileNotFoundError(
                                f"Cached Step 7 input belongs to a different workspace: {cached_label}"
                            )
                    elif expected != self._normalize_path(self.result_dir):
                        raise FileNotFoundError(
                            "Cached Step 7 input has no source workspace metadata. Load Step 7 again."
                        )
                return phot_df
        raise FileNotFoundError("Step 7 extinction input cache not found. Load Step 7 first.")

    def _extract_date_from_file(self, fname: str, hdr=None) -> str:
        """파일명/폴더명 또는 헤더에서 날짜 추출 (YYYY-MM-DD or YYYYMMDD)"""
        import re
        # 1. 파일명에 날짜 폴더가 포함된 경우 (예: "2024-01-15__pp_image.fits")
        if "__" in fname:
            folder_part = fname.split("__")[0]
            # YYYY-MM-DD 패턴
            m = re.match(r"(\d{4}-\d{2}-\d{2})", folder_part)
            if m:
                return m.group(1)
            # YYYYMMDD 패턴
            m = re.match(r"(\d{8})", folder_part)
            if m:
                d = m.group(1)
                return f"{d[:4]}-{d[4:6]}-{d[6:8]}"
        # 2. DATE-OBS 헤더에서 추출
        if hdr is not None:
            date_obs = str(hdr.get("DATE-OBS", ""))
            if date_obs:
                m = re.match(r"(\d{4}-\d{2}-\d{2})", date_obs)
                if m:
                    return m.group(1)
        return "unknown"

    def _build_frame_airmass(self, idx: pd.DataFrame) -> pd.DataFrame:
        P = self.params.P
        lat = self._as_float(getattr(P, "site_lat_deg", 0.0), 0.0)
        lon = self._as_float(getattr(P, "site_lon_deg", 0.0), 0.0)
        alt = self._as_float(getattr(P, "site_alt_m", 0.0), 0.0)
        tz = self._as_float(getattr(P, "site_tz_offset_hours", 0.0), 0.0)

        rows = []
        for _, r in idx.iterrows():
            fname = str(r.get("file", "")).strip()
            if fname == "":
                continue
            fpath = self._resolve_source_fits_path(fname)
            if fpath is None:
                continue
            try:
                hdr = fits.getheader(fpath)
                info = compute_airmass_from_header(hdr, lat, lon, alt, tz)
                filt = str(r.get("filter", hdr.get("FILTER", ""))).strip().lower()
                obs_date = self._extract_date_from_file(fname, hdr)
                rows.append({
                    "file": fname,
                    "filter": filt,
                    "date": obs_date,
                    **info,
                })
            except Exception:
                continue
        df = pd.DataFrame(rows)
        if df.empty:
            df = pd.DataFrame(columns=["file", "filter", "date", "airmass", "airmass_source", "alt_deg", "zenith_deg", "datetime_utc", "datetime_local", "ra_deg", "dec_deg"])
        else:
            df = df.drop_duplicates(subset=["file"], keep="last").reset_index(drop=True)

        out_dir = tool_extinction_dir(self.result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / "frame_airmass.csv"
        df.to_csv(out_path, index=False)
        self._log(f"Saved {out_path.name} | rows={len(df)}")

        try:
            df.to_csv(self.result_dir / "frame_airmass.csv", index=False)
        except Exception:
            pass
        return df

    def run(self):
        """Entry point: dispatch to ensemble or gaia-absolute mode."""
        try:
            if self.mode == "ensemble":
                if self.task == "photometry":
                    phot_df = self.phot_df if isinstance(self.phot_df, pd.DataFrame) and not self.phot_df.empty else self._run_step7_source_load()
                    self._save_photometry(phot_df)
                    self.finished.emit({
                        "phot": phot_df,
                        "mode": "ensemble",
                        "task": "photometry",
                        "source_result_dir": self._normalize_path(self.source_result_dir),
                    })
                elif self.task == "fit":
                    phot_df = self.phot_df if isinstance(self.phot_df, pd.DataFrame) else None
                    if phot_df is None or phot_df.empty:
                        try:
                            phot_df = self._load_photometry(expected_source_dir=self.source_result_dir)
                        except FileNotFoundError as e:
                            self._log(f"{e}. Rebuilding from Step 7 forced photometry.")
                            phot_df = self._run_step7_source_load()
                            self._save_photometry(phot_df)
                    self._run_ensemble_fit(phot_df)
                else:
                    phot_df = self._run_step7_source_load()
                    self._save_photometry(phot_df)
                    self._run_ensemble_fit(phot_df)
            elif self.mode == "per_star":
                if self.task == "photometry":
                    phot_df = self.phot_df if isinstance(self.phot_df, pd.DataFrame) and not self.phot_df.empty else self._run_step7_source_load()
                    self._save_photometry(phot_df)
                    self.finished.emit({
                        "phot": phot_df,
                        "mode": "per_star",
                        "task": "photometry",
                        "source_result_dir": self._normalize_path(self.source_result_dir),
                    })
                elif self.task == "fit":
                    phot_df = self.phot_df if isinstance(self.phot_df, pd.DataFrame) else None
                    if phot_df is None or phot_df.empty:
                        try:
                            phot_df = self._load_photometry(expected_source_dir=self.source_result_dir)
                        except FileNotFoundError as e:
                            self._log(f"{e}. Rebuilding from Step 7 forced photometry.")
                            phot_df = self._run_step7_source_load()
                            self._save_photometry(phot_df)
                    self._run_per_star_fit(phot_df)
                else:
                    phot_df = self._run_step7_source_load()
                    self._save_photometry(phot_df)
                    self._run_per_star_fit(phot_df)
            elif self.mode == "median":
                if self.task == "photometry":
                    phot_df = self.phot_df if isinstance(self.phot_df, pd.DataFrame) and not self.phot_df.empty else self._run_step7_source_load()
                    self._save_photometry(phot_df)
                    self.finished.emit({
                        "phot": phot_df,
                        "mode": "median",
                        "task": "photometry",
                        "source_result_dir": self._normalize_path(self.source_result_dir),
                    })
                elif self.task == "fit":
                    phot_df = self.phot_df if isinstance(self.phot_df, pd.DataFrame) else None
                    if phot_df is None or phot_df.empty:
                        try:
                            phot_df = self._load_photometry(expected_source_dir=self.source_result_dir)
                        except FileNotFoundError as e:
                            self._log(f"{e}. Rebuilding from Step 7 forced photometry.")
                            phot_df = self._run_step7_source_load()
                            self._save_photometry(phot_df)
                    self._run_median_fit(phot_df)
                else:
                    phot_df = self._run_step7_source_load()
                    self._save_photometry(phot_df)
                    self._run_median_fit(phot_df)
            else:
                self._run_gaia_absolute()
        except Exception as e:
            self.error.emit(str(e))

    def _run_gaia_absolute(self):
        """Original Gaia-absolute extinction fit (requires step 9 data)."""
        try:
            P = self.params.P
            result_dir = self.result_dir

            idx_path = step7_forced_phot_dir(result_dir) / "photometry_index.csv"
            if not idx_path.exists():
                raise FileNotFoundError("photometry index csv not found")
            idx = pd.read_csv(idx_path)
            if "file" not in idx.columns:
                for cand in ("fname", "frame", "image", "fits", "name"):
                    if cand in idx.columns:
                        idx = idx.rename(columns={cand: "file"})
                        break
            if "filter" in idx.columns:
                idx["filter"] = idx["filter"].astype(str).str.strip().str.lower()
            elif "FILTER" in idx.columns:
                idx["filter"] = idx["FILTER"].astype(str).str.strip().str.lower()
            else:
                idx["filter"] = ""

            qc_exclude: set[str] = set()
            for qp in [step7_forced_phot_dir(result_dir) / "frame_quality.csv", result_dir / "frame_quality.csv"]:
                if not qp.exists():
                    continue
                try:
                    dfq = pd.read_csv(qp)
                except Exception:
                    continue
                if {"file", "passed"} <= set(dfq.columns):
                    bad = dfq.loc[dfq["passed"] == False, "file"].astype(str)
                    qc_exclude |= set(bad.tolist())
                break
            if qc_exclude:
                before = len(idx)
                idx = idx[~idx["file"].astype(str).isin(qc_exclude)].reset_index(drop=True)
                self._log(f"Frame QC: {before - len(idx)} frame(s) excluded by Step 7 quality flags")

            # Apply step10 manual frame exclusions
            step10_excl = _load_frame_excludes(result_dir)
            if step10_excl and "file" in idx.columns:
                before = len(idx)
                idx = idx[~idx["file"].astype(str).isin(set(step10_excl.keys()))].reset_index(drop=True)
                self._log(f"Frame QC: {before - len(idx)} frame(s) excluded by step10 manual exclusion")

            rows = []
            total = len(idx)
            for i, r in idx.iterrows():
                if self._stop_requested:
                    self.finished.emit({"stopped": True})
                    return
                fname = str(r.get("file", "") or "").strip()
                if not fname:
                    continue
                filt_hint = str(r.get("filter", "") or "").strip().lower()
                dfp = load_frame_photometry(result_dir, fname, filt_hint)
                if dfp is None or dfp.empty:
                    continue
                dfp = dfp.copy()
                mag_col = None
                for cand in ("mag_inst", "mag", "mag_ap", "mag_apcorr"):
                    if cand in dfp.columns:
                        mag_col = cand
                        break
                if mag_col is None:
                    continue
                err_col = None
                for cand in ("mag_err", "emag", "emag_inst", "magerr"):
                    if cand in dfp.columns:
                        err_col = cand
                        break
                if err_col is None:
                    dfp["mag_err"] = np.nan
                    err_col = "mag_err"
                snr_col = "snr" if "snr" in dfp.columns else None
                if "FILTER" not in dfp.columns:
                    if "filter" in dfp.columns:
                        dfp["FILTER"] = dfp["filter"]
                    else:
                        dfp["FILTER"] = filt_hint
                keep_cols = ["ID", "FILTER", mag_col, err_col]
                if "source_id" in dfp.columns:
                    keep_cols.append("source_id")
                if snr_col:
                    keep_cols.append(snr_col)
                if "ID" not in dfp.columns:
                    continue
                tmp = dfp[keep_cols].copy()
                tmp = tmp.rename(columns={mag_col: "mag_inst", err_col: "mag_err"})
                if snr_col is None:
                    tmp["snr"] = np.nan
                else:
                    tmp = tmp.rename(columns={snr_col: "snr"})
                tmp["file"] = fname
                rows.append(tmp)
                self.progress.emit(i + 1, total, fname)

            if not rows:
                raise RuntimeError("No photometry data found")
            all_df = pd.concat(rows, ignore_index=True)
            all_df["FILTER"] = all_df["FILTER"].astype(str).str.strip().str.lower()

            grp = (
                all_df.groupby(["ID", "FILTER"])
                .agg(
                    mag_inst_med=("mag_inst", lambda s: self._robust_median_and_err(s)[0]),
                    n_frames=("mag_inst", lambda s: self._robust_median_and_err(s)[2]),
                    snr_med=("snr", lambda s: float(np.nanmedian(np.asarray(s, float))) if np.isfinite(np.asarray(s, float)).any() else np.nan),
                )
                .reset_index()
            )

            wide_mag = grp.pivot_table(index="ID", columns="FILTER", values="mag_inst_med", aggfunc="median")
            wide_snr = grp.pivot_table(index="ID", columns="FILTER", values="snr_med", aggfunc="median")
            wide_mag.columns = [f"mag_inst_{c}" for c in wide_mag.columns]
            wide_snr.columns = [f"snr_{c}" for c in wide_snr.columns]
            wide = pd.concat([wide_mag, wide_snr], axis=1).reset_index()

            refbuild = step6_refbuild_dir(result_dir)
            per_filter = sorted(refbuild.glob("ref_catalog_*.tsv")) if refbuild.exists() else []
            master_path = next(
                (p for p in [refbuild / "ref_catalog.tsv"] + per_filter
                 + [refbuild / "master_catalog.tsv"] if p.exists()),
                None,
            )
            if master_path is None:
                raise FileNotFoundError("ref_catalog.tsv not found in step6_refbuild/")
            master = pd.read_csv(master_path, sep="\t")
            if "ID" not in master.columns:
                raise RuntimeError("ref_catalog.tsv missing ID column")

            merge_cols = ["ID"]
            if "source_id" in master.columns:
                merge_cols.append("source_id")
            for col in ("gaia_G", "gaia_BP", "gaia_RP", "gmag", "bpmag", "rpmag", "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag"):
                if col in master.columns and col not in merge_cols:
                    merge_cols.append(col)
            df = wide.merge(master[merge_cols], on="ID", how="left")

            g_col = None
            bp_col = None
            rp_col = None
            for cand in ("gaia_G", "gmag", "phot_g_mean_mag"):
                if cand in df.columns:
                    g_col = cand
                    break
            for cand in ("gaia_BP", "bpmag", "phot_bp_mean_mag"):
                if cand in df.columns:
                    bp_col = cand
                    break
            for cand in ("gaia_RP", "rpmag", "phot_rp_mean_mag"):
                if cand in df.columns:
                    rp_col = cand
                    break

            if g_col is None or bp_col is None or rp_col is None:
                if "source_id" not in df.columns:
                    raise RuntimeError("Gaia mags not available and source_id missing")
                src_sid = coerce_int64_source_id(df["source_id"])
                df = df.loc[src_sid.notna()].copy()
                df["source_id"] = src_sid[src_sid.notna()].astype("int64")
                gaia_path = step5_wcs_dir(result_dir) / "gaia_fov.ecsv"
                if not gaia_path.exists():
                    raise RuntimeError("gaia_fov.ecsv not found")
                t_gaia = Table.read(gaia_path, format="ascii.ecsv")
                gaia_df = t_gaia.to_pandas()
                gaia_sid = coerce_int64_source_id(gaia_df["source_id"])
                gaia_df = gaia_df.loc[gaia_sid.notna()].copy()
                gaia_df["source_id"] = gaia_sid[gaia_sid.notna()].astype("int64")
                gaia_cols = ["source_id", "phot_g_mean_mag"]
                if "phot_bp_mean_mag" in gaia_df.columns:
                    gaia_cols.append("phot_bp_mean_mag")
                if "phot_rp_mean_mag" in gaia_df.columns:
                    gaia_cols.append("phot_rp_mean_mag")
                df = df.merge(gaia_df[gaia_cols], on="source_id", how="left")
                g_col = "phot_g_mean_mag"
                bp_col = "phot_bp_mean_mag"
                rp_col = "phot_rp_mean_mag"

            df["gaia_G"] = pd.to_numeric(df[g_col], errors="coerce").astype(float)
            df["gaia_BP"] = pd.to_numeric(df[bp_col], errors="coerce").astype(float)
            df["gaia_RP"] = pd.to_numeric(df[rp_col], errors="coerce").astype(float)
            dfm = df[np.isfinite(df["gaia_G"]) & np.isfinite(df["gaia_BP"]) & np.isfinite(df["gaia_RP"])].copy()
            dfm["gaia_BP_RP"] = dfm["gaia_BP"] - dfm["gaia_RP"]

            out_cal = dfm.copy()
            xcol = out_cal["gaia_BP_RP"].to_numpy(float)
            G = out_cal["gaia_G"].to_numpy(float)

            G_minus = {}
            for band, info in GAIA_TO_SDSS.items():
                lo, hi = info["color_range"]
                m = np.isfinite(xcol) & (xcol >= lo) & (xcol <= hi)
                arr = np.full_like(G, np.nan)
                arr[m] = self._poly_eval(xcol[m], info["coeffs"])
                G_minus[band] = arr

            out_cal["sdss_g_ref"] = G - G_minus["g"]
            out_cal["sdss_r_ref"] = G - G_minus["r"]
            out_cal["sdss_i_ref"] = G - G_minus["i"]

            g_inst = out_cal.get("mag_inst_g", pd.Series(np.full(len(out_cal), np.nan))).to_numpy(float)
            r_inst = out_cal.get("mag_inst_r", pd.Series(np.full(len(out_cal), np.nan))).to_numpy(float)
            i_inst = out_cal.get("mag_inst_i", pd.Series(np.full(len(out_cal), np.nan))).to_numpy(float)

            color_gr = g_inst - r_inst
            color_ri = r_inst - i_inst

            clip_sigma = self._fit_clip_sigma()
            fit_iters = self._fit_iters()
            min_match = self._fit_min_points()

            delta_g = out_cal["sdss_g_ref"].to_numpy(float) - g_inst
            mg = np.isfinite(delta_g) & np.isfinite(color_gr) & np.isfinite(g_inst)
            ct_g = self._robust_linfit(color_gr[mg], delta_g[mg], clip_sigma=clip_sigma, iters=fit_iters, min_n=min_match)[0]
            if not np.isfinite(ct_g):
                ct_g = 0.0

            delta_r = out_cal["sdss_r_ref"].to_numpy(float) - r_inst
            mr = np.isfinite(delta_r) & np.isfinite(color_gr) & np.isfinite(r_inst)
            ct_r = self._robust_linfit(color_gr[mr], delta_r[mr], clip_sigma=clip_sigma, iters=fit_iters, min_n=min_match)[0]
            if not np.isfinite(ct_r):
                ct_r = 0.0

            delta_i = out_cal["sdss_i_ref"].to_numpy(float) - i_inst
            mi = np.isfinite(delta_i) & np.isfinite(color_ri) & np.isfinite(i_inst)
            ct_i = self._robust_linfit(color_ri[mi], delta_i[mi], clip_sigma=clip_sigma, iters=fit_iters, min_n=min_match)[0]
            if not np.isfinite(ct_i):
                ct_i = 0.0

            out_cal["color_gr"] = color_gr
            out_cal["color_ri"] = color_ri

            frame_airmass = self._build_frame_airmass(idx)
            if frame_airmass is None or frame_airmass.empty:
                raise RuntimeError("frame_airmass.csv missing and airmass computation failed")

            cal_cols = ["ID", "sdss_g_ref", "sdss_r_ref", "sdss_i_ref", "color_gr", "color_ri", "gaia_BP_RP"]
            obs = all_df.merge(out_cal[cal_cols], on="ID", how="left")
            # date 컬럼도 함께 merge
            merge_cols = ["file", "filter", "airmass"]
            if "date" in frame_airmass.columns:
                merge_cols.append("date")
            obs = obs.merge(frame_airmass[merge_cols], left_on=["file", "FILTER"], right_on=["file", "filter"], how="left")
            if "date" not in obs.columns:
                obs["date"] = "unknown"

            obs["ref_mag"] = np.nan
            obs.loc[obs["FILTER"] == "g", "ref_mag"] = obs.loc[obs["FILTER"] == "g", "sdss_g_ref"]
            obs.loc[obs["FILTER"] == "r", "ref_mag"] = obs.loc[obs["FILTER"] == "r", "sdss_r_ref"]
            obs.loc[obs["FILTER"] == "i", "ref_mag"] = obs.loc[obs["FILTER"] == "i", "sdss_i_ref"]

            obs["color_term"] = np.nan
            obs.loc[obs["FILTER"] == "g", "color_term"] = ct_g * obs.loc[obs["FILTER"] == "g", "color_gr"]
            obs.loc[obs["FILTER"] == "r", "color_term"] = ct_r * obs.loc[obs["FILTER"] == "r", "color_gr"]
            obs.loc[obs["FILTER"] == "i", "color_term"] = ct_i * obs.loc[obs["FILTER"] == "i", "color_ri"]

            obs["delta"] = obs["ref_mag"] - (obs["mag_inst"] + obs["color_term"])
            obs["cal_ok"] = False
            bp = obs["gaia_BP_RP"].to_numpy(float)
            obs.loc[(obs["FILTER"] == "g") & np.isfinite(bp) & (bp >= 0.3) & (bp <= 3.0), "cal_ok"] = True
            obs.loc[(obs["FILTER"] == "r") & np.isfinite(bp) & (bp >= 0.0) & (bp <= 3.0), "cal_ok"] = True
            obs.loc[(obs["FILTER"] == "i") & np.isfinite(bp) & (bp >= 0.5) & (bp <= 2.0), "cal_ok"] = True

            snr_cut = self._as_float(getattr(P, "cmd_snr_calib_min", 20.0), 20.0)
            obs["snr_ok"] = True
            if "snr" in obs.columns:
                svals = obs["snr"].to_numpy(float)
                obs["snr_ok"] = np.isfinite(svals) & (svals >= snr_cut)

            filters_seen = sorted(set(all_df["FILTER"].dropna().astype(str)))
            out_dir = tool_extinction_dir(result_dir)
            out_dir.mkdir(parents=True, exist_ok=True)
            step11_out = out_dir
            stats_rows = []
            for filt in filters_seen:
                fmask = obs["FILTER"] == filt
                n_total = int(fmask.sum())
                n_ref = int(np.isfinite(obs.loc[fmask, "ref_mag"]).sum())
                n_air = int(np.isfinite(obs.loc[fmask, "airmass"]).sum())
                n_color = int(np.isfinite(obs.loc[fmask, "color_term"]).sum())
                n_snr = int(obs.loc[fmask, "snr_ok"].sum()) if "snr_ok" in obs.columns else n_total
                n_cal = int(obs.loc[fmask, "cal_ok"].sum())
                n_delta = int(np.isfinite(obs.loc[fmask, "delta"]).sum())
                self._log(
                    f"Filter[{filt}] total={n_total}, ref={n_ref}, color={n_color}, "
                    f"airmass={n_air}, snr_ok={n_snr}, cal_ok={n_cal}, delta_ok={n_delta}"
                )
                stats_rows.append({
                    "filter": filt,
                    "n_total": n_total,
                    "n_ref": n_ref,
                    "n_color": n_color,
                    "n_airmass": n_air,
                    "n_snr_ok": n_snr,
                    "n_cal_ok": n_cal,
                    "n_delta_ok": n_delta,
                })

            if stats_rows:
                stats_df = pd.DataFrame(stats_rows)
                stats_path = out_dir / "extinction_fit_filter_stats.csv"
                stats_df.to_csv(stats_path, index=False)
                self._log(f"Saved {stats_path.name} | rows={len(stats_df)}")

            obs = obs[np.isfinite(obs["delta"]) & np.isfinite(obs["airmass"]) & obs["cal_ok"] & obs["snr_ok"]].copy()

            # 날짜별 분류 통계
            dates_seen = sorted(set(obs["date"].dropna().astype(str)))
            self._log(f"Dates found: {dates_seen}")
            if len(obs):
                for (date_val, filt), sub in obs.groupby(["date", "FILTER"]):
                    self._log(f"Fit candidates [{date_val}][{filt}]: {len(sub)} points")

            fit_rows = []
            point_rows = []
            use_quadratic = bool(getattr(P, "extinction_use_quadratic", True))
            use_color_extinction = bool(getattr(P, "extinction_use_color_dependent", True))
            min_quad = max(min_match, int(getattr(P, "extinction_min_points_quadratic", 20)))
            min_color = max(min_match, int(getattr(P, "extinction_min_points_color", 30)))

            # 날짜별 + 필터별 피팅
            for (date_val, filt), sub in obs.groupby(["date", "FILTER"]):
                x = sub["airmass"].to_numpy(float)
                y = sub["delta"].to_numpy(float)

                # 색지수 선택: g/r 필터면 g-r, i 필터면 r-i
                if filt in ("g", "r"):
                    color_col = "color_gr"
                elif filt == "i":
                    color_col = "color_ri"
                else:
                    color_col = "color_gr"  # 기본값

                color = sub[color_col].to_numpy(float) if color_col in sub.columns else np.full(len(x), np.nan)

                # 1차 피팅 (기본)
                k1_lin, zp_lin, scatter_lin, n_ref_lin, out_frac_lin = self._robust_linfit(
                    x, y, clip_sigma=clip_sigma, iters=fit_iters, min_n=min_match
                )

                # 색 의존 소광 피팅: delta = k' * X + k'' * C * X
                k_prime, k_color, zp_col, scatter_col, n_ref_col, out_frac_col, k_col_sig = (
                    np.nan, np.nan, np.nan, np.nan, 0, np.nan, False
                )
                n_valid_color = np.sum(np.isfinite(x) & np.isfinite(color) & np.isfinite(y))
                if use_color_extinction and n_valid_color >= min_color:
                    k_prime, k_color, zp_col, scatter_col, n_ref_col, out_frac_col, k_col_sig = self._robust_color_extinction_fit(
                        x, color, y, clip_sigma=clip_sigma, iters=fit_iters, min_n=min_color
                    )
                    if k_col_sig and np.isfinite(k_color):
                        self._log(f"[{date_val}][{filt}] color extinction fit: k'={k_prime:.4f}, k''={k_color:.4f} (n={n_ref_col})")

                # 2차 피팅 시도 (에어매스 비선형)
                k1_q, k2_q, zp_q, scatter_q, n_ref_q, out_frac_q, k2_sig = (
                    np.nan, np.nan, np.nan, np.nan, 0, np.nan, False
                )
                if use_quadratic and len(x) >= min_quad:
                    k1_q, k2_q, zp_q, scatter_q, n_ref_q, out_frac_q, k2_sig = self._robust_quadfit(
                        x, y, clip_sigma=clip_sigma, iters=fit_iters, min_n=min_quad
                    )

                # 최선의 피팅 선택
                # 우선순위: 색 의존 피팅 > 2차 피팅 > 1차 피팅
                k1, k2, k_color_final, zp = k1_lin, 0.0, 0.0, zp_lin
                scatter, n_ref, out_frac = scatter_lin, n_ref_lin, out_frac_lin
                fit_order = 1

                # 색 의존 피팅이 유의미하면 사용
                if k_col_sig and np.isfinite(k_color) and scatter_col < scatter_lin * 0.90:
                    k1, k2, k_color_final, zp = k_prime, 0.0, k_color, zp_col
                    scatter, n_ref, out_frac = scatter_col, n_ref_col, out_frac_col
                    fit_order = 3  # color-dependent
                    yhat = zp + k1 * x + k_color_final * color * x
                    self._log(f"[{date_val}][{filt}] using COLOR-DEPENDENT fit: k'={k1:.4f}, k''={k_color_final:.4f}")
                # 2차 피팅이 유의미하면 사용
                elif k2_sig and np.isfinite(k2_q) and scatter_q < scatter_lin * 0.95:
                    k1, k2, zp = k1_q, k2_q, zp_q
                    scatter, n_ref, out_frac = scatter_q, n_ref_q, out_frac_q
                    fit_order = 2
                    yhat = zp + k1 * x + k2 * x * x
                    self._log(f"[{date_val}][{filt}] using quadratic fit (k2={k2:.4f})")
                else:
                    yhat = zp + k1 * x

                if not np.isfinite(k1):
                    self._log(f"[{date_val}][{filt}] FAILED: insufficient points (n={len(x)})")
                    continue

                resid = y - yhat
                fit_rows.append({
                    "date": date_val,
                    "filter": filt,
                    "k1": k1,
                    "k2": k2,
                    "k_color": k_color_final,
                    "zp": zp,
                    "scatter": scatter,
                    "n_ref": n_ref,
                    "outlier_fraction": out_frac,
                    "fit_order": fit_order,
                    "method": "gaia",
                })
                for xi, yi, ri, ci in zip(x, y, resid, color):
                    point_rows.append({"date": date_val, "filter": filt, "airmass": xi, "delta": yi, "resid": ri, "color": ci})

            fit_df = pd.DataFrame(fit_rows)
            if fit_df.empty:
                raise RuntimeError("No valid extinction fits produced")

            # Step 11 전용 폴더에 저장
            fit_path = step11_out / "step11_extinction_fit_by_filter.csv"
            fit_df.to_csv(fit_path, index=False)
            self._log(f"Saved {fit_path.name} | rows={len(fit_df)}")

            # 기존 경로에도 저장 (호환성)
            fit_path_old = out_dir / "extinction_fit_by_filter.csv"
            fit_df.to_csv(fit_path_old, index=False)

            for _, row in fit_df.iterrows():
                k2_str = f", k2={row['k2']:.4f}" if row.get('k2', 0) != 0 else ""
                k_color_str = f", k''={row['k_color']:.4f}" if row.get('k_color', 0) != 0 else ""
                fit_type = {1: "linear", 2: "quadratic", 3: "color-dep"}.get(row.get('fit_order', 1), "linear")
                date_str = row.get('date', 'unknown')
                self._log(
                    f"[{date_str}][{row['filter']}] k'={row['k1']:.4f}{k2_str}{k_color_str}, zp={row['zp']:.4f}, "
                    f"scatter={row['scatter']:.4f}, n={int(row['n_ref'])}, type={fit_type}"
                )

            points_df = pd.DataFrame(point_rows)
            points_path = step11_out / "step11_extinction_fit_points.csv"
            points_df.to_csv(points_path, index=False)
            self._log(f"Saved {points_path.name} | rows={len(points_df)}")

            self.finished.emit({"fit": fit_df, "points": points_df, "mode": "gaia"})
        except Exception as e:
            self.error.emit(str(e))


class ExtinctionFitWindow(QWidget):
    """Per-star Bouguer extinction tool with image-based star selection."""

    _DATE_KEY_RE = re.compile(r"(20\d{6})")

    def __init__(self, params, data_dir: Path, result_dir: Path, parent=None):
        super().__init__(parent)
        self.params = params
        self.data_dir = Path(data_dir)
        self.result_dir = Path(result_dir)
        self.worker = None
        self.loaded_source_result_dir: Path | None = None
        self._current_mode = "per_star"

        self.points_df: pd.DataFrame | None = None
        self.fit_df: pd.DataFrame | None = None
        self.phot_df: pd.DataFrame | None = None
        self.selection_stats_df = pd.DataFrame()
        self.selection_state: dict[str, dict[str, set[int]]] = {}
        self._selection_meta_df = pd.DataFrame()
        self._step8_target_by_filter: dict[str, int] = {}
        self._source_file_path_map: dict[str, str] | None = None
        self._source_file_path_map_dir: str | None = None

        self.selection_frame_df = pd.DataFrame()
        self.selection_image_data = None
        self.selection_header = None
        self.selection_selected_source_id: int | None = None
        self._selection_sid_to_row: dict[int, int] = {}
        self._selection_imshow = None
        self._selection_scat_candidate = None
        self._selection_scat_selected = None
        self._selection_scat_rejected = None
        self._selection_scat_highlight = None
        self._selection_panning = False
        self._selection_pan_start = None

        self.setWindowTitle("Atmospheric Extinction Fit Tool")
        self.resize(1220, 860)

        layout = QVBoxLayout(self)

        info = QLabel(
            "Atmospheric extinction coefficient tool.\n"
            "Load Step 7 forced photometry, choose stable stars per date/filter on the image, "
            "and run per-star Bouguer fits to measure k1."
        )
        info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 8px; border-radius: 5px; }")
        layout.addWidget(info)

        settings_group = QGroupBox("Tool")
        settings_layout = QHBoxLayout(settings_group)
        mode_label = QLabel("Mode: Per-star Bouguer only")
        mode_label.setStyleSheet("QLabel { font-weight: bold; }")
        settings_layout.addWidget(mode_label)
        settings_layout.addStretch()
        btn_params = QPushButton("Parameters")
        btn_params.setStyleSheet(
            "QPushButton { background-color: #9C27B0; color: white; "
            "font-weight: bold; padding: 6px 12px; }"
        )
        btn_params.clicked.connect(self.open_parameters_dialog)
        settings_layout.addWidget(btn_params)
        layout.addWidget(settings_group)

        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)

        self.tabs.addTab(self._build_source_tab(), "Step 7 Source")
        self.selection_tab = self._build_selection_tab()
        self.tabs.addTab(self.selection_tab, "Selection")
        self.tabs.addTab(self._build_fit_tab(), "Extinction Fit")

        log_group = QGroupBox("Log")
        log_layout = QVBoxLayout(log_group)
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setStyleSheet("QTextEdit { font-family: monospace; font-size: 9pt; }")
        log_layout.addWidget(self.log_text)
        layout.addWidget(log_group)

    def _build_source_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)

        info = QLabel(
            "Load and cache Step 7 forced photometry + frame airmass from the selected workspace.\n"
            "The cached table is reused by the selection tab and the Bouguer fit."
        )
        info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 8px; border-radius: 5px; }")
        layout.addWidget(info)

        source_layout = QHBoxLayout()
        source_layout.addWidget(QLabel("Source workspace:"))
        self.source_workspace_edit = QLineEdit(str(self.result_dir))
        self.source_workspace_edit.setPlaceholderText("result folder containing step7_forced_phot")
        self.source_workspace_edit.textChanged.connect(self._on_source_workspace_changed)
        source_layout.addWidget(self.source_workspace_edit, 1)

        btn_browse_source = QPushButton("Browse...")
        btn_browse_source.clicked.connect(self._browse_source_workspace)
        source_layout.addWidget(btn_browse_source)

        btn_use_current = QPushButton("Current")
        btn_use_current.clicked.connect(self._use_current_source_workspace)
        source_layout.addWidget(btn_use_current)
        layout.addLayout(source_layout)

        phot_controls = QHBoxLayout()
        self.btn_run_phot = QPushButton("Load Step 7")
        self.btn_run_phot.setStyleSheet(
            "QPushButton { background-color: #4CAF50; color: white; "
            "font-weight: bold; padding: 6px 12px; }"
        )
        self.btn_run_phot.clicked.connect(self.run_photometry)
        phot_controls.addWidget(self.btn_run_phot)

        self.btn_stop = QPushButton("Stop")
        self.btn_stop.setStyleSheet(
            "QPushButton { background-color: #f44336; color: white; "
            "font-weight: bold; padding: 6px 12px; }"
        )
        self.btn_stop.clicked.connect(self.stop_fit)
        self.btn_stop.setEnabled(False)
        phot_controls.addWidget(self.btn_stop)
        phot_controls.addStretch()
        layout.addLayout(phot_controls)

        progress_layout = QHBoxLayout()
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        progress_layout.addWidget(self.progress_bar)
        self.progress_label = QLabel("Ready")
        self.progress_label.setMinimumWidth(320)
        progress_layout.addWidget(self.progress_label)
        layout.addLayout(progress_layout)

        return tab

    def _build_selection_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)

        info = QLabel(
            "Choose extinction stars on a representative frame.\n"
            "State is saved per date/filter as use, reject, or candidate, and those choices feed the fit."
        )
        info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 8px; border-radius: 5px; }")
        layout.addWidget(info)

        group_layout = QHBoxLayout()
        group_layout.addWidget(QLabel("Date:"))
        self.sel_date_combo = QComboBox()
        self.sel_date_combo.currentIndexChanged.connect(self._on_selection_group_changed)
        group_layout.addWidget(self.sel_date_combo)

        group_layout.addWidget(QLabel("Filter:"))
        self.sel_filter_combo = QComboBox()
        self.sel_filter_combo.currentIndexChanged.connect(self._on_selection_group_changed)
        group_layout.addWidget(self.sel_filter_combo)

        group_layout.addWidget(QLabel("Frame:"))
        self.sel_frame_combo = QComboBox()
        self.sel_frame_combo.currentIndexChanged.connect(self._on_selection_frame_changed)
        group_layout.addWidget(self.sel_frame_combo)

        group_layout.addStretch()

        self.chk_exclude_gaia_var = QCheckBox("Exclude Gaia VARIABLE")
        self.chk_exclude_gaia_var.setChecked(True)
        group_layout.addWidget(self.chk_exclude_gaia_var)

        self.chk_exclude_step8_target = QCheckBox("Exclude Step 8 Target")
        self.chk_exclude_step8_target.setChecked(True)
        group_layout.addWidget(self.chk_exclude_step8_target)

        self.chk_show_selection_ids = QCheckBox("Show IDs")
        self.chk_show_selection_ids.setChecked(True)
        self.chk_show_selection_ids.stateChanged.connect(self._update_selection_overlay)
        group_layout.addWidget(self.chk_show_selection_ids)

        layout.addLayout(group_layout)

        action_layout = QHBoxLayout()
        self.btn_sel_auto = QPushButton("Auto Pick")
        self.btn_sel_auto.clicked.connect(self._auto_pick_current_group)
        action_layout.addWidget(self.btn_sel_auto)

        self.btn_sel_use = QPushButton("Use")
        self.btn_sel_use.clicked.connect(lambda: self._apply_state_to_selected("selected"))
        action_layout.addWidget(self.btn_sel_use)

        self.btn_sel_reject = QPushButton("Reject")
        self.btn_sel_reject.clicked.connect(lambda: self._apply_state_to_selected("rejected"))
        action_layout.addWidget(self.btn_sel_reject)

        self.btn_sel_reset = QPushButton("Reset")
        self.btn_sel_reset.clicked.connect(lambda: self._apply_state_to_selected("candidate"))
        action_layout.addWidget(self.btn_sel_reset)

        self.btn_sel_clear = QPushButton("Clear Group")
        self.btn_sel_clear.clicked.connect(self._clear_current_group_states)
        action_layout.addWidget(self.btn_sel_clear)

        self.btn_sel_copy_date = QPushButton("Copy Same Date")
        self.btn_sel_copy_date.clicked.connect(self._copy_current_selection_to_same_date)
        action_layout.addWidget(self.btn_sel_copy_date)

        self.btn_sel_copy_all = QPushButton("Copy All Groups")
        self.btn_sel_copy_all.clicked.connect(self._copy_current_selection_to_all_groups)
        action_layout.addWidget(self.btn_sel_copy_all)

        action_layout.addStretch()
        layout.addLayout(action_layout)

        self.selection_summary_label = QLabel("No Step 7 source table loaded.")
        self.selection_summary_label.setStyleSheet("QLabel { color: #333; font-weight: bold; }")
        layout.addWidget(self.selection_summary_label)

        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter, 1)

        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        left_layout.setContentsMargins(0, 0, 0, 0)

        img_group = QGroupBox("Representative Frame")
        img_layout = QVBoxLayout(img_group)
        self.selection_selected_label = QLabel("Selected: (none)")
        img_layout.addWidget(self.selection_selected_label)

        self.selection_figure = Figure(figsize=(6, 5))
        self.selection_canvas = FigureCanvas(self.selection_figure)
        self.selection_ax = self.selection_figure.add_subplot(111)
        img_layout.addWidget(NavigationToolbar(self.selection_canvas, self))
        img_layout.addWidget(self.selection_canvas)
        self.selection_canvas.mpl_connect("scroll_event", self._on_selection_scroll)
        self.selection_canvas.mpl_connect("button_press_event", self._on_selection_button_press)
        self.selection_canvas.mpl_connect("button_release_event", self._on_selection_button_release)
        self.selection_canvas.mpl_connect("motion_notify_event", self._on_selection_motion)
        self.selection_canvas.mpl_connect("button_press_event", self._on_selection_click)
        left_layout.addWidget(img_group)

        splitter.addWidget(left_widget)

        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        right_layout.setContentsMargins(0, 0, 0, 0)

        table_group = QGroupBox("Extinction Star Table")
        table_layout = QVBoxLayout(table_group)
        self.selection_table = QTableWidget()
        self.selection_table.setColumnCount(11)
        self.selection_table.setHorizontalHeaderLabels(
            ["State", "ID", "source_id", "G", "BP-RP", "Gaia Var", "Nfr", "ΔX", "SNRmed", "RMS", "k_i"]
        )
        self.selection_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.selection_table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.selection_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.selection_table.verticalHeader().setVisible(False)
        header = self.selection_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(4, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(5, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(6, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(7, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(8, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(9, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(10, QHeaderView.Stretch)
        self.selection_table.itemSelectionChanged.connect(self._on_selection_table_changed)
        table_layout.addWidget(self.selection_table)
        right_layout.addWidget(table_group, 3)

        preview_group = QGroupBox("Selected Star Preview")
        preview_layout = QVBoxLayout(preview_group)
        self.selection_preview_figure = Figure(figsize=(4, 2.4))
        self.selection_preview_canvas = FigureCanvas(self.selection_preview_figure)
        preview_layout.addWidget(self.selection_preview_canvas)
        right_layout.addWidget(preview_group, 2)

        splitter.addWidget(right_widget)
        splitter.setSizes([760, 460])

        return tab

    def _build_fit_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)

        info = QLabel(
            "Run per-star Bouguer extinction fitting on the cached Step 7 source table.\n"
            "If Selection tab states exist, rejected stars are dropped and selected stars are preferred."
        )
        info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 8px; border-radius: 5px; }")
        layout.addWidget(info)

        controls = QHBoxLayout()
        self.btn_run_fit = QPushButton("Run Fit")
        self.btn_run_fit.setStyleSheet(
            "QPushButton { background-color: #4CAF50; color: white; "
            "font-weight: bold; padding: 6px 12px; }"
        )
        self.btn_run_fit.clicked.connect(self.run_fit)
        controls.addWidget(self.btn_run_fit)

        self.btn_save = QPushButton("Save Plots")
        self.btn_save.setStyleSheet(
            "QPushButton { background-color: #607D8B; color: white; "
            "font-weight: bold; padding: 6px 12px; }"
        )
        self.btn_save.clicked.connect(self.save_plots)
        self.btn_save.setEnabled(False)
        controls.addWidget(self.btn_save)
        controls.addStretch()
        layout.addLayout(controls)

        filter_layout = QHBoxLayout()
        filter_layout.addWidget(QLabel("Date:"))
        self.date_combo = QComboBox()
        self.date_combo.addItem("All Dates")
        self.date_combo.currentIndexChanged.connect(self._on_filter_changed)
        filter_layout.addWidget(self.date_combo)

        filter_layout.addWidget(QLabel("Filter:"))
        self.filter_combo = QComboBox()
        self.filter_combo.addItem("All Filters")
        self.filter_combo.currentIndexChanged.connect(self._on_filter_changed)
        filter_layout.addWidget(self.filter_combo)

        self.chk_show_fit = QCheckBox("Show Fit Lines")
        self.chk_show_fit.setChecked(True)
        self.chk_show_fit.stateChanged.connect(self._on_filter_changed)
        filter_layout.addWidget(self.chk_show_fit)

        filter_layout.addStretch()
        layout.addLayout(filter_layout)

        plot_group = QGroupBox("Fit Diagnostics")
        plot_layout = QVBoxLayout(plot_group)
        self.figure = Figure(figsize=(7, 5))
        self.canvas = FigureCanvas(self.figure)
        plot_layout.addWidget(NavigationToolbar(self.canvas, self))
        plot_layout.addWidget(self.canvas)
        layout.addWidget(plot_group)

        return tab

    @staticmethod
    def _as_float(value, default: float) -> float:
        try:
            if value is None:
                return float(default)
            out = float(value)
            return out if np.isfinite(out) else float(default)
        except Exception:
            return float(default)

    @staticmethod
    def _as_int(value, default: int) -> int:
        try:
            if value is None:
                return int(default)
            return int(value)
        except Exception:
            return int(default)

    @staticmethod
    def _first_valid_int(series: pd.Series):
        vals = pd.to_numeric(series, errors="coerce").dropna()
        if vals.empty:
            return pd.NA
        try:
            return int(vals.iloc[0])
        except Exception:
            return pd.NA

    @staticmethod
    def _normalize_dir_path(path) -> str | None:
        if path is None:
            return None
        try:
            return str(Path(path).expanduser().resolve())
        except Exception:
            return str(Path(path).expanduser())

    @staticmethod
    def _safe_offsets(x, y):
        if len(x) == 0:
            return np.empty((0, 2))
        return np.column_stack([np.asarray(x, float), np.asarray(y, float)])

    @staticmethod
    def _pick_first_column(columns, candidates):
        colset = set(columns)
        for cand in candidates:
            if cand in colset:
                return cand
        return None

    def _current_source_result_dir(self) -> Path:
        text = self.source_workspace_edit.text().strip() if hasattr(self, "source_workspace_edit") else ""
        return Path(text) if text else Path(self.result_dir)

    def _selection_path(self, source_dir: Path | None = None) -> Path:
        base = Path(source_dir or self.loaded_source_result_dir or self._current_source_result_dir())
        return tool_extinction_dir(base) / "extinction_star_selection.json"

    def _current_selection_date(self) -> str:
        return str(self.sel_date_combo.currentText()).strip()

    def _current_selection_filter(self) -> str:
        return str(self.sel_filter_combo.currentText()).strip().lower()

    def _current_selection_key(self) -> str | None:
        date_val = self._current_selection_date()
        filt = self._current_selection_filter()
        if not date_val or not filt:
            return None
        return _ext_group_key(date_val, filt)

    def _selection_state_for_key(self, key: str | None) -> dict[str, set[int]]:
        if not key:
            return {"selected": set(), "rejected": set()}
        if key not in self.selection_state:
            self.selection_state[key] = {"selected": set(), "rejected": set()}
        return self.selection_state[key]

    def _on_source_workspace_changed(self):
        current = self._normalize_dir_path(self._current_source_result_dir())
        loaded = self._normalize_dir_path(self.loaded_source_result_dir)
        if loaded and current != loaded:
            self.phot_df = None
            self._clear_selection_context()

    def _browse_source_workspace(self):
        start_dir = self._current_source_result_dir()
        start = str(start_dir if start_dir.exists() else Path(self.result_dir).parent)
        path = QFileDialog.getExistingDirectory(self, "Source workspace 선택", start)
        if path:
            self.source_workspace_edit.setText(path)
            self.log(f"Source workspace selected: {path}")

    def _use_current_source_workspace(self):
        self.source_workspace_edit.setText(str(self.result_dir))
        self.log(f"Source workspace reset to current result dir: {self.result_dir}")

    def log(self, message: str):
        timestamp = time.strftime("%H:%M:%S")
        self.log_text.append(f"[{timestamp}] {message}")

    def on_progress(self, current: int, total: int, fname: str):
        pct = int(100 * current / max(total, 1))
        self.progress_bar.setValue(pct)
        self.progress_label.setText(f"{current}/{total} | {fname}")

    def _clear_selection_context(self):
        self.selection_stats_df = pd.DataFrame()
        self.selection_state = {}
        self._selection_meta_df = pd.DataFrame()
        self._step8_target_by_filter = {}
        self.selection_frame_df = pd.DataFrame()
        self.selection_image_data = None
        self.selection_header = None
        self.selection_selected_source_id = None
        self._selection_sid_to_row = {}
        self.sel_date_combo.blockSignals(True)
        self.sel_filter_combo.blockSignals(True)
        self.sel_frame_combo.blockSignals(True)
        self.sel_date_combo.clear()
        self.sel_filter_combo.clear()
        self.sel_frame_combo.clear()
        self.sel_date_combo.blockSignals(False)
        self.sel_filter_combo.blockSignals(False)
        self.sel_frame_combo.blockSignals(False)
        self.selection_table.setRowCount(0)
        self.selection_summary_label.setText("No Step 7 source table loaded.")
        self.selection_selected_label.setText("Selected: (none)")
        self.selection_figure.clear()
        self.selection_canvas.draw()
        self.selection_preview_figure.clear()
        self.selection_preview_canvas.draw()

    def _source_path_map_data(self, source_dir: Path | None = None) -> dict[str, str]:
        source_dir = Path(source_dir or self._current_source_result_dir())
        source_key = self._normalize_dir_path(source_dir)
        if self._source_file_path_map is not None and self._source_file_path_map_dir == source_key:
            return self._source_file_path_map

        raw = None
        try:
            if source_key == self._normalize_dir_path(self.result_dir):
                raw = getattr(self.params.P, "file_path_map", None)
        except Exception:
            raw = None

        if isinstance(raw, dict) and raw:
            self._source_file_path_map = {str(k): str(v) for k, v in raw.items() if k and v}
        else:
            self._source_file_path_map = load_file_path_map(source_dir)
        self._source_file_path_map_dir = source_key
        return self._source_file_path_map

    def _resolve_source_fits_path(self, fname: str, source_dir: Path | None = None) -> Path | None:
        source_dir = Path(source_dir or self._current_source_result_dir())
        fname = str(fname).strip()
        if not fname:
            return None
        mapped = self._source_path_map_data(source_dir).get(fname)
        candidates: list[Path] = []
        if mapped:
            candidates.append(Path(mapped))
        try:
            if self._normalize_dir_path(source_dir) == self._normalize_dir_path(self.result_dir):
                candidates.append(Path(self.params.get_file_path(fname)))
        except Exception:
            pass
        candidates.extend([
            self.data_dir / fname,
            source_dir / fname,
            step2_cropped_dir(source_dir) / fname,
        ])
        seen: set[str] = set()
        for cand in candidates:
            key = str(cand)
            if key in seen:
                continue
            seen.add(key)
            if cand.exists():
                return cand
        return None

    def _load_saved_selection_state(self, source_dir: Path):
        self.selection_state = {}
        path = self._selection_path(source_dir)
        if not path.exists():
            self.log("No saved extinction-star selection found for this source workspace.")
            return
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            groups = payload.get("groups", {}) if isinstance(payload, dict) else {}
            for key, item in groups.items():
                if not isinstance(item, dict):
                    continue
                selected = {int(v) for v in item.get("selected", []) if v is not None}
                rejected = {int(v) for v in item.get("rejected", []) if v is not None}
                if selected or rejected:
                    self.selection_state[str(key)] = {"selected": selected, "rejected": rejected}
            self.log(f"Loaded saved extinction-star selection: {path}")
        except Exception as e:
            self.log(f"[WARN] Failed to load saved selection: {e}")

    def _save_selection_state(self):
        source_dir = self.loaded_source_result_dir or self._current_source_result_dir()
        path = self._selection_path(source_dir)
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            groups = {}
            for key, state in self.selection_state.items():
                selected = sorted(int(v) for v in state.get("selected", set()))
                rejected = sorted(int(v) for v in state.get("rejected", set()))
                if selected or rejected:
                    groups[str(key)] = {"selected": selected, "rejected": rejected}
            payload = {"version": 1, "groups": groups}
            path.write_text(json.dumps(payload, indent=2, ensure_ascii=True), encoding="utf-8")
        except Exception as e:
            self.log(f"[WARN] Failed to save extinction-star selection: {e}")

    def _load_step8_target_hints(self, source_dir: Path):
        self._step8_target_by_filter = {}
        selection_dir = step8_selection_dir(source_dir)
        if not selection_dir.exists():
            return
        for sel_path in sorted(selection_dir.glob("selection_*.json")):
            filt = sel_path.stem.replace("selection_", "").strip().lower()
            if not filt:
                continue
            try:
                payload = json.loads(sel_path.read_text(encoding="utf-8"))
                sid = payload.get("target_source_id")
                if sid is None:
                    continue
                self._step8_target_by_filter[filt] = int(sid)
            except Exception:
                continue
        if self._step8_target_by_filter:
            self.log(f"Step 8 target hints loaded: {sorted(self._step8_target_by_filter.keys())}")

    def _merge_meta_columns(self, base: pd.DataFrame, extra: pd.DataFrame) -> pd.DataFrame:
        if extra is None or extra.empty:
            return base
        out = base.merge(extra, on="source_id", how="left", suffixes=("", "__new"))
        for col in list(out.columns):
            if not col.endswith("__new"):
                continue
            base_col = col[:-5]
            if base_col in out.columns:
                out[base_col] = out[base_col].where(out[base_col].notna(), out[col])
                out = out.drop(columns=[col])
            else:
                out = out.rename(columns={col: base_col})
        return out

    def _load_selection_metadata(self, source_dir: Path, phot_df: pd.DataFrame) -> pd.DataFrame:
        df = phot_df.copy()
        df = df[df["source_id"].notna()].copy()
        if df.empty:
            return pd.DataFrame(columns=["source_id", "ID", "gaia_G", "gaia_BP_RP", "gaia_variable_flag"])

        df["source_id"] = coerce_int64_source_id(df["source_id"]).astype("Int64")
        base = (
            df.groupby("source_id", dropna=False)
            .agg(ID=("ID", self._first_valid_int))
            .reset_index()
        )
        base = base[base["source_id"].notna()].copy()
        base["source_id"] = base["source_id"].astype("int64")

        ref_path = step6_refbuild_dir(source_dir) / "ref_catalog.tsv"
        ref_meta = pd.DataFrame()
        if ref_path.exists():
            try:
                ref_df = read_csv_int64_source_id(ref_path, sep="\t")
                if "source_id" in ref_df.columns:
                    ref_meta = pd.DataFrame({
                        "source_id": coerce_int64_source_id(ref_df["source_id"]).astype("Int64"),
                    })
                    if "ID" in ref_df.columns:
                        ref_meta["ID"] = pd.to_numeric(ref_df["ID"], errors="coerce").astype("Int64")
                    g_col = self._pick_first_column(ref_df.columns, ("gaia_G", "gaia_g", "phot_g_mean_mag"))
                    bp_col = self._pick_first_column(ref_df.columns, ("gaia_BP", "gaia_bp", "phot_bp_mean_mag"))
                    rp_col = self._pick_first_column(ref_df.columns, ("gaia_RP", "gaia_rp", "phot_rp_mean_mag"))
                    color_col = self._pick_first_column(ref_df.columns, ("bp_rp", "BP_RP", "color_gr"))
                    if g_col:
                        ref_meta["gaia_G"] = pd.to_numeric(ref_df[g_col], errors="coerce")
                    if color_col:
                        ref_meta["gaia_BP_RP"] = pd.to_numeric(ref_df[color_col], errors="coerce")
                    elif bp_col and rp_col:
                        bp = pd.to_numeric(ref_df[bp_col], errors="coerce")
                        rp = pd.to_numeric(ref_df[rp_col], errors="coerce")
                        ref_meta["gaia_BP_RP"] = bp - rp
                    ref_meta = ref_meta[ref_meta["source_id"].notna()].copy()
                    ref_meta["source_id"] = ref_meta["source_id"].astype("int64")
            except Exception as e:
                self.log(f"[WARN] Failed to load ref_catalog metadata: {e}")

        gaia_path = step5_wcs_dir(source_dir) / "gaia_fov.ecsv"
        gaia_meta = pd.DataFrame()
        if gaia_path.exists():
            try:
                gaia_df = Table.read(gaia_path, format="ascii.ecsv").to_pandas()
                if "source_id" in gaia_df.columns:
                    gaia_meta = pd.DataFrame({
                        "source_id": coerce_int64_source_id(gaia_df["source_id"]).astype("Int64"),
                        "gaia_G": pd.to_numeric(gaia_df.get("phot_g_mean_mag"), errors="coerce"),
                        "gaia_BP_RP": pd.to_numeric(gaia_df.get("bp_rp"), errors="coerce"),
                        "gaia_variable_flag": gaia_df.get("phot_variable_flag", pd.Series([""] * len(gaia_df))),
                    })
                    if gaia_meta["gaia_BP_RP"].isna().all():
                        bp = pd.to_numeric(gaia_df.get("phot_bp_mean_mag"), errors="coerce")
                        rp = pd.to_numeric(gaia_df.get("phot_rp_mean_mag"), errors="coerce")
                        gaia_meta["gaia_BP_RP"] = bp - rp
                    gaia_meta = gaia_meta[gaia_meta["source_id"].notna()].copy()
                    gaia_meta["source_id"] = gaia_meta["source_id"].astype("int64")
            except Exception as e:
                self.log(f"[WARN] Failed to load Gaia metadata: {e}")

        meta = base.copy()
        meta = self._merge_meta_columns(meta, ref_meta)
        meta = self._merge_meta_columns(meta, gaia_meta)
        if "gaia_G" not in meta.columns:
            meta["gaia_G"] = np.nan
        if "gaia_BP_RP" not in meta.columns:
            meta["gaia_BP_RP"] = np.nan
        if "gaia_variable_flag" not in meta.columns:
            meta["gaia_variable_flag"] = ""
        meta["gaia_variable_flag"] = meta["gaia_variable_flag"].fillna("").astype(str).str.strip().str.upper()
        return meta

    def _robust_weighted_linfit_local(self, x, y, w=None, clip_sigma=3.0, iters=5, min_n=5):
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        if w is not None:
            w = np.asarray(w, float)
            mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
        else:
            mask = np.isfinite(x) & np.isfinite(y)
        x = x[mask]
        y = y[mask]
        w = w[mask] if w is not None else None
        if len(x) < min_n:
            return (np.nan, np.nan, np.nan, 0, np.nan)

        def _fit(xx, yy, ww):
            if ww is None:
                k, zp = np.polyfit(xx, yy, 1)
            else:
                k, zp = np.polyfit(xx, yy, 1, w=ww)
            return k, zp

        k, zp = _fit(x, y, w)
        base_n = len(x)
        for _ in range(int(iters)):
            yhat = zp + k * x
            resid = y - yhat
            med = np.nanmedian(resid)
            mad = np.nanmedian(np.abs(resid - med)) + 1e-12
            sig = MAD_TO_SIGMA * mad
            keep = np.abs(resid - med) <= float(clip_sigma) * sig
            if keep.sum() < min_n or keep.sum() == len(x):
                break
            x = x[keep]
            y = y[keep]
            if w is not None:
                w = w[keep]
            k, zp = _fit(x, y, w)
        yhat = zp + k * x
        scatter = float(np.nanstd(y - yhat)) if len(x) else np.nan
        out_frac = float(1.0 - (len(x) / max(base_n, 1)))
        return (float(k), float(zp), scatter, int(len(x)), out_frac)

    def _rebuild_selection_catalog(self):
        self.selection_stats_df = pd.DataFrame()
        if not isinstance(self.phot_df, pd.DataFrame) or self.phot_df.empty:
            self._clear_selection_context()
            return

        df = self.phot_df.copy()
        df = df[df["source_id"].notna()].copy()
        if df.empty:
            self._clear_selection_context()
            return

        df["source_id"] = coerce_int64_source_id(df["source_id"]).astype("Int64")
        df = df[df["source_id"].notna()].copy()
        df["source_id"] = df["source_id"].astype("int64")

        clip_sigma = self._as_float(getattr(self.params.P, "extfit_clip_sigma", 3.0), 3.0)
        fit_iters = self._as_int(getattr(self.params.P, "extfit_fit_iters", 5), 5)
        use_weights = bool(getattr(self.params.P, "extinction_star_use_weights", True))

        rows = []
        for (date_val, filt, sid), ssub in df.groupby(["date", "filter", "source_id"]):
            x = pd.to_numeric(ssub["airmass"], errors="coerce").to_numpy(float)
            y = pd.to_numeric(ssub["mag_inst"], errors="coerce").to_numpy(float)
            n_frames = int(ssub["file"].nunique())
            snr_med = float(np.nanmedian(pd.to_numeric(ssub.get("snr"), errors="coerce"))) if "snr" in ssub.columns else np.nan
            delta_x = _finite_range(x)
            w = None
            if use_weights and "snr" in ssub.columns:
                snr = pd.to_numeric(ssub["snr"], errors="coerce").to_numpy(float)
                sig_mag = MAG_ERR_COEFF / np.clip(snr, 1e-6, None)
                w = 1.0 / np.clip(sig_mag ** 2, 1e-12, None)
            k_i, m0_i, rms_i, _, _ = self._robust_weighted_linfit_local(
                x, y, w=w, clip_sigma=clip_sigma, iters=fit_iters, min_n=min(5, max(n_frames, 3))
            )
            rows.append({
                "date": str(date_val),
                "filter": str(filt).strip().lower(),
                "source_id": int(sid),
                "ID": self._first_valid_int(ssub["ID"]) if "ID" in ssub.columns else pd.NA,
                "n_frames": n_frames,
                "delta_x": delta_x,
                "snr_med": snr_med,
                "k_i": k_i,
                "m0_i": m0_i,
                "rms_i": rms_i,
            })

        stats_df = pd.DataFrame(rows)
        if not stats_df.empty and not self._selection_meta_df.empty:
            stats_df = self._merge_meta_columns(stats_df, self._selection_meta_df)
        if "gaia_G" not in stats_df.columns:
            stats_df["gaia_G"] = np.nan
        if "gaia_BP_RP" not in stats_df.columns:
            stats_df["gaia_BP_RP"] = np.nan
        if "gaia_variable_flag" not in stats_df.columns:
            stats_df["gaia_variable_flag"] = ""
        stats_df["step8_target"] = stats_df.apply(
            lambda r: self._step8_target_by_filter.get(str(r["filter"]).strip().lower()) == int(r["source_id"]),
            axis=1,
        )
        self.selection_stats_df = stats_df
        self._refresh_selection_group_controls()

    def _refresh_selection_group_controls(self):
        current_date = self._current_selection_date()
        current_filter = self._current_selection_filter()
        if self.selection_stats_df.empty:
            self._clear_selection_context()
            return

        dates = sorted(self.selection_stats_df["date"].dropna().astype(str).unique().tolist())
        self.sel_date_combo.blockSignals(True)
        self.sel_date_combo.clear()
        self.sel_date_combo.addItems(dates)
        if current_date in dates:
            self.sel_date_combo.setCurrentText(current_date)
        elif dates:
            self.sel_date_combo.setCurrentIndex(0)
        self.sel_date_combo.blockSignals(False)

        self._update_selection_filter_options(preferred=current_filter)
        self._update_selection_frame_options()
        self._refresh_selection_table()
        self._load_selection_frame()

    def _update_selection_filter_options(self, preferred: str | None = None):
        date_val = self._current_selection_date()
        if not date_val:
            return
        sub = self.selection_stats_df[self.selection_stats_df["date"] == date_val]
        filters = sorted(sub["filter"].dropna().astype(str).unique().tolist())
        self.sel_filter_combo.blockSignals(True)
        self.sel_filter_combo.clear()
        self.sel_filter_combo.addItems(filters)
        if preferred and preferred in filters:
            self.sel_filter_combo.setCurrentText(preferred)
        elif filters:
            self.sel_filter_combo.setCurrentIndex(0)
        self.sel_filter_combo.blockSignals(False)

    def _default_frame_for_group(self, group_df: pd.DataFrame) -> str | None:
        if group_df.empty:
            return None
        frame_df = (
            group_df.groupby("file")
            .agg(airmass=("airmass", "median"), n_star=("source_id", "nunique"))
            .reset_index()
        )
        if frame_df.empty:
            return None
        med_x = float(np.nanmedian(pd.to_numeric(frame_df["airmass"], errors="coerce")))
        frame_df["sort_dx"] = np.abs(pd.to_numeric(frame_df["airmass"], errors="coerce") - med_x)
        frame_df = frame_df.sort_values(["sort_dx", "n_star", "file"], ascending=[True, False, True])
        return str(frame_df.iloc[0]["file"])

    def _update_selection_frame_options(self):
        date_val = self._current_selection_date()
        filt = self._current_selection_filter()
        current_file = str(self.sel_frame_combo.currentText()).strip()
        if not date_val or not filt or not isinstance(self.phot_df, pd.DataFrame):
            return
        group_df = self.phot_df[
            (self.phot_df["date"].astype(str) == date_val)
            & (self.phot_df["filter"].astype(str).str.strip().str.lower() == filt)
        ].copy()
        files = sorted(group_df["file"].dropna().astype(str).unique().tolist())
        preferred = current_file if current_file in files else self._default_frame_for_group(group_df)
        self.sel_frame_combo.blockSignals(True)
        self.sel_frame_combo.clear()
        self.sel_frame_combo.addItems(files)
        if preferred and preferred in files:
            self.sel_frame_combo.setCurrentText(preferred)
        elif files:
            self.sel_frame_combo.setCurrentIndex(0)
        self.sel_frame_combo.blockSignals(False)

    def _current_group_stats(self) -> pd.DataFrame:
        if self.selection_stats_df.empty:
            return pd.DataFrame()
        date_val = self._current_selection_date()
        filt = self._current_selection_filter()
        if not date_val or not filt:
            return pd.DataFrame()
        return self.selection_stats_df[
            (self.selection_stats_df["date"] == date_val)
            & (self.selection_stats_df["filter"] == filt)
        ].copy()

    def _status_for_source_id(self, sid: int) -> str:
        state = self._selection_state_for_key(self._current_selection_key())
        sid = int(sid)
        if sid in state["selected"]:
            return "use"
        if sid in state["rejected"]:
            return "reject"
        return "candidate"

    def _refresh_selection_table(self):
        sub = self._current_group_stats()
        self.selection_table.setRowCount(0)
        self._selection_sid_to_row = {}
        self._update_selection_summary()
        self._update_selected_star_label()
        self._refresh_selected_star_preview()
        self._update_selection_overlay()
        if sub.empty:
            return

        sub["status"] = sub["source_id"].map(self._status_for_source_id)
        sub["status_rank"] = sub["status"].map({"use": 0, "candidate": 1, "reject": 2}).fillna(9)
        sub["rms_sort"] = pd.to_numeric(sub["rms_i"], errors="coerce").fillna(np.inf)
        sub["snr_sort"] = pd.to_numeric(sub["snr_med"], errors="coerce").fillna(-np.inf)
        sub = sub.sort_values(
            ["status_rank", "step8_target", "rms_sort", "snr_sort", "source_id"],
            ascending=[True, True, True, False, True],
        ).reset_index(drop=True)

        self.selection_table.setRowCount(len(sub))
        for row_idx, row in sub.iterrows():
            sid = int(row["source_id"])
            state_text = str(row["status"])
            display_id = self._display_id_for_sid(sid, fallback=row.get("ID"))

            items = [
                QTableWidgetItem(state_text),
                QTableWidgetItem("" if pd.isna(display_id) else str(display_id)),
                QTableWidgetItem(str(sid)),
                QTableWidgetItem(f"{float(row['gaia_G']):.2f}" if np.isfinite(row["gaia_G"]) else "-"),
                QTableWidgetItem(f"{float(row['gaia_BP_RP']):.2f}" if np.isfinite(row["gaia_BP_RP"]) else "-"),
                QTableWidgetItem("VAR" if str(row.get("gaia_variable_flag", "")) == "VARIABLE" else ""),
                QTableWidgetItem(str(int(row["n_frames"]))),
                QTableWidgetItem(f"{float(row['delta_x']):.3f}" if np.isfinite(row["delta_x"]) else "-"),
                QTableWidgetItem(f"{float(row['snr_med']):.1f}" if np.isfinite(row["snr_med"]) else "-"),
                QTableWidgetItem(f"{float(row['rms_i']):.4f}" if np.isfinite(row["rms_i"]) else "-"),
                QTableWidgetItem(f"{float(row['k_i']):.4f}" if np.isfinite(row["k_i"]) else "-"),
            ]
            for col_idx, item in enumerate(items):
                if col_idx == 0:
                    if state_text == "use":
                        item.setBackground(QColor("#C8E6C9"))
                    elif state_text == "reject":
                        item.setBackground(QColor("#FFCDD2"))
                    else:
                        item.setBackground(QColor("#ECEFF1"))
                if col_idx == 5 and item.text() == "VAR":
                    item.setForeground(QColor("#C62828"))
                if bool(row.get("step8_target", False)):
                    item.setToolTip("Step 8 target hint")
                self.selection_table.setItem(row_idx, col_idx, item)
            self._selection_sid_to_row[sid] = row_idx

        if self.selection_selected_source_id in self._selection_sid_to_row:
            self._select_selection_row_by_sid(self.selection_selected_source_id, keep_focus=False)
        self._update_selection_summary()
        self._update_selected_star_label()
        self._refresh_selected_star_preview()
        self._update_selection_overlay()

    def _update_selection_summary(self):
        sub = self._current_group_stats()
        if sub.empty:
            self.selection_summary_label.setText("No date/filter group loaded.")
            return
        state = self._selection_state_for_key(self._current_selection_key())
        n_total = len(sub)
        n_use = len(state["selected"])
        n_reject = len(state["rejected"])
        n_cand = max(n_total - n_use - n_reject, 0)
        chosen = sub[sub["source_id"].isin(state["selected"])] if state["selected"] else sub
        color_span = _finite_range(chosen["gaia_BP_RP"]) if "gaia_BP_RP" in chosen.columns else np.nan
        dx_span = _finite_range(chosen["delta_x"]) if "delta_x" in chosen.columns else np.nan
        parts = [
            f"{self._current_selection_date()} / {self._current_selection_filter()}",
            f"use={n_use}",
            f"reject={n_reject}",
            f"candidate={n_cand}",
            f"total={n_total}",
        ]
        if np.isfinite(color_span):
            parts.append(f"color span={color_span:.2f}")
        if np.isfinite(dx_span):
            parts.append(f"ΔX span={dx_span:.2f}")
        self.selection_summary_label.setText(" | ".join(parts))

    def _table_selected_source_ids(self) -> set[int]:
        rows = self.selection_table.selectionModel().selectedRows()
        sids = set()
        for row in rows:
            item = self.selection_table.item(row.row(), 2)
            if item is None:
                continue
            try:
                sids.add(int(item.text()))
            except Exception:
                continue
        if not sids and self.selection_selected_source_id is not None:
            sids.add(int(self.selection_selected_source_id))
        return sids

    def _apply_state_to_selected(self, new_state: str):
        key = self._current_selection_key()
        if not key:
            return
        sids = self._table_selected_source_ids()
        if not sids:
            QMessageBox.information(self, "Selection", "Select one or more stars first.")
            return
        state = self._selection_state_for_key(key)
        if new_state == "selected":
            state["selected"] |= sids
            state["rejected"] -= sids
        elif new_state == "rejected":
            state["rejected"] |= sids
            state["selected"] -= sids
        else:
            state["selected"] -= sids
            state["rejected"] -= sids
        self._save_selection_state()
        self._refresh_selection_table()
        self.log(f"{self._current_selection_date()}/{self._current_selection_filter()}: {new_state} -> {sorted(sids)}")

    def _clear_current_group_states(self):
        key = self._current_selection_key()
        if not key:
            return
        self.selection_state[key] = {"selected": set(), "rejected": set()}
        self._save_selection_state()
        self._refresh_selection_table()
        self.log(f"{self._current_selection_date()}/{self._current_selection_filter()}: manual state cleared")

    def _copy_selection_state_to_keys(self, keys: list[str], scope_label: str):
        src_key = self._current_selection_key()
        if not src_key:
            return
        src_state = self._selection_state_for_key(src_key)
        selected = set(src_state.get("selected", set()))
        rejected = set(src_state.get("rejected", set()))
        if not selected and not rejected:
            QMessageBox.information(self, "Copy Selection", "Current group has no manual use/reject state to copy.")
            return

        copied = 0
        for key in keys:
            if not key or key == src_key:
                continue
            self.selection_state[key] = {
                "selected": set(selected),
                "rejected": set(rejected),
            }
            copied += 1

        if copied == 0:
            QMessageBox.information(self, "Copy Selection", "No target groups available for copy.")
            return

        self._save_selection_state()
        self._refresh_selection_table()
        self.log(
            f"{self._current_selection_date()}/{self._current_selection_filter()}: "
            f"copied manual state to {copied} {scope_label} group(s)"
        )

    def _copy_current_selection_to_same_date(self):
        date_val = self._current_selection_date()
        filt = self._current_selection_filter()
        if not date_val or not filt or self.selection_stats_df.empty:
            return
        keys = [
            _ext_group_key(date_val, str(target_filt))
            for target_filt in sorted(
                self.selection_stats_df.loc[
                    self.selection_stats_df["date"] == date_val, "filter"
                ].dropna().astype(str).unique().tolist()
            )
        ]
        self._copy_selection_state_to_keys(keys, "same-date")

    def _copy_current_selection_to_all_groups(self):
        if self.selection_stats_df.empty:
            return
        keys = [
            _ext_group_key(str(row["date"]), str(row["filter"]))
            for _, row in (
                self.selection_stats_df[["date", "filter"]]
                .drop_duplicates()
                .sort_values(["date", "filter"])
                .iterrows()
            )
        ]
        self._copy_selection_state_to_keys(keys, "all")

    def _auto_pick_current_group(self):
        sub = self._current_group_stats()
        key = self._current_selection_key()
        if sub.empty or not key:
            return
        star_min_frames = self._as_int(getattr(self.params.P, "extinction_star_min_frames", 8), 8)
        star_rms_max = self._as_float(getattr(self.params.P, "extinction_star_rms_max", 0.10), 0.10)
        star_snr_min = self._as_float(getattr(self.params.P, "extinction_star_snr_med_min", 10.0), 10.0)
        filt = self._current_selection_filter()

        total = len(sub)
        mask = sub["n_frames"] >= star_min_frames
        n_after_frames = int(mask.sum())
        if star_snr_min > 0:
            mask &= pd.to_numeric(sub["snr_med"], errors="coerce") >= star_snr_min
        n_after_snr = int(mask.sum())
        if star_rms_max > 0:
            mask &= pd.to_numeric(sub["rms_i"], errors="coerce").fillna(np.inf) <= star_rms_max
        n_after_rms = int(mask.sum())
        if self.chk_exclude_gaia_var.isChecked():
            mask &= sub["gaia_variable_flag"].astype(str).str.upper().ne("VARIABLE")
        n_after_var = int(mask.sum())
        if self.chk_exclude_step8_target.isChecked():
            target_sid = self._step8_target_by_filter.get(filt)
            if target_sid is not None:
                mask &= sub["source_id"] != int(target_sid)
        n_after_target = int(mask.sum())

        state = self._selection_state_for_key(key)
        rejected = set(state["rejected"])
        chosen = {int(v) for v in sub.loc[mask, "source_id"].tolist() if int(v) not in rejected}
        state["selected"] = chosen
        state["rejected"] -= chosen
        self._save_selection_state()
        self._refresh_selection_table()
        self.log(
            f"{self._current_selection_date()}/{self._current_selection_filter()}: "
            f"auto-picked {len(chosen)} star(s) "
            f"[total={total}, frames={n_after_frames}, snr={n_after_snr}, "
            f"rms={n_after_rms}, var={n_after_var}, target={n_after_target}]"
        )

    def _on_selection_group_changed(self):
        self._update_selection_filter_options(preferred=self._current_selection_filter())
        self._update_selection_frame_options()
        self.selection_selected_source_id = None
        self._refresh_selection_table()
        self._load_selection_frame()

    def _on_selection_frame_changed(self):
        self._load_selection_frame()

    def _extract_date_key(self, value: str) -> str:
        match = self._DATE_KEY_RE.search(str(value))
        return match.group(1) if match else ""

    def _resolve_idmatch_path(self, source_dir: Path, fname: str) -> Path | None:
        # New pipeline: master_id is a column in forced phot TSV; no separate idmatch file.
        # Check forced_phot dir for per-frame photometry TSV.
        forced_dir = step7_forced_phot_dir(source_dir)
        candidates = [forced_dir / f"photometry_{fname}.tsv"]
        for path in candidates:
            if path.exists():
                return path
        return None

    def _load_selection_frame(self):
        self.selection_frame_df = pd.DataFrame()
        self.selection_image_data = None
        self.selection_header = None

        fname = str(self.sel_frame_combo.currentText()).strip()
        filt = self._current_selection_filter()
        source_dir = self.loaded_source_result_dir or self._current_source_result_dir()
        if not fname or not filt:
            self._update_selection_overlay()
            self.selection_figure.clear()
            self.selection_canvas.draw()
            return

        try:
            fpath = self._resolve_source_fits_path(fname, source_dir)
            if fpath is None or not fpath.exists():
                raise FileNotFoundError(f"Cannot find {fname}")
            with fits.open(fpath, memmap=False) as hdul:
                data = hdul[0].data
                self.selection_header = hdul[0].header
            arr = np.asarray(data)
            while arr.ndim > 2:
                arr = arr[0]
            self.selection_image_data = np.asarray(arr, float)
        except Exception as e:
            self.log(f"[WARN] Failed to load representative frame {fname}: {e}")

        try:
            frame_df = load_frame_photometry(source_dir, fname, filt)
        except Exception:
            frame_df = None
        overlay = pd.DataFrame()
        if frame_df is not None and not frame_df.empty:
            frame_df = frame_df.copy()
            if "source_id" in frame_df.columns:
                frame_df["source_id"] = coerce_int64_source_id(frame_df["source_id"]).astype("Int64")
            if {"x", "y", "source_id"} <= set(frame_df.columns):
                overlay = frame_df[["x", "y", "source_id"]].copy()
                if "ID" in frame_df.columns:
                    overlay["ID"] = pd.to_numeric(frame_df["ID"], errors="coerce").astype("Int64")

        if overlay.empty:
            idm_path = self._resolve_idmatch_path(source_dir, fname)
            if idm_path is not None and idm_path.exists():
                try:
                    idm = read_csv_int64_source_id(idm_path)
                    if {"x", "y", "source_id"} <= set(idm.columns):
                        overlay = idm[["x", "y", "source_id"]].copy()
                        overlay["source_id"] = coerce_int64_source_id(overlay["source_id"]).astype("Int64")
                except Exception as e:
                    self.log(f"[WARN] Failed to load idmatch overlay {idm_path.name}: {e}")

        if not overlay.empty:
            overlay = overlay[overlay["source_id"].notna()].copy()
            overlay["source_id"] = overlay["source_id"].astype("int64")
            if "ID" not in overlay.columns:
                overlay["ID"] = overlay["source_id"].map(
                    dict(zip(self._selection_meta_df["source_id"], self._selection_meta_df["ID"]))
                ).astype("Int64") if not self._selection_meta_df.empty else pd.Series(pd.NA, index=overlay.index, dtype="Int64")
        self.selection_frame_df = overlay
        self._refresh_selection_image()
        self._update_selection_overlay()

    def _stretch_selection_image(self, data):
        arr = np.asarray(data, float)
        finite = arr[np.isfinite(arr)]
        if finite.size == 0:
            return np.zeros_like(arr, dtype=float)
        lo = float(np.nanpercentile(finite, 1.0))
        hi = float(np.nanpercentile(finite, 99.5))
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            lo = float(np.nanmin(finite))
            hi = float(np.nanmax(finite))
        denom = max(hi - lo, 1e-12)
        return np.clip((arr - lo) / denom, 0.0, 1.0)

    def _refresh_selection_image(self):
        self.selection_figure.clear()
        self.selection_ax = self.selection_figure.add_subplot(111)
        self._selection_imshow = None
        self._selection_scat_candidate = None
        self._selection_scat_selected = None
        self._selection_scat_rejected = None
        self._selection_scat_highlight = None

        if self.selection_image_data is None:
            self.selection_ax.text(0.5, 0.5, "No frame loaded", ha="center", va="center", transform=self.selection_ax.transAxes)
            self.selection_canvas.draw()
            return

        stretched = self._stretch_selection_image(self.selection_image_data)
        self._selection_imshow = self.selection_ax.imshow(
            stretched, cmap="gray", origin="lower", vmin=0, vmax=1, interpolation="nearest"
        )
        title = str(self.sel_frame_combo.currentText()).strip()
        self.selection_ax.set_title(title)
        self.selection_ax.set_xlabel("X")
        self.selection_ax.set_ylabel("Y")
        self._selection_scat_candidate = self.selection_ax.scatter(
            [], [], s=26, facecolors="none", edgecolors="#90A4AE", linewidths=0.9, alpha=0.8
        )
        self._selection_scat_selected = self.selection_ax.scatter(
            [], [], s=42, facecolors="none", edgecolors="#2E7D32", linewidths=1.5, alpha=0.95
        )
        self._selection_scat_rejected = self.selection_ax.scatter(
            [], [], s=42, facecolors="none", edgecolors="#C62828", linewidths=1.5, alpha=0.95, marker="x"
        )
        self._selection_scat_highlight = self.selection_ax.scatter(
            [], [], s=80, facecolors="none", edgecolors="#FFD54F", linewidths=1.8, alpha=0.95
        )
        self.selection_canvas.draw()

    def _display_id_for_sid(self, sid: int, fallback=None):
        sid = int(sid)
        if pd.notna(fallback):
            try:
                return int(fallback)
            except Exception:
                pass
        if not self._selection_meta_df.empty:
            sub = self._selection_meta_df[self._selection_meta_df["source_id"] == sid]
            if not sub.empty:
                val = pd.to_numeric(sub.iloc[0].get("ID"), errors="coerce")
                if np.isfinite(val):
                    return int(val)
        return sid

    def _update_selection_overlay(self):
        if self.selection_ax is None:
            return
        if self.selection_frame_df is None or self.selection_frame_df.empty:
            if self._selection_scat_candidate is not None:
                self._selection_scat_candidate.set_offsets(np.empty((0, 2)))
                self._selection_scat_selected.set_offsets(np.empty((0, 2)))
                self._selection_scat_rejected.set_offsets(np.empty((0, 2)))
                self._selection_scat_highlight.set_offsets(np.empty((0, 2)))
                self.selection_canvas.draw_idle()
            return

        for txt in list(self.selection_ax.texts):
            try:
                txt.remove()
            except Exception:
                pass

        x = pd.to_numeric(self.selection_frame_df["x"], errors="coerce").to_numpy(float)
        y = pd.to_numeric(self.selection_frame_df["y"], errors="coerce").to_numpy(float)
        sids = coerce_int64_source_id(self.selection_frame_df["source_id"]).fillna(-1).astype("int64").to_numpy()

        state = self._selection_state_for_key(self._current_selection_key())
        is_selected = np.array([sid in state["selected"] for sid in sids], dtype=bool)
        is_rejected = np.array([sid in state["rejected"] for sid in sids], dtype=bool)
        is_candidate = ~(is_selected | is_rejected)
        is_highlight = sids == int(self.selection_selected_source_id) if self.selection_selected_source_id is not None else np.zeros(len(sids), dtype=bool)

        self._selection_scat_candidate.set_offsets(self._safe_offsets(x[is_candidate], y[is_candidate]))
        self._selection_scat_selected.set_offsets(self._safe_offsets(x[is_selected], y[is_selected]))
        self._selection_scat_rejected.set_offsets(self._safe_offsets(x[is_rejected], y[is_rejected]))
        self._selection_scat_highlight.set_offsets(self._safe_offsets(x[is_highlight], y[is_highlight]))

        if self.chk_show_selection_ids.isChecked():
            label_mask = is_selected | is_rejected | is_highlight
            for xi, yi, sid in zip(x[label_mask], y[label_mask], sids[label_mask]):
                if not np.isfinite(xi) or not np.isfinite(yi):
                    continue
                disp_id = self._display_id_for_sid(int(sid))
                self.selection_ax.text(
                    xi - 4.0, yi + 4.0, str(disp_id),
                    color="#FFD54F", fontsize=8, fontweight="bold",
                    ha="right", va="bottom", alpha=0.98,
                )
        self.selection_canvas.draw_idle()

    def _on_selection_scroll(self, event):
        if event.inaxes != self.selection_ax or event.xdata is None or event.ydata is None:
            return
        scale = 1.2 if event.button == "down" else 1 / 1.2
        xlim = self.selection_ax.get_xlim()
        ylim = self.selection_ax.get_ylim()
        xdata, ydata = event.xdata, event.ydata
        new_width = (xlim[1] - xlim[0]) * scale
        new_height = (ylim[1] - ylim[0]) * scale
        relx = (xlim[1] - xdata) / (xlim[1] - xlim[0])
        rely = (ylim[1] - ydata) / (ylim[1] - ylim[0])
        self.selection_ax.set_xlim([xdata - new_width * (1 - relx), xdata + new_width * relx])
        self.selection_ax.set_ylim([ydata - new_height * (1 - rely), ydata + new_height * rely])
        self.selection_canvas.draw_idle()

    def _on_selection_button_press(self, event):
        if event.inaxes != self.selection_ax:
            return
        if event.button == 3:
            self._selection_panning = True
            self._selection_pan_start = (event.xdata, event.ydata)

    def _on_selection_button_release(self, event):
        if event.button == 3:
            self._selection_panning = False
            self._selection_pan_start = None

    def _on_selection_motion(self, event):
        if not self._selection_panning or event.inaxes != self.selection_ax:
            return
        if self._selection_pan_start is None or event.xdata is None or event.ydata is None:
            return
        dx = self._selection_pan_start[0] - event.xdata
        dy = self._selection_pan_start[1] - event.ydata
        xlim = self.selection_ax.get_xlim()
        ylim = self.selection_ax.get_ylim()
        self.selection_ax.set_xlim([xlim[0] + dx, xlim[1] + dx])
        self.selection_ax.set_ylim([ylim[0] + dy, ylim[1] + dy])
        self.selection_canvas.draw_idle()

    def _on_selection_click(self, event):
        if event.inaxes != self.selection_ax or event.button != 1:
            return
        if self.selection_frame_df is None or self.selection_frame_df.empty:
            return
        x = event.xdata
        y = event.ydata
        if x is None or y is None:
            return
        search_r = float(getattr(self.params.P, "search_radius_px", 7.0))
        dx = pd.to_numeric(self.selection_frame_df["x"], errors="coerce").to_numpy(float) - x
        dy = pd.to_numeric(self.selection_frame_df["y"], errors="coerce").to_numpy(float) - y
        dist2 = dx * dx + dy * dy
        valid = np.isfinite(dist2)
        if dist2.size == 0 or not np.any(valid):
            return
        i = int(np.argmin(np.where(valid, dist2, np.inf)))
        if dist2[i] > search_r * search_r:
            self.selection_selected_source_id = None
            self.selection_table.clearSelection()
            self._update_selected_star_label()
            self._refresh_selected_star_preview()
            self._update_selection_overlay()
            return
        sid = int(coerce_int64_source_id(pd.Series([self.selection_frame_df.iloc[i]["source_id"]])).iloc[0])
        self._select_source_id(sid)

    def _on_selection_table_changed(self):
        rows = self.selection_table.selectionModel().selectedRows()
        if not rows:
            return
        item = self.selection_table.item(rows[0].row(), 2)
        if item is None:
            return
        try:
            sid = int(item.text())
        except Exception:
            return
        self._select_source_id(sid, sync_table=False)

    def _select_selection_row_by_sid(self, sid: int, keep_focus: bool = True):
        row = self._selection_sid_to_row.get(int(sid))
        if row is None:
            return
        self.selection_table.blockSignals(True)
        self.selection_table.selectRow(row)
        self.selection_table.setCurrentCell(row, 0)
        item = self.selection_table.item(row, 0)
        if keep_focus and item is not None:
            self.selection_table.scrollToItem(item, QAbstractItemView.PositionAtCenter)
        self.selection_table.blockSignals(False)

    def _select_source_id(self, sid: int, sync_table: bool = True):
        self.selection_selected_source_id = int(sid)
        if sync_table:
            self._select_selection_row_by_sid(int(sid))
        self._update_selected_star_label()
        self._refresh_selected_star_preview()
        self._update_selection_overlay()

    def _update_selected_star_label(self):
        sid = self.selection_selected_source_id
        if sid is None:
            self.selection_selected_label.setText("Selected: (none)")
            return
        sub = self._current_group_stats()
        row = sub[sub["source_id"] == int(sid)]
        if row.empty:
            self.selection_selected_label.setText(f"Selected: {sid}")
            return
        row = row.iloc[0]
        disp_id = self._display_id_for_sid(int(sid), fallback=row.get("ID"))
        txt = f"Selected: ID {disp_id} | source_id={sid}"
        if np.isfinite(row.get("gaia_G", np.nan)):
            txt += f" | G={float(row['gaia_G']):.2f}"
        if np.isfinite(row.get("gaia_BP_RP", np.nan)):
            txt += f" | BP-RP={float(row['gaia_BP_RP']):.2f}"
        txt += f" | Nfr={int(row['n_frames'])}"
        if np.isfinite(row.get("rms_i", np.nan)):
            txt += f" | RMS={float(row['rms_i']):.4f}"
        self.selection_selected_label.setText(txt)

    def _refresh_selected_star_preview(self):
        self.selection_preview_figure.clear()
        ax = self.selection_preview_figure.add_subplot(111)
        sid = self.selection_selected_source_id
        if sid is None or not isinstance(self.phot_df, pd.DataFrame):
            ax.text(0.5, 0.5, "No star selected", ha="center", va="center", transform=ax.transAxes)
            self.selection_preview_canvas.draw()
            return
        date_val = self._current_selection_date()
        filt = self._current_selection_filter()
        sub = self.phot_df[
            (self.phot_df["date"].astype(str) == date_val)
            & (self.phot_df["filter"].astype(str).str.strip().str.lower() == filt)
            & (coerce_int64_source_id(self.phot_df["source_id"]).fillna(-1).astype("int64") == int(sid))
        ].copy()
        if sub.empty:
            ax.text(0.5, 0.5, "No points for selected star", ha="center", va="center", transform=ax.transAxes)
            self.selection_preview_canvas.draw()
            return

        x = pd.to_numeric(sub["airmass"], errors="coerce").to_numpy(float)
        y = pd.to_numeric(sub["mag_inst"], errors="coerce").to_numpy(float)
        ax.scatter(x, y, s=12, alpha=0.7, color="#1f77b4")

        group_stats = self._current_group_stats()
        row = group_stats[group_stats["source_id"] == int(sid)]
        if not row.empty:
            row = row.iloc[0]
            k_i = float(row.get("k_i", np.nan))
            m0_i = float(row.get("m0_i", np.nan))
            if np.isfinite(k_i) and np.isfinite(m0_i):
                xx = np.linspace(np.nanmin(x), np.nanmax(x), 100)
                yy = m0_i + k_i * xx
                ax.plot(xx, yy, color="#D32F2F", lw=1.5, label=f"k_i={k_i:.4f}")
                ax.legend(loc="best", fontsize=8)

        ax.set_xlabel("Airmass")
        ax.set_ylabel("mag_inst")
        ax.set_title("Bouguer Preview")
        ax.grid(True, alpha=0.3)
        ax.invert_yaxis()
        self.selection_preview_figure.tight_layout()
        self.selection_preview_canvas.draw()

    def run_photometry(self):
        if self.worker and self.worker.isRunning():
            return
        source_dir = self._current_source_result_dir()
        if not source_dir.exists():
            QMessageBox.warning(self, "Missing Source Workspace", f"Source workspace not found:\n{source_dir}")
            return
        idx_path = step7_forced_phot_dir(source_dir) / "photometry_index.csv"
        if not idx_path.exists():
            QMessageBox.warning(self, "Missing Step 7", f"Step 7 forced photometry_index.csv not found in:\n{source_dir}")
            return

        self.log_text.clear()
        self.btn_run_phot.setEnabled(False)
        self.btn_run_fit.setEnabled(False)
        self.btn_stop.setEnabled(True)
        self.btn_save.setEnabled(False)
        self.progress_bar.setValue(0)
        self.progress_label.setText("Loading Step 7 source table...")
        self.log(f"Source workspace: {source_dir}")

        self.worker = ExtinctionFitWorker(
            self.params,
            self.data_dir,
            self.result_dir,
            mode="per_star",
            task="photometry",
            source_result_dir=source_dir,
        )
        self.worker.progress.connect(self.on_progress)
        self.worker.log.connect(self.log)
        self.worker.finished.connect(self.on_photometry_finished)
        self.worker.error.connect(self.on_error)
        self.worker.start()

    def on_photometry_finished(self, results: dict):
        self.btn_run_phot.setEnabled(True)
        self.btn_stop.setEnabled(False)
        if results.get("stopped"):
            self.progress_label.setText("Stopped")
            self.log("Step 7 source load stopped")
            self.btn_run_fit.setEnabled(True)
            return

        phot_df = results.get("phot")
        if isinstance(phot_df, pd.DataFrame):
            self.phot_df = phot_df
        source_dir = results.get("source_result_dir")
        self.loaded_source_result_dir = Path(source_dir) if source_dir else self._current_source_result_dir()
        self._source_file_path_map = None
        self._source_file_path_map_dir = None
        self._load_saved_selection_state(self.loaded_source_result_dir)
        self._load_step8_target_hints(self.loaded_source_result_dir)
        self._selection_meta_df = self._load_selection_metadata(self.loaded_source_result_dir, self.phot_df)
        self._rebuild_selection_catalog()
        self.btn_run_fit.setEnabled(True)
        self.progress_label.setText("Step 7 source table ready")
        self.log("Step 7 source load complete")
        self.tabs.setCurrentWidget(self.selection_tab)

    def open_parameters_dialog(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Extinction Parameters")
        dialog.resize(420, 320)
        layout = QVBoxLayout(dialog)

        var_group = QGroupBox("Variable-star filter")
        var_form = QFormLayout(var_group)

        method_combo = QComboBox()
        method_combo.addItems(["MAD", "STD"])
        current_method = str(getattr(self.params.P, "extinction_varstar_method", "mad") or "mad").strip().lower()
        method_combo.setCurrentIndex(1 if current_method == "std" else 0)
        var_form.addRow("Metric:", method_combo)

        var_sigma = QDoubleSpinBox()
        var_sigma.setRange(0.5, 10.0)
        var_sigma.setSingleStep(0.1)
        var_sigma.setValue(float(getattr(self.params.P, "extinction_varstar_sigma", 3.0) or 3.0))
        var_form.addRow("Sigma:", var_sigma)

        var_min_frames = QSpinBox()
        var_min_frames.setRange(3, 100)
        var_min_frames.setValue(int(getattr(self.params.P, "extinction_varstar_min_frames", 5) or 5))
        var_form.addRow("Min frames:", var_min_frames)
        layout.addWidget(var_group)

        qc_group = QGroupBox("Frame QC")
        qc_form = QFormLayout(qc_group)

        qc_method_combo = QComboBox()
        qc_method_combo.addItems(["MAD", "STD"])
        qc_method = str(getattr(self.params.P, "extinction_frame_qc_method", "mad") or "mad").strip().lower()
        qc_method_combo.setCurrentIndex(1 if qc_method == "std" else 0)
        qc_form.addRow("Metric:", qc_method_combo)

        qc_sigma = QDoubleSpinBox()
        qc_sigma.setRange(0.5, 10.0)
        qc_sigma.setSingleStep(0.1)
        qc_sigma.setValue(float(getattr(self.params.P, "extinction_frame_qc_sigma", 3.0) or 3.0))
        qc_form.addRow("Sigma:", qc_sigma)
        layout.addWidget(qc_group)

        method_group = QGroupBox("Per-star Bouguer")
        method_form = QFormLayout(method_group)

        snr_cut = QDoubleSpinBox()
        snr_cut.setRange(0.0, 100.0)
        snr_cut.setSingleStep(0.5)
        snr_cut.setValue(float(getattr(self.params.P, "extinction_snr_min", 10.0) or 10.0))
        method_form.addRow("SNR min:", snr_cut)

        clip_sigma = QDoubleSpinBox()
        clip_sigma.setRange(0.5, 10.0)
        clip_sigma.setSingleStep(0.1)
        clip_sigma.setValue(float(getattr(self.params.P, "extfit_clip_sigma", 3.0) or 3.0))
        method_form.addRow("Fit clip sigma:", clip_sigma)

        fit_iters = QSpinBox()
        fit_iters.setRange(1, 20)
        fit_iters.setValue(int(getattr(self.params.P, "extfit_fit_iters", 5) or 5))
        method_form.addRow("Fit iterations:", fit_iters)

        min_match = QSpinBox()
        min_match.setRange(3, 200)
        min_match.setValue(int(getattr(self.params.P, "extfit_min_points", 10) or 10))
        method_form.addRow("Min points:", min_match)

        delta_x_enable = QCheckBox("Enable ΔX minimum")
        delta_x_enable.setChecked(bool(getattr(self.params.P, "extinction_delta_x_enable", True)))
        method_form.addRow("ΔX filter:", delta_x_enable)

        delta_x_min = QDoubleSpinBox()
        delta_x_min.setRange(0.0, 2.0)
        delta_x_min.setSingleStep(0.05)
        delta_x_min.setValue(float(getattr(self.params.P, "extinction_delta_x_min", 0.3) or 0.3))
        method_form.addRow("ΔX min:", delta_x_min)

        star_min_frames = QSpinBox()
        star_min_frames.setRange(3, 200)
        star_min_frames.setValue(int(getattr(self.params.P, "extinction_star_min_frames", 8) or 8))
        method_form.addRow("Per-star min frames:", star_min_frames)

        star_rms_max = QDoubleSpinBox()
        star_rms_max.setRange(0.0, 1.0)
        star_rms_max.setSingleStep(0.01)
        star_rms_max.setValue(float(getattr(self.params.P, "extinction_star_rms_max", 0.10) or 0.10))
        method_form.addRow("Per-star RMS max:", star_rms_max)

        star_snr_min = QDoubleSpinBox()
        star_snr_min.setRange(0.0, 200.0)
        star_snr_min.setSingleStep(1.0)
        star_snr_min.setValue(float(getattr(self.params.P, "extinction_star_snr_med_min", 10.0) or 10.0))
        method_form.addRow("Per-star median SNR min:", star_snr_min)

        min_good_stars = QSpinBox()
        min_good_stars.setRange(1, 200)
        min_good_stars.setValue(int(getattr(self.params.P, "extinction_min_good_stars", 3) or 3))
        method_form.addRow("Min good stars:", min_good_stars)

        use_weights = QCheckBox("Use SNR weights (fit)")
        use_weights.setChecked(bool(getattr(self.params.P, "extinction_star_use_weights", True)))
        method_form.addRow("Weights:", use_weights)
        layout.addWidget(method_group)

        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)

        if dialog.exec_() != QDialog.Accepted:
            return

        self.params.P.extinction_varstar_method = "std" if method_combo.currentText() == "STD" else "mad"
        self.params.P.extinction_varstar_sigma = float(var_sigma.value())
        self.params.P.extinction_varstar_min_frames = int(var_min_frames.value())
        self.params.P.extinction_frame_qc_method = "std" if qc_method_combo.currentText() == "STD" else "mad"
        self.params.P.extinction_frame_qc_sigma = float(qc_sigma.value())
        self.params.P.extinction_snr_min = float(snr_cut.value())
        self.params.P.extfit_clip_sigma = float(clip_sigma.value())
        self.params.P.extfit_fit_iters = int(fit_iters.value())
        self.params.P.extfit_min_points = int(min_match.value())
        self.params.P.extinction_delta_x_enable = bool(delta_x_enable.isChecked())
        self.params.P.extinction_delta_x_min = float(delta_x_min.value())
        self.params.P.extinction_star_min_frames = int(star_min_frames.value())
        self.params.P.extinction_star_rms_max = float(star_rms_max.value())
        self.params.P.extinction_star_snr_med_min = float(star_snr_min.value())
        self.params.P.extinction_min_good_stars = int(min_good_stars.value())
        self.params.P.extinction_star_use_weights = bool(use_weights.isChecked())

        saved = False
        try:
            saved = bool(self.params.save_toml())
        except Exception:
            saved = False

        self._rebuild_selection_catalog()
        if saved:
            QMessageBox.information(self, "Saved", "Extinction parameters saved.")
        else:
            QMessageBox.information(self, "Saved", "Parameters updated in memory.")

    def _collect_manual_state_maps(self):
        selected_map = {}
        rejected_map = {}
        for key, state in self.selection_state.items():
            selected = sorted(int(v) for v in state.get("selected", set()))
            rejected = sorted(int(v) for v in state.get("rejected", set()))
            if selected:
                selected_map[key] = selected
            if rejected:
                rejected_map[key] = rejected
        return selected_map, rejected_map

    def run_fit(self):
        if self.worker and self.worker.isRunning():
            return
        source_dir = self._current_source_result_dir()
        if not source_dir.exists():
            QMessageBox.warning(self, "Missing Source Workspace", f"Source workspace not found:\n{source_dir}")
            return
        idx_path = step7_forced_phot_dir(source_dir) / "photometry_index.csv"
        if not idx_path.exists():
            QMessageBox.warning(self, "Missing Step 7", f"Step 7 forced photometry_index.csv not found in:\n{source_dir}")
            return

        self.log_text.clear()
        self.figure.clear()
        self.canvas.draw()
        self.btn_run_fit.setEnabled(False)
        self.btn_run_phot.setEnabled(False)
        self.btn_stop.setEnabled(True)
        self.btn_save.setEnabled(False)
        self.progress_bar.setValue(0)
        self.progress_label.setText("Starting per-star Bouguer fit...")
        self._current_mode = "per_star"

        current_norm = self._normalize_dir_path(source_dir)
        loaded_norm = self._normalize_dir_path(self.loaded_source_result_dir)
        phot_df = self.phot_df if (
            isinstance(self.phot_df, pd.DataFrame)
            and loaded_norm is not None
            and loaded_norm == current_norm
        ) else None
        selected_map, rejected_map = self._collect_manual_state_maps()
        self.log(f"Source workspace: {source_dir}")
        if selected_map or rejected_map:
            self.log(
                f"Manual state applied: selected groups={len(selected_map)}, rejected groups={len(rejected_map)}"
            )

        self.worker = ExtinctionFitWorker(
            self.params,
            self.data_dir,
            self.result_dir,
            mode="per_star",
            task="fit",
            phot_df=phot_df,
            source_result_dir=source_dir,
            selected_star_map=selected_map,
            rejected_star_map=rejected_map,
        )
        self.worker.progress.connect(self.on_progress)
        self.worker.log.connect(self.log)
        self.worker.finished.connect(self.on_fit_finished)
        self.worker.error.connect(self.on_error)
        self.worker.start()

    def stop_fit(self):
        if self.worker and self.worker.isRunning():
            self.worker.stop()
            self.log("Stop requested")

    def on_fit_finished(self, results):
        self.btn_run_fit.setEnabled(True)
        self.btn_run_phot.setEnabled(True)
        self.btn_stop.setEnabled(False)
        if results.get("stopped"):
            self.progress_label.setText("Stopped")
            self.log("Extinction fit stopped")
            return
        self._current_mode = str(results.get("mode", "per_star"))
        self.fit_df = results.get("fit")
        self.points_df = results.get("points")
        if isinstance(self.fit_df, pd.DataFrame) and isinstance(self.points_df, pd.DataFrame):
            self._populate_combos()
            self._plot_results()
            self.btn_save.setEnabled(True)
            self.tabs.setCurrentIndex(2)
        self.log("Extinction fit complete")
        self.progress_label.setText("Done")

    def _populate_combos(self):
        self.date_combo.blockSignals(True)
        self.date_combo.clear()
        self.date_combo.addItem("All Dates")
        if self.fit_df is not None and "date" in self.fit_df.columns:
            for date_val in sorted(set(self.fit_df["date"].dropna().astype(str))):
                self.date_combo.addItem(date_val)
        self.date_combo.blockSignals(False)

        self.filter_combo.blockSignals(True)
        self.filter_combo.clear()
        self.filter_combo.addItem("All Filters")
        if self.fit_df is not None and "filter" in self.fit_df.columns:
            for filt in sorted(set(self.fit_df["filter"].dropna().astype(str))):
                self.filter_combo.addItem(filt)
        self.filter_combo.blockSignals(False)

    def _on_filter_changed(self):
        self._plot_results()

    def on_error(self, message: str):
        self.btn_run_phot.setEnabled(True)
        self.btn_run_fit.setEnabled(True)
        self.btn_stop.setEnabled(False)
        self.progress_label.setText("Error")
        QMessageBox.critical(self, "Error", f"Extinction fit failed:\n{message}")

    def _plot_results(self):
        self.figure.clear()
        ax1, ax2 = self.figure.subplots(2, 1)

        if self.fit_df is None or self.fit_df.empty:
            self.canvas.draw()
            return

        fit_df = self.fit_df.copy()
        points_df = self.points_df.copy() if isinstance(self.points_df, pd.DataFrame) else pd.DataFrame()
        sel_date = self.date_combo.currentText()
        sel_filter = self.filter_combo.currentText()
        show_fit = self.chk_show_fit.isChecked()

        if sel_date != "All Dates" and "date" in fit_df.columns:
            fit_df = fit_df[fit_df["date"] == sel_date]
            if not points_df.empty and "date" in points_df.columns:
                points_df = points_df[points_df["date"] == sel_date]
        if sel_filter != "All Filters" and "filter" in fit_df.columns:
            fit_df = fit_df[fit_df["filter"] == sel_filter]
            if not points_df.empty and "filter" in points_df.columns:
                points_df = points_df[points_df["filter"] == sel_filter]

        colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
        markers = ["o", "s", "^", "v", "D", "p"]
        all_hist = []

        for i, (_, row) in enumerate(fit_df.iterrows()):
            date_val = str(row.get("date", "unknown"))
            filt = str(row["filter"])
            k1 = float(row.get("k1", np.nan))
            color = colors[i % len(colors)]
            marker = markers[i % len(markers)]
            mask = pd.Series(True, index=points_df.index)
            if not points_df.empty:
                if "date" in points_df.columns:
                    mask &= points_df["date"].astype(str) == date_val
                if "filter" in points_df.columns:
                    mask &= points_df["filter"].astype(str) == filt
            sub = points_df[mask].copy() if not points_df.empty else pd.DataFrame()
            label_prefix = f"{date_val}/{filt}" if sel_date == "All Dates" else filt
            if not sub.empty:
                ax1.scatter(
                    sub["airmass"], sub["delta_m"],
                    s=8, alpha=0.35, color=color, marker=marker,
                    label=f"{label_prefix} (n={len(sub)})",
                )
                if show_fit and np.isfinite(k1):
                    xx = np.linspace(float(np.nanmin(sub["airmass"])), float(np.nanmax(sub["airmass"])), 100)
                    yy = k1 * xx
                    ax1.plot(xx, yy, color=color, lw=2, label=f"{label_prefix}: k1={k1:.4f}")
                resid = pd.to_numeric(sub.get("resid"), errors="coerce").dropna()
                if len(resid):
                    all_hist.append((label_prefix, color, resid))

        ax1.set_xlabel("Airmass (X)")
        ax1.set_ylabel("Δm = m_inst - m0_i (mag)")
        title = "Per-star Bouguer Extinction Fit"
        if sel_date != "All Dates":
            title += f" [{sel_date}]"
        if sel_filter != "All Filters":
            title += f" [{sel_filter}]"
        ax1.set_title(title, fontweight="bold")
        if len(fit_df):
            ax1.legend(loc="best", fontsize=7, ncol=min(3, max(len(fit_df), 1)))
        ax1.grid(True, alpha=0.3)

        for label, color, resid in all_hist:
            ax2.hist(resid, bins=30, alpha=0.5, color=color, label=f"{label} (σ={resid.std():.3f})")
        ax2.set_xlabel("Residual (mag)")
        ax2.set_ylabel("Count")
        ax2.set_title("Residual Histogram", fontweight="bold")
        if all_hist:
            ax2.legend(loc="best", fontsize=7)
        ax2.grid(True, alpha=0.3)

        self.figure.tight_layout()
        self.canvas.draw()

    def save_plots(self):
        try:
            out_dir = tool_extinction_dir(self.result_dir)
            out_dir.mkdir(parents=True, exist_ok=True)
            out = out_dir / "per_star_extinction_fit_plot.png"
            self.figure.savefig(out, dpi=150, bbox_inches="tight")
            self.log(f"Saved {out}")

        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to save plot:\n{e}")
