"""
Step 8 - PSF Photometry  (photutils EPSFBuilder + PSFPhotometry, skippable)

Reads:
    step4_detection/detect_{fname}.csv   (det_uid, x, y)
    step7_forced_phot/photometry_{fname}.tsv (initial flux estimate, optional)
    <FITS image>

Writes to cmd_psf/:
    photometry_{fname}.tsv   – det_uid, x_fit, y_fit, mag_psf, mag_psf_err,
                                chi2, iter_found, flags_psf
    epsf_model_{filter}_{frame_stem}.fits – oversampled ePSF model (per-frame)
    residual_iter{N}_{fname}.fits – sky-subtracted residual image (per iteration)
    starsub_iter{N}_{fname}.fits  – raw image with fitted stars removed (per iteration)
    photometry_index.csv     – per-frame summary
    psf_output_signature.json – selected frames, input mtimes, and PSF params

Step can be SKIPPED: clicking "Skip PSF" marks step as complete and
passes control to Step 9 (Master ID Editor). Downstream steps use Step 7 forced aperture
photometry results when PSF outputs are unavailable.
"""
from __future__ import annotations

import json
import hashlib
import traceback
import time
import copy
import threading
from dataclasses import replace
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from threading import Lock

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.nddata import NDData
from astropy.stats import sigma_clipped_stats, mad_std as _mad_std


def _fast_res_std(arr: np.ndarray) -> float:
    """Robust std for residual images: MAD estimator on a 65K-pixel subsample."""
    flat = arr.ravel()
    flat = flat[np.isfinite(flat)]
    if flat.size == 0:
        return 0.0
    stride = max(1, flat.size // 65536)
    return float(_mad_std(flat[::stride]))
from scipy.spatial import cKDTree

from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QGroupBox, QMessageBox,
    QTextEdit, QDialog, QFormLayout, QDialogButtonBox, QProgressBar,
    QCheckBox, QSpinBox, QDoubleSpinBox, QWidget, QTabWidget,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView,
    QSplitter, QComboBox, QListWidget, QScrollArea,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle, Patch

from apex.gui.layout_rules import FittedDialog, prevent_collapse, scroll_wrap, tame_canvas
from apex.gui.workflow.step_window_base import StepWindowBase
from apex.gui.workflow.log_panel import WorkflowLogWindow, WorkerStatusPanel, append_timestamped_log, show_raised
from apex.gui.workflow.ui_helpers import (
    add_parameter_reset_button,
    create_collapsible_section,
    create_output_reuse_checkbox,
    create_parameter_button,
    configure_parameter_dialog,
    set_table_row_background,
    status_row_background,
)
from apex.utils.step_paths_cmd import (
    step2_cropped_dir, step4_dir, step8_psf_dir,
    crop_is_active,
)
from apex.utils.step_paths import step7_forced_phot_dir
from apex.utils.astro_utils import normalize_filter_name
from apex.utils.constants import get_parallel_workers
from apex.utils.noise_params import resolve_effective_noise_params
from apex.utils.qc_utils import filter_files_by_qc, should_use_frame_quality_qc
from apex.utils.psf_core import (
    PSFCoreCut,
    estimate_psf_core_cut,
    psf_core_keep_mask,
    target_pixel_from_wcs,
)
from apex.analysis.psf_policy import (
    estimate_psf_flux_seeds,
    local_group_policy,
    merge_forced_catalog_seeds,
    plan_epsf_stars,
    plan_psf_fit_window,
    psf_symmetric_mask,
    select_epsf_reference_stars,
    select_spatially_balanced,
)
from apex.analysis.psf_iteration import (
    IterationSnapshot,
    PSFFitFlag,
    assess_psf_frame_quality,
    decide_residual_iteration,
    fit_parameters_changed,
    measure_psf_fit_quality,
    qfit_noise_diagnostics,
)
from apex.analysis.psf_diagnostics import (
    draw_psf_final_diagnostics,
    load_psf_final_diagnostic_data,
)
from apex.analysis.psf_flux_scale import (
    PSFApertureScale,
    apply_psf_aperture_scale,
    estimate_psf_aperture_scale,
)


# ── Scalar helpers ────────────────────────────────────────────────────────────

def _to_float(val, default):
    try:
        if val is None:
            return float(default)
        out = float(val)
        return out if np.isfinite(out) else float(default)
    except Exception:
        return float(default)


def _to_int(val, default):
    try:
        if val is None:
            return int(default)
        return int(float(val))
    except Exception:
        return int(default)


def _safe_float(x, default=np.nan):
    try:
        return float(x)
    except Exception:
        return default


def _numeric_series(df: pd.DataFrame, column: str) -> pd.Series:
    if not isinstance(df, pd.DataFrame) or column not in df.columns:
        return pd.Series(dtype=float)
    return pd.to_numeric(df[column], errors="coerce")


def _finite_values(df: pd.DataFrame, column: str) -> np.ndarray:
    values = _numeric_series(df, column).to_numpy(dtype=float)
    return values[np.isfinite(values)]


def _median_value(df: pd.DataFrame, column: str) -> float:
    values = _finite_values(df, column)
    return float(np.median(values)) if values.size else np.nan


def _mean_value(df: pd.DataFrame, column: str) -> float:
    values = _finite_values(df, column)
    return float(np.mean(values)) if values.size else np.nan


def _std_value(df: pd.DataFrame, column: str) -> float:
    values = _finite_values(df, column)
    return float(np.std(values)) if values.size else np.nan


def _first_value(df: pd.DataFrame, column: str, default=np.nan):
    if not isinstance(df, pd.DataFrame) or df.empty or column not in df.columns:
        return default
    values = df[column].dropna()
    return values.iloc[0] if len(values) else default


def _as_bool(value, default: bool = False) -> bool:
    if value is None:
        return bool(default)
    try:
        if pd.isna(value):
            return bool(default)
    except Exception:
        pass
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "no", "n", "off", ""}:
        return False
    return bool(default)


def _filter_key_series(df: pd.DataFrame, column: str) -> pd.Series:
    if not isinstance(df, pd.DataFrame) or column not in df.columns:
        return pd.Series(dtype=str)
    return df[column].map(lambda v: normalize_filter_name(v) if pd.notna(v) else "")


def _string_filter_values(df: pd.DataFrame, column: str) -> set[str]:
    out: set[str] = set()
    for val in _filter_key_series(df, column):
        val = val.strip()
        if val:
            out.add(val)
    return out


def _filter_subset(df: pd.DataFrame, column: str, filt: str | None) -> pd.DataFrame:
    if not isinstance(df, pd.DataFrame) or df.empty:
        return pd.DataFrame(columns=getattr(df, "columns", []))
    if filt is None:
        return df
    if column not in df.columns:
        return df.iloc[0:0].copy()
    filt_key = normalize_filter_name(filt)
    keys = _filter_key_series(df, column)
    return df[keys == filt_key].copy()


def build_ap_psf_comparison(params, result_dir: Path) -> tuple[pd.DataFrame, int]:
    """Merge Step 7 aperture and Step 8 PSF magnitudes on ``det_uid``.

    Reads only from disk, so the window and a headless run share this path.
    Returns ``(merged, n_split_excluded)``.

    Crowd-safe compare: when ``seed_uid`` + ``flux_psf_e`` exist, PSF components
    that were split off one Step-4 seed are summed back to that seed before the
    merge, and (by default) seeds that decomposed into more than one component
    are dropped entirely — an aperture measurement of a blend cannot be
    compared against one of its PSF pieces.
    """
    result_dir = Path(result_dir)
    ap_dir = step7_forced_phot_dir(result_dir)
    psf_dir = step8_psf_dir(result_dir)

    merged_rows: list[pd.DataFrame] = []
    split_excluded_total = 0
    for psf_tsv in sorted(psf_dir.glob("photometry_*.tsv")):
        fname_key = psf_tsv.name[len("photometry_"):]
        ap_tsv = ap_dir / f"photometry_{fname_key}"
        if not ap_tsv.exists():
            continue
        try:
            df_ap = pd.read_csv(ap_tsv, sep="\t")
            df_psf = pd.read_csv(psf_tsv, sep="\t")
        except Exception:
            continue
        if "det_uid" not in df_ap.columns or "det_uid" not in df_psf.columns:
            continue

        # Step 7 writes the aperture magnitude as `mag_inst`; `mag_ap` has never
        # existed in its output. Hard-coding `mag_ap` made every comparison come
        # out empty ("All magnitudes are NaN") — the rest of the codebase already
        # falls back through this list (photometry_loader, step10, extinction_fit).
        ap_mag_col = next(
            (c for c in ("mag_inst", "mag", "mag_ap", "mag_apcorr") if c in df_ap.columns),
            None,
        )
        if ap_mag_col is None:
            continue
        ap_err_col = next(
            (c for c in ("mag_err", "mag_inst_err", "mag_ap_err") if c in df_ap.columns),
            None,
        )
        ap_cols = ["det_uid", ap_mag_col]
        if ap_err_col:
            ap_cols.append(ap_err_col)
        if "r_ap_px" in df_ap.columns:
            ap_cols.append("r_ap_px")
        try:
            if {"seed_uid", "flux_psf_e", "exptime"} <= set(df_psf.columns):
                zp = _to_float(getattr(params.P, "zp_initial", 25.0), 25.0)
                p = df_psf.copy()
                for c in ("seed_uid", "flux_psf_e", "exptime"):
                    p[c] = pd.to_numeric(p[c], errors="coerce")
                p = p[
                    np.isfinite(p["seed_uid"]) & (p["seed_uid"] >= 0)
                    & np.isfinite(p["flux_psf_e"]) & (p["flux_psf_e"] > 0)
                    & np.isfinite(p["exptime"]) & (p["exptime"] > 0)
                ].copy()
                if len(p) == 0:
                    continue
                agg_map = {"flux_psf_e": "sum", "exptime": "median"}
                for c in ("FILTER", "qfit", "qfit_noise_ratio", "iter_found", "snr_psf", "flags_psf"):
                    if c in p.columns:
                        agg_map[c] = "median" if c in {"qfit", "qfit_noise_ratio", "iter_found"} else "first"
                g = p.groupby("seed_uid", as_index=False).agg(agg_map)
                comp = p.groupby("seed_uid", as_index=False).size().rename(columns={"size": "n_comp"})
                g = g.merge(comp, on="seed_uid", how="left")
                if bool(getattr(params.P, "step6_compare_exclude_split", True)):
                    n_before_g = int(len(g))
                    g = g[g["n_comp"] == 1].copy()
                    split_excluded_total += max(0, n_before_g - int(len(g)))
                    if len(g) == 0:
                        continue
                g["det_uid"] = g["seed_uid"].astype(int)
                g["mag_psf"] = zp - 2.5 * np.log10(
                    np.maximum(g["flux_psf_e"].to_numpy(float), 1e-30)
                    / np.maximum(g["exptime"].to_numpy(float), 1e-30)
                )
                psf_cols = [c for c in ("det_uid", "mag_psf", "FILTER", "qfit", "qfit_noise_ratio",
                                        "iter_found", "snr_psf", "flags_psf") if c in g.columns]
                m = df_ap[ap_cols].merge(g[psf_cols], on="det_uid", how="inner")
            else:
                psf_cols = [c for c in ("det_uid", "mag_psf", "mag_psf_err", "FILTER", "qfit",
                                        "qfit_noise_ratio", "iter_found", "snr_psf", "flags_psf")
                            if c in df_psf.columns]
                m = df_ap[ap_cols].merge(df_psf[psf_cols], on="det_uid", how="inner")
            # Downstream (plot + stats) expects `mag_ap`/`mag_ap_err` names.
            rename = {ap_mag_col: "mag_ap"}
            if ap_err_col:
                rename[ap_err_col] = "mag_ap_err"
            m = m.rename(columns=rename)
            # Strip .tsv so FRAME matches apcorr_summary.csv "file" column
            m["FRAME"] = fname_key[:-4] if fname_key.endswith(".tsv") else fname_key
            merged_rows.append(m)
        except Exception:
            continue

    merged = pd.concat(merged_rows, ignore_index=True) if merged_rows else pd.DataFrame()
    return merged, int(split_excluded_total)


def load_psf_qc_inputs(psf_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read the three tables the Step 8 QC products are built from.

    Everything comes off disk (``photometry_index.csv``, the per-frame TSVs and
    ``residual_meta_*.json``), so a headless run can produce the same QC output
    as the window. The GUI calls this too — one code path, not two.

    Returns ``(idx, all_df, meta_df)``; any of them may be empty.
    """
    psf_dir = Path(psf_dir)
    idx_path = psf_dir / "photometry_index.csv"
    idx = pd.read_csv(idx_path) if idx_path.exists() else pd.DataFrame()

    tsv_files = sorted(psf_dir.glob("photometry_*.tsv"))
    all_df = (
        pd.concat([pd.read_csv(f, sep="\t") for f in tsv_files], ignore_index=True)
        if tsv_files else pd.DataFrame()
    )

    meta_rows: list[dict] = []
    for mf in sorted(psf_dir.glob("residual_meta_*.json")):
        try:
            m = json.loads(mf.read_text(encoding="utf-8"))
        except Exception:
            continue
        core = m.get("core_cut", {}) if isinstance(m.get("core_cut", {}), dict) else {}
        base_meta = {
            "file": m.get("file", mf.name.replace("residual_meta_", "", 1).replace(".json", "")),
            "filter": m.get("filter", "?"),
            # Older products predate this field; "?" marks them as unknown
            # rather than silently claiming the current default.
            "psf_build_mode": m.get("psf_build_mode", "?"),
            "psf_fit_engine": m.get("psf_fit_engine", "?"),
            "core_cut_enabled": bool(core.get("enabled", False)),
            "core_cut_x_px": core.get("center_x", np.nan),
            "core_cut_y_px": core.get("center_y", np.nan),
            "core_cut_radius_px": core.get("radius_px", np.nan),
            "core_cut_method": core.get("method", ""),
            "core_cut_reason": core.get("reason", ""),
            "n_core_excluded_init": core.get("n_excluded_init", np.nan),
            "n_core_excluded_redetect": core.get("n_excluded_redetect", np.nan),
            "n_core_excluded_result": core.get("n_excluded_result", np.nan),
        }
        for it in m.get("iters", []):
            meta_rows.append({
                **base_meta,
                "iter": it.get("iter"),
                "phase": it.get("phase", "residual_fit"),
                "residual_std": it.get("residual_std", np.nan),
                "n_fit": it.get("n_fit", np.nan),
                "n_new_raw": it.get("n_new_raw", np.nan),
                "n_new_kept": it.get("n_new_kept", np.nan),
                "n_candidates_raw": it.get("n_candidates_raw", np.nan),
                "n_candidates_unique": it.get("n_candidates_unique", np.nan),
                "n_candidates_accepted": it.get("n_candidates_accepted", np.nan),
                "n_pruned": it.get("n_pruned", np.nan),
                "median_qfit": it.get("median_qfit", np.nan),
                "median_reduced_chi2": it.get("median_reduced_chi2", np.nan),
                "elapsed_s": it.get("elapsed_s", np.nan),
                "stop_reason": it.get("stop_reason", ""),
            })
    meta_df = pd.DataFrame(meta_rows) if meta_rows else pd.DataFrame()
    return idx, all_df, meta_df


def render_psf_final_diagnostics(
    fig,
    params,
    result_dir: Path,
    fname: str,
    *,
    use_cropped: bool = False,
) -> tuple[pd.DataFrame, dict]:
    """Draw the Step 8 final-diagnostic panel for one frame onto ``fig``.

    Assembles every input off disk (residual meta, ePSF model, reference-star
    catalogue, frame FWHM, pixel scale) so the window and a headless run render
    the same figure. Returns ``(data, summary)``.
    """
    result_dir = Path(result_dir)
    psf_dir = step8_psf_dir(result_dir)
    data = load_psf_final_diagnostic_data(result_dir, fname)

    meta_path = psf_dir / f"residual_meta_{fname}.json"
    meta: dict = {}
    if meta_path.exists():
        try:
            loaded = json.loads(meta_path.read_text(encoding="utf-8"))
            meta = loaded if isinstance(loaded, dict) else {}
        except Exception:
            meta = {}
    core = meta.get("core_cut", {}) if isinstance(meta.get("core_cut", {}), dict) else {}

    cache_dir = Path(getattr(params.P, "cache_dir", result_dir))
    fwhm_px = _load_fwhm_from_meta(
        fname, cache_dir, result_dir,
        _to_float(getattr(params.P, "fwhm_pix_guess", 6.0), 6.0),
    )

    # pixel scale — from the frame WCS when it can be found
    pixel_scale = np.nan
    fits_path = None
    if use_cropped and crop_is_active(result_dir):
        cpath = step2_cropped_dir(result_dir) / fname
        if cpath.exists():
            fits_path = cpath
    if fits_path is None:
        cand = Path(getattr(params.P, "data_dir", "")) / fname
        fits_path = cand if cand.exists() else None
    if fits_path is not None:
        try:
            from astropy.wcs import WCS
            from astropy.wcs.utils import proj_plane_pixel_scales
            celestial = WCS(fits.getheader(fits_path)).celestial
            scales = np.asarray(proj_plane_pixel_scales(celestial), dtype=float) * 3600.0
            s = float(np.nanmedian(np.abs(scales)))
            if np.isfinite(s) and s > 0:
                pixel_scale = s
        except Exception:
            pass
    if not np.isfinite(pixel_scale):
        pixel_scale = _to_float(getattr(params.P, "pixel_scale_arcsec", np.nan), np.nan)

    # ePSF model
    filt = str(meta.get("filter", "")).strip()
    stem = Path(fname).stem
    epsf_model, epsf_path = None, None
    cands: list[Path] = []
    if filt:
        for f in (filt, filt.lower(), filt.upper()):
            cands += [psf_dir / f"epsf_model_{f}_{stem}.fits", psf_dir / f"epsf_model_{f}.fits"]
    cands += sorted(psf_dir.glob(f"epsf_model_*_{stem}.fits"))
    for c in cands:
        if c.exists():
            try:
                epsf_model = fits.getdata(c).astype(float)
                epsf_path = c
                break
            except Exception:
                continue

    # ePSF reference-star catalogue
    reference = meta.get("epsf_reference", {})
    cat_name = reference.get("catalog_path", "") if isinstance(reference, dict) else ""
    cat_path = psf_dir / (str(cat_name).strip() or f"epsf_reference_{fname}.csv")
    epsf_reference = pd.DataFrame()
    if cat_path.exists():
        try:
            epsf_reference = pd.read_csv(cat_path)
        except Exception:
            pass

    summary = draw_psf_final_diagnostics(
        fig, data, epsf_model,
        filename=fname,
        fwhm_px=fwhm_px,
        pixel_scale_arcsec=pixel_scale,
        core_center=(_safe_float(core.get("center_x"), np.nan),
                     _safe_float(core.get("center_y"), np.nan)),
        core_radius_px=_safe_float(core.get("radius_px"), np.nan),
        epsf_reference=epsf_reference,
    )
    if epsf_path is not None:
        summary["epsf_file"] = epsf_path.name
    flux_scale = meta.get("flux_scale", {})
    if isinstance(flux_scale, dict):
        summary["psf_aperture_scale"] = _safe_float(flux_scale.get("scale"), 1.0)
        summary["psf_aperture_scale_applied"] = bool(flux_scale.get("applied", False))
        summary["psf_aperture_scale_n"] = _to_int(flux_scale.get("n_used", 0), 0)
        summary["psf_aperture_scale_raw_offset_mag"] = _safe_float(
            flux_scale.get("median_delta_mag_raw"), np.nan
        )
    return data, summary


def export_psf_qc_products(psf_dir: Path, params=None, result_dir: Path | None = None) -> list[Path]:
    """Write the window-independent Step 8 QC products.

    Covers the parts that need no Qt widget: the two QC tables, the
    residual/core overview figure, and — when ``params``/``result_dir`` are
    given — the aperture-vs-PSF comparison table. The window additionally
    renders its own interactive plots from live widget state.
    """
    psf_dir = Path(psf_dir)
    idx, all_df, meta_df = load_psf_qc_inputs(psf_dir)
    if idx.empty and all_df.empty:
        return []

    saved: list[Path] = []

    cmp_df = None
    if params is not None and result_dir is not None:
        try:
            cmp_df, _n_split = build_ap_psf_comparison(params, result_dir)
            if not cmp_df.empty:
                p = psf_dir / "psf_ap_vs_psf.csv"
                cmp_df.to_csv(p, index=False)
                saved.append(p)
                # 창이 이 표를 캐시로 재사용한다(_load_or_build_comparison).
                # 병합만 10초대라 헤드리스가 남겨 두면 창이 즉시 열린다.
                (psf_dir / "psf_ap_vs_psf_meta.json").write_text(
                    json.dumps({"split_excluded_total": int(_n_split)}),
                    encoding="utf-8",
                )
        except Exception:
            cmp_df = None

    summary = _build_psf_qc_summary(idx, all_df, meta_df, cmp_df)
    if not summary.empty:
        p = psf_dir / "psf_qc_summary.csv"
        summary.to_csv(p, index=False)
        saved.append(p)

    frame_qc = _build_psf_frame_qc_table(idx, meta_df)
    if not frame_qc.empty:
        p = psf_dir / "psf_frame_qc.csv"
        frame_qc.to_csv(p, index=False)
        saved.append(p)

        fig = Figure(figsize=(10.5, 6.8), dpi=120)
        if _draw_psf_frame_qc_overview(fig, frame_qc):
            fp = psf_dir / "step8_residual_core_qc.png"
            fig.savefig(fp, dpi=160, bbox_inches="tight")
            saved.append(fp)

    # 최종 진단 — 대표 프레임 한 장. 창은 사용자가 고른 프레임을 그리지만
    # 배치는 고를 사람이 없으므로 인덱스 첫 프레임을 쓴다.
    if params is not None and result_dir is not None and not idx.empty:
        try:
            fname = Path(str(idx["file"].iloc[0])).name
            fig = Figure(figsize=(12.0, 7.5), dpi=120)
            data, summary = render_psf_final_diagnostics(
                fig, params, Path(result_dir), fname
            )
            stem = Path(fname).stem
            fp = psf_dir / f"step8_final_diagnostics_{stem}.png"
            fig.savefig(fp, dpi=160, bbox_inches="tight")
            saved.append(fp)

            sp = psf_dir / f"psf_final_diagnostics_{stem}.json"
            sp.write_text(
                json.dumps(summary, ensure_ascii=False, indent=2,
                           allow_nan=True, default=str),
                encoding="utf-8",
            )
            saved.append(sp)
        except Exception:
            pass
    return saved


def _build_psf_qc_summary(
    idx: pd.DataFrame,
    phot_df: pd.DataFrame,
    meta_df: pd.DataFrame | None = None,
    cmp_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build the compact Step 8 QC table exported beside PSF products."""
    idx = idx.copy() if isinstance(idx, pd.DataFrame) else pd.DataFrame()
    phot_df = phot_df.copy() if isinstance(phot_df, pd.DataFrame) else pd.DataFrame()
    meta_df = meta_df.copy() if isinstance(meta_df, pd.DataFrame) else pd.DataFrame()
    cmp_df = cmp_df.copy() if isinstance(cmp_df, pd.DataFrame) else pd.DataFrame()

    filters = _string_filter_values(idx, "filter") | _string_filter_values(phot_df, "FILTER")
    filters |= _string_filter_values(meta_df, "filter") | _string_filter_values(cmp_df, "FILTER")
    groups: list[tuple[str, str | None]] = [("ALL", None)]
    groups.extend((filt, filt) for filt in sorted(filters))

    if not cmp_df.empty and {"mag_ap", "mag_psf"} <= set(cmp_df.columns):
        cmp_df["mag_ap"] = pd.to_numeric(cmp_df["mag_ap"], errors="coerce")
        cmp_df["mag_psf"] = pd.to_numeric(cmp_df["mag_psf"], errors="coerce")
        cmp_df = cmp_df[np.isfinite(cmp_df["mag_ap"]) & np.isfinite(cmp_df["mag_psf"])].copy()
        cmp_df["delta_ap_minus_psf"] = cmp_df["mag_ap"] - cmp_df["mag_psf"]
    else:
        cmp_df = pd.DataFrame()

    rows = []
    for label, filt in groups:
        frame_sub = _filter_subset(idx, "filter", filt)
        phot_sub = _filter_subset(phot_df, "FILTER", filt)
        meta_sub = _filter_subset(meta_df, "filter", filt)
        cmp_sub = _filter_subset(cmp_df, "FILTER", filt) if not cmp_df.empty else pd.DataFrame()

        if "flags_psf" in phot_sub.columns:
            flags = _numeric_series(phot_sub, "flags_psf")
            clean_mask = np.isfinite(flags.to_numpy(dtype=float)) & (flags.to_numpy(dtype=float) == 0)
            good = phot_sub.loc[clean_mask].copy()
        else:
            good = phot_sub.copy()

        qfit = _finite_values(good, "qfit")
        qfit_gt5_fraction = float(np.mean(qfit > 5.0)) if qfit.size else np.nan
        qfit_gt1_fraction = float(np.mean(qfit > 1.0)) if qfit.size else np.nan
        qfit_noise = _finite_values(good, "qfit_noise_ratio")
        qfit_noise_gt3_fraction = (
            float(np.mean(qfit_noise > 3.0)) if qfit_noise.size else np.nan
        )
        cfit = _finite_values(good, "cfit")
        cfit_abs_gt01_fraction = float(np.mean(np.abs(cfit) > 0.1)) if cfit.size else np.nan
        forced_mask = phot_sub.get(
            "forced_psf", pd.Series(False, index=phot_sub.index)
        ).map(_as_bool).to_numpy(dtype=bool)
        flux_values = pd.to_numeric(
            phot_sub.get("flux_psf_e", pd.Series(np.nan, index=phot_sub.index)),
            errors="coerce",
        ).to_numpy(dtype=float)
        all_flags = pd.to_numeric(
            phot_sub.get("flags_psf", pd.Series(0, index=phot_sub.index)),
            errors="coerce",
        ).fillna(0).to_numpy(dtype=np.int64)

        if not meta_sub.empty and {"iter", "residual_std"} <= set(meta_sub.columns):
            meta_sub["iter"] = pd.to_numeric(meta_sub["iter"], errors="coerce")
            i1 = meta_sub[meta_sub["iter"] == 1]
            i2 = meta_sub[meta_sub["iter"] == 2]
            residual_i1 = _mean_value(i1, "residual_std")
            residual_i2 = _mean_value(i2, "residual_std")
        else:
            residual_i1 = np.nan
            residual_i2 = np.nan

        rows.append({
            "filter": label,
            "n_frames": int(len(frame_sub)),
            "n_psf_rows": int(len(phot_sub)),
            "n_clean": int(len(good)),
            "clean_fraction": float(len(good) / len(phot_sub)) if len(phot_sub) else np.nan,
            "median_n_psf_per_frame": _median_value(frame_sub, "n"),
            "median_n_goodmag_per_frame": _median_value(frame_sub, "n_goodmag"),
            "median_n_fail_per_frame": _median_value(frame_sub, "n_fail"),
            "median_n_new_iter_per_frame": _median_value(frame_sub, "n_new_iter"),
            "n_forced": int(np.sum(forced_mask)),
            "n_forced_negative": int(np.sum(
                forced_mask & np.isfinite(flux_values) & (flux_values <= 0)
            )),
            "n_crowding_unreliable": int(np.sum(
                (all_flags & int(PSFFitFlag.CROWDING_UNRELIABLE)) != 0
            )),
            "median_mag_psf": _median_value(good, "mag_psf"),
            "median_mag_psf_err": _median_value(good, "mag_psf_err"),
            "median_snr_psf": _median_value(good, "snr_psf"),
            "median_qfit": _median_value(good, "qfit"),
            "median_qfit_noise_ratio": _median_value(good, "qfit_noise_ratio"),
            "qfit_noise_gt3_fraction": qfit_noise_gt3_fraction,
            "qfit_gt5_fraction": qfit_gt5_fraction,
            "qfit_gt1_fraction": qfit_gt1_fraction,
            "median_cfit": _median_value(good, "cfit"),
            "cfit_abs_gt01_fraction": cfit_abs_gt01_fraction,
            "median_reduced_chi2": _median_value(good, "reduced_chi2"),
            "n_ap_psf_matches": int(len(cmp_sub)),
            "median_ap_minus_psf": _median_value(cmp_sub, "delta_ap_minus_psf"),
            "std_ap_minus_psf": _std_value(cmp_sub, "delta_ap_minus_psf"),
            "residual_std_iter1_mean": residual_i1,
            "residual_std_iter2_mean": residual_i2,
        })

    return pd.DataFrame(rows)


def _build_psf_frame_qc_table(idx: pd.DataFrame, meta_df: pd.DataFrame | None = None) -> pd.DataFrame:
    """Build one Step 8 QC row per frame from photometry_index and residual metadata."""
    idx = idx.copy() if isinstance(idx, pd.DataFrame) else pd.DataFrame()
    meta_df = meta_df.copy() if isinstance(meta_df, pd.DataFrame) else pd.DataFrame()
    if idx.empty and meta_df.empty:
        return pd.DataFrame()

    files = set()
    if "file" in idx.columns:
        files |= {str(v) for v in idx["file"].dropna().astype(str)}
    if "file" in meta_df.columns:
        files |= {str(v) for v in meta_df["file"].dropna().astype(str)}

    rows = []
    for fname in sorted(files):
        frame_sub = idx[idx["file"].astype(str) == fname].copy() if "file" in idx.columns else pd.DataFrame()
        meta_sub = meta_df[meta_df["file"].astype(str) == fname].copy() if "file" in meta_df.columns else pd.DataFrame()

        filt = _first_value(frame_sub, "filter", _first_value(meta_sub, "filter", ""))
        filt = normalize_filter_name(filt)
        # Which PSF model produced this frame's magnitudes. The per-frame QC
        # table is where a reader compares frames, so a mixed-mode result set
        # has to be visible here rather than only inside each residual JSON.
        build_mode = _first_value(meta_sub, "psf_build_mode", "?")
        fit_engine = _first_value(meta_sub, "psf_fit_engine", "?")

        if "iter" in meta_sub.columns:
            meta_sub["iter"] = pd.to_numeric(meta_sub["iter"], errors="coerce")
            i1 = meta_sub[meta_sub["iter"] == 1]
            i2 = meta_sub[meta_sub["iter"] == 2]
            finite_iter = meta_sub[np.isfinite(meta_sub["iter"].to_numpy(dtype=float))]
            if not finite_iter.empty:
                final_iter = int(np.nanmax(finite_iter["iter"].to_numpy(dtype=float)))
                ifinal = meta_sub[meta_sub["iter"] == final_iter]
            else:
                final_iter = np.nan
                ifinal = pd.DataFrame()
        else:
            i1 = pd.DataFrame()
            i2 = pd.DataFrame()
            ifinal = pd.DataFrame()
            final_iter = np.nan

        r1 = _mean_value(i1, "residual_std")
        r2 = _mean_value(i2, "residual_std")
        rfinal = _mean_value(ifinal, "residual_std")
        if np.isfinite(r1) and np.isfinite(rfinal):
            residual_delta = float(rfinal - r1)
            residual_frac = float((rfinal - r1) / r1) if r1 != 0 else np.nan
        else:
            residual_delta = np.nan
            residual_frac = np.nan

        core_enabled = _as_bool(
            _first_value(meta_sub, "core_cut_enabled", _first_value(frame_sub, "core_cut_enabled", False))
        )

        rows.append({
            "file": fname,
            "filter": filt,
            "psf_build_mode": build_mode,
            "psf_fit_engine": fit_engine,
            "frame_fwhm_px": _first_value(frame_sub, "frame_fwhm_px", np.nan),
            "frame_fwhm_arcsec": _first_value(frame_sub, "frame_fwhm_arcsec", np.nan),
            "psf_qc_status": _first_value(frame_sub, "psf_qc_status", ""),
            "psf_qc_score": _first_value(frame_sub, "psf_qc_score", np.nan),
            "psf_qc_reasons": _first_value(frame_sub, "psf_qc_reasons", ""),
            "psf_clean_fraction": _first_value(frame_sub, "psf_clean_fraction", np.nan),
            "psf_fit_failure_fraction": _first_value(
                frame_sub, "psf_fit_failure_fraction", np.nan
            ),
            "psf_crowding_unreliable_fraction": _first_value(
                frame_sub, "psf_crowding_unreliable_fraction", np.nan
            ),
            "frame_total_elapsed_s": _first_value(
                frame_sub, "frame_total_elapsed_s", np.nan
            ),
            "fit_elapsed_s": _first_value(frame_sub, "fit_elapsed_s", np.nan),
            "epsf_elapsed_s": _first_value(frame_sub, "epsf_elapsed_s", np.nan),
            "n_psf": _first_value(frame_sub, "n", np.nan),
            "n_goodmag": _first_value(frame_sub, "n_goodmag", np.nan),
            "n_fail": _first_value(frame_sub, "n_fail", np.nan),
            "n_new_iter": _first_value(frame_sub, "n_new_iter", np.nan),
            "n_forced": _first_value(frame_sub, "n_forced", np.nan),
            "n_forced_negative": _first_value(frame_sub, "n_forced_negative", np.nan),
            "n_crowding_unreliable": _first_value(frame_sub, "n_crowding_unreliable", np.nan),
            "median_qfit_noise_ratio": _first_value(
                frame_sub, "median_qfit_noise_ratio", np.nan
            ),
            "fit_window_mode": _first_value(frame_sub, "fit_window_mode", ""),
            "fit_window_px": _first_value(frame_sub, "fit_window_px", np.nan),
            "fit_window_energy": _first_value(
                frame_sub, "fit_window_energy", np.nan
            ),
            "psf_nea_px": _first_value(frame_sub, "psf_nea_px", np.nan),
            "core_cut_enabled": core_enabled,
            "core_cut_method": _first_value(meta_sub, "core_cut_method", ""),
            "core_cut_reason": _first_value(meta_sub, "core_cut_reason", ""),
            "core_cut_x_px": _first_value(meta_sub, "core_cut_x_px", _first_value(frame_sub, "core_cut_x_px", np.nan)),
            "core_cut_y_px": _first_value(meta_sub, "core_cut_y_px", _first_value(frame_sub, "core_cut_y_px", np.nan)),
            "core_cut_radius_px": _first_value(
                meta_sub, "core_cut_radius_px", _first_value(frame_sub, "core_cut_radius_px", np.nan)
            ),
            "n_core_excluded_init": _first_value(
                meta_sub, "n_core_excluded_init", _first_value(frame_sub, "n_core_excluded_init", np.nan)
            ),
            "n_core_excluded_redetect": _first_value(
                meta_sub, "n_core_excluded_redetect", _first_value(frame_sub, "n_core_excluded_redetect", np.nan)
            ),
            "n_core_excluded_result": _first_value(
                meta_sub, "n_core_excluded_result", _first_value(frame_sub, "n_core_excluded_result", np.nan)
            ),
            "final_iter": final_iter,
            "residual_std_iter1": r1,
            "residual_std_iter2": r2,
            "residual_std_final": rfinal,
            "residual_std_final_minus_iter1": residual_delta,
            "residual_std_frac_change": residual_frac,
            "n_fit_iter1": _first_value(i1, "n_fit", np.nan),
            "n_fit_iter2": _first_value(i2, "n_fit", np.nan),
            "n_new_raw_iter2": _first_value(i2, "n_new_raw", np.nan),
            "n_new_kept_iter2": _first_value(i2, "n_new_kept", np.nan),
        })

    return pd.DataFrame(rows)


def _draw_psf_frame_qc_overview(fig: Figure, frame_qc: pd.DataFrame) -> bool:
    if not isinstance(frame_qc, pd.DataFrame) or frame_qc.empty:
        return False
    fig.clear()
    ax_res = fig.add_subplot(211)
    ax_core = fig.add_subplot(212)

    work = frame_qc.copy().reset_index(drop=True)
    x = np.arange(len(work), dtype=float)
    labels = work["file"].astype(str).tolist() if "file" in work.columns else [str(i) for i in range(len(work))]

    r1 = _numeric_series(work, "residual_std_iter1").to_numpy(dtype=float)
    rf = _numeric_series(work, "residual_std_final").to_numpy(dtype=float)
    ok1 = np.isfinite(r1)
    okf = np.isfinite(rf)
    if np.any(ok1):
        ax_res.plot(x[ok1], r1[ok1], "o-", ms=3.0, lw=0.9, color="#4C78A8", label="iter1")
    if np.any(okf):
        ax_res.plot(x[okf], rf[okf], "o-", ms=3.0, lw=0.9, color="#F58518", label="final")
    ax_res.set_ylabel("residual std (ADU)")
    ax_res.set_title("Step 8 PSF Residual QC")
    ax_res.grid(True, alpha=0.2)
    if ax_res.get_legend_handles_labels()[0]:
        ax_res.legend(loc="best", fontsize=8, frameon=False)

    excluded = _numeric_series(work, "n_core_excluded_result").fillna(0.0).to_numpy(dtype=float)
    crowding = _numeric_series(work, "n_crowding_unreliable").fillna(0.0).to_numpy(dtype=float)
    negative = _numeric_series(work, "n_forced_negative").fillna(0.0).to_numpy(dtype=float)
    enabled = work.get("core_cut_enabled", pd.Series([False] * len(work))).map(_as_bool).to_numpy(dtype=bool)
    width = 0.26
    ax_core.bar(x - width, crowding, color="#E69F00", width=width, label="unresolved")
    ax_core.bar(x, negative, color="#0072B2", width=width, label="forced flux <= 0")
    ax_core.bar(x + width, excluded, color="#CC79A7", width=width, label="hard-core excluded")
    ax_core.set_ylabel("flagged / excluded sources")
    ax_core.set_xlabel("frame")
    ax_core.grid(True, axis="y", alpha=0.2)
    ax_core.set_xlim(-0.6, max(0.6, len(work) - 0.4))
    ax_core.legend(loc="best", fontsize=8, frameon=False)

    if len(labels) <= 24:
        ax_core.set_xticks(x)
        ax_core.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    else:
        ax_core.set_xticks([])

    n_enabled = int(np.sum(enabled))
    fig.suptitle(f"Step 8 Frame QC | frames={len(work)} hard_core_cut={n_enabled}", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return True


def _odd_int(value: float, min_value: int = 3, max_value: int | None = None) -> int:
    """Convert to odd integer within optional bounds."""
    try:
        v = int(round(float(value)))
    except Exception:
        v = int(min_value)
    v = max(int(min_value), v)
    if max_value is not None:
        v = min(int(max_value), v)
    if v % 2 == 0:
        v += 1
    if max_value is not None and v > int(max_value):
        v = int(max_value) - 1
        if v < int(min_value):
            v = int(min_value)
        if v % 2 == 0:
            v = max(int(min_value), v - 1)
    return int(v)


# Korean → ASCII translation for matplotlib (matplotlib default font lacks Korean glyphs)
_KO_TO_ASCII = {
    "신규검출 (step4 미검출)": "New (not in step4)",
    "재검출 (step4 기검출)": "Re-detected (step4)",
    "경계소스": "Edge",
}

_PSF_SIGNATURE_FILE = "psf_output_signature.json"
_PSF_SIGNATURE_VERSION = 3

_PSF_SIGNATURE_PARAMS = (
    "phot_use_qc_pass_only",
    "gain_e_per_adu",
    "zp_initial",
    "rdnoise_e",
    "noise_use_fits_header",
    "noise_reference_binning",
    "noise_scale_by_binning",
    "saturation_adu",
    "min_snr_for_mag",
    "fwhm_pix_guess",
    "psf_mode",
    "psf_model_mode",
    "psf_fit_engine",
    "psf_build_mode",
    "psf_parallel_workers",
    "psf_epsf_oversampling",
    "psf_epsf_maxiters",
    "psf_epsf_size_px",
    "psf_epsf_size_fwhm_mult",
    "psf_n_stars_max",
    "psf_isolation_fwhm_mult",
    "psf_epsf_contamination_filter",
    "psf_flux_scale_correction",
    "psf_flux_scale_min_snr",
    "psf_flux_scale_min_stars",
    "psf_flux_scale_min_neighbor_fwhm",
    "psf_flux_scale_max_scatter_mag",
    "psf_flux_percentile_lo",
    "psf_flux_percentile_hi",
    "psf_fit_shape_px",
    "psf_fit_shape_fwhm_mult",
    "psf_fit_window_mode",
    "psf_fit_encircled_energy",
    "psf_max_iter",
    "psf_fitter_max_iter",
    "psf_fit_mode",
    "psf_redetect_sigma",
    "psf_redetect_sigma_g",
    "psf_redetect_sigma_r",
    "psf_redetect_sigma_i",
    "psf_epsf_sharp_lo",
    "psf_epsf_sharp_hi",
    "psf_epsf_round_abs_max",
    "psf_epsf_elong_max",
    "psf_duplicate_radius_fwhm_mult",
    "psf_duplicate_radius_px",
    "psf_new_sources_cap_per_iter",
    "psf_new_sources_cap_frac",
    "psf_fit_init_max_sources",
    "psf_core_cut_enable",
    "psf_core_cut_center_mode",
    "psf_core_cut_x_px",
    "psf_core_cut_y_px",
    "psf_core_cut_radius_px",
    "psf_core_cut_radius_fwhm_mult",
    "psf_core_cut_auto_cell_fwhm_mult",
    "psf_core_cut_auto_min_density_ratio",
    "psf_core_cut_auto_min_sources",
    "psf_core_cut_max_exclude_frac",
    "psf_substar_iters",
    "psf_substar_neighbor_r_fwhm_mult",
    "psf_substar_max_sources",
    "psf_conv_new_frac",
    "psf_postfit_snr_min",
    "psf_postfit_qfit_max",
    "psf_postfit_reduced_chi2_max",
    "psf_blend_residual_ratio",
    "psf_flux_conv_threshold",
    "psf_use_grouper",
    "psf_grouper_max_size",
    "psf_grouper_radius_fwhm",
    "psf_forced_match_radius_fwhm",
    "psf_use_error_image",
    "psf_shared_filter_epsf",
    "psf_min_epsf_stars",
    "psf_redetect_sharp_lo",
    "psf_redetect_sharp_hi",
    "psf_redetect_round_abs_max",
    "psf_save_residuals",
    "psf_save_all_iter_residuals",
)


def _psf_signature_value(value):
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (bool, int, str)) or value is None:
        return value
    if isinstance(value, (list, tuple, set)):
        return [_psf_signature_value(item) for item in value]
    if isinstance(value, dict):
        return {
            str(key): _psf_signature_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass
    return str(value)


def _psf_file_signature(path: Path | None) -> dict | None:
    if path is None:
        return None
    try:
        candidate = Path(path)
        if not candidate.exists():
            return None
        stat = candidate.stat()
        try:
            path_text = str(candidate.resolve())
        except Exception:
            path_text = str(candidate)
        return {
            "path": path_text,
            "size": int(stat.st_size),
            "mtime_ns": int(stat.st_mtime_ns),
        }
    except Exception:
        return None


def _first_psf_input(candidates: list[Path], *, newest: bool = False) -> Path | None:
    existing = []
    for candidate in candidates:
        try:
            if candidate.exists() and candidate.stat().st_size > 0:
                existing.append(candidate)
        except Exception:
            continue
    if not existing:
        return None
    if newest:
        return max(existing, key=lambda path: path.stat().st_mtime_ns)
    return existing[0]


def build_psf_output_signature(
    params,
    frames: list[str],
    *,
    use_cropped: bool = False,
    cache_dir: Path | None = None,
) -> dict:
    """Build the Step 8 completion signature for GUI and headless runs."""
    P = params.P
    result_dir = Path(P.result_dir)
    cache_path = Path(cache_dir or getattr(P, "cache_dir", result_dir / "cache"))
    if not cache_path.is_absolute():
        cache_path = result_dir / cache_path
    step4_path = step4_dir(result_dir)
    step7_path = step7_forced_phot_dir(result_dir)
    cropped_dir = step2_cropped_dir(result_dir)
    data_dir = Path(P.data_dir)

    frame_inputs = []
    for frame in frames:
        detect_csv = _first_psf_input([
            cache_path / f"detect_{frame}.csv",
            step4_path / f"detect_{frame}.csv",
        ])
        detect_json = _first_psf_input(
            [
                cache_path / f"detect_{frame}.json",
                step4_path / f"detect_{frame}.json",
            ],
            newest=True,
        )
        fits_path = cropped_dir / frame if use_cropped else data_dir / frame
        frame_inputs.append({
            "file": str(frame),
            "fits": _psf_file_signature(fits_path),
            "detect_csv": _psf_file_signature(detect_csv),
            "detect_json": _psf_file_signature(detect_json),
            "step7_tsv": _psf_file_signature(step7_path / f"photometry_{frame}.tsv"),
        })

    payload = {
        "signature_version": _PSF_SIGNATURE_VERSION,
        "step": "cmd_step8_psf_photometry",
        "frames": [str(frame) for frame in frames],
        "use_cropped": bool(use_cropped),
        "params": {
            key: _psf_signature_value(getattr(P, key, None))
            for key in _PSF_SIGNATURE_PARAMS
        },
        "inputs": {
            "step7_index": _psf_file_signature(step7_path / "photometry_index.csv"),
            "step7_apcorr_summary": _psf_file_signature(step7_path / "apcorr_summary.csv"),
            "frames": frame_inputs,
        },
    }
    encoded = json.dumps(payload, sort_keys=True, ensure_ascii=False, allow_nan=False)
    payload["signature_hash"] = hashlib.sha1(encoded.encode("utf-8")).hexdigest()
    return payload


def write_psf_output_signature(result_dir: Path | str, signature: dict) -> Path:
    output_dir = step8_psf_dir(result_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    signature_path = output_dir / _PSF_SIGNATURE_FILE
    signature_path.write_text(
        json.dumps(signature, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False),
        encoding="utf-8",
    )
    return signature_path


_PSF_MODE_PRESETS = {
    "normal": {
        "psf_n_stars_max": 0,
        "psf_isolation_fwhm_mult": 3.0,
        "psf_epsf_contamination_filter": True,
        "psf_flux_scale_correction": False,
        "psf_fit_shape_fwhm_mult": 2.4,
        "psf_fit_window_mode": "auto",
        "psf_fit_encircled_energy": 0.90,
        "psf_max_iter": 2,
        "psf_fitter_max_iter": 6,
        "psf_redetect_sigma": 4.0,
        "psf_duplicate_radius_fwhm_mult": 0.8,
        "psf_new_sources_cap_per_iter": 70,
        "psf_new_sources_cap_frac": 0.02,
        "psf_fit_init_max_sources": 0,
        "psf_core_cut_enable": False,
        "psf_core_cut_radius_px": 0.0,
        "psf_core_cut_radius_fwhm_mult": 20.0,
        "psf_core_cut_auto_min_density_ratio": 1.5,
        "psf_substar_neighbor_r_fwhm_mult": 8.0,
        "psf_substar_max_sources": 1500,
        "psf_substar_iters": 1,
        "psf_conv_new_frac": 0.02,
        "psf_postfit_snr_min": 3.0,
        "psf_postfit_qfit_max": 3.0,
        "psf_postfit_reduced_chi2_max": 25.0,
        "psf_blend_residual_ratio": 0.3,
        "psf_flux_conv_threshold": 0.01,
        "psf_use_grouper": False,
        "psf_grouper_radius_fwhm": 1.5,
        "psf_forced_match_radius_fwhm": 1.25,
        "psf_redetect_sharp_lo": 0.15,
        "psf_redetect_sharp_hi": 0.95,
        "psf_redetect_round_abs_max": 0.8,
    },
    "crowded": {
        "psf_n_stars_max": 0,
        "psf_isolation_fwhm_mult": 2.0,
        "psf_epsf_contamination_filter": True,
        "psf_flux_scale_correction": False,
        "psf_fit_shape_fwhm_mult": 2.4,
        "psf_fit_window_mode": "auto",
        "psf_fit_encircled_energy": 0.90,
        "psf_max_iter": 2,
        "psf_fitter_max_iter": 8,
        "psf_redetect_sigma": 4.5,
        "psf_duplicate_radius_fwhm_mult": 0.4,
        "psf_new_sources_cap_per_iter": 50,
        "psf_new_sources_cap_frac": 0.015,
        "psf_fit_init_max_sources": 3000,
        "psf_core_cut_enable": False,
        "psf_core_cut_radius_px": 0.0,
        "psf_core_cut_radius_fwhm_mult": 20.0,
        "psf_core_cut_auto_min_density_ratio": 1.5,
        "psf_substar_neighbor_r_fwhm_mult": 5.0,
        "psf_substar_max_sources": 1000,
        "psf_substar_iters": 1,
        "psf_conv_new_frac": 0.02,
        "psf_postfit_snr_min": 3.0,
        "psf_postfit_qfit_max": 3.0,
        "psf_postfit_reduced_chi2_max": 25.0,
        "psf_blend_residual_ratio": 0.3,
        "psf_flux_conv_threshold": 0.01,
        "psf_use_grouper": False,
        "psf_grouper_radius_fwhm": 1.5,
        "psf_forced_match_radius_fwhm": 1.25,
        "psf_redetect_sharp_lo": 0.2,
        "psf_redetect_sharp_hi": 0.9,
        "psf_redetect_round_abs_max": 0.6,
    },
    "faint": {
        "psf_n_stars_max": 0,
        "psf_isolation_fwhm_mult": 2.5,
        "psf_epsf_contamination_filter": True,
        "psf_flux_scale_correction": False,
        "psf_fit_shape_fwhm_mult": 2.4,
        "psf_fit_window_mode": "auto",
        "psf_fit_encircled_energy": 0.90,
        "psf_max_iter": 3,
        "psf_fitter_max_iter": 8,
        "psf_redetect_sigma": 3.0,
        "psf_duplicate_radius_fwhm_mult": 1.0,
        "psf_new_sources_cap_per_iter": 100,
        "psf_new_sources_cap_frac": 0.05,
        "psf_fit_init_max_sources": 0,
        "psf_core_cut_enable": False,
        "psf_core_cut_radius_px": 0.0,
        "psf_core_cut_radius_fwhm_mult": 20.0,
        "psf_core_cut_auto_min_density_ratio": 1.5,
        "psf_substar_neighbor_r_fwhm_mult": 8.0,
        "psf_substar_max_sources": 1500,
        "psf_substar_iters": 1,
        "psf_conv_new_frac": 0.03,
        "psf_postfit_snr_min": 3.0,
        "psf_postfit_qfit_max": 3.0,
        "psf_postfit_reduced_chi2_max": 25.0,
        "psf_blend_residual_ratio": 0.25,
        "psf_flux_conv_threshold": 0.01,
        "psf_use_grouper": False,
        "psf_grouper_radius_fwhm": 1.5,
        "psf_forced_match_radius_fwhm": 1.25,
        "psf_redetect_sharp_lo": 0.1,
        "psf_redetect_sharp_hi": 0.95,
        "psf_redetect_round_abs_max": 0.9,
    },
}


def _clone_psf_model(model):
    """Return a per-frame copy of PSF model to avoid thread-shared mutation."""
    try:
        return model.copy()
    except Exception:
        try:
            return copy.deepcopy(model)
        except Exception:
            return model


# ── FITS helpers ───────────────────────────────────────────────────────────────

def _get_filter_lower(fits_path: Path) -> str:
    try:
        h = fits.getheader(fits_path)
        f = h.get("FILTER", None)
        if f is None:
            return "unknown"
        return normalize_filter_name(f)
    except Exception:
        return "unknown"


def _get_exptime(fits_path: Path, default=1.0) -> float:
    try:
        h = fits.getheader(fits_path)
        for k in ("EXPTIME", "EXPOSURE", "ITIME", "ELAPTIME"):
            if k in h:
                v = float(h[k])
                if np.isfinite(v) and v > 0:
                    return v
    except Exception:
        pass
    return float(default)


# ── Detect helpers ────────────────────────────────────────────────────────────

def _load_detect_positions(fname: str, cache_dir: Path, result_dir: Path):
    candidates = [
        cache_dir / f"detect_{fname}.csv",
        step4_dir(result_dir) / f"detect_{fname}.csv",
    ]
    for p in candidates:
        if p.exists() and p.stat().st_size > 0:
            try:
                df = pd.read_csv(p)
                x_col = "x" if "x" in df.columns else ("xcenter" if "xcenter" in df.columns else None)
                y_col = "y" if "y" in df.columns else ("ycenter" if "ycenter" in df.columns else None)
                if x_col is None or y_col is None:
                    continue
                out = pd.DataFrame({"x": pd.to_numeric(df[x_col], errors="coerce"),
                                    "y": pd.to_numeric(df[y_col], errors="coerce")})
                if "det_uid" in df.columns:
                    det_uid_raw = pd.to_numeric(df["det_uid"], errors="coerce")
                    if det_uid_raw.notna().any():
                        missing = det_uid_raw.isna()
                        if missing.any():
                            fallback = np.arange(len(det_uid_raw), dtype=float)
                            det_uid_raw.loc[missing] = fallback[missing.to_numpy()]
                        out["det_uid"] = det_uid_raw.to_numpy(dtype=np.int64, copy=False)
                    else:
                        out["det_uid"] = np.arange(len(df), dtype=np.int64)
                else:
                    out["det_uid"] = np.arange(len(df), dtype=np.int64)
                for flux_col in ("flux_for_quality", "dao_flux", "peak_adu", "dao_peak", "flux", "peak", "amplitude"):
                    if flux_col in df.columns:
                        out["flux_init"] = pd.to_numeric(df[flux_col], errors="coerce")
                        break
                # Pass through morphology quality metrics for EPSF star selection
                for _src, _dst in (
                    ("sharpness", "sharpness"),
                    ("roundness",  "roundness"),
                    ("roundness1", "roundness"),   # DAOStarFinder; prefer plain "roundness"
                    ("elongation", "elong"),
                    ("elong",      "elong"),
                ):
                    if _src in df.columns and _dst not in out.columns:
                        out[_dst] = pd.to_numeric(df[_src], errors="coerce")
                for col in (
                    "quality_score", "nearest_neighbor_px", "nearest_neighbor_fwhm",
                    "edge_margin_px", "fwhm_ratio_to_frame", "flux_percentile",
                    "peak_fraction_to_sat",
                ):
                    if col in df.columns:
                        out[col] = pd.to_numeric(df[col], errors="coerce")
                for col in ("quality_flags",):
                    if col in df.columns:
                        out[col] = df[col].astype(str).reset_index(drop=True)
                for col in ("anchor_candidate", "apcorr_candidate", "epsf_candidate", "psf_seed_candidate"):
                    if col in df.columns:
                        out[col] = (
                            df[col]
                            .astype(str)
                            .str.strip()
                            .str.lower()
                            .isin({"1", "true", "t", "yes", "y"})
                            .reset_index(drop=True)
                        )
                out = out.dropna(subset=["x", "y"]).reset_index(drop=True)
                return out
            except Exception:
                continue
    return None


def _load_fwhm_from_meta(fname: str, cache_dir: Path, result_dir: Path,
                          params_fwhm_guess=6.0) -> float:
    candidates = [
        cache_dir / f"detect_{fname}.json",
        step4_dir(result_dir) / f"detect_{fname}.json",
    ]
    for p in sorted([c for c in candidates if c.exists()],
                    key=lambda q: q.stat().st_mtime_ns, reverse=True):
        try:
            meta = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            continue
        # Prefer explicit FWHM(px); keep radius-derived values as fallback metadata.
        for k in ("fwhm_med_px", "fwhm_px", "fwhm_med", "fwhm_med_rad_px"):
            v = meta.get(k, None)
            if v is not None:
                try:
                    v = float(v)
                    if np.isfinite(v) and v > 0:
                        return v
                except Exception:
                    continue
    return float(params_fwhm_guess)


# ── Moffat PSF builder ────────────────────────────────────────────────────────

def _build_moffat_psf(img_sub: np.ndarray, xy_stars: np.ndarray,
                      fwhm_safe: float, cutout_size: int,
                      log_fn=None):
    """Build normalized Moffat2D PSF from isolated star cutouts.
    Returns (moffat_model, n_good).  Model has x_0=y_0=0, amplitude=1/integral.
    """
    from astropy.modeling.models import Moffat2D
    from astropy.modeling.fitting import LevMarLSQFitter

    h, w = img_sub.shape
    half = cutout_size // 2
    fitter = LevMarLSQFitter()
    yy_c, xx_c = np.mgrid[-half:half + 1, -half:half + 1].astype(float)

    gammas, alphas = [], []
    for xi_f, yi_f in np.asarray(xy_stars, dtype=float):
        xi, yi = int(round(xi_f)), int(round(yi_f))
        if xi - half < 0 or xi + half + 1 > w or yi - half < 0 or yi + half + 1 > h:
            continue
        cutout = img_sub[yi - half:yi + half + 1, xi - half:xi + half + 1].copy()
        if cutout.shape != yy_c.shape:
            continue
        peak = float(np.nanmax(cutout))
        if peak <= 0:
            continue
        alpha_init = 2.5
        gamma_init = max(0.5, fwhm_safe / (2.0 * np.sqrt(2.0 ** (1.0 / alpha_init) - 1.0)))
        try:
            fitted = fitter(
                Moffat2D(amplitude=peak, x_0=0.0, y_0=0.0, gamma=gamma_init, alpha=alpha_init),
                xx_c, yy_c, cutout, maxiter=50,
            )
            g, a = float(fitted.gamma.value), float(fitted.alpha.value)
            if 0.3 < g < 25.0 and 0.5 < a < 10.0:
                gammas.append(g)
                alphas.append(a)
        except Exception:
            continue

    n_good = len(gammas)
    if n_good < 2:
        alpha_med = 2.5
        gamma_med = max(0.5, fwhm_safe / (2.0 * np.sqrt(2.0 ** (1.0 / alpha_med) - 1.0)))
        if log_fn:
            log_fn(f"[MOFFAT] {n_good} good fits; using FWHM estimate gamma={gamma_med:.2f} alpha={alpha_med:.2f}")
    else:
        gamma_med = float(np.median(gammas))
        alpha_med = float(np.median(alphas))
        if log_fn:
            fwhm_est = 2.0 * gamma_med * np.sqrt(2.0 ** (1.0 / alpha_med) - 1.0)
            log_fn(f"[MOFFAT] {n_good} stars: gamma={gamma_med:.3f} alpha={alpha_med:.3f} FWHM≈{fwhm_est:.2f}px")

    # Normalize: integral of Moffat2D = pi*gamma^2/(alpha-1) for alpha>1
    if alpha_med > 1.0:
        integral = np.pi * gamma_med ** 2 / (alpha_med - 1.0)
    else:
        sz = max(cutout_size, 51)
        _h2 = sz // 2
        _yy, _xx = np.mgrid[-_h2:_h2 + 1, -_h2:_h2 + 1].astype(float)
        integral = float(Moffat2D(1.0, 0.0, 0.0, gamma_med, alpha_med)(_xx, _yy).sum())
        if integral <= 0:
            integral = 1.0

    model = Moffat2D(amplitude=1.0 / integral, x_0=0.0, y_0=0.0,
                     gamma=gamma_med, alpha=alpha_med)
    return model, n_good


class MoffatHybridPSF:
    """An analytic Moffat plus the grid of what it missed — DAOPHOT's shape.

    The empirical ePSF and the analytic Moffat are the two extremes APEX
    already had, and DAOPHOT beats both by sitting between them: five fitted
    numbers carry the bulk of the light, and a look-up table carries only the
    leftover (0.4 % of the model's flux on the M13 frame measured 2026-08-14).

    Decomposing does not by itself reduce noise — averaging N stars gives the
    same σ/√N either way. What it removes is *interpolation* error. A star's
    fit samples the model at sub-pixel offsets, and APEX interpolates its ePSF
    linearly; on a core that falls by roughly a factor of two per pixel that
    interpolation is itself a source of error. Here the steep part is evaluated
    from the closed-form Moffat at the exact offset, and only the small, smooth
    residual is interpolated.

    ``residual`` is stored on the same oversampled grid as the ePSF it came
    from, already normalised the same way, so the two pieces add directly.
    """

    def __init__(self, analytic, residual: np.ndarray, oversampling: int):
        self.analytic = analytic
        self.residual = np.asarray(residual, dtype=float)
        self.oversampling = max(1, int(oversampling))

    def _grid_offsets(self) -> tuple[np.ndarray, np.ndarray]:
        ny, nx = self.residual.shape
        yy, xx = np.mgrid[:ny, :nx].astype(float)
        os = self.oversampling
        return ((xx - (nx - 1) / 2.0) / os, (yy - (ny - 1) / 2.0) / os)

    @property
    def data(self) -> np.ndarray:
        """Rendered model on the oversampled grid, for QC and file output.

        Both pieces live on the evaluator's normalisation (values that sum to
        ~1 when sampled at native spacing), so they add with no extra factor.
        """
        dx, dy = self._grid_offsets()
        return np.asarray(self.analytic(dx, dy), dtype=float) + self.residual


def build_moffat_hybrid_psf(epsf_model, analytic, oversampling: int) -> MoffatHybridPSF:
    """Split an ePSF into an analytic core and the residual it leaves behind.

    The ePSF is the empirical average of the reference stars, so subtracting
    the analytic evaluated on the same grid gives exactly "what the analytic
    missed", averaged over the same stars — no separate stacking pass and no
    new star-selection policy.
    """
    grid = np.asarray(epsf_model.data, dtype=float)
    os = max(1, int(oversampling))
    ny, nx = grid.shape
    yy, xx = np.mgrid[:ny, :nx].astype(float)
    dx = (xx - (nx - 1) / 2.0) / os
    dy = (yy - (ny - 1) / 2.0) / os

    # Put the ePSF on the evaluator's normalisation — the same `sum()/os**2`
    # the ePSF evaluator divides by — so that sampling it at native spacing
    # sums to 1. The analytic already integrates to 1 per native pixel, so the
    # two are directly comparable and the residual is their difference.
    norm = grid.sum() / os ** 2
    if not np.isfinite(norm) or norm <= 0:
        norm = 1.0
    residual = grid / norm - np.asarray(analytic(dx, dy), dtype=float)
    return MoffatHybridPSF(analytic, residual, os)


# ── Unified PSF evaluator (EPSF, Moffat, or the hybrid) ──────────────────────

def _make_psf_evaluator(psf_model, psf_type: str, oversampling: int = 2,
                        interp_order: int = 1):
    """Return eval_fn(dx_2d, dy_2d) -> normalized PSF values.

    dx, dy are pixel offsets from star centre.  Output sums to ≈ 1 / pixel².

    ``interp_order`` is the spline order used to sample the oversampled grid.
    It is a knob rather than a constant because the grid is what a sub-pixel
    fit interpolates, and on a core that falls by roughly a factor of two per
    pixel the interpolation is itself a source of error — the same error the
    hybrid model avoids by evaluating its analytic part in closed form. Raising
    the order is the cheap control that separates "the hybrid helped because it
    is analytic" from "it helped because linear interpolation was the problem".
    """
    order = int(np.clip(int(interp_order), 1, 5))
    if psf_type == 'moffat_hybrid':
        from scipy.ndimage import map_coordinates as _mc
        res = psf_model.residual
        os = max(1, int(psf_model.oversampling))
        cy, cx = res.shape[0] // 2, res.shape[1] // 2

        def _eval_hybrid(dx, dy):
            dx_a = np.asarray(dx, dtype=float)
            dy_a = np.asarray(dy, dtype=float)
            # Steep part: closed form at the exact offset, no interpolation.
            core = np.asarray(psf_model.analytic(dx_a, dy_a), dtype=float)
            # Leftover: small and smooth, so linear interpolation is cheap here.
            # Already on the evaluator's normalisation — no os**2 factor.
            vals = _mc(res, [dy_a.ravel() * os + cy, dx_a.ravel() * os + cx],
                       order=order, mode='constant', cval=0.0)
            return core + vals.reshape(dx_a.shape)
        return _eval_hybrid

    if psf_type == 'moffat':
        def _eval_moffat(dx, dy):
            return np.asarray(psf_model(np.asarray(dx, dtype=float),
                                        np.asarray(dy, dtype=float)), dtype=float)
        return _eval_moffat

    # EPSF
    from scipy.ndimage import map_coordinates as _mc
    psf_data = np.asarray(psf_model.data, dtype=float)
    os = max(1, int(oversampling))
    cy = psf_data.shape[0] // 2
    cx = psf_data.shape[1] // 2
    norm = psf_data.sum() / os ** 2
    if norm <= 0:
        norm = 1.0

    def _eval_epsf(dx, dy):
        dx_a = np.asarray(dx, dtype=float)
        dy_a = np.asarray(dy, dtype=float)
        coords_y = dy_a.ravel() * os + cy
        coords_x = dx_a.ravel() * os + cx
        vals = _mc(psf_data, [coords_y, coords_x], order=order, mode='constant', cval=0.0)
        return (vals / norm).reshape(dx_a.shape)

    return _eval_epsf


def _sample_native_psf(eval_psf, support_size: int) -> np.ndarray:
    """Sample a centered native-pixel PSF for fit-window and noise policy."""
    size = max(3, int(support_size))
    if size % 2 == 0:
        size += 1
    half = size // 2
    yy, xx = np.mgrid[-half:half + 1, -half:half + 1]
    return np.asarray(eval_psf(xx, yy), dtype=float)


# ── APEX iterative engine (ALLSTAR-inspired) ──────────────────────────────────

_NEWTON_GRAD_DELTA = 0.5  # sub-pixel step for numerical PSF gradient (pixels)
# Independent real-image injections in M13 and M3 show that single-source fits
# inside 1.5 FWHM can retain acceptable qfit/chi2 while absorbing neighbour flux.
_UNRESOLVED_NEIGHBOR_FWHM = 1.5


def _allstar_newton_one(cleaned_patch: np.ndarray,
                        x0: float, y0: float,
                        patch_y0: int, patch_x0: int,
                        flux0: float, eval_psf,
                        max_shift: float = 2.0,
                        weights: np.ndarray | None = None,
                        position_fixed: bool = False,
                        allow_negative_flux: bool = False):
    """Single linearized Newton step for one star (true DAOPHOT ALLSTAR style).

    Solves weighted normal equations for flux, position, and local sky.
    One matrix solve replaces 20-40 LM iterations.
    Returns (x_new, y_new, flux_new, chi2, ok).
    """
    ny, nx = cleaned_patch.shape
    if ny < 3 or nx < 3:
        flux_out = float(flux0) if allow_negative_flux and np.isfinite(flux0) else max(flux0, 1.0)
        return x0, y0, flux_out, np.nan, False

    xc = float(x0)
    yc = float(y0)
    if allow_negative_flux:
        flux_safe = float(flux0) if np.isfinite(flux0) else 0.0
    else:
        flux_safe = max(float(flux0), 1.0)
    d = _NEWTON_GRAD_DELTA

    yy = np.arange(ny, dtype=float) + patch_y0
    xx = np.arange(nx, dtype=float) + patch_x0
    YY, XX = np.meshgrid(yy, xx, indexing='ij')

    # PSF value and partial derivatives at current star centre
    c_f = eval_psf(XX - xc, YY - yc)
    c_x = (eval_psf(XX - xc - d, YY - yc) - eval_psf(XX - xc + d, YY - yc)) / (2.0 * d)
    c_y = (eval_psf(XX - xc, YY - yc - d) - eval_psf(XX - xc, YY - yc + d)) / (2.0 * d)

    residual = (cleaned_patch - c_f * flux_safe).ravel()
    Cf = c_f.ravel()
    Cx = (flux_safe * c_x).ravel()
    Cy = (flux_safe * c_y).ravel()

    columns = [Cf]
    if not position_fixed:
        columns.extend([Cx, Cy])
    columns.append(np.ones_like(Cf))
    design = np.column_stack(columns)
    if weights is None:
        sqrt_weight = np.ones_like(residual)
    else:
        weight_arr = np.asarray(weights, dtype=float).ravel()
        sqrt_weight = np.sqrt(
            np.where(np.isfinite(weight_arr) & (weight_arr > 0), weight_arr, 0.0)
        )
    design_weighted = design * sqrt_weight[:, None]
    residual_weighted = residual * sqrt_weight

    try:
        params, _, rank, _ = np.linalg.lstsq(
            design_weighted,
            residual_weighted,
            rcond=None,
        )
        if rank < design.shape[1]:
            raise np.linalg.LinAlgError("rank-deficient PSF fit")
    except (np.linalg.LinAlgError, ValueError):
        return x0, y0, flux_safe, np.nan, False

    dflux = float(params[0])
    if position_fixed:
        dx = 0.0
        dy = 0.0
    else:
        dx = float(params[1])
        dy = float(params[2])

    if abs(dx) > max_shift or abs(dy) > max_shift:
        return x0, y0, flux_safe, np.nan, False

    flux_new = flux_safe + dflux
    if not allow_negative_flux and flux_new < flux_safe * 0.1:
        flux_new = flux_safe * 0.5

    local_sky = float(params[-1])
    model_new = eval_psf(XX - (xc + dx), YY - (yc + dy)) * flux_new
    res_new = (cleaned_patch - model_new - local_sky).ravel()
    qfit = float(np.sum(np.abs(res_new))) / max(abs(float(flux_new)), 1e-20)

    return xc + dx, yc + dy, flux_new, qfit, True


def _allstar_newton_group(cleaned_patch: np.ndarray,
                          group_info: list,
                          patch_y0: int, patch_x0: int,
                          eval_psf, max_shift: float = 2.0,
                          weights: np.ndarray | None = None,
                          position_fixed: bool | np.ndarray = False,
                          allow_negative_flux: bool | np.ndarray = False):
    """Single Newton step for N close stars simultaneously (3N×3N normal equations).

    group_info: list of (x, y, flux) — absolute image coordinates.
    Returns list of (x_new, y_new, flux_new, chi2, ok).
    """
    N = len(group_info)
    if N == 0:
        return []
    fixed_mask = np.asarray(position_fixed, dtype=bool)
    if fixed_mask.ndim == 0:
        fixed_mask = np.full(N, bool(fixed_mask), dtype=bool)
    if fixed_mask.shape != (N,):
        raise ValueError("position_fixed must be scalar or have one value per source")
    signed_mask = np.asarray(allow_negative_flux, dtype=bool)
    if signed_mask.ndim == 0:
        signed_mask = np.full(N, bool(signed_mask), dtype=bool)
    if signed_mask.shape != (N,):
        raise ValueError("allow_negative_flux must be scalar or have one value per source")
    if N == 1:
        x0, y0, fl0 = group_info[0]
        return [_allstar_newton_one(cleaned_patch, x0, y0, patch_y0, patch_x0,
                                    fl0, eval_psf, max_shift, weights,
                                    bool(fixed_mask[0]), bool(signed_mask[0]))]

    ny, nx = cleaned_patch.shape
    if ny < 3 or nx < 3:
        return [(x, y, max(fl, 1.0), np.nan, False) for x, y, fl in group_info]

    yy_abs = np.arange(ny, dtype=float) + patch_y0
    xx_abs = np.arange(nx, dtype=float) + patch_x0
    YY, XX = np.meshgrid(yy_abs, xx_abs, indexing='ij')
    n_pix = ny * nx
    d = _NEWTON_GRAD_DELTA

    xc_arr = np.array([float(x) for x, y, fl in group_info])
    yc_arr = np.array([float(y) for x, y, fl in group_info])
    fl_arr = np.array([
        float(fl) if signed and np.isfinite(float(fl)) else max(float(fl), 1.0)
        for (_, _, fl), signed in zip(group_info, signed_mask)
    ])

    # Every source has a flux delta. Only non-forced sources add dx/dy, so a
    # local group may safely mix catalog-anchored and freely-centroided stars.
    parameter_columns: list[tuple[int, int | None, int | None]] = []
    n_source_parameters = 0
    for fixed in fixed_mask:
        flux_column = n_source_parameters
        if fixed:
            parameter_columns.append((flux_column, None, None))
            n_source_parameters += 1
        else:
            parameter_columns.append((flux_column, flux_column + 1, flux_column + 2))
            n_source_parameters += 3
    # One final column solves a local constant sky for the whole group patch.
    A_mat = np.zeros((n_pix, n_source_parameters + 1), dtype=float)
    model = np.zeros(n_pix, dtype=float)

    for n in range(N):
        xc, yc, fl = xc_arr[n], yc_arr[n], fl_arr[n]
        c_f = eval_psf(XX - xc, YY - yc).ravel()
        flux_column, dx_column, dy_column = parameter_columns[n]
        A_mat[:, flux_column] = c_f
        if dx_column is not None and dy_column is not None:
            c_x = ((eval_psf(XX - xc - d, YY - yc) -
                    eval_psf(XX - xc + d, YY - yc)) / (2.0 * d)).ravel()
            c_y = ((eval_psf(XX - xc, YY - yc - d) -
                    eval_psf(XX - xc, YY - yc + d)) / (2.0 * d)).ravel()
            A_mat[:, dx_column] = fl * c_x
            A_mat[:, dy_column] = fl * c_y
        model += c_f * fl
    A_mat[:, -1] = 1.0

    residual = cleaned_patch.ravel() - model
    if weights is None:
        sqrt_weight = np.ones_like(residual)
    else:
        weight_arr = np.asarray(weights, dtype=float).ravel()
        sqrt_weight = np.sqrt(
            np.where(np.isfinite(weight_arr) & (weight_arr > 0), weight_arr, 0.0)
        )
    design_weighted = A_mat * sqrt_weight[:, None]
    residual_weighted = residual * sqrt_weight

    try:
        if A_mat.shape[1] >= 13:
            from scipy.sparse import csr_matrix
            from scipy.sparse.linalg import lsqr

            sparse_result = lsqr(
                csr_matrix(design_weighted),
                residual_weighted,
                atol=1e-6,
                btol=1e-6,
                iter_lim=max(50, 4 * A_mat.shape[1]),
            )
            params = np.asarray(sparse_result[0], dtype=float)
            if int(sparse_result[1]) not in {0, 1, 2} or not np.all(np.isfinite(params)):
                raise np.linalg.LinAlgError("sparse group fit did not converge")
        else:
            params, _, rank, _ = np.linalg.lstsq(
                design_weighted,
                residual_weighted,
                rcond=None,
            )
            if rank < A_mat.shape[1]:
                raise np.linalg.LinAlgError("rank-deficient group fit")
    except (np.linalg.LinAlgError, ValueError):
        return [(x, y, max(fl, 1.0), np.nan, False) for x, y, fl in group_info]

    results = []
    for n in range(N):
        flux_column, dx_column, dy_column = parameter_columns[n]
        dflux = float(params[flux_column])
        if fixed_mask[n]:
            dx = 0.0
            dy = 0.0
        else:
            dx = float(params[dx_column])
            dy = float(params[dy_column])
        x0, y0, fl0 = group_info[n]
        if abs(dx) > max_shift or abs(dy) > max_shift:
            fallback_flux = float(fl0) if signed_mask[n] else max(float(fl0), 1.0)
            results.append((x0, y0, fallback_flux, np.nan, False))
            continue
        flux_new = fl_arr[n] + dflux
        if not signed_mask[n] and flux_new < fl_arr[n] * 0.1:
            flux_new = fl_arr[n] * 0.5
        group_model = np.zeros(n_pix, dtype=float)
        for source_index in range(N):
            source_flux_column, source_dx_column, source_dy_column = parameter_columns[source_index]
            source_flux = fl_arr[source_index] + float(params[source_flux_column])
            if fixed_mask[source_index]:
                source_dx = source_dy = 0.0
            else:
                source_dx = float(params[source_dx_column])
                source_dy = float(params[source_dy_column])
            group_model += eval_psf(
                XX - (xc_arr[source_index] + source_dx),
                YY - (yc_arr[source_index] + source_dy),
            ).ravel() * source_flux
        group_residual = residual + model - group_model - float(params[-1])
        qfit = float(np.sum(np.abs(group_residual))) / max(abs(float(flux_new)), 1e-20)
        results.append((xc_arr[n] + dx, yc_arr[n] + dy, flux_new, qfit, True))

    return results

def _build_groups(x_arr: np.ndarray, y_arr: np.ndarray, f_arr: np.ndarray,
                  radius: float, max_size: int,
                  max_grouped_sources: int = 0) -> list:
    """Greedy brightest-first group assignment for simultaneous PSF fitting.
    Stars within `radius` pixels of each other are grouped together (max `max_size`).
    Each star appears in exactly one group.  Returns list of index-lists.
    """
    N = len(x_arr)
    if N == 0 or radius <= 0 or max_size <= 1:
        return [[i] for i in range(N)]

    xy = np.column_stack([np.asarray(x_arr, dtype=float),
                          np.asarray(y_arr, dtype=float)])
    tree = cKDTree(xy)
    assigned = np.zeros(N, dtype=bool)
    groups: list = []

    for i in np.argsort(f_arr)[::-1]:  # brightest first
        if assigned[i]:
            continue
        raw = tree.query_ball_point([float(x_arr[i]), float(y_arr[i])], r=radius)
        neighbors = np.array([j for j in raw if not assigned[j]], dtype=int)
        if len(neighbors) > max_size:
            dists = np.hypot(x_arr[neighbors] - x_arr[i],
                             y_arr[neighbors] - y_arr[i])
            neighbors = neighbors[np.argsort(dists)[:max_size]]
        for j in neighbors:
            assigned[j] = True
        groups.append(neighbors.tolist())

    if max_grouped_sources <= 0:
        return groups

    limited_groups: list[list[int]] = []
    grouped_count = 0
    for group in groups:
        if len(group) <= 1:
            limited_groups.append(group)
            continue
        if grouped_count + len(group) <= max_grouped_sources:
            limited_groups.append(group)
            grouped_count += len(group)
        else:
            limited_groups.extend([[index] for index in group])
    return limited_groups


def _allstar_fit_group(cleaned_patch: np.ndarray,
                       group_info: list,
                       patch_y0: int, patch_x0: int,
                       eval_psf, max_shift: float):
    """Simultaneously fit N stars on a pre-neighbour-subtracted patch.

    group_info : list of (x, y, flux) — absolute image coordinates.
    Returns    : list of (x_fit, y_fit, flux_fit, chi2, ok) per star.

    Parameters are subpixel offsets (dx, dy) + log-flux for numerical stability.
    Each dx/dy is relative to the integer-rounded star centre, keeping values near 0.
    log-flux prevents negative-flux blow-up during LM iterations.
    """
    from scipy.optimize import least_squares

    N = len(group_info)
    if N == 0:
        return []
    if N == 1:
        x0, y0, fl0 = group_info[0]
        return [_allstar_fit_one(cleaned_patch, x0, y0, patch_y0, patch_x0,
                                 fl0, eval_psf, max_shift)]

    ny, nx = cleaned_patch.shape
    if ny < 3 or nx < 3:
        return [(x, y, max(fl, 1.0), np.nan, False) for x, y, fl in group_info]

    yy_abs = np.arange(ny, dtype=float) + patch_y0
    xx_abs = np.arange(nx, dtype=float) + patch_x0
    YY, XX = np.meshgrid(yy_abs, xx_abs, indexing='ij')

    # Integer reference centres (keep dx/dy near zero for good LM conditioning)
    xi_refs = np.array([int(round(float(x))) for x, y, fl in group_info], dtype=int)
    yi_refs = np.array([int(round(float(y))) for x, y, fl in group_info], dtype=int)

    # p = [dx1, dy1, log_fl1, dx2, dy2, log_fl2, ...]
    # Using log-flux so flux stays positive and scale is comparable to dx/dy
    p0 = []
    fl_refs = []
    for n, (x, y, fl) in enumerate(group_info):
        fl_safe = max(float(fl), 1.0)
        fl_refs.append(fl_safe)
        p0.extend([float(x) - xi_refs[n], float(y) - yi_refs[n],
                   float(np.log(fl_safe))])
    p0 = np.array(p0, dtype=float)

    def _res(p):
        model = np.zeros((ny, nx), dtype=float)
        for n in range(N):
            dx, dy = p[3 * n], p[3 * n + 1]
            fl = float(np.exp(min(p[3 * n + 2], 30.0)))  # cap prevents overflow
            xc = xi_refs[n] + dx
            yc = yi_refs[n] + dy
            model += eval_psf(XX - xc, YY - yc) * fl
        diff = (cleaned_patch - model).ravel()
        return np.where(np.isfinite(diff), diff, 0.0)

    try:
        r = least_squares(_res, p0, method='lm',
                          ftol=1e-4, xtol=0.01, gtol=1e-6,
                          max_nfev=40 * N)
        chi2_base = float(np.mean(r.fun ** 2))
        out = []
        for n in range(N):
            dx_f, dy_f = r.x[3 * n], r.x[3 * n + 1]
            fl_f = float(np.exp(min(r.x[3 * n + 2], 30.0)))
            x0, y0, fl0 = group_info[n]
            xc_f = xi_refs[n] + dx_f
            yc_f = yi_refs[n] + dy_f
            if abs(dx_f) > max_shift or abs(dy_f) > max_shift or fl_f <= 0:
                out.append((x0, y0, max(fl0, 1.0), np.nan, False))
            else:
                out.append((xc_f, yc_f, fl_f,
                             chi2_base / max(fl_f ** 2 * 1e-8, 1e-20), True))
        return out
    except Exception:
        return [(x, y, max(fl, 1.0), np.nan, False) for x, y, fl in group_info]

def _allstar_stamp(img_shape, x_cen: float, y_cen: float, flux: float,
                   eval_psf, stamp_size: int):
    """Return (slice_y, slice_x, stamp_array) or (None, None, None)."""
    h, w = img_shape
    half = stamp_size // 2
    xi, yi = int(round(x_cen)), int(round(y_cen))
    y_lo, y_hi = max(0, yi - half), min(h, yi + half + 1)
    x_lo, x_hi = max(0, xi - half), min(w, xi + half + 1)
    if y_hi <= y_lo or x_hi <= x_lo:
        return None, None, None
    yy = np.arange(y_lo, y_hi, dtype=float) - yi
    xx = np.arange(x_lo, x_hi, dtype=float) - xi
    YY, XX = np.meshgrid(yy, xx, indexing='ij')
    dx_sub = x_cen - xi
    dy_sub = y_cen - yi
    stamp = (eval_psf(XX - dx_sub, YY - dy_sub) * flux).astype(np.float32)
    return slice(y_lo, y_hi), slice(x_lo, x_hi), stamp


def _allstar_build_model(img_shape, x_arr, y_arr, f_arr, eval_psf, stamp_size: int):
    model = np.zeros(img_shape, dtype=np.float32)
    _allstar_apply_model_inplace(
        model, x_arr, y_arr, f_arr, eval_psf, stamp_size, subtract=False
    )
    return model


def _allstar_apply_model_inplace(
    image: np.ndarray,
    x_arr,
    y_arr,
    f_arr,
    eval_psf,
    stamp_size: int,
    *,
    subtract: bool,
) -> None:
    """Accumulate source stamps without allocating a full-frame model."""
    for xi, yi, fi in zip(x_arr, y_arr, f_arr):
        sy, sx, stamp = _allstar_stamp(
            image.shape, xi, yi, fi, eval_psf, stamp_size
        )
        if sy is not None:
            if subtract:
                np.subtract(image[sy, sx], stamp, out=image[sy, sx])
            else:
                np.add(image[sy, sx], stamp, out=image[sy, sx])


def _float32_difference(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Return ``left - right`` with one full-frame float32 allocation."""
    if np.shape(left) != np.shape(right):
        raise ValueError("difference operands must have the same shape")
    result = np.empty(np.shape(left), dtype=np.float32)
    np.subtract(left, right, out=result, casting="unsafe")
    return result


def _allstar_fit_one(cleaned_patch: np.ndarray,
                     x0: float, y0: float,
                     patch_y0: int, patch_x0: int,
                     flux0: float, eval_psf,
                     max_shift: float = 2.0):
    """Fit one source on a pre-cleaned local patch.

    cleaned_patch : 2D array, already neighbour-subtracted, positioned at
                    image coords [patch_y0:patch_y0+ny, patch_x0:patch_x0+nx].
    Returns (x_fit, y_fit, flux_fit, chi2, ok).
    """
    from scipy.optimize import least_squares
    ny, nx = cleaned_patch.shape
    if ny < 3 or nx < 3:
        return x0, y0, flux0, np.nan, False

    xi, yi = int(round(x0)), int(round(y0))
    # Pixel offset grids (relative to integer star centre)
    yy = np.arange(ny, dtype=float) + patch_y0 - yi
    xx = np.arange(nx, dtype=float) + patch_x0 - xi
    YY, XX = np.meshgrid(yy, xx, indexing='ij')

    dx0, dy0 = x0 - xi, y0 - yi
    flux_safe = max(float(flux0) if np.isfinite(flux0) else 1.0, 1.0)

    # Cluster cores retain unresolved stellar light after global sky removal.
    # Fit a local constant term so that diffuse core light is not forced into
    # the target star's PSF flux.
    edge_mask = np.zeros(cleaned_patch.shape, dtype=bool)
    edge_mask[0, :] = True
    edge_mask[-1, :] = True
    edge_mask[:, 0] = True
    edge_mask[:, -1] = True
    edge_vals = cleaned_patch[edge_mask]
    edge_vals = edge_vals[np.isfinite(edge_vals)]
    if edge_vals.size:
        bg0 = float(np.nanmedian(edge_vals))
        bg_scale = float(_mad_std(edge_vals)) if edge_vals.size > 3 else float(np.nanstd(edge_vals))
    else:
        bg0 = 0.0
        bg_scale = 0.0
    if not np.isfinite(bg0):
        bg0 = 0.0
    if not np.isfinite(bg_scale) or bg_scale <= 0:
        bg_scale = max(1.0, abs(bg0) * 0.1)

    def _res(p):
        dx, dy, fl, bg = p
        diff = (cleaned_patch - (eval_psf(XX - dx, YY - dy) * fl + bg)).ravel()
        return np.where(np.isfinite(diff), diff, 0.0)

    try:
        bg_pad = max(5.0 * bg_scale, abs(bg0) + 10.0)
        r = least_squares(
            _res,
            [dx0, dy0, flux_safe, bg0],
            bounds=(
                [-max_shift, -max_shift, 1e-12, bg0 - bg_pad],
                [ max_shift,  max_shift, np.inf, bg0 + bg_pad],
            ),
            method='trf',
            ftol=1e-4,
            xtol=0.01,
            gtol=1e-6,
            max_nfev=80,
        )
        dx_f, dy_f, fl_f, _bg_f = r.x
        if abs(dx_f) > max_shift or abs(dy_f) > max_shift or fl_f <= 0:
            return x0, y0, flux_safe, np.nan, False
        chi2 = float(np.mean(r.fun ** 2)) / max(fl_f ** 2 * 1e-8, 1e-20)
        return xi + dx_f, yi + dy_f, fl_f, chi2, True
    except Exception:
        return x0, y0, flux_safe, np.nan, False


def _allstar_fit(img_sub: np.ndarray, positions: np.ndarray, fluxes: np.ndarray,
                 eval_psf, fit_shape: int, stamp_size: int,
                 max_iter: int, flux_conv: float,
                 max_shift: float = 2.5,
                 group_radius: float = 0.0,
                 max_group_size: int = 1,
                 max_grouped_sources: int = 0,
                 background_rms: float = 1.0,
                 gain: float = 1.0,
                 initial_positions: np.ndarray | None = None,
                 initial_fit_valid: np.ndarray | None = None,
                 position_bound: float | None = None,
                 position_fixed: bool = False,
                 position_fixed_mask: np.ndarray | None = None,
                 allow_negative_flux_mask: np.ndarray | None = None,
                 fit_active_mask: np.ndarray | None = None,
                 log_fn=None, stop_fn=None):
    """DAOPHOT ALLSTAR-style iterative PSF fitting.

    Each iteration: for every star, subtract neighbours, fit star plus a local
    constant background, then update the model with a delta stamp.
    O(N) per iteration vs O(N_group × size²).

    Returns astropy Table with columns matching photutils output format.
    """
    from astropy.table import Table as _Tab
    N = len(positions)
    if N == 0:
        return _Tab({'x_fit': np.array([]), 'y_fit': np.array([]),
                     'flux_fit': np.array([]), 'flux_err': np.array([]),
                     'qfit': np.array([]), 'cfit': np.array([]),
                     'reduced_chi2': np.array([]), 'flags': np.array([], dtype=int),
                     'n_pixels_fit': np.array([], dtype=int),
                     'iter_detected': np.array([], dtype=int)})

    h, w = img_sub.shape
    fit_half = fit_shape // 2
    x = positions[:, 0].copy().astype(float)
    y = positions[:, 1].copy().astype(float)
    fixed_mask = np.full(N, bool(position_fixed), dtype=bool)
    if position_fixed_mask is not None:
        supplied_fixed = np.asarray(position_fixed_mask, dtype=bool)
        if supplied_fixed.shape == (N,):
            fixed_mask |= supplied_fixed
    signed_flux_mask = np.zeros(N, dtype=bool)
    if allow_negative_flux_mask is not None:
        supplied_signed = np.asarray(allow_negative_flux_mask, dtype=bool)
        if supplied_signed.shape == (N,):
            signed_flux_mask = supplied_signed & fixed_mask
    active_mask = np.ones(N, dtype=bool)
    if fit_active_mask is not None:
        supplied_active = np.asarray(fit_active_mask, dtype=bool)
        if supplied_active.shape != (N,):
            raise ValueError("fit_active_mask must have shape (N,)")
        active_mask = supplied_active.copy()
    flux_arr = np.asarray(fluxes, dtype=float)
    f = np.where(
        signed_flux_mask,
        np.where(np.isfinite(flux_arr), flux_arr, 0.0),
        np.where(np.isfinite(flux_arr) & (flux_arr > 0), flux_arr, 1.0),
    ).astype(float, copy=True)
    chi2 = np.full(N, np.nan)
    fit_ok = (
        np.asarray(initial_fit_valid, dtype=bool).copy()
        if initial_fit_valid is not None
        and np.asarray(initial_fit_valid).shape == (N,)
        else np.zeros(N, dtype=bool)
    )
    fit_flags = np.zeros(N, dtype=np.int32)
    anchors = np.asarray(
        initial_positions if initial_positions is not None else positions,
        dtype=float,
    )
    if anchors.shape != (N, 2):
        anchors = np.asarray(positions, dtype=float).copy()
    bound = float(position_bound) if position_bound is not None else float(max_shift)
    bound = max(0.0, bound)
    last_changed = np.zeros(N, dtype=bool)
    exhausted_iterations = True

    use_groups = (group_radius > 0 and max_group_size > 1)
    model_img = _allstar_build_model(img_sub.shape, x, y, f, eval_psf, stamp_size)
    if log_fn:
        log_fn(
            f"  [APEX] N={N} fit_shape={fit_shape} stamp={stamp_size} "
            f"max_iter={max_iter} active={int(np.sum(active_mask))} "
            f"grouping={'on r=%.1fpx max=%d budget=%d' % (group_radius, max_group_size, max_grouped_sources) if use_groups else 'off'}"
        )

    for it in range(max_iter):
        if stop_fn and stop_fn():
            break
        n_changed = 0
        max_df = 0.0
        changed_this_iter = np.zeros(N, dtype=bool)

        # Build groups for this iteration (brightest-first)
        active_indices = np.flatnonzero(active_mask)
        if use_groups and active_indices.size:
            local_groups = _build_groups(
                x[active_indices],
                y[active_indices],
                f[active_indices],
                group_radius,
                max_group_size,
                max_grouped_sources=max_grouped_sources,
            )
            groups = [
                active_indices[np.asarray(group, dtype=int)].tolist()
                for group in local_groups
            ]
        else:
            ordered = active_indices[np.argsort(f[active_indices])[::-1]]
            groups = [[int(i)] for i in ordered]  # brightest-first, singles only

        for group in groups:
            if stop_fn and stop_fn():
                break

            if len(group) == 1:
                # ── Single-star fit (unchanged) ───────────────────────────
                i = group[0]
                xi, yi = int(round(x[i])), int(round(y[i]))
                fy_lo = max(0, yi - fit_half)
                fy_hi = min(h, yi + fit_half + 1)
                fx_lo = max(0, xi - fit_half)
                fx_hi = min(w, xi + fit_half + 1)
                if fy_hi - fy_lo < 3 or fx_hi - fx_lo < 3:
                    continue
                sy_old, sx_old, stamp_old = _allstar_stamp(
                    img_sub.shape, x[i], y[i], f[i], eval_psf, stamp_size)
                if sy_old is None:
                    continue
                fit_raw   = img_sub  [fy_lo:fy_hi, fx_lo:fx_hi].copy()
                fit_model = model_img[fy_lo:fy_hi, fx_lo:fx_hi].copy()
                y0_s, y1_s = sy_old.start, sy_old.stop
                x0_s, x1_s = sx_old.start, sx_old.stop
                oy_lo, oy_hi = max(fy_lo, y0_s), min(fy_hi, y1_s)
                ox_lo, ox_hi = max(fx_lo, x0_s), min(fx_hi, x1_s)
                if oy_hi > oy_lo and ox_hi > ox_lo:
                    fit_model[oy_lo - fy_lo:oy_hi - fy_lo,
                              ox_lo - fx_lo:ox_hi - fx_lo] -= \
                        stamp_old[oy_lo - y0_s:oy_hi - y0_s,
                                  ox_lo - x0_s:ox_hi - x0_s]
                cleaned = (fit_raw - fit_model).astype(np.float32, copy=False)
                variance = (
                    max(float(background_rms), 1e-6) ** 2
                    + np.clip(fit_raw, 0.0, None) / max(float(gain), 1e-6)
                )
                weights = np.where(np.isfinite(variance) & (variance > 0), 1.0 / variance, 0.0)
                x_new, y_new, f_new, chi2_i, ok = _allstar_newton_one(
                    cleaned,
                    x[i],
                    y[i],
                    fy_lo,
                    fx_lo,
                    f[i],
                    eval_psf,
                    max_shift,
                    weights,
                    bool(fixed_mask[i]),
                    bool(signed_flux_mask[i]),
                )
                if ok and not fixed_mask[i] and bound > 0:
                    if (
                        abs(float(x_new) - float(anchors[i, 0])) > bound
                        or abs(float(y_new) - float(anchors[i, 1])) > bound
                    ):
                        fit_flags[i] |= int(PSFFitFlag.NEAR_BOUND | PSFFitFlag.NONCONVERGENCE)
                        ok = False
                if ok:
                    sy_new, sx_new, stamp_new = _allstar_stamp(
                        img_sub.shape, x_new, y_new, f_new, eval_psf, stamp_size)
                    model_img[sy_old, sx_old] -= stamp_old
                    if sy_new is not None:
                        model_img[sy_new, sx_new] += stamp_new
                    df = abs(f_new - f[i]) / max(abs(f[i]), 1e-10)
                    changed = fit_parameters_changed(
                        x[i],
                        y[i],
                        f[i],
                        x_new,
                        y_new,
                        f_new,
                        flux_fraction=flux_conv,
                    )
                    if changed:
                        n_changed += 1
                        changed_this_iter[i] = True
                    max_df = max(max_df, df)
                    x[i], y[i], f[i], chi2[i] = x_new, y_new, f_new, chi2_i
                    fit_ok[i] = True

            else:
                # ── Multi-star simultaneous fit ───────────────────────────
                xi_arr = [int(round(x[k])) for k in group]
                yi_arr = [int(round(y[k])) for k in group]
                # Group patch: union of all fit windows
                gy_lo = max(0, min(yi_arr) - fit_half)
                gy_hi = min(h, max(yi_arr) + fit_half + 1)
                gx_lo = max(0, min(xi_arr) - fit_half)
                gx_hi = min(w, max(xi_arr) + fit_half + 1)
                if gy_hi - gy_lo < 3 or gx_hi - gx_lo < 3:
                    continue

                # Cache old stamps
                old_stamps = [_allstar_stamp(img_sub.shape, x[k], y[k], f[k],
                                             eval_psf, stamp_size) for k in group]

                # Build cleaned group patch
                fit_raw   = img_sub  [gy_lo:gy_hi, gx_lo:gx_hi].copy()
                fit_model = model_img[gy_lo:gy_hi, gx_lo:gx_hi].copy()
                for sy_k, sx_k, stamp_k in old_stamps:
                    if sy_k is None:
                        continue
                    y0_s, y1_s = sy_k.start, sy_k.stop
                    x0_s, x1_s = sx_k.start, sx_k.stop
                    oy_lo2, oy_hi2 = max(gy_lo, y0_s), min(gy_hi, y1_s)
                    ox_lo2, ox_hi2 = max(gx_lo, x0_s), min(gx_hi, x1_s)
                    if oy_hi2 > oy_lo2 and ox_hi2 > ox_lo2:
                        fit_model[oy_lo2 - gy_lo:oy_hi2 - gy_lo,
                                  ox_lo2 - gx_lo:ox_hi2 - gx_lo] -= \
                            stamp_k[oy_lo2 - y0_s:oy_hi2 - y0_s,
                                    ox_lo2 - x0_s:ox_hi2 - x0_s]
                cleaned = (fit_raw - fit_model).astype(np.float32, copy=False)
                variance = (
                    max(float(background_rms), 1e-6) ** 2
                    + np.clip(fit_raw, 0.0, None) / max(float(gain), 1e-6)
                )
                weights = np.where(np.isfinite(variance) & (variance > 0), 1.0 / variance, 0.0)

                group_info = [(x[k], y[k], f[k]) for k in group]
                results = _allstar_newton_group(cleaned, group_info,
                                                gy_lo, gx_lo, eval_psf, max_shift,
                                                weights, fixed_mask[group],
                                                signed_flux_mask[group])

                # Apply results and update model
                for idx, k in enumerate(group):
                    x_new, y_new, f_new, chi2_k, ok = results[idx]
                    sy_old_k, sx_old_k, stamp_old_k = old_stamps[idx]
                    if ok and not fixed_mask[k] and bound > 0:
                        if (
                            abs(float(x_new) - float(anchors[k, 0])) > bound
                            or abs(float(y_new) - float(anchors[k, 1])) > bound
                        ):
                            fit_flags[k] |= int(PSFFitFlag.NEAR_BOUND | PSFFitFlag.NONCONVERGENCE)
                            ok = False
                    if ok:
                        sy_new_k, sx_new_k, stamp_new_k = _allstar_stamp(
                            img_sub.shape, x_new, y_new, f_new, eval_psf, stamp_size)
                        if sy_old_k is not None:
                            model_img[sy_old_k, sx_old_k] -= stamp_old_k
                        if sy_new_k is not None:
                            model_img[sy_new_k, sx_new_k] += stamp_new_k
                        df = abs(f_new - f[k]) / max(abs(f[k]), 1e-10)
                        changed = fit_parameters_changed(
                            x[k],
                            y[k],
                            f[k],
                            x_new,
                            y_new,
                            f_new,
                            flux_fraction=flux_conv,
                        )
                        if changed:
                            n_changed += 1
                            changed_this_iter[k] = True
                        max_df = max(max_df, df)
                        x[k], y[k], f[k], chi2[k] = x_new, y_new, f_new, chi2_k
                        fit_ok[k] = True

        last_changed = changed_this_iter
        if log_fn:
            log_fn(f"  [APEX] iter {it + 1}/{max_iter} | changed={n_changed} max_dflux={max_df:.4f}")
        if n_changed == 0 and it > 0:
            exhausted_iterations = False
            if log_fn:
                log_fn(f"  [APEX] converged at iter {it + 1}")
            break

    if exhausted_iterations and np.any(last_changed & active_mask & ~fixed_mask):
        fit_flags[last_changed & active_mask & ~fixed_mask] |= int(PSFFitFlag.NONCONVERGENCE)

    final_model = _allstar_build_model(img_sub.shape, x, y, f, eval_psf, stamp_size)
    quality = measure_psf_fit_quality(
        img_sub,
        final_model,
        x,
        y,
        f,
        eval_psf,
        fit_shape=fit_shape,
        background_rms=background_rms,
        gain=gain,
        fit_ok=fit_ok,
        initial_xy=anchors,
        xy_bound=bound if np.any(~fixed_mask) else None,
    )
    output_flags = np.asarray(quality.flags, dtype=np.int32) | fit_flags

    return _Tab({
        'x_fit': x, 'y_fit': y, 'flux_fit': f,
        'flux_err': quality.flux_err,
        'qfit': quality.qfit,
        'cfit': quality.cfit,
        'reduced_chi2': quality.reduced_chi2,
        'flags': output_flags,
        'n_pixels_fit': quality.n_pixels_fit,
        'iter_detected': np.ones(N, dtype=int),
    })


# ── PSF Worker ────────────────────────────────────────────────────────────────

class Step6PSFWorker(QThread):
    """Per-frame EPSFBuilder + PSFPhotometry worker.

    Algorithm per frame:
    1. Load detected positions from detect_{fname}.csv
    2. Select bright isolated stars for EPSF building
    3. Build oversampled EPSF with EPSFBuilder
    4. Run iterative PSFPhotometry: fit → residual → detect new → re-fit
    5. Save residual FITS and epsf_model FITS
    6. Emit per-frame result
    """
    progress = pyqtSignal(int, int, str)
    worker_status = pyqtSignal(int, str, str, int)  # worker_id, frame, stage, progress(0-100)
    frame_done = pyqtSignal(str, dict)
    epsf_ready = pyqtSignal(str, str, object)      # display_key, frame_name, epsf_array (numpy)
    residual_ready = pyqtSignal(str, object, object)  # fname, residual_meta(dict), new_xy (or None)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str, str)
    log = pyqtSignal(str)

    FLAG_SAT = int(PSFFitFlag.SATURATED)
    FLAG_EDGE = int(PSFFitFlag.INCOMPLETE_REGION)
    FLAG_FIT_FAIL = int(PSFFitFlag.NONCONVERGENCE)

    def __init__(self, file_list, params, data_dir, result_dir, cache_dir, use_cropped=False):
        super().__init__()
        self.file_list = list(file_list)
        self.params = params
        self.data_dir = Path(data_dir)
        self.result_dir = Path(result_dir)
        self.cache_dir = Path(cache_dir)
        self.use_cropped = use_cropped
        w_override = _to_int(getattr(self.params.P, "psf_parallel_workers", 0), 0)
        self.max_workers = max(1, w_override) if w_override > 0 else get_parallel_workers(params)
        self._workers_override = (w_override > 0)
        self._executor = None
        self._stop_requested = False

    def stop(self):
        self._stop_requested = True
        self._log("Stop requested — finishing running frames, cancelling queued frames.")

    def _log(self, msg):
        self.log.emit(msg)

    def _resolve_fits_path(self, fname: str) -> Path | None:
        if self.use_cropped and crop_is_active(self.result_dir):
            cdir = step2_cropped_dir(self.result_dir)
            cpath = cdir / fname
            if cpath.exists():
                return cpath
        fpath = self.data_dir / fname
        return fpath if fpath.exists() else None

    def run(self):  # noqa: C901
        try:
            try:
                from photutils.psf import EPSFBuilder, extract_stars, PSFPhotometry
                from photutils.detection import DAOStarFinder
                from photutils.background import LocalBackground, MMMBackground, Background2D, MedianBackground
                try:
                    from photutils.psf import SourceGrouper
                    _has_grouper = True
                except ImportError:
                    _has_grouper = False
                from astropy.table import Table, vstack as astropy_vstack
                import photutils as _pu
                self._log(f"photutils version: {_pu.__version__}")
            except ImportError as e:
                self.error.emit("IMPORT", f"photutils required: {e}")
                self.finished.emit({})
                return

            P = self.params.P
            output_dir = step8_psf_dir(self.result_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            GAIN = _to_float(getattr(P, "gain_e_per_adu", 1.0), 1.0)
            ZP = _to_float(getattr(P, "zp_initial", 25.0), 25.0)
            rn_e = _to_float(getattr(P, "rdnoise_e", 7.5), 7.5)
            sat_adu = _to_float(getattr(P, "saturation_adu", 60000.0), 60000.0)
            min_snr = _to_float(getattr(P, "min_snr_for_mag", 3.0), 3.0)
            fwhm_guess = _to_float(getattr(P, "fwhm_pix_guess", 6.0), 6.0)

            oversampling = _to_int(getattr(P, "psf_epsf_oversampling", 2), 2)
            epsf_size_fwhm_mult = _to_float(getattr(P, "psf_epsf_size_fwhm_mult", 4.0), 4.0)
            n_stars_max = max(0, _to_int(getattr(P, "psf_n_stars_max", 0), 0))
            isolation_mult = _to_float(getattr(P, "psf_isolation_fwhm_mult", 3.0), 3.0)
            epsf_contamination_filter = _as_bool(
                getattr(P, "psf_epsf_contamination_filter", True),
                True,
            )
            flux_scale_correction = _as_bool(
                getattr(P, "psf_flux_scale_correction", False),
                True,
            )
            flux_scale_min_snr = max(
                0.0, _to_float(getattr(P, "psf_flux_scale_min_snr", 50.0), 50.0)
            )
            flux_scale_min_stars = max(
                3, _to_int(getattr(P, "psf_flux_scale_min_stars", 8), 8)
            )
            flux_scale_min_neighbor_fwhm = max(
                0.0,
                _to_float(getattr(P, "psf_flux_scale_min_neighbor_fwhm", 4.0), 4.0),
            )
            flux_scale_max_scatter_mag = max(
                0.0,
                _to_float(getattr(P, "psf_flux_scale_max_scatter_mag", 0.10), 0.10),
            )
            flux_pct_lo = _to_float(getattr(P, "psf_flux_percentile_lo", 75.0), 75.0)
            flux_pct_hi = _to_float(getattr(P, "psf_flux_percentile_hi", 95.0), 95.0)
            fit_shape_fwhm_mult = _to_float(
                getattr(P, "psf_fit_shape_fwhm_mult", 2.4), 2.4
            )
            fit_window_mode = str(
                getattr(P, "psf_fit_window_mode", "auto")
            ).strip().lower()
            if fit_window_mode not in {"auto", "manual"}:
                fit_window_mode = "auto"
            fit_encircled_energy = min(
                0.995,
                max(
                    0.50,
                    _to_float(getattr(P, "psf_fit_encircled_energy", 0.90), 0.90),
                ),
            )
            max_iter = _to_int(getattr(P, "psf_max_iter", 2), 2)
            fitter_max_iter = max(1, _to_int(getattr(P, "psf_fitter_max_iter", 6), 6))
            redetect_sigma = _to_float(getattr(P, "psf_redetect_sigma", 3.5), 3.5)
            # EPSF star selection quality cuts (tighter than re-detection cuts)
            epsf_sharp_lo       = _to_float(getattr(P, "psf_epsf_sharp_lo",      0.3), 0.3)
            epsf_sharp_hi       = _to_float(getattr(P, "psf_epsf_sharp_hi",      0.8), 0.8)
            epsf_round_abs_max  = _to_float(getattr(P, "psf_epsf_round_abs_max", 0.5), 0.5)
            epsf_elong_max      = _to_float(getattr(P, "psf_epsf_elong_max",     1.3), 1.3)
            # IterativePSFPhotometry iteration mode: "new" (fast) or "all" (accurate, slow)
            fit_mode_cfg = str(getattr(P, "psf_fit_mode", "new")).strip().lower()
            if fit_mode_cfg not in ("new", "all"):
                fit_mode_cfg = "new"
            redetect_sharp_lo = _to_float(getattr(P, "psf_redetect_sharp_lo", 0.15), 0.15)
            redetect_sharp_hi = _to_float(getattr(P, "psf_redetect_sharp_hi", 0.95), 0.95)
            redetect_round_abs_max = _to_float(getattr(P, "psf_redetect_round_abs_max", 0.8), 0.8)
            duplicate_radius_px_cfg = _to_float(getattr(P, "psf_duplicate_radius_px", np.nan), np.nan)
            duplicate_radius_mult = _to_float(getattr(P, "psf_duplicate_radius_fwhm_mult", 0.8), 0.8)
            new_sources_cap_per_iter = _to_int(getattr(P, "psf_new_sources_cap_per_iter", 70), 70)
            new_sources_cap_frac = _to_float(getattr(P, "psf_new_sources_cap_frac", 0.02), 0.02)
            conv_new_frac = _to_float(getattr(P, "psf_conv_new_frac", 0.02), 0.02)
            flux_conv_threshold = _to_float(getattr(P, "psf_flux_conv_threshold", 0.01), 0.01)
            postfit_snr_min = max(0.0, _to_float(getattr(P, "psf_postfit_snr_min", 3.0), 3.0))
            postfit_qfit_max = max(
                0.0,
                _to_float(getattr(P, "psf_postfit_qfit_max", 3.0), 3.0),
            )
            postfit_reduced_chi2_max = max(
                0.0,
                _to_float(getattr(P, "psf_postfit_reduced_chi2_max", 25.0), 25.0),
            )
            blend_residual_ratio = min(
                1.0,
                max(0.0, _to_float(getattr(P, "psf_blend_residual_ratio", 0.3), 0.3)),
            )
            fit_init_max_sources = _to_int(getattr(P, "psf_fit_init_max_sources", 0), 0)
            core_cut_enable = bool(getattr(P, "psf_core_cut_enable", False))
            core_cut_center_mode = str(getattr(P, "psf_core_cut_center_mode", "auto")).strip().lower() or "auto"
            if core_cut_center_mode not in ("auto", "image", "manual"):
                core_cut_center_mode = "auto"
            core_cut_x_px = _to_float(getattr(P, "psf_core_cut_x_px", np.nan), np.nan)
            core_cut_y_px = _to_float(getattr(P, "psf_core_cut_y_px", np.nan), np.nan)
            core_cut_radius_px = _to_float(getattr(P, "psf_core_cut_radius_px", 0.0), 0.0)
            core_cut_radius_fwhm_mult = _to_float(getattr(P, "psf_core_cut_radius_fwhm_mult", 20.0), 20.0)
            core_cut_auto_cell_fwhm_mult = _to_float(getattr(P, "psf_core_cut_auto_cell_fwhm_mult", 8.0), 8.0)
            core_cut_auto_min_density_ratio = _to_float(getattr(P, "psf_core_cut_auto_min_density_ratio", 1.5), 1.5)
            core_cut_auto_min_sources = _to_int(getattr(P, "psf_core_cut_auto_min_sources", 50), 50)
            core_cut_max_exclude_frac = _to_float(getattr(P, "psf_core_cut_max_exclude_frac", 0.70), 0.70)
            use_error_image = bool(getattr(P, "psf_use_error_image", True))
            use_grouper = bool(getattr(P, "psf_use_grouper", True))
            grouper_max_size = _to_int(getattr(P, "psf_grouper_max_size", 3), 3)
            grouper_max_size = min(25, max(1, grouper_max_size))
            grouper_radius_fwhm = min(
                5.0,
                max(0.5, _to_float(getattr(P, "psf_grouper_radius_fwhm", 1.5), 1.5)),
            )
            forced_match_radius_fwhm = min(
                3.0,
                max(
                    0.1,
                    _to_float(
                        getattr(P, "psf_forced_match_radius_fwhm", 1.25),
                        1.25,
                    ),
                ),
            )
            save_all_iter_residuals = bool(getattr(P, "psf_save_all_iter_residuals", False))
            model_mode = str(getattr(P, "psf_model_mode", "per_frame")).strip().lower()
            max_workers = max(1, int(self.max_workers))

            if model_mode != "per_frame":
                self._log(f"PSF mode '{model_mode}' is disabled; forcing per_frame")
                model_mode = "per_frame"
            use_shared_filter_epsf = bool(getattr(P, "psf_shared_filter_epsf", False))
            min_epsf_stars = max(1, _to_int(getattr(P, "psf_min_epsf_stars", 10), 10))
            psf_fit_engine_cfg = str(
                getattr(P, "psf_fit_engine", "apex_iterative")
            ).strip().lower()
            if psf_fit_engine_cfg == "allstar":
                psf_fit_engine_cfg = "apex_iterative"
            if psf_fit_engine_cfg not in ("photutils", "apex_iterative"):
                psf_fit_engine_cfg = "apex_iterative"
            psf_build_mode_cfg = str(getattr(P, "psf_build_mode", "epsf")).strip().lower()
            if psf_build_mode_cfg not in ("epsf", "moffat", "moffat_hybrid"):
                self._log(
                    f"PSF build mode '{psf_build_mode_cfg}' is unknown; using epsf"
                )
                psf_build_mode_cfg = "epsf"
            # How the oversampled grid is sampled at sub-pixel offsets. Linear
            # has been the only behaviour; higher orders are the control for
            # whether interpolation error was costing accuracy.
            psf_interp_order = int(np.clip(
                _to_int(getattr(P, "psf_interp_order", 1), 1), 1, 5))
            if psf_interp_order != 1:
                self._log(f"PSF grid interpolation order = {psf_interp_order}")
            if psf_build_mode_cfg == "moffat_hybrid":
                self._log(
                    "PSF build mode 'moffat_hybrid': analytic Moffat evaluated in "
                    "closed form plus an interpolated residual grid — the steep "
                    "core avoids interpolation error, the leftover keeps the shape."
                )
            if psf_build_mode_cfg == "moffat":
                # An analytic Moffat is smooth by construction, so it does not
                # carry the reference stars' pixel noise into every fit the way
                # the empirical grid does. It is also circular: astropy's
                # Moffat2D has no axis ratio, and this frame's bright stars run
                # 8-12 % elongated (p90 ~1.45), which the empirical ePSF
                # reproduces and this model cannot.
                self._log(
                    "PSF build mode 'moffat': analytic model, no empirical "
                    "residual grid. Circular by construction — check elongation."
                )
            self._log(
                f"PSF engine={psf_fit_engine_cfg} | build={psf_build_mode_cfg}"
            )

            redetect_sigma = max(1.0, redetect_sigma)
            new_sources_cap_per_iter = max(0, new_sources_cap_per_iter)
            new_sources_cap_frac = min(max(0.0, new_sources_cap_frac), 1.0)
            conv_new_frac = min(max(0.0, conv_new_frac), 1.0)
            core_cut_radius_px = max(0.0, core_cut_radius_px) if np.isfinite(core_cut_radius_px) else 0.0
            core_cut_radius_fwhm_mult = max(1.0, core_cut_radius_fwhm_mult)
            core_cut_auto_cell_fwhm_mult = max(2.0, core_cut_auto_cell_fwhm_mult)
            core_cut_auto_min_density_ratio = max(1.0, core_cut_auto_min_density_ratio)
            core_cut_auto_min_sources = max(5, core_cut_auto_min_sources)
            core_cut_max_exclude_frac = min(max(0.05, core_cut_max_exclude_frac), 0.95)
            duplicate_radius_mult = max(0.0, duplicate_radius_mult)
            if np.isfinite(duplicate_radius_px_cfg):
                duplicate_radius_px_cfg = max(0.0, float(duplicate_radius_px_cfg))
            dedup_enabled = bool(
                (np.isfinite(duplicate_radius_px_cfg) and duplicate_radius_px_cfg > 0.0)
                or (duplicate_radius_mult > 0.0)
            )
            # Outdated sentinel values (-999/999) effectively disable morphology cuts and
            # can explode residual re-detections in crowded fields.
            if redetect_sharp_lo <= -900.0 and redetect_sharp_hi >= 900.0:
                redetect_sharp_lo, redetect_sharp_hi = 0.15, 0.95
            if redetect_round_abs_max >= 9.0:
                redetect_round_abs_max = 0.8
            # If user/state has extremely loose residual cuts (e.g. sharp=[0,1], round=2),
            # tighten them to suppress ring/halo false detections.
            if redetect_sharp_lo <= 0.01 and redetect_sharp_hi >= 0.99 and redetect_round_abs_max >= 1.5:
                redetect_sharp_lo, redetect_sharp_hi, redetect_round_abs_max = 0.15, 0.95, 0.8

            self._log(
                "PSF settings | "
                f"model_mode={model_mode} | fit_mode={fit_mode_cfg} | "
                f"max_iter={max_iter} | redetect_sigma={redetect_sigma:.2f} | "
                f"cap_iter={new_sources_cap_per_iter} | cap_frac={new_sources_cap_frac:.3f} | "
                f"use_error_image={'on' if use_error_image else 'off'} | "
                f"use_grouper={'on' if use_grouper else 'off'}"
            )
            self._log(
                f"PSF redetect cuts | sharp=[{redetect_sharp_lo:.2f},{redetect_sharp_hi:.2f}] "
                f"| |round|<={redetect_round_abs_max:.2f}"
            )
            if np.isfinite(duplicate_radius_px_cfg):
                self._log(f"PSF dedup radius: {duplicate_radius_px_cfg:.2f}px (absolute)")
            else:
                self._log(f"PSF dedup radius: {duplicate_radius_mult:.2f}xFWHM")
            if not use_grouper:
                self._log("PSF fit mode: iterative 'new' (grouper off; photutils 2.3 requires grouper for mode='all')")
            self._log(
                "PSF core cut | "
                f"{'on' if core_cut_enable else 'off'} | center={core_cut_center_mode} | "
                f"radius_px={core_cut_radius_px:.1f} "
                f"fallback={core_cut_radius_fwhm_mult:.1f}xFWHM | "
                f"density_ratio>={core_cut_auto_min_density_ratio:.2f}"
            )
            self._log(
                f"PSF scales | epsf_cutout={epsf_size_fwhm_mult:.2f}xFWHM | "
                f"fit_window={fit_shape_fwhm_mult:.2f}xFWHM | "
                "subtract_window≈2xEPSF"
            )

            self._log(
                "EPSF contamination-aware references | "
                f"{'on' if epsf_contamination_filter else 'off'}"
            )

            frames = list(self.file_list)

            use_qc = should_use_frame_quality_qc(
                self.result_dir,
                self.params.P,
                "phot_use_qc_pass_only",
                default=False,
            )
            frames, qc_info = filter_files_by_qc(self.result_dir, frames, require_qc=use_qc)
            if use_qc:
                if qc_info.get("applied"):
                    self._log(f"Step4 QC: {qc_info['kept']}/{qc_info['total']} frame(s) kept.")
                elif qc_info.get("path") is None:
                    self._log("Step4 QC: frame_quality.csv not found; using all frames.")
                else:
                    self._log(f"Step4 QC: frame_quality.csv ignored ({qc_info['reason']}); using all frames.")
            if not frames:
                raise RuntimeError("No frames remain after Step 4 QC filtering.")

            total = len(frames)
            index_rows = []
            counters = {"processed": 0, "no_detect": 0, "no_fits": 0, "stopped": 0}
            completed = [0]
            epsf_cache: dict[str, object] = {}  # filter → epsf model
            epsf_cache_lock = Lock()
            if self._workers_override:
                self._log(f"PSF parallel workers={max_workers} (Step6 override)")
            else:
                self._log(f"PSF parallel workers={max_workers}")

            run_t0 = time.time()
            last_hb = 0.0
            last_stall_log = 0.0
            last_done_count = -1

            def _fmt_eta(sec: float) -> str:
                s = int(max(0, round(float(sec))))
                h, rem = divmod(s, 3600)
                m, ss = divmod(rem, 60)
                if h > 0:
                    return f"{h:d}:{m:02d}:{ss:02d}"
                return f"{m:02d}:{ss:02d}"

            def process_single_frame(fname: str):
                if self._stop_requested:
                    self.worker_status.emit(0, fname, "Stopped", 100)
                    return {"file": fname, "status": "stopped"}

                wid = int(threading.get_ident() % 10000)
                _t_frame = time.time()
                _t: dict[str, float] = {"start": _t_frame}
                self.progress.emit(completed[0], total, f"RUN | {fname}")
                self.worker_status.emit(wid, fname, "Load", 5)
                img_path = self._resolve_fits_path(fname)
                if img_path is None:
                    self.worker_status.emit(wid, fname, "No FITS", 100)
                    return {"file": fname, "status": "no_fits", "reason": "no FITS"}

                det_df = _load_detect_positions(fname, self.cache_dir, self.result_dir)
                if det_df is None or len(det_df) == 0:
                    self.worker_status.emit(wid, fname, "No detect", 100)
                    return {"file": fname, "status": "no_detect", "reason": "no detect csv"}

                try:
                    img = fits.getdata(img_path).astype(np.float32, copy=False)
                except Exception as e:
                    self.worker_status.emit(wid, fname, "FITS error", 100)
                    return {"file": fname, "status": "no_fits", "reason": f"FITS read: {e}"}

                try:
                    header = fits.getheader(img_path)
                except Exception:
                    header = None
                noise = resolve_effective_noise_params(P, header)
                GAIN = noise.gain_e_per_adu
                rn_e = noise.rdnoise_e
                noise_info = {
                    "gain_e_per_adu": float(noise.gain_e_per_adu),
                    "rdnoise_e": float(noise.rdnoise_e),
                    "binning_x": int(noise.bin_x),
                    "binning_y": int(noise.bin_y),
                    "gain_source": noise.gain_source,
                    "rdnoise_source": noise.rdnoise_source,
                }

                this_filter = _get_filter_lower(img_path)
                exptime = _get_exptime(img_path, default=1.0)
                fwhm_med = _load_fwhm_from_meta(fname, self.cache_dir, self.result_dir, fwhm_guess)
                fwhm_safe = max(float(fwhm_med), 1.0)
                pixel_scale_arcsec = _to_float(
                    getattr(P, "pixel_scale_arcsec", np.nan), np.nan
                )
                fwhm_arcsec = (
                    float(fwhm_med) * pixel_scale_arcsec
                    if np.isfinite(pixel_scale_arcsec) and pixel_scale_arcsec > 0
                    else np.nan
                )
                fwhm_qc_max_px = _to_float(
                    getattr(P, "fwhm_px_max", np.nan), np.nan
                )
                # Size controls are driven primarily by FWHM multipliers for per-frame adaptation.
                epsf_size_frame = _odd_int(
                    float(epsf_size_fwhm_mult) * fwhm_safe,
                    min_value=25,
                    max_value=101,
                )
                _epsf_desired = int(round(float(epsf_size_fwhm_mult) * fwhm_safe))
                if _epsf_desired > 101:
                    self._log(
                        f"  [PSF] epsf_size capped to max: desired={_epsf_desired}px → {epsf_size_frame}px "
                        f"(fwhm={fwhm_safe:.1f}px, mult={epsf_size_fwhm_mult:.2f}x)"
                    )
                fit_shape_frame = _odd_int(
                    float(fit_shape_fwhm_mult) * fwhm_safe,
                    min_value=9,
                    max_value=31,
                )
                if fit_shape_frame >= epsf_size_frame:
                    fit_shape_frame = _odd_int(max(3, epsf_size_frame - 4), min_value=3, max_value=31)
                render_shape_frame = _odd_int(
                    max(float(epsf_size_frame) * 2.0, float(fit_shape_frame)),
                    min_value=11,
                    max_value=201,
                )

                epsf_cache_key = this_filter if use_shared_filter_epsf else f"{this_filter}:{fname}"
                h, w = img.shape
                det_xy_all_for_core = det_df[["x", "y"]].to_numpy(float)
                core_center_mode_for_estimate = core_cut_center_mode
                core_manual_center = (core_cut_x_px, core_cut_y_px)
                core_center_method = ""
                if str(core_cut_center_mode).strip().lower() == "auto":
                    target_pixel = target_pixel_from_wcs(
                        header,
                        getattr(P, "target_ra_deg", None),
                        getattr(P, "target_dec_deg", None),
                        img.shape,
                    )
                    if target_pixel is not None:
                        core_center_mode_for_estimate = "manual"
                        core_manual_center = target_pixel
                        core_center_method = "target_wcs"
                core_diagnostic = estimate_psf_core_cut(
                    det_xy_all_for_core,
                    img.shape,
                    fwhm_safe,
                    enabled=True,
                    center_mode=core_center_mode_for_estimate,
                    manual_center=core_manual_center,
                    radius_px=core_cut_radius_px,
                    radius_fwhm_mult=core_cut_radius_fwhm_mult,
                    auto_cell_fwhm_mult=core_cut_auto_cell_fwhm_mult,
                    auto_min_density_ratio=core_cut_auto_min_density_ratio,
                    auto_min_sources=core_cut_auto_min_sources,
                    max_exclude_frac=core_cut_max_exclude_frac,
                    image=img,
                )
                if core_center_method:
                    core_diagnostic = replace(
                        core_diagnostic,
                        method=core_diagnostic.method.replace("manual", core_center_method, 1),
                    )
                if core_cut_enable:
                    core_cut = core_diagnostic
                else:
                    core_cut = PSFCoreCut(
                        False,
                        center_x=core_diagnostic.center_x,
                        center_y=core_diagnostic.center_y,
                        radius_px=core_diagnostic.radius_px,
                        method=core_diagnostic.method,
                        n_total=core_diagnostic.n_total,
                        n_excluded=core_diagnostic.n_excluded,
                        n_kept=core_diagnostic.n_kept,
                        density_ratio=core_diagnostic.density_ratio,
                        reason="disabled",
                    )
                n_core_excluded_init = 0
                n_core_excluded_redetect = 0
                n_core_excluded_result = 0
                # EPSF 품질 지표 기본값 — 모델 재사용 등으로 품질 검사 블록을
                # 건너뛰는 경로에서도 residual_meta 조립이 참조할 수 있게 한다.
                epsf_quality_n_blobs = 0
                epsf_quality_max_quadrant_frac = float("nan")

                def _core_keep(xy_like) -> np.ndarray:
                    return psf_core_keep_mask(np.asarray(xy_like, dtype=float), core_cut)

                if core_cut.enabled:
                    self._log(
                        f"  [CORE] enabled | center=({core_cut.center_x:.1f},{core_cut.center_y:.1f}) "
                        f"r={core_cut.radius_px:.1f}px | method={core_cut.method} | "
                        f"exclude={core_cut.n_excluded}/{core_cut.n_total}"
                    )
                elif core_cut_enable:
                    self._log(
                        f"  [CORE] auto cut off for {fname}: {core_cut.reason or 'not_applicable'}"
                    )

                try:
                    self.worker_status.emit(wid, fname, "Background", 20)
                    _t["bkg"] = time.time()
                    _box = max(32, min(128, h // 16, w // 16))
                    try:
                        from astropy.stats import SigmaClip as _SigmaClip
                        _bkg2d = Background2D(img, (_box, _box), filter_size=(3, 3),
                                              sigma_clip=_SigmaClip(sigma=3.0, maxiters=3),
                                              bkg_estimator=MedianBackground())
                        bkg_map = np.asarray(_bkg2d.background, dtype=np.float32)
                        bkg_rms_scalar = float(_bkg2d.background_rms_median)
                        bkg_med = float(_bkg2d.background_median)
                        bkg_std = float(bkg_rms_scalar)
                        img_sub = (img - bkg_map).astype(np.float32, copy=False)
                        del bkg_map
                        self._log(f"  [BKG] Background2D | box={_box}px | "
                                  f"bkg_med={bkg_med:.2f} rms={bkg_std:.3f}")
                    except Exception as _bkg_e:
                        # Fallback to scalar sigma-clipped stats
                        self._log(f"  [BKG] Background2D failed ({_bkg_e}); using scalar median")
                        _, bkg_med, bkg_std = sigma_clipped_stats(img, sigma=3.0, maxiters=3)
                        bkg_rms_scalar = float(bkg_std)
                        img_sub = (img - float(bkg_med)).astype(np.float32, copy=False)

                    if self._stop_requested:
                        self.worker_status.emit(wid, fname, "Stopped", 100)
                        return {"file": fname, "status": "stopped", "reason": "stop requested"}

                    _t["bkg_done"] = time.time()
                    epsf_emit_arr = None
                    n_epsf_detected = 0
                    n_epsf_candidates = 0
                    n_epsf_candidates_pre_morph = 0
                    n_epsf_candidates_post_morph = 0
                    n_epsf_selected = 0
                    n_iso = 0
                    epsf_plan_target = 0
                    epsf_grid_size = 1
                    n_epsf_low_contamination = 0
                    n_epsf_core_rejected = 0
                    n_epsf_fallback_selected = 0
                    n_epsf_morphology_relaxed_selected = 0
                    epsf_selected_median_contamination = np.nan
                    epsf_reference_catalog_name = ""
                    psf_type_built = psf_build_mode_cfg  # 'epsf' or 'moffat'
                    with epsf_cache_lock:
                        epsf_model = epsf_cache.get(epsf_cache_key)
                    if epsf_model is None:
                        self.worker_status.emit(wid, fname, "PSF build", 40)
                        _t["epsf"] = time.time()
                        xy_all = det_df[["x", "y"]].to_numpy(float)
                        finite_xy = np.isfinite(xy_all[:, 0]) & np.isfinite(xy_all[:, 1])
                        xy_all = xy_all[finite_xy]
                        if len(xy_all) < 5:
                            raise RuntimeError("Too few detected sources for EPSF building")

                        if "flux_init" in det_df.columns:
                            fluxes = det_df["flux_init"].to_numpy(float)[finite_xy]
                        else:
                            xi = xy_all[:, 0].astype(int).clip(0, w - 1)
                            yi = xy_all[:, 1].astype(int).clip(0, h - 1)
                            fluxes = img_sub[yi, xi]

                        valid_flux = np.isfinite(fluxes)
                        lo = np.nanpercentile(fluxes[valid_flux], flux_pct_lo) if np.any(valid_flux) else -np.inf
                        hi = np.nanpercentile(fluxes[valid_flux], flux_pct_hi) if np.any(valid_flux) else np.inf
                        in_range = valid_flux & (fluxes >= lo) & (fluxes <= hi)

                        peak_vals = img[
                            xy_all[:, 1].astype(int).clip(0, h - 1),
                            xy_all[:, 0].astype(int).clip(0, w - 1),
                        ]
                        not_sat = peak_vals < sat_adu
                        in_range = in_range & not_sat

                        epsf_half = epsf_size_frame // 2 + 5
                        not_edge = (
                            (xy_all[:, 0] >= epsf_half) & (xy_all[:, 0] <= w - 1 - epsf_half) &
                            (xy_all[:, 1] >= epsf_half) & (xy_all[:, 1] <= h - 1 - epsf_half)
                        )
                        in_range = in_range & not_edge
                        if core_cut.enabled:
                            epsf_core_keep = _core_keep(xy_all)
                            n_epsf_core_drop = int(np.sum(in_range & ~epsf_core_keep))
                            if n_epsf_core_drop > 0:
                                self._log(
                                    f"[EPSF] core cut removed {n_epsf_core_drop} candidates "
                                    f"(r<{core_cut.radius_px:.1f}px)"
                                )
                            in_range = in_range & epsf_core_keep

                        if "epsf_candidate" in det_df.columns:
                            epsf_candidate = det_df["epsf_candidate"].to_numpy(bool)[finite_xy]
                            cand_range = in_range & epsf_candidate
                            n_before = int(np.sum(in_range))
                            n_after = int(np.sum(cand_range))
                            if n_after >= min_epsf_stars:
                                in_range = cand_range
                                self._log(
                                    f"[EPSF] Step4 epsf_candidate filter: {n_before} -> {n_after}"
                                )
                            else:
                                self._log(
                                    f"[WARN][EPSF] Step4 epsf_candidate left {n_after} stars; "
                                    f"fallback to local EPSF cuts ({n_before})"
                                )

                        # Treat morphology as a preference so a narrow strict pool
                        # cannot force contaminated stars into the ePSF model.
                        _in_range_pre_morph = in_range.copy()
                        _morphology_ok = np.ones(len(xy_all), dtype=bool)
                        _morph_applied = False
                        if "sharpness" in det_df.columns:
                            _sharp = det_df["sharpness"].to_numpy(float)[finite_xy]
                            _morphology_ok &= np.isfinite(_sharp) & (_sharp >= epsf_sharp_lo) & (_sharp <= epsf_sharp_hi)
                            _morph_applied = True
                        if "roundness" in det_df.columns:
                            _round = det_df["roundness"].to_numpy(float)[finite_xy]
                            _morphology_ok &= np.isfinite(_round) & (np.abs(_round) <= epsf_round_abs_max)
                            _morph_applied = True
                        if "elong" in det_df.columns:
                            _elong = det_df["elong"].to_numpy(float)[finite_xy]
                            _morphology_ok &= np.isfinite(_elong) & (_elong <= epsf_elong_max)
                            _morph_applied = True
                        n_epsf_detected = len(xy_all)
                        n_epsf_candidates_pre_morph = int(np.sum(_in_range_pre_morph))
                        n_epsf_candidates_post_morph = int(
                            np.sum(_in_range_pre_morph & _morphology_ok)
                        )
                        epsf_plan = plan_epsf_stars(
                            n_epsf_detected,
                            n_epsf_candidates_pre_morph,
                            user_cap=n_stars_max,
                        )
                        epsf_plan_target = epsf_plan.target
                        epsf_grid_size = epsf_plan.grid_size
                        # 완화하더라도 절대 넘기면 안 되는 최소 방어선: ePSF 참조별의
                        # FWHM 은 프레임 대표 FWHM 의 절반은 되어야 한다. 점광원 잡음
                        # (우주선·핫픽셀)은 FWHM 이 1 px 근처라 여기서 걸린다.
                        # 이 하한이 없으면, 형태 컷이 목표 수를 못 채웠을 때 컷을 통째로
                        # 버리고 「밝고 고립」 기준만 남는데 — 그 기준의 최적해가 바로
                        # 우주선이다(플럭스가 1 px 에 몰려 피크 최대, 이웃 없음).
                        # 2026-07-29 M67/QHY600 에서 실제로 FWHM 1.00 px 소스가 참조별로
                        # 뽑혀 ePSF 가 2.75배 좁아졌고 PSF 플럭스가 구경의 32% 로 떨어졌다.
                        _fwhm_floor_ok = np.ones(len(xy_all), dtype=bool)
                        _fwhm_floor = 0.5 * float(fwhm_safe)
                        if "fwhm_px" in det_df.columns:
                            _fw = det_df["fwhm_px"].to_numpy(float)[finite_xy]
                            # 측정 실패(NaN)는 통과시킨다 — 하한은 「점광원임이 확인된
                            # 것」만 거르기 위한 것이지 미측정 별을 버리기 위한 게 아니다.
                            _fwhm_floor_ok = ~np.isfinite(_fw) | (_fw >= _fwhm_floor)

                        # FWHM 은 표본(measure_max)에만 측정되므로 대부분 NaN 이다.
                        # 픽셀만으로 판정하는 PSF 대칭 검사를 함께 건다 — 별은 등방이라
                        # 좌우 또는 상하 **양쪽** 이웃이 피크의 일정 비율 이상이지만,
                        # 고립 스파이크·2픽셀 쌍·L자/대각 클러스터는 그렇지 못하다.
                        # 실측(M67/QHY600, CR 미제거 프레임): astroscrappy 가 CR 로 지목한
                        # 45,433곳 중 67.3% 를 걸러내면서 검출된 별 736개는 100% 보존.
                        _psf_sym_ok = psf_symmetric_mask(
                            img_sub, xy_all, background=0.0, neighbor_frac=0.3
                        )
                        _fwhm_floor_ok = _fwhm_floor_ok & _psf_sym_ok

                        if _morph_applied:
                            if n_epsf_candidates_post_morph < epsf_plan.target:
                                in_range = _in_range_pre_morph & _fwhm_floor_ok
                                n_floor_drop = int(
                                    np.sum(_in_range_pre_morph & ~_fwhm_floor_ok)
                                )
                                _morphology_ok_cand = _morphology_ok[in_range]
                                self._log(
                                    f"[EPSF] morphology filter: {n_epsf_candidates_pre_morph} -> "
                                    f"{n_epsf_candidates_post_morph}; target={epsf_plan.target} "
                                    f"-> relaxed to pre-morph pool "
                                    f"(FWHM>={_fwhm_floor:.2f}px floor dropped {n_floor_drop})"
                                )
                            else:
                                in_range = _in_range_pre_morph & _morphology_ok
                                _morphology_ok_cand = np.ones(
                                    int(np.sum(in_range)), dtype=bool
                                )
                                self._log(
                                    f"[EPSF] morphology filter: {n_epsf_candidates_pre_morph} -> "
                                    f"{n_epsf_candidates_post_morph} "
                                    f"(sharp=[{epsf_sharp_lo:.2f},{epsf_sharp_hi:.2f}] "
                                    f"|round|<={epsf_round_abs_max:.2f} "
                                    f"elong<={epsf_elong_max:.2f})"
                                )
                        else:
                            _morphology_ok_cand = np.ones(
                                n_epsf_candidates_pre_morph, dtype=bool
                            )
                        xy_cand = xy_all[in_range]
                        fluxes_cand = fluxes[in_range]
                        if len(xy_cand) < 3:
                            raise RuntimeError("Too few candidates after flux/sat/edge filter")

                        n_epsf_candidates = len(xy_cand)

                        if len(xy_cand) >= 2:
                            tree = cKDTree(xy_cand)
                            nn_dists, _ = tree.query(xy_cand, k=min(2, len(xy_cand)), workers=1)
                            nn_d = nn_dists[:, 1] if nn_dists.ndim > 1 else nn_dists
                            isolated = nn_d > isolation_mult * fwhm_med
                            n_iso = int(np.count_nonzero(isolated))
                            if np.any(isolated):
                                xy_iso = xy_cand[isolated]
                                fl_iso = fluxes_cand[isolated]
                                if not epsf_contamination_filter:
                                    self._log(
                                        f"[EPSF] isolate pass | frame={fname} | cand={len(xy_cand)} | "
                                        f"isolated={n_iso} (thr={isolation_mult:.2f}xFWHM)"
                                    )
                            else:
                                xy_iso = xy_cand
                                fl_iso = fluxes_cand
                                # P4-9: isolation fallback → WARN level (EPSF quality degraded)
                                if not epsf_contamination_filter:
                                    self._log(
                                        f"[WARN][EPSF] isolated=0 for {fname} | "
                                        f"falling back to {len(xy_cand)} non-isolated candidates. "
                                        f"EPSF may be contaminated by neighbours. "
                                        f"Consider lowering isolation_fwhm_mult (current={isolation_mult:.1f}) "
                                        f"or using a less crowded frame."
                                    )
                                    self.log.emit(
                                        f"⚠ EPSF isolation fallback [{fname}]: "
                                        f"no isolated stars (thr={isolation_mult:.1f}×FWHM). "
                                        f"PSF model may be degraded — check log."
                                    )
                        else:
                            xy_iso = xy_cand
                            fl_iso = fluxes_cand
                            n_iso = len(xy_cand)

                        flux_order = np.argsort(fluxes_cand, kind="stable")
                        flux_rank = np.empty(len(fluxes_cand), dtype=float)
                        flux_rank[flux_order] = np.linspace(0.0, 1.0, len(fluxes_cand))
                        separation_score = np.clip(
                            nn_d / max(isolation_mult * fwhm_med, 1.0),
                            0.0,
                            2.0,
                        )
                        quality_score = 4.0 * separation_score + flux_rank

                        selected_indices: list[int] = []
                        isolated_indices = np.flatnonzero(isolated)
                        if isolated_indices.size:
                            selected_iso = select_spatially_balanced(
                                xy_cand[isolated_indices],
                                quality_score[isolated_indices],
                                target=min(epsf_plan.target, isolated_indices.size),
                                image_shape=img.shape,
                                grid_size=epsf_plan.grid_size,
                            )
                            selected_indices.extend(isolated_indices[selected_iso].tolist())

                        n_supplement = epsf_plan.target - len(selected_indices)
                        if n_supplement > 0:
                            remaining = np.setdiff1d(
                                np.arange(len(xy_cand), dtype=int),
                                np.asarray(selected_indices, dtype=int),
                                assume_unique=False,
                            )
                            selected_extra = select_spatially_balanced(
                                xy_cand[remaining],
                                quality_score[remaining],
                                target=n_supplement,
                                image_shape=img.shape,
                                grid_size=epsf_plan.grid_size,
                            )
                            selected_indices.extend(remaining[selected_extra].tolist())

                        selected = np.asarray(selected_indices, dtype=int)
                        contamination_score = np.full(len(xy_cand), np.nan, dtype=float)
                        low_contamination = np.ones(len(xy_cand), dtype=bool)
                        core_safe = np.ones(len(xy_cand), dtype=bool)
                        core_distance = np.full(len(xy_cand), np.nan, dtype=float)
                        selection_tier = np.full(len(xy_cand), -1, dtype=int)
                        selection_tier[selected] = 0
                        n_epsf_low_contamination = len(xy_cand)
                        if epsf_contamination_filter:
                            reference_selection = select_epsf_reference_stars(
                                xy_cand,
                                fluxes_cand,
                                xy_all,
                                img_sub,
                                target=epsf_plan.target,
                                image_shape=img.shape,
                                grid_size=epsf_plan.grid_size,
                                fwhm_px=fwhm_med,
                                isolation_fwhm_mult=isolation_mult,
                                background_rms=bkg_rms_scalar,
                                core_center=(core_diagnostic.center_x, core_diagnostic.center_y),
                                core_radius_px=core_diagnostic.radius_px,
                                minimum_required=min_epsf_stars,
                                morphology_ok=_morphology_ok_cand,
                            )
                            selected = reference_selection.selected_indices
                            nn_d = reference_selection.nearest_neighbor_px
                            isolated = reference_selection.isolated
                            quality_score = reference_selection.quality_score
                            contamination_score = reference_selection.contamination_score
                            low_contamination = reference_selection.low_contamination
                            core_safe = reference_selection.core_safe
                            core_distance = reference_selection.core_distance_px
                            selection_tier = reference_selection.selection_tier
                            n_iso = reference_selection.n_isolated
                            n_epsf_low_contamination = reference_selection.n_low_contamination
                            n_epsf_core_rejected = reference_selection.n_core_rejected
                            n_epsf_fallback_selected = reference_selection.n_fallback_selected
                            n_epsf_morphology_relaxed_selected = (
                                reference_selection.n_morphology_relaxed_selected
                            )
                            selected_contamination = contamination_score[selected]
                            selected_contamination = selected_contamination[
                                np.isfinite(selected_contamination)
                            ]
                            if selected_contamination.size:
                                epsf_selected_median_contamination = float(
                                    np.median(selected_contamination)
                                )
                            contamination_level = (
                                "WARN][EPSF" if n_epsf_fallback_selected > 0 else "EPSF"
                            )
                            self._log(
                                f"[{contamination_level}] contamination filter | frame={fname} "
                                f"isolated(all detections)={n_iso}/{len(xy_cand)} "
                                f"low_contam={n_epsf_low_contamination} "
                                f"inside_core={n_epsf_core_rejected} "
                                f"fallback_selected={n_epsf_fallback_selected}"
                            )
                        xy_iso = xy_cand[selected]
                        fl_iso = fluxes_cand[selected]
                        n_epsf_selected = len(xy_iso)
                        n_epsf_morphology_relaxed_selected = int(
                            np.count_nonzero(~_morphology_ok_cand[selected])
                        )
                        selected_mask = np.zeros(len(xy_cand), dtype=bool)
                        selected_mask[selected] = True
                        epsf_reference_catalog_name = f"epsf_reference_{fname}.csv"
                        pd.DataFrame({
                            "x": xy_cand[:, 0],
                            "y": xy_cand[:, 1],
                            "flux": fluxes_cand,
                            "morphology_ok": _morphology_ok_cand,
                            "morphology_relaxed_selected": (
                                (~_morphology_ok_cand) & selected_mask
                            ),
                            "nearest_neighbor_px": nn_d,
                            "nearest_neighbor_fwhm": nn_d / max(fwhm_med, 1.0),
                            "contamination_score": contamination_score,
                            "core_distance_px": core_distance,
                            "core_safe": core_safe,
                            "isolated": isolated,
                            "low_contamination": low_contamination,
                            "quality_score": quality_score,
                            "selected": selected_mask,
                            "selection_tier": selection_tier,
                        }).to_csv(
                            output_dir / epsf_reference_catalog_name,
                            index=False,
                        )
                        cap_label = str(n_stars_max) if n_stars_max > 0 else "auto"
                        self._log(
                            f"[EPSF] reference plan | detected={n_epsf_detected} "
                            f"cand={n_epsf_candidates} pre={n_epsf_candidates_pre_morph} "
                            f"post={n_epsf_candidates_post_morph} isolated={n_iso} "
                            f"morph_relaxed={n_epsf_morphology_relaxed_selected} "
                            f"selected={n_epsf_selected}/{epsf_plan.target} "
                            f"grid={epsf_plan.grid_size}x{epsf_plan.grid_size} cap={cap_label}"
                        )
                        if n_iso < min_epsf_stars:
                            self._log(
                                f"[WARN][EPSF] only {n_iso} strictly isolated stars; "
                                "supplemented with the best separated candidates after morphology cuts"
                            )

                        from astropy.table import Table as AstropyTable
                        star_table = AstropyTable({"x": xy_iso[:, 0], "y": xy_iso[:, 1]})
                        nddata = NDData(data=img_sub)
                        stars_extracted = extract_stars(nddata, star_table, size=epsf_size_frame)
                        if len(stars_extracted) < 3:
                            raise RuntimeError(f"Only {len(stars_extracted)} stars extracted; need ≥3")

                        epsf_maxiters = _to_int(getattr(P, "psf_epsf_maxiters", 5), 5)
                        builder = EPSFBuilder(
                            oversampling=oversampling,
                            maxiters=max(3, epsf_maxiters),
                            progress_bar=False,
                            smoothing_kernel="quadratic",
                        )
                        if self._stop_requested:
                            self.worker_status.emit(wid, fname, "Stopped", 100)
                            return {"file": fname, "status": "stopped", "reason": "stop requested"}

                        if psf_build_mode_cfg == 'moffat':
                            epsf, n_moffat_good = _build_moffat_psf(
                                img_sub, xy_iso, fwhm_safe, epsf_size_frame, self._log)
                            psf_type_built = 'moffat'
                            self._log(
                                f"[PSF] Moffat build | filter={this_filter} "
                                f"n_stars={n_moffat_good} fwhm_guess={fwhm_safe:.2f}px"
                            )
                        elif psf_build_mode_cfg == 'moffat_hybrid':
                            # Both pieces come from the same reference stars: the
                            # ePSF is their empirical average, the Moffat is fitted
                            # to the same cutouts, and the residual is the
                            # difference. No extra star-selection policy enters.
                            _epsf_emp, _ = builder(stars_extracted)
                            _analytic, n_moffat_good = _build_moffat_psf(
                                img_sub, xy_iso, fwhm_safe, epsf_size_frame, self._log)
                            epsf = build_moffat_hybrid_psf(
                                _epsf_emp, _analytic, oversampling)
                            psf_type_built = 'moffat_hybrid'
                            _res_frac = (
                                np.abs(epsf.residual).sum()
                                / max(np.abs(epsf.data).sum(), 1e-30)
                            )
                            self._log(
                                f"[PSF] Moffat+residual build | filter={this_filter} "
                                f"n_stars={n_moffat_good} residual={_res_frac * 100:.1f}% "
                                f"of model"
                            )
                        else:
                            epsf, _ = builder(stars_extracted)
                            psf_type_built = 'epsf'

                        # ── P4-10: EPSF quality check ─────────────────────────────────
                        # 경고를 로그로만 남기지 않고 residual_meta 에도 영속화한다
                        # (M13 20250308 0002-R: 이중 피크 프레임이 로그로만 경고돼
                        # 다운스트림이 기계적으로 제외할 수 없었다 — 2026-07-29).
                        epsf_quality_n_blobs = 0
                        epsf_quality_max_quadrant_frac = float("nan")
                        try:
                            _ed = np.asarray(epsf.data, dtype=float)
                            _ed_pos = np.where(_ed > 0, _ed, 0.0)
                            _peak = float(np.nanmax(_ed_pos)) if _ed_pos.size else 0.0
                            if _peak > 0:
                                _norm = _ed_pos / _peak
                                # Double-peak check: count pixels > 50% of peak
                                # For a clean PSF, these should form one connected blob
                                _high = (_norm > 0.5).astype(float)
                                from scipy.ndimage import label as _label
                                _labeled, _n_blobs = _label(_high)
                                epsf_quality_n_blobs = int(_n_blobs)
                                if _n_blobs > 1:
                                    self._log(
                                        f"[WARN][EPSF] {fname}: possible double-peak PSF "
                                        f"({_n_blobs} blobs above 50% peak). "
                                        f"Check focus/tracking."
                                    )
                                # Asymmetry check: compare quadrant sums
                                _cy, _cx = np.array(_ed.shape) // 2
                                _q1 = float(_ed[:_cy, :_cx].sum())
                                _q2 = float(_ed[:_cy, _cx:].sum())
                                _q3 = float(_ed[_cy:, :_cx].sum())
                                _q4 = float(_ed[_cy:, _cx:].sum())
                                _qtot = _q1 + _q2 + _q3 + _q4
                                if _qtot > 0:
                                    _qmax = max(_q1, _q2, _q3, _q4) / _qtot
                                    epsf_quality_max_quadrant_frac = float(_qmax)
                                    if _qmax > 0.45:  # >45% in one quadrant = asymmetric
                                        self._log(
                                            f"[WARN][EPSF] {fname}: asymmetric PSF "
                                            f"(max quadrant fraction={_qmax:.2f}). "
                                            f"Possible tracking drift or coma."
                                        )
                        except Exception:
                            pass
                        # ─────────────────────────────────────────────────────────────

                        # ── IRAF SUBSTAR-style iterative PSF rebuild ──────────────────
                        # For each rebuild pass:
                        #   1. Render model of ALL detected sources (rough flux)
                        #   2. Render model of PSF-selection stars only
                        #   3. cleaned = img_sub - all_model + psf_star_model
                        #      → each PSF star cutout is free of neighbours
                        #   4. Re-extract cutouts from cleaned image → rebuild EPSF
                        n_substar_iters = _to_int(getattr(P, "psf_substar_iters", 1), 1)
                        substar_neighbor_r_mult = _to_float(getattr(P, "psf_substar_neighbor_r_fwhm_mult", 8.0), 8.0)
                        substar_max_sources = _to_int(getattr(P, "psf_substar_max_sources", 1500), 1500)
                        _t["substar"] = time.time()
                        _cleaned_img = None
                        _nddata_clean = None
                        _stars_clean = None
                        if n_substar_iters > 0 and len(xy_all) > 1:
                            try:
                                from astropy.nddata import NDData as _NDData

                                # Use the detection-stage flux already computed above
                                # (from det_df["flux_init"] or peak-pixel fallback).
                                # This avoids re-running aperture photometry and is
                                # consistent with the positions in xy_all.
                                _all_flux = np.where(
                                    np.isfinite(fluxes) & (fluxes > 0),
                                    fluxes,
                                    0.0,
                                )
                                # Speed optimization:
                                # substar neighbor-cleaning needs sources affecting PSF stars,
                                # not necessarily every detection in the frame.
                                _neighbor_r = max(
                                    float(substar_neighbor_r_mult) * float(fwhm_safe),
                                    float(epsf_size_frame),
                                )
                                _src_tree = cKDTree(np.asarray(xy_all, dtype=float))
                                _neighbor_set = set()
                                for _px, _py in np.asarray(xy_iso, dtype=float):
                                    _hits = _src_tree.query_ball_point([float(_px), float(_py)], r=float(_neighbor_r))
                                    _neighbor_set.update(int(h) for h in _hits)
                                if _neighbor_set:
                                    _idx_nei = np.array(sorted(_neighbor_set), dtype=int)
                                else:
                                    _idx_nei = np.arange(len(xy_all), dtype=int)

                                if substar_max_sources > 0 and len(_idx_nei) > substar_max_sources:
                                    _fsel = np.asarray(_all_flux[_idx_nei], dtype=float)
                                    _fsel = np.where(np.isfinite(_fsel), _fsel, -np.inf)
                                    _ord = np.argsort(_fsel)[::-1][:int(substar_max_sources)]
                                    _idx_nei = _idx_nei[_ord]

                                _xy_sub = np.asarray(xy_all[_idx_nei], dtype=float)
                                _all_flux_sub = np.asarray(_all_flux[_idx_nei], dtype=float)
                                _psf_nn_tree = cKDTree(_xy_sub) if len(_xy_sub) else None
                                self._log(
                                    f"[EPSF] substar sources | frame={fname} | "
                                    f"all={len(xy_all)} near_psf={len(_xy_sub)} "
                                    f"(r={_neighbor_r:.1f}px, cap={substar_max_sources})"
                                )

                                _render_sz = int(epsf_size_frame)
                                _rough_epsf = epsf

                                for _si in range(n_substar_iters):
                                    # Full source model (all detected)
                                    _rough_eval = _make_psf_evaluator(
                                        _rough_epsf, psf_type_built, oversampling,
                                        psf_interp_order
                                    )

                                    # PSF-star-only model (add back after subtraction)
                                    _psf_flux = np.zeros(len(xy_iso), dtype=float)
                                    if _psf_nn_tree is not None and len(xy_iso):
                                        _d_psf, _i_psf = _psf_nn_tree.query(
                                            np.asarray(xy_iso, dtype=float), k=1, workers=1
                                        )
                                        _psf_flux = _all_flux_sub[np.asarray(_i_psf, dtype=int)]
                                    _cleaned_img = np.array(
                                        img_sub, dtype=np.float32, copy=True
                                    )
                                    _allstar_apply_model_inplace(
                                        _cleaned_img,
                                        _xy_sub[:, 0],
                                        _xy_sub[:, 1],
                                        _all_flux_sub,
                                        _rough_eval,
                                        _render_sz,
                                        subtract=True,
                                    )
                                    _allstar_apply_model_inplace(
                                        _cleaned_img,
                                        xy_iso[:, 0],
                                        xy_iso[:, 1],
                                        _psf_flux,
                                        _rough_eval,
                                        _render_sz,
                                        subtract=False,
                                    )

                                    # Neighbour-cleaned image:
                                    # img_sub - all_model + psf_only_model
                                    # ≡ img_sub - neighbour_model
                                    _nddata_clean = _NDData(data=_cleaned_img)
                                    _stars_clean = extract_stars(
                                        _nddata_clean, star_table, size=epsf_size_frame
                                    )
                                    if len(_stars_clean) < 3:
                                        self._log(
                                            f"[EPSF] substar {_si+1}/{n_substar_iters}"
                                            f" | too few clean stars ({len(_stars_clean)})"
                                            " → stop"
                                        )
                                        break
                                    _rough_epsf, _ = builder(_stars_clean)
                                    if psf_type_built == 'moffat_hybrid':
                                        # The rebuild produces a plain ePSF; keep
                                        # the model the chosen mode promised by
                                        # re-splitting it against the analytic.
                                        _rough_epsf = build_moffat_hybrid_psf(
                                            _rough_epsf, epsf.analytic, oversampling)
                                    elif psf_type_built == 'moffat':
                                        # Pure-analytic mode has nothing to refine
                                        # from a re-stack; keep the fitted Moffat.
                                        _rough_epsf = epsf
                                    self._log(
                                        f"[EPSF] substar {_si+1}/{n_substar_iters}"
                                        f" | n_psf={len(_stars_clean)}"
                                        f" | neighbours from {len(_xy_sub)} sources"
                                    )
                                    _stars_clean = None
                                    _nddata_clean = None
                                    _cleaned_img = None

                                epsf = _rough_epsf

                            except Exception as _se:
                                self._log(
                                    f"[EPSF] substar rebuild error: {_se}"
                                    " | using initial EPSF"
                                )

                        _stars_clean = None
                        _nddata_clean = None
                        _cleaned_img = None

                        _t["epsf_done"] = time.time()
                        if psf_build_mode_cfg == 'moffat':
                            # Render Moffat to native-scale 2D array for PSF tab display
                            _disp_half = max(25, int(fwhm_safe * 4))
                            _yy_d, _xx_d = np.mgrid[
                                -_disp_half:_disp_half + 1,
                                -_disp_half:_disp_half + 1,
                            ].astype(float)
                            _moffat_disp = epsf(_xx_d, _yy_d)
                            epsf_emit_arr = np.asarray(_moffat_disp, dtype=np.float32)
                        else:
                            epsf_emit_arr = epsf.data.copy()

                        if use_shared_filter_epsf:
                            epsf_path = output_dir / f"epsf_model_{this_filter}.fits"
                        else:
                            epsf_path = output_dir / f"epsf_model_{this_filter}_{Path(fname).stem}.fits"
                        hdr = fits.Header()
                        hdr["FILTER"] = this_filter
                        hdr["OVERSAMPL"] = oversampling
                        hdr["NSTARS"] = len(stars_extracted)
                        hdr["NDETECT"] = int(n_epsf_detected)
                        hdr["NCAND"] = int(n_epsf_candidates)
                        hdr["NCPRE"] = int(n_epsf_candidates_pre_morph)
                        hdr["NCPOST"] = int(n_epsf_candidates_post_morph)
                        hdr["NISOL"] = int(n_iso)
                        hdr["NSELECT"] = int(n_epsf_selected)
                        hdr["NMRPH"] = int(n_epsf_morphology_relaxed_selected)
                        hdr["PLANTRG"] = int(epsf_plan_target)
                        hdr["GRID"] = int(epsf_grid_size)
                        hdr["CTMAWARE"] = bool(epsf_contamination_filter)
                        hdr["NLOWCONT"] = int(n_epsf_low_contamination)
                        hdr["NCOREREJ"] = int(n_epsf_core_rejected)
                        hdr["NFALLBK"] = int(n_epsf_fallback_selected)
                        if np.isfinite(epsf_selected_median_contamination):
                            hdr["MEDCONT"] = float(epsf_selected_median_contamination)
                        hdr["EPSFSIZE"] = int(epsf_size_frame)
                        fits.writeto(str(epsf_path), epsf.data.astype(np.float32), hdr, overwrite=True)
                        self._log(
                            f"[EPSF] filter={this_filter} | "
                            f"n_stars={len(stars_extracted)} | oversampling={oversampling} | "
                            f"epsf_size={epsf_size_frame} | fit_shape={fit_shape_frame}"
                        )
                        with epsf_cache_lock:
                            _has_cached = epsf_cache_key in epsf_cache
                            _enough_stars = n_iso >= min_epsf_stars
                            if use_shared_filter_epsf and not _enough_stars and _has_cached:
                                # Too few isolated stars — reuse existing shared ePSF
                                self._log(
                                    f"[EPSF] {fname}: only {n_iso} isolated stars "
                                    f"(min={min_epsf_stars}) → reusing cached ePSF for filter={this_filter}"
                                )
                                self.log.emit(
                                    f"⚠ EPSF [{fname}]: {n_iso} isolated stars < {min_epsf_stars} → using shared ePSF"
                                )
                            else:
                                if not _has_cached:
                                    epsf_cache[epsf_cache_key] = epsf
                                    if use_shared_filter_epsf and not _enough_stars:
                                        self._log(
                                            f"[WARN][EPSF] {fname}: only {n_iso} isolated stars "
                                            f"(min={min_epsf_stars}), no cached ePSF yet → using this frame's ePSF"
                                        )
                            epsf_model = epsf_cache[epsf_cache_key]

                    try:
                        _policy_evaluator = _make_psf_evaluator(
                            epsf_model, psf_type_built, oversampling,
                            psf_interp_order
                        )
                        _native_psf_policy = _sample_native_psf(
                            _policy_evaluator, epsf_size_frame
                        )
                    except Exception as _fit_policy_error:
                        self._log(
                            f"[WARN][FIT] PSF sampling failed: {_fit_policy_error}; "
                            "using manual footprint fallback"
                        )
                        _native_psf_policy = np.zeros((3, 3), dtype=float)
                    fit_window_plan = plan_psf_fit_window(
                        _native_psf_policy,
                        fwhm_safe,
                        mode=fit_window_mode,
                        manual_fwhm_mult=fit_shape_fwhm_mult,
                        target_energy_fraction=fit_encircled_energy,
                        minimum_fwhm_mult=2.0,
                        maximum_size_px=min(31, max(9, epsf_size_frame - 4)),
                    )
                    fit_shape_frame = int(fit_window_plan.shape_px)
                    render_shape_frame = _odd_int(
                        max(float(epsf_size_frame) * 2.0, float(fit_shape_frame)),
                        min_value=11,
                        max_value=201,
                    )
                    psf_nea_frame = float(fit_window_plan.noise_equivalent_area_px)
                    self._log(
                        "  [FIT] window policy | "
                        f"mode={fit_window_plan.mode} shape={fit_shape_frame}px "
                        f"energy={fit_window_plan.energy_fraction:.3f}/"
                        f"{fit_window_plan.target_energy_fraction:.3f} "
                        f"NEA={psf_nea_frame:.1f}px reason={fit_window_plan.reason}"
                    )

                    from astropy.table import Table as AstropyTable
                    xy_det = det_df[["x", "y"]].to_numpy(float)
                    finite_xy = np.isfinite(xy_det[:, 0]) & np.isfinite(xy_det[:, 1])
                    xy_det = xy_det[finite_xy]
                    det_uids = det_df["det_uid"].to_numpy(int)[finite_xy]
                    init_forced_positions = np.zeros(len(xy_det), dtype=bool)

                    # Remove edge detections that cannot support fit window.
                    edge_init = fit_shape_frame // 2 + 2
                    valid_init = (
                        (xy_det[:, 0] >= edge_init) & (xy_det[:, 0] < (w - edge_init)) &
                        (xy_det[:, 1] >= edge_init) & (xy_det[:, 1] < (h - edge_init))
                    )
                    n_init_drop = int(np.count_nonzero(~valid_init))
                    xy_det = xy_det[valid_init]
                    det_uids = det_uids[valid_init]
                    init_forced_positions = init_forced_positions[valid_init]
                    if len(xy_det) == 0:
                        return {
                            "file": fname,
                            "status": "no_valid_init",
                            "reason": f"all detections near edge for fit_shape={fit_shape_frame}",
                        }

                    # Exclude saturated sources from PSF fitting.
                    # EPSF cannot model saturated profiles; including them degrades
                    # the fit for nearby unsaturated sources as well.
                    xi_init = xy_det[:, 0].astype(int).clip(0, w - 1)
                    yi_init = xy_det[:, 1].astype(int).clip(0, h - 1)
                    not_sat_init = img[yi_init, xi_init] < sat_adu
                    n_sat_drop = int(np.count_nonzero(~not_sat_init))
                    if n_sat_drop > 0:
                        self._log(
                            f"  [init_params] excluded {n_sat_drop} saturated sources "
                            f"(peak ≥ {sat_adu:.0f} ADU) from PSF fitting"
                        )
                    xy_det = xy_det[not_sat_init]
                    det_uids = det_uids[not_sat_init]
                    init_forced_positions = init_forced_positions[not_sat_init]
                    if len(xy_det) == 0:
                        return {
                            "file": fname,
                            "status": "no_valid_init",
                            "reason": f"all detections saturated (sat_adu={sat_adu:.0f})",
                        }

                    ap_tsv = step7_forced_phot_dir(self.result_dir) / f"photometry_{fname}.tsv"
                    flux_init_map = {}
                    df_ap = pd.DataFrame()
                    if ap_tsv.exists():
                        try:
                            df_ap = pd.read_csv(ap_tsv, sep="\t")
                            _ap_cols = set(df_ap.columns)
                            if "det_uid" in _ap_cols:
                                _uid = pd.to_numeric(df_ap["det_uid"], errors="coerce")
                                # Use ADU flux (matches PSF fitting image units).
                                # flux_net_adu is the sky-subtracted aperture flux in ADU.
                                # flux_e is in electrons = flux_net_adu × GAIN (10× smaller
                                # for gain=0.1); using electrons as flux_0 shifts the LM
                                # optimizer 10× from the true minimum and causes flux
                                # redistribution errors in crowded group fits.
                                if "flux_net_adu" in _ap_cols:
                                    _flx = pd.to_numeric(df_ap["flux_net_adu"], errors="coerce")
                                elif "flux_e" in _ap_cols:
                                    _flx = pd.to_numeric(df_ap["flux_e"], errors="coerce") / max(GAIN, 1e-6)
                                else:
                                    _flx = None
                                if _flx is not None:
                                    _ok = _uid.notna() & _flx.notna() & (_flx > 0)
                                    if _ok.any():
                                        for _u, _v in zip(
                                            _uid.loc[_ok].to_numpy(dtype=np.int64, copy=False),
                                            _flx.loc[_ok].to_numpy(dtype=float, copy=False),
                                        ):
                                            flux_init_map[int(_u)] = float(_v)
                            if psf_fit_engine_cfg == "apex_iterative" and {"x_fit", "y_fit"} <= _ap_cols:
                                if "flux_net_adu" in _ap_cols:
                                    _seed_flux = pd.to_numeric(df_ap["flux_net_adu"], errors="coerce")
                                elif "flux_e" in _ap_cols:
                                    _seed_flux = pd.to_numeric(df_ap["flux_e"], errors="coerce") / max(GAIN, 1e-6)
                                else:
                                    _seed_flux = pd.Series(np.nan, index=df_ap.index)

                                _seed_x = pd.to_numeric(df_ap["x_fit"], errors="coerce")
                                _seed_y = pd.to_numeric(df_ap["y_fit"], errors="coerce")
                                _seed_uid = (
                                    pd.to_numeric(df_ap["det_uid"], errors="coerce")
                                    if "det_uid" in _ap_cols
                                    else pd.Series(np.nan, index=df_ap.index)
                                )

                                def _bool_series(col: str, default: bool) -> pd.Series:
                                    if col not in _ap_cols:
                                        return pd.Series(default, index=df_ap.index)
                                    raw = df_ap[col]
                                    if raw.dtype == bool:
                                        return raw.fillna(default)
                                    text = raw.astype(str).str.strip().str.lower()
                                    true_vals = {"1", "true", "t", "yes", "y"}
                                    false_vals = {"0", "false", "f", "no", "n", ""}
                                    out = text.map(
                                        lambda v: True if v in true_vals else (False if v in false_vals else default)
                                    )
                                    return out.astype(bool)

                                _forced_like = (
                                    _bool_series("forced_flag", False)
                                    | (~_bool_series("detected_flag", True))
                                    | (_seed_uid.fillna(-1) < 0)
                                )
                                _phot_ok = (
                                    ~_bool_series("off_frame_flag", False)
                                    & ~_bool_series("is_saturated", False)
                                    & ~_bool_series("is_nonlinear", False)
                                )

                                _edge = fit_shape_frame // 2 + 2
                                _seed_ok = (
                                    _forced_like
                                    & _phot_ok
                                    & _seed_x.notna()
                                    & _seed_y.notna()
                                    & _seed_flux.notna()
                                    & (_seed_x >= _edge)
                                    & (_seed_x < (w - _edge))
                                    & (_seed_y >= _edge)
                                    & (_seed_y < (h - _edge))
                                )
                                if _seed_ok.any():
                                    _sx = _seed_x.loc[_seed_ok].to_numpy(dtype=float, copy=False)
                                    _sy = _seed_y.loc[_seed_ok].to_numpy(dtype=float, copy=False)
                                    _sf = _seed_flux.loc[_seed_ok].to_numpy(dtype=float, copy=False)
                                    _smid = (
                                        pd.to_numeric(df_ap.loc[_seed_ok, "master_id"], errors="coerce").to_numpy(dtype=float, copy=False)
                                        if "master_id" in _ap_cols
                                        else np.full(len(_sx), np.nan, dtype=float)
                                    )

                                    _xi = np.rint(_sx).astype(int).clip(0, w - 1)
                                    _yi = np.rint(_sy).astype(int).clip(0, h - 1)
                                    _not_sat = img[_yi, _xi] < sat_adu
                                    _sx, _sy, _sf, _smid = _sx[_not_sat], _sy[_not_sat], _sf[_not_sat], _smid[_not_sat]

                                    _forced_match_radius_px = max(
                                        1.0,
                                        forced_match_radius_fwhm * float(fwhm_safe),
                                    )
                                    _merge = merge_forced_catalog_seeds(
                                        xy_det,
                                        det_uids,
                                        init_forced_positions,
                                        np.column_stack([_sx, _sy]),
                                        _sf,
                                        _smid,
                                        match_radius_px=_forced_match_radius_px,
                                    )
                                    xy_det = _merge.xy
                                    det_uids = _merge.det_uids
                                    init_forced_positions = _merge.forced_mask
                                    flux_init_map.update(_merge.flux_by_uid)
                                    self._log(
                                        "  [INIT] Step7 forced catalog | "
                                        f"matched={_merge.n_matched} added={_merge.n_added} "
                                        f"radius={_forced_match_radius_px:.2f}px "
                                        f"({forced_match_radius_fwhm:.2f}xFWHM)"
                                    )
                        except Exception as exc:
                            self._log(f"  [WARN] Step7 flux/forced seed load failed: {exc}")

                    if core_cut.enabled and len(xy_det):
                        _keep_core_init = _core_keep(xy_det)
                        n_core_excluded_init = int(np.sum(~_keep_core_init))
                        if n_core_excluded_init > 0:
                            xy_det = xy_det[_keep_core_init]
                            det_uids = det_uids[_keep_core_init]
                            init_forced_positions = init_forced_positions[_keep_core_init]
                            self._log(
                                f"  [CORE] initial PSF seeds excluded: {n_core_excluded_init} "
                                f"(r<{core_cut.radius_px:.1f}px)"
                            )
                        if len(xy_det) == 0:
                            return {
                                "file": fname,
                                "status": "no_valid_init",
                                "reason": f"all detections inside PSF core cut r<{core_cut.radius_px:.1f}px",
                            }

                    default_flux = max(1.0, float(bkg_std) * 10.0)
                    init_flux_list = []
                    init_flux_from_aperture = []
                    for _seed_index, (_uid, (x0, y0)) in enumerate(zip(det_uids, xy_det)):
                        v = flux_init_map.get(int(_uid), np.nan)
                        if np.isfinite(v) and (
                            float(v) > 0 or bool(init_forced_positions[_seed_index])
                        ):
                            init_flux_list.append(float(v))
                            init_flux_from_aperture.append(True)
                            continue
                        xi0 = int(np.clip(round(float(x0)), 0, w - 1))
                        yi0 = int(np.clip(round(float(y0)), 0, h - 1))
                        pv = _safe_float(img_sub[yi0, xi0], np.nan)
                        if not np.isfinite(pv):
                            pv = default_flux
                        init_flux_list.append(max(default_flux, float(pv)))
                        init_flux_from_aperture.append(False)
                    init_flux = np.asarray(init_flux_list, dtype=float)
                    init_flux_from_aperture = np.asarray(init_flux_from_aperture, dtype=bool)
                    if fit_init_max_sources > 0 and len(xy_det) > fit_init_max_sources:
                        _ord_fit = np.argsort(np.where(np.isfinite(init_flux), init_flux, -np.inf))[::-1][:fit_init_max_sources]
                        xy_det = xy_det[_ord_fit]
                        det_uids = det_uids[_ord_fit]
                        init_flux = init_flux[_ord_fit]
                        init_flux_from_aperture = init_flux_from_aperture[_ord_fit]
                        init_forced_positions = init_forced_positions[_ord_fit]
                        self._log(
                            f"  [INIT] capped initial fit sources: kept={len(xy_det)} "
                            f"(psf_fit_init_max_sources={fit_init_max_sources})"
                        )
                    init_params = AstropyTable({"x_0": xy_det[:, 0], "y_0": xy_det[:, 1], "flux_0": init_flux})

                    # ── IterativePSFPhotometry  (Stetson 1987 / DAOPHOT style) ──────────
                    # localbkg_estimator=None: background already removed by Background2D.
                    # SourceGrouper(2.5×FWHM): Stetson's critical separation — sources
                    #   within this radius are fitted SIMULTANEOUSLY, correctly accounting
                    #   for mutual flux contamination (crowded cluster requirement).
                    # mode='all': every iteration refits ALL sources on the original data
                    #   (not just the residual), allowing later-found faint stars to improve
                    #   the fit of already-found bright neighbors.
                    # Note: photutils 2.3.0 introduced a 'flat model' that eliminates
                    #   the compound-model recursion crash seen in 2.2.0 for large groups.
                    #   If a RecursionError occurs, we fall back to no grouper.
                    # ─────────────────────────────────────────────────────────────────────

                    # Error image for photon-noise-correct flux_err
                    if use_error_image:
                        try:
                            from photutils.utils import calc_total_error
                            error_img = calc_total_error(img_sub, bkg_rms_scalar, GAIN)
                        except Exception:
                            error_img = None
                    else:
                        error_img = None

                    # Re-detection finder (used internally by IterativePSFPhotometry)
                    # Per-filter sigma overrides default redetect_sigma when specified.
                    _sigma_key = f"psf_redetect_sigma_{this_filter}"
                    _sigma_override = _to_float(getattr(P, _sigma_key, float("nan")), float("nan"))
                    if np.isfinite(_sigma_override) and _sigma_override > 0:
                        redetect_sigma_eff = _sigma_override
                    else:
                        redetect_sigma_eff = float(redetect_sigma)

                    dao_redetect_finder = DAOStarFinder(
                        fwhm=fwhm_safe,
                        threshold=redetect_sigma_eff * bkg_std,
                        peakmax=sat_adu,
                        sharplo=redetect_sharp_lo,
                        sharphi=redetect_sharp_hi,
                        roundlo=-redetect_round_abs_max,
                        roundhi=redetect_round_abs_max,
                    )
                    if core_cut.enabled:
                        def redetect_finder(data):
                            nonlocal n_core_excluded_redetect
                            tbl = dao_redetect_finder(data)
                            if tbl is None or len(tbl) == 0:
                                return tbl
                            try:
                                xy_new = np.column_stack([
                                    np.asarray(tbl["xcentroid"], dtype=float),
                                    np.asarray(tbl["ycentroid"], dtype=float),
                                ])
                                keep_core = _core_keep(xy_new)
                                n_drop_core = int(np.sum(~keep_core))
                                if n_drop_core > 0:
                                    n_core_excluded_redetect += n_drop_core
                                return tbl[keep_core]
                            except Exception:
                                return tbl
                    else:
                        redetect_finder = dao_redetect_finder
                    if redetect_sigma_eff != redetect_sigma:
                        self._log(
                            f"  [REDETECT] filter={this_filter} sigma override: {redetect_sigma:.2f} -> {redetect_sigma_eff:.2f}"
                        )

                    ap_rad = max(int(round(fwhm_safe * 2.0)), fit_shape_frame // 2 + 1)

                    def _build_iterative_phot(with_grouper: bool, n_seed: int):
                        from photutils.psf import IterativePSFPhotometry
                        import inspect as _ins
                        psf_m = _clone_psf_model(epsf_model)
                        kw: dict = dict(
                            psf_model=psf_m,
                            fit_shape=fit_shape_frame,
                            finder=redetect_finder,
                            aperture_radius=ap_rad,
                            localbkg_estimator=None,
                        )
                        sig = _ins.signature(IterativePSFPhotometry).parameters
                        if "maxiters" in sig:
                            kw["maxiters"] = max_iter
                        if "mode" in sig:
                            # mode='new' (default): iter1 fits all, iter2+ only new sources — fast
                            # mode='all': every iteration refits ALL sources — accurate but O(n×iter)
                            #   → can be slow for large fields; a performance warning is logged
                            if fit_mode_cfg == "all" and n_seed > 800:
                                self._log(
                                    f"  [PSF] fit_mode='all' | {n_seed} sources "
                                    "— expect significantly slower fitting"
                                )
                            kw["mode"] = fit_mode_cfg
                        if (
                            with_grouper
                            and _has_grouper
                            and grouper_max_size > 1
                            and "grouper" in sig
                        ):
                            _grouper_kw: dict = {"min_separation": 2.5 * fwhm_safe}
                            _sg_sig = _ins.signature(SourceGrouper).parameters
                            if "max_group_size" not in _sg_sig:
                                self._log(
                                    "  [PSF] SourceGrouper disabled: installed photutils "
                                    "cannot enforce the 3-source CPU limit"
                                )
                            else:
                                _grouper_kw["max_group_size"] = grouper_max_size
                                kw["grouper"] = SourceGrouper(**_grouper_kw)
                        self._log(
                            f"  [PSF] IterativePSFPhotometry | mode={kw.get('mode', 'N/A')} "
                            f"maxiters={kw.get('maxiters', 'N/A')} "
                            f"grouper={'on' if 'grouper' in kw else 'off'} "
                            f"n_seed={n_seed}"
                        )
                        return IterativePSFPhotometry(**kw)

                    def _results_to_init_params(results_tbl, photometry_obj=None):
                        if results_tbl is None or len(results_tbl) == 0:
                            return None
                        if photometry_obj is not None and hasattr(photometry_obj, "results_to_init_params"):
                            try:
                                tbl = photometry_obj.results_to_init_params()
                                if tbl is not None and len(tbl) > 0:
                                    return tbl
                            except Exception as _ri:
                                self._log(f"  [PSF] results_to_init_params fallback: {_ri}")
                        try:
                            cols = list(results_tbl.colnames)
                            x_col = "x_fit" if "x_fit" in cols else ("x_0" if "x_0" in cols else None)
                            y_col = "y_fit" if "y_fit" in cols else ("y_0" if "y_0" in cols else None)
                            f_col = next((c for c in ("flux_fit", "flux", "flux_0") if c in cols), None)
                            if x_col is None or y_col is None or f_col is None:
                                return None
                            x_arr = np.asarray(results_tbl[x_col], dtype=float)
                            y_arr = np.asarray(results_tbl[y_col], dtype=float)
                            f_arr = np.asarray(results_tbl[f_col], dtype=float)
                            keep = (
                                np.isfinite(x_arr) &
                                np.isfinite(y_arr) &
                                np.isfinite(f_arr) &
                                (f_arr > 0)
                            )
                            if not np.any(keep):
                                return None
                            return AstropyTable({
                                "x_0": x_arr[keep],
                                "y_0": y_arr[keep],
                                "flux_0": f_arr[keep],
                            })
                        except Exception:
                            return None

                    def _run_iterative_fit(seed_params, stage_label: str):
                        fit_reason = None
                        fit_photometry = None
                        fit_results = None
                        for _attempt, _use_grouper in enumerate(attempt_plan):
                            if self._stop_requested:
                                return None, None, "stopped"
                            try:
                                fit_photometry = _build_iterative_phot(
                                    with_grouper=_use_grouper,
                                    n_seed=len(seed_params),
                                )
                                if _attempt == 1:
                                    self._log(
                                        f"  [PSF] {stage_label} retry without SourceGrouper (fallback)"
                                    )
                                call_kw = {"init_params": seed_params}
                                if error_img is not None:
                                    call_kw["error"] = error_img
                                fit_results = fit_photometry(img_sub, **call_kw)
                                fit_reason = None
                                break
                            except RecursionError as _re:
                                self._log(
                                    f"  [PSF] {stage_label} RecursionError with grouper "
                                    f"(photutils<2.3 compound-model bug): {_re}. Retrying without grouper."
                                )
                                fit_reason = str(_re)
                                if _attempt + 1 < len(attempt_plan):
                                    continue
                            except Exception as _fe:
                                fit_reason = str(_fe)
                                self._log(
                                    f"  [PSF] {stage_label} fit failed (attempt {_attempt+1}): {fit_reason}"
                                )
                                if _attempt + 1 < len(attempt_plan):
                                    continue
                                break
                        return fit_photometry, fit_results, fit_reason

                    def _render_model_from_results(results_tbl, photometry_obj=None):
                        if results_tbl is None or len(results_tbl) == 0:
                            return None
                        try:
                            cols = list(results_tbl.colnames)
                            x_col = "x_fit" if "x_fit" in cols else ("x_0" if "x_0" in cols else None)
                            y_col = "y_fit" if "y_fit" in cols else ("y_0" if "y_0" in cols else None)
                            f_col = next((c for c in ("flux_fit", "flux", "flux_0") if c in cols), None)
                            if x_col is None or y_col is None or f_col is None:
                                return None
                            x_arr = np.asarray(results_tbl[x_col], dtype=float)
                            y_arr = np.asarray(results_tbl[y_col], dtype=float)
                            f_arr = np.asarray(results_tbl[f_col], dtype=float)
                            keep = (
                                np.isfinite(x_arr) &
                                np.isfinite(y_arr) &
                                np.isfinite(f_arr) &
                                (f_arr > 0)
                            )
                            if not np.any(keep):
                                return None
                            from photutils.datasets import make_model_image as _make_model_image
                            pt = AstropyTable()
                            pt["x_0"] = np.asarray(x_arr[keep], dtype=float)
                            pt["y_0"] = np.asarray(y_arr[keep], dtype=float)
                            pt["flux"] = np.asarray(f_arr[keep], dtype=float)
                            out = _make_model_image(
                                img_sub.shape,
                                _clone_psf_model(epsf_model),
                                pt,
                                model_shape=(int(render_shape_frame), int(render_shape_frame)),
                                x_name="x_0",
                                y_name="y_0",
                            )
                            return np.asarray(out, dtype=np.float32)
                        except Exception as _re:
                            self._log(f"  [DIAG] wide model render failed: {_re}")
                            if photometry_obj is not None:
                                try:
                                    out = photometry_obj.make_model_image(
                                        img_sub.shape,
                                        psf_shape=(int(render_shape_frame), int(render_shape_frame)),
                                    )
                                    return np.asarray(out, dtype=np.float32)
                                except Exception as _pe:
                                    self._log(f"  [DIAG] make_model_image fallback failed: {_pe}")
                        return None

                    fit_fail_reason = None
                    phot_result = None
                    photometry = None
                    model_img = None

                    _t["fit1"] = time.time()
                    self.progress.emit(completed[0], total, f"FIT | {fname}")

                    _psf_evaluator = None  # set in APEX branch; used for final residual rendering
                    _engine_snapshots: list[dict] = []
                    _engine_stop_reason = ""

                    if psf_fit_engine_cfg == 'apex_iterative':
                        self.worker_status.emit(wid, fname, "APEX iterative fit", 70)
                        _psf_evaluator = _make_psf_evaluator(
                            epsf_model, psf_type_built, oversampling, psf_interp_order)

                        # Outer loop: ALLSTAR fit → residual → re-detect → add → repeat
                        # max_iter controls number of find+fit cycles (DAOPHOT style)
                        _INNER_ITERS = fitter_max_iter
                        _max_grp_size, _group_budget = local_group_policy(
                            len(init_params),
                            enabled=use_grouper,
                            requested_max_size=grouper_max_size,
                            hard_max_size=25,
                            max_fraction=0.10,
                            absolute_cap=200,
                        )
                        # Respect the shared grouper switch for the APEX engine too.
                        # With the core excluded, single-star neighbour subtraction is
                        # usually the practical fast path.
                        if use_grouper and _max_grp_size > 1:
                            _group_radius = fwhm_safe * grouper_radius_fwhm
                        else:
                            _group_radius = 0.0
                            _max_grp_size = 1
                            _group_budget = 0
                        self._log(
                            f"  [APEX] local grouping | max_size={_max_grp_size} "
                            f"radius={grouper_radius_fwhm:.2f}xFWHM "
                            f"budget={_group_budget}/{len(init_params)}"
                        )

                        _cur_xy  = np.column_stack([
                            np.asarray(init_params["x_0"], dtype=float),
                            np.asarray(init_params["y_0"], dtype=float),
                        ])
                        _cur_flux  = np.asarray(init_params["flux_0"], dtype=float)
                        _needs_psf_flux_seed = ~np.asarray(
                            init_flux_from_aperture, dtype=bool
                        )
                        if np.any(_needs_psf_flux_seed):
                            _estimated_initial_flux = estimate_psf_flux_seeds(
                                img_sub,
                                _cur_xy[_needs_psf_flux_seed],
                                _psf_evaluator,
                                fit_shape=fit_shape_frame,
                                fallback=_cur_flux[_needs_psf_flux_seed],
                            )
                            _replace_flux_seed = (
                                np.isfinite(_estimated_initial_flux)
                                & (_estimated_initial_flux > 0)
                            )
                            _fallback_indices = np.flatnonzero(_needs_psf_flux_seed)
                            _cur_flux[_fallback_indices[_replace_flux_seed]] = (
                                _estimated_initial_flux[_replace_flux_seed]
                            )
                            self._log(
                                "  [APEX] ePSF flux seeds | "
                                f"estimated={int(np.sum(_replace_flux_seed))}/"
                                f"{int(np.sum(_needs_psf_flux_seed))} "
                                "(Step7 aperture seed unavailable)"
                            )
                        _cur_idet  = np.ones(len(_cur_xy), dtype=int)  # iter_detected
                        _cur_forced = np.asarray(init_forced_positions, dtype=bool).copy()
                        _cur_anchor = _cur_xy.copy()
                        _cur_position_flags = np.zeros(len(_cur_xy), dtype=np.int32)
                        _cur_fit_valid = np.zeros(len(_cur_xy), dtype=bool)
                        _n_initial = len(_cur_xy)

                        _dedup_r = (float(duplicate_radius_px_cfg)
                                    if np.isfinite(duplicate_radius_px_cfg) and duplicate_radius_px_cfg > 0
                                    else float(max(0.0, duplicate_radius_mult * fwhm_safe)))
                        _dedup_r = max(_dedup_r, 1.0)

                        _outer_fit_result = None
                        _stop_after_refit = False
                        _pending_stop_reason = ""

                        def _save_actual_snapshot(
                            *,
                            sequence: int,
                            phase: str,
                            fit_result,
                            residual_image: np.ndarray,
                            model_image: np.ndarray,
                            candidate_xy: np.ndarray,
                            candidate_counts: tuple[int, int, int],
                            n_pruned: int,
                            stop_reason: str,
                            elapsed_s: float,
                        ) -> None:
                            fit_xy_snapshot = np.column_stack([
                                np.asarray(fit_result["x_fit"], dtype=float),
                                np.asarray(fit_result["y_fit"], dtype=float),
                            ])
                            new_mask_snapshot = (
                                np.zeros(len(_cur_idet), dtype=bool)
                                if phase == "final_flux"
                                else _cur_idet == sequence
                            )
                            model_xy_snapshot = (
                                fit_xy_snapshot[new_mask_snapshot]
                                if len(new_mask_snapshot) == len(fit_xy_snapshot)
                                else np.zeros((0, 2), dtype=float)
                            )
                            applied_xy_snapshot = (
                                fit_xy_snapshot[~new_mask_snapshot]
                                if len(new_mask_snapshot) == len(fit_xy_snapshot)
                                else fit_xy_snapshot
                            )
                            residual_path = output_dir / f"residual_iter{sequence}_{fname}"
                            model_path = output_dir / f"model_iter{sequence}_{fname}"
                            starsub_path = output_dir / f"starsub_iter{sequence}_{fname}"
                            fitxy_path = output_dir / f"fitxy_iter{sequence}_{fname}.npy"
                            modelxy_path = output_dir / f"modelxy_iter{sequence}_{fname}.npy"
                            detxy_path = output_dir / f"detxy_iter{sequence}_{fname}.npy"
                            candidatexy_path = output_dir / f"candidatexy_iter{sequence}_{fname}.npy"
                            appliedxy_path = output_dir / f"appliedxy_iter{sequence}_{fname}.npy"
                            boxxy_path = output_dir / f"boxxy_iter{sequence}_{fname}.npy"

                            header_snapshot = fits.Header()
                            header_snapshot["FILTER"] = this_filter
                            header_snapshot["BKGMED"] = float(bkg_med)
                            header_snapshot["ITER"] = int(sequence)
                            header_snapshot["PHASE"] = str(phase)[:16]
                            fits.CompImageHDU(
                                np.asarray(residual_image, dtype=np.float32),
                                header_snapshot,
                                compression_type="RICE_1",
                                quantize_level=16.0,
                            ).writeto(str(residual_path), overwrite=True)
                            fits.CompImageHDU(
                                np.asarray(model_image, dtype=np.float32),
                                header_snapshot,
                                compression_type="RICE_1",
                                quantize_level=16.0,
                            ).writeto(str(model_path), overwrite=True)
                            starsub_name = None
                            if save_all_iter_residuals:
                                fits.CompImageHDU(
                                    np.asarray(img - model_image, dtype=np.float32),
                                    header_snapshot,
                                    compression_type="RICE_1",
                                    quantize_level=16.0,
                                ).writeto(str(starsub_path), overwrite=True)
                                starsub_name = starsub_path.name

                            np.save(str(fitxy_path), np.asarray(fit_xy_snapshot, dtype=np.float32))
                            np.save(str(modelxy_path), np.asarray(model_xy_snapshot, dtype=np.float32))
                            np.save(str(detxy_path), np.asarray(model_xy_snapshot, dtype=np.float32))
                            np.save(str(candidatexy_path), np.asarray(candidate_xy, dtype=np.float32))
                            np.save(str(appliedxy_path), np.asarray(applied_xy_snapshot, dtype=np.float32))
                            np.save(str(boxxy_path), np.asarray(model_xy_snapshot, dtype=np.float32))

                            raw_count, unique_count, accepted_count = candidate_counts
                            qfit_values = np.asarray(fit_result["qfit"], dtype=float)
                            redchi_values = np.asarray(fit_result["reduced_chi2"], dtype=float)
                            summary = IterationSnapshot(
                                iteration=sequence,
                                n_fit=len(fit_result),
                                n_candidates_raw=raw_count,
                                n_candidates_unique=unique_count,
                                n_candidates_accepted=accepted_count,
                                residual_std=float(_fast_res_std(residual_image)),
                                median_qfit=(
                                    float(np.nanmedian(qfit_values))
                                    if np.any(np.isfinite(qfit_values)) else np.nan
                                ),
                                median_reduced_chi2=(
                                    float(np.nanmedian(redchi_values))
                                    if np.any(np.isfinite(redchi_values)) else np.nan
                                ),
                                elapsed_s=float(elapsed_s),
                                stop_reason=str(stop_reason),
                            ).to_dict()
                            summary.update({
                                "iter": int(sequence),
                                "phase": str(phase),
                                "fit_shape_px": int(fit_shape_frame),
                                "epsf_size_px": int(epsf_size_frame),
                                "n_new_raw": int(np.sum(new_mask_snapshot)) if sequence > 1 else 0,
                                "n_new_kept": int(np.sum(new_mask_snapshot)) if sequence > 1 else 0,
                                "n_pruned": int(n_pruned),
                                "n_applied_prev": int(len(applied_xy_snapshot)),
                                "residual_path": residual_path.name,
                                "model_path": model_path.name,
                                "starsub_path": starsub_name,
                                "fitxy_path": fitxy_path.name,
                                "modelxy_path": modelxy_path.name,
                                "detxy_path": detxy_path.name,
                                "candidatexy_path": candidatexy_path.name,
                                "appliedxy_path": appliedxy_path.name,
                                "boxxy_path": boxxy_path.name,
                            })
                            _engine_snapshots.append(summary)

                        for _outer in range(max_iter):
                            if self._stop_requested:
                                break
                            _outer_started = time.perf_counter()
                            self._log(f"  [APEX] outer {_outer+1}/{max_iter} | n={len(_cur_xy)}")
                            _outer_active_mask = None
                            if _outer > 0:
                                _new_generation_before_fit = _cur_idet == (_outer + 1)
                                _retry_previous_fit = (
                                    _cur_position_flags
                                    & int(PSFFitFlag.NONCONVERGENCE)
                                ) != 0
                                _outer_active_mask = (
                                    _new_generation_before_fit | _retry_previous_fit
                                )
                                if np.any(_outer_active_mask):
                                    _local_refit_radius = max(
                                        min(
                                            2.5,
                                            max(1.5, float(grouper_radius_fwhm)),
                                        ) * float(fwhm_safe),
                                        float(fit_shape_frame // 2),
                                        float(_group_radius) if _max_grp_size > 1 else 0.0,
                                    )
                                    _local_tree = cKDTree(_cur_xy)
                                    _local_neighbors = _local_tree.query_ball_point(
                                        _cur_xy[_outer_active_mask],
                                        r=_local_refit_radius,
                                        workers=1,
                                    )
                                    for _indices in _local_neighbors:
                                        _outer_active_mask[
                                            np.asarray(_indices, dtype=int)
                                        ] = True
                                    self._log(
                                        "  [APEX] local residual refit | "
                                        f"new={int(np.sum(_new_generation_before_fit))} "
                                        f"retry={int(np.sum(_retry_previous_fit))} "
                                        f"active={int(np.sum(_outer_active_mask))}/{len(_cur_xy)} "
                                        f"radius={_local_refit_radius:.1f}px"
                                    )
                            _previous_position_flags = _cur_position_flags.copy()
                            _outer_fit_result = _allstar_fit(
                                img_sub, _cur_xy, _cur_flux,
                                _psf_evaluator,
                                fit_shape=fit_shape_frame,
                                stamp_size=render_shape_frame,
                                max_iter=_INNER_ITERS,
                                flux_conv=flux_conv_threshold,
                                max_shift=float(fit_shape_frame // 2),
                                group_radius=_group_radius,
                                max_group_size=_max_grp_size,
                                max_grouped_sources=_group_budget,
                                background_rms=float(bkg_std),
                                gain=float(GAIN),
                                initial_positions=_cur_anchor,
                                initial_fit_valid=_cur_fit_valid,
                                position_bound=float(fit_shape_frame // 2),
                                position_fixed_mask=_cur_forced,
                                allow_negative_flux_mask=_cur_forced,
                                fit_active_mask=_outer_active_mask,
                                log_fn=self._log,
                                stop_fn=lambda: self._stop_requested,
                            )

                            if _outer_active_mask is not None:
                                _local_flags = np.asarray(
                                    _outer_fit_result["flags"], dtype=np.int32
                                )
                                _local_flags[~_outer_active_mask] |= (
                                    _previous_position_flags[~_outer_active_mask]
                                )
                                _outer_fit_result["flags"] = _local_flags

                            # Update positions/fluxes from fit
                            _cur_xy   = np.column_stack([
                                np.asarray(_outer_fit_result["x_fit"],   dtype=float),
                                np.asarray(_outer_fit_result["y_fit"],   dtype=float),
                            ])
                            _cur_flux = np.asarray(_outer_fit_result["flux_fit"], dtype=float)
                            _cur_position_flags = np.asarray(
                                _outer_fit_result["flags"], dtype=np.int32
                            )
                            _fit_invalid_mask = int(
                                PSFFitFlag.NONCONVERGENCE
                                | PSFFitFlag.NO_OVERLAP
                                | PSFFitFlag.NONFINITE_POSITION
                                | PSFFitFlag.NONFINITE_FLUX
                            )
                            _cur_fit_valid = (
                                np.isfinite(_cur_xy[:, 0])
                                & np.isfinite(_cur_xy[:, 1])
                                & np.isfinite(_cur_flux)
                                & ((_cur_flux > 0) | _cur_forced)
                                & ((_cur_position_flags & _fit_invalid_mask) == 0)
                            )

                            _n_pruned = 0
                            _new_generation = _cur_idet == (_outer + 1)
                            if _outer > 0 and np.any(_new_generation):
                                _fit_flags = np.asarray(_outer_fit_result["flags"], dtype=np.int64)
                                _fit_err = np.asarray(_outer_fit_result["flux_err"], dtype=float)
                                _fit_snr = np.divide(
                                    _cur_flux,
                                    _fit_err,
                                    out=np.full_like(_cur_flux, np.nan),
                                    where=np.isfinite(_fit_err) & (_fit_err > 0),
                                )
                                _fit_qfit = np.asarray(_outer_fit_result["qfit"], dtype=float)
                                _fit_npix = np.asarray(
                                    _outer_fit_result["n_pixels_fit"], dtype=float
                                )
                                _fit_qfit_expected, _fit_qfit_ratio = qfit_noise_diagnostics(
                                    _fit_qfit,
                                    _fit_npix,
                                    _fit_snr,
                                    psf_nea_frame,
                                )
                                _fit_redchi = np.asarray(
                                    _outer_fit_result["reduced_chi2"], dtype=float
                                )
                                _severe_flags = int(
                                    PSFFitFlag.NONPOSITIVE_FLUX
                                    | PSFFitFlag.NONCONVERGENCE
                                    | PSFFitFlag.NONFINITE_POSITION
                                    | PSFFitFlag.NONFINITE_FLUX
                                    | PSFFitFlag.NO_OVERLAP
                                )
                                _keep_new_fit = (
                                    np.isfinite(_cur_flux)
                                    & (_cur_flux > 0)
                                    & ((_fit_flags & _severe_flags) == 0)
                                    & (np.isfinite(_fit_snr) & (_fit_snr >= postfit_snr_min))
                                    & (
                                        (postfit_qfit_max <= 0)
                                        | (
                                            np.isfinite(_fit_qfit_ratio)
                                            & (_fit_qfit_ratio <= postfit_qfit_max)
                                        )
                                    )
                                    & (
                                        (postfit_reduced_chi2_max <= 0)
                                        | (
                                            np.isfinite(_fit_redchi)
                                            & (_fit_redchi <= postfit_reduced_chi2_max)
                                        )
                                    )
                                )
                                _keep_fit = (~_new_generation) | _keep_new_fit
                                _n_pruned = int(np.sum(~_keep_fit))
                                if _n_pruned > 0:
                                    self._log(
                                        f"  [APEX] post-fit prune | generation={_outer + 1} "
                                        f"removed={_n_pruned} snr_min={postfit_snr_min:.1f} "
                                        f"qfit/noise_max={postfit_qfit_max:.2f} "
                                        f"redchi_max={postfit_reduced_chi2_max:.1f}"
                                    )
                                    _outer_fit_result = _outer_fit_result[_keep_fit]
                                    _cur_xy = _cur_xy[_keep_fit]
                                    _cur_flux = _cur_flux[_keep_fit]
                                    _cur_idet = _cur_idet[_keep_fit]
                                    _cur_forced = _cur_forced[_keep_fit]
                                    _cur_anchor = _cur_anchor[_keep_fit]
                                    _cur_position_flags = _cur_position_flags[_keep_fit]
                                    _cur_fit_valid = _cur_fit_valid[_keep_fit]

                            # Build residual image
                            _model_temp = _allstar_build_model(
                                img_sub.shape,
                                _cur_xy[:, 0], _cur_xy[:, 1], _cur_flux,
                                _psf_evaluator, render_shape_frame,
                            )
                            _resid_temp = _float32_difference(img_sub, _model_temp)

                            _candidate_raw = 0
                            _candidate_unique = 0
                            _candidate_accepted = 0
                            _candidate_xy = np.zeros((0, 2), dtype=float)
                            _new_flux_u = np.zeros(0, dtype=float)
                            _stop_reason = ""

                            if _stop_after_refit:
                                _stop_reason = _pending_stop_reason or "candidate_fraction"
                            elif _outer > 0 and _n_pruned > 0 and not np.any(_cur_idet == (_outer + 1)):
                                _stop_reason = "postfit_pruned_all"
                            elif _outer == max_iter - 1:
                                _stop_reason = "max_residual_passes"
                            else:
                                try:
                                    _new_tbl = redetect_finder(_resid_temp)
                                except Exception as _rd_e:
                                    self._log(f"  [APEX] residual detect error: {_rd_e}")
                                    _new_tbl = None
                                    _stop_reason = "detection_error"

                                if not _stop_reason:
                                    _candidate_raw = int(len(_new_tbl)) if _new_tbl is not None else 0
                                    if _candidate_raw > 0:
                                        _new_x = np.asarray(_new_tbl["xcentroid"], dtype=float)
                                        _new_y = np.asarray(_new_tbl["ycentroid"], dtype=float)
                                        _new_pk = np.asarray(_new_tbl["peak"], dtype=float)
                                        _new_xy_all = np.column_stack([_new_x, _new_y])
                                        _tree_cur = cKDTree(_cur_xy)
                                        _d_cur, _ = _tree_cur.query(_new_xy_all, k=1, workers=1)
                                        _unique_mask = np.asarray(_d_cur, dtype=float) > _dedup_r
                                        _new_xy_u = _new_xy_all[_unique_mask]
                                        _new_peak_u = np.maximum(_new_pk[_unique_mask], 1.0)

                                        if blend_residual_ratio > 0 and len(_new_xy_u):
                                            _cx = np.rint(_new_xy_u[:, 0]).astype(int).clip(0, w - 1)
                                            _cy = np.rint(_new_xy_u[:, 1]).astype(int).clip(0, h - 1)
                                            _model_level = np.maximum(
                                                np.abs(_model_temp[_cy, _cx]),
                                                float(bkg_std),
                                            )
                                            _blend_keep = (
                                                _new_peak_u / _model_level
                                            ) >= blend_residual_ratio
                                            _new_xy_u = _new_xy_u[_blend_keep]
                                            _new_peak_u = _new_peak_u[_blend_keep]

                                        _candidate_unique = int(len(_new_xy_u))
                                        _new_flux_u = estimate_psf_flux_seeds(
                                            _resid_temp,
                                            _new_xy_u,
                                            _psf_evaluator,
                                            fit_shape=fit_shape_frame,
                                            fallback=_new_peak_u,
                                        )

                                        _cap = (
                                            new_sources_cap_per_iter
                                            if new_sources_cap_per_iter > 0 else len(_new_xy_u)
                                        )
                                        _cap_f = (
                                            max(1, int(np.floor(new_sources_cap_frac * _n_initial)))
                                            if new_sources_cap_frac > 0 and len(_new_xy_u) > 0
                                            else len(_new_xy_u)
                                        )
                                        _cap = min(_cap, _cap_f)
                                        if _cap < len(_new_xy_u):
                                            _order = np.argsort(_new_flux_u)[::-1][:_cap]
                                            _new_xy_u = _new_xy_u[_order]
                                            _new_flux_u = _new_flux_u[_order]
                                        _candidate_xy = np.asarray(_new_xy_u, dtype=float)
                                        _candidate_accepted = int(len(_candidate_xy))

                                    _decision = decide_residual_iteration(
                                        n_candidates_raw=_candidate_raw,
                                        n_candidates_unique=_candidate_unique,
                                        n_candidates_accepted=_candidate_accepted,
                                        n_current=len(_cur_xy),
                                        convergence_fraction=conv_new_frac,
                                    )
                                    if _decision.stop_now:
                                        _stop_reason = _decision.reason
                                    elif _decision.stop_after_refit:
                                        _stop_after_refit = True
                                        _pending_stop_reason = _decision.reason
                                        self._log(
                                            "  [APEX] convergence requested after fitting accepted sources "
                                            f"(pre-cap candidate_frac={_decision.candidate_fraction:.4f})"
                                        )

                            _save_actual_snapshot(
                                sequence=_outer + 1,
                                phase="residual_fit",
                                fit_result=_outer_fit_result,
                                residual_image=_resid_temp,
                                model_image=_model_temp,
                                candidate_xy=_candidate_xy,
                                candidate_counts=(
                                    _candidate_raw,
                                    _candidate_unique,
                                    _candidate_accepted,
                                ),
                                n_pruned=_n_pruned,
                                stop_reason=_stop_reason,
                                elapsed_s=time.perf_counter() - _outer_started,
                            )

                            if _stop_reason:
                                _engine_stop_reason = _stop_reason
                                self._log(f"  [APEX] stopped: {_stop_reason}")
                                del _model_temp, _resid_temp
                                break

                            self._log(
                                f"  [APEX] residual candidates | raw={_candidate_raw} "
                                f"unique={_candidate_unique} accepted={_candidate_accepted}"
                            )
                            _new_idet = np.full(_candidate_accepted, _outer + 2, dtype=int)
                            _cur_xy = np.vstack([_cur_xy, _candidate_xy])
                            _cur_flux = np.concatenate([_cur_flux, _new_flux_u])
                            _cur_idet = np.concatenate([_cur_idet, _new_idet])
                            _cur_forced = np.concatenate([
                                _cur_forced,
                                np.zeros(_candidate_accepted, dtype=bool),
                            ])
                            _cur_anchor = np.vstack([_cur_anchor, _candidate_xy])
                            _cur_position_flags = np.concatenate([
                                _cur_position_flags,
                                np.zeros(_candidate_accepted, dtype=np.int32),
                            ])
                            _cur_fit_valid = np.concatenate([
                                _cur_fit_valid,
                                np.zeros(_candidate_accepted, dtype=bool),
                            ])
                            del _model_temp, _resid_temp

                        # Final fixed-position pass stabilizes flux after source discovery.
                        _final_started = time.perf_counter()
                        phot_result = _allstar_fit(
                            img_sub,
                            _cur_xy,
                            _cur_flux,
                            _psf_evaluator,
                            fit_shape=fit_shape_frame,
                            stamp_size=render_shape_frame,
                            max_iter=max(1, min(2, _INNER_ITERS)),
                            flux_conv=flux_conv_threshold,
                            max_shift=float(fit_shape_frame // 2),
                            group_radius=_group_radius,
                            max_group_size=_max_grp_size,
                            max_grouped_sources=_group_budget,
                            background_rms=float(bkg_std),
                            gain=float(GAIN),
                            initial_positions=_cur_xy,
                            initial_fit_valid=_cur_fit_valid,
                            position_bound=float(fit_shape_frame // 2),
                            position_fixed=True,
                            allow_negative_flux_mask=_cur_forced,
                            log_fn=self._log,
                            stop_fn=lambda: self._stop_requested,
                        )
                        if len(phot_result) == len(_cur_idet):
                            phot_result["iter_detected"] = np.asarray(_cur_idet, dtype=int)
                            phot_result["forced_psf"] = np.asarray(_cur_forced, dtype=bool)
                            phot_result["flags"] = (
                                np.asarray(phot_result["flags"], dtype=np.int32)
                                | np.asarray(_cur_position_flags, dtype=np.int32)
                            )
                        _cur_xy = np.column_stack([
                            np.asarray(phot_result["x_fit"], dtype=float),
                            np.asarray(phot_result["y_fit"], dtype=float),
                        ])
                        _cur_flux = np.asarray(phot_result["flux_fit"], dtype=float)
                        _final_model = _allstar_build_model(
                            img_sub.shape,
                            _cur_xy[:, 0],
                            _cur_xy[:, 1],
                            _cur_flux,
                            _psf_evaluator,
                            render_shape_frame,
                        )
                        _final_residual = _float32_difference(img_sub, _final_model)
                        _save_actual_snapshot(
                            sequence=len(_engine_snapshots) + 1,
                            phase="final_flux",
                            fit_result=phot_result,
                            residual_image=_final_residual,
                            model_image=_final_model,
                            candidate_xy=np.zeros((0, 2), dtype=float),
                            candidate_counts=(0, 0, 0),
                            n_pruned=0,
                            stop_reason=_engine_stop_reason or "final_flux_complete",
                            elapsed_s=time.perf_counter() - _final_started,
                        )
                        model_img = _final_model
                        residual = _final_residual

                        photometry = None
                        fit_fail_reason = None
                        _t["fit1_done"] = time.time()
                        self._log(
                            f"  [TIME] {fname} APEX={_t['fit1_done'] - _t['fit1']:.1f}s "
                            f"n_fit={len(phot_result)} "
                            f"n_new={int(np.sum(np.asarray(phot_result['iter_detected'], dtype=int) > 1))}"
                        )

                        if self._stop_requested:
                            return {"file": fname, "status": "stopped"}
                    else:
                        self.worker_status.emit(wid, fname, "PSF fit", 70)
                        attempt_plan = [False]
                        if use_grouper and _has_grouper:
                            attempt_plan = [True, False]
                        photometry, phot_result, fit_fail_reason = _run_iterative_fit(init_params, "pass1")
                        _t["fit1_done"] = time.time()
                        self._log(f"  [TIME] {fname} pass1={_t['fit1_done'] - _t['fit1']:.1f}s")
                    if psf_fit_engine_cfg != 'apex_iterative' and fit_fail_reason == "stopped":
                        return {"file": fname, "status": "stopped"}
                    refine_pass_max_sources = 2500
                    if psf_fit_engine_cfg != 'apex_iterative' and phot_result is not None and len(phot_result) > 0:
                        refine_init = _results_to_init_params(phot_result, photometry_obj=photometry)
                        if refine_init is not None and len(refine_init) > 0:
                            if len(refine_init) <= refine_pass_max_sources:
                                _skip_pass2 = False
                                if "iter_detected" in phot_result.colnames:
                                    _it_p1 = np.asarray(phot_result["iter_detected"], dtype=float)
                                    _it_p1 = np.where(np.isfinite(_it_p1), _it_p1, 1.0).astype(int)
                                    _n_new_p1 = int(np.sum(_it_p1 > 1))
                                    if _n_new_p1 == 0:
                                        _skip_pass2 = True
                                        self._log("  [PSF] pass2 skipped: no new sources in pass1 (converged)")
                                    elif conv_new_frac > 0:
                                        _new_frac_p1 = float(_n_new_p1) / max(1, len(phot_result))
                                        if _new_frac_p1 <= conv_new_frac:
                                            _skip_pass2 = True
                                            self._log(
                                                f"  [PSF] pass2 skipped: converged "
                                                f"(new_frac={_new_frac_p1:.3f} <= conv_new_frac={conv_new_frac:.3f})"
                                            )
                                if not _skip_pass2:
                                    _t["fit2"] = time.time()
                                    photometry_refine, phot_result_refine, refine_reason = _run_iterative_fit(
                                        refine_init,
                                        "pass2",
                                    )
                                    _t["fit2_done"] = time.time()
                                    self._log(f"  [TIME] {fname} pass2={_t['fit2_done'] - _t['fit2']:.1f}s")
                                    if refine_reason == "stopped":
                                        return {"file": fname, "status": "stopped"}
                                    if phot_result_refine is not None and len(phot_result_refine) > 0:
                                        self._log(
                                            f"  [PSF] refine pass accepted | seed={len(refine_init)} "
                                            f"fit={len(phot_result_refine)}"
                                        )
                                        photometry = photometry_refine
                                        phot_result = phot_result_refine
                                        fit_fail_reason = None
                                    elif refine_reason:
                                        self._log(
                                            f"  [PSF] refine pass failed; keeping pass1 solution | {refine_reason}"
                                        )
                            else:
                                self._log(
                                    f"  [PSF] refine pass skipped | seed={len(refine_init)} "
                                    f"> {refine_pass_max_sources}"
                                )

                    if core_cut.enabled and phot_result is not None and len(phot_result) > 0:
                        try:
                            _x_core = np.asarray(phot_result["x_fit"], dtype=float)
                            _y_core = np.asarray(phot_result["y_fit"], dtype=float)
                            _keep_core_result = _core_keep(np.column_stack([_x_core, _y_core]))
                            n_core_excluded_result = int(np.sum(~_keep_core_result))
                            if n_core_excluded_result > 0:
                                phot_result = phot_result[_keep_core_result]
                                self._log(
                                    f"  [CORE] fitted result rows excluded: {n_core_excluded_result} "
                                    f"(r<{core_cut.radius_px:.1f}px)"
                                )
                        except Exception as _core_e:
                            self._log(f"  [CORE] result filter skipped: {_core_e}")

                    raw_iter_counts: dict[int, int] = {}
                    n_new_raw_total = 0
                    n_new_kept_total = 0
                    raw_new_xy = np.zeros((0, 2), dtype=float)

                    if phot_result is not None and len(phot_result) > 0 and "iter_detected" in phot_result.colnames:
                        try:
                            _x0 = np.asarray(phot_result["x_fit"], dtype=float)
                            _y0 = np.asarray(phot_result["y_fit"], dtype=float)
                            _it_raw0 = np.asarray(phot_result["iter_detected"], dtype=float)
                            _it0 = np.where(np.isfinite(_it_raw0), _it_raw0, 1.0).astype(int)
                            _uniq_raw, _cnt_raw = np.unique(_it0[_it0 > 1], return_counts=True) if np.any(_it0 > 1) else ([], [])
                            raw_iter_counts = {int(i): int(c) for i, c in zip(_uniq_raw, _cnt_raw)}
                            n_new_raw_total = int(np.sum(_it0 > 1))
                            n_new_kept_total = n_new_raw_total
                            _m_raw_xy = np.isfinite(_x0) & np.isfinite(_y0) & (_it0 > 1)
                            if np.any(_m_raw_xy):
                                raw_new_xy = np.column_stack([_x0[_m_raw_xy], _y0[_m_raw_xy]])
                                _n_show = min(6, len(raw_new_xy))
                                _pts = ", ".join(
                                    [f"({raw_new_xy[i,0]:.2f},{raw_new_xy[i,1]:.2f})" for i in range(_n_show)]
                                )
                                self._log(
                                    f"  [RAWXY] iter>1 raw first={_n_show}/{len(raw_new_xy)} | {_pts}"
                                )
                                try:
                                    _tree_seed = cKDTree(np.asarray(xy_det, dtype=float))
                                    _d_seed, _ = _tree_seed.query(raw_new_xy, k=1, workers=1)
                                    _seed_tol = 1.0  # px
                                    _n_near = int(np.sum(np.asarray(_d_seed, dtype=float) <= _seed_tol))
                                    _n_far = int(len(_d_seed) - _n_near)
                                    self._log(
                                        f"  [RAWXY] vs Step4 seed | near<=1.00px={_n_near} | far={_n_far}"
                                    )
                                except Exception:
                                    pass
                        except Exception:
                            pass

                    if psf_fit_engine_cfg == "apex_iterative" and _engine_snapshots:
                        n_new_raw_total = int(sum(
                            int(record.get("n_candidates_raw", 0))
                            for record in _engine_snapshots
                            if record.get("phase") == "residual_fit"
                        ))

                    # De-duplicate residual re-detections against iter1 fitted sources first.
                    # ALLSTAR outer loop already deduped at detection time; skip here.
                    if (
                        psf_fit_engine_cfg != 'apex_iterative'
                        and phot_result is not None
                        and len(phot_result) > 0
                        and "iter_detected" in phot_result.colnames
                        and dedup_enabled
                    ):
                        try:
                            _x = np.asarray(phot_result["x_fit"], dtype=float)
                            _y = np.asarray(phot_result["y_fit"], dtype=float)
                            _it_raw = np.asarray(phot_result["iter_detected"], dtype=float)
                            _it = np.where(np.isfinite(_it_raw), _it_raw, 1.0).astype(int)
                            _finite = np.isfinite(_x) & np.isfinite(_y)
                            _idx_base = np.where(_finite & (_it <= 1))[0]
                            _idx_new = np.where(_finite & (_it > 1))[0]
                            if len(_idx_base) and len(_idx_new):
                                _xy_base = np.column_stack([_x[_idx_base], _y[_idx_base]])
                                _xy_new = np.column_stack([_x[_idx_new], _y[_idx_new]])
                                _tree = cKDTree(_xy_base)
                                _dnn, _ = _tree.query(_xy_new, k=1, workers=1)
                                if len(_dnn):
                                    try:
                                        _dnn_arr = np.asarray(_dnn, dtype=float)
                                        self._log(
                                            "  [DEDUP] d_nn(iter2->iter1) px | "
                                            f"p50={np.nanpercentile(_dnn_arr, 50):.2f} "
                                            f"p90={np.nanpercentile(_dnn_arr, 90):.2f} "
                                            f"p99={np.nanpercentile(_dnn_arr, 99):.2f}"
                                        )
                                    except Exception:
                                        pass
                                if np.isfinite(duplicate_radius_px_cfg):
                                    _dup_r_px = float(duplicate_radius_px_cfg)
                                else:
                                    _dup_r_px = float(max(0.0, duplicate_radius_mult * fwhm_safe))
                                _keep_new = np.asarray(_dnn, dtype=float) > _dup_r_px
                                if np.any(~_keep_new):
                                    _drop_n = int(np.sum(~_keep_new))
                                    _keep_mask = np.ones(len(phot_result), dtype=bool)
                                    _keep_mask[_idx_new[~_keep_new]] = False
                                    phot_result = phot_result[_keep_mask]
                                    _it2_raw = np.asarray(phot_result["iter_detected"], dtype=float)
                                    _it2 = np.where(np.isfinite(_it2_raw), _it2_raw, 1.0).astype(int)
                                    n_new_kept_total = int(np.sum(_it2 > 1))
                                    self._log(
                                        f"  [DEDUP] dropped near-duplicate iter>1 sources: {_drop_n} "
                                        f"(r<{_dup_r_px:.2f}px)"
                                    )
                        except Exception as _de:
                            self._log(f"  [DEDUP] skipped: {_de}")

                    # Apply residual new-source cap.
                    # ALLSTAR outer loop already applied per-cycle caps; skip here.
                    if (
                        psf_fit_engine_cfg != 'apex_iterative'
                        and phot_result is not None
                        and len(phot_result) > 0
                        and "iter_detected" in phot_result.colnames
                    ):
                        try:
                            _cols_cap = list(phot_result.colnames)
                            _it_raw = np.asarray(phot_result["iter_detected"], dtype=float)
                            _it = np.where(np.isfinite(_it_raw), _it_raw, 1.0).astype(int)
                            _new_now = int(np.sum(_it > 1))

                            _cap_abs = int(new_sources_cap_per_iter) if int(new_sources_cap_per_iter) > 0 else None
                            _cap_frac_n = None
                            if float(new_sources_cap_frac) > 0:
                                _cap_frac_n = int(np.floor(float(new_sources_cap_frac) * max(1, len(init_params))))

                            if _cap_abs is not None and _cap_frac_n is not None:
                                _cap_new = min(_cap_abs, _cap_frac_n)
                            elif _cap_abs is not None:
                                _cap_new = _cap_abs
                            else:
                                _cap_new = _cap_frac_n

                            if _cap_new is not None:
                                _cap_new = max(0, int(_cap_new))
                                if _new_now > _cap_new:
                                    _ff_col_cap = next((c for c in ("flux_fit", "flux") if c in _cols_cap), None)
                                    _flux_all = (
                                        np.asarray(phot_result[_ff_col_cap], dtype=float)
                                        if _ff_col_cap is not None else
                                        np.full(len(phot_result), np.nan, dtype=float)
                                    )
                                    _idx_new = np.where(_it > 1)[0]
                                    if len(_idx_new):
                                        _m = _flux_all[_idx_new]
                                        _m = np.where(np.isfinite(_m), _m, -np.inf)
                                        _order = np.argsort(_m)[::-1]
                                        _keep_new_idx = _idx_new[_order[:_cap_new]]
                                        _keep_mask = (_it <= 1)
                                        _keep_mask[_keep_new_idx] = True
                                        phot_result = phot_result[_keep_mask]
                                        _it_kept_raw = np.asarray(phot_result["iter_detected"], dtype=float)
                                        _it_kept = np.where(np.isfinite(_it_kept_raw), _it_kept_raw, 1.0).astype(int)
                                        _new_after_cap = int(np.sum(_it_kept > 1))
                                        self._log(
                                            f"  [CAP] new sources capped | raw={n_new_raw_total} dedup={_new_now} kept={_new_after_cap} "
                                            f"(cap={_cap_new}, abs={new_sources_cap_per_iter}, frac={new_sources_cap_frac:.3f})"
                                        )
                                        n_new_kept_total = _new_after_cap
                                else:
                                    n_new_kept_total = _new_now
                            else:
                                n_new_kept_total = _new_now
                        except Exception as _ce:
                            self._log(f"  [CAP] cap logic skipped: {_ce}")

                    # ── Diagnostics ──────────────────────────────────────────────────────
                    if phot_result is not None and len(phot_result) > 0:
                        try:
                            _cols = list(phot_result.colnames)
                            _ff_col = next((c for c in ("flux_fit", "flux") if c in _cols), None)
                            if _ff_col:
                                _ff = np.asarray(phot_result[_ff_col], dtype=float)
                                _ff_pos = _ff[np.isfinite(_ff) & (_ff > 0)]
                                self._log(
                                    f"  [DIAG] {fname} | n_fit={len(_ff)} | "
                                    f"flux_fit: n>0={len(_ff_pos)} "
                                    f"med={np.nanmedian(_ff):.2f} max={np.nanmax(_ff):.2f} | "
                                    f"img_sub peak={float(np.nanmax(img_sub)):.2f} bkg_std={float(bkg_std):.3f}"
                                )
                            if "group_size" in _cols:
                                _gs = np.asarray(phot_result["group_size"], dtype=int)
                                self._log(
                                    f"  [DIAG] group_size: max={_gs.max()} "
                                    f"med={np.median(_gs):.0f} "
                                    f"n_groups={len(np.unique(phot_result['group_id']))}"
                                )
                            if "iter_detected" in _cols:
                                _idet_raw = np.asarray(phot_result["iter_detected"], dtype=float)
                                _idet = np.where(np.isfinite(_idet_raw), _idet_raw, 1.0).astype(int)
                                self._log(
                                    f"  [DIAG] iters used={_idet.max()} | "
                                    f"new sources (iter>1)={int(np.sum(_idet > 1))}"
                                )
                        except Exception as _de:
                            self._log(f"  [DIAG] diag error: {_de}")

                    # ── Model image & residual ────────────────────────────────────────────
                    if model_img is None:
                        if phot_result is not None and len(phot_result) > 0:
                            if psf_fit_engine_cfg == 'apex_iterative' and _psf_evaluator is not None:
                                _x_f = np.asarray(phot_result["x_fit"],   dtype=float)
                                _y_f = np.asarray(phot_result["y_fit"],   dtype=float)
                                _fl_f = np.asarray(phot_result["flux_fit"], dtype=float)
                                _v = np.isfinite(_x_f) & np.isfinite(_y_f) & np.isfinite(_fl_f) & (_fl_f > 0)
                                model_img = _allstar_build_model(
                                    img_sub.shape, _x_f[_v], _y_f[_v], _fl_f[_v],
                                    _psf_evaluator, render_shape_frame,
                                )
                            else:
                                model_img = _render_model_from_results(
                                    phot_result, photometry_obj=photometry
                                )
                        residual = (
                            _float32_difference(img_sub, model_img)
                            if model_img is not None
                            else img_sub.copy()
                        )
                    if model_img is not None:
                        self._log(
                            f"  [DIAG] model_img sum={float(np.nansum(model_img)):.2f} "
                            f"peak={float(np.nanmax(model_img)):.2f} | "
                            f"img_sub peak={float(np.nanmax(img_sub)):.2f} | "
                            f"subtract_shape={render_shape_frame}"
                        )

                    res_std = _fast_res_std(residual)

                    # n_new_total: kept sources first detected in iteration > 1
                    n_new_total = int(n_new_kept_total)
                    if n_new_total <= 0 and phot_result is not None and "iter_detected" in phot_result.colnames:
                        _iter_raw = np.asarray(phot_result["iter_detected"], dtype=float)
                        _iter_safe = np.where(np.isfinite(_iter_raw), _iter_raw, 1.0).astype(int)
                        n_new_total = int(np.sum(_iter_safe > 1))

                    # ── Save starsub / residual for cutout viewer ─────────────────────────
                    fit_xy = np.zeros((0, 2), dtype=float)
                    fit_flux = np.zeros((0,), dtype=float)
                    fit_iter = np.zeros((0,), dtype=int)
                    n_fit = 0
                    if phot_result is not None and len(phot_result) > 0:
                        try:
                            x_it = np.asarray(phot_result["x_fit"], dtype=float)
                            y_it = np.asarray(phot_result["y_fit"], dtype=float)
                            if "flux_fit" in phot_result.colnames:
                                f_it = np.asarray(phot_result["flux_fit"], dtype=float)
                            elif "flux" in phot_result.colnames:
                                f_it = np.asarray(phot_result["flux"], dtype=float)
                            else:
                                f_it = np.full(len(x_it), np.nan, dtype=float)
                            if "iter_detected" in phot_result.colnames:
                                it_raw = np.asarray(phot_result["iter_detected"], dtype=float)
                                it_safe = np.where(np.isfinite(it_raw), it_raw, 1.0).astype(int)
                            else:
                                it_safe = np.ones(len(x_it), dtype=int)
                            valid_it = np.isfinite(x_it) & np.isfinite(y_it) & np.isfinite(f_it)
                            if np.any(valid_it):
                                fit_xy = np.column_stack([x_it[valid_it], y_it[valid_it]])
                                fit_flux = f_it[valid_it]
                                fit_iter = it_safe[valid_it]
                            n_fit = int(len(fit_xy))
                        except Exception:
                            pass

                    starsub_raw = None
                    init_xy_ui = np.column_stack(
                        [np.asarray(init_params["x_0"], dtype=float), np.asarray(init_params["y_0"], dtype=float)]
                    ) if len(init_params) > 0 else np.zeros((0, 2), dtype=float)
                    iter_max_used = int(np.max(fit_iter)) if len(fit_iter) else 1
                    iter_max_used = max(1, min(iter_max_used, max(1, int(max_iter))))

                    iter_records = (
                        list(_engine_snapshots)
                        if psf_fit_engine_cfg == "apex_iterative" and _engine_snapshots
                        else []
                    )
                    if not iter_records:
                        starsub_raw = (
                            _float32_difference(img, model_img)
                            if model_img is not None
                            else np.asarray(img_sub + float(bkg_med), dtype=np.float32)
                        )
                    if len(fit_iter):
                        try:
                            _uniq, _cnt = np.unique(fit_iter, return_counts=True)
                            _iter_counts = ", ".join([f"i{int(i)}={int(c)}" for i, c in zip(_uniq, _cnt)])
                            self._log(f"  [DIAG] iter source counts: {_iter_counts}")
                        except Exception:
                            pass

                    def _render_model_subset(xy_sub: np.ndarray, flux_sub: np.ndarray) -> np.ndarray:
                        if len(xy_sub) == 0 or len(flux_sub) == 0:
                            return np.zeros_like(img_sub, dtype=np.float32)
                        try:
                            from photutils.datasets import make_model_image as _make_model_image
                            xy_sub = np.asarray(xy_sub, dtype=float)
                            flux_sub = np.asarray(flux_sub, dtype=float)
                            valid = (
                                np.isfinite(xy_sub[:, 0]) &
                                np.isfinite(xy_sub[:, 1]) &
                                np.isfinite(flux_sub) &
                                (flux_sub > 0)
                            )
                            if not np.any(valid):
                                return np.zeros_like(img_sub, dtype=np.float32)
                            xy_sub = xy_sub[valid]
                            flux_sub = flux_sub[valid]

                            # Use a wider rendering footprint than fit window.
                            # fit_shape covers only the PSF core; residual subtraction
                            # needs the wings too.  Use 2× epsf_size so the rendered
                            # stamp captures flux out to ~4×FWHM from each source.
                            pt = AstropyTable()
                            pt["x_0"] = np.asarray(xy_sub[:, 0], dtype=float)
                            pt["y_0"] = np.asarray(xy_sub[:, 1], dtype=float)
                            pt["flux"] = np.asarray(flux_sub, dtype=float)
                            mod = _clone_psf_model(epsf_model)
                            out = _make_model_image(
                                img_sub.shape,
                                mod,
                                pt,
                                model_shape=(int(render_shape_frame), int(render_shape_frame)),
                                x_name="x_0",
                                y_name="y_0",
                            )
                            return np.asarray(out, dtype=np.float32)
                        except Exception as _re:
                            self._log(f"  [DIAG] iter model render failed: {_re}")
                            # Never return full final model here; this helper is for
                            # subset-by-iteration rendering used by diagnostics.
                            return np.zeros_like(img_sub, dtype=np.float32)

                    _reconstructed_iters = (
                        [] if iter_records else range(1, iter_max_used + 1)
                    )
                    for it_no in _reconstructed_iters:
                        m_le = fit_iter <= it_no if len(fit_iter) else np.zeros((0,), dtype=bool)
                        m_eq = fit_iter == it_no if len(fit_iter) else np.zeros((0,), dtype=bool)
                        fit_xy_i = fit_xy[m_le] if len(fit_xy) else np.zeros((0, 2), dtype=float)
                        fit_flux_i = fit_flux[m_le] if len(fit_flux) else np.zeros((0,), dtype=float)
                        if it_no <= 1:
                            applied_xy_i = init_xy_ui
                            det_xy_i = np.zeros((0, 2), dtype=float)
                        else:
                            applied_xy_i = fit_xy[fit_iter < it_no] if len(fit_iter) else np.zeros((0, 2), dtype=float)
                            det_xy_i = fit_xy[m_eq] if len(fit_xy) else np.zeros((0, 2), dtype=float)

                        if len(applied_xy_i) and len(det_xy_i):
                            box_xy_i = np.vstack([applied_xy_i, det_xy_i])
                        elif len(applied_xy_i):
                            box_xy_i = applied_xy_i
                        elif len(det_xy_i):
                            box_xy_i = det_xy_i
                        else:
                            box_xy_i = fit_xy_i

                        # Both residual_i and starsub_i are derived from the EXACT
                        # photometry.make_model_image() result where available.
                        #
                        # Last iter  → use exact model_img directly.
                        # Earlier iter → exact_full_model minus later-iter contributions
                        #   (add back later-iter _render_model_subset so only iter1..N remain).
                        # Fallback (model_img is None) → _render_model_subset only.
                        if model_img is not None:
                            if it_no == iter_max_used:
                                # Exact: photometry.make_model_image() covers all fitted sources
                                residual_i = np.asarray(residual, dtype=np.float32)
                                starsub_i = np.asarray(starsub_raw, dtype=np.float32)
                                res_std_i = float(res_std)
                            else:
                                later_mask = (
                                    (fit_iter > it_no) if len(fit_iter) > 0
                                    else np.zeros(0, dtype=bool)
                                )
                                if np.any(later_mask):
                                    later_contrib = _render_model_subset(
                                        fit_xy[later_mask], fit_flux[later_mask]
                                    )
                                    residual_i = np.asarray(
                                        residual + later_contrib, dtype=np.float32
                                    )
                                    starsub_i = np.asarray(
                                        starsub_raw + later_contrib, dtype=np.float32
                                    )
                                else:
                                    residual_i = np.asarray(residual, dtype=np.float32)
                                    starsub_i = np.asarray(starsub_raw, dtype=np.float32)
                                res_std_i = _fast_res_std(residual_i)
                        else:
                            model_i = _render_model_subset(fit_xy_i, fit_flux_i)
                            residual_i = (img_sub - model_i).astype(np.float32, copy=False)
                            starsub_i = (img - model_i).astype(np.float32, copy=False)
                            res_std_i = _fast_res_std(residual_i)

                        hdr_it = fits.Header()
                        hdr_it["FILTER"] = this_filter
                        hdr_it["BKGMED"] = float(bkg_med)
                        hdr_it["ITER"] = int(it_no)
                        residual_iter_path = output_dir / f"residual_iter{it_no}_{fname}"
                        starsub_iter_path = output_dir / f"starsub_iter{it_no}_{fname}"
                        fitxy_iter_path = output_dir / f"fitxy_iter{it_no}_{fname}.npy"
                        modelxy_iter_path = output_dir / f"modelxy_iter{it_no}_{fname}.npy"
                        appliedxy_iter_path = output_dir / f"appliedxy_iter{it_no}_{fname}.npy"
                        detxy_iter_path = output_dir / f"detxy_iter{it_no}_{fname}.npy"
                        boxxy_iter_path = output_dir / f"boxxy_iter{it_no}_{fname}.npy"
                        _is_final_iter = (it_no == iter_max_used)
                        _is_first_iter = (it_no == 1)
                        _write_fits = save_all_iter_residuals or _is_final_iter or _is_first_iter
                        if _write_fits:
                            # Rice-compressed FITS — 3-5× smaller, transparent to readers
                            _chdu_res = fits.CompImageHDU(
                                residual_i.astype(np.float32), hdr_it,
                                compression_type='RICE_1', quantize_level=16.0,
                            )
                            _chdu_res.writeto(str(residual_iter_path), overwrite=True)
                            # starsub only for final (or all-iters mode) — saves ~50% disk
                            if _is_final_iter or save_all_iter_residuals:
                                _chdu_sub = fits.CompImageHDU(
                                    starsub_i.astype(np.float32), hdr_it,
                                    compression_type='RICE_1', quantize_level=16.0,
                                )
                                _chdu_sub.writeto(str(starsub_iter_path), overwrite=True)
                        np.save(str(fitxy_iter_path), np.asarray(fit_xy_i, dtype=np.float32))
                        # modelxy: sources first detected at THIS iter only (not cumulative).
                        # iter1 → initial seeds (1023), iter2 → new residual detections (70).
                        fit_xy_this = fit_xy[m_eq] if len(fit_xy) else np.zeros((0, 2), dtype=float)
                        fit_flux_this = fit_flux[m_eq] if len(fit_flux) else np.zeros((0,), dtype=float)
                        if it_no == 1:
                            # iter1: "new" mask is empty; use all iter1 sources
                            fit_xy_this = fit_xy[fit_iter == 1] if len(fit_iter) else np.zeros((0, 2), dtype=float)
                            fit_flux_this = fit_flux[fit_iter == 1] if len(fit_iter) else np.zeros((0,), dtype=float)
                        _m_model = (
                            np.isfinite(fit_xy_this[:, 0]) &
                            np.isfinite(fit_xy_this[:, 1]) &
                            np.isfinite(fit_flux_this) &
                            (fit_flux_this > 0)
                        ) if len(fit_xy_this) else np.zeros((0,), dtype=bool)
                        model_xy_i = fit_xy_this[_m_model] if len(fit_xy_this) else np.zeros((0, 2), dtype=float)
                        np.save(str(modelxy_iter_path), np.asarray(model_xy_i, dtype=np.float32))
                        np.save(str(appliedxy_iter_path), np.asarray(applied_xy_i, dtype=np.float32))
                        np.save(str(detxy_iter_path), np.asarray(det_xy_i, dtype=np.float32))
                        np.save(str(boxxy_iter_path), np.asarray(box_xy_i, dtype=np.float32))

                        n_new_i_kept = int(np.sum(m_eq)) if it_no > 1 else 0
                        n_new_i_raw = int(raw_iter_counts.get(int(it_no), n_new_i_kept)) if it_no > 1 else 0
                        iter_records.append({
                            "iter": int(it_no),
                            "fit_shape_px": _to_int(fit_shape_frame, 9),
                            "epsf_size_px": _to_int(epsf_size_frame, 25),
                            "n_fit": int(len(fit_xy_i)),
                            "residual_std": float(res_std_i),
                            "n_new_raw": int(n_new_i_raw),
                            "n_new_kept": int(n_new_i_kept),
                            "n_applied_prev": int(len(applied_xy_i)),
                            "residual_path": residual_iter_path.name if _write_fits else None,
                            "starsub_path": starsub_iter_path.name if (_is_final_iter or save_all_iter_residuals) else None,
                            "fitxy_path": fitxy_iter_path.name,
                            "modelxy_path": modelxy_iter_path.name,
                            "detxy_path": detxy_iter_path.name,
                            "appliedxy_path": appliedxy_iter_path.name,
                            "boxxy_path": boxxy_iter_path.name,
                        })
                    _t["done"] = time.time()
                    _t_bkg   = _t.get("bkg_done", _t.get("bkg", _t["done"])) - _t.get("bkg", _t["done"])
                    _t_epsf  = _t.get("epsf_done", _t.get("epsf", _t["done"])) - _t.get("epsf", _t.get("bkg_done", _t["done"]))
                    _t_sub   = _t.get("fit1", _t["done"]) - _t.get("substar", _t.get("fit1", _t["done"]))
                    _t_p1    = _t.get("fit1_done", _t["done"]) - _t.get("fit1", _t["done"])
                    _t_p2    = _t.get("fit2_done", _t.get("fit2", _t["done"])) - _t.get("fit2", _t.get("fit2_done", _t["done"]))
                    _t_total = _t["done"] - _t["start"]
                    self._log(
                        f"  [TIME] {fname} total={_t_total:.1f}s | "
                        f"bkg={_t_bkg:.1f}s epsf={_t_epsf:.1f}s substar={_t_sub:.1f}s "
                        f"pass1={_t_p1:.1f}s pass2={_t_p2:.1f}s"
                    )
                    self._log(
                        f"  fit done | n_fit={n_fit} | n_new={n_new_total} | "
                        f"residual_std={res_std:.4f}"
                    )

                    if (phot_result is None) or (len(phot_result) == 0):
                        reason = fit_fail_reason or "no fitted sources"
                        if n_init_drop > 0:
                            reason = f"{reason} | dropped_edge_init={n_init_drop}"
                        self.worker_status.emit(wid, fname, "Fit failed", 100)
                        return {
                            "file": fname,
                            "status": "fit_failed",
                            "reason": reason,
                        }

                    phot_rows = []
                    if phot_result is not None:
                        x_fit = np.array(phot_result["x_fit"])
                        y_fit = np.array(phot_result["y_fit"])
                        _ff_col_main = next(
                            (c for c in ("flux_fit", "flux", "flux_0") if c in phot_result.colnames), None
                        )
                        flux_fit = (
                            np.array(phot_result[_ff_col_main], dtype=float)
                            if _ff_col_main is not None
                            else np.full(len(x_fit), np.nan, dtype=float)
                        )
                        flux_err = (np.array(phot_result["flux_err"]) if "flux_err" in phot_result.colnames else np.full(len(x_fit), np.nan))
                        qfit_col = (
                            np.array(phot_result["qfit"], dtype=float)
                            if "qfit" in phot_result.colnames else np.full(len(x_fit), np.nan)
                        )
                        cfit_col = (
                            np.array(phot_result["cfit"], dtype=float)
                            if "cfit" in phot_result.colnames else np.full(len(x_fit), np.nan)
                        )
                        redchi_col = (
                            np.array(phot_result["reduced_chi2"], dtype=float)
                            if "reduced_chi2" in phot_result.colnames else np.full(len(x_fit), np.nan)
                        )
                        n_pixels_col = (
                            np.array(phot_result["n_pixels_fit"], dtype=int)
                            if "n_pixels_fit" in phot_result.colnames else np.zeros(len(x_fit), dtype=int)
                        )
                        flags_col = (
                            np.array(phot_result["flags"], dtype=int)
                            if "flags" in phot_result.colnames else np.zeros(len(x_fit), dtype=int)
                        )
                        forced_col = (
                            np.array(phot_result["forced_psf"], dtype=bool)
                            if "forced_psf" in phot_result.colnames
                            else np.zeros(len(x_fit), dtype=bool)
                        )

                        fit_neighbor_dist = np.full(len(x_fit), np.inf, dtype=float)
                        finite_fit_xy = np.isfinite(x_fit) & np.isfinite(y_fit)
                        finite_indices = np.flatnonzero(finite_fit_xy)
                        if len(finite_indices) >= 2:
                            fit_tree = cKDTree(np.column_stack([x_fit[finite_fit_xy], y_fit[finite_fit_xy]]))
                            fit_distances, _ = fit_tree.query(
                                np.column_stack([x_fit[finite_fit_xy], y_fit[finite_fit_xy]]),
                                k=2,
                                workers=1,
                            )
                            fit_neighbor_dist[finite_indices] = np.asarray(fit_distances[:, 1], dtype=float)
                        unresolved_threshold_px = _UNRESOLVED_NEIGHBOR_FWHM * float(fwhm_safe)
                        crowding_unreliable = (
                            np.isfinite(fit_neighbor_dist)
                            & (fit_neighbor_dist < unresolved_threshold_px)
                        )
                        flags_col[crowding_unreliable] |= int(PSFFitFlag.CROWDING_UNRELIABLE)

                        valid_fit_xy = np.isfinite(x_fit) & np.isfinite(y_fit)
                        for k in np.where(valid_fit_xy)[0]:
                            xi = int(round(float(x_fit[k])))
                            yi = int(round(float(y_fit[k])))
                            if 0 <= xi < w and 0 <= yi < h and img[yi, xi] >= sat_adu:
                                flags_col[k] |= self.FLAG_SAT
                            edge_m = fit_shape_frame // 2 + 1
                            if xi < edge_m or xi >= w - edge_m or yi < edge_m or yi >= h - edge_m:
                                flags_col[k] |= self.FLAG_EDGE

                        if len(xy_det) and len(x_fit):
                            src_xy = np.column_stack([x_fit, y_fit])
                            tree_ref = cKDTree(xy_det)
                            matched_det_uids = np.full(len(x_fit), -1, dtype=int)
                            nn_dists = np.full(len(x_fit), np.inf, dtype=float)
                            valid_src = np.isfinite(src_xy[:, 0]) & np.isfinite(src_xy[:, 1])
                            if np.any(valid_src):
                                q_dists, q_idx = tree_ref.query(src_xy[valid_src], k=1, workers=1)
                                matched_det_uids[valid_src] = det_uids[q_idx]
                                nn_dists[valid_src] = q_dists
                            match_tol = 2.0 * fwhm_med
                        else:
                            matched_det_uids = np.arange(len(x_fit), dtype=int)
                            nn_dists = np.zeros(len(x_fit))
                            match_tol = np.inf

                        _psf_only_uid = -1  # counts down for sources with no step4 match
                        _used_det_uids = set()
                        _uid_collision = 0
                        for k in range(len(x_fit)):
                            xk = float(x_fit[k]) if np.isfinite(x_fit[k]) else np.nan
                            yk = float(y_fit[k]) if np.isfinite(y_fit[k]) else np.nan
                            if not (np.isfinite(xk) and np.isfinite(yk)):
                                continue
                            fe = float(flux_fit[k]) * GAIN  # ADU → electrons (same as step5)
                            se = float(flux_err[k]) * GAIN if np.isfinite(flux_err[k]) else np.nan
                            snr = fe / se if (np.isfinite(se) and se > 0) else np.nan
                            _qfit_expected, _qfit_noise_ratio = qfit_noise_diagnostics(
                                float(qfit_col[k]),
                                float(n_pixels_col[k]),
                                float(snr),
                                float(psf_nea_frame),
                            )
                            qfit_expected = float(_qfit_expected)
                            qfit_noise_ratio = float(_qfit_noise_ratio)
                            if np.isfinite(snr) and snr >= min_snr and fe > 0:
                                mag_psf = ZP - 2.5 * np.log10(max(fe, 1e-30) / exptime)
                                mag_psf_err = (2.5 / np.log(10) * se / fe if (np.isfinite(se) and fe > 0) else np.nan)
                            else:
                                mag_psf = np.nan
                                mag_psf_err = np.nan
                            # Assign unique negative UIDs for PSF-only new detections
                            # (no matching step4 source within match_tol).
                            # Negative UIDs are excluded by downstream steps that join
                            # on step4 det_uid; they are preserved for traceability.
                            if nn_dists[k] <= match_tol:
                                cand_uid = int(matched_det_uids[k])
                                # Keep det_uid unique per frame: when two PSF components
                                # map to the same Step4 seed, keep first as seed UID and
                                # force others to PSF-only negative UIDs.
                                if cand_uid not in _used_det_uids:
                                    det_uid = cand_uid
                                    _used_det_uids.add(cand_uid)
                                else:
                                    _uid_collision += 1
                                    det_uid = _psf_only_uid
                                    _psf_only_uid -= 1
                            else:
                                det_uid = _psf_only_uid
                                _psf_only_uid -= 1
                                cand_uid = -1
                            if "iter_detected" in phot_result.colnames:
                                iter_val = _safe_float(phot_result["iter_detected"][k], np.nan)
                                iter_found = int(iter_val) if np.isfinite(iter_val) and iter_val > 0 else 1
                            else:
                                iter_found = 1
                            r_core_px = (
                                float(np.hypot(xk - core_cut.center_x, yk - core_cut.center_y))
                                if np.isfinite(core_cut.center_x) and np.isfinite(core_cut.center_y)
                                else np.nan
                            )
                            phot_rows.append({
                                "det_uid": det_uid,
                                "seed_uid": int(cand_uid) if np.isfinite(cand_uid) else -1,
                                "x_fit": round(xk, 4),
                                "y_fit": round(yk, 4),
                                "FILTER": this_filter,
                                "flux_psf_e": round(fe, 4) if np.isfinite(fe) else np.nan,
                                "flux_psf_err_e": round(float(se), 4) if np.isfinite(se) else np.nan,
                                "gain_e_per_adu": round(float(GAIN), 8),
                                "rdnoise_e": round(float(rn_e), 6),
                                "binning_x": int(noise.bin_x),
                                "binning_y": int(noise.bin_y),
                                "gain_source": noise.gain_source,
                                "rdnoise_source": noise.rdnoise_source,
                                "mag_psf": round(mag_psf, 6) if np.isfinite(mag_psf) else np.nan,
                                "mag_psf_err": round(mag_psf_err, 6) if np.isfinite(mag_psf_err) else np.nan,
                                "snr_psf": round(float(snr), 3) if np.isfinite(snr) else np.nan,
                                "qfit": round(float(qfit_col[k]), 6) if np.isfinite(qfit_col[k]) else np.nan,
                                "qfit_noise_expected": (
                                    round(float(qfit_expected), 6)
                                    if np.isfinite(qfit_expected) else np.nan
                                ),
                                "qfit_noise_ratio": (
                                    round(float(qfit_noise_ratio), 6)
                                    if np.isfinite(qfit_noise_ratio) else np.nan
                                ),
                                "cfit": round(float(cfit_col[k]), 6) if np.isfinite(cfit_col[k]) else np.nan,
                                "reduced_chi2": (
                                    round(float(redchi_col[k]), 6)
                                    if np.isfinite(redchi_col[k]) else np.nan
                                ),
                                "n_pixels_fit": int(n_pixels_col[k]),
                                "psf_nea_px": (
                                    round(float(psf_nea_frame), 4)
                                    if np.isfinite(psf_nea_frame) else np.nan
                                ),
                                "fit_window_px": int(fit_shape_frame),
                                "fit_window_energy": (
                                    round(float(fit_window_plan.energy_fraction), 6)
                                    if np.isfinite(fit_window_plan.energy_fraction)
                                    else np.nan
                                ),
                                "iter_found": iter_found,
                                "forced_psf": bool(forced_col[k]),
                                "neighbor_dist_px": (
                                    round(float(fit_neighbor_dist[k]), 4)
                                    if np.isfinite(fit_neighbor_dist[k]) else np.nan
                                ),
                                "neighbor_dist_fwhm": (
                                    round(float(fit_neighbor_dist[k]) / max(float(fwhm_safe), 1e-6), 4)
                                    if np.isfinite(fit_neighbor_dist[k]) else np.nan
                                ),
                                "crowding_unreliable_psf": bool(crowding_unreliable[k]),
                                "flags_psf": int(flags_col[k]),
                                "saturated_psf": bool(int(flags_col[k]) & self.FLAG_SAT),
                                "edge_psf": bool(int(flags_col[k]) & self.FLAG_EDGE),
                                "psf_core_r_px": round(r_core_px, 3) if np.isfinite(r_core_px) else np.nan,
                                "psf_core_cut_px": (
                                    round(float(core_cut.radius_px), 3)
                                    if np.isfinite(core_cut.radius_px)
                                    else np.nan
                                ),
                                "exptime": round(exptime, 4),
                            })
                        if _uid_collision > 0:
                            self._log(
                                f"  [UID] det_uid collision resolved: {_uid_collision} "
                                f"(assigned PSF-only negative det_uid)"
                            )

                    df_out = pd.DataFrame(phot_rows)

                    scale_result = PSFApertureScale(
                        scale=1.0,
                        applied=False,
                        n_matched=0,
                        n_candidates=0,
                        n_used=0,
                        median_delta_mag_raw=np.nan,
                        scatter_mag=np.nan,
                        reason="disabled",
                    )
                    scale_references = pd.DataFrame()
                    flux_scale_reference_name = ""
                    if flux_scale_correction and not df_out.empty:
                        scale_result, scale_references = estimate_psf_aperture_scale(
                            df_out,
                            df_ap,
                            min_snr=flux_scale_min_snr,
                            min_stars=flux_scale_min_stars,
                            min_neighbor_fwhm=flux_scale_min_neighbor_fwhm,
                            max_scatter_mag=flux_scale_max_scatter_mag,
                        )
                    if not df_out.empty:
                        df_out = apply_psf_aperture_scale(
                            df_out,
                            scale_result,
                            zeropoint=ZP,
                            exptime=exptime,
                        )
                    if not scale_references.empty:
                        flux_scale_reference_name = f"psf_flux_scale_reference_{fname}.csv"
                        scale_references.to_csv(
                            output_dir / flux_scale_reference_name,
                            index=False,
                        )
                    scale_level = (
                        "WARN][PSF-SCALE"
                        if flux_scale_correction and not scale_result.applied
                        else "PSF-SCALE"
                    )
                    self._log(
                        f"[{scale_level}] frame={fname} scale={scale_result.scale:.6f} "
                        f"refs={scale_result.n_used}/{scale_result.n_candidates} "
                        f"scatter={scale_result.scatter_mag:.4f} mag "
                        f"status={scale_result.reason}"
                    )

                    # ── Flux unit sanity check (P1-2) ────────────────────────
                    # PSF fitting runs on img_sub (ADU); flux_fit is in ADU.
                    # flux_psf_e = flux_fit * GAIN.  If the ratio PSF/aperture
                    # deviates far from 1.0 across bright sources, GAIN or the
                    # aperture data may be in the wrong unit.
                    if flux_init_map and len(df_out) > 5:
                        try:
                            _psf_e = pd.to_numeric(df_out["flux_psf_e"], errors="coerce")
                            _det_uid_col = pd.to_numeric(df_out["det_uid"], errors="coerce")
                            _ap_e_vals = np.array([
                                flux_init_map.get(int(u), np.nan) * GAIN
                                for u in _det_uid_col
                            ], dtype=float)
                            _ratio = _psf_e.to_numpy(float) / _ap_e_vals
                            _ratio_ok = _ratio[np.isfinite(_ratio) & (_ratio > 0)]
                            if len(_ratio_ok) >= 5:
                                med_ratio = float(np.median(_ratio_ok))
                                if not (0.3 < med_ratio < 3.0):
                                    self._log(
                                        f"  [WARN] flux unit mismatch? "
                                        f"median(psf_e/ap_e)={med_ratio:.3f} "
                                        f"(expected ~1.0). Check GAIN setting."
                                    )
                                else:
                                    self._log(
                                        f"  [UNIT] flux sanity OK: "
                                        f"median(psf_e/ap_e)={med_ratio:.3f} n={len(_ratio_ok)}"
                                    )
                        except Exception:
                            pass
                    # ─────────────────────────────────────────────────────────

                    out_tsv = output_dir / f"photometry_{fname}.tsv"
                    df_out.to_csv(out_tsv, sep="\t", index=False, encoding="utf-8-sig")
                    # Save step4 seed positions so the UI can tag iter2+ detections
                    # as "신규검출 (step4 미검출)" vs "재검출 (step4 기검출)".
                    seed_xy_path = output_dir / f"seed_xy_{fname}.npy"
                    np.save(str(seed_xy_path), init_xy_ui.astype(np.float32))

                    residual_meta = {
                        "file": fname,
                        "filter": this_filter,
                        # Which PSF model produced these magnitudes. Without it
                        # an ePSF product and a Moffat product are
                        # indistinguishable after the fact, and the two differ
                        # by a median 8 mmag on the same frame.
                        "psf_build_mode": str(psf_type_built),
                        "psf_fit_engine": str(psf_fit_engine_cfg),
                        "bkg_med": float(bkg_med),
                        "timing": {
                            "total_s": float(_t_total),
                            "background_s": float(_t_bkg),
                            "epsf_s": float(_t_epsf),
                            "substar_s": float(_t_sub),
                            "fit_s": float(_t_p1),
                            "second_fit_s": float(_t_p2),
                        },
                        "n_new_raw": int(n_new_raw_total),
                        "rawxy_iter2_path": f"rawxy_iter2_{fname}.npy",
                        "seedxy_path": seed_xy_path.name,
                        "fit_window": {
                            "mode": fit_window_plan.mode,
                            "shape_px": int(fit_shape_frame),
                            "energy_fraction": (
                                float(fit_window_plan.energy_fraction)
                                if np.isfinite(fit_window_plan.energy_fraction)
                                else None
                            ),
                            "target_energy_fraction": float(
                                fit_window_plan.target_energy_fraction
                            ),
                            "noise_equivalent_area_px": (
                                float(psf_nea_frame)
                                if np.isfinite(psf_nea_frame) else None
                            ),
                            "reason": fit_window_plan.reason,
                        },
                        "epsf_reference": {
                            "n_detected": int(n_epsf_detected),
                            "n_candidates": int(n_epsf_candidates),
                            "n_candidates_pre_morph": int(n_epsf_candidates_pre_morph),
                            "n_candidates_post_morph": int(n_epsf_candidates_post_morph),
                            "n_isolated": int(n_iso),
                            "n_selected": int(n_epsf_selected),
                            "n_morphology_relaxed_selected": int(
                                n_epsf_morphology_relaxed_selected
                            ),
                            "target": int(epsf_plan_target),
                            "grid_size": int(epsf_grid_size),
                            "contamination_aware": bool(epsf_contamination_filter),
                            "n_low_contamination": int(n_epsf_low_contamination),
                            "n_core_rejected": int(n_epsf_core_rejected),
                            "n_fallback_selected": int(n_epsf_fallback_selected),
                            "selected_median_contamination": (
                                float(epsf_selected_median_contamination)
                                if np.isfinite(epsf_selected_median_contamination)
                                else None
                            ),
                            "catalog_path": epsf_reference_catalog_name,
                            # EPSF 품질 검사 (로그 경고의 기계 판독 가능한 영속본)
                            "quality_n_blobs": int(epsf_quality_n_blobs),
                            "quality_double_peak": bool(epsf_quality_n_blobs > 1),
                            "quality_max_quadrant_frac": (
                                float(epsf_quality_max_quadrant_frac)
                                if np.isfinite(epsf_quality_max_quadrant_frac)
                                else None
                            ),
                            "quality_asymmetric": bool(
                                np.isfinite(epsf_quality_max_quadrant_frac)
                                and epsf_quality_max_quadrant_frac > 0.45
                            ),
                        },
                        "flux_scale": {
                            "enabled": bool(flux_scale_correction),
                            "applied": bool(scale_result.applied),
                            "scale": float(scale_result.scale),
                            "n_matched": int(scale_result.n_matched),
                            "n_candidates": int(scale_result.n_candidates),
                            "n_used": int(scale_result.n_used),
                            "median_delta_mag_raw": (
                                float(scale_result.median_delta_mag_raw)
                                if np.isfinite(scale_result.median_delta_mag_raw)
                                else None
                            ),
                            "scatter_mag": (
                                float(scale_result.scatter_mag)
                                if np.isfinite(scale_result.scatter_mag)
                                else None
                            ),
                            "reason": scale_result.reason,
                            "catalog_path": flux_scale_reference_name,
                        },
                        "iters": iter_records,
                        "core_cut": {
                            "enabled": bool(core_cut.enabled),
                            "center_x": float(core_cut.center_x) if np.isfinite(core_cut.center_x) else None,
                            "center_y": float(core_cut.center_y) if np.isfinite(core_cut.center_y) else None,
                            "radius_px": float(core_cut.radius_px) if np.isfinite(core_cut.radius_px) else None,
                            "method": core_cut.method,
                            "reason": core_cut.reason,
                            "n_excluded_init": int(n_core_excluded_init),
                            "n_excluded_redetect": int(n_core_excluded_redetect),
                            "n_excluded_result": int(n_core_excluded_result),
                        },
                    }
                    self.worker_status.emit(wid, fname, "Save", 95)
                    # Keep final products and metadata for UI reload/QA.
                    res_path = output_dir / f"residual_{fname}"
                    starsub_path = output_dir / f"starsub_{fname}"
                    hdr_res = fits.Header()
                    hdr_res["FILTER"] = this_filter
                    hdr_res["BKGMED"] = float(bkg_med)
                    hdr_res["FITWIN"] = int(fit_shape_frame)
                    if np.isfinite(psf_nea_frame):
                        hdr_res["PSFNEA"] = float(psf_nea_frame)
                    residual_out = np.asarray(residual, dtype=np.float32)
                    fits.writeto(str(res_path), residual_out, hdr_res, overwrite=True)
                    starsub_out = np.empty_like(residual_out)
                    np.add(
                        residual_out,
                        np.float32(bkg_med),
                        out=starsub_out,
                        casting="unsafe",
                    )
                    fits.writeto(str(starsub_path), starsub_out, hdr_res, overwrite=True)
                    del starsub_out
                    meta_path = output_dir / f"residual_meta_{fname}.json"
                    meta_path.write_text(json.dumps(residual_meta, ensure_ascii=False, indent=2), encoding="utf-8")
                    rawxy_iter2_path = output_dir / f"rawxy_iter2_{fname}.npy"
                    np.save(str(rawxy_iter2_path), np.asarray(raw_new_xy, dtype=np.float32))

                    merged_new_xy = None
                    if phot_result is not None and len(phot_result) > 0:
                        try:
                            if "iter_detected" in phot_result.colnames:
                                x_all = np.asarray(phot_result["x_fit"], dtype=float)
                                y_all = np.asarray(phot_result["y_fit"], dtype=float)
                                it_raw = np.asarray(phot_result["iter_detected"], dtype=float)
                                it_all = np.where(np.isfinite(it_raw), it_raw, 1.0).astype(int)
                                m_new = np.isfinite(x_all) & np.isfinite(y_all) & (it_all > 1)
                                if np.any(m_new):
                                    merged_new_xy = np.column_stack([x_all[m_new], y_all[m_new]])
                        except Exception:
                            merged_new_xy = None
                    n_rows = len(phot_rows)
                    if not df_out.empty:
                        _clean_output = (
                            df_out["mag_psf"].notna()
                            & (pd.to_numeric(df_out["flags_psf"], errors="coerce").fillna(-1) == 0)
                        )
                        n_good = int(np.sum(_clean_output))
                    else:
                        _clean_output = pd.Series(dtype=bool)
                        n_good = 0
                    if not df_out.empty:
                        _forced_output = df_out.get(
                            "forced_psf", pd.Series(False, index=df_out.index)
                        ).map(_as_bool)
                        _flux_output = pd.to_numeric(df_out["flux_psf_e"], errors="coerce")
                        _flags_output = pd.to_numeric(
                            df_out["flags_psf"], errors="coerce"
                        ).fillna(0).astype(np.int64)
                        n_forced = int(np.sum(_forced_output))
                        n_forced_negative = int(np.sum(
                            _forced_output & np.isfinite(_flux_output) & (_flux_output <= 0)
                        ))
                        n_crowding_unreliable = int(np.sum(
                            (_flags_output & int(PSFFitFlag.CROWDING_UNRELIABLE)) != 0
                        ))
                    else:
                        n_forced = 0
                        n_forced_negative = 0
                        n_crowding_unreliable = 0
                    median_qfit_noise_ratio = (
                        _median_value(df_out.loc[_clean_output], "qfit_noise_ratio")
                        if n_good else np.nan
                    )
                    frame_assessment = assess_psf_frame_quality(
                        n_sources=n_rows,
                        n_good=n_good,
                        n_crowding_unreliable=n_crowding_unreliable,
                        median_qfit_noise_ratio=median_qfit_noise_ratio,
                        epsf_n_selected=n_epsf_selected,
                        epsf_median_contamination=epsf_selected_median_contamination,
                        frame_fwhm_px=fwhm_med,
                        frame_fwhm_max_px=fwhm_qc_max_px,
                    )
                    idx_row = {
                        "file": fname,
                        "filter": this_filter,
                        "frame_fwhm_px": float(fwhm_med),
                        "frame_fwhm_arcsec": fwhm_arcsec,
                        "frame_fwhm_qc_max_px": fwhm_qc_max_px,
                        "frame_total_elapsed_s": float(_t_total),
                        "background_elapsed_s": float(_t_bkg),
                        "epsf_elapsed_s": float(_t_epsf),
                        "substar_elapsed_s": float(_t_sub),
                        "fit_elapsed_s": float(_t_p1),
                        "second_fit_elapsed_s": float(_t_p2),
                        "n": n_rows,
                        "n_goodmag": n_good,
                        "n_fail": n_rows - n_good,
                        "psf_clean_fraction": frame_assessment.clean_fraction,
                        "psf_fit_failure_fraction": frame_assessment.fit_failure_fraction,
                        "psf_crowding_unreliable_fraction": (
                            frame_assessment.crowding_unreliable_fraction
                        ),
                        "psf_qc_status": frame_assessment.status,
                        "psf_qc_score": frame_assessment.score,
                        "psf_qc_reasons": ",".join(frame_assessment.reasons),
                        "n_new_iter": int(n_new_total),
                        "n_forced": n_forced,
                        "n_forced_negative": n_forced_negative,
                        "n_crowding_unreliable": n_crowding_unreliable,
                        "median_qfit": _median_value(df_out.loc[_clean_output], "qfit") if n_good else np.nan,
                        "median_qfit_noise_ratio": median_qfit_noise_ratio,
                        "median_cfit": _median_value(df_out.loc[_clean_output], "cfit") if n_good else np.nan,
                        "median_reduced_chi2": (
                            _median_value(df_out.loc[_clean_output], "reduced_chi2") if n_good else np.nan
                        ),
                        "stop_reason": _engine_stop_reason if psf_fit_engine_cfg == "apex_iterative" else "",
                        "epsf_n_detected": int(n_epsf_detected),
                        "epsf_n_candidates": int(n_epsf_candidates),
                        "epsf_n_candidates_pre_morph": int(n_epsf_candidates_pre_morph),
                        "epsf_n_candidates_post_morph": int(n_epsf_candidates_post_morph),
                        "epsf_n_isolated": int(n_iso),
                        "epsf_n_selected": int(n_epsf_selected),
                        "epsf_n_morphology_relaxed_selected": int(
                            n_epsf_morphology_relaxed_selected
                        ),
                        "epsf_target": int(epsf_plan_target),
                        "epsf_grid_size": int(epsf_grid_size),
                        "epsf_contamination_aware": bool(epsf_contamination_filter),
                        "epsf_n_low_contamination": int(n_epsf_low_contamination),
                        "epsf_n_core_rejected": int(n_epsf_core_rejected),
                        "epsf_n_fallback_selected": int(n_epsf_fallback_selected),
                        "epsf_median_contamination": epsf_selected_median_contamination,
                        "fit_window_mode": fit_window_plan.mode,
                        "fit_window_px": int(fit_shape_frame),
                        "fit_window_energy": fit_window_plan.energy_fraction,
                        "fit_window_target_energy": fit_window_plan.target_energy_fraction,
                        "psf_nea_px": psf_nea_frame,
                        "psf_aperture_scale_enabled": bool(flux_scale_correction),
                        "psf_aperture_scale_applied": bool(scale_result.applied),
                        "psf_aperture_scale": float(scale_result.scale),
                        "psf_aperture_scale_n": int(scale_result.n_used),
                        "psf_aperture_scale_scatter_mag": float(scale_result.scatter_mag),
                        "psf_aperture_scale_reason": scale_result.reason,
                        "core_cut_enabled": bool(core_cut.enabled),
                        "core_cut_x_px": round(float(core_cut.center_x), 3) if np.isfinite(core_cut.center_x) else np.nan,
                        "core_cut_y_px": round(float(core_cut.center_y), 3) if np.isfinite(core_cut.center_y) else np.nan,
                        "core_cut_radius_px": round(float(core_cut.radius_px), 3) if np.isfinite(core_cut.radius_px) else np.nan,
                        "n_core_excluded_init": int(n_core_excluded_init),
                        "n_core_excluded_redetect": int(n_core_excluded_redetect),
                        "n_core_excluded_result": int(n_core_excluded_result),
                    }
                    idx_row.update(noise_info)
                    self.worker_status.emit(wid, fname, "Done", 100)
                    return {
                        "file": fname,
                        "status": "processed",
                        "idx_row": idx_row,
                        "epsf_key": (f"[{psf_type_built.upper()}] {this_filter} | {fname}"
                                    if epsf_emit_arr is not None else None),
                        "epsf_frame": fname if epsf_emit_arr is not None else None,
                        "epsf_arr": epsf_emit_arr,
                        "residual_meta": residual_meta,
                        "new_xy": merged_new_xy,
                    }
                except Exception as frame_e:
                    self.worker_status.emit(wid, fname, "Error", 100)
                    return {"file": fname, "status": "error", "reason": f"{frame_e}\n{traceback.format_exc()}"}

            ex = ThreadPoolExecutor(max_workers=max_workers)
            self._executor = ex
            future_map: dict = {}
            next_idx = 0

            def _submit_next():
                nonlocal next_idx
                if next_idx >= total:
                    return False
                fname_n = frames[next_idx]
                future_map[ex.submit(process_single_frame, fname_n)] = fname_n
                next_idx += 1
                return True

            for _ in range(min(max_workers, total)):
                _submit_next()

            try:
                while future_map:
                    # Stop mode: cancel queued (not-started) futures and do not submit new ones.
                    if self._stop_requested:
                        n_cancel = 0
                        for fut, fname_c in list(future_map.items()):
                            if fut.cancel():
                                del future_map[fut]
                                completed[0] += 1
                                counters["stopped"] += 1
                                n_cancel += 1
                                self.progress.emit(completed[0], total, fname_c)
                                self._log(f"[{completed[0]}/{total}] STOP {fname_c} | cancelled")
                        if n_cancel > 0:
                            self._log(f"Stop requested | cancelled pending={n_cancel}")
                        if not future_map:
                            break

                    done, _ = wait(tuple(future_map.keys()), timeout=0.2, return_when=FIRST_COMPLETED)
                    now = time.time()
                    n_done = int(completed[0])
                    n_running = int(len(future_map))
                    n_queued = int(max(0, total - next_idx))
                    progress_changed = (n_done != last_done_count)
                    if progress_changed:
                        last_done_count = n_done

                    if (now - last_hb) >= 8.0:
                        eta_txt = "--:--"
                        if n_done > 0:
                            elapsed = max(0.0, now - run_t0)
                            eta_txt = _fmt_eta((elapsed / n_done) * max(0, total - n_done))

                        if progress_changed:
                            self._log(
                                f"[PROG] done={n_done}/{total} | running={n_running} | queued={n_queued} | ETA~{eta_txt}"
                            )
                            self.progress.emit(n_done, total, f"RUN={n_running} QUEUE={n_queued} ETA~{eta_txt}")
                            last_hb = now
                            last_stall_log = now
                        else:
                            # Long fit phases can run for minutes; avoid spamming identical lines.
                            if (now - last_stall_log) >= 30.0:
                                active_names = list(future_map.values())
                                active_txt = ", ".join(active_names[:3]) if active_names else "-"
                                self._log(
                                    f"[PROG] waiting | done={n_done}/{total} | running={n_running} | "
                                    f"queued={n_queued} | active={active_txt} | ETA~{eta_txt}"
                                )
                                self.progress.emit(n_done, total, f"RUN={n_running} QUEUE={n_queued} ETA~{eta_txt}")
                                last_stall_log = now
                                last_hb = now
                    if not done:
                        continue

                    for fut in done:
                        fname = future_map.pop(fut, None)
                        if fname is None:
                            continue

                        if fut.cancelled():
                            completed[0] += 1
                            counters["stopped"] += 1
                            self.progress.emit(completed[0], total, fname)
                            self._log(f"[{completed[0]}/{total}] STOP {fname} | cancelled")
                            continue

                        try:
                            out = fut.result()
                        except Exception as e:
                            out = {"file": fname, "status": "error", "reason": str(e)}

                        completed[0] += 1
                        self.progress.emit(completed[0], total, out.get("file", fname))
                        status = out.get("status", "error")

                        if status == "processed":
                            idx_row = out.get("idx_row", {})
                            if idx_row:
                                index_rows.append(idx_row)
                                self.frame_done.emit(out["file"], idx_row)
                            if out.get("epsf_key") and out.get("epsf_arr") is not None:
                                self.epsf_ready.emit(out["epsf_key"], out.get("epsf_frame", out["file"]), out["epsf_arr"])
                            self.residual_ready.emit(out["file"], out.get("residual_meta"), out.get("new_xy"))
                            counters["processed"] += 1
                            self._log(
                                f"[{completed[0]}/{total}] OK {out['file']} | "
                                f"f={idx_row.get('filter', '?')} n={idx_row.get('n', 0)} "
                                f"good={idx_row.get('n_goodmag', 0)} new_iter={idx_row.get('n_new_iter', 0)}"
                            )
                        elif status == "no_detect":
                            counters["no_detect"] += 1
                            self._log(f"[{completed[0]}/{total}] SKIP {out['file']} | reason={out.get('reason', status)}")
                        elif status == "no_fits":
                            counters["no_fits"] += 1
                            self._log(f"[{completed[0]}/{total}] SKIP {out['file']} | reason={out.get('reason', status)}")
                        elif status == "stopped":
                            counters["stopped"] += 1
                            self._log(f"[{completed[0]}/{total}] STOP {out['file']} | reason={out.get('reason', status)}")
                        elif status == "fit_failed":
                            self._log(f"[{completed[0]}/{total}] FAIL {out['file']} | reason={out.get('reason', status)}")
                        elif status == "no_valid_init":
                            self._log(f"[{completed[0]}/{total}] SKIP {out['file']} | reason={out.get('reason', status)}")
                        else:
                            self._log(f"[{completed[0]}/{total}] ERROR {out['file']} | {out.get('reason', 'unknown')}")

                        # Keep pipeline fed only while not stopping.
                        if not self._stop_requested:
                            _submit_next()
            finally:
                remaining_unscheduled = max(0, total - next_idx)
                if self._stop_requested and remaining_unscheduled > 0:
                    counters["stopped"] += remaining_unscheduled
                    self._log(f"Stop requested | not_submitted={remaining_unscheduled}")
                try:
                    ex.shutdown(wait=False, cancel_futures=True)
                except Exception:
                    pass
                self._executor = None

            if index_rows:
                pd.DataFrame(index_rows).to_csv(output_dir / "photometry_index.csv", index=False)

            self._log(
                f"Done | processed={counters['processed']} | "
                f"no_detect={counters['no_detect']} | no_fits={counters['no_fits']} | "
                f"stopped={counters['stopped']}"
            )
            self.finished.emit({"frames": total, **counters})

        except Exception as e:
            self.error.emit("PSF_WORKER", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")
            self.finished.emit({})


# ── PSF Photometry Window ─────────────────────────────────────────────────────

class PSFPhotometryWindow(StepWindowBase):
    """Step 8 - PSF Photometry (skippable).

    If skipped, Step 9 Master ID Editor falls back to Step 7 forced photometry results.
    """

    def __init__(self, params, file_manager, project_state, main_window):
        self.file_manager = file_manager
        self.workflow_mode = str(getattr(main_window, "mode", "cmd")).lower()
        self.downstream_name = (
            "Target/Comparison Selection"
            if self.workflow_mode == "lc"
            else "Master ID Editor"
        )
        self.worker = None
        self.file_list = []
        self.use_cropped = False
        self.log_window = None
        self._skip_psf = False
        # In-memory cache (lost on window close → reloaded from disk in restore_state)
        self._last_epsf: dict[str, np.ndarray] = {}          # display_key → epsf array
        self._residual_meta: dict[str, dict] = {}            # fname  → residual metadata + iter records
        self._last_new_xy: dict[str, np.ndarray | None] = {} # fname  → new-detect XY or None
        self._cutout_idx: int = 0  # current star index in cutout viewer
        self._run_started_ts: float | None = None
        self._log_worker_frame: dict[int, str] = {}    # worker_id → current frame name
        self._current_psf_run_frames: list[str] = []
        self._current_psf_run_signature: dict | None = None
        self._psf_cache_validation_key: str | None = None
        self._psf_cache_validation_result: tuple[bool, str] = (False, "not checked")
        self._final_diag_data = pd.DataFrame()
        self._final_diag_summary: dict[str, object] = {}

        super().__init__(
            step_index=7,
            step_name="PSF Photometry",
            params=params,
            project_state=project_state,
            main_window=main_window,
        )
        self.setup_step_ui()
        self.restore_state()

    def setup_step_ui(self):
        info = QLabel(
            "Optional PSF photometry using photutils EPSFBuilder.\n"
            f"Click Skip PSF to continue to {self.downstream_name}; downstream "
            "steps will use Step 7 forced aperture photometry."
        )
        info.setWordWrap(True)
        info.setStyleSheet("QLabel { background-color: #E8F5E9; padding: 8px; margin-bottom: 6px; }")
        self.content_layout.addWidget(info)

        # ── Controls ──────────────────────────────────────────────────────────
        ctrl = QHBoxLayout()
        self.btn_params = create_parameter_button("PSF Parameters")
        self.btn_params.clicked.connect(self.open_parameters_dialog)
        ctrl.addWidget(self.btn_params)

        self.btn_skip = QPushButton("Skip PSF →")
        self.btn_skip.setStyleSheet(
            "QPushButton { background-color: #FF7043; color: white; font-weight: bold; padding: 8px 20px; }"
        )
        self.btn_skip.setToolTip(
            f"Skip PSF photometry; {self.downstream_name} will use Step 7 forced aperture results."
        )
        self.btn_skip.clicked.connect(self.skip_psf)
        ctrl.addWidget(self.btn_skip)

        self.chk_use_existing_output = create_output_reuse_checkbox(
            True,
            "Load Step 8 PSF outputs when photometry_index.csv, per-frame TSVs, "
            "ePSF/residual files, and the saved PSF signature all match the current run. "
            "Disable to force recomputation.",
        )
        ctrl.addWidget(self.chk_use_existing_output)

        ctrl.addStretch()

        self.btn_run = QPushButton("Run PSF")
        self.btn_run.setStyleSheet(
            "QPushButton { background-color: #388E3C; color: white; font-weight: bold; padding: 8px 24px; }"
        )
        self.btn_run.clicked.connect(self.run_psf)
        ctrl.addWidget(self.btn_run)

        self.btn_stop = QPushButton("Stop")
        self.btn_stop.setEnabled(False)
        self.btn_stop.clicked.connect(self.stop_psf)
        ctrl.addWidget(self.btn_stop)

        self.btn_log = QPushButton("Log")
        self.btn_log.clicked.connect(self.show_log_window)
        ctrl.addWidget(self.btn_log)

        self.content_layout.addLayout(ctrl)

        # ── Progress ──────────────────────────────────────────────────────────
        prog = QHBoxLayout()
        self.progress_bar = QProgressBar()
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        prog.addWidget(self.progress_bar, 1)
        self.progress_label = QLabel("Ready")
        self.progress_label.setMinimumWidth(420)
        prog.addWidget(self.progress_label)
        self.content_layout.addLayout(prog)


        # ── Skip status label ─────────────────────────────────────────────────
        self.skip_label = QLabel("")
        self.skip_label.setStyleSheet("QLabel { color: #FF7043; font-weight: bold; padding: 4px; }")
        self.content_layout.addWidget(self.skip_label)

        # ── Diagnostic tabs ───────────────────────────────────────────────────
        self.main_tabs = QTabWidget()
        self.content_layout.addWidget(self.main_tabs, 1)

        # Tab 0: EPSF Model
        epsf_tab = QWidget()
        epsf_layout = QVBoxLayout(epsf_tab)
        epsf_top = QHBoxLayout()
        epsf_top.addWidget(QLabel("Model:"))
        self.epsf_filter_combo = QComboBox()
        self.epsf_filter_combo.currentTextChanged.connect(self._on_epsf_filter_changed)
        epsf_top.addWidget(self.epsf_filter_combo)
        epsf_top.addStretch()
        epsf_layout.addLayout(epsf_top)

        self.epsf_fig = Figure(figsize=(8, 4))
        self.epsf_canvas = FigureCanvas(self.epsf_fig)
        self.epsf_toolbar = NavigationToolbar(self.epsf_canvas, self)
        epsf_layout.addWidget(self.epsf_toolbar)
        epsf_layout.addWidget(tame_canvas(self.epsf_canvas), 1)
        self.main_tabs.addTab(epsf_tab, "PSF Model")

        # Tab 1: Cutout viewer – Raw | Star-subtracted per star, per iter
        res_tab = QWidget()
        res_layout = QVBoxLayout(res_tab)

        # Row 1: frame / iter / mode selectors
        res_top = QHBoxLayout()
        res_top.addWidget(QLabel("Frame:"))
        self.res_file_combo = QComboBox()
        self.res_file_combo.currentTextChanged.connect(self._on_residual_frame_changed)
        res_top.addWidget(self.res_file_combo)
        res_top.addWidget(QLabel("Iter:"))
        self.res_iter_combo = QComboBox()
        self.res_iter_combo.currentTextChanged.connect(self._on_residual_iter_changed)
        res_top.addWidget(self.res_iter_combo)
        res_top.addStretch()
        res_layout.addLayout(res_top)

        # Row 2: star navigation (◀ / label / ▶)
        res_nav = QHBoxLayout()
        self.res_prev_btn = QPushButton("◀")
        self.res_prev_btn.setFixedWidth(36)
        self.res_prev_btn.clicked.connect(self._on_cutout_prev)
        res_nav.addWidget(self.res_prev_btn)
        self.res_star_label = QLabel("—")
        self.res_star_label.setMinimumWidth(70)
        self.res_star_label.setAlignment(Qt.AlignCenter)
        res_nav.addWidget(self.res_star_label)
        self.res_next_btn = QPushButton("▶")
        self.res_next_btn.setFixedWidth(36)
        self.res_next_btn.clicked.connect(self._on_cutout_next)
        res_nav.addWidget(self.res_next_btn)
        res_nav.addStretch()
        self.res_info_label = QLabel("")
        res_nav.addWidget(self.res_info_label)
        res_layout.addLayout(res_nav)

        self.res_fig = Figure(figsize=(8, 4))
        self.res_canvas = FigureCanvas(self.res_fig)
        self.res_toolbar = NavigationToolbar(self.res_canvas, self)
        res_layout.addWidget(self.res_toolbar)
        res_layout.addWidget(tame_canvas(self.res_canvas), 1)
        self.main_tabs.addTab(res_tab, "Residuals")

        # Tab 2: Photometry Table
        phot_tab = QWidget()
        phot_layout = QVBoxLayout(phot_tab)
        self.frame_table = QTableWidget()
        self.frame_table.setColumnCount(9)
        self.frame_table.setHorizontalHeaderLabels(
            ["Frame", "Filter", "N_psf", "N_goodmag", "N_fail", "N_new_iter",
             "Forced %", "PSF QC", "Time"]
        )
        # Forced % 는 PSF 결과를 읽을 때 반드시 함께 봐야 하는 값이다. 강제 측광
        # 위치(그 프레임에서 검출되지 않은 어두운 별)는 구경이 거의 0 을 재므로
        # PSF/구경 비교가 그쪽에서 폭주한다. 이 값이 높은 프레임의 PSF vs 구경
        # 통계는 검출된 별로만 걸러서 봐야 한다.
        _ft_tips = [
            "FITS 파일명", "필터명",
            "PSF 적합 대상 소스 수",
            "유효 등급을 얻은 수",
            "적합 실패 수",
            "잔차 재검출로 추가된 수",
            "강제 측광 비율 = n_forced / N_psf.\n"
            "그 프레임에서 검출되지 않아 마스터 위치로 강제 측광한 별의 비율.\n"
            "노출이 짧거나 청색 필터일수록 높다(실측: 같은 밤 B 75% vs R 50%).\n"
            "높으면 PSF/구경 플럭스 비교는 검출된 별로만 해야 한다.",
            "PSF QC 판정", "프레임 처리 시간",
        ]
        for _c, _tip in enumerate(_ft_tips):
            _it = self.frame_table.horizontalHeaderItem(_c)
            if _it is not None:
                _it.setToolTip(_tip)
        self.frame_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        for c in range(1, 9):
            self.frame_table.horizontalHeader().setSectionResizeMode(c, QHeaderView.ResizeToContents)
        self.frame_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        phot_layout.addWidget(self.frame_table)
        self.main_tabs.addTab(phot_tab, "Photometry")

        # Tab 3: QC Report (PSF statistics + Ap vs PSF comparison)
        qc_tab = QWidget()
        qc_outer = QVBoxLayout(qc_tab)

        qc_top_bar = QHBoxLayout()
        qc_refresh_btn = QPushButton("Refresh QC")
        qc_refresh_btn.clicked.connect(self._refresh_qc)
        qc_top_bar.addWidget(qc_refresh_btn)
        qc_top_bar.addStretch()
        qc_outer.addLayout(qc_top_bar)

        qc_splitter = QSplitter(Qt.Vertical)
        prevent_collapse(qc_splitter)
        qc_outer.addWidget(qc_splitter, 1)

        # Top: PSF statistics summary text
        self.qc_text = QTextEdit()
        self.qc_text.setReadOnly(True)
        self.qc_text.setStyleSheet("QTextEdit { font-family: monospace; font-size: 9pt; }")
        self.qc_text.setMaximumHeight(300)
        qc_splitter.addWidget(self.qc_text)

        # Bottom: Ap vs PSF matplotlib plot
        cmp_widget = QWidget()
        cmp_layout = QVBoxLayout(cmp_widget)
        cmp_top = QHBoxLayout()
        cmp_refresh_btn = QPushButton("Refresh Plot")
        cmp_refresh_btn.clicked.connect(self._plot_mag_comparison)
        cmp_top.addWidget(cmp_refresh_btn)
        cmp_top.addWidget(QLabel("Filter:"))
        self.cmp_filter_combo = QComboBox()
        self.cmp_filter_combo.addItem("all")
        self.cmp_filter_combo.currentTextChanged.connect(self._plot_mag_comparison)
        cmp_top.addWidget(self.cmp_filter_combo)
        cmp_top.addWidget(QLabel("Frame:"))
        self.cmp_frame_combo = QComboBox()
        self.cmp_frame_combo.addItem("all")
        self.cmp_frame_combo.currentTextChanged.connect(self._plot_mag_comparison)
        cmp_top.addWidget(self.cmp_frame_combo)
        self.cmp_flags0_only = QCheckBox("flags=0 only")
        self.cmp_flags0_only.setChecked(False)
        self.cmp_flags0_only.toggled.connect(self._plot_mag_comparison)
        cmp_top.addWidget(self.cmp_flags0_only)
        cmp_top.addWidget(QLabel("SNR ≥"))
        self.cmp_snr_min = QDoubleSpinBox()
        self.cmp_snr_min.setRange(0.0, 200.0)
        self.cmp_snr_min.setSingleStep(1.0)
        self.cmp_snr_min.setValue(0.0)
        self.cmp_snr_min.setDecimals(1)
        self.cmp_snr_min.setToolTip("0 = off")
        self.cmp_snr_min.valueChanged.connect(self._plot_mag_comparison)
        cmp_top.addWidget(self.cmp_snr_min)
        cmp_top.addWidget(QLabel("qfit/noise ≤"))
        self.cmp_qfit_max = QDoubleSpinBox()
        self.cmp_qfit_max.setRange(0.0, 10.0)
        self.cmp_qfit_max.setSingleStep(0.05)
        self.cmp_qfit_max.setValue(0.0)
        self.cmp_qfit_max.setDecimals(3)
        self.cmp_qfit_max.setToolTip("0 = off")
        self.cmp_qfit_max.valueChanged.connect(self._plot_mag_comparison)
        cmp_top.addWidget(self.cmp_qfit_max)
        cmp_top.addWidget(QLabel("|Δmag| ≤"))
        self.cmp_dmag_clip = QDoubleSpinBox()
        self.cmp_dmag_clip.setRange(0.0, 5.0)
        self.cmp_dmag_clip.setSingleStep(0.05)
        self.cmp_dmag_clip.setValue(0.0)
        self.cmp_dmag_clip.setDecimals(3)
        self.cmp_dmag_clip.setToolTip("0 = off")
        self.cmp_dmag_clip.valueChanged.connect(self._plot_mag_comparison)
        cmp_top.addWidget(self.cmp_dmag_clip)
        self.cmp_stats_label = QLabel("")
        self.cmp_stats_label.setWordWrap(True)
        cmp_top.addWidget(self.cmp_stats_label, 1)
        cmp_layout.addLayout(cmp_top)
        self.cmp_fig = Figure(figsize=(10, 4))
        self.cmp_canvas = FigureCanvas(self.cmp_fig)
        self.cmp_toolbar = NavigationToolbar(self.cmp_canvas, self)
        cmp_layout.addWidget(self.cmp_toolbar)
        cmp_layout.addWidget(tame_canvas(self.cmp_canvas), 1)
        qc_splitter.addWidget(cmp_widget)

        # Tallest page (478 px): scroll it so it does not set the window's
        # minimum height — see layout_rules.scroll_wrap.
        self.main_tabs.addTab(scroll_wrap(qc_tab), "QC")

        # Tab 4: per-frame final PSF diagnostics.
        final_tab = QWidget()
        final_layout = QVBoxLayout(final_tab)
        final_top = QHBoxLayout()
        final_top.addWidget(QLabel("Frame:"))
        self.final_diag_frame_combo = QComboBox()
        self.final_diag_frame_combo.setMinimumWidth(280)
        self.final_diag_frame_combo.currentTextChanged.connect(self._plot_final_diagnostics)
        final_top.addWidget(self.final_diag_frame_combo)
        final_refresh_btn = QPushButton("Refresh Diagnostics")
        final_refresh_btn.clicked.connect(self._refresh_final_diagnostics)
        final_top.addWidget(final_refresh_btn)
        self.final_diag_status = QLabel("Run Step 8 to generate diagnostics.")
        self.final_diag_status.setWordWrap(True)
        final_top.addWidget(self.final_diag_status, 1)
        final_layout.addLayout(final_top)

        self.final_diag_fig = Figure(figsize=(12, 7.2))
        self.final_diag_canvas = FigureCanvas(self.final_diag_fig)
        self.final_diag_toolbar = NavigationToolbar(self.final_diag_canvas, self)
        final_layout.addWidget(self.final_diag_toolbar)
        final_layout.addWidget(tame_canvas(self.final_diag_canvas), 1)
        self.main_tabs.addTab(final_tab, "Final Diagnostics")

        self.main_tabs.setCurrentIndex(0)

        # ── Log window ────────────────────────────────────────────────────────
        _log_worker_group = QGroupBox("Workers")
        _log_worker_group_layout = QVBoxLayout(_log_worker_group)
        _log_worker_group_layout.setContentsMargins(5, 5, 5, 5)
        self._worker_panel = WorkerStatusPanel(_log_worker_group)
        _log_worker_group_layout.addWidget(self._worker_panel)

        self.log_window = WorkflowLogWindow(
            self, "PSF Photometry Log & Workers",
            width=900, height=500,
            side_widget=_log_worker_group,
        )
        self.log_text = self.log_window.log_text

        # Keyboard shortcuts: ← → navigate cutout stars
        from PyQt5.QtWidgets import QShortcut
        from PyQt5.QtGui import QKeySequence
        QShortcut(QKeySequence(Qt.Key_Left),  self).activated.connect(self._on_cutout_prev)
        QShortcut(QKeySequence(Qt.Key_Right), self).activated.connect(self._on_cutout_next)

        self.populate_file_list()
        self.update_frame_table()
        self._update_skip_label()
        self._refresh_final_diagnostics()

    # ── File list ─────────────────────────────────────────────────────────────

    def populate_file_list(self):
        crop_active = crop_is_active(self.params.P.result_dir)
        cropped_dir = step2_cropped_dir(self.params.P.result_dir)
        if crop_active and cropped_dir.exists() and list(cropped_dir.glob("*.fit*")):
            files = sorted([f.name for f in cropped_dir.glob("*.fit*")])
            self.use_cropped = True
        else:
            if not self.file_manager.filenames:
                try:
                    self.file_manager.scan_files()
                except Exception:
                    pass
            files = self.file_manager.filenames
            self.use_cropped = False
        files = list(files)

        # Hard gate: downstream should skip frames where apcorr was not applied.
        apcorr_sum = step7_forced_phot_dir(self.params.P.result_dir) / "apcorr_summary.csv"
        if apcorr_sum.exists():
            try:
                df_apc = pd.read_csv(apcorr_sum)
                if (not df_apc.empty) and {"file", "apply"} <= set(df_apc.columns):
                    ok_vals = df_apc["apply"].astype(str).str.strip().str.lower().isin(
                        {"true", "1", "yes", "y", "on"}
                    )
                    ok_set = set(df_apc.loc[ok_vals, "file"].astype(str).map(lambda s: Path(str(s)).name))
                    before_n = len(files)
                    files = [f for f in files if Path(str(f)).name in ok_set]
                    self.log(f"[APCORR] apply=True frame filter: {len(files)}/{before_n} kept")
            except Exception as e:
                self.log(f"[APCORR] frame filter skipped ({e})")

        self.file_list = list(files)

    def _cache_dir_path(self) -> Path:
        return Path(getattr(self.params.P, "cache_dir", self.params.P.result_dir))

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
            return [PSFPhotometryWindow._signature_value(v) for v in value]
        if isinstance(value, dict):
            return {
                str(k): PSFPhotometryWindow._signature_value(v)
                for k, v in sorted(value.items(), key=lambda item: str(item[0]))
            }
        try:
            if pd.isna(value):
                return None
        except Exception:
            pass
        return str(value)

    @staticmethod
    def _first_existing_path(candidates: list[Path]) -> Path | None:
        for p in candidates:
            try:
                if p.exists() and p.stat().st_size > 0:
                    return p
            except Exception:
                continue
        return None

    @staticmethod
    def _newest_existing_path(candidates: list[Path]) -> Path | None:
        found = []
        for p in candidates:
            try:
                if p.exists() and p.stat().st_size > 0:
                    found.append(p)
            except Exception:
                continue
        if not found:
            return None
        return max(found, key=lambda p: p.stat().st_mtime_ns)

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

    def _psf_frames_after_qc(self) -> tuple[list[str], dict]:
        use_qc = should_use_frame_quality_qc(
            Path(self.params.P.result_dir),
            self.params.P,
            "phot_use_qc_pass_only",
            default=False,
        )
        frames, qc_info = filter_files_by_qc(
            Path(self.params.P.result_dir),
            list(self.file_list),
            require_qc=use_qc,
        )
        return list(frames), dict(qc_info or {})

    def _log_step4_qc_selection(self, qc_info: dict):
        if not should_use_frame_quality_qc(
            Path(self.params.P.result_dir),
            self.params.P,
            "phot_use_qc_pass_only",
            default=False,
        ):
            return
        if qc_info.get("applied"):
            self.log(f"Step4 QC: {qc_info.get('kept', 0)}/{qc_info.get('total', 0)} frame(s) kept.")
        elif qc_info.get("path") is None:
            self.log("Step4 QC: frame_quality.csv not found; using all frames.")
        else:
            self.log(f"Step4 QC: frame_quality.csv ignored ({qc_info.get('reason', 'unknown')}); using all frames.")

    def _build_psf_output_signature(self, frames: list[str]) -> dict:
        return build_psf_output_signature(
            self.params,
            frames,
            use_cropped=self.use_cropped,
            cache_dir=self._cache_dir_path(),
        )

    def _stored_psf_signature(self) -> dict | None:
        sig_path = step8_psf_dir(self.params.P.result_dir) / _PSF_SIGNATURE_FILE
        if not sig_path.exists():
            return None
        try:
            data = json.loads(sig_path.read_text(encoding="utf-8"))
            return data if isinstance(data, dict) else None
        except Exception:
            return None

    def _write_psf_output_signature(self, signature: dict):
        write_psf_output_signature(self.params.P.result_dir, signature)

    def _psf_signature_matches(self, signature: dict) -> tuple[bool, str]:
        stored = self._stored_psf_signature()
        if not stored:
            return False, "missing signature"
        if stored.get("signature_version") != _PSF_SIGNATURE_VERSION:
            return False, "signature version mismatch"
        if stored.get("signature_hash") != signature.get("signature_hash"):
            return False, "signature hash mismatch"
        return True, "ok"

    @staticmethod
    def _path_from_meta(out_dir: Path, name) -> Path | None:
        if name is None:
            return None
        text = str(name).strip()
        if not text:
            return None
        return out_dir / text

    def _existing_psf_output_covers(self, frames: list[str], signature: dict) -> tuple[bool, str]:
        ok, reason = self._psf_signature_matches(signature)
        if not ok:
            return False, reason

        out_dir = step8_psf_dir(self.params.P.result_dir)
        idx_path = out_dir / "photometry_index.csv"
        if not idx_path.exists():
            return False, "missing photometry_index.csv"
        try:
            idx = pd.read_csv(idx_path)
        except Exception as exc:
            return False, f"cannot read photometry_index.csv: {exc}"
        required_idx_cols = {"file", "filter"}
        if not required_idx_cols <= set(idx.columns):
            return False, "photometry_index.csv missing file/filter columns"

        expected_frames = [Path(str(f)).name for f in frames]
        expected_set = set(expected_frames)
        idx_files = [Path(str(f)).name for f in idx["file"].astype(str).tolist()]
        if len(idx_files) != len(expected_frames) or set(idx_files) != expected_set:
            return False, "photometry_index.csv frame set mismatch"

        expected_tsv_names = {f"photometry_{fname}.tsv" for fname in expected_frames}
        actual_tsv_names = {p.name for p in out_dir.glob("photometry_*.tsv")}
        if actual_tsv_names != expected_tsv_names:
            return False, "per-frame photometry TSV set mismatch"

        expected_epsf_paths: set[Path] = set()
        shared_epsf = bool(getattr(self.params.P, "psf_shared_filter_epsf", False))
        required_tsv_cols = {"det_uid", "x_fit", "y_fit", "mag_psf", "flags_psf"}
        if bool(getattr(self.params.P, "psf_flux_scale_correction", False)):
            required_tsv_cols.update({
                "flux_psf_raw_e",
                "psf_aperture_scale",
                "psf_aperture_scale_applied",
            })

        for fname in expected_frames:
            rows = idx[idx["file"].astype(str).map(lambda s: Path(s).name) == fname]
            if rows.empty:
                return False, f"missing index row for {fname}"
            filt = str(rows.iloc[0].get("filter", "")).strip()
            if not filt:
                return False, f"missing filter for {fname}"

            tsv_path = out_dir / f"photometry_{fname}.tsv"
            if not tsv_path.exists() or tsv_path.stat().st_size <= 0:
                return False, f"missing TSV for {fname}"
            try:
                tsv_head = pd.read_csv(tsv_path, sep="\t", nrows=5)
            except Exception as exc:
                return False, f"cannot read TSV for {fname}: {exc}"
            if not required_tsv_cols <= set(tsv_head.columns):
                return False, f"TSV columns incomplete for {fname}"

            for product_name in (f"residual_{fname}", f"starsub_{fname}"):
                product_path = out_dir / product_name
                if not product_path.exists() or product_path.stat().st_size <= 0:
                    return False, f"missing {product_name}"

            meta_path = out_dir / f"residual_meta_{fname}.json"
            if not meta_path.exists() or meta_path.stat().st_size <= 0:
                return False, f"missing residual metadata for {fname}"
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
            except Exception as exc:
                return False, f"cannot read residual metadata for {fname}: {exc}"
            if not isinstance(meta, dict):
                return False, f"invalid residual metadata for {fname}"
            epsf_reference = meta.get("epsf_reference", {})
            if (
                isinstance(epsf_reference, dict)
                and bool(epsf_reference.get("contamination_aware", False))
            ):
                reference_path = self._path_from_meta(
                    out_dir,
                    epsf_reference.get("catalog_path"),
                )
                if (
                    reference_path is None
                    or not reference_path.exists()
                    or reference_path.stat().st_size <= 0
                ):
                    return False, f"missing ePSF reference catalogue for {fname}"
            flux_scale = meta.get("flux_scale", {})
            if bool(getattr(self.params.P, "psf_flux_scale_correction", False)):
                if not isinstance(flux_scale, dict) or not bool(flux_scale.get("enabled", False)):
                    return False, f"missing PSF aperture-scale metadata for {fname}"
                scale_catalog = self._path_from_meta(out_dir, flux_scale.get("catalog_path"))
                if scale_catalog is not None and (
                    not scale_catalog.exists() or scale_catalog.stat().st_size <= 0
                ):
                    return False, f"missing PSF aperture-scale catalogue for {fname}"
            iters = meta.get("iters", [])
            if not isinstance(iters, list) or len(iters) == 0:
                return False, f"missing iteration metadata for {fname}"
            for key in ("seedxy_path", "rawxy_iter2_path"):
                p = self._path_from_meta(out_dir, meta.get(key))
                if p is None or not p.exists() or p.stat().st_size <= 0:
                    return False, f"missing {key} for {fname}"
            for rec in iters:
                if not isinstance(rec, dict):
                    return False, f"invalid iteration record for {fname}"
                for key in ("fitxy_path", "modelxy_path", "detxy_path", "appliedxy_path", "boxxy_path"):
                    p = self._path_from_meta(out_dir, rec.get(key))
                    if p is None or not p.exists() or p.stat().st_size <= 0:
                        return False, f"missing {key} for {fname}"
                for key in ("residual_path", "starsub_path"):
                    p = self._path_from_meta(out_dir, rec.get(key))
                    if p is not None and (not p.exists() or p.stat().st_size <= 0):
                        return False, f"missing {key} for {fname}"

            if shared_epsf:
                expected_epsf_paths.add(out_dir / f"epsf_model_{filt}.fits")
            else:
                expected_epsf_paths.add(out_dir / f"epsf_model_{filt}_{Path(fname).stem}.fits")

        for epsf_path in expected_epsf_paths:
            if not epsf_path.exists() or epsf_path.stat().st_size <= 0:
                return False, f"missing {epsf_path.name}"

        return True, "ok"

    def _current_psf_cache_status(self) -> tuple[bool, str]:
        frames, _ = self._psf_frames_after_qc()
        if not frames:
            return False, "no current frames"
        signature = self._build_psf_output_signature(frames)
        key = str(signature.get("signature_hash", ""))
        if key and key == self._psf_cache_validation_key:
            return self._psf_cache_validation_result
        result = self._existing_psf_output_covers(frames, signature)
        self._psf_cache_validation_key = key
        self._psf_cache_validation_result = result
        return result

    def _clear_psf_outputs(self) -> int:
        out_dir = step8_psf_dir(self.params.P.result_dir)
        if not out_dir.exists():
            return 0
        patterns = [
            _PSF_SIGNATURE_FILE,
            "photometry_index.csv",
            "photometry_*.tsv",
            "epsf_reference_*.csv",
            "psf_flux_scale_reference_*.csv",
            "epsf_model_*.fits",
            "residual_*",
            "starsub_*",
            "fitxy_iter*.npy",
            "modelxy_iter*.npy",
            "appliedxy_iter*.npy",
            "detxy_iter*.npy",
            "boxxy_iter*.npy",
            "seed_xy_*.npy",
            "rawxy_iter*.npy",
        ]
        removed = 0
        seen: set[Path] = set()
        for pat in patterns:
            for p in out_dir.glob(pat):
                if p in seen or not p.is_file():
                    continue
                seen.add(p)
                try:
                    p.unlink()
                    removed += 1
                except Exception:
                    pass
        return removed

    # ── Actions ───────────────────────────────────────────────────────────────

    def skip_psf(self):
        self._skip_psf = True
        self.save_state()
        self._update_skip_label()
        self.update_navigation_buttons()
        self.log(
            f"PSF skipped; {self.downstream_name} will use Step 7 forced aperture results."
        )

    def run_psf(self):
        if not self.file_list:
            QMessageBox.warning(self, "Warning", "No files to process")
            return
        if self.worker and self.worker.isRunning():
            return
        if not (step7_forced_phot_dir(self.params.P.result_dir) / "photometry_index.csv").exists():
            QMessageBox.warning(
                self, "Prerequisite",
                "Step 7 Forced Aperture Photometry must be completed first."
            )
            return

        frames_for_run, qc_info = self._psf_frames_after_qc()
        if not frames_for_run:
            QMessageBox.warning(
                self,
                "No frames",
                "No frames remain after Step 4 QC / Step 7 apcorr filtering.",
            )
            return
        signature = self._build_psf_output_signature(frames_for_run)
        self._psf_cache_validation_key = None
        self._current_psf_run_frames = list(frames_for_run)
        self._current_psf_run_signature = signature

        self._skip_psf = False
        self.log_text.clear()
        self.frame_table.setRowCount(0)
        self._last_epsf.clear()
        self._residual_meta.clear()
        self._last_new_xy.clear()
        self.epsf_filter_combo.clear()
        self.res_file_combo.clear()
        self.res_iter_combo.clear()
        # Keep frame selector populated during run so UI is not visually blank
        # before first residual metadata arrives.
        self.res_file_combo.addItems(frames_for_run)
        if self.res_file_combo.count() > 0:
            self.res_file_combo.setCurrentIndex(0)
        self._log_step4_qc_selection(qc_info)

        if getattr(self, "chk_use_existing_output", None) is not None and self.chk_use_existing_output.isChecked():
            cache_ok, cache_reason = self._existing_psf_output_covers(frames_for_run, signature)
            if cache_ok:
                self._psf_cache_validation_key = str(signature.get("signature_hash", ""))
                self._psf_cache_validation_result = (True, "ok")
                self.log(
                    f"[PSF][CACHE] Existing Step 8 output is complete "
                    f"({len(frames_for_run)} frame(s)); loading from disk."
                )
                self.progress_bar.setMaximum(len(frames_for_run))
                self.progress_bar.setValue(len(frames_for_run))
                self.progress_label.setText(
                    f"Cached Step 8 PSF output loaded ({len(frames_for_run)} frame(s))"
                )
                self._load_from_disk()
                self.update_frame_table()
                self._refresh_qc()
                self.save_state()
                self.update_navigation_buttons()
                self._current_psf_run_frames = []
                self._current_psf_run_signature = None
                return
            self.log(f"[PSF][CACHE] Existing output not reusable: {cache_reason}")

        removed = self._clear_psf_outputs()
        if removed:
            self.log(f"[PSF][CACHE] Removed {removed} stale Step 8 output file(s) before recompute.")

        # Clear log window worker bars from previous run
        self._log_worker_frame.clear()
        if hasattr(self, "_worker_panel") and self._worker_panel is not None:
            self._worker_panel.clear()

        self.log(f"Start PSF photometry | {len(frames_for_run)} frames")
        self._run_started_ts = time.time()

        self.worker = Step6PSFWorker(
            frames_for_run, self.params,
            self.params.P.data_dir, self.params.P.result_dir,
            self.params.P.cache_dir, self.use_cropped,
        )
        self.worker.progress.connect(self.on_progress)
        self.worker.worker_status.connect(self.on_worker_status)
        self.worker.frame_done.connect(self.on_frame_done)
        self.worker.epsf_ready.connect(self.on_epsf_ready)
        self.worker.residual_ready.connect(self.on_residual_ready)
        self.worker.finished.connect(self.on_finished)
        self.worker.error.connect(self.on_error)
        self.worker.log.connect(self.log)

        self.btn_run.setEnabled(False)
        self.btn_stop.setEnabled(True)
        self.btn_skip.setEnabled(False)
        self.progress_bar.setMaximum(len(frames_for_run))
        self.progress_bar.setValue(0)
        self.progress_label.setText(f"0/{len(frames_for_run)} | ETA --:-- | Starting...")
        self.worker.start()
        self.show_log_window()

    def stop_psf(self):
        if self.worker and self.worker.isRunning():
            self.worker.stop()
            self.btn_stop.setEnabled(False)
            self.progress_label.setText("Stopping... (running frames will finish current fit)")

    # ── Callbacks ─────────────────────────────────────────────────────────────

    def on_worker_status(self, worker_id: int, frame: str, stage: str, pct: int):
        self._log_worker_frame[int(worker_id)] = frame
        if hasattr(self, "_worker_panel") and self._worker_panel is not None:
            self._worker_panel.update_worker(worker_id, frame, stage, pct)

    def on_progress(self, current, total, filename):
        self.progress_bar.setValue(current)
        eta_txt = "pending"
        if self._run_started_ts is not None and current > 0 and total > 0:
            elapsed = max(0.0, time.time() - float(self._run_started_ts))
            per_frame = elapsed / float(current)
            eta_sec = max(0.0, per_frame * float(max(0, total - current)))
            eta_txt = self._fmt_duration(eta_sec)
        self.progress_label.setText(f"{current}/{total} | ETA {eta_txt} | {filename}")

    def on_frame_done(self, filename, result):
        r = self.frame_table.rowCount()
        self.frame_table.insertRow(r)
        self.frame_table.setItem(r, 0, QTableWidgetItem(filename))
        self.frame_table.setItem(r, 1, QTableWidgetItem(str(result.get("filter", ""))))
        self.frame_table.setItem(r, 2, QTableWidgetItem(str(result.get("n", 0))))
        self.frame_table.setItem(r, 3, QTableWidgetItem(str(result.get("n_goodmag", 0))))
        self.frame_table.setItem(r, 4, QTableWidgetItem(str(result.get("n_fail", 0))))
        self.frame_table.setItem(r, 5, QTableWidgetItem(str(result.get("n_new_iter", 0))))
        qc_status = str(result.get("psf_qc_status", "") or "")
        qc_item = QTableWidgetItem(qc_status)
        qc_item.setToolTip(str(result.get("psf_qc_reasons", "") or ""))
        self.frame_table.setItem(r, 6, qc_item)
        elapsed = _safe_float(result.get("frame_total_elapsed_s", np.nan), np.nan)
        self.frame_table.setItem(
            r, 7, QTableWidgetItem(f"{elapsed:.1f} s" if np.isfinite(elapsed) else "")
        )
        try:
            has_good_phot = int(result.get("n_goodmag", 0) or 0) > 0
        except (TypeError, ValueError):
            has_good_phot = False
        if qc_status == "FAIL":
            row_background = status_row_background(False)
        elif qc_status == "REVIEW":
            row_background = status_row_background(True, warning=True)
        else:
            row_background = status_row_background(has_good_phot)
        set_table_row_background(self.frame_table, r, row_background)
        self.frame_table.scrollToBottom()
        # Mark log window worker bar as done
        for w_key, fname in self._log_worker_frame.items():
            if fname == filename:
                if hasattr(self, "_worker_panel") and self._worker_panel is not None:
                    self._worker_panel.update_worker(w_key, fname, "Done", 100)
                break

    def on_epsf_ready(self, display_key: str, _frame_name: str, epsf_arr: np.ndarray):
        self._last_epsf[display_key] = epsf_arr
        current = self.epsf_filter_combo.currentText()
        self.epsf_filter_combo.blockSignals(True)
        self.epsf_filter_combo.clear()
        self.epsf_filter_combo.addItems(sorted(self._last_epsf.keys()))
        if current in self._last_epsf:
            self.epsf_filter_combo.setCurrentText(current)
        else:
            try:
                self.epsf_filter_combo.setCurrentText(display_key)
            except Exception:
                self.epsf_filter_combo.setCurrentIndex(0)
        self.epsf_filter_combo.blockSignals(False)
        self._plot_epsf(display_key)

    def on_residual_ready(self, fname: str, residual_meta: dict, new_xy):
        if isinstance(residual_meta, dict):
            self._residual_meta[fname] = residual_meta
        self._last_new_xy[fname] = new_xy  # ndarray or None
        current = self.res_file_combo.currentText()
        self.res_file_combo.blockSignals(True)
        self.res_file_combo.clear()
        self.res_file_combo.addItems(sorted(self._residual_meta.keys()))
        if current in self._residual_meta:
            self.res_file_combo.setCurrentText(current)
        else:
            self.res_file_combo.setCurrentIndex(self.res_file_combo.count() - 1)
        self.res_file_combo.blockSignals(False)
        self._cutout_idx = 0
        self._refresh_residual_iter_combo(fname)
        self._plot_cutout(fname)

    def on_error(self, src, err):
        self.log(f"ERROR {src}: {err}")

    def on_finished(self, summary):
        self.btn_run.setEnabled(True)
        self.btn_stop.setEnabled(False)
        self.btn_skip.setEnabled(True)
        elapsed_txt = ""
        if self._run_started_ts is not None:
            elapsed_txt = f" | elapsed {self._fmt_duration(time.time() - float(self._run_started_ts))}"
        self.progress_label.setText(f"Done{elapsed_txt}")
        self._run_started_ts = None
        self.log(f"PSF done: {summary}")
        if isinstance(summary, dict) and self._current_psf_run_signature:
            expected = len(self._current_psf_run_frames)
            processed = _to_int(summary.get("processed", 0), 0)
            stopped = _to_int(summary.get("stopped", 0), 0)
            if expected > 0 and stopped == 0 and processed == expected:
                self._write_psf_output_signature(self._current_psf_run_signature)
                cache_ok, cache_reason = self._existing_psf_output_covers(
                    self._current_psf_run_frames,
                    self._current_psf_run_signature,
                )
                if cache_ok:
                    self.log("[PSF][CACHE] Output signature saved for future reuse.")
                else:
                    sig_path = step8_psf_dir(self.params.P.result_dir) / _PSF_SIGNATURE_FILE
                    try:
                        sig_path.unlink()
                    except Exception:
                        pass
                    self.log(f"[PSF][CACHE] Signature not saved: output incomplete ({cache_reason}).")
            else:
                self.log(
                    "[PSF][CACHE] Output reuse disabled for this run: "
                    f"processed={processed}/{expected}, stopped={stopped}."
                )
        self._cleanup_worker()
        self._update_skip_label()
        self.update_frame_table()  # refresh Photometry tab from disk
        self._refresh_qc()           # refresh QC tab (stats + Ap vs PSF plot)
        self.save_state()
        self._psf_cache_validation_key = None
        self.update_navigation_buttons()
        self._current_psf_run_frames = []
        self._current_psf_run_signature = None

    # ── EPSF plot ─────────────────────────────────────────────────────────────

    def _on_epsf_filter_changed(self, display_key: str):
        self._plot_epsf(display_key)

    @staticmethod
    def _fmt_duration(seconds: float) -> str:
        s = int(max(0, round(float(seconds))))
        h, rem = divmod(s, 3600)
        m, sec = divmod(rem, 60)
        if h > 0:
            return f"{h:d}:{m:02d}:{sec:02d}"
        return f"{m:02d}:{sec:02d}"

    def _plot_epsf(self, display_key: str):
        if display_key not in self._last_epsf:
            return
        epsf_arr = self._last_epsf[display_key]

        is_moffat = display_key.startswith("[MOFFAT]")
        is_epsf   = display_key.startswith("[EPSF]")
        psf_label = "Moffat PSF" if is_moffat else "ePSF"
        px_label  = "px (native)" if is_moffat else "px (oversampled)"

        self.epsf_fig.clf()
        ax2d  = self.epsf_fig.add_subplot(121)
        ax_rad = self.epsf_fig.add_subplot(122)

        vmax = np.nanpercentile(epsf_arr, 99)
        im = ax2d.imshow(epsf_arr, origin="lower", cmap="viridis",
                         norm=mcolors.PowerNorm(gamma=0.5, vmin=0, vmax=max(vmax, 1e-10)))
        self.epsf_fig.colorbar(im, ax=ax2d, fraction=0.046, pad=0.04)
        ax2d.set_title(f"{psf_label} — {display_key}", fontsize=9)
        ax2d.set_xlabel(px_label, fontsize=8)
        ax2d.set_ylabel(px_label, fontsize=8)

        cy, cx = np.array(epsf_arr.shape) / 2.0
        yy, xx = np.mgrid[0:epsf_arr.shape[0], 0:epsf_arr.shape[1]]
        rr = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
        r_flat, v_flat = rr.ravel(), epsf_arr.ravel()
        order = np.argsort(r_flat)
        ax_rad.plot(r_flat[order], v_flat[order], ".", markersize=1, alpha=0.3, color="#1565C0")
        ax_rad.set_xlabel(f"Radius ({px_label})", fontsize=8)
        ax_rad.set_ylabel("PSF value", fontsize=8)
        ax_rad.set_title("Radial profile", fontsize=9)
        ax_rad.grid(True, alpha=0.3)

        self.epsf_fig.tight_layout()
        self.epsf_canvas.draw_idle()

    # ── Residual plot ─────────────────────────────────────────────────────────

    def _on_residual_frame_changed(self, fname: str):
        self._cutout_idx = 0
        self._refresh_residual_iter_combo(fname)
        self._plot_cutout(fname)

    def _on_residual_iter_changed(self, _iter_label: str):
        self._cutout_idx = 0
        fname = self.res_file_combo.currentText().strip()
        if fname:
            self._plot_cutout(fname)

    def _on_cutout_prev(self):
        if self._cutout_idx > 0:
            self._cutout_idx -= 1
            fname = self.res_file_combo.currentText().strip()
            if fname:
                self._plot_cutout(fname)

    def _on_cutout_next(self):
        self._cutout_idx += 1
        fname = self.res_file_combo.currentText().strip()
        if fname:
            self._plot_cutout(fname)

    def _get_iter_records(self, fname: str) -> list[dict]:
        meta = self._residual_meta.get(fname, {})
        recs = meta.get("iters", []) if isinstance(meta, dict) else []
        if not isinstance(recs, list):
            return []
        out = [r for r in recs if isinstance(r, dict)]
        out.sort(key=lambda d: int(d.get("iter", 0)))
        return out

    def _get_iter_record_by_no(self, fname: str, iter_no: int) -> dict | None:
        for rec in self._get_iter_records(fname):
            if int(rec.get("iter", -1)) == int(iter_no):
                return rec
        return None

    def _refresh_residual_iter_combo(self, fname: str):
        recs = self._get_iter_records(fname)
        self.res_iter_combo.blockSignals(True)
        self.res_iter_combo.clear()
        for rec in recs:
            i = int(rec.get("iter", 0))
            phase = str(rec.get("phase", "residual_fit"))
            label = f"{i} final flux" if phase == "final_flux" else f"{i} residual"
            self.res_iter_combo.addItem(label, i)
        if self.res_iter_combo.count() > 0:
            self.res_iter_combo.setCurrentIndex(self.res_iter_combo.count() - 1)
        self.res_iter_combo.blockSignals(False)

    def _load_starsub_for_iter(self, fname: str, rec: dict) -> np.ndarray | None:
        out_dir = step8_psf_dir(self.params.P.result_dir)
        meta = self._residual_meta.get(fname, {})
        bkg_med = float(meta.get("bkg_med", 0.0)) if isinstance(meta, dict) else 0.0

        starsub_name = str(rec.get("starsub_path", "")).strip()
        if starsub_name:
            p = out_dir / starsub_name
            if p.exists():
                try:
                    return fits.getdata(str(p)).astype(float)
                except Exception:
                    pass

        residual_name = str(rec.get("residual_path", "")).strip()
        if residual_name:
            p = out_dir / residual_name
            if p.exists():
                try:
                    res = fits.getdata(str(p)).astype(float)
                    return res + bkg_med
                except Exception:
                    pass

        return None

    def _load_snapshot_image(self, rec: dict, key: str) -> np.ndarray | None:
        out_dir = step8_psf_dir(self.params.P.result_dir)
        image_name = str(rec.get(key, "")).strip()
        if not image_name:
            return None
        path = out_dir / image_name
        if not path.exists():
            return None
        try:
            return fits.getdata(str(path)).astype(float)
        except Exception:
            return None

    def _load_xy_npy_for_iter(self, rec: dict, key: str, max_points: int = 500) -> np.ndarray:
        out_dir = step8_psf_dir(self.params.P.result_dir)
        arr_name = str(rec.get(key, "")).strip()
        if not arr_name:
            return np.zeros((0, 2), dtype=float)
        p = out_dir / arr_name
        if not p.exists():
            return np.zeros((0, 2), dtype=float)
        try:
            arr = np.load(str(p), allow_pickle=False)
            arr = np.asarray(arr, dtype=float)
            if arr.ndim != 2 or arr.shape[1] < 2:
                return np.zeros((0, 2), dtype=float)
            arr = arr[:, :2]
            finite = np.isfinite(arr[:, 0]) & np.isfinite(arr[:, 1])
            arr = arr[finite]
            if int(max_points) > 0:
                return arr[:max(0, int(max_points))]
            return arr
        except Exception:
            return np.zeros((0, 2), dtype=float)

    def _load_boxxy_for_iter(self, rec: dict, max_boxes: int = 500) -> np.ndarray:
        # Preferred: delta boxes (applied-from-previous + detected-this-iter).
        arr = self._load_xy_npy_for_iter(rec, "boxxy_path", max_points=max_boxes)
        if len(arr):
            return arr
        # Fallback 1: compose from separate arrays if present.
        arr_applied = self._load_xy_npy_for_iter(rec, "appliedxy_path", max_points=0)
        arr_detected = self._load_xy_npy_for_iter(rec, "detxy_path", max_points=0)
        if len(arr_applied) and len(arr_detected):
            arr = np.vstack([arr_applied, arr_detected])
        elif len(arr_applied):
            arr = arr_applied
        elif len(arr_detected):
            arr = arr_detected
        else:
            # Fallback 2: all fitted stars in this iteration.
            arr = self._load_xy_npy_for_iter(rec, "fitxy_path", max_points=0)
        if int(max_boxes) > 0:
            return arr[:max(0, int(max_boxes))]
        return arr

    def _load_modelxy_for_iter(self, rec: dict, max_points: int = 500) -> np.ndarray:
        arr = self._load_xy_npy_for_iter(rec, "modelxy_path", max_points=max_points)
        if len(arr):
            return arr
        # Backward compatibility with runs before modelxy_path existed:
        # iter>=2 should map to "new in this iter" rather than cumulative fit list.
        iter_no = int(rec.get("iter", 1))
        if iter_no > 1:
            arr = self._load_xy_npy_for_iter(rec, "detxy_path", max_points=max_points)
            if len(arr):
                return arr
        return self._load_xy_npy_for_iter(rec, "fitxy_path", max_points=max_points)

    def _resolve_fits_path_window(self, fname: str) -> Path | None:
        if self.use_cropped and crop_is_active(self.params.P.result_dir):
            cdir = step2_cropped_dir(self.params.P.result_dir)
            cpath = cdir / fname
            if cpath.exists():
                return cpath
        fpath = Path(self.params.P.data_dir) / fname
        return fpath if fpath.exists() else None

    # ── Cutout viewer ─────────────────────────────────────────────────────────

    def _plot_cutout(self, fname: str):  # noqa: C901
        """Show Raw | Star-subtracted cutouts for the selected star."""
        if fname not in self._residual_meta:
            self.res_fig.clf()
            ax = self.res_fig.add_subplot(111)
            ax.text(
                0.5, 0.5,
                "No residual result yet for this frame.\n(Still processing or frame skipped/failed)",
                transform=ax.transAxes,
                ha="center", va="center", fontsize=10, color="gray",
            )
            ax.set_title(fname, fontsize=9)
            self.res_star_label.setText("0/0")
            self.res_info_label.setText("waiting for residual_meta...")
            self.res_canvas.draw_idle()
            return
        recs = self._get_iter_records(fname)
        if not recs:
            return

        iter_no = self.res_iter_combo.currentData()
        selected = recs[-1]
        if iter_no is not None:
            for rec in recs:
                if int(rec.get("iter", -1)) == int(iter_no):
                    selected = rec
                    break

        # Fixed semantics requested by user:
        # - iter1: detected/fitted stars in iter1
        # - iter>=2: stars detected from residual(iter-1)
        iter_val = int(selected.get("iter", 0))
        phase = str(selected.get("phase", "residual_fit"))
        det_xy = self._load_xy_npy_for_iter(selected, "detxy_path", max_points=0)
        model_xy = self._load_modelxy_for_iter(selected, max_points=0)
        if phase == "final_flux":
            xy_list = self._load_xy_npy_for_iter(selected, "fitxy_path", max_points=0)
            mode_label = "fixed-position final flux"
        elif iter_val <= 1:
            xy_list = model_xy
            mode_label = "iter1 fitted stars"
        else:
            xy_list = det_xy if len(det_xy) > 0 else model_xy
            mode_label = f"iter{iter_val} detected-from-residual"

        res_std = float(selected.get("residual_std", np.nan))
        n_new_raw = int(selected.get("n_new_raw", 0))
        n_new_kept = int(selected.get("n_new_kept", 0))
        n_candidate_raw = int(selected.get("n_candidates_raw", 0))
        n_candidate_accepted = int(selected.get("n_candidates_accepted", 0))
        median_qfit = float(selected.get("median_qfit", np.nan))
        median_redchi = float(selected.get("median_reduced_chi2", np.nan))
        stop_reason = str(selected.get("stop_reason", ""))

        # Cutout half-size: driven by epsf_size_px so the PSF footprint is visible
        epsf_sz = int(selected.get("epsf_size_px", 25))
        half = max(epsf_sz // 2 + 4, 10)

        # Load raw FITS
        raw_img = None
        try:
            p = self._resolve_fits_path_window(fname)
            if p is not None:
                raw_img = fits.getdata(str(p)).astype(float)
        except Exception:
            pass

        model_img_snapshot = self._load_snapshot_image(selected, "model_path")
        residual_img_snapshot = self._load_snapshot_image(selected, "residual_path")

        full_cut_sz = 2 * half + 1

        def _cut_at(img, x_val: float, y_val: float):
            if img is None:
                return None
            nr, nc = img.shape
            cx = int(round(float(x_val)))
            cy = int(round(float(y_val)))
            x0 = max(0, cx - half)
            x1 = min(nc, cx + half + 1)
            y0 = max(0, cy - half)
            y1 = min(nr, cy + half + 1)
            cut = img[y0:y1, x0:x1]
            if cut.size == 0:
                return None
            if cut.shape == (full_cut_sz, full_cut_sz):
                return cut
            # Edge source: pad with NaN so the star stays centred in the panel.
            # NaN renders as the colormap bad-colour (neutral), making the
            # padding region visually distinct from real background.
            padded = np.full((full_cut_sz, full_cut_sz), np.nan, dtype=np.float64)
            dst_y = max(0, half - cy)
            dst_x = max(0, half - cx)
            padded[dst_y:dst_y + (y1 - y0), dst_x:dst_x + (x1 - x0)] = cut
            return padded

        # ── Filter xy_list to photometry-successful sources only ─────────────
        # Load the Step 8 PSF photometry TSV and keep only positions where
        # mag_psf is finite and FLAG_SAT is not set.  Saturated / edge /
        # fit-fail sources have NaN mag_psf or FLAG_SAT=1 — showing their
        # cutouts is misleading (over-subtraction rings, clipped PSF, etc.).
        if len(xy_list) > 0:
            try:
                _psf_tsv = step8_psf_dir(self.params.P.result_dir) / f"photometry_{fname}.tsv"
                if _psf_tsv.exists():
                    _df_phot = pd.read_csv(_psf_tsv, sep="\t")
                    if "saturated_psf" in _df_phot.columns:
                        _not_saturated = ~_df_phot["saturated_psf"].astype(str).str.lower().isin(
                            {"true", "1", "yes"}
                        )
                    else:
                        # Compatibility with catalogs written before standard fit flags.
                        _not_saturated = (
                            pd.to_numeric(
                                _df_phot.get("flags_psf", pd.Series(0, index=_df_phot.index)),
                                errors="coerce",
                            ).fillna(0).astype(int) & 1
                        ) == 0
                    _good = (
                        pd.to_numeric(_df_phot.get("mag_psf", pd.Series(dtype=float)), errors="coerce").notna() &
                        _not_saturated
                    )
                    _good_xy = _df_phot.loc[_good, ["x_fit", "y_fit"]].to_numpy(dtype=float)
                    if len(_good_xy) > 0:
                        _tree_good = cKDTree(_good_xy)
                        _d, _ = _tree_good.query(xy_list, k=1, workers=1)
                        _match_r = 1.5  # px
                        _keep = np.asarray(_d, dtype=float) <= _match_r
                        if np.any(_keep):
                            xy_list = xy_list[_keep]
            except Exception:
                pass  # on any error fall through and show all sources

        # Move edge sources to the end so idx=0 always shows a well-centred star.
        if raw_img is not None and len(xy_list) > 1:
            _nr, _nc = raw_img.shape
            _not_edge = (
                (xy_list[:, 0] >= half) & (xy_list[:, 0] < _nc - half) &
                (xy_list[:, 1] >= half) & (xy_list[:, 1] < _nr - half)
            )
            if np.any(_not_edge) and not np.all(_not_edge):
                _order = np.concatenate([np.where(_not_edge)[0], np.where(~_not_edge)[0]])
                xy_list = xy_list[_order]

        n = int(len(xy_list))
        idx = max(0, min(self._cutout_idx, n - 1)) if n > 0 else 0
        self._cutout_idx = idx
        self.res_star_label.setText(f"{idx + 1}/{n}" if n > 0 else "0/0")

        self.res_fig.clf()
        if n == 0:
            self.res_info_label.setText(
                f"iter {iter_val} | stars={n} | "
                f"new(raw/used)={n_new_raw}/{n_new_kept} | res_std={res_std:.3f}"
            )
            ax = self.res_fig.add_subplot(111)
            ax.text(
                0.5, 0.5,
                f"No sources for {mode_label}.",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=10, color="gray",
            )
            ax.set_title(f"{fname}  iter {iter_val}", fontsize=9)
            self.res_canvas.draw_idle()
            return

        x_c, y_c = float(xy_list[idx, 0]), float(xy_list[idx, 1])

        # ── Classify iter2+ detection: "신규검출" vs "재검출" ──────────────────
        # Compare the star's fitted position against step4 seed positions.
        # Within 2px → the seed existed in step4 (iter1 fit it but subtraction was poor).
        # Beyond 2px → genuinely new source not in step4.
        seed_tag = ""
        if iter_val >= 2:
            out_dir_ui = step8_psf_dir(self.params.P.result_dir)
            meta_ui = self._residual_meta.get(fname, {})
            seedxy_name = meta_ui.get("seedxy_path", "") if isinstance(meta_ui, dict) else ""
            if seedxy_name:
                seed_path_ui = out_dir_ui / seedxy_name
                if seed_path_ui.exists():
                    try:
                        seed_xy_ui = np.load(str(seed_path_ui)).astype(float)
                        if len(seed_xy_ui) > 0:
                            d_seed = np.hypot(seed_xy_ui[:, 0] - x_c, seed_xy_ui[:, 1] - y_c)
                            if np.min(d_seed) <= 2.0:
                                seed_tag = "재검출 (step4 기검출)"
                            else:
                                seed_tag = "신규검출 (step4 미검출)"
                    except Exception:
                        pass

        # Edge tag: shown when the star is within `half` pixels of the image boundary.
        edge_tag = ""
        if raw_img is not None:
            _nr, _nc = raw_img.shape
            if not (half <= x_c < _nc - half and half <= y_c < _nr - half):
                edge_tag = "경계소스"

        tags = "  " + "  ".join(f"[{t}]" for t in [seed_tag, edge_tag] if t) if (seed_tag or edge_tag) else ""
        self.res_info_label.setText(
            f"pass {iter_val} {phase} | stars={n} | new fit={n_new_kept} | "
            f"candidates={n_candidate_raw}/{n_candidate_accepted} | "
            f"res_std={res_std:.3f} qfit50={median_qfit:.3f} "
            f"chi2r50={median_redchi:.2f} | stop={stop_reason or '-'} | "
            f"xy=({x_c:.2f},{y_c:.2f}){tags}"
        )

        def _cut(img):
            return _cut_at(img, x_c, y_c)

        cut_raw = _cut(raw_img)
        cut_model = _cut(model_img_snapshot)
        cut_residual = _cut(residual_img_snapshot)
        panels: list[dict] = [
            {"img": cut_raw, "title": "Raw", "mark_detect": False, "residual": False},
            {"img": cut_model, "title": "PSF model", "mark_detect": False, "residual": False},
            {
                "img": cut_residual,
                "title": "Sky-sub residual",
                "mark_detect": phase != "final_flux" and iter_val > 1,
                "residual": True,
            },
        ]

        for i, p in enumerate(panels):
            cut = p.get("img", None)
            title = str(p.get("title", ""))
            mark_detect = bool(p.get("mark_detect", False))
            ax = self.res_fig.add_subplot(1, len(panels), i + 1)
            if cut is not None and cut.size > 0:
                if bool(p.get("residual", False)):
                    vmax = float(np.nanpercentile(np.abs(cut), 99))
                    panel_vmin, panel_vmax = -max(vmax, 1e-10), max(vmax, 1e-10)
                    panel_cmap = "coolwarm"
                else:
                    panel_vmin, panel_vmax = np.nanpercentile(cut, [1, 99])
                    panel_cmap = "gray"
                im = ax.imshow(
                    cut,
                    origin="lower",
                    cmap=panel_cmap,
                    vmin=panel_vmin,
                    vmax=panel_vmax,
                    interpolation="nearest",
                )
                self.res_fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                # Crosshair at star center
                cy_c = (cut.shape[0] - 1) / 2.0
                cx_c = (cut.shape[1] - 1) / 2.0
                ax.axhline(cy_c, color="#FF4444", lw=0.8, alpha=0.55, ls="--")
                ax.axvline(cx_c, color="#FF4444", lw=0.8, alpha=0.55, ls="--")
                if mark_detect:
                    ax.plot(
                        [cx_c], [cy_c],
                        marker="o", markersize=8,
                        markerfacecolor="none",
                        markeredgecolor="#FF5555",
                        markeredgewidth=1.2,
                    )
            else:
                ax.text(0.5, 0.5, "N/A", transform=ax.transAxes,
                        ha="center", va="center", fontsize=10, color="gray")
            ax.set_title(title, fontsize=8)
            ax.set_xlabel("Δx (px)", fontsize=7)
            ax.set_ylabel("Δy (px)", fontsize=7)

        suptitle_tags = "  " + "  ".join(f"[{_KO_TO_ASCII.get(t, t)}]" for t in [seed_tag, edge_tag] if t) if (seed_tag or edge_tag) else ""
        self.res_fig.suptitle(
            f"{fname}  |  {mode_label} #{idx + 1}/{n}{suptitle_tags}",
            fontsize=8, y=1.01,
        )
        self.res_fig.tight_layout()
        self.res_canvas.draw_idle()

    # ── Frame table refresh ───────────────────────────────────────────────────

    def update_frame_table(self):
        idx_path = step8_psf_dir(self.params.P.result_dir) / "photometry_index.csv"
        if not idx_path.exists() or not hasattr(self, "frame_table"):
            return
        try:
            idx = pd.read_csv(idx_path)
        except Exception:
            return
        self.frame_table.setRowCount(len(idx))
        for r, row in enumerate(idx.itertuples(index=False)):
            self.frame_table.setItem(r, 0, QTableWidgetItem(str(getattr(row, "file", ""))))
            self.frame_table.setItem(r, 1, QTableWidgetItem(str(getattr(row, "filter", ""))))
            self.frame_table.setItem(r, 2, QTableWidgetItem(str(int(_safe_float(getattr(row, "n", 0), 0)))))
            self.frame_table.setItem(r, 3, QTableWidgetItem(str(int(_safe_float(getattr(row, "n_goodmag", 0), 0)))))
            self.frame_table.setItem(r, 4, QTableWidgetItem(str(int(_safe_float(getattr(row, "n_fail", 0), 0)))))
            self.frame_table.setItem(r, 5, QTableWidgetItem(str(int(_safe_float(getattr(row, "n_new_iter", 0), 0)))))
            n_psf = _safe_float(getattr(row, "n", 0), 0)
            n_forced = _safe_float(getattr(row, "n_forced", np.nan), np.nan)
            if np.isfinite(n_forced) and n_psf > 0:
                frac = 100.0 * n_forced / n_psf
                fitem = QTableWidgetItem(f"{frac:.0f}%")
                fitem.setToolTip(
                    f"{int(n_forced)} / {int(n_psf)} — 그 프레임에서 검출되지 않아 "
                    "마스터 위치로 강제 측광한 별.\n높으면 PSF vs 구경 비교를 "
                    "검출된 별로만 해야 한다(강제 위치는 구경이 거의 0)."
                )
            else:
                fitem = QTableWidgetItem("")
            self.frame_table.setItem(r, 6, fitem)
            qc_status = str(getattr(row, "psf_qc_status", "") or "")
            qc_item = QTableWidgetItem(qc_status)
            qc_item.setToolTip(str(getattr(row, "psf_qc_reasons", "") or ""))
            self.frame_table.setItem(r, 7, qc_item)
            elapsed = _safe_float(getattr(row, "frame_total_elapsed_s", np.nan), np.nan)
            self.frame_table.setItem(
                r, 8, QTableWidgetItem(f"{elapsed:.1f} s" if np.isfinite(elapsed) else "")
            )
            if qc_status == "FAIL":
                set_table_row_background(self.frame_table, r, status_row_background(False))
            elif qc_status == "REVIEW":
                set_table_row_background(
                    self.frame_table, r, status_row_background(True, warning=True)
                )
            else:
                has_good_phot = int(
                    _safe_float(getattr(row, "n_goodmag", 0), 0)
                ) > 0
                set_table_row_background(
                    self.frame_table, r, status_row_background(has_good_phot)
                )

    # ── QC Report ─────────────────────────────────────────────────────────────

    def _refresh_qc(self):
        """Compute PSF QC statistics and update the QC tab text + Ap vs PSF plot."""
        if not hasattr(self, "qc_text"):
            return
        psf_dir = step8_psf_dir(self.params.P.result_dir)
        idx_path = psf_dir / "photometry_index.csv"
        export_inputs = None
        if not idx_path.exists():
            self.qc_text.setPlainText("photometry_index.csv not found.\nRun Step 8 first.")
            self._cmp_merged_df = None
            self._plot_mag_comparison()
            self._refresh_final_diagnostics()
            return
        try:
            # 헤드리스 러너와 **같은 함수**로 읽는다 — 창과 배치가 서로 다른
            # 코드로 QC 를 만들면 둘의 산출물이 조용히 갈라진다.
            idx, all_df, meta_df = load_psf_qc_inputs(psf_dir)
            if all_df.empty:
                self.qc_text.setPlainText("No photometry TSV files found.")
                self._refresh_final_diagnostics()
                return

            good = all_df[all_df["flags_psf"] == 0].copy() if "flags_psf" in all_df.columns else all_df.copy()
            filters = sorted(all_df["FILTER"].dropna().unique().tolist()) if "FILTER" in all_df.columns else []
            n_total = len(all_df)
            n_clean = len(good)

            W = 60
            lines = []
            lines.append("─" * W)
            lines.append("  PSF Photometry QC Report")
            lines.append("─" * W)
            filt_counts = "  ".join(
                f"{f}:{(idx['filter'] == f).sum()}" for f in filters if "filter" in idx.columns
            )
            lines.append(f"  총 프레임    : {len(idx)}  ({filt_counts})")
            lines.append(f"  총 검출 소스 : {n_total:,}")
            lines.append(f"  flags=0      : {n_clean:,} / {n_total:,} = {100 * n_clean / max(n_total, 1):.1f}%")
            lines.append("")

            # 1. Filter stats
            lines.append("  1. 필터별 검출 통계")
            lines.append(f"  {'필터':^4}  {'프레임':^6}  {'평균검출':^8}  {'성공률':^7}  {'실패율':^6}  {'iter2추가':^9}")
            lines.append("  " + "─" * 52)
            for filt in filters:
                si = idx[idx["filter"] == filt] if "filter" in idx.columns else pd.DataFrame()
                if si.empty:
                    continue
                avg_n = si["n"].mean() if "n" in si.columns else 0
                avg_g = si["n_goodmag"].mean() if "n_goodmag" in si.columns else avg_n
                ok_pct = 100 * avg_g / max(avg_n, 1)
                avg_new = si["n_new_iter"].mean() if "n_new_iter" in si.columns else 0
                lines.append(
                    f"  {filt:^4}  {len(si):^6}  {avg_n:^8.0f}  {ok_pct:^6.1f}%  "
                    f"{100 - ok_pct:^5.1f}%  avg {avg_new:.1f}"
                )
            lines.append("")

            # 2. Mag range & error
            lines.append("  2. 등급 범위 (mag_psf, flags=0)")
            lines.append(f"  {'필터':^4}  {'범위':^15}  {'평균':^6}  {'중앙값':^6}  {'σ':^5}  {'err중앙값':^9}")
            lines.append("  " + "─" * 58)
            for filt in filters:
                sub = good[good["FILTER"] == filt] if "FILTER" in good.columns else pd.DataFrame()
                mag = sub["mag_psf"].dropna() if "mag_psf" in sub.columns else pd.Series()
                err = sub["mag_psf_err"].dropna() if "mag_psf_err" in sub.columns else pd.Series()
                if mag.empty:
                    continue
                lines.append(
                    f"  {filt:^4}  {mag.min():.2f} ~ {mag.max():.2f}  "
                    f"{mag.mean():^6.2f}  {mag.median():^6.2f}  {mag.std():^5.2f}  "
                    f"{err.median():^9.4f}" if not err.empty else
                    f"  {filt:^4}  {mag.min():.2f} ~ {mag.max():.2f}  "
                    f"{mag.mean():^6.2f}  {mag.median():^6.2f}  {mag.std():^5.2f}  {'N/A':^9}"
                )
            lines.append("")

            # 3. SNR
            if "snr_psf" in good.columns:
                lines.append("  3. SNR 분포 (flags=0)")
                lines.append(f"  {'필터':^4}  {'10%':^8}  {'median':^8}  {'90%':^8}")
                lines.append("  " + "─" * 36)
                for filt in filters:
                    sub = good[good["FILTER"] == filt]["snr_psf"].dropna() if "FILTER" in good.columns else pd.Series()
                    if sub.empty:
                        continue
                    lines.append(
                        f"  {filt:^4}  {np.percentile(sub, 10):^8.1f}  "
                        f"{sub.median():^8.1f}  {np.percentile(sub, 90):^8.1f}"
                    )
                lines.append("")

            # 4. qfit
            if "qfit" in good.columns:
                lines.append("  4. PSF 적합 품질 (qfit, flags=0)")
                lines.append(f"  {'필터':^4}  {'중앙값':^8}  {'>5 비율':^8}")
                lines.append("  " + "─" * 28)
                for filt in filters:
                    sub = good[good["FILTER"] == filt]["qfit"].dropna() if "FILTER" in good.columns else pd.Series()
                    if sub.empty:
                        continue
                    bad_pct = 100 * (sub > 5).sum() / max(len(sub), 1)
                    warn = " ⚠" if bad_pct > 5 else ""
                    lines.append(f"  {filt:^4}  {sub.median():^8.3f}  {bad_pct:^7.1f}%{warn}")
                lines.append("")

            # 5. Residual STD
            if not meta_df.empty and "residual_std" in meta_df.columns:
                lines.append("  5. Residual STD (ADU, per frame mean)")
                lines.append(f"  {'필터':^4}  {'iter1':^10}  {'iter2':^10}")
                lines.append("  " + "─" * 30)
                for filt in filters:
                    i1 = meta_df[(meta_df["filter"] == filt) & (meta_df["iter"] == 1)]["residual_std"]
                    i2 = meta_df[(meta_df["filter"] == filt) & (meta_df["iter"] == 2)]["residual_std"]
                    if i1.empty:
                        continue
                    i2_mean = f"{i2.mean():.2f}" if not i2.empty else "N/A"
                    lines.append(f"  {filt:^4}  {i1.mean():^10.2f}  {i2_mean:^10}")
                lines.append("")

            self.qc_text.setPlainText("\n".join(lines))
            export_inputs = (idx, all_df, meta_df)
        except Exception as e:
            self.qc_text.setPlainText(f"QC 생성 오류: {e}")

        self._cmp_merged_df = None
        self._plot_mag_comparison()
        self._refresh_final_diagnostics()
        if export_inputs is not None:
            self._export_psf_qc_products(*export_inputs)

    def _refresh_final_diagnostics(self, *_args) -> None:
        """Refresh the frame selector and redraw the selected final diagnostic."""
        if not hasattr(self, "final_diag_frame_combo"):
            return
        psf_dir = step8_psf_dir(self.params.P.result_dir)
        frames = [
            path.name[len("photometry_"):-len(".tsv")]
            for path in sorted(psf_dir.glob("photometry_*.tsv"))
        ]
        current = self.final_diag_frame_combo.currentText()
        self.final_diag_frame_combo.blockSignals(True)
        self.final_diag_frame_combo.clear()
        self.final_diag_frame_combo.addItems(frames)
        if current in frames:
            self.final_diag_frame_combo.setCurrentText(current)
        elif frames:
            self.final_diag_frame_combo.setCurrentIndex(0)
        self.final_diag_frame_combo.blockSignals(False)
        self._plot_final_diagnostics(self.final_diag_frame_combo.currentText())

    def _show_final_diagnostic_message(self, message: str) -> None:
        self.final_diag_fig.clear()
        ax = self.final_diag_fig.add_subplot(111)
        ax.set_axis_off()
        ax.text(0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes, color="0.4")
        self.final_diag_status.setText(message)
        self.final_diag_status.setStyleSheet("QLabel { color: #616161; }")
        self.final_diag_canvas.draw_idle()

    def _final_diagnostic_meta(self, fname: str) -> dict:
        meta = self._residual_meta.get(fname, {})
        if isinstance(meta, dict) and meta:
            return meta
        path = step8_psf_dir(self.params.P.result_dir) / f"residual_meta_{fname}.json"
        if not path.exists():
            return {}
        try:
            meta = json.loads(path.read_text(encoding="utf-8"))
            return meta if isinstance(meta, dict) else {}
        except Exception:
            return {}

    def _final_diagnostic_epsf(self, fname: str, meta: dict) -> tuple[np.ndarray | None, Path | None]:
        psf_dir = step8_psf_dir(self.params.P.result_dir)
        filt = str(meta.get("filter", "")).strip()
        stem = Path(fname).stem
        candidates = []
        if filt:
            candidates.extend(
                [
                    psf_dir / f"epsf_model_{filt}_{stem}.fits",
                    psf_dir / f"epsf_model_{filt.lower()}_{stem}.fits",
                    psf_dir / f"epsf_model_{filt.upper()}_{stem}.fits",
                    psf_dir / f"epsf_model_{filt}.fits",
                    psf_dir / f"epsf_model_{filt.lower()}.fits",
                    psf_dir / f"epsf_model_{filt.upper()}.fits",
                ]
            )
        candidates.extend(sorted(psf_dir.glob(f"epsf_model_*_{stem}.fits")))
        seen: set[Path] = set()
        for path in candidates:
            if path in seen or not path.exists():
                continue
            seen.add(path)
            try:
                return np.asarray(fits.getdata(path), dtype=float), path
            except Exception:
                continue
        return None, None

    def _final_diagnostic_reference_catalog(self, fname: str, meta: dict) -> pd.DataFrame:
        reference = meta.get("epsf_reference", {})
        catalog_name = reference.get("catalog_path", "") if isinstance(reference, dict) else ""
        path = step8_psf_dir(self.params.P.result_dir) / (
            str(catalog_name).strip() or f"epsf_reference_{fname}.csv"
        )
        if not path.exists():
            return pd.DataFrame()
        try:
            return pd.read_csv(path)
        except Exception:
            return pd.DataFrame()

    def _final_diagnostic_pixel_scale(self, fname: str) -> float:
        fits_path = self._resolve_fits_path_window(fname)
        if fits_path is not None:
            try:
                from astropy.wcs import WCS
                from astropy.wcs.utils import proj_plane_pixel_scales

                celestial = WCS(fits.getheader(fits_path)).celestial
                scales = np.asarray(proj_plane_pixel_scales(celestial), dtype=float) * 3600.0
                scale = float(np.nanmedian(np.abs(scales)))
                if np.isfinite(scale) and scale > 0:
                    return scale
            except Exception:
                pass

        for path in (
            self._cache_dir_path() / f"detect_{fname}.json",
            step4_dir(self.params.P.result_dir) / f"detect_{fname}.json",
        ):
            if not path.exists():
                continue
            try:
                meta = json.loads(path.read_text(encoding="utf-8"))
                fwhm_arcsec = _safe_float(meta.get("fwhm_arcsec"), np.nan)
                fwhm_px = _safe_float(
                    meta.get("fwhm_px", meta.get("fwhm_med_px", meta.get("fwhm_med"))),
                    np.nan,
                )
                scale = fwhm_arcsec / fwhm_px
                if np.isfinite(scale) and scale > 0:
                    return float(scale)
            except Exception:
                continue
        return np.nan

    def _plot_final_diagnostics(self, fname: str | None = None) -> None:
        """Draw the six-panel final diagnostic for one Step 8 frame."""
        if not hasattr(self, "final_diag_fig"):
            return
        if not isinstance(fname, str) or not fname:
            fname = self.final_diag_frame_combo.currentText()
        if not fname:
            self._final_diag_data = pd.DataFrame()
            self._final_diag_summary = {}
            self._show_final_diagnostic_message("No PSF result is available.")
            return

        try:
            # 헤드리스와 **같은 함수**로 그린다 — 조립을 창이 따로 하면 두
            # 산출물이 갈라진다.
            data, summary = render_psf_final_diagnostics(
                self.final_diag_fig,
                self.params,
                Path(self.params.P.result_dir),
                fname,
                use_cropped=bool(getattr(self, "use_cropped", False)),
            )
            self._final_diag_data = data
            self._final_diag_summary = summary
            status = str(summary.get("status", "CHECK"))
            status_parts = [
                status,
                f"N={int(summary.get('n_matched', 0))}",
                f"offset={_safe_float(summary.get('high_snr_reference_offset_mag')):+.3f} mag",
                f"low-SNR={_safe_float(summary.get('low_snr_5_10_median_centered_mag')):+.3f} mag",
                f"high-SNR scatter={_safe_float(summary.get('high_snr_robust_scatter_mag')):.3f} mag",
                f"e={_safe_float(summary.get('epsf_ellipticity')):.3f}",
                f"A180={_safe_float(summary.get('epsf_rotation_asymmetry')):.3f}",
                f"refs={int(summary.get('epsf_reference_n', 0))}",
            ]
            if bool(summary.get("psf_aperture_scale_applied", False)):
                status_parts.append(
                    f"scale={_safe_float(summary.get('psf_aperture_scale'), 1.0):.4f} "
                    f"(N={_to_int(summary.get('psf_aperture_scale_n', 0), 0)})"
                )
            warnings = summary.get("warnings", [])
            if isinstance(warnings, list) and warnings:
                status_parts.append("; ".join(str(item) for item in warnings))
            self.final_diag_status.setText(" | ".join(status_parts))
            color = "#2E7D32" if status == "OK" else "#E65100"
            self.final_diag_status.setStyleSheet(f"QLabel {{ color: {color}; font-weight: bold; }}")
            self.final_diag_canvas.draw_idle()
        except Exception as exc:
            self._final_diag_data = pd.DataFrame()
            self._final_diag_summary = {"file": fname, "status": "ERROR", "error": str(exc)}
            self._show_final_diagnostic_message(f"Final diagnostics failed for {fname}:\n{exc}")
            try:
                self.log(f"Final diagnostics failed for {fname}: {exc}")
            except Exception:
                pass

    def _export_psf_qc_products(
        self,
        idx: pd.DataFrame,
        all_df: pd.DataFrame,
        meta_df: pd.DataFrame,
    ) -> list[Path]:
        """Export reproducible Step 8 QC products for papers and run audits."""
        psf_dir = step8_psf_dir(self.params.P.result_dir)
        psf_dir.mkdir(parents=True, exist_ok=True)

        saved: list[Path] = []
        summary = _build_psf_qc_summary(
            idx,
            all_df,
            meta_df,
            getattr(self, "_cmp_merged_df", None),
        )
        if not summary.empty:
            summary_path = psf_dir / "psf_qc_summary.csv"
            summary.to_csv(summary_path, index=False)
            saved.append(summary_path)

        frame_qc = _build_psf_frame_qc_table(idx, meta_df)
        if not frame_qc.empty:
            frame_qc_path = psf_dir / "psf_frame_qc.csv"
            frame_qc.to_csv(frame_qc_path, index=False)
            saved.append(frame_qc_path)

            frame_fig = Figure(figsize=(10.5, 6.8), dpi=120)
            if _draw_psf_frame_qc_overview(frame_fig, frame_qc):
                frame_fig_path = psf_dir / "step8_residual_core_qc.png"
                frame_fig.savefig(frame_fig_path, dpi=160, bbox_inches="tight")
                saved.append(frame_fig_path)

        report = self.qc_text.toPlainText() if hasattr(self, "qc_text") else ""
        if report.strip():
            report_path = psf_dir / "psf_qc_report.txt"
            report_path.write_text(report, encoding="utf-8")
            saved.append(report_path)

        if hasattr(self, "cmp_fig"):
            fig_path = psf_dir / "step8_ap_vs_psf_comparison.png"
            self.cmp_fig.savefig(fig_path, dpi=160, bbox_inches="tight")
            saved.append(fig_path)

        final_summary = getattr(self, "_final_diag_summary", {})
        final_data = getattr(self, "_final_diag_data", pd.DataFrame())
        if isinstance(final_summary, dict) and final_summary.get("file"):
            stem = Path(str(final_summary["file"])).stem
            final_fig_path = psf_dir / f"step8_final_diagnostics_{stem}.png"
            self.final_diag_fig.savefig(final_fig_path, dpi=160, bbox_inches="tight")
            saved.append(final_fig_path)

            summary_path = psf_dir / f"psf_final_diagnostics_{stem}.json"
            summary_path.write_text(
                json.dumps(final_summary, ensure_ascii=False, indent=2, allow_nan=True),
                encoding="utf-8",
            )
            saved.append(summary_path)
            if isinstance(final_data, pd.DataFrame) and not final_data.empty:
                data_path = psf_dir / f"psf_final_diagnostics_{stem}.csv"
                final_data.to_csv(data_path, index=False)
                saved.append(data_path)

        if saved:
            try:
                names = ", ".join(path.name for path in saved)
                self.log(f"Step8 QC products exported: {names}")
            except Exception:
                pass
        return saved

    # ── Aperture vs PSF magnitude comparison ──────────────────────────────────

    def _load_or_build_comparison(self) -> tuple[pd.DataFrame, int]:
        """구경 vs PSF 병합표를 디스크 캐시에서 읽고, 없거나 낡았으면 다시 만든다.

        이 병합은 프레임 수에 비례해 무겁다 — M13 15프레임에 **10.9초**이고,
        창을 열 때마다 `_refresh_qc()` 가 캐시를 버리고 처음부터 다시 만든다.
        그래서 Step 8 창 로드가 17.8초였다. Step 8 산출물이 그대로면 결과도
        같으므로, `photometry_index.csv` 보다 새 캐시가 있으면 그것을 쓴다.
        """
        psf_dir = step8_psf_dir(self.params.P.result_dir)
        cache_path = psf_dir / "psf_ap_vs_psf.csv"
        meta_path = psf_dir / "psf_ap_vs_psf_meta.json"
        index_path = psf_dir / "photometry_index.csv"
        try:
            if (
                cache_path.exists()
                and index_path.exists()
                and cache_path.stat().st_mtime >= index_path.stat().st_mtime
            ):
                cached = pd.read_csv(cache_path)
                n_split = 0
                if meta_path.exists():
                    n_split = int(
                        json.loads(meta_path.read_text(encoding="utf-8")).get(
                            "split_excluded_total", 0
                        )
                    )
                return cached, n_split
        except Exception:
            pass

        merged, split_excluded_total = build_ap_psf_comparison(
            self.params, self.params.P.result_dir
        )
        try:
            psf_dir.mkdir(parents=True, exist_ok=True)
            merged.to_csv(cache_path, index=False)
            meta_path.write_text(
                json.dumps({"split_excluded_total": int(split_excluded_total)}),
                encoding="utf-8",
            )
        except Exception:
            pass
        return merged, int(split_excluded_total)

    def _plot_mag_comparison(self):  # noqa: C901
        """Scatter: mag_ap (Step5) vs mag_psf (Step6), merged on det_uid."""
        if not hasattr(self, "cmp_fig"):
            return

        _FILT_COLORS = {
            "u": "#9467bd", "g": "#2ca02c", "r": "#d62728",
            "i": "#ff7f0e", "z": "#8c564b", "b": "#1f77b4",
            "v": "#bcbd22", "ha": "#e377c2",
        }

        # 병합은 헤드리스와 **같은 함수**로 한다 — 창과 배치가 각자 병합하면
        # 두 산출물이 조용히 갈라진다. 캐시(_cmp_merged_df)는 창 쪽 사정이다.
        if not hasattr(self, "_cmp_merged_df") or self._cmp_merged_df is None:
            merged, split_excluded_total = self._load_or_build_comparison()
            self._cmp_merged_df = merged
            self._cmp_split_excluded_total = int(split_excluded_total)

        self.cmp_fig.clf()

        def _empty(msg):
            ax = self.cmp_fig.add_subplot(111)
            ax.text(0.5, 0.5, msg, transform=ax.transAxes,
                    ha="center", va="center", fontsize=10, color="gray")
            self.cmp_canvas.draw_idle()
            self.cmp_stats_label.setText(msg)

        if self._cmp_merged_df.empty:
            _empty("No matched data.\nRun Step 7 and Step 8 first.")
            return

        df = self._cmp_merged_df.copy()

        # Refresh filter/frame selectors from available merged data.
        if hasattr(self, "cmp_filter_combo"):
            try:
                _fvals = sorted(df["FILTER"].dropna().astype(str).unique().tolist()) if "FILTER" in df.columns else []
                _cur = self.cmp_filter_combo.currentText().strip() or "all"
                self.cmp_filter_combo.blockSignals(True)
                self.cmp_filter_combo.clear()
                self.cmp_filter_combo.addItem("all")
                for _v in _fvals:
                    self.cmp_filter_combo.addItem(_v)
                self.cmp_filter_combo.setCurrentText(_cur if _cur in (["all"] + _fvals) else "all")
                self.cmp_filter_combo.blockSignals(False)
            except Exception:
                pass
        if hasattr(self, "cmp_frame_combo"):
            try:
                _frames = sorted(df["FRAME"].dropna().astype(str).unique().tolist()) if "FRAME" in df.columns else []
                _cur = self.cmp_frame_combo.currentText().strip() or "all"
                self.cmp_frame_combo.blockSignals(True)
                self.cmp_frame_combo.clear()
                self.cmp_frame_combo.addItem("all")
                for _v in _frames:
                    self.cmp_frame_combo.addItem(_v)
                self.cmp_frame_combo.setCurrentText(_cur if _cur in (["all"] + _frames) else "all")
                self.cmp_frame_combo.blockSignals(False)
            except Exception:
                pass

        df["mag_ap"] = pd.to_numeric(df.get("mag_ap"), errors="coerce")
        df["mag_psf"] = pd.to_numeric(df.get("mag_psf"), errors="coerce")
        df = df[np.isfinite(df["mag_ap"]) & np.isfinite(df["mag_psf"])].copy()

        if len(df) == 0:
            _empty("All magnitudes are NaN.\nCheck Step 7 forced photometry and Step 8 PSF outputs.")
            return

        df["delta"] = df["mag_ap"] - df["mag_psf"]
        n_before = int(len(df))

        # Pre-convert numeric filter columns once
        if "flags_psf" in df.columns:
            df["flags_psf"] = pd.to_numeric(df["flags_psf"], errors="coerce")
        if "snr_psf" in df.columns:
            df["snr_psf"] = pd.to_numeric(df["snr_psf"], errors="coerce")
        if "qfit" in df.columns:
            df["qfit"] = pd.to_numeric(df["qfit"], errors="coerce")
        if "qfit_noise_ratio" in df.columns:
            df["qfit_noise_ratio"] = pd.to_numeric(
                df["qfit_noise_ratio"], errors="coerce"
            )

        # Selector filters (filter/frame)
        if hasattr(self, "cmp_filter_combo") and "FILTER" in df.columns:
            _fsel = str(self.cmp_filter_combo.currentText()).strip()
            if _fsel and _fsel.lower() != "all":
                _fkey = normalize_filter_name(_fsel)
                df = df[df["FILTER"].astype(str).map(normalize_filter_name) == _fkey].copy()
        if hasattr(self, "cmp_frame_combo") and "FRAME" in df.columns:
            _rsel = str(self.cmp_frame_combo.currentText()).strip()
            if _rsel and _rsel.lower() != "all":
                df = df[df["FRAME"].astype(str) == _rsel].copy()

        # User filters
        if getattr(self, "cmp_flags0_only", None) is not None and self.cmp_flags0_only.isChecked() and "flags_psf" in df.columns:
            df = df[np.isfinite(df["flags_psf"]) & (df["flags_psf"] == 0)].copy()
        if getattr(self, "cmp_snr_min", None) is not None:
            _snr_min = float(self.cmp_snr_min.value())
            if _snr_min > 0 and "snr_psf" in df.columns:
                df = df[np.isfinite(df["snr_psf"]) & (df["snr_psf"] >= _snr_min)].copy()
        if getattr(self, "cmp_qfit_max", None) is not None:
            _qmax = float(self.cmp_qfit_max.value())
            _qfit_column = (
                "qfit_noise_ratio"
                if "qfit_noise_ratio" in df.columns
                and np.isfinite(df["qfit_noise_ratio"]).any()
                else "qfit"
            )
            if _qmax > 0 and _qfit_column in df.columns:
                df = df[
                    np.isfinite(df[_qfit_column]) & (df[_qfit_column] <= _qmax)
                ].copy()
        if getattr(self, "cmp_dmag_clip", None) is not None:
            _dclip = float(self.cmp_dmag_clip.value())
            if _dclip > 0:
                df = df[np.isfinite(df["delta"]) & (np.abs(df["delta"]) <= _dclip)].copy()

        if len(df) == 0:
            _empty("No data after filters.\nRelax SNR/qfit/|Δmag| settings.")
            return

        filt_col = "FILTER" if "FILTER" in df.columns else None

        ax1 = self.cmp_fig.add_subplot(121)
        ax2 = self.cmp_fig.add_subplot(122)

        stats_parts = []
        groups = df.groupby(filt_col) if filt_col else [("all", df)]
        for filt, sub in groups:
            color = _FILT_COLORS.get(str(filt).lower(), "#999999")
            ax1.scatter(sub["mag_ap"], sub["mag_psf"],
                        s=4, alpha=0.35, color=color, label=str(filt), rasterized=True)
            ax2.scatter(sub["mag_ap"], sub["delta"],
                        s=4, alpha=0.35, color=color, label=str(filt), rasterized=True)
            n = len(sub)
            med = float(np.nanmedian(sub["delta"]))
            std = float(np.nanstd(sub["delta"]))
            stats_parts.append(f"{filt}: N={n}  Δmed={med:+.3f}  σ={std:.3f}")

        # 1:1 reference line (ax1)
        all_mag = np.concatenate([df["mag_ap"].values, df["mag_psf"].values])
        lo, hi = np.nanmin(all_mag) - 0.2, np.nanmax(all_mag) + 0.2
        ax1.plot([lo, hi], [lo, hi], "k--", lw=0.8, alpha=0.5, zorder=0)
        ax1.set_xlim(lo, hi)
        ax1.set_ylim(lo, hi)
        ax1.set_xlabel("mag_ap", fontsize=9)
        ax1.set_ylabel("mag_psf", fontsize=9)
        ax1.set_title("Aperture vs PSF magnitude", fontsize=9)
        ax1.legend(fontsize=7, markerscale=2, loc="upper left")

        # Zero ± 0.05 reference lines (ax2) + robust median guide
        dmed_all = float(np.nanmedian(df["delta"]))
        ax2.axhline(
            dmed_all,
            color="#D62728",
            lw=2.0,
            ls="-",
            alpha=0.95,
            zorder=1,
            label=f"Δmag median {dmed_all:+.3f}",
        )
        ax2.axhline(0.0,  color="k",    lw=0.8, ls="--", alpha=0.6, zorder=0)
        ax2.axhline(+0.05, color="gray", lw=0.5, ls=":",  alpha=0.5, zorder=0)
        ax2.axhline(-0.05, color="gray", lw=0.5, ls=":",  alpha=0.5, zorder=0)
        ax2.set_xlabel("mag_ap", fontsize=9)
        ax2.set_ylabel("Δmag  (Ap − PSF)", fontsize=9)
        ax2.set_title("Δmag vs mag_ap", fontsize=9)
        ax2.legend(fontsize=7, markerscale=2, loc="upper left")

        self.cmp_fig.tight_layout()
        self.cmp_canvas.draw_idle()
        self.cmp_stats_label.setText(
            f"N={len(df)}/{n_before}  |  split_excluded={int(getattr(self, '_cmp_split_excluded_total', 0))}  |  "
            + "  |  ".join(stats_parts)
        )

    # ── Load existing results from disk (called on restore_state) ─────────────

    def _load_from_disk(self):
        """Reload EPSF models and residual images from disk into memory caches."""
        out_dir = step8_psf_dir(self.params.P.result_dir)
        if not out_dir.exists():
            return
        self._last_epsf.clear()
        self._residual_meta.clear()
        self._last_new_xy.clear()
        self.epsf_filter_combo.clear()
        self.res_file_combo.clear()
        self.res_iter_combo.clear()

        def _epsf_display_key_from_path(epsf_path: Path) -> str:
            stem = epsf_path.stem  # epsf_model_{filter}_{frame_stem} or epsf_model_{filter}
            body = stem.replace("epsf_model_", "", 1)
            if "_" not in body:
                return body
            filt, frame_stem = body.split("_", 1)
            return f"{frame_stem} | {filt}"

        # Load EPSF FITS files
        for epsf_path in out_dir.glob("epsf_model_*.fits"):
            try:
                display_key = _epsf_display_key_from_path(epsf_path)
                if not display_key:
                    continue
                arr = fits.getdata(str(epsf_path)).astype(float)
                self._last_epsf[display_key] = arr
            except Exception:
                pass

        if self._last_epsf:
            self.epsf_filter_combo.blockSignals(True)
            self.epsf_filter_combo.clear()
            self.epsf_filter_combo.addItems(sorted(self._last_epsf.keys()))
            self.epsf_filter_combo.setCurrentIndex(0)
            self.epsf_filter_combo.blockSignals(False)
            first_filter = self.epsf_filter_combo.currentText()
            if first_filter:
                self._plot_epsf(first_filter)

        # Load residual metadata (preferred, supports iteration-wise view).
        for meta_path in sorted(out_dir.glob("residual_meta_*.json")):
            try:
                name = meta_path.name
                if not name.startswith("residual_meta_") or not name.endswith(".json"):
                    continue
                fname = name[len("residual_meta_"):-len(".json")]
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                if isinstance(meta, dict):
                    self._residual_meta[fname] = meta
                    self._last_new_xy.setdefault(fname, None)
            except Exception:
                pass

        if self._residual_meta:
            self.res_file_combo.blockSignals(True)
            self.res_file_combo.clear()
            self.res_file_combo.addItems(sorted(self._residual_meta.keys()))
            self.res_file_combo.setCurrentIndex(0)
            self.res_file_combo.blockSignals(False)
            first_fname = self.res_file_combo.currentText()
            if first_fname:
                self._refresh_residual_iter_combo(first_fname)
                self._plot_cutout(first_fname)

        self._refresh_qc()  # refresh QC tab (stats + Ap vs PSF plot)

    # ── Parameters dialog ─────────────────────────────────────────────────────

    def open_parameters_dialog(self):
        dialog = FittedDialog(self)
        configure_parameter_dialog(dialog, "Step 8 PSF Parameters", 620, 720)
        layout = QVBoxLayout(dialog)

        scroll = QScrollArea(dialog)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        content = QWidget()
        body = QVBoxLayout(content)
        body.setContentsMargins(4, 4, 4, 4)
        body.setSpacing(8)
        scroll.setWidget(content)
        layout.addWidget(scroll, 1)

        _info = QLabel("Adjust PSF photometry parameters. Changes apply to the next run.")
        _info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 10px; margin-bottom: 10px; }")
        _info.setWordWrap(True)
        body.addWidget(_info)

        def _add_group(title: str, *, expanded: bool = False) -> QFormLayout:
            group, container = create_collapsible_section(title, initial_expanded=expanded)
            form = QFormLayout(container)
            form.setLabelAlignment(Qt.AlignRight)
            form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
            body.addWidget(group)
            return form

        mode_form = _add_group("Mode", expanded=True)
        epsf_form = _add_group("ePSF Model", expanded=True)
        scale_form = _add_group("PSF Flux Scale")
        fit_form = _add_group("PSF Fit")
        core_form = _add_group("Crowded Core Cut")
        redetect_form = _add_group("Residual Re-detection")
        output_form = _add_group("Output")

        # ── Field mode preset ──────────────────────────────────────────────
        mode_combo = QComboBox()
        mode_combo.addItem("Normal (일반)", "normal")
        mode_combo.addItem("Crowded (구상성단/혼잡장)", "crowded")
        mode_combo.addItem("Faint (희미한 필드)", "faint")
        mode_combo.addItem("Custom (수동)", "custom")
        _saved_mode = str(getattr(self.params.P, "psf_mode", "normal"))
        _mi = mode_combo.findData(_saved_mode)
        mode_combo.setCurrentIndex(_mi if _mi >= 0 else 0)
        mode_form.addRow("Field mode:", mode_combo)

        self.p_model_mode = QComboBox()
        self.p_model_mode.addItems(["per_frame"])
        self.p_model_mode.setCurrentText("per_frame")

        self.p_fit_engine = QComboBox()
        self.p_fit_engine.addItem("APEX iterative - CPU, recommended", "apex_iterative")
        self.p_fit_engine.addItem("photutils - validation, slower", "photutils")
        _eng = str(
            getattr(self.params.P, "psf_fit_engine", "apex_iterative")
        ).strip().lower()
        if _eng == "allstar":
            _eng = "apex_iterative"
        _ei = self.p_fit_engine.findData(_eng)
        self.p_fit_engine.setCurrentIndex(_ei if _ei >= 0 else 0)
        mode_form.addRow("Fit engine:", self.p_fit_engine)

        self.p_build_mode = QComboBox()
        self.p_build_mode.addItem("epsf  —  EPSFBuilder 경험적 PSF", "epsf")
        self.p_build_mode.addItem("moffat  —  해석 Moffat (γ·β 적합)", "moffat")
        self.p_build_mode.addItem("moffat+residual  —  해석 + 잔차격자", "moffat_hybrid")
        self.p_build_mode.setToolTip(
            "epsf: 별 이미지를 그대로 평균낸 경험적 모형. 어떤 모양이든 담지만"
            " 가파른 핵심까지 보간해야 한다.\n"
            "moffat: γ·β 를 적합한 해석 모형. 핵심을 정확히 계산하지만"
            " 원형이라 늘어난 별은 못 맞춘다.\n"
            "moffat+residual: 핵심은 해석식으로 정확히, 남은 모양은 격자로."
            " DAOPHOT 과 같은 구성."
        )
        _bm = str(getattr(self.params.P, "psf_build_mode", "epsf")).strip().lower()
        _bi = self.p_build_mode.findData(_bm)
        self.p_build_mode.setCurrentIndex(_bi if _bi >= 0 else 0)

        self.p_workers = QSpinBox()
        self.p_workers.setRange(0, 64)
        self.p_workers.setValue(_to_int(getattr(self.params.P, "psf_parallel_workers", 0), 0))
        self.p_workers.setToolTip("0 = auto/global parallel workers")
        mode_form.addRow("PSF workers (0=auto):", self.p_workers)

        self.p_oversampling = QSpinBox()
        self.p_oversampling.setRange(1, 8)
        self.p_oversampling.setValue(_to_int(getattr(self.params.P, "psf_epsf_oversampling", 2), 2))
        epsf_form.addRow("EPSF oversampling:", self.p_oversampling)

        self.p_epsf_mult = QDoubleSpinBox()
        self.p_epsf_mult.setRange(1.0, 10.0)
        self.p_epsf_mult.setSingleStep(0.1)
        self.p_epsf_mult.setValue(_to_float(getattr(self.params.P, "psf_epsf_size_fwhm_mult", 4.0), 4.0))
        epsf_form.addRow("EPSF cutout (×FWHM):", self.p_epsf_mult)

        self.p_n_stars = QSpinBox()
        self.p_n_stars.setRange(0, 500)
        self.p_n_stars.setSpecialValueText("auto")
        self.p_n_stars.setValue(_to_int(getattr(self.params.P, "psf_n_stars_max", 0), 0))
        self.p_n_stars.setToolTip(
            "0 = automatic per-frame budget; positive values are a CPU/memory cap (minimum 3)"
        )
        epsf_form.addRow("PSF star cap:", self.p_n_stars)

        self.p_isolation = QDoubleSpinBox()
        self.p_isolation.setRange(1.0, 10.0)
        self.p_isolation.setSingleStep(0.5)
        self.p_isolation.setValue(_to_float(getattr(self.params.P, "psf_isolation_fwhm_mult", 3.0), 3.0))
        epsf_form.addRow("Isolation (×FWHM):", self.p_isolation)

        self.p_epsf_contamination_filter = QCheckBox(
            "Reject locally contaminated/core ePSF reference stars"
        )
        self.p_epsf_contamination_filter.setChecked(
            _as_bool(getattr(self.params.P, "psf_epsf_contamination_filter", True), True)
        )
        self.p_epsf_contamination_filter.setToolTip(
            "Uses all Step 4 detections, local annulus residuals, and the automatic "
            "cluster-core estimate only for ePSF reference-star selection."
        )
        epsf_form.addRow("", self.p_epsf_contamination_filter)

        self.p_flux_scale_correction = QCheckBox(
            "Anchor PSF fluxes to clean Step 7 aperture references"
        )
        self.p_flux_scale_correction.setChecked(
            _as_bool(getattr(self.params.P, "psf_flux_scale_correction", False), False)
        )
        self.p_flux_scale_correction.setToolTip(
            "Applies one robust per-frame multiplicative scale after PSF fitting. "
            "Raw PSF fluxes are preserved in separate output columns."
        )
        scale_form.addRow("", self.p_flux_scale_correction)

        self.p_flux_scale_min_snr = QDoubleSpinBox()
        self.p_flux_scale_min_snr.setRange(5.0, 1000.0)
        self.p_flux_scale_min_snr.setSingleStep(10.0)
        self.p_flux_scale_min_snr.setValue(
            _to_float(getattr(self.params.P, "psf_flux_scale_min_snr", 50.0), 50.0)
        )
        scale_form.addRow("Minimum aperture SNR:", self.p_flux_scale_min_snr)

        self.p_flux_scale_min_stars = QSpinBox()
        self.p_flux_scale_min_stars.setRange(3, 500)
        self.p_flux_scale_min_stars.setValue(
            _to_int(getattr(self.params.P, "psf_flux_scale_min_stars", 8), 8)
        )
        scale_form.addRow("Minimum references:", self.p_flux_scale_min_stars)

        self.p_flux_scale_min_neighbor = QDoubleSpinBox()
        self.p_flux_scale_min_neighbor.setRange(0.0, 20.0)
        self.p_flux_scale_min_neighbor.setSingleStep(0.5)
        self.p_flux_scale_min_neighbor.setValue(
            _to_float(
                getattr(self.params.P, "psf_flux_scale_min_neighbor_fwhm", 4.0),
                4.0,
            )
        )
        scale_form.addRow("Minimum neighbor distance (xFWHM):", self.p_flux_scale_min_neighbor)

        self.p_flux_scale_max_scatter = QDoubleSpinBox()
        self.p_flux_scale_max_scatter.setRange(0.01, 1.0)
        self.p_flux_scale_max_scatter.setSingleStep(0.01)
        self.p_flux_scale_max_scatter.setDecimals(3)
        self.p_flux_scale_max_scatter.setValue(
            _to_float(
                getattr(self.params.P, "psf_flux_scale_max_scatter_mag", 0.10),
                0.10,
            )
        )
        scale_form.addRow("Maximum reference scatter (mag):", self.p_flux_scale_max_scatter)

        self.p_fit_window_mode = QComboBox()
        self.p_fit_window_mode.addItem("Auto (PSF energy)", "auto")
        self.p_fit_window_mode.addItem("Manual (FWHM multiplier)", "manual")
        _fit_window_mode = str(
            getattr(self.params.P, "psf_fit_window_mode", "auto")
        ).strip().lower()
        _fit_window_index = self.p_fit_window_mode.findData(_fit_window_mode)
        self.p_fit_window_mode.setCurrentIndex(
            _fit_window_index if _fit_window_index >= 0 else 0
        )
        fit_form.addRow("Fit window mode:", self.p_fit_window_mode)

        self.p_fit_energy = QDoubleSpinBox()
        self.p_fit_energy.setRange(0.50, 0.995)
        self.p_fit_energy.setSingleStep(0.01)
        self.p_fit_energy.setDecimals(3)
        self.p_fit_energy.setValue(
            _to_float(
                getattr(self.params.P, "psf_fit_encircled_energy", 0.90), 0.90
            )
        )
        fit_form.addRow("Target PSF energy:", self.p_fit_energy)

        self.p_fit_mult = QDoubleSpinBox()
        self.p_fit_mult.setRange(0.5, 5.0)
        self.p_fit_mult.setSingleStep(0.1)
        self.p_fit_mult.setValue(
            _to_float(getattr(self.params.P, "psf_fit_shape_fwhm_mult", 2.4), 2.4)
        )
        fit_form.addRow("Manual fit window (xFWHM):", self.p_fit_mult)

        def _sync_fit_window_controls() -> None:
            automatic = self.p_fit_window_mode.currentData() == "auto"
            self.p_fit_energy.setEnabled(automatic)
            self.p_fit_mult.setEnabled(not automatic)

        self.p_fit_window_mode.currentIndexChanged.connect(
            lambda _index: _sync_fit_window_controls()
        )
        _sync_fit_window_controls()

        self.p_max_iter = QSpinBox()
        self.p_max_iter.setRange(1, 3)
        self.p_max_iter.setValue(_to_int(getattr(self.params.P, "psf_max_iter", 2), 2))
        fit_form.addRow("Residual passes:", self.p_max_iter)

        self.p_fitter_max_iter = QSpinBox()
        self.p_fitter_max_iter.setRange(1, 10)
        self.p_fitter_max_iter.setValue(
            _to_int(getattr(self.params.P, "psf_fitter_max_iter", 6), 6)
        )
        self.p_fitter_max_iter.setToolTip(
            "Maximum weighted Newton updates inside each residual pass"
        )
        fit_form.addRow("Fitter updates/pass:", self.p_fitter_max_iter)

        self.p_redetect = QDoubleSpinBox()
        self.p_redetect.setRange(1.0, 10.0)
        self.p_redetect.setSingleStep(0.5)
        self.p_redetect.setValue(_to_float(getattr(self.params.P, "psf_redetect_sigma", 4.0), 4.0))
        redetect_form.addRow("Re-detect sigma (base):", self.p_redetect)

        def _make_filter_sigma_spin(attr, label):
            sp = QDoubleSpinBox()
            sp.setRange(0.0, 10.0)
            sp.setSingleStep(0.5)
            sp.setDecimals(1)
            sp.setSpecialValueText("base")
            _v = _to_float(getattr(self.params.P, attr, float("nan")), float("nan"))
            sp.setValue(0.0 if not np.isfinite(_v) else float(_v))
            sp.setToolTip("0 = use base sigma")
            redetect_form.addRow(label, sp)
            return sp

        self.p_redetect_g = _make_filter_sigma_spin("psf_redetect_sigma_g", "g-band override:")
        self.p_redetect_r = _make_filter_sigma_spin("psf_redetect_sigma_r", "r-band override:")
        self.p_redetect_i = _make_filter_sigma_spin("psf_redetect_sigma_i", "i-band override:")

        self.p_dup_mult = QDoubleSpinBox()
        self.p_dup_mult.setRange(0.0, 5.0)
        self.p_dup_mult.setSingleStep(0.1)
        self.p_dup_mult.setValue(_to_float(getattr(self.params.P, "psf_duplicate_radius_fwhm_mult", 0.8), 0.8))
        redetect_form.addRow("Duplicate radius (×FWHM):", self.p_dup_mult)

        self.p_dup_px = QDoubleSpinBox()
        self.p_dup_px.setRange(0.0, 50.0)
        self.p_dup_px.setSingleStep(0.1)
        self.p_dup_px.setDecimals(2)
        _dup_px = _to_float(getattr(self.params.P, "psf_duplicate_radius_px", np.nan), np.nan)
        self.p_dup_px.setValue(0.0 if not np.isfinite(_dup_px) else float(_dup_px))
        self.p_dup_px.setToolTip("0이면 비활성(×FWHM 값 사용), >0이면 절대 px 반경 사용")
        redetect_form.addRow("Duplicate radius (px override):", self.p_dup_px)

        self.p_cap_per_iter = QSpinBox()
        self.p_cap_per_iter.setRange(0, 50000)
        self.p_cap_per_iter.setSingleStep(50)
        self.p_cap_per_iter.setValue(_to_int(getattr(self.params.P, "psf_new_sources_cap_per_iter", 70), 70))
        redetect_form.addRow("Max new/iter (abs):", self.p_cap_per_iter)

        self.p_cap_frac = QDoubleSpinBox()
        self.p_cap_frac.setRange(0.0, 1.0)
        self.p_cap_frac.setSingleStep(0.01)
        self.p_cap_frac.setValue(_to_float(getattr(self.params.P, "psf_new_sources_cap_frac", 0.02), 0.02))
        redetect_form.addRow("Max new/iter (frac):", self.p_cap_frac)

        self.p_blend_ratio = QDoubleSpinBox()
        self.p_blend_ratio.setRange(0.0, 1.0)
        self.p_blend_ratio.setDecimals(2)
        self.p_blend_ratio.setSingleStep(0.05)
        self.p_blend_ratio.setValue(
            _to_float(getattr(self.params.P, "psf_blend_residual_ratio", 0.3), 0.3)
        )
        self.p_blend_ratio.setToolTip(
            "Reject residual peaks that are too weak relative to the current source model; 0 disables"
        )
        redetect_form.addRow("Residual/model minimum:", self.p_blend_ratio)

        self.p_postfit_snr = QDoubleSpinBox()
        self.p_postfit_snr.setRange(0.0, 100.0)
        self.p_postfit_snr.setDecimals(1)
        self.p_postfit_snr.setSingleStep(0.5)
        self.p_postfit_snr.setValue(
            _to_float(getattr(self.params.P, "psf_postfit_snr_min", 3.0), 3.0)
        )
        self.p_postfit_snr.setToolTip(
            "New residual detections below this fitted S/N are removed; initial Step 4 sources are retained"
        )
        redetect_form.addRow("Post-fit S/N minimum:", self.p_postfit_snr)

        self.p_postfit_qfit = QDoubleSpinBox()
        self.p_postfit_qfit.setRange(0.0, 100.0)
        self.p_postfit_qfit.setDecimals(2)
        self.p_postfit_qfit.setSingleStep(0.1)
        self.p_postfit_qfit.setValue(
            _to_float(getattr(self.params.P, "psf_postfit_qfit_max", 3.0), 3.0)
        )
        self.p_postfit_qfit.setToolTip(
            "Remove new residual sources above qfit / expected-noise qfit; "
            "0 disables. Initial Step 4 sources are retained."
        )
        redetect_form.addRow("Post-fit qfit/noise maximum:", self.p_postfit_qfit)

        self.p_postfit_redchi = QDoubleSpinBox()
        self.p_postfit_redchi.setRange(0.0, 100000.0)
        self.p_postfit_redchi.setDecimals(1)
        self.p_postfit_redchi.setSingleStep(5.0)
        self.p_postfit_redchi.setValue(
            _to_float(
                getattr(self.params.P, "psf_postfit_reduced_chi2_max", 25.0),
                25.0,
            )
        )
        self.p_postfit_redchi.setToolTip(
            "Remove newly detected sources above this reduced chi-square; 0 disables."
        )
        redetect_form.addRow("Post-fit reduced chi2 maximum:", self.p_postfit_redchi)

        self.p_fit_init_max = QSpinBox()
        self.p_fit_init_max.setRange(0, 200000)
        self.p_fit_init_max.setSingleStep(100)
        self.p_fit_init_max.setValue(_to_int(getattr(self.params.P, "psf_fit_init_max_sources", 0), 0))
        self.p_fit_init_max.setToolTip("0이면 초기 피팅 소스 무제한")
        fit_form.addRow("Initial fit source cap (0=off):", self.p_fit_init_max)

        self.p_core_enable = QCheckBox("Hard-exclude crowded core during PSF fit")
        self.p_core_enable.setChecked(bool(getattr(self.params.P, "psf_core_cut_enable", False)))
        self.p_core_enable.setToolTip(
            "Optional CPU/quality safeguard. Leave off to fit the full field; unresolved pairs are retained with a crowding flag."
        )
        core_form.addRow("", self.p_core_enable)

        self.p_core_center_mode = QComboBox()
        self.p_core_center_mode.addItem("Auto density peak", "auto")
        self.p_core_center_mode.addItem("Image center", "image")
        self.p_core_center_mode.addItem("Manual x/y", "manual")
        _core_mode = str(getattr(self.params.P, "psf_core_cut_center_mode", "auto")).strip().lower()
        _core_mode_i = self.p_core_center_mode.findData(_core_mode)
        self.p_core_center_mode.setCurrentIndex(_core_mode_i if _core_mode_i >= 0 else 0)
        core_form.addRow("Center:", self.p_core_center_mode)

        self.p_core_x = QDoubleSpinBox()
        self.p_core_x.setRange(0.0, 200000.0)
        self.p_core_x.setDecimals(1)
        self.p_core_x.setSingleStep(10.0)
        _cx = _to_float(getattr(self.params.P, "psf_core_cut_x_px", 0.0), 0.0)
        self.p_core_x.setValue(0.0 if not np.isfinite(_cx) else float(_cx))
        core_form.addRow("Manual center x (px):", self.p_core_x)

        self.p_core_y = QDoubleSpinBox()
        self.p_core_y.setRange(0.0, 200000.0)
        self.p_core_y.setDecimals(1)
        self.p_core_y.setSingleStep(10.0)
        _cy = _to_float(getattr(self.params.P, "psf_core_cut_y_px", 0.0), 0.0)
        self.p_core_y.setValue(0.0 if not np.isfinite(_cy) else float(_cy))
        core_form.addRow("Manual center y (px):", self.p_core_y)

        self.p_core_radius_px = QDoubleSpinBox()
        self.p_core_radius_px.setRange(0.0, 200000.0)
        self.p_core_radius_px.setDecimals(1)
        self.p_core_radius_px.setSingleStep(10.0)
        self.p_core_radius_px.setSpecialValueText("auto")
        self.p_core_radius_px.setValue(_to_float(getattr(self.params.P, "psf_core_cut_radius_px", 0.0), 0.0))
        self.p_core_radius_px.setToolTip("0 = estimate radius from the detection-density profile")
        core_form.addRow("Cut radius (px):", self.p_core_radius_px)

        self.p_core_radius_mult = QDoubleSpinBox()
        self.p_core_radius_mult.setRange(1.0, 200.0)
        self.p_core_radius_mult.setDecimals(1)
        self.p_core_radius_mult.setSingleStep(1.0)
        self.p_core_radius_mult.setValue(_to_float(getattr(self.params.P, "psf_core_cut_radius_fwhm_mult", 20.0), 20.0))
        self.p_core_radius_mult.setToolTip("Fallback and safety cap for the automatic core radius")
        core_form.addRow("Auto radius cap (xFWHM):", self.p_core_radius_mult)

        self.p_core_density_ratio = QDoubleSpinBox()
        self.p_core_density_ratio.setRange(1.0, 20.0)
        self.p_core_density_ratio.setDecimals(2)
        self.p_core_density_ratio.setSingleStep(0.1)
        self.p_core_density_ratio.setValue(
            _to_float(getattr(self.params.P, "psf_core_cut_auto_min_density_ratio", 1.5), 1.5)
        )
        self.p_core_density_ratio.setToolTip("Auto center/cut is disabled when the density peak is below this contrast")
        core_form.addRow("Min density contrast:", self.p_core_density_ratio)

        self.p_substar_nei_mult = QDoubleSpinBox()
        self.p_substar_nei_mult.setRange(2.0, 30.0)
        self.p_substar_nei_mult.setSingleStep(0.5)
        self.p_substar_nei_mult.setValue(_to_float(getattr(self.params.P, "psf_substar_neighbor_r_fwhm_mult", 8.0), 8.0))
        fit_form.addRow("Substar neighbor radius (×FWHM):", self.p_substar_nei_mult)

        self.p_substar_iters = QSpinBox()
        self.p_substar_iters.setRange(0, 2)
        self.p_substar_iters.setValue(_to_int(getattr(self.params.P, "psf_substar_iters", 1), 1))
        self.p_substar_iters.setToolTip("0 disables neighbour cleaning; 1 is the recommended CPU default")
        fit_form.addRow("Substar passes:", self.p_substar_iters)

        self.p_substar_max_src = QSpinBox()
        self.p_substar_max_src.setRange(0, 200000)
        self.p_substar_max_src.setSingleStep(100)
        self.p_substar_max_src.setValue(_to_int(getattr(self.params.P, "psf_substar_max_sources", 1500), 1500))
        self.p_substar_max_src.setToolTip("0이면 substar 이웃 소스 캡 무제한")
        fit_form.addRow("Substar max neighbor sources:", self.p_substar_max_src)

        self.p_conv_new = QDoubleSpinBox()
        self.p_conv_new.setRange(0.0, 1.0)
        self.p_conv_new.setSingleStep(0.005)
        self.p_conv_new.setValue(_to_float(getattr(self.params.P, "psf_conv_new_frac", 0.02), 0.02))
        self.p_conv_new.setToolTip("Uses unique candidates before the CPU source cap")
        redetect_form.addRow("Converge candidate frac <", self.p_conv_new)

        self.p_conv_flux = QDoubleSpinBox()
        self.p_conv_flux.setRange(0.0, 1.0)
        self.p_conv_flux.setSingleStep(0.001)
        self.p_conv_flux.setValue(_to_float(getattr(self.params.P, "psf_flux_conv_threshold", 0.01), 0.01))
        fit_form.addRow("Flux convergence fraction:", self.p_conv_flux)

        self.p_use_grouper = QCheckBox("Fit close neighbours together (CPU-limited)")
        self.p_use_grouper.setChecked(bool(getattr(self.params.P, "psf_use_grouper", False)))
        fit_form.addRow("", self.p_use_grouper)

        self.p_grouper_max_size = QSpinBox()
        self.p_grouper_max_size.setRange(1, 25)
        self.p_grouper_max_size.setValue(_to_int(getattr(self.params.P, "psf_grouper_max_size", 3), 3))
        self.p_grouper_max_size.setToolTip(
            "1 disables grouping; 2-3 is the CPU default. Groups above 4 use sparse LSQR."
        )
        fit_form.addRow("Group max size:", self.p_grouper_max_size)

        self.p_grouper_radius = QDoubleSpinBox()
        self.p_grouper_radius.setRange(0.5, 5.0)
        self.p_grouper_radius.setSingleStep(0.25)
        self.p_grouper_radius.setSuffix(" FWHM")
        self.p_grouper_radius.setValue(
            _to_float(getattr(self.params.P, "psf_grouper_radius_fwhm", 1.5), 1.5)
        )
        self.p_grouper_radius.setToolTip(
            "Neighbours inside this separation are fit simultaneously. Larger values cost more CPU."
        )
        fit_form.addRow("Group radius:", self.p_grouper_radius)

        self.p_forced_match_radius = QDoubleSpinBox()
        self.p_forced_match_radius.setRange(0.1, 3.0)
        self.p_forced_match_radius.setSingleStep(0.05)
        self.p_forced_match_radius.setSuffix(" FWHM")
        self.p_forced_match_radius.setValue(
            _to_float(
                getattr(self.params.P, "psf_forced_match_radius_fwhm", 1.25),
                1.25,
            )
        )
        self.p_forced_match_radius.setToolTip(
            "Step 4 detections inside this radius are anchored to their Step 7 catalog positions."
        )
        fit_form.addRow("Forced-catalog match:", self.p_forced_match_radius)

        self.p_use_error_img = QCheckBox("Use error image (slower, higher RAM)")
        self.p_use_error_img.setChecked(bool(getattr(self.params.P, "psf_use_error_image", True)))
        fit_form.addRow("", self.p_use_error_img)

        self.p_shared_filter_epsf = QCheckBox(
            "Share EPSF per filter (faster; disable if seeing varies >1px across frames)"
        )
        self.p_shared_filter_epsf.setChecked(
            bool(getattr(self.params.P, "psf_shared_filter_epsf", False))
        )
        epsf_form.addRow("", self.p_shared_filter_epsf)

        self.p_min_epsf_stars = QSpinBox()
        self.p_min_epsf_stars.setRange(1, 200)
        self.p_min_epsf_stars.setSingleStep(1)
        self.p_min_epsf_stars.setValue(_to_int(getattr(self.params.P, "psf_min_epsf_stars", 10), 10))
        self.p_min_epsf_stars.setToolTip(
            "Min isolated PSF stars required to build/cache a new EPSF.\n"
            "With 'Share EPSF' ON: frames below this threshold reuse the cached filter EPSF.\n"
            "Raise to avoid bad EPSF from crowded/trailed frames (e.g. 10–20)."
        )
        epsf_form.addRow("Min isolated PSF stars:", self.p_min_epsf_stars)

        self.p_sharp_lo = QDoubleSpinBox()
        self.p_sharp_lo.setRange(0.0, 1.0)
        self.p_sharp_lo.setSingleStep(0.05)
        self.p_sharp_lo.setDecimals(2)
        self.p_sharp_lo.setValue(_to_float(getattr(self.params.P, "psf_redetect_sharp_lo", 0.15), 0.15))
        redetect_form.addRow("Re-detect sharpness min:", self.p_sharp_lo)

        self.p_sharp_hi = QDoubleSpinBox()
        self.p_sharp_hi.setRange(0.0, 1.0)
        self.p_sharp_hi.setSingleStep(0.05)
        self.p_sharp_hi.setDecimals(2)
        self.p_sharp_hi.setValue(_to_float(getattr(self.params.P, "psf_redetect_sharp_hi", 0.95), 0.95))
        redetect_form.addRow("Re-detect sharpness max:", self.p_sharp_hi)

        self.p_round_max = QDoubleSpinBox()
        self.p_round_max.setRange(0.0, 2.0)
        self.p_round_max.setSingleStep(0.05)
        self.p_round_max.setDecimals(2)
        self.p_round_max.setValue(_to_float(getattr(self.params.P, "psf_redetect_round_abs_max", 0.8), 0.8))
        redetect_form.addRow("Re-detect |roundness| max:", self.p_round_max)

        self.p_save_residuals = QCheckBox("Save residual FITS (required for iter viewer)")
        self.p_save_residuals.setChecked(True)
        self.p_save_residuals.setEnabled(False)
        output_form.addRow("", self.p_save_residuals)

        self.p_save_all_iter_residuals = QCheckBox(
            "Also save background-restored star-subtracted images for every pass"
        )
        self.p_save_all_iter_residuals.setChecked(
            bool(getattr(self.params.P, "psf_save_all_iter_residuals", False))
        )
        output_form.addRow("", self.p_save_all_iter_residuals)

        # ── mode logic ────────────────────────────────────────────────────
        _manual_widgets = [
            self.p_n_stars, self.p_isolation, self.p_epsf_contamination_filter,
            self.p_flux_scale_correction, self.p_flux_scale_min_snr,
            self.p_flux_scale_min_stars, self.p_flux_scale_min_neighbor,
            self.p_flux_scale_max_scatter,
            self.p_fit_window_mode, self.p_fit_energy, self.p_fit_mult,
            self.p_max_iter,
            self.p_fitter_max_iter,
            self.p_redetect, self.p_dup_mult, self.p_dup_px,
            self.p_cap_per_iter, self.p_cap_frac, self.p_blend_ratio,
            self.p_postfit_snr, self.p_postfit_qfit, self.p_postfit_redchi,
            self.p_fit_init_max,
            self.p_core_enable, self.p_core_center_mode, self.p_core_x, self.p_core_y,
            self.p_core_radius_px, self.p_core_radius_mult, self.p_core_density_ratio,
            self.p_substar_iters, self.p_substar_nei_mult, self.p_substar_max_src,
            self.p_conv_new, self.p_conv_flux, self.p_use_grouper,
            self.p_grouper_max_size, self.p_grouper_radius,
            self.p_forced_match_radius,
            self.p_sharp_lo, self.p_sharp_hi, self.p_round_max,
        ]

        def _apply_mode_to_widgets(mode_key):
            p = _PSF_MODE_PRESETS.get(mode_key, _PSF_MODE_PRESETS["normal"])
            self.p_n_stars.setValue(p["psf_n_stars_max"])
            self.p_isolation.setValue(p["psf_isolation_fwhm_mult"])
            self.p_epsf_contamination_filter.setChecked(
                bool(p["psf_epsf_contamination_filter"])
            )
            self.p_flux_scale_correction.setChecked(
                bool(p["psf_flux_scale_correction"])
            )
            self.p_fit_window_mode.setCurrentIndex(
                max(0, self.p_fit_window_mode.findData(p["psf_fit_window_mode"]))
            )
            self.p_fit_energy.setValue(p["psf_fit_encircled_energy"])
            self.p_fit_mult.setValue(p["psf_fit_shape_fwhm_mult"])
            self.p_max_iter.setValue(p["psf_max_iter"])
            self.p_fitter_max_iter.setValue(p["psf_fitter_max_iter"])
            self.p_redetect.setValue(p["psf_redetect_sigma"])
            self.p_dup_mult.setValue(p["psf_duplicate_radius_fwhm_mult"])
            self.p_dup_px.setValue(0.0)
            self.p_cap_per_iter.setValue(p["psf_new_sources_cap_per_iter"])
            self.p_cap_frac.setValue(p["psf_new_sources_cap_frac"])
            self.p_blend_ratio.setValue(p["psf_blend_residual_ratio"])
            self.p_postfit_snr.setValue(p["psf_postfit_snr_min"])
            self.p_postfit_qfit.setValue(p["psf_postfit_qfit_max"])
            self.p_postfit_redchi.setValue(p["psf_postfit_reduced_chi2_max"])
            self.p_fit_init_max.setValue(p["psf_fit_init_max_sources"])
            self.p_core_enable.setChecked(bool(p["psf_core_cut_enable"]))
            self.p_core_center_mode.setCurrentIndex(max(0, self.p_core_center_mode.findData("auto")))
            self.p_core_x.setValue(0.0)
            self.p_core_y.setValue(0.0)
            self.p_core_radius_px.setValue(p["psf_core_cut_radius_px"])
            self.p_core_radius_mult.setValue(p["psf_core_cut_radius_fwhm_mult"])
            self.p_core_density_ratio.setValue(p["psf_core_cut_auto_min_density_ratio"])
            self.p_substar_iters.setValue(p["psf_substar_iters"])
            self.p_substar_nei_mult.setValue(p["psf_substar_neighbor_r_fwhm_mult"])
            self.p_substar_max_src.setValue(p["psf_substar_max_sources"])
            self.p_conv_new.setValue(p["psf_conv_new_frac"])
            self.p_conv_flux.setValue(p["psf_flux_conv_threshold"])
            self.p_use_grouper.setChecked(p["psf_use_grouper"])
            self.p_grouper_radius.setValue(p["psf_grouper_radius_fwhm"])
            self.p_forced_match_radius.setValue(p["psf_forced_match_radius_fwhm"])
            self.p_sharp_lo.setValue(p["psf_redetect_sharp_lo"])
            self.p_sharp_hi.setValue(p["psf_redetect_sharp_hi"])
            self.p_round_max.setValue(p["psf_redetect_round_abs_max"])

        _epsf_only_widgets = [
            self.p_oversampling, self.p_shared_filter_epsf, self.p_min_epsf_stars,
        ]

        def _refresh_controls():
            engine = self.p_fit_engine.currentData()
            is_moffat  = (self.p_build_mode.currentData() == "moffat")

            for w in _manual_widgets:
                w.setEnabled(True)
            for w in _epsf_only_widgets:
                w.setEnabled(not is_moffat)
            self.p_use_error_img.setEnabled(engine == "photutils")
            self.p_grouper_max_size.setEnabled(self.p_use_grouper.isChecked())
            self.p_grouper_radius.setEnabled(self.p_use_grouper.isChecked())
            scale_enabled = self.p_flux_scale_correction.isChecked()
            for widget in (
                self.p_flux_scale_min_snr,
                self.p_flux_scale_min_stars,
                self.p_flux_scale_min_neighbor,
                self.p_flux_scale_max_scatter,
            ):
                widget.setEnabled(scale_enabled)
            _sync_fit_window_controls()

        def _on_mode_changed():
            mode_key = mode_combo.currentData()
            if mode_key != "custom":
                _apply_mode_to_widgets(mode_key)
            _refresh_controls()

        mode_combo.currentIndexChanged.connect(lambda *_: _on_mode_changed())
        self.p_fit_engine.currentIndexChanged.connect(lambda *_: _refresh_controls())
        self.p_build_mode.currentIndexChanged.connect(lambda *_: _refresh_controls())
        self.p_use_grouper.toggled.connect(lambda *_: _refresh_controls())
        self.p_flux_scale_correction.toggled.connect(lambda *_: _refresh_controls())
        _refresh_controls()

        btns = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        def _after_psf_reset():
            _refresh_controls()

        add_parameter_reset_button(
            btns,
            [
                (mode_combo, "crowded"),
                (self.p_fit_engine, "apex_iterative"),
                (self.p_build_mode, "epsf"),
                (self.p_workers, 2),
                (self.p_oversampling, 2),
                (self.p_epsf_mult, 4.0),
                (self.p_n_stars, 0),
                (self.p_isolation, 2.0),
                (self.p_epsf_contamination_filter, True),
                (self.p_flux_scale_correction, False),
                (self.p_flux_scale_min_snr, 50.0),
                (self.p_flux_scale_min_stars, 8),
                (self.p_flux_scale_min_neighbor, 4.0),
                (self.p_flux_scale_max_scatter, 0.10),
                (self.p_fit_window_mode, "auto"),
                (self.p_fit_energy, 0.90),
                (self.p_fit_mult, 2.4),
                (self.p_max_iter, 2),
                (self.p_fitter_max_iter, 8),
                (self.p_redetect, 4.5),
                (self.p_redetect_g, 4.0),
                (self.p_redetect_r, 4.0),
                (self.p_redetect_i, 4.5),
                (self.p_dup_mult, 0.4),
                (self.p_dup_px, 0.0),
                (self.p_cap_per_iter, 50),
                (self.p_cap_frac, 0.01),
                (self.p_blend_ratio, 0.3),
                (self.p_postfit_snr, 3.0),
                (self.p_postfit_qfit, 3.0),
                (self.p_postfit_redchi, 25.0),
                (self.p_fit_init_max, 3000),
                (self.p_core_enable, False),
                (self.p_core_center_mode, "auto"),
                (self.p_core_x, 0.0),
                (self.p_core_y, 0.0),
                (self.p_core_radius_px, 0.0),
                (self.p_core_radius_mult, 20.0),
                (self.p_core_density_ratio, 1.5),
                (self.p_substar_iters, 1),
                (self.p_substar_nei_mult, 5.0),
                (self.p_substar_max_src, 1000),
                (self.p_conv_new, 0.02),
                (self.p_conv_flux, 0.01),
                (self.p_use_grouper, False),
                (self.p_grouper_max_size, 3),
                (self.p_grouper_radius, 1.5),
                (self.p_forced_match_radius, 1.25),
                (self.p_use_error_img, False),
                (self.p_shared_filter_epsf, False),
                (self.p_min_epsf_stars, 10),
                (self.p_sharp_lo, 0.2),
                (self.p_sharp_hi, 0.9),
                (self.p_round_max, 0.6),
                (self.p_save_residuals, True),
                (self.p_save_all_iter_residuals, False),
            ],
            on_reset=_after_psf_reset,
        )
        btns.accepted.connect(lambda: self._save_params(dialog, mode_combo.currentData()))
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)
        dialog.exec_()


    def _save_params(self, dialog, mode_key="normal"):
        self.params.P.psf_mode = mode_key
        self.params.P.psf_model_mode = "per_frame"
        self.params.P.psf_fit_engine = self.p_fit_engine.currentData()
        self.params.P.psf_build_mode = self.p_build_mode.currentData()
        self.params.P.psf_parallel_workers = self.p_workers.value()
        self.params.P.psf_epsf_oversampling = self.p_oversampling.value()
        self.params.P.psf_epsf_size_fwhm_mult = self.p_epsf_mult.value()
        self.params.P.psf_n_stars_max = self.p_n_stars.value()
        self.params.P.psf_isolation_fwhm_mult = self.p_isolation.value()
        self.params.P.psf_epsf_contamination_filter = (
            self.p_epsf_contamination_filter.isChecked()
        )
        self.params.P.psf_flux_scale_correction = self.p_flux_scale_correction.isChecked()
        self.params.P.psf_flux_scale_min_snr = self.p_flux_scale_min_snr.value()
        self.params.P.psf_flux_scale_min_stars = self.p_flux_scale_min_stars.value()
        self.params.P.psf_flux_scale_min_neighbor_fwhm = self.p_flux_scale_min_neighbor.value()
        self.params.P.psf_flux_scale_max_scatter_mag = self.p_flux_scale_max_scatter.value()
        self.params.P.psf_fit_window_mode = self.p_fit_window_mode.currentData()
        self.params.P.psf_fit_encircled_energy = self.p_fit_energy.value()
        self.params.P.psf_fit_shape_fwhm_mult = self.p_fit_mult.value()
        self.params.P.psf_max_iter = self.p_max_iter.value()
        self.params.P.psf_fitter_max_iter = self.p_fitter_max_iter.value()
        self.params.P.psf_redetect_sigma = self.p_redetect.value()
        def _spin_to_sigma(sp):
            v = sp.value()
            return float("nan") if v <= 0.0 else v
        self.params.P.psf_redetect_sigma_g = _spin_to_sigma(self.p_redetect_g)
        self.params.P.psf_redetect_sigma_r = _spin_to_sigma(self.p_redetect_r)
        self.params.P.psf_redetect_sigma_i = _spin_to_sigma(self.p_redetect_i)
        self.params.P.psf_duplicate_radius_fwhm_mult = self.p_dup_mult.value()
        self.params.P.psf_duplicate_radius_px = self.p_dup_px.value() if self.p_dup_px.value() > 0 else np.nan
        self.params.P.psf_new_sources_cap_per_iter = self.p_cap_per_iter.value()
        self.params.P.psf_new_sources_cap_frac = self.p_cap_frac.value()
        self.params.P.psf_blend_residual_ratio = self.p_blend_ratio.value()
        self.params.P.psf_postfit_snr_min = self.p_postfit_snr.value()
        self.params.P.psf_postfit_qfit_max = self.p_postfit_qfit.value()
        self.params.P.psf_postfit_reduced_chi2_max = self.p_postfit_redchi.value()
        self.params.P.psf_fit_init_max_sources = self.p_fit_init_max.value()
        self.params.P.psf_core_cut_enable = self.p_core_enable.isChecked()
        self.params.P.psf_core_cut_center_mode = self.p_core_center_mode.currentData()
        self.params.P.psf_core_cut_x_px = self.p_core_x.value()
        self.params.P.psf_core_cut_y_px = self.p_core_y.value()
        self.params.P.psf_core_cut_radius_px = self.p_core_radius_px.value()
        self.params.P.psf_core_cut_radius_fwhm_mult = self.p_core_radius_mult.value()
        self.params.P.psf_core_cut_auto_min_density_ratio = self.p_core_density_ratio.value()
        self.params.P.psf_substar_iters = self.p_substar_iters.value()
        self.params.P.psf_substar_neighbor_r_fwhm_mult = self.p_substar_nei_mult.value()
        self.params.P.psf_substar_max_sources = self.p_substar_max_src.value()
        self.params.P.psf_conv_new_frac = self.p_conv_new.value()
        self.params.P.psf_flux_conv_threshold = self.p_conv_flux.value()
        self.params.P.psf_use_grouper = self.p_use_grouper.isChecked()
        self.params.P.psf_grouper_max_size = self.p_grouper_max_size.value()
        self.params.P.psf_grouper_radius_fwhm = self.p_grouper_radius.value()
        self.params.P.psf_forced_match_radius_fwhm = self.p_forced_match_radius.value()
        self.params.P.psf_use_error_image = self.p_use_error_img.isChecked()
        self.params.P.psf_shared_filter_epsf = self.p_shared_filter_epsf.isChecked()
        self.params.P.psf_min_epsf_stars = self.p_min_epsf_stars.value()
        self.params.P.psf_save_all_iter_residuals = self.p_save_all_iter_residuals.isChecked()
        self.params.P.psf_redetect_sharp_lo = self.p_sharp_lo.value()
        self.params.P.psf_redetect_sharp_hi = self.p_sharp_hi.value()
        self.params.P.psf_redetect_round_abs_max = self.p_round_max.value()
        self.params.P.psf_save_residuals = self.p_save_residuals.isChecked()
        self.save_state()
        self.persist_params()
        QMessageBox.information(dialog, "Saved", "Parameters saved.")
        dialog.accept()

    # ── Thread cleanup ────────────────────────────────────────────────────────

    def _cleanup_worker(self, timeout_ms=5000):
        if not self.worker:
            return
        if self.worker.isRunning():
            self.worker.stop()
            self.worker.quit()
            self.worker.wait(timeout_ms)
        try:
            self.worker.deleteLater()
        except Exception:
            pass
        self.worker = None

    # ── Log ───────────────────────────────────────────────────────────────────

    def log(self, message: str):
        append_timestamped_log(self.log_text, message)

    def show_log_window(self):
        show_raised(self.log_window)

    # ── Skip label ────────────────────────────────────────────────────────────

    def _update_skip_label(self):
        if not hasattr(self, "skip_label"):
            return
        if self._skip_psf:
            self.skip_label.setText(
                f"PSF SKIPPED — {self.downstream_name} will use Step 7 forced aperture results."
            )
        else:
            psf_idx = step8_psf_dir(self.params.P.result_dir) / "photometry_index.csv"
            if psf_idx.exists():
                self.skip_label.setText("PSF photometry results available.")
            else:
                self.skip_label.setText("")

    # ── Validation / State ────────────────────────────────────────────────────

    def validate_step(self) -> bool:
        """Step 8 is always valid: either PSF was run or it was skipped."""
        if self._skip_psf:
            return True
        valid, _ = self._current_psf_cache_status()
        return valid

    def save_state(self):
        self.project_state.store_step_data("psf_photometry", {
            "skip_psf": self._skip_psf,
            "use_existing_psf_output": (
                self.chk_use_existing_output.isChecked()
                if hasattr(self, "chk_use_existing_output")
                else True
            ),
            "psf_mode": getattr(self.params.P, "psf_mode", "normal"),
            "psf_model_mode": getattr(self.params.P, "psf_model_mode", "per_frame"),
            "psf_fit_engine": getattr(self.params.P, "psf_fit_engine", "apex_iterative"),
            "psf_build_mode": getattr(self.params.P, "psf_build_mode", "epsf"),
            "psf_parallel_workers": getattr(self.params.P, "psf_parallel_workers", 0),
            "psf_epsf_oversampling": getattr(self.params.P, "psf_epsf_oversampling", 2),
            "psf_epsf_size_px": getattr(self.params.P, "psf_epsf_size_px", 25),
            "psf_epsf_size_fwhm_mult": getattr(self.params.P, "psf_epsf_size_fwhm_mult", 4.0),
            "psf_n_stars_max": getattr(self.params.P, "psf_n_stars_max", 0),
            "psf_isolation_fwhm_mult": getattr(self.params.P, "psf_isolation_fwhm_mult", 3.0),
            "psf_epsf_contamination_filter": getattr(
                self.params.P,
                "psf_epsf_contamination_filter",
                True,
            ),
            "psf_flux_scale_correction": getattr(
                self.params.P, "psf_flux_scale_correction", False
            ),
            "psf_flux_scale_min_snr": getattr(
                self.params.P, "psf_flux_scale_min_snr", 50.0
            ),
            "psf_flux_scale_min_stars": getattr(
                self.params.P, "psf_flux_scale_min_stars", 8
            ),
            "psf_flux_scale_min_neighbor_fwhm": getattr(
                self.params.P, "psf_flux_scale_min_neighbor_fwhm", 4.0
            ),
            "psf_flux_scale_max_scatter_mag": getattr(
                self.params.P, "psf_flux_scale_max_scatter_mag", 0.10
            ),
            "psf_fit_shape_px": getattr(self.params.P, "psf_fit_shape_px", 5),
            "psf_fit_shape_fwhm_mult": getattr(self.params.P, "psf_fit_shape_fwhm_mult", 2.4),
            "psf_fit_window_mode": getattr(
                self.params.P, "psf_fit_window_mode", "auto"
            ),
            "psf_fit_encircled_energy": getattr(
                self.params.P, "psf_fit_encircled_energy", 0.90
            ),
            "psf_use_grouper": getattr(self.params.P, "psf_use_grouper", False),
            "psf_max_iter": getattr(self.params.P, "psf_max_iter", 2),
            "psf_fitter_max_iter": getattr(self.params.P, "psf_fitter_max_iter", 6),
            "psf_redetect_sigma": getattr(self.params.P, "psf_redetect_sigma", 4.0),
            "psf_redetect_sigma_g": getattr(self.params.P, "psf_redetect_sigma_g", float("nan")),
            "psf_redetect_sigma_r": getattr(self.params.P, "psf_redetect_sigma_r", float("nan")),
            "psf_redetect_sigma_i": getattr(self.params.P, "psf_redetect_sigma_i", float("nan")),
            "psf_duplicate_radius_fwhm_mult": getattr(self.params.P, "psf_duplicate_radius_fwhm_mult", 0.8),
            "psf_duplicate_radius_px": getattr(self.params.P, "psf_duplicate_radius_px", np.nan),
            "psf_new_sources_cap_per_iter": getattr(self.params.P, "psf_new_sources_cap_per_iter", 70),
            "psf_new_sources_cap_frac": getattr(self.params.P, "psf_new_sources_cap_frac", 0.02),
            "psf_blend_residual_ratio": getattr(self.params.P, "psf_blend_residual_ratio", 0.3),
            "psf_postfit_snr_min": getattr(self.params.P, "psf_postfit_snr_min", 3.0),
            "psf_postfit_qfit_max": getattr(self.params.P, "psf_postfit_qfit_max", 3.0),
            "psf_postfit_reduced_chi2_max": getattr(
                self.params.P, "psf_postfit_reduced_chi2_max", 25.0
            ),
            "psf_fit_init_max_sources": getattr(self.params.P, "psf_fit_init_max_sources", 0),
            "psf_core_cut_enable": getattr(self.params.P, "psf_core_cut_enable", False),
            "psf_core_cut_center_mode": getattr(self.params.P, "psf_core_cut_center_mode", "auto"),
            "psf_core_cut_x_px": getattr(self.params.P, "psf_core_cut_x_px", 0.0),
            "psf_core_cut_y_px": getattr(self.params.P, "psf_core_cut_y_px", 0.0),
            "psf_core_cut_radius_px": getattr(self.params.P, "psf_core_cut_radius_px", 0.0),
            "psf_core_cut_radius_fwhm_mult": getattr(self.params.P, "psf_core_cut_radius_fwhm_mult", 20.0),
            "psf_core_cut_auto_min_density_ratio": getattr(self.params.P, "psf_core_cut_auto_min_density_ratio", 1.5),
            "psf_substar_iters": getattr(self.params.P, "psf_substar_iters", 1),
            "psf_substar_neighbor_r_fwhm_mult": getattr(self.params.P, "psf_substar_neighbor_r_fwhm_mult", 8.0),
            "psf_substar_max_sources": getattr(self.params.P, "psf_substar_max_sources", 1500),
            "psf_conv_new_frac": getattr(self.params.P, "psf_conv_new_frac", 0.02),
            "psf_flux_conv_threshold": getattr(self.params.P, "psf_flux_conv_threshold", 0.01),
            "psf_use_error_image": getattr(self.params.P, "psf_use_error_image", False),
            "psf_shared_filter_epsf": getattr(self.params.P, "psf_shared_filter_epsf", False),
            "psf_grouper_max_size": getattr(self.params.P, "psf_grouper_max_size", 3),
            "psf_grouper_radius_fwhm": getattr(
                self.params.P, "psf_grouper_radius_fwhm", 1.5
            ),
            "psf_forced_match_radius_fwhm": getattr(
                self.params.P, "psf_forced_match_radius_fwhm", 1.25
            ),
            "psf_min_epsf_stars": getattr(self.params.P, "psf_min_epsf_stars", 10),
            "psf_save_all_iter_residuals": getattr(self.params.P, "psf_save_all_iter_residuals", False),
            "psf_redetect_sharp_lo": getattr(self.params.P, "psf_redetect_sharp_lo", 0.15),
            "psf_redetect_sharp_hi": getattr(self.params.P, "psf_redetect_sharp_hi", 0.95),
            "psf_redetect_round_abs_max": getattr(self.params.P, "psf_redetect_round_abs_max", 0.8),
            "psf_save_residuals": getattr(self.params.P, "psf_save_residuals", True),
        })

    def restore_state(self):
        state = self.project_state.get_step_data("psf_photometry")
        if state:
            self._skip_psf = bool(state.get("skip_psf", False))
            if hasattr(self, "chk_use_existing_output"):
                self.chk_use_existing_output.setChecked(
                    bool(state.get("use_existing_psf_output", True))
                )
            for k, v in state.items():
                if k not in {"skip_psf", "use_existing_psf_output"} and hasattr(self.params.P, k):
                    setattr(self.params.P, k, v)
        if str(getattr(self.params.P, "psf_model_mode", "per_frame")).strip().lower() != "per_frame":
            self.params.P.psf_model_mode = "per_frame"
        if str(getattr(self.params.P, "psf_fit_engine", "apex_iterative")).strip().lower() == "allstar":
            self.params.P.psf_fit_engine = "apex_iterative"
        _mode = str(getattr(self.params.P, "psf_mode", "normal")).strip().lower()
        _star_cap = _to_int(getattr(self.params.P, "psf_n_stars_max", 0), 0)
        if _mode != "custom" and _star_cap in {30, 40, 50}:
            self.params.P.psf_n_stars_max = 0
        self.params.P.psf_grouper_max_size = min(
            25,
            max(1, _to_int(getattr(self.params.P, "psf_grouper_max_size", 3), 3)),
        )
        self.params.P.psf_grouper_radius_fwhm = min(
            5.0,
            max(
                0.5,
                _to_float(getattr(self.params.P, "psf_grouper_radius_fwhm", 1.5), 1.5),
            ),
        )
        self.params.P.psf_forced_match_radius_fwhm = min(
            3.0,
            max(
                0.1,
                _to_float(
                    getattr(self.params.P, "psf_forced_match_radius_fwhm", 1.25),
                    1.25,
                ),
            ),
        )
        # Clamp fit_shape_fwhm_mult to a sensible minimum (< 1.0 is unusable).
        _fmult = _to_float(getattr(self.params.P, "psf_fit_shape_fwhm_mult", 1.5), 1.5)
        if _fmult < 1.0:
            self.params.P.psf_fit_shape_fwhm_mult = 1.5
        # Migrate broad defaults to tuned defaults unless user explicitly changed them.
        _rsig = _to_float(getattr(self.params.P, "psf_redetect_sigma", 4.0), 4.0)
        if abs(_rsig - 6.0) < 1e-6 or abs(_rsig - 7.5) < 1e-6:
            # 6.0 and 7.5 were old defaults; migrate to current default.
            self.params.P.psf_redetect_sigma = 4.0
        _cap_abs = _to_int(getattr(self.params.P, "psf_new_sources_cap_per_iter", 70), 70)
        if _cap_abs == 100:
            self.params.P.psf_new_sources_cap_per_iter = 70
        _cap_frac = _to_float(getattr(self.params.P, "psf_new_sources_cap_frac", 0.02), 0.02)
        if abs(_cap_frac - 0.04) < 1e-6:
            self.params.P.psf_new_sources_cap_frac = 0.02
        _slo = _to_float(getattr(self.params.P, "psf_redetect_sharp_lo", 0.15), 0.15)
        _shi = _to_float(getattr(self.params.P, "psf_redetect_sharp_hi", 0.95), 0.95)
        _rnd = _to_float(getattr(self.params.P, "psf_redetect_round_abs_max", 0.8), 0.8)
        if _slo <= -900.0 and _shi >= 900.0:
            self.params.P.psf_redetect_sharp_lo = 0.15
            self.params.P.psf_redetect_sharp_hi = 0.95
        if _rnd >= 9.0:
            self.params.P.psf_redetect_round_abs_max = 0.8
        # Migrate overly-loose UI values.
        if _slo <= 0.01 and _shi >= 0.99 and _rnd >= 1.5:
            self.params.P.psf_redetect_sharp_lo = 0.15
            self.params.P.psf_redetect_sharp_hi = 0.95
            self.params.P.psf_redetect_round_abs_max = 0.8
        self._update_skip_label()
        self.update_frame_table()
        idx_path = step8_psf_dir(self.params.P.result_dir) / "photometry_index.csv"
        if (not self._skip_psf) and idx_path.exists():
            valid, reason = self._current_psf_cache_status()
            if valid:
                try:
                    idx = pd.read_csv(idx_path)
                    n_frames = len(idx)
                    self.progress_bar.setMaximum(max(1, n_frames))
                    self.progress_bar.setValue(n_frames)
                    self.progress_label.setText(f"Loaded previous PSF output ({n_frames} frames)")
                    self.log(f"[PSF][CACHE] Loaded previous Step 8 PSF output from disk ({n_frames} frames).")
                    self._load_from_disk()
                except Exception as exc:
                    self.log(f"[PSF][CACHE] Previous output could not be loaded: {exc}")
            else:
                self.log(f"[PSF][CACHE] Previous output not restored ({reason}).")

    def closeEvent(self, event):
        if self.worker and self.worker.isRunning():
            self.stop_psf()
        self._cleanup_worker(timeout_ms=10000)
        super().closeEvent(event)
