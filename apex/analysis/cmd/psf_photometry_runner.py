"""Step 8 PSF photometry — the computation, with no Qt in it.

Lifted out of `apex/gui/workflow/cmd/step8_psf_photometry.py`. The headless
pipeline used to import that GUI module for this worker, so PSF photometry was
impossible on a Qt-free install: the step reported NOT_IMPLEMENTED and stopped.
The worker never called a single `QThread` method — `QThread` was only the base
class, and the eight `pyqtSignal`s only a way to report progress.

The body is unchanged apart from where it announces: `self.progress.emit(...)`
became `self.on_progress.send(...)`, and the GUI subclass subscribes its Qt
signals to those channels. The window and the script therefore run the same
object, which is what makes "the app and the batch pipeline agree" a fact about
identity rather than a claim to be tested.
"""

from __future__ import annotations

from apex.analysis.psf_diagnostics import draw_psf_final_diagnostics, load_psf_final_diagnostic_data
from apex.analysis.psf_flux_scale import PSFApertureScale, apply_psf_aperture_scale, estimate_psf_aperture_scale
from apex.analysis.psf_iteration import IterationSnapshot, PSFFitFlag, assess_psf_frame_quality, decide_residual_iteration, fit_parameters_changed, measure_psf_fit_quality, qfit_noise_diagnostics
from apex.analysis.psf_policy import estimate_psf_flux_seeds, local_group_policy, merge_forced_catalog_seeds, plan_epsf_stars, plan_psf_fit_window, psf_symmetric_mask, select_epsf_reference_stars, select_spatially_balanced
from apex.utils.astro_utils import normalize_filter_name
from apex.utils.constants import get_parallel_workers
from apex.utils.noise_params import resolve_effective_noise_params
from apex.utils.psf_core import PSFCoreCut, estimate_psf_core_cut, psf_core_keep_mask, target_pixel_from_wcs
from apex.utils.qc_utils import filter_files_by_qc, should_use_frame_quality_qc
from apex.utils.step_paths import step7_forced_phot_dir
from apex.utils.step_paths_cmd import step2_cropped_dir, step4_dir, step8_psf_dir, crop_is_active
from astropy.io import fits
from astropy.nddata import NDData
from astropy.stats import sigma_clipped_stats, mad_std as _mad_std
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from dataclasses import replace
from matplotlib.figure import Figure
from pathlib import Path
from scipy.spatial import cKDTree
from threading import Lock
import copy
import hashlib
import json
import numpy as np
import pandas as pd
import threading
import time
import traceback
from apex.analysis.worker_signals import ReportsProgress


def _fast_res_std(arr: np.ndarray) -> float:
    """Robust std for residual images: MAD estimator on a 65K-pixel subsample."""
    flat = arr.ravel()
    flat = flat[np.isfinite(flat)]
    if flat.size == 0:
        return 0.0
    stride = max(1, flat.size // 65536)
    return float(_mad_std(flat[::stride]))

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

    # 구경 vs PSF 비교 그림. 창에서만 나오던 것이라 헤드리스 실행에는
    # 없었다 — 두 측광이 서로 맞는지 말해주는 바로 그 그림이다.
    if cmp_df is not None and not getattr(cmp_df, "empty", True):
        try:
            df_cut, n_before = filter_ap_vs_psf(cmp_df)
            fig = Figure(figsize=(10.0, 4.0), dpi=120)
            draw_ap_vs_psf(fig, df_cut, n_before)
            fp = psf_dir / "step8_ap_vs_psf_comparison.png"
            fig.savefig(fp, dpi=160, bbox_inches="tight")
            saved.append(fp)
        except Exception:  # noqa: BLE001 - a figure must not fail the export
            pass

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

def _clone_psf_model(model):
    """Return a per-frame copy of PSF model to avoid thread-shared mutation."""
    try:
        return model.copy()
    except Exception:
        try:
            return copy.deepcopy(model)
        except Exception:
            return model

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

_NEWTON_GRAD_DELTA = 0.5  # sub-pixel step for numerical PSF gradient (pixels)

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

def _fit_variance(counts: np.ndarray, background_rms: float, gain: float,
                  profile_error_frac: float = 0.0) -> np.ndarray:
    """Per-pixel variance for the fit weights.

    Background plus photon noise is the whole story only if the PSF model is
    exact. DAOPHOT does not assume that: its `proferr` (5 % by default) adds a
    term proportional to the model itself, which dominates wherever the model
    is bright. The consequence is a different fit, not a different error bar —
    with pure photon weighting the core carries the most weight because it has
    the most signal; with a profile-error term the core is deliberately
    distrusted, because that is exactly where a slightly wrong PSF is most
    wrong, and the fit leans on the wings instead.

    Measured on this frame's numbers (2026-08-14): for a star peaking at
    20,000 ADU the central pixel carries 66x more weight under photon-only
    weighting than under DAOPHOT's, and the fraction of total weight inside one
    FWHM is 11.8 % versus 2.0 %. Faint stars are barely affected (1.26x). In a
    tight blend the core is where the two stars overlap, so this is where the
    two engines' answers diverge most.

    `profile_error_frac` of 0 reproduces the previous behaviour exactly.
    """
    var = (max(float(background_rms), 1e-6) ** 2
           + np.clip(counts, 0.0, None) / max(float(gain), 1e-6))
    frac = float(profile_error_frac)
    if frac > 0.0:
        var = var + (frac * np.clip(counts, 0.0, None)) ** 2
    return var

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

def _allstar_fit(img_sub: np.ndarray, positions: np.ndarray, fluxes: np.ndarray,
                 eval_psf, fit_shape: int, stamp_size: int,
                 max_iter: int, flux_conv: float,
                 max_shift: float = 2.5,
                 group_radius: float = 0.0,
                 max_group_size: int = 1,
                 max_grouped_sources: int = 0,
                 background_rms: float = 1.0,
                 gain: float = 1.0,
                 profile_error_frac: float = 0.0,
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
                variance = _fit_variance(fit_raw, background_rms, gain,
                                         profile_error_frac)
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
                variance = _fit_variance(fit_raw, background_rms, gain,
                                         profile_error_frac)
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

class PsfPhotometryRunner(ReportsProgress):
    """Per-frame EPSFBuilder + PSFPhotometry worker.

    Algorithm per frame:
    1. Load detected positions from detect_{fname}.csv
    2. Select bright isolated stars for EPSF building
    3. Build oversampled EPSF with EPSFBuilder
    4. Run iterative PSFPhotometry: fit → residual → detect new → re-fit
    5. Save residual FITS and epsf_model FITS
    6. Emit per-frame result
    """

    _CHANNELS = ('progress', 'worker_status', 'frame_done', 'epsf_ready', 'residual_ready', 'finished', 'error', 'log')

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
        self.on_log.send(msg)

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
                self.on_error.send("IMPORT", f"photutils required: {e}")
                self.on_finished.send({})
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
            # DAOPHOT's proferr, in fraction rather than percent. 0 keeps the
            # photon-only weighting APEX has always used; 0.05 is IRAF's default.
            profile_error_frac = max(0.0, min(0.5, _to_float(
                getattr(P, "psf_profile_error_frac", 0.0), 0.0)))
            # Ceiling for the final fixed-position pass; 2 was the hard-coded
            # value. ALLSTAR's equivalent is 50.
            final_pass_max_iter = max(
                1, _to_int(getattr(P, "psf_final_pass_max_iter", 2), 2)
            )
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
            # The ceiling exists because a group is solved as one linear system,
            # so cost grows faster than the star count; 25 was chosen for a
            # laptop. It is a resource limit, not a science one, and clamping a
            # requested 60 down to 25 silently made a grouping experiment
            # unrunnable (2026-08-14). Keep a bound, but a generous one, and say
            # so when the request is cut.
            _requested_group_size = grouper_max_size
            grouper_max_size = min(200, max(1, grouper_max_size))
            if grouper_max_size != _requested_group_size:
                self._log(
                    f"[APEX] grouper max_size {_requested_group_size} -> "
                    f"{grouper_max_size} (bound 1-200)"
                )
            # How much of a frame may be solved as joint groups. 0.10 keeps the
            # old default; 1.0 is the ALLSTAR-like setting where every star is
            # eligible. Cost grows with it, so it stays a user dial.
            grouper_budget_frac = min(
                1.0, max(0.0, _to_float(getattr(P, "psf_grouper_budget_frac", 0.10), 0.10))
            )
            grouper_budget_cap = max(
                0, _to_int(getattr(P, "psf_grouper_budget_cap", 200), 200)
            )
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
                    self.on_worker_status.send(0, fname, "Stopped", 100)
                    return {"file": fname, "status": "stopped"}

                wid = int(threading.get_ident() % 10000)
                _t_frame = time.time()
                _t: dict[str, float] = {"start": _t_frame}
                self.on_progress.send(completed[0], total, f"RUN | {fname}")
                self.on_worker_status.send(wid, fname, "Load", 5)
                img_path = self._resolve_fits_path(fname)
                if img_path is None:
                    self.on_worker_status.send(wid, fname, "No FITS", 100)
                    return {"file": fname, "status": "no_fits", "reason": "no FITS"}

                det_df = _load_detect_positions(fname, self.cache_dir, self.result_dir)
                if det_df is None or len(det_df) == 0:
                    self.on_worker_status.send(wid, fname, "No detect", 100)
                    return {"file": fname, "status": "no_detect", "reason": "no detect csv"}

                try:
                    img = fits.getdata(img_path).astype(np.float32, copy=False)
                except Exception as e:
                    self.on_worker_status.send(wid, fname, "FITS error", 100)
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
                    self.on_worker_status.send(wid, fname, "Background", 20)
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
                        self.on_worker_status.send(wid, fname, "Stopped", 100)
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
                        self.on_worker_status.send(wid, fname, "PSF build", 40)
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
                                    self.on_log.send(
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
                            self.on_worker_status.send(wid, fname, "Stopped", 100)
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
                                self.on_log.send(
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
                                # Detected catalog stars, as context for the merge.
                                # Without them a forced star inside the match
                                # radius of a *different, detected* star claims
                                # that star's detection and the snap deletes its
                                # seed — the tight-blend under-subtraction found
                                # 2026-08-14. They join the one-to-one matching
                                # so each detected star defends its own
                                # detection, but they are never snapped or added.
                                _ctx_ok = (
                                    (~_forced_like)
                                    & _phot_ok
                                    & _seed_x.notna()
                                    & _seed_y.notna()
                                    & (_seed_x >= _edge)
                                    & (_seed_x < (w - _edge))
                                    & (_seed_y >= _edge)
                                    & (_seed_y < (h - _edge))
                                )
                                _ctx_xy = np.column_stack([
                                    _seed_x.loc[_ctx_ok].to_numpy(dtype=float, copy=False),
                                    _seed_y.loc[_ctx_ok].to_numpy(dtype=float, copy=False),
                                ]) if _ctx_ok.any() else np.zeros((0, 2), dtype=float)
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
                                        context_xy=_ctx_xy,
                                    )
                                    xy_det = _merge.xy
                                    det_uids = _merge.det_uids
                                    init_forced_positions = _merge.forced_mask
                                    flux_init_map.update(_merge.flux_by_uid)
                                    self._log(
                                        "  [INIT] Step7 forced catalog | "
                                        f"matched={_merge.n_matched} added={_merge.n_added} "
                                        f"ctx={len(_ctx_xy)} "
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
                    self.on_progress.send(completed[0], total, f"FIT | {fname}")

                    _psf_evaluator = None  # set in APEX branch; used for final residual rendering
                    _engine_snapshots: list[dict] = []
                    _engine_stop_reason = ""

                    if psf_fit_engine_cfg == 'apex_iterative':
                        self.on_worker_status.send(wid, fname, "APEX iterative fit", 70)
                        _psf_evaluator = _make_psf_evaluator(
                            epsf_model, psf_type_built, oversampling, psf_interp_order)

                        # Outer loop: ALLSTAR fit → residual → re-detect → add → repeat
                        # max_iter controls number of find+fit cycles (DAOPHOT style)
                        _INNER_ITERS = fitter_max_iter
                        _max_grp_size, _group_budget = local_group_policy(
                            len(init_params),
                            enabled=use_grouper,
                            requested_max_size=grouper_max_size,
                            # Both bounds used to be written here as literals, so
                            # a configured max_size of 60 silently became 25 and
                            # only a tenth of the frame was ever solved jointly.
                            # ALLSTAR groups every star it fits; leaving the
                            # fraction unreachable made that comparison
                            # impossible to set up (2026-08-14).
                            hard_max_size=grouper_max_size,
                            max_fraction=grouper_budget_frac,
                            absolute_cap=grouper_budget_cap,
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
                        # Locking catalog positions is right for an isolated
                        # star — the catalog knows better than one frame. In a
                        # blend it is not symmetric: the free neighbour carries
                        # an extra degree of freedom the locked star does not,
                        # so it absorbs whatever the pair's model cannot
                        # explain. Measured on the M13 blends (2026-08-14): the
                        # pair's summed flux is right (1.063) while the locked
                        # star gets 1.036 — 2.6 %, or 25 mmag, handed to the
                        # neighbour; and the neighbour's drift direction swings
                        # the locked star by +18 / -30 mmag either way.
                        # `psf_forced_position_lock` releases the lock so both
                        # members of a blend are fitted on equal terms.
                        _forced_lock = str(
                            getattr(P, "psf_forced_position_lock", "always")
                        ).strip().lower()
                        if _forced_lock not in ("always", "never"):
                            self._log(
                                f"  [APEX] unknown psf_forced_position_lock "
                                f"'{_forced_lock}'; using 'always'"
                            )
                            _forced_lock = "always"
                        if _forced_lock != "always":
                            self._log(
                                f"  [APEX] forced position lock: {_forced_lock} "
                                f"({int(np.sum(_cur_forced))} sources released)"
                            )
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
                                profile_error_frac=profile_error_frac,
                                gain=float(GAIN),
                                initial_positions=_cur_anchor,
                                initial_fit_valid=_cur_fit_valid,
                                position_bound=float(fit_shape_frame // 2),
                                position_fixed_mask=(
                                    _cur_forced if _forced_lock == "always"
                                    else np.zeros(len(_cur_forced), dtype=bool)
                                ),
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
                                    # `_cur_fit_valid` is the one array here not
                                    # rebuilt from the fit result each pass, so a
                                    # generation that adds sources can leave it a
                                    # step behind. Locked catalog positions hid
                                    # this: freeing them (2026-08-14) made the
                                    # residual pass accept enough new sources to
                                    # expose it as an IndexError. Rebuild from the
                                    # arrays that are authoritative rather than
                                    # index with a stale length.
                                    if len(_cur_fit_valid) == len(_keep_fit):
                                        _cur_fit_valid = _cur_fit_valid[_keep_fit]
                                    else:
                                        self._log(
                                            "  [APEX] fit_valid out of step "
                                            f"({len(_cur_fit_valid)} vs {len(_keep_fit)}); rebuilt"
                                        )
                                        _cur_fit_valid = (
                                            np.isfinite(_cur_xy[:, 0])
                                            & np.isfinite(_cur_xy[:, 1])
                                            & np.isfinite(_cur_flux)
                                            & ((_cur_flux > 0) | _cur_forced)
                                            & ((_cur_position_flags & _fit_invalid_mask) == 0)
                                        )

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
                        # This pass sets every published flux, and it was capped at two
                        # Newton steps regardless of the configured iteration count. An
                        # isolated star converges in one; a blended group is a coupled
                        # system that does not, and ALLSTAR gives the same solve up to
                        # 50. That is the shape of the measured deficit — APEX wins on
                        # the most isolated stars (0.033 vs 0.041) and loses on the most
                        # blended ones (2026-08-14). `flux_conv` still exits early, so a
                        # larger ceiling costs nothing where two steps were enough.
                        _final_started = time.perf_counter()
                        phot_result = _allstar_fit(
                            img_sub,
                            _cur_xy,
                            _cur_flux,
                            _psf_evaluator,
                            fit_shape=fit_shape_frame,
                            stamp_size=render_shape_frame,
                            max_iter=final_pass_max_iter,
                            flux_conv=flux_conv_threshold,
                            max_shift=float(fit_shape_frame // 2),
                            group_radius=_group_radius,
                            max_group_size=_max_grp_size,
                            max_grouped_sources=_group_budget,
                            background_rms=float(bkg_std),
                            profile_error_frac=profile_error_frac,
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
                        self.on_worker_status.send(wid, fname, "PSF fit", 70)
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
                        self.on_worker_status.send(wid, fname, "Fit failed", 100)
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
                    self.on_worker_status.send(wid, fname, "Save", 95)
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
                    self.on_worker_status.send(wid, fname, "Done", 100)
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
                    self.on_worker_status.send(wid, fname, "Error", 100)
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
                                self.on_progress.send(completed[0], total, fname_c)
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
                            self.on_progress.send(n_done, total, f"RUN={n_running} QUEUE={n_queued} ETA~{eta_txt}")
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
                                self.on_progress.send(n_done, total, f"RUN={n_running} QUEUE={n_queued} ETA~{eta_txt}")
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
                            self.on_progress.send(completed[0], total, fname)
                            self._log(f"[{completed[0]}/{total}] STOP {fname} | cancelled")
                            continue

                        try:
                            out = fut.result()
                        except Exception as e:
                            out = {"file": fname, "status": "error", "reason": str(e)}

                        completed[0] += 1
                        self.on_progress.send(completed[0], total, out.get("file", fname))
                        status = out.get("status", "error")

                        if status == "processed":
                            idx_row = out.get("idx_row", {})
                            if idx_row:
                                index_rows.append(idx_row)
                                self.on_frame_done.send(out["file"], idx_row)
                            if out.get("epsf_key") and out.get("epsf_arr") is not None:
                                self.on_epsf_ready.send(out["epsf_key"], out.get("epsf_frame", out["file"]), out["epsf_arr"])
                            self.on_residual_ready.send(out["file"], out.get("residual_meta"), out.get("new_xy"))
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
            self.on_finished.send({"frames": total, **counters})

        except Exception as e:
            self.on_error.send("PSF_WORKER", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")
            self.on_finished.send({})



# ── Aperture vs PSF comparison figure ────────────────────────────────────────
#
# This lived in the Step 8 window as a 180-line method wired to three spin
# boxes, two combo boxes, a canvas and a status label. The consequence was that
# `step8_ap_vs_psf_comparison.png` only existed if somebody had opened the
# window and looked at the tab — a headless run produced every other Step 8
# product and not this one, which is the figure that says whether the two
# photometries agree.
#
# The window keeps the widgets; this keeps the drawing. Both go through here, so
# the batch figure and the on-screen figure cannot drift apart.

_AP_VS_PSF_COLORS = {
    "u": "#9467bd", "g": "#2ca02c", "r": "#d62728",
    "i": "#ff7f0e", "z": "#8c564b", "b": "#1f77b4",
    "v": "#bcbd22", "ha": "#e377c2",
}


def filter_ap_vs_psf(merged, *, flags0_only: bool = False, snr_min: float = 0.0,
                     qfit_max: float = 0.0, dmag_clip: float = 0.0,
                     filter_name: str = "", frame_name: str = ""):
    """The cuts, separated from the drawing so both callers apply the same ones.

    Returns (df, n_before). `df` carries a `delta` column (mag_ap − mag_psf) and
    is empty when nothing survives — the caller decides what to say about that.
    """
    import numpy as np
    import pandas as pd

    if merged is None or getattr(merged, "empty", True):
        return pd.DataFrame(), 0

    df = merged.copy()
    df["mag_ap"] = pd.to_numeric(df.get("mag_ap"), errors="coerce")
    df["mag_psf"] = pd.to_numeric(df.get("mag_psf"), errors="coerce")
    df = df[np.isfinite(df["mag_ap"]) & np.isfinite(df["mag_psf"])].copy()
    if df.empty:
        return df, 0

    df["delta"] = df["mag_ap"] - df["mag_psf"]
    n_before = int(len(df))

    for col in ("flags_psf", "snr_psf", "qfit", "qfit_noise_ratio"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    if filter_name and filter_name.lower() != "all" and "FILTER" in df.columns:
        from apex.utils.astro_utils import normalize_filter_name
        key = normalize_filter_name(filter_name)
        df = df[df["FILTER"].astype(str).map(normalize_filter_name) == key].copy()
    if frame_name and frame_name.lower() != "all" and "FRAME" in df.columns:
        df = df[df["FRAME"].astype(str) == frame_name].copy()

    if flags0_only and "flags_psf" in df.columns:
        df = df[np.isfinite(df["flags_psf"]) & (df["flags_psf"] == 0)].copy()
    if snr_min > 0 and "snr_psf" in df.columns:
        df = df[np.isfinite(df["snr_psf"]) & (df["snr_psf"] >= snr_min)].copy()
    if qfit_max > 0:
        # `qfit_noise_ratio` when it exists and is populated — it is the
        # normalised version, so the same threshold means the same thing across
        # frames of different depth.
        col = ("qfit_noise_ratio"
               if "qfit_noise_ratio" in df.columns
               and np.isfinite(df["qfit_noise_ratio"]).any() else "qfit")
        if col in df.columns:
            df = df[np.isfinite(df[col]) & (df[col] <= qfit_max)].copy()
    if dmag_clip > 0:
        df = df[np.isfinite(df["delta"]) & (np.abs(df["delta"]) <= dmag_clip)].copy()

    return df, n_before


def draw_ap_vs_psf(fig, df, n_before: int, *, split_excluded: int = 0) -> str:
    """Draw the two panels onto `fig`. Returns the one-line statistics text.

    An empty `df` gets a figure that says so rather than an empty axis — a blank
    panel and "no data survived the cuts" look identical on screen and only one
    of them is informative.
    """
    import numpy as np

    fig.clf()
    if df is None or df.empty:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No data after filters.\nRelax SNR/qfit/|Δmag| settings.",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=10, color="gray")
        ax.set_xticks([]); ax.set_yticks([])
        return "No data after filters."

    filt_col = "FILTER" if "FILTER" in df.columns else None
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    stats_parts = []
    groups = df.groupby(filt_col) if filt_col else [("all", df)]
    for filt, sub in groups:
        color = _AP_VS_PSF_COLORS.get(str(filt).lower(), "#999999")
        ax1.scatter(sub["mag_ap"], sub["mag_psf"], s=4, alpha=0.35,
                    color=color, label=str(filt), rasterized=True)
        ax2.scatter(sub["mag_ap"], sub["delta"], s=4, alpha=0.35,
                    color=color, label=str(filt), rasterized=True)
        stats_parts.append(
            f"{filt}: N={len(sub)}  Δmed={float(np.nanmedian(sub['delta'])):+.3f}"
            f"  σ={float(np.nanstd(sub['delta'])):.3f}")

    all_mag = np.concatenate([df["mag_ap"].values, df["mag_psf"].values])
    lo, hi = np.nanmin(all_mag) - 0.2, np.nanmax(all_mag) + 0.2
    ax1.plot([lo, hi], [lo, hi], "k--", lw=0.8, alpha=0.5, zorder=0)
    ax1.set_xlim(lo, hi); ax1.set_ylim(lo, hi)
    ax1.set_xlabel("mag_ap", fontsize=9)
    ax1.set_ylabel("mag_psf", fontsize=9)
    ax1.set_title("Aperture vs PSF magnitude", fontsize=9)
    ax1.legend(fontsize=7, markerscale=2, loc="upper left")

    dmed_all = float(np.nanmedian(df["delta"]))
    ax2.axhline(dmed_all, color="#D62728", lw=2.0, ls="-", alpha=0.95, zorder=1,
                label=f"Δmag median {dmed_all:+.3f}")
    ax2.axhline(0.0, color="k", lw=0.8, ls="--", alpha=0.6, zorder=0)
    ax2.axhline(+0.05, color="gray", lw=0.5, ls=":", alpha=0.5, zorder=0)
    ax2.axhline(-0.05, color="gray", lw=0.5, ls=":", alpha=0.5, zorder=0)
    ax2.set_xlabel("mag_ap", fontsize=9)
    ax2.set_ylabel("Δmag  (Ap − PSF)", fontsize=9)
    ax2.set_title("Δmag vs mag_ap", fontsize=9)
    ax2.legend(fontsize=7, markerscale=2, loc="upper left")

    fig.tight_layout()
    return (f"N={len(df)}/{n_before}  |  split_excluded={int(split_excluded)}  |  "
            + "  |  ".join(stats_parts))
