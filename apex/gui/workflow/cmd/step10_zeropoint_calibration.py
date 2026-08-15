"""
Step 10: Zeropoint & Standardization
"""

from __future__ import annotations

import hashlib
import json
import time
import traceback
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import cKDTree

from apex.utils.constants import MAG_ERR_COEFF, MAD_TO_SIGMA

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib as mpl

from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QGroupBox, QMessageBox,
    QTextEdit, QDialog, QFormLayout, QDialogButtonBox, QDoubleSpinBox,
    QSpinBox, QCheckBox, QComboBox, QWidget, QTabWidget, QFileDialog, QLineEdit
)

from apex.gui.workflow.step_window_base import StepWindowBase
from apex.analysis.light_curve.photometry_source_service import (
    resolve_lightcurve_photometry_source,
)
from apex.gui.theme import mono_note_style


def _set_label_role(label, prop: str, value) -> None:
    """Swap a theme role property at runtime and repolish so the QSS re-runs.

    Roles are static in QSS; a label whose meaning changes (ROI set vs unset)
    must clear the other property or both selectors stay live.
    """
    for name in ("role", "status", "banner"):
        label.setProperty(name, value if name == prop else None)
    style = label.style()
    style.unpolish(label)
    style.polish(label)


from apex.gui.workflow.run_control import RunControlBar
from apex.gui.workflow.log_panel import WorkflowLogWindow, WorkerStatusPanel, append_timestamped_log, show_raised
from apex.gui.workflow.ui_helpers import (
    add_parameter_reset_button,
    build_scroll_param_dialog,
    create_collapsible_section,
    create_parameter_button,
    install_parameter_wheel_guard,
)
from apex.utils.astro_utils import normalize_filter_name
from apex.utils.common_helpers import format_cmd_title, photometric_system_label
from apex.utils.cmd_gaia_enrichment import (
    load_master_table as _load_master_table,
    merge_gaia_columns_from_catalog as _merge_gaia_columns_from_catalog,
)
from apex.utils.step_paths import (
    step2_cropped_dir,
    crop_is_active,
    step7_forced_phot_dir,
    step5_wcs_dir,
    tool_extinction_dir,
)
from apex.utils.step_paths_cmd import step8_psf_dir, step9_selection_dir, step10_zp_dir
from apex.utils.io_utils import parse_int64_series, read_ecsv_int64_source_id
from apex.utils.gaia_quality import (
    gaia_corrected_excess_factor,
    gaia_cstar_sigma,
    gaia_quality_mask,
    gaia_quality_report,
)
from apex.utils.qc_utils import filter_frame_df_by_qc, should_use_frame_quality_qc
from apex.utils.photometry_provenance import (
    build_photometry_provenance,
    collapse_provenance_values,
    format_photometry_provenance,
    summarize_photometry_table,
)


from apex.utils.gaia_transforms import (
    GAIA_TO_BAND       as _GAIA_TO_BAND,
    FILTER_COLOR_PREF  as _FILTER_COLOR_PREF,
    BAND_ALIASES       as _BAND_ALIASES,
    build_color_pairs  as _build_color_pairs,
    teff_from_color    as _teff_from_color,
    TEFF_COLOR_ANCHORS as _TEFF_COLOR_ANCHORS,
    filter_bands_from_columns as _filter_bands_from_columns,
)

_ZP_SIGNATURE_FILE = "zeropoint_signature.json"
_ZP_SIGNATURE_VERSION = 3
_ZP_SIGNATURE_PARAMS = (
    "match_tol_px",
    "min_master_gaia_matches",
    "cmd_snr_calib_min",
    "frame_zp_min_n",
    "cmd_apply_extinction",
    "cmd_extinction_mode",
    "zp_clip_sigma",
    "zp_fit_iters",
    "zp_slope_absmax",
    "gaia_snr_calib_min",
    "gaia_gi_min",
    "gaia_gi_max",
    "gaia_zp_slope_absmax",
    "gaia_color_slope_absmax",
    "min_snr_for_mag",
    "phot_ref_apcorr_min_keep",
    "phot_ref_require_apcorr_candidate",
    "phot_use_qc_pass_only",
    "ref_frame",
    "site_lat_deg",
    "site_lon_deg",
    "site_alt_m",
    "site_tz_offset_hours",
)


def _cmd_photometry_index_candidates(result_dir: Path | str) -> list[Path]:
    root = Path(result_dir)
    aperture_dir = step7_forced_phot_dir(root)
    candidates = [
        aperture_dir / "photometry_index.csv",
        root / "photometry_index.csv",
        root / "phot_index.csv",
        root / "phot" / "phot_index.csv",
        root / "phot" / "photometry_index.csv",
    ]
    if aperture_dir.exists():
        candidates.extend(sorted(aperture_dir.glob("*phot*index*.csv")))
    candidates.extend(sorted(root.glob("*phot*index*.csv")))
    if (root / "phot").exists():
        candidates.extend(sorted((root / "phot").glob("*phot*index*.csv")))
    return candidates


def resolve_cmd_photometry_provenance(result_dir: Path | str) -> dict[str, str]:
    """Describe the last CMD product, or the input that Step 10 would use."""
    root = Path(result_dir)
    output_dir = step10_zp_dir(root)
    for name in (
        "median_by_ID_filter_wide_cmd.csv",
        "median_by_ID_filter_wide.csv",
        "median_by_ID_filter_wide_raw.csv",
    ):
        path = output_dir / name
        if not path.exists() or path.stat().st_size <= 0:
            continue
        try:
            info = summarize_photometry_table(pd.read_csv(path, nrows=500))
            if info["source"] != "unknown":
                return info
        except Exception:
            continue

    source_info = resolve_lightcurve_photometry_source(root)
    index_path = Path(source_info["index_path"])
    if index_path.exists() and index_path.stat().st_size > 0:
        return build_photometry_provenance(
            source_info.get("source"),
            source_info.get("mag_column"),
            source_info.get("mag_error_column"),
        )

    index_path = next(
        (
            path
            for path in _cmd_photometry_index_candidates(root)
            if path.exists() and path.stat().st_size > 0
        ),
        None,
    )
    return build_photometry_provenance("aperture")


def resolve_cmd_photometry_input(
    result_dir: Path | str,
    project_state=None,
) -> dict:
    """Select complete PSF output unless Step 8 was skipped or is stale."""
    return resolve_lightcurve_photometry_source(result_dir, project_state)


def _zp_numeric(df: pd.DataFrame, column: str) -> pd.Series:
    if not isinstance(df, pd.DataFrame) or column not in df.columns:
        return pd.Series(dtype=float)
    return pd.to_numeric(df[column], errors="coerce")


def _zp_finite(df: pd.DataFrame, column: str) -> np.ndarray:
    values = _zp_numeric(df, column).to_numpy(dtype=float)
    return values[np.isfinite(values)]


def _zp_median(df: pd.DataFrame, column: str) -> float:
    values = _zp_finite(df, column)
    return float(np.median(values)) if values.size else np.nan


def _zp_mad_sigma(values) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.nan
    med = float(np.median(arr))
    return float(MAD_TO_SIGMA * np.median(np.abs(arr - med)))


def _zp_filter_values(*dfs: pd.DataFrame) -> list[str]:
    out: set[str] = set()
    for df in dfs:
        if not isinstance(df, pd.DataFrame) or "filter" not in df.columns:
            continue
        for val in df["filter"].dropna():
            filt = normalize_filter_name(val)
            if filt:
                out.add(filt)
    return sorted(out)


def _zp_filter_subset(df: pd.DataFrame, filt: str) -> pd.DataFrame:
    if not isinstance(df, pd.DataFrame) or df.empty or "filter" not in df.columns:
        return pd.DataFrame(columns=getattr(df, "columns", []))
    keys = df["filter"].map(normalize_filter_name)
    return df[keys == normalize_filter_name(filt)].copy()


# Minimum sigma-clipped calibrators before the quadratic color term is
# attempted; below this the linear fit is kept (a poorly constrained curvature
# does more harm than the ±0.02-0.03 mag it corrects on rich fields).
_QUAD_MIN_CALIBRATORS = 60


def robust_weighted_polyfit(
    x, y, w=None, degree: int = 2, clip_sigma: float = 3.0, iters: int = 5, min_n: int = 10
):
    """Sigma-clipping weighted polynomial fit.

    Returns (coeffs highest-power-first, n_inlier, mad_scatter) or
    (None, 0, nan) when the fit is not possible. Mirrors the clipping scheme
    of ``_robust_linfit`` so the quadratic refinement rejects the same kind of
    outliers the linear fit does.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if w is not None:
        w = np.asarray(w, dtype=float)
        m &= np.isfinite(w) & (w > 0)
    if int(m.sum()) < max(min_n, degree + 2):
        return None, 0, float("nan")
    coeffs = None
    scatter = float("nan")
    for _ in range(max(1, int(iters))):
        try:
            coeffs = np.polyfit(x[m], y[m], degree, w=np.sqrt(w[m]) if w is not None else None)
        except (np.linalg.LinAlgError, ValueError):
            return None, 0, float("nan")
        r = y - np.polyval(coeffs, x)
        med = float(np.nanmedian(r[m]))
        scatter = float(MAD_TO_SIGMA * np.nanmedian(np.abs(r[m] - med)))
        m_new = m & (np.abs(r - med) <= clip_sigma * max(scatter, 1e-6))
        if int(m_new.sum()) < max(min_n, degree + 2) or m_new.sum() == m.sum():
            break
        m = m_new
    return coeffs, int(m.sum()), scatter


def solve_standard_colors(
    inst_mags: dict[str, np.ndarray],
    fit_params: dict[str, dict],
    iters: int = 6,
) -> dict[str, np.ndarray]:
    """Solve the *standard* color indices from instrumental magnitudes alone.

    The transformation model per filter f is ``std_f = inst_f + zp_f +
    ct_f * C(color_col_f)`` where ``C`` is a **standard** color (difference of
    two standard magnitudes). Substituting the model into each color gives a
    small linear system across the filter chain (e.g. B,V use B-V while R uses
    V-R); it is solved by fixed-point iteration, which contracts by a factor
    ~max|ct| per pass (|ct| <~ 0.15 in practice, so ``iters=6`` converges to
    <1e-5 mag). No external catalog is used at application time — the only
    inputs are the star's own instrumental magnitudes and the already-fitted
    constants, so faint-star catalog systematics cannot leak in, and (unlike
    feeding the raw instrumental color through the color term) the applied
    color sits on the same scale the coefficients were fitted on.

    Parameters
    ----------
    inst_mags : mapping of filter name -> instrumental magnitude array.
    fit_params : mapping of filter name -> {"zp", "ct", "color_col"} as built
        by the ZP fit ("color_col" like ``"B_V"`` or ``"none"``). Optional
        ``"color_min"``/``"color_max"`` give the color range the coefficients
        were fitted over; when present the applied color is clamped to it.
    iters : fixed-point iterations.

    Returns
    -------
    dict of color-column name (e.g. ``"B_V"``) -> standard-color array.
    Stars lacking a needed instrumental magnitude get NaN for that color, as do
    stars whose iteration did not converge.

    Convergence
    -----------
    The contraction argument above holds for the *linear* model. Once the
    quadratic term was added the per-pass factor became ``|ct + 2*ct2*C|``,
    which depends on the star, so bounding the coefficients is not enough: a
    star far outside the calibrator color range can leave the basin and run
    away. M13 B (ct = +0.307, ct2 = -0.170, so the basin is
    ``-2.1 < C < 3.8``) produced exactly that — one star started at C = -1.82
    and reached C = -70 in six passes, for ``mag_std_B = -844``.

    Two guards prevent it. Clamping C to the fitted color range makes the map
    bounded, so the iteration cannot diverge and the quadratic is never
    extrapolated beyond the calibrators that constrain it. Because callers may
    supply no range, the iteration is also checked for convergence and
    non-converged stars are returned as NaN rather than as a large number.
    """
    pairs: dict[str, tuple[str, str]] = {}
    for fp in fit_params.values():
        name = str(fp.get("color_col", "none"))
        if name != "none" and "_" in name and name not in pairs:
            fa, fb = name.split("_", 1)
            # Both bands must have instrumental mags AND fitted coefficients;
            # otherwise the pair cannot be placed on the standard scale and the
            # caller must keep its legacy instrumental color.
            if fa in inst_mags and fb in inst_mags and fa in fit_params and fb in fit_params:
                pairs[name] = (fa, fb)
    if not pairs:
        return {}

    def _coef(f: str, key: str) -> float:
        try:
            v = float(fit_params[f][key])
            return v if np.isfinite(v) else 0.0
        except (KeyError, TypeError, ValueError):
            return 0.0

    # Initial guess: instrumental color + zeropoint difference (exact when all
    # color terms are zero).
    colors = {
        name: (
            np.asarray(inst_mags[fa], dtype=float)
            - np.asarray(inst_mags[fb], dtype=float)
            + (_coef(fa, "zp") - _coef(fb, "zp"))
        )
        for name, (fa, fb) in pairs.items()
    }

    def _bounds(f: str) -> tuple[float, float]:
        """Color range the filter's coefficients were fitted over."""
        lo = _coef(f, "color_min") if "color_min" in fit_params.get(f, {}) else -np.inf
        hi = _coef(f, "color_max") if "color_max" in fit_params.get(f, {}) else np.inf
        return (lo, hi) if lo < hi else (-np.inf, np.inf)

    def _std_mag(f: str) -> np.ndarray:
        base = np.asarray(inst_mags[f], dtype=float) + _coef(f, "zp")
        ccol = str(fit_params.get(f, {}).get("color_col", "none"))
        if ccol in colors:
            # Applying the color term outside the calibrator color range both
            # extrapolates the quadratic and can break the fixed point.
            c = np.clip(colors[ccol], *_bounds(f))
            # Optional quadratic color term (ct2 defaults to 0 for legacy fits).
            return base + _coef(f, "ct") * c + _coef(f, "ct2") * c * c
        return base

    previous = dict(colors)
    for _ in range(max(1, int(iters))):
        previous = colors
        colors = {name: _std_mag(fa) - _std_mag(fb) for name, (fa, fb) in pairs.items()}

    # A converged star moves by ~1e-5 mag on the final pass; anything still
    # moving by more than a milli-magnitude is diverging and has no calibration.
    for name in colors:
        step = np.abs(colors[name] - previous[name])
        colors[name] = np.where(
            np.isfinite(colors[name]) & np.isfinite(step) & (step <= 1e-3),
            colors[name], np.nan,
        )
    return colors


def _first_text(df: pd.DataFrame, column: str) -> str:
    if not isinstance(df, pd.DataFrame) or df.empty or column not in df.columns:
        return ""
    vals = df[column].dropna().astype(str)
    return vals.iloc[0] if len(vals) else ""


def build_zp_qc_summary(
    coeff_df: pd.DataFrame,
    frame_df: pd.DataFrame | None = None,
    cut_df: pd.DataFrame | None = None,
    reject_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build a compact per-filter Step 10 calibration QC table."""
    coeff_df = coeff_df.copy() if isinstance(coeff_df, pd.DataFrame) else pd.DataFrame()
    frame_df = frame_df.copy() if isinstance(frame_df, pd.DataFrame) else pd.DataFrame()
    cut_df = cut_df.copy() if isinstance(cut_df, pd.DataFrame) else pd.DataFrame()
    reject_df = reject_df.copy() if isinstance(reject_df, pd.DataFrame) else pd.DataFrame()

    filters = _zp_filter_values(coeff_df, frame_df, cut_df, reject_df)
    rows = []
    for filt in filters:
        coeff = _zp_filter_subset(coeff_df, filt)
        frame = _zp_filter_subset(frame_df, filt)
        cuts = _zp_filter_subset(cut_df, filt)
        rejects = _zp_filter_subset(reject_df, filt)

        zp_vals = _zp_finite(frame, "zp_frame")
        rows.append({
            "filter": filt,
            "global_zp": _zp_median(coeff, "zp"),
            "color_term": _zp_median(coeff, "ct"),
            "fit_scatter_rms": _zp_median(coeff, "scatter_rms"),
            "n_fit_calibrators": int(_zp_median(coeff, "N")) if np.isfinite(_zp_median(coeff, "N")) else 0,
            "color_col": _first_text(coeff, "color_col"),
            "ref_source": _first_text(coeff, "ref_source"),
            "n_frame_zp": int(len(frame)),
            "frame_zp_median": float(np.median(zp_vals)) if zp_vals.size else np.nan,
            "frame_zp_sigma_mad": _zp_mad_sigma(zp_vals),
            "median_frame_zp_scatter": _zp_median(frame, "zp_scatter"),
            "median_n_ref_per_frame": _zp_median(frame, "n_ref"),
            "median_outlier_fraction": _zp_median(frame, "outlier_fraction"),
            "median_snr_ref": _zp_median(frame, "snr_med"),
            "n_rejected_frames": int(len(rejects)),
            "n_total_measurements": int(_zp_median(cuts, "n_total")) if np.isfinite(_zp_median(cuts, "n_total")) else 0,
            "n_kept_measurements": int(_zp_median(cuts, "n_kept")) if np.isfinite(_zp_median(cuts, "n_kept")) else 0,
        })

    summary = pd.DataFrame(rows)
    if not summary.empty:
        total = pd.to_numeric(summary["n_total_measurements"], errors="coerce").replace(0, np.nan)
        kept = pd.to_numeric(summary["n_kept_measurements"], errors="coerce")
        summary["kept_measurement_fraction"] = kept / total
    return summary


def draw_zp_qc_overview(
    fig: Figure,
    coeff_df: pd.DataFrame,
    frame_df: pd.DataFrame,
    summary_df: pd.DataFrame,
) -> bool:
    """Draw the static Step 10 paper/audit QC overview."""
    if not isinstance(summary_df, pd.DataFrame) or summary_df.empty:
        return False
    fig.clear()
    axes = fig.subplots(2, 2, squeeze=False)
    ax_frame, ax_scatter = axes[0]
    ax_nref, ax_kept = axes[1]

    frame_df = frame_df.copy() if isinstance(frame_df, pd.DataFrame) else pd.DataFrame()
    coeff_df = coeff_df.copy() if isinstance(coeff_df, pd.DataFrame) else pd.DataFrame()
    filters = summary_df["filter"].astype(str).tolist()
    cmap = mpl.colormaps.get_cmap("tab10")
    color_map = {filt: cmap(i % 10) for i, filt in enumerate(filters)}

    has_frame = False
    if not frame_df.empty and {"filter", "zp_frame"} <= set(frame_df.columns):
        for filt in filters:
            sub = _zp_filter_subset(frame_df, filt).reset_index(drop=True)
            vals = _zp_numeric(sub, "zp_frame").to_numpy(dtype=float)
            ok = np.isfinite(vals)
            if not np.any(ok):
                continue
            has_frame = True
            x = np.arange(len(sub), dtype=float)
            ax_frame.plot(x[ok], vals[ok], "o-", ms=3.0, lw=0.9, color=color_map[filt], label=filt)
            med = float(np.median(vals[ok]))
            sig = _zp_mad_sigma(vals[ok])
            ax_frame.axhline(med, color=color_map[filt], lw=0.8, ls="--", alpha=0.7)
            if np.isfinite(sig) and sig > 0:
                ax_frame.fill_between(
                    [np.nanmin(x[ok]), np.nanmax(x[ok])],
                    med - sig,
                    med + sig,
                    color=color_map[filt],
                    alpha=0.08,
                )
    ax_frame.set_title("Per-Frame Zeropoint")
    ax_frame.set_xlabel("frame index")
    ax_frame.set_ylabel("ZP (mag)")
    ax_frame.grid(True, alpha=0.25)
    if has_frame:
        ax_frame.legend(fontsize=8, frameon=False)
    else:
        ax_frame.text(0.5, 0.5, "No frame ZP data", ha="center", va="center", transform=ax_frame.transAxes)

    x = np.arange(len(filters), dtype=float)
    fit_scatter = pd.to_numeric(summary_df.get("fit_scatter_rms"), errors="coerce").to_numpy(dtype=float)
    frame_sigma = pd.to_numeric(summary_df.get("frame_zp_sigma_mad"), errors="coerce").to_numpy(dtype=float)
    width = 0.36
    ax_scatter.bar(x - width / 2, fit_scatter, width=width, color="#4C78A8", label="fit residual")
    ax_scatter.bar(x + width / 2, frame_sigma, width=width, color="#F58518", label="frame ZP sigma")
    ax_scatter.set_xticks(x)
    ax_scatter.set_xticklabels(filters)
    ax_scatter.set_ylabel("mag")
    ax_scatter.set_title("Calibration Scatter")
    ax_scatter.grid(True, axis="y", alpha=0.25)
    ax_scatter.legend(fontsize=8, frameon=False)

    n_fit = pd.to_numeric(summary_df.get("n_fit_calibrators"), errors="coerce").fillna(0).to_numpy(dtype=float)
    n_frame = pd.to_numeric(summary_df.get("n_frame_zp"), errors="coerce").fillna(0).to_numpy(dtype=float)
    n_rej = pd.to_numeric(summary_df.get("n_rejected_frames"), errors="coerce").fillna(0).to_numpy(dtype=float)
    ax_nref.bar(x - width / 2, n_fit, width=width, color="#54A24B", label="fit stars")
    ax_nref.bar(x + width / 2, n_frame, width=width, color="#B279A2", label="frame ZPs")
    if np.any(n_rej > 0):
        ax_nref.scatter(x + width / 2, n_frame + n_rej, marker="x", color="#D62728", label="rejected frames")
    ax_nref.set_xticks(x)
    ax_nref.set_xticklabels(filters)
    ax_nref.set_ylabel("count")
    ax_nref.set_title("Calibration Sample Size")
    ax_nref.grid(True, axis="y", alpha=0.25)
    ax_nref.legend(fontsize=8, frameon=False)

    kept = pd.to_numeric(summary_df.get("kept_measurement_fraction"), errors="coerce").to_numpy(dtype=float)
    med_nref = pd.to_numeric(summary_df.get("median_n_ref_per_frame"), errors="coerce").to_numpy(dtype=float)
    ax_kept.bar(x, kept * 100.0, color="#E45756", alpha=0.75, label="kept measurements")
    ax_kept.set_ylim(0, 105)
    ax_kept.set_xticks(x)
    ax_kept.set_xticklabels(filters)
    ax_kept.set_ylabel("kept (%)")
    ax_kept.set_title("Frame-ZP Reference Cuts")
    ax_kept.grid(True, axis="y", alpha=0.25)
    ax2 = ax_kept.twinx()
    ax2.plot(x, med_nref, "o-", color="#2F4B7C", lw=1.2, ms=4.0, label="median n_ref/frame")
    ax2.set_ylabel("median n_ref/frame")
    ax2.grid(False)

    handles1, labels1 = ax_kept.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax_kept.legend(handles1 + handles2, labels1 + labels2, fontsize=8, frameon=False, loc="best")

    fig.suptitle("Step 10 Zeropoint Calibration QC", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return True


def export_zp_qc_products(output_dir: Path, log_func=None) -> list[Path]:
    """Export Step 10 QC summary and overview figure from existing calibration outputs."""
    output_dir = Path(output_dir)
    coeff_path = output_dir / "zp_fit_coefficients.csv"
    if not coeff_path.exists() or coeff_path.stat().st_size == 0:
        return []
    try:
        coeff_df = pd.read_csv(coeff_path)
    except Exception:
        return []

    def _read_optional(name: str) -> pd.DataFrame:
        path = output_dir / name
        if not path.exists() or path.stat().st_size == 0:
            return pd.DataFrame()
        try:
            return pd.read_csv(path)
        except Exception:
            return pd.DataFrame()

    frame_df = _read_optional("frame_zeropoint.csv")
    cut_df = _read_optional("frame_zeropoint_cut_summary.csv")
    reject_df = _read_optional("frame_zeropoint_rejects.csv")

    summary = build_zp_qc_summary(coeff_df, frame_df, cut_df, reject_df)
    if summary.empty:
        return []

    saved: list[Path] = []
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "zp_qc_summary.csv"
    summary.to_csv(summary_path, index=False)
    saved.append(summary_path)

    fig = Figure(figsize=(11.0, 7.6), dpi=120)
    if draw_zp_qc_overview(fig, coeff_df, frame_df, summary):
        fig_path = output_dir / "step10_zp_qc_overview.png"
        fig.savefig(fig_path, dpi=160, bbox_inches="tight")
        saved.append(fig_path)

    if saved and log_func is not None:
        try:
            log_func("Step10 QC products exported: " + ", ".join(path.name for path in saved))
        except Exception:
            pass
    return saved


def _cmd_mag_col(system: str, band: str) -> str:
    return f"mag_{system}_{band}"


def _cmd_err_col(system: str, band: str) -> str:
    return f"mag_{system}_err_{band}"


def select_cmd_qc_axes(df: pd.DataFrame, system: str = "std") -> dict | None:
    """Choose the calibrated CMD axes with the most finite stars."""
    if not isinstance(df, pd.DataFrame) or df.empty:
        return None
    bands = _filter_bands_from_columns(df.columns, f"mag_{system}_")
    pairs = _build_color_pairs(bands, adjacent_only=True)
    best = None
    for a, b in pairs:
        ca = _cmd_mag_col(system, a)
        cb = _cmd_mag_col(system, b)
        cy = _cmd_mag_col(system, b)
        if ca not in df.columns or cb not in df.columns or cy not in df.columns:
            continue
        ma = pd.to_numeric(df[ca], errors="coerce").to_numpy(dtype=float)
        mb = pd.to_numeric(df[cb], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(df[cy], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(ma) & np.isfinite(mb) & np.isfinite(y)
        n = int(np.sum(mask))
        if best is None or n > int(best["n"]):
            best = {
                "system": system,
                "color_a": a,
                "color_b": b,
                "y_band": b,
                "n": n,
            }
    return best


def build_cmd_qc_summary(df: pd.DataFrame, system: str = "std") -> pd.DataFrame:
    """Summarize final wide CMD finite counts and magnitude-error behavior."""
    if not isinstance(df, pd.DataFrame) or df.empty:
        return pd.DataFrame()
    bands = _filter_bands_from_columns(df.columns, f"mag_{system}_")
    axes = select_cmd_qc_axes(df, system=system)
    rows = []
    for band in bands:
        mag_col = _cmd_mag_col(system, band)
        err_col = _cmd_err_col(system, band)
        mag = pd.to_numeric(df.get(mag_col, pd.Series(dtype=float)), errors="coerce").to_numpy(dtype=float)
        err = pd.to_numeric(df.get(err_col, pd.Series(dtype=float)), errors="coerce").to_numpy(dtype=float)
        snr = pd.to_numeric(df.get(f"snr_{band}", pd.Series(dtype=float)), errors="coerce").to_numpy(dtype=float)
        finite_mag = mag[np.isfinite(mag)]
        finite_err = err[np.isfinite(err) & (err >= 0)]
        finite_snr = snr[np.isfinite(snr)]
        rows.append({
            "filter": band,
            "system": system,
            "n_sources_total": int(len(df)),
            "n_finite_mag": int(finite_mag.size),
            "finite_mag_fraction": float(finite_mag.size / max(len(df), 1)),
            "mag_p05": float(np.nanpercentile(finite_mag, 5)) if finite_mag.size else np.nan,
            "mag_p50": float(np.nanpercentile(finite_mag, 50)) if finite_mag.size else np.nan,
            "mag_p95": float(np.nanpercentile(finite_mag, 95)) if finite_mag.size else np.nan,
            "median_mag_err": float(np.nanmedian(finite_err)) if finite_err.size else np.nan,
            "p90_mag_err": float(np.nanpercentile(finite_err, 90)) if finite_err.size else np.nan,
            "median_snr": float(np.nanmedian(finite_snr)) if finite_snr.size else np.nan,
            "cmd_color_a": axes.get("color_a", "") if axes else "",
            "cmd_color_b": axes.get("color_b", "") if axes else "",
            "cmd_y_band": axes.get("y_band", "") if axes else "",
            "cmd_n_points": int(axes.get("n", 0)) if axes else 0,
        })
    return pd.DataFrame(rows)


def _cmd_error_curve(df: pd.DataFrame, system: str, band: str, max_bins: int = 12) -> pd.DataFrame:
    mag_col = _cmd_mag_col(system, band)
    err_col = _cmd_err_col(system, band)
    if mag_col not in df.columns or err_col not in df.columns:
        return pd.DataFrame()
    mag = pd.to_numeric(df[mag_col], errors="coerce").to_numpy(dtype=float)
    err = pd.to_numeric(df[err_col], errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(mag) & np.isfinite(err) & (err >= 0)
    mag = mag[mask]
    err = err[mask]
    if len(mag) < 6:
        return pd.DataFrame()
    order = np.argsort(mag)
    mag = mag[order]
    err = err[order]
    n_bins = max(3, min(int(max_bins), max(3, len(mag) // 8)))
    chunks = np.array_split(np.arange(len(mag)), n_bins)
    rows = []
    for chunk in chunks:
        if len(chunk) == 0:
            continue
        rows.append({
            "filter": band,
            "mag": float(np.nanmedian(mag[chunk])),
            "mag_err": float(np.nanmedian(err[chunk])),
            "n": int(len(chunk)),
        })
    return pd.DataFrame(rows)


def draw_cmd_qc_overview(fig: Figure, df: pd.DataFrame, summary_df: pd.DataFrame, system: str = "std") -> bool:
    """Draw a final CMD plus photometric-error overview."""
    axes = select_cmd_qc_axes(df, system=system)
    if axes is None or int(axes.get("n", 0)) <= 0:
        return False

    a = axes["color_a"]
    b = axes["color_b"]
    y_band = axes["y_band"]
    ca = _cmd_mag_col(system, a)
    cb = _cmd_mag_col(system, b)
    cy = _cmd_mag_col(system, y_band)
    ma = pd.to_numeric(df[ca], errors="coerce").to_numpy(dtype=float)
    mb = pd.to_numeric(df[cb], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[cy], errors="coerce").to_numpy(dtype=float)
    color = ma - mb
    mask = np.isfinite(color) & np.isfinite(y)
    if not np.any(mask):
        return False

    fig.clear()
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0], hspace=0.32, wspace=0.28)
    ax_cmd = fig.add_subplot(gs[:, 0])
    ax_err = fig.add_subplot(gs[0, 1])
    ax_count = fig.add_subplot(gs[1, 1])

    sc = ax_cmd.scatter(
        color[mask],
        y[mask],
        s=5,
        c=y[mask],
        cmap="viridis_r",
        alpha=0.65,
        linewidths=0,
        rasterized=True,
    )
    ax_cmd.set_xlabel(f"{a} - {b}")
    ax_cmd.set_ylabel(f"{y_band} ({system})")
    ax_cmd.set_title(f"Final CMD | N={int(np.sum(mask))}")
    ax_cmd.invert_yaxis()
    ax_cmd.grid(True, alpha=0.18)
    fig.colorbar(sc, ax=ax_cmd, label=f"{y_band} mag", fraction=0.046, pad=0.04)

    bands = summary_df["filter"].astype(str).tolist() if isinstance(summary_df, pd.DataFrame) and not summary_df.empty else []
    cmap = mpl.colormaps.get_cmap("tab10")
    for i, band in enumerate(bands[:6]):
        curve = _cmd_error_curve(df, system, band)
        if curve.empty:
            continue
        ax_err.plot(
            curve["mag"].to_numpy(float),
            curve["mag_err"].to_numpy(float),
            "o-",
            ms=3.0,
            lw=1.0,
            color=cmap(i % 10),
            label=band,
        )
    ax_err.set_xlabel("magnitude")
    ax_err.set_ylabel("median mag err")
    ax_err.set_title("Photometric Error Curve")
    ax_err.grid(True, alpha=0.25)
    if ax_err.get_lines():
        ax_err.legend(fontsize=8, frameon=False)

    if isinstance(summary_df, pd.DataFrame) and not summary_df.empty:
        x = np.arange(len(summary_df), dtype=float)
        counts = pd.to_numeric(summary_df["n_finite_mag"], errors="coerce").fillna(0).to_numpy(dtype=float)
        frac = pd.to_numeric(summary_df["finite_mag_fraction"], errors="coerce").fillna(0).to_numpy(dtype=float)
        ax_count.bar(x, counts, color="#4C78A8", alpha=0.78)
        ax_count.set_xticks(x)
        ax_count.set_xticklabels(summary_df["filter"].astype(str).tolist())
        ax_count.set_ylabel("finite mag count")
        ax_count.set_title("Final Catalog Coverage")
        ax_count.grid(True, axis="y", alpha=0.25)
        ax2 = ax_count.twinx()
        ax2.plot(x, frac * 100.0, "o-", color="#E45756", lw=1.2, ms=4.0)
        ax2.set_ylabel("finite (%)")
        ax2.set_ylim(0, 105)

    fig.suptitle("Step 10 Final CMD QC", fontsize=12)
    fig.subplots_adjust(left=0.07, right=0.93, bottom=0.08, top=0.92, wspace=0.32, hspace=0.34)
    return True


def export_cmd_qc_products(output_dir: Path, log_func=None) -> list[Path]:
    """Export final calibrated CMD QC products from Step 10 wide output."""
    output_dir = Path(output_dir)
    cmd_path = output_dir / "median_by_ID_filter_wide_cmd.csv"
    if not cmd_path.exists() or cmd_path.stat().st_size == 0:
        return []
    try:
        df = pd.read_csv(cmd_path)
    except Exception:
        return []

    summary = build_cmd_qc_summary(df, system="std")
    if summary.empty:
        return []

    saved: list[Path] = []
    summary_path = output_dir / "cmd_qc_summary.csv"
    summary.to_csv(summary_path, index=False)
    saved.append(summary_path)

    fig = Figure(figsize=(11.2, 7.2), dpi=120)
    if draw_cmd_qc_overview(fig, df, summary, system="std"):
        fig_path = output_dir / "step10_cmd_qc_overview.png"
        fig.savefig(fig_path, dpi=160, bbox_inches="tight")
        saved.append(fig_path)

    if saved and log_func is not None:
        try:
            log_func("Step10 CMD QC products exported: " + ", ".join(path.name for path in saved))
        except Exception:
            pass
    return saved


def _gaia_observed_color(df: pd.DataFrame) -> np.ndarray:
    if "gaia_BP_RP" in df.columns:
        return pd.to_numeric(df["gaia_BP_RP"], errors="coerce").to_numpy(dtype=float)
    if {"gaia_BP", "gaia_RP"} <= set(df.columns):
        bp = pd.to_numeric(df["gaia_BP"], errors="coerce").to_numpy(dtype=float)
        rp = pd.to_numeric(df["gaia_RP"], errors="coerce").to_numpy(dtype=float)
        return bp - rp
    return np.full(len(df), np.nan, dtype=float)


def _poly_eval_array(x, coeffs) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.zeros_like(x, dtype=float)
    xp = np.ones_like(x, dtype=float)
    for c in coeffs:
        y += float(c) * xp
        xp *= x
    return y


def _gaia_reference_mag_for_band(df: pd.DataFrame, band: str) -> np.ndarray:
    """Transform Gaia observed G/BP/RP into the requested native filter."""
    if "gaia_G" not in df.columns:
        return np.full(len(df), np.nan, dtype=float)
    key = _BAND_ALIASES.get(band, band)
    if key not in _GAIA_TO_BAND:
        return np.full(len(df), np.nan, dtype=float)
    coeffs, lo, hi, _, _ = _GAIA_TO_BAND[key]
    g = pd.to_numeric(df["gaia_G"], errors="coerce").to_numpy(dtype=float)
    color = _gaia_observed_color(df)
    out = np.full(len(df), np.nan, dtype=float)
    ok = np.isfinite(g) & np.isfinite(color) & (color >= float(lo)) & (color <= float(hi))
    if np.any(ok):
        out[ok] = g[ok] - _poly_eval_array(color[ok], coeffs)
    return out


def _synthetic_gaia_cmd_arrays(df: pd.DataFrame) -> dict | None:
    need = {"gaia_G", "gaia_G_syn", "gaia_BP_RP_syn"}
    if not need <= set(df.columns):
        return None
    gaia_mag = pd.to_numeric(df["gaia_G"], errors="coerce").to_numpy(dtype=float)
    gaia_color = _gaia_observed_color(df)
    apex_mag = pd.to_numeric(df["gaia_G_syn"], errors="coerce").to_numpy(dtype=float)
    apex_color = pd.to_numeric(df["gaia_BP_RP_syn"], errors="coerce").to_numpy(dtype=float)
    matched = np.isfinite(gaia_mag) & np.isfinite(gaia_color) & np.isfinite(apex_mag) & np.isfinite(apex_color)
    if not np.any(matched):
        return None
    snr_bands: list[str] = []
    if {"mag_std_g", "mag_std_i"} <= set(df.columns):
        snr_bands = ["g", "i"]
    elif {"mag_std_V", "mag_std_I"} <= set(df.columns):
        snr_bands = ["V", "I"]
    return {
        "comparison": "apex_synthetic_gaia_minus_gaia_observed",
        "mode": "synthetic_gaia",
        "color_label": "BP - RP",
        "apex_color_label": "BP - RP synthetic",
        "mag_label": "G",
        "apex_mag_label": "G synthetic",
        "gaia_color": gaia_color,
        "gaia_mag": gaia_mag,
        "apex_color": apex_color,
        "apex_mag": apex_mag,
        "matched": matched,
        "snr_bands": snr_bands,
    }


def _native_gaia_transformed_cmd_arrays(df: pd.DataFrame) -> dict | None:
    axes = select_cmd_qc_axes(df, system="std")
    if axes is None:
        return None
    a = str(axes["color_a"])
    b = str(axes["color_b"])
    y_band = str(axes["y_band"])
    required = [_cmd_mag_col("std", a), _cmd_mag_col("std", b), _cmd_mag_col("std", y_band)]
    if not all(col in df.columns for col in required):
        return None

    apex_a = pd.to_numeric(df[_cmd_mag_col("std", a)], errors="coerce").to_numpy(dtype=float)
    apex_b = pd.to_numeric(df[_cmd_mag_col("std", b)], errors="coerce").to_numpy(dtype=float)
    apex_y = pd.to_numeric(df[_cmd_mag_col("std", y_band)], errors="coerce").to_numpy(dtype=float)

    gaia_a = _gaia_reference_mag_for_band(df, a)
    gaia_b = _gaia_reference_mag_for_band(df, b)
    gaia_y = _gaia_reference_mag_for_band(df, y_band)

    apex_color = apex_a - apex_b
    gaia_color = gaia_a - gaia_b
    matched = np.isfinite(gaia_color) & np.isfinite(gaia_y) & np.isfinite(apex_color) & np.isfinite(apex_y)
    if not np.any(matched):
        return None
    return {
        "comparison": "apex_standard_minus_gaia_transformed_standard",
        "mode": "native_standard",
        "color_label": f"{a} - {b}",
        "apex_color_label": f"{a} - {b}",
        "mag_label": y_band,
        "apex_mag_label": f"{y_band} std",
        "gaia_color": gaia_color,
        "gaia_mag": gaia_y,
        "apex_color": apex_color,
        "apex_mag": apex_y,
        "matched": matched,
        "snr_bands": list(dict.fromkeys([a, b, y_band])),
    }


def build_gaia_cmd_comparison(df: pd.DataFrame) -> pd.DataFrame:
    """Compare final APEX CMD with Gaia on the same matched-star CMD axes."""
    if not isinstance(df, pd.DataFrame) or df.empty:
        return pd.DataFrame()
    arrs = _synthetic_gaia_cmd_arrays(df)
    if arrs is None:
        arrs = _native_gaia_transformed_cmd_arrays(df)
    if arrs is None:
        return pd.DataFrame()

    gaia_mag = arrs["gaia_mag"]
    gaia_color = arrs["gaia_color"]
    apex_mag = arrs["apex_mag"]
    apex_color = arrs["apex_color"]
    matched = arrs["matched"]
    gaia_mask = np.isfinite(gaia_mag) & np.isfinite(gaia_color)
    apex_mask = np.isfinite(apex_mag) & np.isfinite(apex_color)

    d_mag = apex_mag[matched] - gaia_mag[matched]
    d_color = apex_color[matched] - gaia_color[matched]
    rows = [{
        "comparison": arrs["comparison"],
        "mode": arrs["mode"],
        "cmd_color": arrs["color_label"],
        "cmd_y": arrs["mag_label"],
        "n_total_sources": int(len(df)),
        "n_gaia_cmd": int(np.sum(gaia_mask)),
        "n_apex_cmd": int(np.sum(apex_mask)),
        "n_apex_synthetic_gaia_cmd": int(np.sum(apex_mask)) if arrs["mode"] == "synthetic_gaia" else 0,
        "n_matched_cmd": int(np.sum(matched)),
        "matched_fraction_of_gaia": float(np.sum(matched) / max(np.sum(gaia_mask), 1)),
        "median_delta_mag": float(np.nanmedian(d_mag)),
        "sigma_mad_delta_mag": _zp_mad_sigma(d_mag),
        "p90_abs_delta_mag": float(np.nanpercentile(np.abs(d_mag), 90)),
        "median_delta_color": float(np.nanmedian(d_color)),
        "sigma_mad_delta_color": _zp_mad_sigma(d_color),
        "p90_abs_delta_color": float(np.nanpercentile(np.abs(d_color), 90)),
        "median_delta_G": float(np.nanmedian(d_mag)) if arrs["mode"] == "synthetic_gaia" else np.nan,
        "sigma_mad_delta_G": _zp_mad_sigma(d_mag) if arrs["mode"] == "synthetic_gaia" else np.nan,
        "p90_abs_delta_G": float(np.nanpercentile(np.abs(d_mag), 90)) if arrs["mode"] == "synthetic_gaia" else np.nan,
        "median_delta_BP_RP": float(np.nanmedian(d_color)) if arrs["mode"] == "synthetic_gaia" else np.nan,
        "sigma_mad_delta_BP_RP": _zp_mad_sigma(d_color) if arrs["mode"] == "synthetic_gaia" else np.nan,
        "p90_abs_delta_BP_RP": float(np.nanpercentile(np.abs(d_color), 90)) if arrs["mode"] == "synthetic_gaia" else np.nan,
        "median_gaia_mag": float(np.nanmedian(gaia_mag[matched])),
        "median_apex_mag": float(np.nanmedian(apex_mag[matched])),
        "median_gaia_G": float(np.nanmedian(gaia_mag[matched])) if arrs["mode"] == "synthetic_gaia" else np.nan,
        "median_apex_G_syn": float(np.nanmedian(apex_mag[matched])) if arrs["mode"] == "synthetic_gaia" else np.nan,
        "basis": "same matched IDs with finite Gaia-derived and APEX calibrated CMD magnitudes",
    }]
    return pd.DataFrame(rows)


def build_gaia_cmd_drift_table(df: pd.DataFrame, bin_width: float = 0.5) -> pd.DataFrame:
    """Binned-median Delta(mag)/Delta(color) vs Gaia reference magnitude.

    This automates the bright->faint drift diagnostic that exposed the
    NGC 6811 B-band faint bias: a magnitude-dependent median offset between
    the APEX calibrated CMD and the Gaia-transformed reference is invisible
    in the global medians but jumps out of this table. The final row is a
    ``DRIFT`` summary: median over the faint quintile minus median over the
    bright quintile (by matched-star magnitude), per axis.
    """
    if not isinstance(df, pd.DataFrame) or df.empty:
        return pd.DataFrame()
    arrs = _synthetic_gaia_cmd_arrays(df)
    if arrs is None:
        arrs = _native_gaia_transformed_cmd_arrays(df)
    if arrs is None:
        return pd.DataFrame()

    matched = arrs["matched"]
    gaia_mag = np.asarray(arrs["gaia_mag"], dtype=float)[matched]
    d_mag = np.asarray(arrs["apex_mag"], dtype=float)[matched] - gaia_mag
    d_color = (
        np.asarray(arrs["apex_color"], dtype=float)[matched]
        - np.asarray(arrs["gaia_color"], dtype=float)[matched]
    )
    ok = np.isfinite(gaia_mag) & np.isfinite(d_mag) & np.isfinite(d_color)
    if int(ok.sum()) < 30:
        return pd.DataFrame()
    gaia_mag, d_mag, d_color = gaia_mag[ok], d_mag[ok], d_color[ok]

    rows: list[dict] = []
    lo = np.floor(np.nanmin(gaia_mag) / bin_width) * bin_width
    hi = np.ceil(np.nanmax(gaia_mag) / bin_width) * bin_width
    edges = np.arange(lo, hi + 0.5 * bin_width, bin_width)
    for b_lo, b_hi in zip(edges[:-1], edges[1:]):
        k = (gaia_mag >= b_lo) & (gaia_mag < b_hi)
        if int(k.sum()) < 15:
            continue
        rows.append({
            "kind": "bin",
            "mag_lo": float(b_lo),
            "mag_hi": float(b_hi),
            "n": int(k.sum()),
            "median_delta_mag": float(np.nanmedian(d_mag[k])),
            "sigma_mad_delta_mag": _zp_mad_sigma(d_mag[k]),
            "median_delta_color": float(np.nanmedian(d_color[k])),
            "sigma_mad_delta_color": _zp_mad_sigma(d_color[k]),
        })

    # DRIFT summary: faint quintile minus bright quintile.
    q20, q80 = np.nanpercentile(gaia_mag, [20.0, 80.0])
    bright, faint = gaia_mag <= q20, gaia_mag >= q80
    if int(bright.sum()) >= 15 and int(faint.sum()) >= 15:
        rows.append({
            "kind": "DRIFT",
            "mag_lo": float(q20),   # bright quintile upper edge
            "mag_hi": float(q80),   # faint quintile lower edge
            "n": int(bright.sum() + faint.sum()),
            "median_delta_mag": float(np.nanmedian(d_mag[faint]) - np.nanmedian(d_mag[bright])),
            "sigma_mad_delta_mag": np.nan,
            "median_delta_color": float(np.nanmedian(d_color[faint]) - np.nanmedian(d_color[bright])),
            "sigma_mad_delta_color": np.nan,
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out.insert(0, "cmd_y", arrs["mag_label"])
        out.insert(0, "cmd_color", arrs["color_label"])
    return out


def _robust_line_fit(x, y) -> tuple[float, float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    if len(x) < 5:
        return np.nan, np.nan, np.nan, int(len(x))
    for _ in range(5):
        slope, intercept = np.polyfit(x, y, 1)
        resid = y - (slope * x + intercept)
        med = float(np.nanmedian(resid))
        sig = _zp_mad_sigma(resid)
        if not np.isfinite(sig) or sig <= 0:
            break
        keep = np.abs(resid - med) <= 3.0 * sig
        if int(np.sum(keep)) == len(x):
            break
        x = x[keep]
        y = y[keep]
        if len(x) < 5:
            break
    if len(x) < 5:
        return np.nan, np.nan, np.nan, int(len(x))
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    return float(slope), float(intercept), _zp_mad_sigma(resid), int(len(x))


def _snr_cut_mask(df: pd.DataFrame, bands: list[str], threshold: float) -> tuple[np.ndarray, list[str]]:
    mask = np.ones(len(df), dtype=bool)
    used: list[str] = []
    if threshold <= 0:
        return mask, used
    for band in bands:
        col = f"snr_{band}"
        if col not in df.columns:
            continue
        vals = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
        mask &= np.isfinite(vals) & (vals >= float(threshold))
        used.append(col)
    return mask, used


def build_gaia_cmd_snr_sweep(
    df: pd.DataFrame,
    snr_cuts: tuple[float, ...] = (5, 10, 20, 50, 100),
) -> pd.DataFrame:
    """Measure Gaia/APEX CMD residual sensitivity to the adopted SNR cut."""
    if not isinstance(df, pd.DataFrame) or df.empty:
        return pd.DataFrame()
    arrs = _synthetic_gaia_cmd_arrays(df)
    if arrs is None:
        arrs = _native_gaia_transformed_cmd_arrays(df)
    if arrs is None:
        return pd.DataFrame()

    gaia_mag = arrs["gaia_mag"]
    gaia_color = arrs["gaia_color"]
    apex_mag = arrs["apex_mag"]
    apex_color = arrs["apex_color"]
    base = np.asarray(arrs["matched"], dtype=bool)
    d_mag = apex_mag - gaia_mag
    d_color = apex_color - gaia_color
    snr_bands = list(arrs.get("snr_bands", []))

    rows = []
    for cut in snr_cuts:
        snr_mask, used_cols = _snr_cut_mask(df, snr_bands, float(cut))
        use = base & snr_mask
        n = int(np.sum(use))
        if n > 0:
            slope_mag, _, scatter_mag_fit, nfit_mag = _robust_line_fit(gaia_mag[use], d_mag[use])
            slope_color, _, scatter_color_fit, nfit_color = _robust_line_fit(gaia_color[use], d_color[use])
            dm = d_mag[use]
            dc = d_color[use]
            row = {
                "snr_cut": float(cut),
                "snr_columns": ",".join(used_cols),
                "comparison": arrs["comparison"],
                "mode": arrs["mode"],
                "cmd_color": arrs["color_label"],
                "cmd_y": arrs["mag_label"],
                "n_matched": int(np.sum(base)),
                "n_used": n,
                "used_fraction": float(n / max(int(np.sum(base)), 1)),
                "median_delta_mag": float(np.nanmedian(dm)),
                "sigma_mad_delta_mag": _zp_mad_sigma(dm),
                "p90_abs_delta_mag": float(np.nanpercentile(np.abs(dm), 90)),
                "median_delta_color": float(np.nanmedian(dc)),
                "sigma_mad_delta_color": _zp_mad_sigma(dc),
                "p90_abs_delta_color": float(np.nanpercentile(np.abs(dc), 90)),
                "slope_delta_mag_vs_mag": slope_mag,
                "fit_scatter_delta_mag": scatter_mag_fit,
                "nfit_delta_mag": nfit_mag,
                "slope_delta_color_vs_color": slope_color,
                "fit_scatter_delta_color": scatter_color_fit,
                "nfit_delta_color": nfit_color,
            }
        else:
            row = {
                "snr_cut": float(cut),
                "snr_columns": ",".join(used_cols),
                "comparison": arrs["comparison"],
                "mode": arrs["mode"],
                "cmd_color": arrs["color_label"],
                "cmd_y": arrs["mag_label"],
                "n_matched": int(np.sum(base)),
                "n_used": 0,
                "used_fraction": 0.0,
                "median_delta_mag": np.nan,
                "sigma_mad_delta_mag": np.nan,
                "p90_abs_delta_mag": np.nan,
                "median_delta_color": np.nan,
                "sigma_mad_delta_color": np.nan,
                "p90_abs_delta_color": np.nan,
                "slope_delta_mag_vs_mag": np.nan,
                "fit_scatter_delta_mag": np.nan,
                "nfit_delta_mag": 0,
                "slope_delta_color_vs_color": np.nan,
                "fit_scatter_delta_color": np.nan,
                "nfit_delta_color": 0,
            }
        rows.append(row)
    return pd.DataFrame(rows)


def draw_gaia_cmd_snr_sweep(fig: Figure, sweep_df: pd.DataFrame) -> bool:
    if not isinstance(sweep_df, pd.DataFrame) or sweep_df.empty:
        return False
    fig.clear()
    axes = fig.subplots(2, 2, squeeze=False)
    ax_n, ax_med = axes[0]
    ax_scatter, ax_slope = axes[1]

    x = pd.to_numeric(sweep_df["snr_cut"], errors="coerce").to_numpy(dtype=float)
    order = np.argsort(x)
    x = x[order]
    work = sweep_df.iloc[order].reset_index(drop=True)

    n_used = pd.to_numeric(work["n_used"], errors="coerce").to_numpy(dtype=float)
    used_frac = pd.to_numeric(work["used_fraction"], errors="coerce").to_numpy(dtype=float)
    med_mag = pd.to_numeric(work["median_delta_mag"], errors="coerce").to_numpy(dtype=float)
    med_color = pd.to_numeric(work["median_delta_color"], errors="coerce").to_numpy(dtype=float)
    sig_mag = pd.to_numeric(work["sigma_mad_delta_mag"], errors="coerce").to_numpy(dtype=float)
    sig_color = pd.to_numeric(work["sigma_mad_delta_color"], errors="coerce").to_numpy(dtype=float)
    slope_mag = pd.to_numeric(work["slope_delta_mag_vs_mag"], errors="coerce").to_numpy(dtype=float)
    slope_color = pd.to_numeric(work["slope_delta_color_vs_color"], errors="coerce").to_numpy(dtype=float)

    ax_n.plot(x, n_used, "o-", color="#4C78A8", label="N used")
    ax_n.set_xlabel("SNR cut")
    ax_n.set_ylabel("matched stars")
    ax_n.grid(True, alpha=0.25)
    ax_n2 = ax_n.twinx()
    ax_n2.plot(x, used_frac * 100.0, "s--", color="#E45756", label="used %")
    ax_n2.set_ylabel("used (%)")
    ax_n2.set_ylim(0, 105)
    h1, l1 = ax_n.get_legend_handles_labels()
    h2, l2 = ax_n2.get_legend_handles_labels()
    ax_n.legend(h1 + h2, l1 + l2, fontsize=8, frameon=False)
    ax_n.set_title("Sample Retention")

    ax_med.axhline(0.0, color="#222222", lw=0.8, ls="--", alpha=0.6)
    ax_med.plot(x, med_mag, "o-", color="#54A24B", label="median dMag")
    ax_med.plot(x, med_color, "o-", color="#B279A2", label="median dColor")
    ax_med.set_xlabel("SNR cut")
    ax_med.set_ylabel("median APEX - Gaia")
    ax_med.set_title("Median Residual")
    ax_med.grid(True, alpha=0.25)
    ax_med.legend(fontsize=8, frameon=False)

    ax_scatter.plot(x, sig_mag, "o-", color="#54A24B", label="dMag MAD sigma")
    ax_scatter.plot(x, sig_color, "o-", color="#B279A2", label="dColor MAD sigma")
    ax_scatter.set_xlabel("SNR cut")
    ax_scatter.set_ylabel("robust scatter")
    ax_scatter.set_title("Residual Scatter")
    ax_scatter.grid(True, alpha=0.25)
    ax_scatter.legend(fontsize=8, frameon=False)

    ax_slope.axhline(0.0, color="#222222", lw=0.8, ls="--", alpha=0.6)
    ax_slope.plot(x, slope_mag, "o-", color="#54A24B", label="dMag vs mag slope")
    ax_slope.plot(x, slope_color, "o-", color="#B279A2", label="dColor vs color slope")
    ax_slope.set_xlabel("SNR cut")
    ax_slope.set_ylabel("slope")
    ax_slope.set_title("Residual Trend")
    ax_slope.grid(True, alpha=0.25)
    ax_slope.legend(fontsize=8, frameon=False)

    first = work.iloc[0]
    fig.suptitle(
        f"Gaia/APEX CMD SNR-Cut Sensitivity | {first.get('mode', '')} | "
        f"{first.get('cmd_color', '')} vs {first.get('cmd_y', '')}",
        fontsize=12,
    )
    fig.subplots_adjust(left=0.08, right=0.95, bottom=0.08, top=0.90, wspace=0.46, hspace=0.34)
    return True


def draw_gaia_cmd_comparison(fig: Figure, df: pd.DataFrame, summary_df: pd.DataFrame | None = None) -> bool:
    """Draw Gaia reference CMD beside APEX calibrated CMD on matched-star axes."""
    if not isinstance(df, pd.DataFrame) or df.empty:
        return False
    arrs = _synthetic_gaia_cmd_arrays(df)
    if arrs is None:
        arrs = _native_gaia_transformed_cmd_arrays(df)
    if arrs is None:
        return False

    gaia_mag = arrs["gaia_mag"]
    gaia_color = arrs["gaia_color"]
    apex_mag = arrs["apex_mag"]
    apex_color = arrs["apex_color"]
    matched = arrs["matched"]
    d_mag = apex_mag[matched] - gaia_mag[matched]
    d_color = apex_color[matched] - gaia_color[matched]

    fig.clear()
    axes = fig.subplots(2, 2, squeeze=False)
    ax_gaia, ax_apex = axes[0]
    ax_dmag, ax_dc = axes[1]

    all_color = np.concatenate([gaia_color[matched], apex_color[matched]])
    all_mag = np.concatenate([gaia_mag[matched], apex_mag[matched]])
    c_lo, c_hi = np.nanpercentile(all_color, [1, 99])
    m_lo, m_hi = np.nanpercentile(all_mag, [1, 99])
    c_pad = max(0.05, 0.08 * (c_hi - c_lo))
    m_pad = max(0.05, 0.08 * (m_hi - m_lo))

    common_scatter = {
        "s": 5,
        "alpha": 0.55,
        "linewidths": 0,
        "rasterized": True,
    }
    ax_gaia.scatter(gaia_color[matched], gaia_mag[matched], color="#4C78A8", **common_scatter)
    gaia_title = "Gaia Observed CMD" if arrs["mode"] == "synthetic_gaia" else "Gaia-Transformed CMD"
    ax_gaia.set_title(gaia_title)
    ax_gaia.set_xlabel(arrs["color_label"])
    ax_gaia.set_ylabel(arrs["mag_label"])
    ax_gaia.set_xlim(c_lo - c_pad, c_hi + c_pad)
    ax_gaia.set_ylim(m_hi + m_pad, m_lo - m_pad)
    ax_gaia.grid(True, alpha=0.2)

    ax_apex.scatter(apex_color[matched], apex_mag[matched], color="#F58518", **common_scatter)
    apex_title = "APEX Calibrated CMD in Gaia Space" if arrs["mode"] == "synthetic_gaia" else "APEX Calibrated CMD"
    ax_apex.set_title(apex_title)
    ax_apex.set_xlabel(arrs["apex_color_label"])
    ax_apex.set_ylabel(arrs["apex_mag_label"])
    ax_apex.set_xlim(c_lo - c_pad, c_hi + c_pad)
    ax_apex.set_ylim(m_hi + m_pad, m_lo - m_pad)
    ax_apex.grid(True, alpha=0.2)

    ax_dmag.scatter(gaia_mag[matched], d_mag, color="#54A24B", **common_scatter)
    dmag_med = float(np.nanmedian(d_mag))
    dmag_sig = _zp_mad_sigma(d_mag)
    ax_dmag.axhline(dmag_med, color="#222222", lw=1.0, ls="--", label=f"median {dmag_med:+.3f}")
    if np.isfinite(dmag_sig) and dmag_sig > 0:
        ax_dmag.axhspan(dmag_med - dmag_sig, dmag_med + dmag_sig, color="#54A24B", alpha=0.12)
    ax_dmag.set_xlabel(f"Gaia reference {arrs['mag_label']}")
    ax_dmag.set_ylabel(f"APEX - Gaia ({arrs['mag_label']})")
    ax_dmag.set_title("Magnitude Residual")
    ax_dmag.grid(True, alpha=0.2)
    ax_dmag.legend(fontsize=8, frameon=False)

    ax_dc.scatter(gaia_color[matched], d_color, color="#B279A2", **common_scatter)
    dc_med = float(np.nanmedian(d_color))
    dc_sig = _zp_mad_sigma(d_color)
    ax_dc.axhline(dc_med, color="#222222", lw=1.0, ls="--", label=f"median {dc_med:+.3f}")
    if np.isfinite(dc_sig) and dc_sig > 0:
        ax_dc.axhspan(dc_med - dc_sig, dc_med + dc_sig, color="#B279A2", alpha=0.12)
    ax_dc.set_xlabel(f"Gaia reference {arrs['color_label']}")
    ax_dc.set_ylabel(f"APEX - Gaia ({arrs['color_label']})")
    ax_dc.set_title("Color Residual")
    ax_dc.grid(True, alpha=0.2)
    ax_dc.legend(fontsize=8, frameon=False)

    n = int(np.sum(matched))
    title_extra = ""
    if isinstance(summary_df, pd.DataFrame) and not summary_df.empty:
        row = summary_df.iloc[0]
        title_extra = (
            f" | median dMag={float(row.get('median_delta_mag', np.nan)):+.3f}, "
            f"median dColor={float(row.get('median_delta_color', np.nan)):+.3f}"
        )
    fig.suptitle(f"Gaia vs APEX CMD Comparison | {arrs['mode']} | matched N={n}{title_extra}", fontsize=12)
    fig.subplots_adjust(left=0.08, right=0.97, bottom=0.08, top=0.91, wspace=0.28, hspace=0.34)
    return True


def export_gaia_cmd_comparison_products(output_dir: Path, log_func=None) -> list[Path]:
    """Export Gaia-observed vs APEX-synthetic CMD comparison products."""
    output_dir = Path(output_dir)
    cmd_path = output_dir / "median_by_ID_filter_wide_cmd.csv"
    if not cmd_path.exists() or cmd_path.stat().st_size == 0:
        return []
    try:
        df = pd.read_csv(cmd_path)
    except Exception:
        return []

    summary = build_gaia_cmd_comparison(df)
    if summary.empty:
        return []

    saved: list[Path] = []
    summary_path = output_dir / "gaia_cmd_comparison_summary.csv"
    summary.to_csv(summary_path, index=False)
    saved.append(summary_path)

    fig = Figure(figsize=(11.2, 8.0), dpi=120)
    if draw_gaia_cmd_comparison(fig, df, summary):
        fig_path = output_dir / "step10_gaia_cmd_comparison.png"
        fig.savefig(fig_path, dpi=160, bbox_inches="tight")
        saved.append(fig_path)

    sweep = build_gaia_cmd_snr_sweep(df)
    if not sweep.empty:
        sweep_path = output_dir / "gaia_cmd_snr_sweep.csv"
        sweep.to_csv(sweep_path, index=False)
        saved.append(sweep_path)
        sweep_fig = Figure(figsize=(10.8, 7.6), dpi=120)
        if draw_gaia_cmd_snr_sweep(sweep_fig, sweep):
            sweep_fig_path = output_dir / "step10_gaia_cmd_snr_sweep.png"
            sweep_fig.savefig(sweep_fig_path, dpi=160, bbox_inches="tight")
            saved.append(sweep_fig_path)

    drift = build_gaia_cmd_drift_table(df)
    if not drift.empty:
        drift_path = output_dir / "gaia_cmd_drift_by_mag.csv"
        drift.to_csv(drift_path, index=False)
        saved.append(drift_path)
        if log_func is not None:
            summary_row = drift[drift["kind"] == "DRIFT"]
            if len(summary_row):
                dm = float(summary_row["median_delta_mag"].iloc[0])
                dc = float(summary_row["median_delta_color"].iloc[0])
                try:
                    log_func(
                        f"[Gaia QC] bright->faint drift: dMag={dm:+.4f}, dColor={dc:+.4f} mag "
                        "(|drift|>~0.02 = magnitude-dependent calibration systematic)"
                    )
                except Exception:
                    pass

    if saved and log_func is not None:
        try:
            log_func("Step10 Gaia CMD comparison exported: " + ", ".join(path.name for path in saved))
        except Exception:
            pass
    return saved


class ZeropointCalibrationWorker(QThread):
    progress = pyqtSignal(int, int, str)
    log = pyqtSignal(str)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)

    def __init__(self, params, data_dir: Path, result_dir: Path, cache_dir: Path):
        super().__init__()
        self.params = params
        self.data_dir = Path(data_dir)
        self.result_dir = Path(result_dir)
        self.cache_dir = Path(cache_dir)
        self._stop_requested = False
        self.last_summary: dict = {}
        self.last_error = ""

    def stop(self):
        self._stop_requested = True

    def _log(self, msg: str):
        self.log.emit(msg)

    def _pick_col(self, cols, cands):
        for c in cands:
            if c in cols:
                return c
        return None

    def _resolve_path(self, p):
        p = str(p) if p is not None else ""
        if p.strip() == "":
            return None
        p0 = Path(p)
        if p0.is_absolute() and p0.exists():
            return p0
        # Support Windows absolute paths when running under POSIX-like envs.
        if len(p) >= 3 and p[1] == ":" and (p[2] == "\\" or p[2] == "/"):
            drv = p[0].lower()
            rest = p[2:].replace("\\", "/").lstrip("/")
            p_win = Path(f"/mnt/{drv}/{rest}")
            if p_win.exists():
                return p_win
        for base in (
            step7_forced_phot_dir(self.result_dir),
            step8_psf_dir(self.result_dir),
            step9_selection_dir(self.result_dir),
            self.result_dir,
            self.result_dir / "phot",
            self.result_dir / "photometry",
            self.result_dir / "result",
        ):
            p1 = base / p
            if p1.exists():
                return p1
        return None

    def _inject_psf_master_identity(
        self,
        table: pd.DataFrame,
        filename: str,
    ) -> pd.DataFrame:
        """Attach Step 7 master IDs to a Step 8 PSF table."""
        step7_path = step7_forced_phot_dir(self.result_dir) / f"photometry_{filename}.tsv"
        if not step7_path.exists():
            return table
        try:
            master = pd.read_csv(step7_path, sep="\t")
        except Exception:
            return table
        if "ID" not in master.columns:
            return table

        output = table.copy()
        matched_id = pd.Series(np.nan, index=output.index, dtype=float)
        match_distance = pd.Series(np.inf, index=output.index, dtype=float)
        position_matches = 0
        uid_matches = 0

        master_x_col = "x_fit" if "x_fit" in master.columns else "x"
        master_y_col = "y_fit" if "y_fit" in master.columns else "y"
        psf_x_col = "x_fit" if "x_fit" in output.columns else "x"
        psf_y_col = "y_fit" if "y_fit" in output.columns else "y"
        if (
            {master_x_col, master_y_col} <= set(master.columns)
            and {psf_x_col, psf_y_col} <= set(output.columns)
        ):
            master_xy = master[[master_x_col, master_y_col]].apply(
                pd.to_numeric, errors="coerce"
            ).to_numpy(float)
            master_ids = pd.to_numeric(master["ID"], errors="coerce").to_numpy(float)
            psf_xy = output[[psf_x_col, psf_y_col]].apply(
                pd.to_numeric, errors="coerce"
            ).to_numpy(float)
            valid_master = np.isfinite(master_xy).all(axis=1) & np.isfinite(master_ids)
            valid_psf = np.isfinite(psf_xy).all(axis=1)
            if valid_master.any() and valid_psf.any():
                radius = max(
                    0.1,
                    float(getattr(self.params.P, "psf_cmd_match_radius_px", 1.0)),
                )
                tree = cKDTree(master_xy[valid_master])
                distances, indices = tree.query(
                    psf_xy[valid_psf],
                    distance_upper_bound=radius,
                )
                hit = np.isfinite(distances)
                psf_rows = output.index.to_numpy()[valid_psf][hit]
                selected_ids = master_ids[valid_master][indices[hit]]
                matched_id.loc[psf_rows] = selected_ids
                match_distance.loc[psf_rows] = distances[hit]
                position_matches = int(hit.sum())

        seed_column = "seed_uid" if "seed_uid" in output.columns else "det_uid"
        if seed_column in output.columns and "det_uid" in master.columns:
            master_uid = pd.to_numeric(master["det_uid"], errors="coerce").astype("Int64")
            master_id = pd.to_numeric(master["ID"], errors="coerce")
            uid_map = (
                pd.DataFrame({"uid": master_uid, "ID": master_id})
                .dropna(subset=["uid", "ID"])
                .drop_duplicates("uid", keep="first")
                .set_index("uid")["ID"]
            )
            seeds = pd.to_numeric(output[seed_column], errors="coerce").astype("Int64")
            uid_ids = seeds.map(uid_map)
            # ``seed_uid == -1`` is Step 8's unmatched-residual sentinel, not
            # a stable Step 7 identity. Other negative values are valid forced
            # catalog UIDs and remain eligible for the fallback.
            use_uid = matched_id.isna() & uid_ids.notna() & seeds.ne(-1)
            matched_id.loc[use_uid] = uid_ids.loc[use_uid].astype(float)
            uid_matches = int(use_uid.sum())

        output["ID"] = matched_id
        output = output[output["ID"].notna()].copy()
        output["_psf_id_match_distance"] = match_distance.loc[output.index]
        if "mag_psf_err" in output.columns:
            output["_psf_sort_error"] = pd.to_numeric(
                output["mag_psf_err"], errors="coerce"
            ).fillna(np.inf)
        else:
            output["_psf_sort_error"] = np.inf
        before_dedup = len(output)
        output = output.sort_values(
            ["_psf_id_match_distance", "_psf_sort_error"],
            kind="stable",
        ).drop_duplicates("ID", keep="first")
        duplicates = before_dedup - len(output)
        output = output.drop(
            columns=["_psf_id_match_distance", "_psf_sort_error"],
            errors="ignore",
        )
        self._log(
            f"[ZP][PSF ID] {filename}: matched={len(output)}/{len(table)} "
            f"(position={position_matches}, uid_fallback={uid_matches}, "
            f"duplicate={duplicates})"
        )
        return output

    @staticmethod
    def _robust_median_and_err(arr, per_measurement_err=None):
        """Median plus its uncertainty, from the scatter between measurements.

        ``per_measurement_err`` is the fallback for when that scatter cannot be
        measured. One frame has no scatter, so ``MAD_TO_SIGMA * 0 / 1`` is
        exactly zero — and a reported error of zero claims the magnitude is
        known perfectly. Across the five Phase-3 clusters that was 0.4-1.0 % of
        entries (12-42 stars each), every one of them reporting 0.000. A star
        that dropped from five frames to one after the Step 8 seed fix then
        carried a 2.7 mag outlier at zero uncertainty (M5 ID 986, 2026-08-15).
        Identical values across frames give MAD = 0 the same way.

        Falling back to the star's own photometric error says what is actually
        known. With nothing to fall back on, NaN — unknown, not perfect.
        """
        x = np.asarray(arr, float)
        finite = np.isfinite(x)
        x = x[finite]
        if len(x) == 0:
            return (np.nan, np.nan, 0)
        med = float(np.median(x))
        mad = float(np.median(np.abs(x - med)))
        err = float(MAD_TO_SIGMA * mad / np.sqrt(max(len(x), 1)))
        if not (err > 0):
            err = np.nan
            if per_measurement_err is not None:
                own = np.asarray(per_measurement_err, float)
                if own.shape == finite.shape:
                    own = own[finite]
                own = own[np.isfinite(own) & (own > 0)]
                if own.size:
                    # Median of the frames that went in, divided the same way,
                    # so the one-frame case reads on the same scale as the rest.
                    err = float(np.median(own) / np.sqrt(len(x)))
        return (med, err, int(len(x)))

    @staticmethod
    def _weighted_mean_mag(mag, mag_err, clip_sigma=3.0, iters=4):
        mags = np.asarray(mag, float)
        errs = np.asarray(mag_err, float)
        mask = np.isfinite(mags) & np.isfinite(errs) & (errs > 0)
        if mask.sum() == 0:
            return (np.nan, np.nan, 0)
        x = mags[mask]
        e = errs[mask]
        for _ in range(int(iters)):
            med = np.nanmedian(x)
            mad = np.nanmedian(np.abs(x - med))
            sig = MAD_TO_SIGMA * mad if mad > 0 else np.nanstd(x)
            if not np.isfinite(sig) or sig <= 0:
                break
            keep = np.abs(x - med) <= float(clip_sigma) * sig
            if keep.sum() == len(x):
                break
            x = x[keep]
            e = e[keep]
            if len(x) == 0:
                break
        n = int(len(x))
        if n == 0:
            return (np.nan, np.nan, 0)

        flux = 10.0 ** (-0.4 * x)
        sigma_flux = flux * (np.log(10.0) / 2.5) * e
        w = np.where(np.isfinite(sigma_flux) & (sigma_flux > 0), 1.0 / (sigma_flux ** 2), 0.0)
        wsum = float(np.nansum(w))
        if not np.isfinite(wsum) or wsum <= 0:
            return (np.nan, np.nan, n)
        flux_w = float(np.nansum(w * flux) / wsum)
        if not np.isfinite(flux_w) or flux_w <= 0:
            return (np.nan, np.nan, n)
        sigma_flux_w = float(np.sqrt(1.0 / wsum))
        mag_w = float(-2.5 * np.log10(flux_w))
        mag_w_err = float(MAG_ERR_COEFF * (sigma_flux_w / flux_w))
        return (mag_w, mag_w_err, n)

    @staticmethod
    def _robust_clip(x, clip_sigma=3.0, iters=5):
        x = np.asarray(x, float)
        x = x[np.isfinite(x)]
        if len(x) == 0:
            return x
        for _ in range(int(iters)):
            med = np.nanmedian(x)
            mad = np.nanmedian(np.abs(x - med))
            sig = MAD_TO_SIGMA * mad if mad > 0 else np.nanstd(x)
            if not np.isfinite(sig) or sig <= 0:
                break
            keep = np.abs(x - med) <= float(clip_sigma) * sig
            if keep.sum() == len(x):
                break
            x = x[keep]
            if len(x) == 0:
                break
        return x

    def _robust_location(self, arr, clip_sigma=3.0, iters=5):
        x = self._robust_clip(arr, clip_sigma=clip_sigma, iters=iters)
        if len(x) == 0:
            return (np.nan, np.nan, 0, np.nan)
        med = float(np.nanmedian(x))
        std = float(np.nanstd(x))
        n = int(len(x))
        outlier_frac = np.nan
        try:
            outlier_frac = float(1.0 - (len(x) / max(len(arr), 1)))
        except Exception:
            outlier_frac = np.nan
        return (med, std, n, outlier_frac)

    def _has_wcs(self, header):
        try:
            w0 = self._wcs_from_header(header)
            return bool(w0.has_celestial)
        except Exception:
            return False

    @staticmethod
    def _wcs_from_header(header):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            return WCS(header)

    def _list_frames(self):
        crop_active = crop_is_active(self.result_dir)
        cropped_dir = step2_cropped_dir(self.result_dir)
        if crop_active and cropped_dir.exists() and list(cropped_dir.glob("*.fit*")):
            return sorted(cropped_dir.glob("*.fit*"))
        return sorted(self.data_dir.glob("*.fit*"))

    def _find_idmatch_csv(self, fname: str, result_dir: Path):
        """Find per-frame photometry TSV from forced phot (replaces old idmatch CSV)."""
        forced_dir = step7_forced_phot_dir(result_dir)
        direct = forced_dir / f"photometry_{fname}.tsv"
        if direct.exists():
            return direct
        return None

    def _load_ref_wcs(self):
        P = self.params.P
        files = self._list_frames()

        ref_val = getattr(P, "ref_frame", None)
        if ref_val is not None and str(ref_val).strip() != "":
            ref_txt = str(ref_val).strip()
            if ref_txt.isdigit():
                idx = int(ref_txt)
                if files and 0 <= idx < len(files):
                    fp = files[idx]
                    hdr = fits.getheader(fp)
                    if self._has_wcs(hdr):
                        self._log(f"WCS from ref_frame index: {fp.name}")
                        return self._wcs_from_header(hdr), fp
            else:
                fp = Path(ref_txt)
                if not fp.is_absolute():
                    cropped_dir = step2_cropped_dir(self.result_dir)
                    search_dirs = [self.result_dir, self.data_dir]
                    if crop_is_active(self.result_dir):
                        if cropped_dir.exists():
                            search_dirs.append(cropped_dir)
                    for base in search_dirs:
                        cand = base / fp
                        if cand.exists():
                            fp = cand
                            break
                if fp.exists():
                    hdr = fits.getheader(fp)
                    if self._has_wcs(hdr):
                        self._log(f"WCS from ref_frame: {fp.name}")
                        return self._wcs_from_header(hdr), fp

        # Check cropped directory FIRST (coordinates are from cropped images)
        patterns = ["ref*.fit*", "rc_*.fit*", "Crop_*.fit*", "crop_*.fit*", "*.fit*"]
        cropped_dir = step2_cropped_dir(self.result_dir)
        search_dirs = [self.result_dir, self.data_dir]
        if crop_is_active(self.result_dir):
            if cropped_dir.exists():
                search_dirs.insert(0, cropped_dir)
        for base in search_dirs:
            if not base.exists():
                continue
            for pat in patterns:
                for fp in sorted(base.glob(pat)):
                    try:
                        hdr = fits.getheader(fp)
                        if self._has_wcs(hdr):
                            self._log(f"WCS auto-detected: {fp.name} (from {base.name})")
                            return self._wcs_from_header(hdr), fp
                    except Exception:
                        continue
        return None, None

    def _build_frame_airmass(self, idx: pd.DataFrame) -> pd.DataFrame:
        from apex.utils.astro_utils import compute_airmass_from_header
        P = self.params.P
        lat = float(getattr(P, "site_lat_deg", 0.0))
        lon = float(getattr(P, "site_lon_deg", 0.0))
        alt = float(getattr(P, "site_alt_m", 0.0))
        tz = float(getattr(P, "site_tz_offset_hours", 0.0))

        rows = []
        for _, r in idx.iterrows():
            fname = str(r.get("file", "")).strip()
            if fname == "":
                continue
            fpath = self.data_dir / fname
            if not fpath.exists():
                if crop_is_active(self.result_dir):
                    cand = step2_cropped_dir(self.result_dir) / fname
                    if cand.exists():
                        fpath = cand
            if not fpath.exists():
                continue
            try:
                hdr = fits.getheader(fpath)
                info = compute_airmass_from_header(hdr, lat, lon, alt, tz)
                filt = normalize_filter_name(r.get("filter", hdr.get("FILTER", "")))
                rows.append({
                    "file": fname,
                    "filter": filt,
                    **info,
                })
            except Exception:
                continue
        df = pd.DataFrame(rows)
        if df.empty:
            df = pd.DataFrame(columns=["file", "filter", "airmass", "airmass_source", "alt_deg", "zenith_deg", "datetime_utc", "datetime_local", "ra_deg", "dec_deg"])
        if len(df):
            output_dir = step10_zp_dir(self.result_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            out_path = output_dir / "frame_airmass.csv"
            df.to_csv(out_path, index=False)
            self._log(f"Saved {out_path.name} | rows={len(df)}")
        return df

    def _write_nonlinearity_diag(self, obs: pd.DataFrame, output_dir) -> None:
        """Bridge-star detector-linearity diagnostic for multi-exposure sets.

        For stars measured at >= 2 distinct exposure times, compare the
        calibrated magnitude between exposure levels. A linear detector gives a
        flat delta-vs-magnitude relation; a brightness-dependent slope flags
        non-linearity that the (per-frame-constant) zeropoint cannot absorb.

        Writes ``nonlinearity_by_exposure.csv`` (per star/filter/level median),
        ``nonlinearity_summary.csv`` (per filter/level-pair slope+offset+scatter)
        and ``nonlinearity_check.png``. Quietly skips single-exposure data.
        """
        try:
            need = {"ID", "FILTER", "mag_cal", "exptime"}
            if not need <= set(obs.columns):
                return
            d = obs[["ID", "FILTER", "mag_cal", "exptime"]].copy()
            d["mag_cal"] = pd.to_numeric(d["mag_cal"], errors="coerce")
            d["exptime"] = pd.to_numeric(d["exptime"], errors="coerce")
            d["snr"] = pd.to_numeric(obs["snr"], errors="coerce") if "snr" in obs.columns else np.nan
            d = d[np.isfinite(d["mag_cal"]) & np.isfinite(d["exptime"]) & (d["exptime"] > 0)]
            if d.empty:
                return
            # Exposure level = rounded exptime (absorbs sub-second jitter).
            d["exp_level"] = d["exptime"].round(1)
            if d["exp_level"].nunique() < 2:
                self._log("[ZP][NLIN] single exposure level; linearity check skipped.")
                return

            grp = (
                d.groupby(["ID", "FILTER", "exp_level"], as_index=False)
                 .agg(mag=("mag_cal", "median"), n=("mag_cal", "size"), snr=("snr", "median"))
            )
            grp.to_csv(output_dir / "nonlinearity_by_exposure.csv", index=False, na_rep="NaN")

            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            filters = sorted(grp["FILTER"].dropna().astype(str).unique().tolist())
            params_holder = self.__dict__.get("params")
            params_obj = getattr(params_holder, "P", None)
            try:
                cmd_snr_min = getattr(params_obj, "cmd_snr_calib_min", 20.0)
                snr_min = float(
                    getattr(params_obj, "gaia_snr_calib_min", cmd_snr_min)
                )
            except (TypeError, ValueError):
                snr_min = 20.0
            snr_min = max(0.0, snr_min)
            summary = []
            panels = []  # (filter, short, long, x, delta, slope, intercept)
            for filt in filters:
                sub = grp[grp["FILTER"].astype(str) == filt]
                flevels = sorted(sub["exp_level"].unique().tolist())
                for short_lv, long_lv in zip(flevels[:-1], flevels[1:]):
                    a = (
                        sub[sub["exp_level"] == short_lv][["ID", "mag", "snr"]]
                        .rename(columns={"mag": "mag_s", "snr": "snr_s"})
                    )
                    b = (
                        sub[sub["exp_level"] == long_lv][["ID", "mag", "snr"]]
                        .rename(columns={"mag": "mag_l", "snr": "snr_l"})
                    )
                    m = a.merge(b, on="ID", how="inner")
                    x = m["mag_l"].to_numpy(float)
                    delta = (m["mag_s"].to_numpy(float) - x)
                    ok = np.isfinite(x) & np.isfinite(delta)
                    snr_s = m["snr_s"].to_numpy(float)
                    snr_l = m["snr_l"].to_numpy(float)
                    if np.isfinite(snr_s).any() or np.isfinite(snr_l).any():
                        ok &= (
                            np.isfinite(snr_s)
                            & np.isfinite(snr_l)
                            & (snr_s >= snr_min)
                            & (snr_l >= snr_min)
                        )
                    x, delta = x[ok], delta[ok]
                    if len(delta) < 5:
                        continue

                    # Reject large fit residuals without flattening a real
                    # brightness-dependent trend.
                    for _ in range(3):
                        slope, intercept = (float(v) for v in np.polyfit(x, delta, 1))
                        residual = delta - (slope * x + intercept)
                        residual_med = float(np.median(residual))
                        residual_sigma = float(
                            MAD_TO_SIGMA * np.median(np.abs(residual - residual_med))
                        )
                        if not np.isfinite(residual_sigma) or residual_sigma <= 1e-12:
                            break
                        keep = np.abs(residual - residual_med) <= 3.0 * residual_sigma
                        if keep.all() or int(keep.sum()) < 5:
                            break
                        x, delta = x[keep], delta[keep]
                    slope, intercept = (float(v) for v in np.polyfit(x, delta, 1))
                    med = float(np.median(delta))
                    scatter = float(MAD_TO_SIGMA * np.median(np.abs(delta - med)))
                    flagged = abs(slope) > 0.02
                    summary.append({
                        "FILTER": filt, "exp_short_s": short_lv, "exp_long_s": long_lv,
                        "n_bridge": int(len(delta)), "delta_median": med,
                        "slope_mag_per_mag": slope, "scatter_mag": scatter,
                        "snr_min": snr_min,
                        "nonlinearity_flag": bool(flagged),
                    })
                    panels.append((filt, short_lv, long_lv, x, delta, slope, intercept))
                    self._log(
                        f"[ZP][NLIN] {filt} {short_lv:g}s-{long_lv:g}s: n={len(delta)} "
                        f"SNR>={snr_min:g} slope={slope:+.4f} mag/mag "
                        f"offset={med:+.4f} scatter={scatter:.4f}"
                        + ("  *** slope flags non-linearity ***" if flagged else "")
                    )

            if summary:
                pd.DataFrame(summary).to_csv(output_dir / "nonlinearity_summary.csv", index=False, na_rep="NaN")

            if panels:
                ncol = min(3, len(panels))
                nrow = int(np.ceil(len(panels) / ncol))
                fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.4 * nrow), squeeze=False)
                for ax, (filt, s_lv, l_lv, x, delta, slope, intercept) in zip(axes.ravel(), panels):
                    ax.axhline(0.0, color="0.6", lw=0.8, ls="--")
                    ax.scatter(x, delta, s=10, alpha=0.5, color="#1f77b4")
                    xs = np.linspace(np.min(x), np.max(x), 50)
                    ax.plot(xs, slope * xs + intercept, color="#d62728", lw=1.5,
                            label=f"slope={slope:+.3f}")
                    ax.set_title(f"{filt}: {s_lv:g}s - {l_lv:g}s", fontsize=9)
                    ax.set_xlabel("mag_cal (long)", fontsize=8)
                    ax.set_ylabel("Δmag (short - long)", fontsize=8)
                    ax.legend(fontsize=7, loc="best")
                for ax in axes.ravel()[len(panels):]:
                    ax.set_visible(False)
                fig.suptitle("Bridge-star linearity check (flat = linear)", fontsize=10)
                fig.tight_layout()
                fig.savefig(output_dir / "nonlinearity_check.png", dpi=110)
                plt.close(fig)
                self._log(f"Saved nonlinearity_check.png | {len(panels)} exposure-pair panel(s)")
        except Exception as e:
            self._log(f"[ZP][NLIN] linearity diagnostic failed: {e}")

    def _poly_eval(self, x, coeffs):
        x = np.asarray(x, float)
        y = np.zeros_like(x, dtype=float)
        p = np.ones_like(x, dtype=float)
        for a in coeffs:
            y += a * p
            p *= x
        return y

    def _robust_linfit(self, x, y, w=None, clip_sigma=3.0, iters=5, slope_absmax=1.0, min_n=10):
        """Robust sigma-clipping linear fit.  w = weights (e.g. 1/err²); if None, OLS.

        Returns (zp, ct, n_inlier, scatter_rms) where scatter_rms is the MAD-based
        scatter of residuals from the final fit.
        """
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        m0 = np.isfinite(x) & np.isfinite(y)
        if w is not None:
            w = np.asarray(w, float)
            m0 &= np.isfinite(w) & (w > 0)
        x, y = x[m0], y[m0]
        if w is not None:
            w = w[m0]

        def _scatter(yf, zp_, ct_, xf):
            r = yf - (zp_ + ct_ * xf)
            return float(MAD_TO_SIGMA * np.nanmedian(np.abs(r - np.nanmedian(r)))) if len(r) else np.nan

        if len(x) < min_n:
            med = float(np.nanmedian(y)) if len(y) else np.nan
            return med, 0.0, int(len(x)), _scatter(y, med, 0.0, x)

        def _fit(xf, yf, wf):
            if wf is not None:
                sw = np.sqrt(wf)
                A = np.column_stack([xf * sw, sw])
                result, _, _, _ = np.linalg.lstsq(A, yf * sw, rcond=None)
                return result[1], result[0]   # zp, ct
            else:
                ct_, zp_ = np.polyfit(xf, yf, 1)
                return zp_, ct_

        zp, ct = _fit(x, y, w)

        for _ in range(int(iters)):
            yhat = zp + ct * x
            r = y - yhat
            med = np.nanmedian(r)
            mad = np.nanmedian(np.abs(r - med)) + 1e-12
            sig = MAD_TO_SIGMA * mad
            keep = np.abs(r - med) <= float(clip_sigma) * sig
            if keep.sum() < min_n:
                break
            zp, ct = _fit(x[keep], y[keep], w[keep] if w is not None else None)
            x, y = x[keep], y[keep]
            if w is not None:
                w = w[keep]

        scatter = _scatter(y, zp, ct, x)

        if abs(ct) > float(slope_absmax):
            self._log(f"[ZP] color term {ct:+.4f} exceeds slope_absmax={slope_absmax:.2f}; using median ZP, ct=0")
            return float(np.nanmedian(y)), 0.0, int(len(x)), scatter

        return float(zp), float(ct), int(len(x)), scatter

    def _apply_standard_anchor(self, df_out: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
        """Re-anchor mag_std_* on an external (Gaia-independent) standard
        catalog when [cmd.standard_anchor] is enabled.  Validated on M67:
        the Gaia 'approx' U reference was off by −0.13 mag and drove the
        isochrone [M/H] to a confident wrong value; anchoring U/B/V on
        MMJ93 (VizieR J/AJ/106/181, Landolt-tied) restored the literature
        solution.  Any failure here must not sink Step 10 — it degrades to
        the un-anchored table with a warning."""
        P = self.params.P
        if not getattr(P, "std_anchor_enable", False):
            return df_out
        catalog = str(getattr(P, "std_anchor_catalog", "") or "").strip()
        if not catalog:
            self._log("[ANCHOR] enabled but no catalog id set — skipped "
                      "(cmd.standard_anchor.catalog, e.g. \"J/AJ/106/181\")")
            # Help the user: search VizieR for candidates and log them.
            try:
                from apex.analysis.cmd.standard_anchor import discover_standard_catalogs

                name = str(getattr(P, "target_name", "") or "").strip()
                ra = getattr(P, "target_ra_deg", None) or getattr(P, "ra_deg", None)
                dec = getattr(P, "target_dec_deg", None) or getattr(P, "dec_deg", None)
                if name and ra is not None and dec is not None:
                    cands = discover_standard_catalogs(
                        name, float(ra), float(dec),
                        ["U", "B", "V", "R", "I", "g", "r", "i"], log=self._log)
                    for c in cands[:3]:
                        self._log(f"[ANCHOR] 후보: {c.catalog_id} "
                                  f"({'/'.join(c.bands)}, "
                                  f"{'시야 내' if c.in_field else '시야 밖?'}) — "
                                  f"{c.description[:50]}")
                    if not cands:
                        self._log("[ANCHOR] VizieR에서 이 시야의 표준 측광 카탈로그를 "
                                  "찾지 못함 — Gaia 참조 유지")
            except Exception as exc:
                self._log(f"[ANCHOR] 후보 탐색 실패(무시): {type(exc).__name__}: {exc}")
            return df_out
        try:
            from apex.analysis.cmd.standard_anchor import (
                anchor_qc_frame,
                apply_anchor,
                compute_anchor,
                fetch_standard_catalog,
            )

            cache_dir = Path(self.result_dir) / "cache" / "standard_anchor"
            std = fetch_standard_catalog(catalog, cache_dir, log=self._log)
            result = compute_anchor(
                df_out, std, catalog,
                match_radius_arcsec=float(getattr(P, "std_anchor_match_radius", 1.5)),
                min_stars=int(getattr(P, "std_anchor_min_stars", 20)),
                log=self._log,
            )
            for msg in result.warnings:
                self._log(f"[ANCHOR][warn] {msg}")
            qc = anchor_qc_frame(result)
            if not qc.empty:
                qc.to_csv(output_dir / "standard_anchor_offsets.csv", index=False)
            if result.residuals is not None and not result.residuals.empty:
                result.residuals.to_csv(
                    output_dir / "standard_anchor_residuals.csv", index=False)
            if not result.applied_bands:
                self._log("[ANCHOR] no band anchored — table left on the Gaia reference")
                return df_out
            self._log(f"[ANCHOR] applied to mag_std_{{{','.join(result.applied_bands)}}} "
                      f"(catalog {catalog}, {result.n_matched} matches, "
                      f"sep median {result.sep_median_arcsec:.2f}\")")
            return apply_anchor(df_out, result)
        except Exception as exc:
            self._log(f"[ANCHOR][warn] standard anchor failed — continuing "
                      f"un-anchored: {type(exc).__name__}: {exc}")
            return df_out

    def run(self):
        try:
            P = self.params.P
            result_dir = self.result_dir
            output_dir = step10_zp_dir(result_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            source_info = resolve_cmd_photometry_input(result_dir)
            idx_path = Path(source_info["index_path"])
            if not idx_path.exists() or idx_path.stat().st_size <= 0:
                idx_path = next(
                    (
                        path
                        for path in _cmd_photometry_index_candidates(result_dir)
                        if path.exists() and path.stat().st_size > 0
                    ),
                    None,
                )
            if idx_path is None:
                raise FileNotFoundError("photometry index csv not found (or all candidates are empty)")

            try:
                idx = pd.read_csv(idx_path)
            except pd.errors.EmptyDataError:
                raise FileNotFoundError(f"photometry index csv is empty: {idx_path.name}")
            self._log(
                f"Photometry = {str(source_info['source']).upper()} | "
                f"MAG = {source_info['mag_column']} | {source_info['reason']}"
            )
            self._log(f"Index = {idx_path.parent.name}/{idx_path.name} | rows={len(idx)}")

            if "path" not in idx.columns:
                for cand in ("phot_tsv", "tsv", "out", "output"):
                    if cand in idx.columns:
                        idx = idx.rename(columns={cand: "path"})
                        break

            # Fallback: construct path from file column (photometry_{fname}.tsv in same dir)
            if "path" not in idx.columns and "file" in idx.columns:
                phot_dir = idx_path.parent
                idx["path"] = idx["file"].apply(
                    lambda f: str(phot_dir / f"photometry_{f}.tsv")
                )
                self._log(f"[ZP] path column missing — constructed from file column ({phot_dir.name}/)")

            if "file" not in idx.columns:
                c_file = self._pick_col(idx.columns, ["fname", "frame", "image", "fits", "name"])
                if c_file:
                    idx = idx.rename(columns={c_file: "file"})

            if "filter" in idx.columns:
                idx["filter"] = idx["filter"].map(normalize_filter_name)
            elif "FILTER" in idx.columns:
                idx["filter"] = idx["FILTER"].map(normalize_filter_name)
            else:
                idx["filter"] = "unknown"

            use_qc = should_use_frame_quality_qc(
                result_dir,
                P,
                "phot_use_qc_pass_only",
                default=False,
            )
            idx, qc_info = filter_frame_df_by_qc(result_dir, idx, file_col="file", require_qc=use_qc)
            if use_qc:
                if qc_info.get("applied"):
                    self._log(f"Step4 QC passed only: {qc_info['total']} -> {qc_info['kept']}")
                elif qc_info.get("path") is None:
                    self._log("Step4 QC: frame_quality.csv not found; using all frames.")
                else:
                    self._log(f"Step4 QC: frame_quality.csv ignored ({qc_info['reason']}); using all frames.")

            # Canonical default is 3.0 everywhere (schema ge=1.0, step8, TOML loaders);
            # a 0.0 fallback here would silently disable the per-measurement SNR floor
            # that guards the CMD faint end against positive-flux noise fluctuations
            # (worse with union/forced photometry of undetected stars in short frames).
            min_snr_for_mag = float(getattr(P, "min_snr_for_mag", 3.0))
            # Restrict calibration measurements to Step 4 apcorr-quality refs
            # (isolated, unsaturated, high-flux). Per-measurement: a star is
            # kept only on frames where it was an apcorr_candidate. Falls back
            # to all measurements if too few survive (see below).
            require_apcorr = bool(getattr(P, "phot_ref_require_apcorr_candidate", True))
            apcorr_min_keep = int(getattr(P, "phot_ref_apcorr_min_keep", 8))

            # Pre-load sourceid_to_ID for fallback ID injection (det_uid → source_id → ID)
            _sid_map = None
            for _cand_dir in (step9_selection_dir(result_dir),):
                _sid_csv = _cand_dir / "sourceid_to_ID.csv"
                if _sid_csv.exists():
                    try:
                        _sid_df = pd.read_csv(_sid_csv)
                        if "source_id" in _sid_df.columns and "ID" in _sid_df.columns:
                            _sid_map = _sid_df.set_index("source_id")["ID"].to_dict()
                            self._log(f"[ZP] sourceid_to_ID: {len(_sid_map)} entries ({_sid_csv.parent.name}/)")
                    except Exception:
                        pass
                    break

            rows = []
            n_missing = 0
            missing_examples = []
            total = len(idx)
            for i, (_, r) in enumerate(idx.iterrows(), start=1):
                if self._stop_requested:
                    self.finished.emit({"stopped": True})
                    return
                p = self._resolve_path(r.get("path", ""))
                if p is None or (not p.exists()):
                    n_missing += 1
                    if len(missing_examples) < 5:
                        missing_examples.append(str(r.get("path", "")))
                    continue

                try:
                    dfp = pd.read_csv(p, sep="\t")
                except Exception:
                    dfp = pd.read_csv(p)

                if source_info.get("source") == "psf" and "flags_psf" in dfp.columns:
                    psf_flags = pd.to_numeric(dfp["flags_psf"], errors="coerce")
                    dfp = dfp[psf_flags.eq(0)].copy()
                if source_info.get("source") == "psf" and "mag_psf" in dfp.columns:
                    psf_magnitude = pd.to_numeric(dfp["mag_psf"], errors="coerce")
                    dfp = dfp[np.isfinite(psf_magnitude)].copy()

                if "is_saturated" in dfp.columns:
                    dfp = dfp[~dfp["is_saturated"].fillna(False).astype(bool)]
                if "is_nonlinear" in dfp.columns:
                    dfp = dfp[~dfp["is_nonlinear"].fillna(False).astype(bool)]
                if "centroid_outlier" in dfp.columns:
                    dfp = dfp[~dfp["centroid_outlier"].fillna(False).astype(bool)]
                if "recenter_capped" in dfp.columns:
                    dfp = dfp[~dfp["recenter_capped"].fillna(False).astype(bool)]

                if "ID" not in dfp.columns and source_info.get("source") == "psf":
                    dfp = self._inject_psf_master_identity(
                        dfp,
                        str(r.get("file", p.name)),
                    )

                if "ID" not in dfp.columns:
                    # Fallback: inject ID via idmatch join (det_uid → source_id → ID)
                    if "det_uid" in dfp.columns and _sid_map is not None:
                        _fname = str(r.get("file", "")) if "file" in idx.columns else ""
                        if not _fname:
                            _stem = p.stem
                            _fname = _stem[len("photometry_"):] if _stem.startswith("photometry_") else _stem
                        _idmatch_csv = self._find_idmatch_csv(_fname, result_dir)
                        if _idmatch_csv is not None:
                            try:
                                sep = "\t" if _idmatch_csv.suffix.lower() == ".tsv" else ","
                                _im = pd.read_csv(_idmatch_csv, sep=sep)
                                if "det_idx" in _im.columns:
                                    _im = _im.rename(columns={"det_idx": "det_uid"})
                                if "det_uid" in _im.columns and "source_id" in _im.columns:
                                    _had_sid = "source_id" in dfp.columns
                                    _im_sub = _im[["det_uid", "source_id"]].drop_duplicates("det_uid")
                                    dfp = dfp.merge(_im_sub, on="det_uid", how="left")
                                    dfp["ID"] = dfp["source_id"].map(_sid_map)
                                    if not _had_sid:
                                        dfp = dfp.drop(columns=["source_id"], errors="ignore")
                            except Exception as _inj_e:
                                self._log(f"[ZP] ID injection failed for {p.name}: {_inj_e}")
                    if "ID" not in dfp.columns:
                        raise RuntimeError(f"{p.name}: ID column missing")

                if "FILTER" in dfp.columns:
                    dfp["FILTER"] = dfp["FILTER"].map(normalize_filter_name)
                else:
                    dfp["FILTER"] = normalize_filter_name(r.get("filter", "unknown"))

                mag_col = None
                for cand in ("mag_psf", "mag_inst", "mag", "mag_ap", "mag_apcorr"):
                    if cand in dfp.columns:
                        mag_col = cand
                        break
                if mag_col is None:
                    raise RuntimeError(f"{p.name}: mag column missing")

                err_col = None
                for cand in ("mag_psf_err", "mag_err", "emag", "emag_inst", "magerr"):
                    if cand in dfp.columns:
                        err_col = cand
                        break
                if err_col is None:
                    dfp["mag_err"] = np.nan
                    err_col = "mag_err"

                provenance = build_photometry_provenance(
                    mag_column=mag_col,
                    mag_error_column=err_col,
                )

                snr_col = "snr" if "snr" in dfp.columns else ("snr_psf" if "snr_psf" in dfp.columns else None)

                _keep_cols = ["ID", "FILTER", mag_col, err_col] + ([snr_col] if snr_col else [])
                _has_apcorr_col = "step4_apcorr_candidate" in dfp.columns
                if _has_apcorr_col:
                    _keep_cols.append("step4_apcorr_candidate")
                _has_exptime_col = "exptime" in dfp.columns
                if _has_exptime_col:
                    _keep_cols.append("exptime")
                tmp = dfp[_keep_cols].copy()
                if not _has_exptime_col:
                    tmp["exptime"] = np.nan
                tmp = tmp.rename(columns={mag_col: "mag_inst", err_col: "mag_err"})
                tmp["photometry_source"] = provenance["source"]
                tmp["mag_input_column"] = provenance["mag_column"]
                tmp["mag_error_input_column"] = provenance["mag_error_column"]
                if _has_apcorr_col:
                    tmp["step4_apcorr_candidate"] = (
                        tmp["step4_apcorr_candidate"].astype(str).str.strip()
                        .str.lower().isin({"1", "true", "t", "yes", "y"})
                    )
                if snr_col is None:
                    tmp["snr"] = np.nan
                else:
                    tmp = tmp.rename(columns={snr_col: "snr"})
                if np.isnan(tmp["mag_err"].to_numpy(float)).all():
                    snr_vals = tmp["snr"].to_numpy(float)
                    if np.isfinite(snr_vals).any():
                        _snr_mask = np.isfinite(snr_vals) & (snr_vals > 0)
                        tmp.loc[_snr_mask, "mag_err"] = MAG_ERR_COEFF / snr_vals[_snr_mask]

                tmp["file"] = str(r.get("file", "")) if ("file" in idx.columns) else ""

                if min_snr_for_mag > 0:
                    m = np.isfinite(tmp["snr"].to_numpy(float))
                    tmp = tmp[(~m) | (tmp["snr"].to_numpy(float) >= min_snr_for_mag)].copy()

                rows.append(tmp)
                self.progress.emit(i, total, str(r.get("file", "")))

            self._log(f"Read frames: {len(rows)} | missing paths: {n_missing}")
            if missing_examples:
                self._log(f"Missing path samples: {missing_examples}")

            all_df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["ID", "FILTER", "mag_inst", "mag_err", "snr", "file"])
            all_df["FILTER"] = all_df["FILTER"].map(normalize_filter_name)
            if all_df.empty:
                raise RuntimeError(
                    "No readable photometry rows. "
                    f"missing_paths={n_missing}/{total}. "
                    "Check photometry_index.csv path column and Step5/Step6 photometry outputs."
                )

            # NOTE: the apcorr-quality reference filter is applied ONLY to the
            # frame-zeropoint fit reference set (see below, where `obs` is
            # restricted to calibrators). It must NOT filter `all_df` here,
            # because `all_df` also builds the per-star instrumental table that
            # the CMD is plotted from — filtering it would drop faint stars from
            # the CMD entirely (apcorr_candidate requires flux_pct >= 60%).

            def _combine_group_raw(g):
                med, med_err, n_med = self._robust_median_and_err(
                    g["mag_inst"], g["mag_err"])
                wmean, werr, _ = self._weighted_mean_mag(g["mag_inst"], g["mag_err"])
                # snr_psf (PSF photometry) or snr (aperture photometry)
                _snr_col = "snr" if "snr" in g.columns else ("snr_psf" if "snr_psf" in g.columns else None)
                snr_vals = np.asarray(g[_snr_col], float) if _snr_col else np.array([np.nan])
                snr_med = float(np.nanmedian(snr_vals)) if np.isfinite(snr_vals).any() else np.nan
                return pd.Series({
                    "mag_inst_med": med,
                    "mag_inst_med_err": med_err,
                    "mag_inst_wmean": wmean,
                    "mag_inst_werr": werr,
                    "n_frames": n_med,
                    "snr_med": snr_med,
                    "photometry_source": collapse_provenance_values(g["photometry_source"]),
                    "mag_input_column": collapse_provenance_values(g["mag_input_column"]),
                    "mag_error_input_column": collapse_provenance_values(g["mag_error_input_column"]),
                })

            grp_raw = all_df.groupby(["ID", "FILTER"], as_index=False).apply(_combine_group_raw)
            provenance_by_filter = {
                filt: build_photometry_provenance(
                    collapse_provenance_values(group["photometry_source"]),
                    collapse_provenance_values(group["mag_input_column"]),
                    collapse_provenance_values(group["mag_error_input_column"]),
                )
                for filt, group in all_df.groupby("FILTER")
            }

            grp_raw_path = output_dir / "median_by_ID_filter_raw.csv"
            grp_raw.to_csv(grp_raw_path, index=False, na_rep="NaN")
            self._log(f"Saved {grp_raw_path.name} | rows={len(grp_raw)}")

            wide_raw_mag_w = grp_raw.pivot_table(index="ID", columns="FILTER", values="mag_inst_wmean", aggfunc="median")
            wide_raw_err_w = grp_raw.pivot_table(index="ID", columns="FILTER", values="mag_inst_werr", aggfunc="median")
            wide_raw_mag_med = grp_raw.pivot_table(index="ID", columns="FILTER", values="mag_inst_med", aggfunc="median")
            wide_raw_err_med = grp_raw.pivot_table(index="ID", columns="FILTER", values="mag_inst_med_err", aggfunc="median")
            wide_raw_mag = wide_raw_mag_w.combine_first(wide_raw_mag_med)
            wide_raw_err = wide_raw_err_w.combine_first(wide_raw_err_med)
            wide_raw_snr = grp_raw.pivot_table(index="ID", columns="FILTER", values="snr_med", aggfunc="median")
            wide_raw_source = grp_raw.pivot_table(index="ID", columns="FILTER", values="photometry_source", aggfunc="first")
            wide_raw_mag_input = grp_raw.pivot_table(index="ID", columns="FILTER", values="mag_input_column", aggfunc="first")
            wide_raw_err_input = grp_raw.pivot_table(index="ID", columns="FILTER", values="mag_error_input_column", aggfunc="first")

            wide_raw_mag.columns = [f"mag_inst_{c}" for c in wide_raw_mag.columns]
            wide_raw_err.columns = [f"mag_inst_err_{c}" for c in wide_raw_err.columns]
            wide_raw_snr.columns = [f"snr_{c}" for c in wide_raw_snr.columns]
            wide_raw_mag_w.columns = [f"mag_inst_wmean_{c}" for c in wide_raw_mag_w.columns]
            wide_raw_err_w.columns = [f"mag_inst_werr_{c}" for c in wide_raw_err_w.columns]
            wide_raw_mag_med.columns = [f"mag_inst_med_{c}" for c in wide_raw_mag_med.columns]
            wide_raw_err_med.columns = [f"mag_inst_med_err_{c}" for c in wide_raw_err_med.columns]
            wide_raw_source.columns = [f"photometry_source_{c}" for c in wide_raw_source.columns]
            wide_raw_mag_input.columns = [f"mag_input_column_{c}" for c in wide_raw_mag_input.columns]
            wide_raw_err_input.columns = [f"mag_error_input_column_{c}" for c in wide_raw_err_input.columns]

            wide_raw = pd.concat(
                [
                    wide_raw_mag,
                    wide_raw_err,
                    wide_raw_mag_w,
                    wide_raw_err_w,
                    wide_raw_mag_med,
                    wide_raw_err_med,
                    wide_raw_snr,
                    wide_raw_source,
                    wide_raw_mag_input,
                    wide_raw_err_input,
                ],
                axis=1,
            ).reset_index()

            wide_raw_path = output_dir / "median_by_ID_filter_wide_raw.csv"
            wide_raw.to_csv(wide_raw_path, index=False, na_rep="NaN")
            self._log(f"Saved {wide_raw_path.name} | rows={len(wide_raw)}")

            master, master_source, master_path = _load_master_table(result_dir)
            self._log(f"Master table: {master_path.name} ({master_source})")
            if "ID" not in master.columns and "source_id" in master.columns and _sid_map:
                sid_s = parse_int64_series(master["source_id"]).astype("Int64")
                id_vals = sid_s.map(_sid_map)
                master["ID"] = pd.to_numeric(id_vals, errors="coerce").astype("Int64")
                master = master[master["ID"].notna()].copy()
            if "ID" not in master.columns:
                raise RuntimeError(f"{master_path.name} missing ID column")

            master = _merge_gaia_columns_from_catalog(master, result_dir)
            if "parallax" in master.columns:
                n_plx = int(pd.to_numeric(master["parallax"], errors="coerce").notna().sum())
                if n_plx > 0:
                    self._log(f"Gaia astrometry attached: {n_plx} / {len(master)}")

            # Merge wide with master to get Gaia mags from master_catalog
            merge_cols = ["ID"]
            if "source_id" in master.columns:
                merge_cols.append("source_id")
            col_xm = self._pick_col(master.columns, ["x_ref", "x", "x_pix", "x_center", "x_med"])
            col_ym = self._pick_col(master.columns, ["y_ref", "y", "y_pix", "y_center", "y_med"])
            if col_xm:
                merge_cols.append(col_xm)
            if col_ym:
                merge_cols.append(col_ym)
            for col in ("ra_deg", "dec_deg",
                        "gaia_G", "gaia_BP", "gaia_RP", "gmag", "bpmag", "rpmag", "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag",
                        "parallax", "parallax_error", "pmra", "pmdec", "pmra_error", "pmdec_error",
                        "ruwe", "phot_bp_rp_excess_factor"):
                if col in master.columns and col not in merge_cols:
                    merge_cols.append(col)

            df = wide_raw.merge(master[merge_cols], on="ID", how="left")
            if col_xm:
                df = df.rename(columns={col_xm: "x_pix"})
            if col_ym:
                df = df.rename(columns={col_ym: "y_pix"})
            g_col = self._pick_col(df.columns, ["gaia_G", "gmag", "phot_g_mean_mag"])
            bp_col = self._pick_col(df.columns, ["gaia_BP", "bpmag", "phot_bp_mean_mag"])
            rp_col = self._pick_col(df.columns, ["gaia_RP", "rpmag", "phot_rp_mean_mag"])

            gaia_from_master = True
            gaia_join_name = ""
            if g_col is None or bp_col is None or rp_col is None:
                gaia_from_master = False
                if "source_id" not in df.columns:
                    raise RuntimeError("master_catalog missing Gaia mags and source_id for Gaia join")
                gaia_candidates = [
                    step5_wcs_dir(result_dir) / "gaia_fov.ecsv",
                    step5_wcs_dir(result_dir) / "gaia_derived.csv",
                    result_dir / "gaia_derived.csv",
                    result_dir / "gaia_fov.ecsv",
                ]
                gaia_df = None
                gaia_path = None
                for cand in gaia_candidates:
                    if not cand.exists():
                        continue
                    try:
                        if cand.suffix.lower() == ".ecsv":
                            gaia_df = read_ecsv_int64_source_id(cand)
                        else:
                            gaia_df = pd.read_csv(cand)
                    except Exception:
                        gaia_df = None
                        continue
                    if gaia_df is None or gaia_df.empty:
                        gaia_df = None
                        continue
                    gaia_path = cand
                    break
                if gaia_df is None or gaia_path is None:
                    raise RuntimeError("master_catalog missing Gaia mags and step5_wcs/gaia_derived.csv not found")
                gaia_join_name = gaia_path.name
                if "source_id" in gaia_df.columns:
                    gaia_df["source_id"] = parse_int64_series(gaia_df["source_id"]).astype("Int64")
                    df["source_id"] = parse_int64_series(df["source_id"]).astype("Int64")
                gaia_cols = ["source_id", "phot_g_mean_mag"]
                if "phot_bp_mean_mag" in gaia_df.columns:
                    gaia_cols.append("phot_bp_mean_mag")
                if "phot_rp_mean_mag" in gaia_df.columns:
                    gaia_cols.append("phot_rp_mean_mag")
                df = df.merge(gaia_df[gaia_cols], on="source_id", how="left")
                g_col = self._pick_col(df.columns, ["phot_g_mean_mag"])
                bp_col = self._pick_col(df.columns, ["phot_bp_mean_mag"])
                rp_col = self._pick_col(df.columns, ["phot_rp_mean_mag"])

            if g_col is None or bp_col is None or rp_col is None:
                raise RuntimeError("Gaia magnitude columns not available (gaia_G/gaia_BP/gaia_RP or Gaia ECSV)")

            df["gaia_G"] = pd.to_numeric(df[g_col], errors="coerce")
            df["gaia_BP"] = pd.to_numeric(df[bp_col], errors="coerce")
            df["gaia_RP"] = pd.to_numeric(df[rp_col], errors="coerce")

            dfm = df[np.isfinite(df["gaia_G"]) & np.isfinite(df["gaia_BP"]) & np.isfinite(df["gaia_RP"])].copy()
            dfm["gaia_BP_RP"] = dfm["gaia_BP"] - dfm["gaia_RP"]
            src_note = str(master_source) if gaia_from_master else gaia_join_name
            self._log(f"Gaia mags from {src_note}: {len(dfm)} / {len(df)}")

            min_match = int(getattr(P, "min_master_gaia_matches", 10))
            if len(dfm) < min_match:
                dfm.to_csv(output_dir / "gaia_sdss_calibrator_by_ID.csv", index=False)
                raise RuntimeError("Not enough Gaia matches for calibration")

            out_cal = dfm.copy()

            xcol = out_cal["gaia_BP_RP"].to_numpy(float)
            G = out_cal["gaia_G"].to_numpy(float)

            # User-configurable global BP-RP pre-filter
            bpRP_lo = float(getattr(P, "gaia_gi_min", -0.5))
            bpRP_hi = float(getattr(P, "gaia_gi_max", 3.5))
            m_bpRP = np.isfinite(xcol) & (xcol >= bpRP_lo) & (xcol <= bpRP_hi)
            self._log(f"BP-RP pre-filter [{bpRP_lo:.2f}, {bpRP_hi:.2f}]: kept {m_bpRP.sum()}/{len(xcol)} stars")

            clip_sigma   = float(getattr(P, "zp_clip_sigma", 3.0))
            fit_iters    = int(getattr(P, "zp_fit_iters", 5))
            slope_absmax = float(getattr(P, "zp_slope_absmax", 1.0))
            snr_cut      = float(getattr(P, "gaia_snr_calib_min", getattr(P, "cmd_snr_calib_min", 20.0)))

            def _arr(col):
                return out_cal[col].to_numpy(float) if col in out_cal.columns else np.full(len(out_cal), np.nan)

            def _wls_weights(err_col):
                e = _arr(err_col)
                return np.where(np.isfinite(e) & (e > 0), 1.0 / e**2, np.nan)

            def _color_pair(fa, fb):
                ca, cb = f"mag_inst_{fa}", f"mag_inst_{fb}"
                if ca in out_cal.columns and cb in out_cal.columns:
                    return out_cal[ca].to_numpy(float) - out_cal[cb].to_numpy(float)
                return np.full(len(out_cal), np.nan)

            # Detect filters present in photometry data
            data_filters = sorted(all_df["FILTER"].dropna().unique())
            self._log(f"Photometric filters in data: {data_filters}")

            # ── Gaia → filter reference magnitude for each detected filter ──────
            self._log("=== Gaia → Filter Transformations ===")
            ref_col_map: dict[str, str] = {}  # filt -> column name in out_cal
            ref_source_map: dict[str, str] = {}
            for filt in data_filters:
                key = _BAND_ALIASES.get(filt, filt)
                if key not in _GAIA_TO_BAND:
                    self._log(f"[ZP][{filt}] No Gaia transformation available — skipping")
                    continue
                coeffs, lo, hi, source, sig_approx = _GAIA_TO_BAND[key]
                warn = " [WARNING: σ≈{:.2f} mag — use with caution]".format(sig_approx) if sig_approx >= 0.10 else ""
                self._log(f"[ZP][{filt}] {source}  G-{filt}=poly(BP-RP)  σ≈{sig_approx:.3f}{warn}")
                if sig_approx >= 0.10 and not getattr(self.params.P, "std_anchor_enable", False):
                    # Validated on M67: the approx U reference was off −0.13 mag
                    # and railed the isochrone [M/H]. The cure is an external,
                    # Gaia-independent standard catalog (Parameters > External
                    # Standard Anchor), not more colours.
                    self._log(f"[ZP][{filt}] 권고: 이 밴드의 Gaia 참조는 근사(σ≥0.1)입니다 — "
                              f"Parameters의 'External Standard Anchor'로 표준성 카탈로그에 "
                              f"재앵커하세요 (M67 실측: U 영점 -0.13 mag 편차가 이소크론 "
                              f"[M/H]를 rail시킴)")
                m_filt = m_bpRP & (xcol >= lo) & (xcol <= hi)
                G_minus_filt = np.full_like(G, np.nan)
                G_minus_filt[m_filt] = self._poly_eval(xcol[m_filt], coeffs)
                col_name = f"ref_{filt}"
                out_cal[col_name] = G - G_minus_filt
                ref_col_map[filt] = col_name
                ref_source_map[filt] = source
                self._log(f"[ZP][{filt}] ref_mag valid: {np.isfinite(out_cal[col_name]).sum()}/{len(out_cal)}")

            if not ref_col_map:
                raise RuntimeError("No supported filters found in photometry data for Gaia calibration")

            # ── ZP + color-term fit per filter ────────────────────────────────
            coeff_rows: list[dict] = []
            fit_params: dict[str, dict] = {}

            # Gaia calibrator quality: RUWE <= 1.4 + Riello+2021 |C*| <= 3sigma
            # (BP/RP contamination in crowded fields biases the transformed
            # reference mags of faint stars). Permissive when the columns are
            # absent (older master catalogs behave as before).
            # The C* cut (Riello+2021) could not run until 2026-08-11 because
            # no query fetched `phot_bp_rp_excess_factor`. Now that it can, it
            # is left OFF by default: measured on real catalogues it rejects
            # 3.6 % of Gaia references in M67 but 49.7 % in M13, because that
            # is exactly the crowding it detects. Turning it on is a science
            # decision that changes every globular-cluster zero point, so it is
            # the user's to make — set gaia.cstar_cut = true.
            _cstar_on = bool(getattr(P, "gaia_cstar_cut", False))
            m_gaia_qual, qual_report = gaia_quality_report(
                out_cal, cstar_nsigma=3.0 if _cstar_on else None)
            n_qual_cut = int(len(out_cal) - int(m_gaia_qual.sum()))
            if n_qual_cut:
                self._log(
                    f"[ZP] Gaia quality cut (RUWE<=1.4, |C*|<=3sig): removed "
                    f"{n_qual_cut}/{len(out_cal)} calibrator candidates"
                )
            # Say out loud which cuts ran. A skipped cut is not a warning-free
            # state: on M67 the RUWE cut ran in one run and not the next
            # (ESA TAP timed out; APEX's VizieR query SELECTs `ruwe` and its
            # ESA query does not), which moved the calibrator count by ~10 %
            # with no other visible symptom. Write it to disk so the number is
            # explainable from the outputs alone.
            for _cut, _info in qual_report["cuts"].items():
                if not _info.get("applied"):
                    self._log(
                        f"[ZP] Gaia quality cut '{_cut}' NOT APPLIED — "
                        f"{_info.get('reason')}. Calibrator counts are not "
                        f"comparable with runs where it did apply."
                    )
            try:
                # Which server answered decides which columns exist, so record
                # it next to the cuts rather than leaving them to be correlated
                # by hand with a file three directories away.
                _meta = output_dir.parent / "step5_wcs" / "gaia_fov_meta.json"
                qual_report["gaia_source"] = (
                    json.loads(_meta.read_text(encoding="utf-8")).get(
                        "gaia_source", "unknown")
                    if _meta.exists() else "unknown")
                (output_dir / "gaia_quality_report.json").write_text(
                    json.dumps(qual_report, indent=1), encoding="utf-8")
            except Exception as exc:      # diagnostics must never break the fit
                self._log(f"[ZP] could not write gaia_quality_report.json: {exc}")

            for filt in data_filters:
                if filt not in ref_col_map:
                    continue
                inst_col = f"mag_inst_{filt}"
                if inst_col not in out_cal.columns:
                    self._log(f"[ZP][{filt}] mag_inst_{filt} missing in calibrator table — skipping fit")
                    continue

                ref_mag_arr = _arr(ref_col_map[filt])
                inst_arr    = _arr(inst_col)
                delta       = ref_mag_arr - inst_arr

                # SNR mask for this filter
                snr_col_f = f"snr_{filt}"
                if snr_col_f in out_cal.columns:
                    sv = out_cal[snr_col_f].to_numpy(float)
                    m_snr_f = np.isfinite(sv) & (sv >= snr_cut)
                else:
                    m_snr_f = np.ones(len(out_cal), dtype=bool)

                # Find the best available color index. Prefer the REFERENCE
                # (standard) color when both bands have Gaia transformations:
                # the fit axis then matches the definition of the color term
                # and carries no instrumental noise correlated with delta
                # (the faint-end B-V coupling diagnosed on NGC 6811, 2026-07).
                key = _BAND_ALIASES.get(filt, filt)
                color_prefs = _FILTER_COLOR_PREF.get(key, _FILTER_COLOR_PREF.get(filt, []))
                color_x = np.full(len(out_cal), np.nan)
                color_col_name = "none"
                color_axis = "instrumental"
                for (ca, cb) in color_prefs:
                    ref_a, ref_b = ref_col_map.get(ca), ref_col_map.get(cb)
                    if ref_a and ref_b:
                        cidx = _arr(ref_a) - _arr(ref_b)
                        axis = "standard"
                    else:
                        cidx = _color_pair(ca, cb)
                        axis = "instrumental"
                    if np.isfinite(cidx).sum() >= min_match:
                        color_x = cidx
                        color_col_name = f"{ca}_{cb}"
                        color_axis = axis
                        break

                w_filt  = _wls_weights(f"mag_inst_err_{filt}")
                # |ct| must stay below the fixed-point contraction bound of
                # solve_standard_colors (error shrinks ~|ct| per pass, so
                # |ct| >= 1 diverges).  A former U-band exemption up to 3.0
                # let a wild fit (ct = -2.94, U vs approx Gaia reference)
                # pass and blow mag_std_U up to +-1000 mag.
                s_max   = min(slope_absmax, 0.8)

                # When no color index is available, fall back to ZP-only fit (CT forced to 0)
                if color_col_name == "none":
                    self._log(f"[ZP][{filt}] No instrumental color index available — fitting ZP only (CT=0)")
                    color_x = np.zeros(len(out_cal))

                m_fit   = np.isfinite(delta) & np.isfinite(color_x) & np.isfinite(inst_arr) & m_snr_f & m_gaia_qual

                if m_fit.sum() < min_match:
                    self._log(f"[ZP][{filt}] Only {m_fit.sum()} calibrators (need {min_match}) — skipping fit")
                    continue

                zp_f, ct_f, Nf, sc_f = self._robust_linfit(
                    color_x[m_fit], delta[m_fit], w=w_filt[m_fit],
                    clip_sigma=clip_sigma, iters=fit_iters, slope_absmax=s_max, min_n=min_match,
                )

                # Quadratic color-term refinement: the linear term leaves a
                # ±0.02-0.03 mag curvature vs color on rich fields (measured on
                # NGC 6811 B). Adopt the quadratic only when well constrained
                # (enough calibrators, sane coefficients, scatter not worse).
                ct2_f = 0.0
                if color_col_name != "none" and int(m_fit.sum()) >= _QUAD_MIN_CALIBRATORS:
                    q_coeffs, Nq, sq = robust_weighted_polyfit(
                        color_x[m_fit], delta[m_fit], w=w_filt[m_fit], degree=2,
                        clip_sigma=clip_sigma, iters=fit_iters, min_n=min_match,
                    )
                    # |ct2| <= 0.25: broadband transformation curvatures are
                    # ~0.05 (NGC 6811 B: -0.07); larger values are unphysical
                    # AND would threaten the fixed-point contraction in
                    # solve_standard_colors (derivative ~ ct + 2*ct2*C).
                    if (
                        q_coeffs is not None
                        and np.all(np.isfinite(q_coeffs))
                        and abs(float(q_coeffs[0])) <= 0.25
                        and abs(float(q_coeffs[1])) <= s_max
                        and np.isfinite(sq)
                        and (not np.isfinite(sc_f) or sq <= sc_f + 1e-6)
                    ):
                        ct2_f = float(q_coeffs[0])
                        ct_f = float(q_coeffs[1])
                        zp_f = float(q_coeffs[2])
                        Nf, sc_f = int(Nq), float(sq)

                clabel = color_col_name.replace("_", "-")
                axis_tag = "_std" if color_axis == "standard" else "_inst"
                quad_str = f" {ct2_f:+.4f}*({clabel})^2" if ct2_f else ""
                self._log(f"{filt}_std = {filt}_inst + {zp_f:+.4f} + {ct_f:+.4f}*({clabel}){axis_tag}{quad_str}  N={Nf}  scatter={sc_f:.4f}")

                # Store delta and color columns for CSV. color_<pair> is the
                # FIT axis (standard color when available); the raw
                # instrumental pair is kept alongside for diagnostics.
                out_cal[f"delta_{filt}"] = delta
                if color_col_name != "none":
                    ccol = f"color_{color_col_name}"
                    if ccol not in out_cal.columns:
                        out_cal[ccol] = color_x
                    ccol_inst = f"color_{color_col_name}_inst"
                    if ccol_inst not in out_cal.columns:
                        ca, cb = color_col_name.split("_", 1)
                        out_cal[ccol_inst] = _color_pair(ca, cb)

                filter_provenance = provenance_by_filter.get(
                    filt, build_photometry_provenance()
                )
                coeff_rows.append({"filter": filt, "zp": zp_f, "ct": ct_f, "ct2": ct2_f, "N": Nf,
                                   "scatter_rms": sc_f, "color_col": color_col_name,
                                   "color_axis": color_axis,
                                   "ref_source": ref_source_map.get(filt, ""),
                                   "photometry_source": filter_provenance["source"],
                                   "mag_input_column": filter_provenance["mag_column"],
                                   "mag_error_input_column": filter_provenance["mag_error_column"]})
                # The color support of the calibrators. solve_standard_colors
                # clamps to it so the (quadratic) color term is never evaluated
                # where no calibrator constrained it.
                _cfit = color_x[m_fit]
                _cfit = _cfit[np.isfinite(_cfit)]
                fit_params[filt] = {"zp": zp_f, "ct": ct_f, "ct2": ct2_f, "scatter_rms": sc_f,
                                    "color_min": float(_cfit.min()) if _cfit.size else float("nan"),
                                    "color_max": float(_cfit.max()) if _cfit.size else float("nan"),
                                    "color_col": color_col_name,
                                    "color_axis": color_axis}

            if not fit_params:
                raise RuntimeError("ZP fit failed for all detected filters")

            coeff_df = pd.DataFrame(coeff_rows)
            coeff_df.to_csv(output_dir / "zp_fit_coefficients.csv", index=False)
            self._log("Saved zp_fit_coefficients.csv")
            apply_ext = bool(getattr(P, "cmd_apply_extinction", False))
            ext_mode = str(getattr(P, "cmd_extinction_mode", "absorb")).strip().lower()
            if ext_mode not in ("absorb", "two_step"):
                ext_mode = "absorb"
            min_frame_refs = int(getattr(P, "frame_zp_min_n", 5))

            _am_empty = pd.DataFrame(columns=["file", "filter", "airmass"])
            if apply_ext:
                frame_airmass = self._build_frame_airmass(idx)
                if frame_airmass is None or frame_airmass.empty:
                    frame_airmass = _am_empty
            else:
                frame_airmass = _am_empty

            ext_df = None
            ext_map = {}
            if apply_ext:
                ext_dir = tool_extinction_dir(result_dir)
                ext_path = ext_dir / "extinction_fit_by_filter.csv"
                if not ext_path.exists():
                    ext_path = result_dir / "extinction_fit_by_filter.csv"
                if ext_path.exists():
                    ext_df = pd.read_csv(ext_path)
                    if {"filter", "k"}.issubset(ext_df.columns):
                        for _, er in ext_df.iterrows():
                            ext_map[normalize_filter_name(er["filter"])] = float(er["k"])

            out_cal_path = output_dir / "gaia_sdss_calibrator_by_ID.csv"
            out_cal.to_csv(out_cal_path, index=False)
            self._log(f"Saved {out_cal_path.name} | rows={len(out_cal)}")

            # Color index DataFrame for all stars (used for color-term
            # application and per-frame ZP). The applied color is the exact
            # STANDARD color solved from each star's own instrumental
            # magnitudes + the fitted constants (solve_standard_colors) — not
            # the raw instrumental color, and not a catalog color, so no
            # external faint-end systematics can leak in.
            inst_mags_all = {
                f: wide_raw[f"mag_inst_{f}"].to_numpy(float)
                for f in fit_params
                if f"mag_inst_{f}" in wide_raw.columns
            }
            solved_colors = solve_standard_colors(inst_mags_all, fit_params)
            color_df = wide_raw[["ID"]].copy()
            for filt, fp in fit_params.items():
                ccol_name = fp["color_col"]
                if ccol_name == "none":
                    continue
                fa, fb = ccol_name.split("_", 1)
                ca, cb = f"mag_inst_{fa}", f"mag_inst_{fb}"
                col_out = f"color_{ccol_name}"
                if col_out in color_df.columns:
                    continue
                if ccol_name in solved_colors:
                    color_df[col_out] = solved_colors[ccol_name]
                    self._log(f"[ZP][{ccol_name}] applied color = standard color (joint solve from instrumental mags)")
                elif ca in wide_raw.columns and cb in wide_raw.columns:
                    # One band of the pair has no fit — legacy instrumental color.
                    color_df[col_out] = wide_raw[ca].to_numpy(float) - wide_raw[cb].to_numpy(float)
                    self._log(f"[ZP][{ccol_name}] applied color = instrumental (pair not fully fitted)")

            # Merge calibrator columns into per-observation table
            cal_merge_cols = ["ID"]
            for filt in fit_params:
                rc = ref_col_map.get(filt)
                if rc and rc in out_cal.columns and rc not in cal_merge_cols:
                    cal_merge_cols.append(rc)
            for ccol in [c for c in out_cal.columns if c.startswith("color_")]:
                if ccol not in cal_merge_cols:
                    cal_merge_cols.append(ccol)
            if "gaia_BP_RP" in out_cal.columns:
                cal_merge_cols.append("gaia_BP_RP")

            obs = all_df.merge(out_cal[cal_merge_cols], on="ID", how="left")
            obs = obs.merge(frame_airmass[["file", "filter", "airmass"]], left_on=["file", "FILTER"], right_on=["file", "filter"], how="left")

            obs["ref_mag"] = np.nan
            obs["color_term"] = np.nan
            for filt, fp in fit_params.items():
                m_f = obs["FILTER"] == filt
                rc = ref_col_map.get(filt)
                if rc and rc in obs.columns:
                    obs.loc[m_f, "ref_mag"] = obs.loc[m_f, rc]
                ccol_name = fp["color_col"]
                if ccol_name != "none":
                    ccol = f"color_{ccol_name}"
                    if ccol in obs.columns:
                        _cvals = obs.loc[m_f, ccol].to_numpy(float)
                        obs.loc[m_f, "color_term"] = (
                            fp["ct"] * _cvals + float(fp.get("ct2", 0.0)) * _cvals * _cvals
                        )
                    else:
                        obs.loc[m_f, "color_term"] = 0.0
                else:
                    obs.loc[m_f, "color_term"] = 0.0

            obs["k_term"] = 0.0
            if apply_ext and ext_map:
                for f, k in ext_map.items():
                    m = obs["FILTER"] == f
                    obs.loc[m, "k_term"] = k * obs.loc[m, "airmass"].to_numpy(float)

            obs["delta"] = obs["ref_mag"] - (obs["mag_inst"] + obs["color_term"])
            if apply_ext and ext_mode == "two_step":
                obs["delta"] = obs["ref_mag"] - (obs["mag_inst"] + obs["color_term"] + obs["k_term"])

            obs["snr_ok"] = True
            if "snr" in obs.columns:
                svals = obs["snr"].to_numpy(float)
                obs["snr_ok"] = np.isfinite(svals) & (svals >= snr_cut)

            obs["cal_ok"] = False
            bp = obs["gaia_BP_RP"].to_numpy(float) if "gaia_BP_RP" in obs.columns else np.full(len(obs), np.nan)
            m_bp_global = np.isfinite(bp) & (bp >= bpRP_lo) & (bp <= bpRP_hi)
            for filt in fit_params:
                key = _BAND_ALIASES.get(filt, filt)
                entry = _GAIA_TO_BAND.get(key)
                if entry:
                    _, lo, hi, _, _ = entry
                    obs.loc[(obs["FILTER"] == filt) & m_bp_global & (bp >= lo) & (bp <= hi), "cal_ok"] = True

            obs_all = obs.copy()
            obs = obs_all[np.isfinite(obs_all["delta"]) & obs_all["snr_ok"] & obs_all["cal_ok"]].copy()

            if len(obs_all):
                summary_rows = []
                for filt, sub in obs_all.groupby("FILTER"):
                    ref_ok = np.isfinite(sub["ref_mag"].to_numpy(float))
                    delta_ok = np.isfinite(sub["delta"].to_numpy(float))
                    snr_ok = sub["snr_ok"].to_numpy(bool)
                    cal_ok = sub["cal_ok"].to_numpy(bool)
                    kept = delta_ok & snr_ok & cal_ok
                    summary_rows.append({
                        "filter": filt,
                        "n_total": int(len(sub)),
                        "n_ref_ok": int(np.sum(ref_ok)),
                        "n_delta_ok": int(np.sum(delta_ok)),
                        "n_snr_ok": int(np.sum(snr_ok)),
                        "n_cal_ok": int(np.sum(cal_ok)),
                        "n_kept": int(np.sum(kept)),
                    })
                cut_df = pd.DataFrame(summary_rows)
                cut_path = output_dir / "frame_zeropoint_cut_summary.csv"
                cut_df.to_csv(cut_path, index=False)
                self._log(f"Saved {cut_path.name} | rows={len(cut_df)}")

            frame_rows = []
            reject_rows = []
            for (fname, filt), sub in obs.groupby(["file", "FILTER"]):
                med, std, n, out_frac = self._robust_location(sub["delta"].to_numpy(float), clip_sigma=clip_sigma, iters=fit_iters)
                if n < min_frame_refs:
                    reject_rows.append({
                        "file": fname,
                        "filter": filt,
                        "n_ref": int(n),
                        "min_required": int(min_frame_refs),
                        "reason": "n_ref_below_min",
                    })
                    continue
                frame_rows.append({
                    "file": fname,
                    "filter": filt,
                    "zp_frame": med,
                    "zp_scatter": std,
                    "n_ref": n,
                    "outlier_fraction": out_frac,
                    "snr_med": float(np.nanmedian(sub["snr"].to_numpy(float))) if "snr" in sub.columns else np.nan,
                    "photometry_source": collapse_provenance_values(sub["photometry_source"]),
                    "mag_input_column": collapse_provenance_values(sub["mag_input_column"]),
                    "mag_error_input_column": collapse_provenance_values(sub["mag_error_input_column"]),
                })

            frame_df = pd.DataFrame(frame_rows)
            if frame_df.empty:
                frame_df = pd.DataFrame(columns=[
                    "file", "filter", "zp_frame", "zp_scatter", "n_ref",
                    "outlier_fraction", "snr_med", "photometry_source",
                    "mag_input_column", "mag_error_input_column",
                ])
                self._log("No per-frame ZP points; falling back to global ZP by filter.")
            if len(frame_df):
                frame_df = frame_df.merge(frame_airmass, on=["file", "filter"], how="left")
                frame_zp_path = output_dir / "frame_zeropoint.csv"
                frame_df.to_csv(frame_zp_path, index=False)
                self._log(f"Saved {frame_zp_path.name} | rows={len(frame_df)}")
            if reject_rows:
                reject_df = pd.DataFrame(reject_rows)
                reject_path = output_dir / "frame_zeropoint_rejects.csv"
                reject_df.to_csv(reject_path, index=False)
                self._log(f"Saved {reject_path.name} | rows={len(reject_df)}")

            obs = all_df.merge(frame_df[["file", "filter", "zp_frame"]], left_on=["file", "FILTER"], right_on=["file", "filter"], how="left")
            zp_map = {filt: fp["zp"] for filt, fp in fit_params.items()}
            obs["zp_frame"] = obs["zp_frame"].fillna(obs["FILTER"].map(zp_map))

            color_merge_cols = ["ID"] + [c for c in color_df.columns if c != "ID"]
            obs = obs.merge(color_df[color_merge_cols], on="ID", how="left")
            obs = obs.merge(frame_airmass[["file", "filter", "airmass"]], left_on=["file", "FILTER"], right_on=["file", "filter"], how="left")

            obs["color_term"] = np.nan
            for filt, fp in fit_params.items():
                m_f = obs["FILTER"] == filt
                ccol_name = fp["color_col"]
                if ccol_name != "none":
                    ccol = f"color_{ccol_name}"
                    if ccol in obs.columns:
                        _cvals = obs.loc[m_f, ccol].to_numpy(float)
                        obs.loc[m_f, "color_term"] = (
                            fp["ct"] * _cvals + float(fp.get("ct2", 0.0)) * _cvals * _cvals
                        )
                    else:
                        obs.loc[m_f, "color_term"] = 0.0
                else:
                    obs.loc[m_f, "color_term"] = 0.0

            obs["k_term"] = 0.0
            if apply_ext and ext_map:
                for f, k in ext_map.items():
                    m = obs["FILTER"] == f
                    obs.loc[m, "k_term"] = k * obs.loc[m, "airmass"].to_numpy(float)

            obs["mag_cal"] = obs["mag_inst"] + obs["zp_frame"] + obs["color_term"]
            if apply_ext and ext_mode == "two_step":
                obs["mag_cal"] = obs["mag_cal"] + obs["k_term"]

            def _combine_group_cal(g):
                med, med_err, n_med = self._robust_median_and_err(
                    g["mag_cal"], g["mag_err"])
                wmean, werr, _ = self._weighted_mean_mag(g["mag_cal"], g["mag_err"])
                _snr_col = "snr" if "snr" in g.columns else ("snr_psf" if "snr_psf" in g.columns else None)
                snr_vals = np.asarray(g[_snr_col], float) if _snr_col else np.array([np.nan])
                snr_med = float(np.nanmedian(snr_vals)) if np.isfinite(snr_vals).any() else np.nan
                return pd.Series({
                    "mag_cal_med": med,
                    "mag_cal_med_err": med_err,
                    "mag_cal_wmean": wmean,
                    "mag_cal_werr": werr,
                    "n_frames": n_med,
                    "snr_med": snr_med,
                    "photometry_source": collapse_provenance_values(g["photometry_source"]),
                    "mag_input_column": collapse_provenance_values(g["mag_input_column"]),
                    "mag_error_input_column": collapse_provenance_values(g["mag_error_input_column"]),
                })

            grp_cal = obs.groupby(["ID", "FILTER"], as_index=False).apply(_combine_group_cal)

            grp_path = output_dir / "median_by_ID_filter.csv"
            grp_cal.to_csv(grp_path, index=False, na_rep="NaN")
            self._log(f"Saved {grp_path.name} | rows={len(grp_cal)}")

            # Multi-exposure detector-linearity diagnostic (no-op for single exposure).
            self._write_nonlinearity_diag(obs, output_dir)

            wide_mag_w = grp_cal.pivot_table(index="ID", columns="FILTER", values="mag_cal_wmean", aggfunc="median")
            wide_err_w = grp_cal.pivot_table(index="ID", columns="FILTER", values="mag_cal_werr", aggfunc="median")
            wide_mag_med = grp_cal.pivot_table(index="ID", columns="FILTER", values="mag_cal_med", aggfunc="median")
            wide_err_med = grp_cal.pivot_table(index="ID", columns="FILTER", values="mag_cal_med_err", aggfunc="median")
            wide_mag = wide_mag_w.combine_first(wide_mag_med)
            wide_err = wide_err_w.combine_first(wide_err_med)
            wide_snr = grp_cal.pivot_table(index="ID", columns="FILTER", values="snr_med", aggfunc="median")
            wide_source = grp_cal.pivot_table(index="ID", columns="FILTER", values="photometry_source", aggfunc="first")
            wide_mag_input = grp_cal.pivot_table(index="ID", columns="FILTER", values="mag_input_column", aggfunc="first")
            wide_err_input = grp_cal.pivot_table(index="ID", columns="FILTER", values="mag_error_input_column", aggfunc="first")

            wide_mag.columns = [f"mag_cal_{c}" for c in wide_mag.columns]
            wide_err.columns = [f"mag_cal_err_{c}" for c in wide_err.columns]
            wide_snr.columns = [f"snr_{c}" for c in wide_snr.columns]
            wide_mag_w.columns = [f"mag_cal_wmean_{c}" for c in wide_mag_w.columns]
            wide_err_w.columns = [f"mag_cal_werr_{c}" for c in wide_err_w.columns]
            wide_mag_med.columns = [f"mag_cal_med_{c}" for c in wide_mag_med.columns]
            wide_err_med.columns = [f"mag_cal_med_err_{c}" for c in wide_err_med.columns]
            wide_source.columns = [f"photometry_source_{c}" for c in wide_source.columns]
            wide_mag_input.columns = [f"mag_input_column_{c}" for c in wide_mag_input.columns]
            wide_err_input.columns = [f"mag_error_input_column_{c}" for c in wide_err_input.columns]

            wide = pd.concat(
                [
                    wide_mag,
                    wide_err,
                    wide_mag_w,
                    wide_err_w,
                    wide_mag_med,
                    wide_err_med,
                    wide_snr,
                    wide_source,
                    wide_mag_input,
                    wide_err_input,
                ],
                axis=1,
            ).reset_index()

            wide_path = output_dir / "median_by_ID_filter_wide.csv"
            wide.to_csv(wide_path, index=False, na_rep="NaN")
            self._log(f"Saved {wide_path.name} | rows={len(wide)}")

            df_out = wide_raw.merge(master[merge_cols], on="ID", how="left")
            if col_xm:
                df_out = df_out.rename(columns={col_xm: "x_pix"})
            if col_ym:
                df_out = df_out.rename(columns={col_ym: "y_pix"})

            for col in ("gaia_G", "gaia_BP", "gaia_RP"):
                if col not in df_out.columns and col in df.columns:
                    df_out[col] = df[col]

            wide_by_id = wide.drop_duplicates(subset=["ID"], keep="first").set_index("ID", drop=False)
            cal_cols_added = []
            std_cols_added = []
            for filt in fit_params:
                c_cal = f"mag_cal_{filt}"
                c_cal_err = f"mag_cal_err_{filt}"
                c_std = f"mag_std_{filt}"
                c_std_err = f"mag_std_err_{filt}"
                if c_cal in wide.columns:
                    cal_values = pd.to_numeric(
                        wide_by_id.reindex(df_out["ID"])[c_cal],
                        errors="coerce",
                    ).to_numpy(float)
                    df_out[c_cal] = cal_values
                    df_out[c_std] = cal_values
                else:
                    df_out[c_cal] = np.nan
                    df_out[c_std] = np.nan
                if c_cal_err in wide.columns:
                    cal_err_values = pd.to_numeric(
                        wide_by_id.reindex(df_out["ID"])[c_cal_err],
                        errors="coerce",
                    ).to_numpy(float)
                    df_out[c_cal_err] = cal_err_values
                    df_out[c_std_err] = cal_err_values
                else:
                    df_out[c_cal_err] = np.nan
                    df_out[c_std_err] = np.nan
                cal_cols_added.append(c_cal)
                std_cols_added.append(c_std)

            if std_cols_added:
                n_missing = int(df_out[std_cols_added].isna().all(axis=1).sum())
                if n_missing:
                    self._log(f"CMD export: {n_missing} IDs missing all calibrated magnitudes; keeping mag_cal_* as NaN.")
                # Per-band finite counts. A CMD needs BOTH bands of its color
                # pair calibrated, so an all-NaN single band silently empties
                # the Std CMD even when n_missing is small (that count only
                # flags rows where EVERY band is NaN). Log each band so an empty
                # one is obvious.
                _band_counts = []
                for c_cal in cal_cols_added:
                    _nf = int(np.isfinite(df_out[c_cal].to_numpy(float)).sum())
                    _band_counts.append(f"{c_cal}={_nf}/{len(df_out)}")
                self._log("CMD export: calibrated-band finite counts | " + ", ".join(_band_counts))

            # External standard anchor must shift mag_std_* BEFORE anything is
            # derived from them — the synthetic Gaia columns below feed the
            # Gaia CMD QC, and computing them from pre-anchor magnitudes would
            # bias the QC drift metric by exactly the anchor offset.
            df_out = self._apply_standard_anchor(df_out, output_dir)

            # Synthetic Gaia magnitudes: prefer SDSS g-i, fall back to Johnson V-I
            gaia_G_syn     = np.full(len(df_out), np.nan)
            gaia_BP_RP_syn = np.full(len(df_out), np.nan)

            if "mag_std_g" in df_out.columns and "mag_std_i" in df_out.columns:
                gi_std  = df_out["mag_std_g"].to_numpy(float) - df_out["mag_std_i"].to_numpy(float)
                m_gi    = np.isfinite(gi_std) & (gi_std >= 1.0) & (gi_std <= 9.0)
                gaia_G_syn[m_gi]     = df_out["mag_std_g"].to_numpy(float)[m_gi] + self._poly_eval(gi_std[m_gi], [-0.1064, -0.4964, -0.09339, 0.004444])
                gaia_BP_RP_syn[m_gi] = self._poly_eval(gi_std[m_gi], [0.3971, 0.777, -0.04164, 0.008237])
            elif "mag_std_V" in df_out.columns and "mag_std_I" in df_out.columns:
                # Johnson V-I: GBP-GRP ≈ f(V-I) from Riello+2021 inverse
                vi_std  = df_out["mag_std_V"].to_numpy(float) - df_out["mag_std_I"].to_numpy(float)
                m_vi    = np.isfinite(vi_std) & (vi_std >= -0.5) & (vi_std <= 5.0)
                bprp    = self._poly_eval(vi_std[m_vi], [-0.033, 1.259, -0.128, 0.016])
                gm_v    = self._poly_eval(bprp, [-0.02704, 0.01424, -0.2156, 0.01426])
                gaia_G_syn[m_vi]     = df_out["mag_std_V"].to_numpy(float)[m_vi] + gm_v
                gaia_BP_RP_syn[m_vi] = bprp

            df_out["gaia_G_syn"]     = gaia_G_syn
            df_out["gaia_BP_RP_syn"] = gaia_BP_RP_syn

            out_cmd_path = output_dir / "median_by_ID_filter_wide_cmd.csv"
            df_out.to_csv(out_cmd_path, index=False, na_rep="NaN")
            self._log(f"Saved {out_cmd_path.name} | rows={len(df_out)}")
            export_zp_qc_products(output_dir, self._log)
            export_cmd_qc_products(output_dir, self._log)
            export_gaia_cmd_comparison_products(output_dir, self._log)

            summary = {
                "ok": True,
                "wide": str(wide_path),
                "cmd": str(out_cmd_path),
                "frame_airmass": str((output_dir / "frame_airmass.csv")) if (output_dir / "frame_airmass.csv").exists() else "",
                "frame_zeropoint": str((output_dir / "frame_zeropoint.csv")) if (output_dir / "frame_zeropoint.csv").exists() else "",
            }
            summary.update(summarize_photometry_table(df_out))
            self.last_summary = dict(summary)
            self.last_error = ""
            self.finished.emit(summary)
        except Exception as e:
            error_msg = f"{type(e).__name__}: {str(e)}\n{traceback.format_exc()}"
            self.last_error = error_msg
            self.last_summary = {}
            self.error.emit(error_msg)


class CmdViewerWindow(QWidget):
    """Interactive CMD viewer (Qt)."""

    def __init__(self, df: pd.DataFrame, result_dir: Path, parent=None, embedded: bool = False, params=None):
        super().__init__(parent)
        self.df = self._with_calibrated_aliases(df)
        self.photometry_provenance = summarize_photometry_table(self.df)
        self.result_dir = Path(result_dir)
        self.params = params

        self.setWindowTitle("CMD Viewer")
        if embedded:
            self.setWindowFlags(Qt.Widget)
            # Embedded in Step 11, this is a *preferred* size, not a floor:
            # a 600 px minimum exceeded the room Step 11 had, so the host
            # scroll showed only ~150 px of the plot at a time. Keep it low
            # enough to fit a laptop screen and let the canvas expand.
            self.setMinimumSize(760, 420)
        else:
            self.setWindowFlag(Qt.Window, True)
            self.resize(1200, 900)
            self.setMinimumSize(1000, 720)

        # View mode is selected after available magnitude products are detected.
        self.view_mode = 0

        self.inst_bands = _filter_bands_from_columns(self.df.columns, "mag_inst_")
        self.std_bands  = _filter_bands_from_columns(self.df.columns, "mag_std_")
        std_value_cols = [f"mag_std_{band}" for band in self.std_bands if f"mag_std_{band}" in self.df.columns]
        self.has_std = bool(std_value_cols) and np.isfinite(self.df[std_value_cols].to_numpy(float)).any()

        self.all_bands = sorted(set(self.inst_bands) | set(self.std_bands))

        # X axis: adjacent color pairs only (e.g. B-V, V-R — standard CMD indices)
        # Y axis: scalar magnitudes only (CMD viewer convention)
        axis_bands = self.inst_bands or self.std_bands
        x_allowed = _build_color_pairs(axis_bands, adjacent_only=True)
        self.x_allowed         = x_allowed
        self.y_allowed_scalars = axis_bands  # already wavelength-sorted
        self.y_allowed_colors  = []

        self.x_pairs       = x_allowed
        self.y_scalar_opts = axis_bands
        self.y_color_pairs = []

        self.snr_cols = [c for c in self.df.columns if c.startswith("snr_")]
        self.has_snr = len(self.snr_cols) > 0

        self.has_gaia_inst = (
            {"gaia_G_inst", "gaia_BP_RP_inst"}.issubset(df.columns)
            and np.isfinite(df["gaia_G_inst"].to_numpy(float)).any()
            and np.isfinite(df["gaia_BP_RP_inst"].to_numpy(float)).any()
        )
        self.has_gaia_syn = (
            {"gaia_G_syn", "gaia_BP_RP_syn"}.issubset(df.columns)
            and np.isfinite(df["gaia_G_syn"].to_numpy(float)).any()
            and np.isfinite(df["gaia_BP_RP_syn"].to_numpy(float)).any()
        )
        # Gaia CMD is a diagnostic-only view.  Keep Gaia columns available for
        # membership and click details, but do not add a third CMD panel.
        self.gaia_mode = None

        self.teff_vmin = 2400.0
        self.teff_vmax = 40000.0
        self.ob_norm = Normalize(vmin=self.teff_vmin, vmax=self.teff_vmax, clip=True)

        anchors = [
            (2400, "#E53935"),
            (3200, "#FF6A3D"),
            (4500, "#FFB84D"),
            (5800, "#FFE36A"),
            (6500, "#FFF6C7"),
            (8000, "#FFFFFF"),
            (10000, "#FFFFFF"),
            (20000, "#2D5BFF"),
            (40000, "#7A3CFF"),
        ]
        anchors = sorted(anchors, key=lambda x: x[0])
        pos = [(t - self.teff_vmin) / (self.teff_vmax - self.teff_vmin) for t, _ in anchors]
        pos[0] = 0.0
        pos[-1] = 1.0

        self.ob_cmap = LinearSegmentedColormap.from_list(
            "obafgkm_like",
            list(zip(pos, [c for _, c in anchors])),
            N=256
        )
        self.ob_cmap.set_bad("#777777")

        self.color_anchors = _TEFF_COLOR_ANCHORS

        # Determine available views
        self.available_views = []
        if self.inst_bands:
            self.available_views.append("inst")
        if self.has_std:
            self.available_views.append("std")
        if self.gaia_mode is not None:
            self.available_views.append("gaia")
        if len(self.available_views) > 1:
            self.available_views.append("all")
        if not self.available_views:
            self.available_views = ["inst"]
        self.view_mode = self.available_views.index("std") if "std" in self.available_views else 0

        self._plot_cache = {}
        self.last_pick_info = None
        self.pick_log = []
        self._membership_prob = None
        self._membership_source = "none"
        self._membership_note = ""
        self._membership_ready = False
        self._parallax_range_initialized = False
        self._roi_data: dict | None = None
        self._build_ui()
        self._update_view_label()
        self._load_roi()
        self._initialize_parallax_range(force=True)
        self._build_figure()
        self._redraw()
        self.setFocusPolicy(Qt.StrongFocus)

    @staticmethod
    def _with_calibrated_aliases(df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        for band in _filter_bands_from_columns(out.columns, "mag_cal_"):
            cal_col = f"mag_cal_{band}"
            std_col = f"mag_std_{band}"
            if std_col not in out.columns and cal_col in out.columns:
                out[std_col] = out[cal_col]
            cal_err_col = f"mag_cal_err_{band}"
            std_err_col = f"mag_std_err_{band}"
            if std_err_col not in out.columns and cal_err_col in out.columns:
                out[std_err_col] = out[cal_err_col]
        return out

    def _update_view_label(self) -> None:
        view_name = self.available_views[self.view_mode] if self.available_views else "inst"
        view_labels = {"inst": "Instrumental", "std": "Calibrated", "gaia": "Gaia", "all": "All CMDs"}
        if hasattr(self, "view_label"):
            self.view_label.setText(f"View: {view_labels.get(view_name, view_name)}")

    def _build_ui(self):
        layout = QVBoxLayout(self)

        # Top controls split across two rows so labels don't get clipped
        # in narrow windows.  Row 1 = selection axes & filters; Row 2 =
        # export buttons + current view label.
        controls = QHBoxLayout()
        controls.addWidget(QLabel("X(color):"))
        self.x_combo = QComboBox()
        self.x_combo.addItems([f"{a}-{b}" for (a, b) in self.x_pairs] or ["(none)"])
        controls.addWidget(self.x_combo)

        controls.addWidget(QLabel("Y:"))
        y_opts = self.y_scalar_opts + [f"{a}-{b}" for (a, b) in self.y_color_pairs]
        self.y_combo = QComboBox()
        self.y_combo.addItems(y_opts or ["(none)"])
        controls.addWidget(self.y_combo)

        controls.addWidget(QLabel("SNR >="))
        self.snr_spin = QSpinBox()
        self.snr_spin.setRange(0, 100)
        self.snr_spin.setValue(20)
        controls.addWidget(self.snr_spin)

        self.invert_y = QCheckBox("Invert Y")
        self.invert_y.setChecked(True)
        controls.addWidget(self.invert_y)

        # Optional extra ZP nudge for the Instrumental view. Instrumental
        # magnitudes already carry the IRAF Z=25 convention baked in at Step 7
        # (mag_inst = 25 - 2.5*log10(flux_e/exptime)), so they read in the
        # usual positive range without any shift — this control defaults to 0
        # and is only for manual fine-tuning. Colors (X = a-b) are unaffected
        # because a constant ZP cancels in a difference.
        self.manual_zp_check = QCheckBox("Manual ZP")
        self.manual_zp_check.setToolTip(
            "Add an extra constant zeropoint to Instrumental magnitudes for display only.\n"
            "mag_inst already includes the IRAF Z=25 convention, so leave at 0 normally.\n"
            "Colors are unchanged."
        )
        controls.addWidget(self.manual_zp_check)
        self.manual_zp_spin = QDoubleSpinBox()
        self.manual_zp_spin.setRange(0.0, 50.0)
        self.manual_zp_spin.setDecimals(3)
        self.manual_zp_spin.setSingleStep(0.1)
        self.manual_zp_spin.setValue(0.0)
        self.manual_zp_spin.setToolTip("Extra Instrumental-view zeropoint added to Y (display only).")
        controls.addWidget(self.manual_zp_spin)

        controls.addSpacing(8)
        controls.addWidget(QLabel("Membership:"))
        self.member_mode_combo = QComboBox()
        self.member_mode_combo.addItems([
            "Off",
            "Loose (P>=0.30)",
            "Normal (P>=0.50)",
            "Strict (P>=0.80)",
        ])
        mode_raw = "off"
        if self.params is not None and hasattr(self.params, "P"):
            mode_raw = str(getattr(self.params.P, "cmd_membership_mode", "off")).strip().lower()
        mode_to_idx = {"off": 0, "loose": 1, "normal": 2, "strict": 3}
        self.member_mode_combo.setCurrentIndex(mode_to_idx.get(mode_raw, 0))
        controls.addWidget(self.member_mode_combo)

        self.member_compare = QCheckBox("Compare")
        cmp_default = True
        if self.params is not None and hasattr(self.params, "P"):
            cmp_default = bool(getattr(self.params.P, "cmd_membership_compare", True))
        self.member_compare.setChecked(cmp_default)
        controls.addWidget(self.member_compare)
        controls.addStretch()
        layout.addLayout(controls)

        controls_row2 = QHBoxLayout()
        self.btn_save_membership = QPushButton("Save Pmem CSV")
        controls_row2.addWidget(self.btn_save_membership)

        self.btn_reset_filters = QPushButton("Reset View Filters")
        self.btn_reset_filters.setToolTip("Restore CMD viewer filters to the project defaults.")
        controls_row2.addWidget(self.btn_reset_filters)

        self.save_btn = QPushButton("Save PNG")
        controls_row2.addWidget(self.save_btn)

        self.xerr_check = QCheckBox("X err")
        self.xerr_check.setToolTip("Show color-index error bars for foreground CMD points.")
        self.xerr_check.setChecked(False)
        controls_row2.addWidget(self.xerr_check)

        self.yerr_check = QCheckBox("Y err")
        self.yerr_check.setToolTip("Show magnitude error bars for foreground CMD points.")
        self.yerr_check.setChecked(False)
        controls_row2.addWidget(self.yerr_check)

        controls_row2.addStretch()
        self.photometry_source_label = QLabel(
            format_photometry_provenance(self.photometry_provenance)
        )
        self.photometry_source_label.setProperty("role", "caption")
        _f = self.photometry_source_label.font(); _f.setBold(True)
        self.photometry_source_label.setFont(_f)
        controls_row2.addWidget(self.photometry_source_label)

        self.view_label = QLabel("View: Instrumental")
        self.view_label.setProperty("role", "info")
        controls_row2.addWidget(self.view_label)
        layout.addLayout(controls_row2)

        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setFixedHeight(90)
        self.info_text.setObjectName("Log")     # themed mono surface
        layout.addWidget(self.info_text)

        # 10×6 instead of 12×6 — the wider aspect made the CMD look stretched
        # on 1400px+ windows.
        self.figure = Figure(figsize=(10, 6), dpi=100)
        self.figure.subplots_adjust(bottom=0.14)
        self.canvas = FigureCanvas(self.figure)
        # Canvas is the only stretch=1 widget below; let it absorb ALL
        # spare vertical space.  (A maxHeight here caused leftover space
        # to be distributed as ugly gaps between the control rows in the
        # standalone CMD+Isochrone tool.)
        self.canvas.setMinimumSize(640, 300)

        # Prev/Next View are mounted next to the matplotlib navigation
        # toolbar (above the canvas) so they cannot visually collide with
        # the X-axis label below.
        self.btn_prev_view = QPushButton("Prev View")
        self.btn_next_view = QPushButton("Next View")
        self.btn_prev_view.setFixedHeight(28)
        self.btn_next_view.setFixedHeight(28)
        self.toolbar = NavigationToolbar(self.canvas, self)
        toolbar_row = QHBoxLayout()
        toolbar_row.setContentsMargins(0, 0, 0, 0)
        toolbar_row.addWidget(self.toolbar)
        toolbar_row.addStretch()
        toolbar_row.addWidget(self.btn_prev_view)
        toolbar_row.addWidget(self.btn_next_view)
        layout.addLayout(toolbar_row)
        layout.addWidget(self.canvas, stretch=1)

        self.x_combo.currentTextChanged.connect(self._redraw)
        self.y_combo.currentTextChanged.connect(self._redraw)
        self.snr_spin.valueChanged.connect(self._redraw)
        self.invert_y.stateChanged.connect(self._redraw)
        self.manual_zp_check.stateChanged.connect(self._redraw)
        self.manual_zp_spin.valueChanged.connect(self._redraw)
        self.xerr_check.stateChanged.connect(self._redraw)
        self.yerr_check.stateChanged.connect(self._redraw)
        self.member_mode_combo.currentIndexChanged.connect(self._on_membership_ui_changed)
        self.member_compare.stateChanged.connect(self._on_membership_ui_changed)
        self.btn_save_membership.clicked.connect(self._save_membership_csv)
        self.btn_reset_filters.clicked.connect(self._reset_view_filters)
        self.save_btn.clicked.connect(self._save_png)
        self.btn_prev_view.clicked.connect(lambda: self._switch_view(-1))
        self.btn_next_view.clicked.connect(lambda: self._switch_view(1))
        self.canvas.mpl_connect("button_press_event", self._on_plot_click)

        controls2 = QHBoxLayout()
        self.plx_check = QCheckBox("Parallax filter")
        self.plx_check.setChecked(False)
        controls2.addWidget(self.plx_check)
        controls2.addWidget(QLabel("min:"))
        self.plx_min_spin = QDoubleSpinBox()
        self.plx_min_spin.setRange(-5.0, 20.0)
        self.plx_min_spin.setDecimals(3)
        self.plx_min_spin.setSingleStep(0.05)
        self.plx_min_spin.setValue(-0.5)
        self.plx_min_spin.setSuffix(" mas")
        controls2.addWidget(self.plx_min_spin)
        controls2.addWidget(QLabel("max:"))
        self.plx_max_spin = QDoubleSpinBox()
        self.plx_max_spin.setRange(-5.0, 20.0)
        self.plx_max_spin.setDecimals(3)
        self.plx_max_spin.setSingleStep(0.05)
        self.plx_max_spin.setValue(0.5)
        self.plx_max_spin.setSuffix(" mas")
        controls2.addWidget(self.plx_max_spin)
        controls2.addSpacing(16)
        self.roi_check = QCheckBox("ROI filter")
        self.roi_check.setChecked(False)
        self.roi_check.setToolTip("Filter CMD sources by the spatial ROI circle set in Step 9.\nDoes not affect ZP calibration.")
        controls2.addWidget(self.roi_check)
        self.roi_info_label = QLabel("(no ROI)")
        self.roi_info_label.setProperty("role", "caption")
        controls2.addWidget(self.roi_info_label)
        self.btn_reload_roi = QPushButton("Reload")
        # sizeHint, not a hard 56 px: with the theme's button padding that
        # width clipped the label to "elo:".
        self.btn_reload_roi.setFixedWidth(self.btn_reload_roi.sizeHint().width())
        self.btn_reload_roi.setToolTip("Re-read cmd_roi.json from Step 9 output directory")
        controls2.addWidget(self.btn_reload_roi)
        controls2.addSpacing(16)

        # Gaia astrometric/photometric quality, as *display* filters. These are
        # the same two cuts step10 can apply when choosing zero-point
        # calibrators, offered here so the CMD can be inspected with and
        # without them — in a globular the C* cut alone removes about half the
        # Gaia-matched sources, and that is worth seeing rather than assuming.
        # Filtering here never touches the calibration; it only changes what
        # is drawn.
        self.ruwe_check = QCheckBox("RUWE ≤")
        self.ruwe_check.setChecked(False)
        self.ruwe_check.setToolTip(
            "Hide sources whose Gaia astrometric fit is poor (RUWE above the\n"
            "threshold) — typically unresolved binaries and blends.\n"
            "Display only; does not affect ZP calibration.")
        controls2.addWidget(self.ruwe_check)
        self.ruwe_spin = QDoubleSpinBox()
        self.ruwe_spin.setRange(1.0, 10.0)
        self.ruwe_spin.setDecimals(2)
        self.ruwe_spin.setSingleStep(0.1)
        self.ruwe_spin.setValue(1.4)
        self.ruwe_spin.setFixedWidth(self.ruwe_spin.sizeHint().width())
        controls2.addWidget(self.ruwe_spin)

        self.cstar_check = QCheckBox("C* ≤")
        self.cstar_check.setChecked(False)
        self.cstar_check.setToolTip(
            "Hide sources whose BP/RP flux excess is inconsistent with G by\n"
            "more than N sigma (Riello+2021) — BP/RP window contamination in\n"
            "crowded fields. Needs phot_bp_rp_excess_factor in the CMD table.\n"
            "Display only; does not affect ZP calibration.")
        controls2.addWidget(self.cstar_check)
        self.cstar_spin = QDoubleSpinBox()
        self.cstar_spin.setRange(1.0, 10.0)
        self.cstar_spin.setDecimals(1)
        self.cstar_spin.setSingleStep(0.5)
        self.cstar_spin.setValue(3.0)
        self.cstar_spin.setSuffix(" σ")
        self.cstar_spin.setFixedWidth(self.cstar_spin.sizeHint().width())
        controls2.addWidget(self.cstar_spin)

        self.quality_info_label = QLabel("")
        self.quality_info_label.setProperty("role", "caption")
        controls2.addWidget(self.quality_info_label)

        controls2.addStretch()
        layout.addLayout(controls2)

        self.plx_check.stateChanged.connect(self._on_plx_filter_changed)
        self.plx_min_spin.valueChanged.connect(self._redraw)
        self.plx_max_spin.valueChanged.connect(self._redraw)
        self.roi_check.stateChanged.connect(self._redraw)
        self.ruwe_check.stateChanged.connect(self._on_quality_filter_changed)
        self.ruwe_spin.valueChanged.connect(self._redraw)
        self.cstar_check.stateChanged.connect(self._on_quality_filter_changed)
        self.cstar_spin.valueChanged.connect(self._redraw)
        self.btn_reload_roi.clicked.connect(self._on_reload_roi)
        install_parameter_wheel_guard(self)

    def _reset_view_filters(self):
        for widget in (
            self.snr_spin,
            self.invert_y,
            self.member_mode_combo,
            self.member_compare,
            self.xerr_check,
            self.yerr_check,
            self.plx_check,
            self.plx_min_spin,
            self.plx_max_spin,
            self.roi_check,
        ):
            widget.blockSignals(True)

        self.snr_spin.setValue(20)
        self.invert_y.setChecked(True)
        self.member_mode_combo.setCurrentIndex(2)  # Normal (P>=0.50)
        self.member_compare.setChecked(True)
        self.xerr_check.setChecked(False)
        self.yerr_check.setChecked(False)
        self.plx_check.setChecked(False)
        self.plx_min_spin.setValue(-0.5)
        self.plx_max_spin.setValue(0.5)
        self._initialize_parallax_range(force=True)
        self.roi_check.setChecked(False)

        for widget in (
            self.snr_spin,
            self.invert_y,
            self.member_mode_combo,
            self.member_compare,
            self.xerr_check,
            self.yerr_check,
            self.plx_check,
            self.plx_min_spin,
            self.plx_max_spin,
            self.roi_check,
        ):
            widget.blockSignals(False)

        if self.params is not None and hasattr(self.params, "P"):
            self.params.P.cmd_membership_mode = "normal"
            self.params.P.cmd_membership_compare = True
            if hasattr(self.params, "save_toml"):
                try:
                    self.params.save_toml()
                except Exception:
                    pass
        self._redraw()

    def _load_roi(self):
        """Load cmd_roi.json from step8 output directory and update UI."""
        roi_path = step9_selection_dir(self.result_dir) / "cmd_roi.json"
        try:
            if roi_path.exists():
                self._roi_data = json.loads(roi_path.read_text())
            else:
                self._roi_data = None
        except Exception:
            self._roi_data = None
        if hasattr(self, "roi_check"):
            self.roi_check.setEnabled(self._roi_data is not None)
        if hasattr(self, "roi_info_label"):
            if self._roi_data:
                ra = self._roi_data.get("ra_deg", 0.0)
                dec = self._roi_data.get("dec_deg", 0.0)
                r = self._roi_data.get("radius_arcsec", 0.0)
                self.roi_info_label.setText(f"RA={ra:.4f} Dec={dec:.4f}  r={r:.0f}\"")
                _set_label_role(self.roi_info_label, "status", "ok")
            else:
                self.roi_info_label.setText("(no ROI)")
                _set_label_role(self.roi_info_label, "role", "caption")

    def _on_reload_roi(self):
        self._load_roi()
        self._redraw()

    def _parallax_values(self):
        if "parallax" not in self.df.columns:
            self._ensure_membership_columns_from_master()
        if "parallax" not in self.df.columns:
            return None
        return pd.to_numeric(self.df["parallax"], errors="coerce").to_numpy(float)

    def _set_parallax_range(self, plx_min: float, plx_max: float):
        if not np.isfinite(plx_min) or not np.isfinite(plx_max):
            return
        if plx_min > plx_max:
            plx_min, plx_max = plx_max, plx_min
        lo = max(float(self.plx_min_spin.minimum()), float(plx_min))
        hi = min(float(self.plx_max_spin.maximum()), float(plx_max))
        if lo >= hi:
            pad = max(0.05, abs(float(plx_min)) * 0.05)
            lo = max(float(self.plx_min_spin.minimum()), float(plx_min) - pad)
            hi = min(float(self.plx_max_spin.maximum()), float(plx_max) + pad)
        for widget in (self.plx_min_spin, self.plx_max_spin):
            widget.blockSignals(True)
        self.plx_min_spin.setValue(lo)
        self.plx_max_spin.setValue(hi)
        for widget in (self.plx_min_spin, self.plx_max_spin):
            widget.blockSignals(False)

    def _auto_set_parallax_range(self, preferred_mask=None) -> bool:
        plx = self._parallax_values()
        if plx is None:
            return False
        finite = np.isfinite(plx)
        base = finite.copy()
        if preferred_mask is not None and len(preferred_mask) == len(plx):
            preferred = np.asarray(preferred_mask, bool) & finite
            if int(preferred.sum()) >= 10:
                base = preferred
        if int(base.sum()) == 0:
            return False

        vals = plx[base]
        center = float(np.nanmedian(vals))
        mad = float(np.nanmedian(np.abs(vals - center)))
        robust_sigma = MAD_TO_SIGMA * mad if np.isfinite(mad) and mad > 0 else 0.0
        half_width = max(0.5, 4.0 * robust_sigma)
        half_width = min(5.0, half_width)
        self._set_parallax_range(center - half_width, center + half_width)
        return True

    def _initialize_parallax_range(self, force: bool = False) -> bool:
        if self._parallax_range_initialized and not force:
            return False
        if self._auto_set_parallax_range():
            self._parallax_range_initialized = True
            return True
        return False

    def _roi_mask(self):
        """Returns boolean mask selecting sources inside the CMD ROI circle (sky coords), or None if disabled."""
        if not self.roi_check.isChecked() or self._roi_data is None:
            return None
        roi_ra = float(self._roi_data["ra_deg"])
        roi_dec = float(self._roi_data["dec_deg"])
        roi_r_arcsec = float(self._roi_data["radius_arcsec"])
        # Prefer RA/Dec angular distance (correct across frames)
        if "ra_deg" in self.df.columns and "dec_deg" in self.df.columns:
            ra = pd.to_numeric(self.df["ra_deg"], errors="coerce").to_numpy(float)
            dec = pd.to_numeric(self.df["dec_deg"], errors="coerce").to_numpy(float)
            valid = np.isfinite(ra) & np.isfinite(dec)
            # Small-field approximation (accurate to ~0.01% within 1 deg)
            cos_dec = np.cos(np.radians(roi_dec))
            d_ra = (ra - roi_ra) * cos_dec * 3600.0   # arcsec
            d_dec = (dec - roi_dec) * 3600.0           # arcsec
            return valid & (d_ra ** 2 + d_dec ** 2 <= roi_r_arcsec ** 2)
        return None

    def _on_plx_filter_changed(self):
        if self.plx_check.isChecked():
            plx = self._parallax_values()
            if plx is None:
                self.plx_check.blockSignals(True)
                self.plx_check.setChecked(False)
                self.plx_check.blockSignals(False)
                from PyQt5.QtWidgets import QMessageBox
                QMessageBox.warning(self, "Parallax Unavailable",
                    "parallax column not found in CMD data or master_catalog.\n"
                    "Rerun Step 6 (Master Catalog Build) and Step 10 (ZP Calibration).")
                return
            plx_min = float(self.plx_min_spin.value())
            plx_max = float(self.plx_max_spin.value())
            finite = np.isfinite(plx)
            selected = finite & (plx >= plx_min) & (plx <= plx_max)
            if int(finite.sum()) > 0 and int(selected.sum()) < 5:
                member_mask, _, _ = self._current_membership_mask()
                self._auto_set_parallax_range(member_mask)
        self._redraw()

    def _parallax_mask(self):
        """Returns boolean mask for parallax range filter, or None if disabled/unavailable."""
        if not self.plx_check.isChecked():
            return None
        plx = self._parallax_values()
        if plx is None:
            return None
        plx_min = float(self.plx_min_spin.value())
        plx_max = float(self.plx_max_spin.value())
        mask = np.isfinite(plx) & (plx >= plx_min) & (plx <= plx_max)
        n_finite = int(np.isfinite(plx).sum())
        if n_finite == 0:
            return None
        return mask

    def _quality_column(self, name: str):
        """A Gaia quality column from the CMD table, as floats, or None."""
        df = getattr(self, "df", None)
        if df is None or name not in df.columns:
            return None
        vals = pd.to_numeric(df[name], errors="coerce").to_numpy(dtype=float)
        return vals if np.isfinite(vals).any() else None

    def _quality_mask(self):
        """Combined RUWE / C* display filter, or None when neither is on.

        Missing values are KEPT, matching `gaia_quality_mask`: a star with no
        RUWE has not failed the cut, it was never measured. Rejecting those
        would quietly drop every non-Gaia source from the plot.
        """
        mask = None
        n_before = None
        parts: list[str] = []

        if self.ruwe_check.isChecked():
            ruwe = self._quality_column("ruwe")
            if ruwe is not None:
                thr = float(self.ruwe_spin.value())
                keep = ~(np.isfinite(ruwe) & (ruwe > thr))
                n_before = keep.size
                mask = keep if mask is None else (mask & keep)
                parts.append(f"RUWE −{int((~keep).sum())}")
            else:
                parts.append("RUWE n/a")

        if self.cstar_check.isChecked():
            excess = self._quality_column("phot_bp_rp_excess_factor")
            bp_rp = self._quality_column("gaia_BP_RP")
            if bp_rp is None:
                bp = self._quality_column("gaia_BP")
                rp = self._quality_column("gaia_RP")
                bp_rp = (bp - rp) if (bp is not None and rp is not None) else None
            gmag = self._quality_column("gaia_G")
            if excess is not None and bp_rp is not None and gmag is not None:
                cstar = gaia_corrected_excess_factor(bp_rp, excess)
                sigma = gaia_cstar_sigma(gmag)
                nsig = float(self.cstar_spin.value())
                bad = (np.isfinite(cstar) & np.isfinite(sigma)
                       & (np.abs(cstar) > nsig * sigma))
                keep = ~bad
                n_before = keep.size
                mask = keep if mask is None else (mask & keep)
                parts.append(f"C* −{int(bad.sum())}")
            else:
                parts.append("C* n/a (no excess factor)")

        label = getattr(self, "quality_info_label", None)
        if label is not None:
            if mask is not None and n_before:
                label.setText(f"({' · '.join(parts)} → {int(mask.sum())}/{n_before})")
            else:
                label.setText(f"({' · '.join(parts)})" if parts else "")
        return mask

    def _on_quality_filter_changed(self):
        """Warn once when a cut is asked for but its column is not there."""
        if self.cstar_check.isChecked() and \
                self._quality_column("phot_bp_rp_excess_factor") is None:
            self.cstar_check.blockSignals(True)
            self.cstar_check.setChecked(False)
            self.cstar_check.blockSignals(False)
            from PyQt5.QtWidgets import QMessageBox
            QMessageBox.information(
                self, "C* Unavailable",
                "phot_bp_rp_excess_factor is not in this CMD table.\n\n"
                "Catalogues fetched before 2026-08-11 did not request it. "
                "Re-run Step 5 (Gaia query) and Step 6 (Master Catalog Build) "
                "to make the C* filter available.")
        self._redraw()

    def _membership_mode_key(self) -> str:
        idx = int(self.member_mode_combo.currentIndex())
        if idx == 1:
            return "loose"
        if idx == 2:
            return "normal"
        if idx == 3:
            return "strict"
        return "off"

    def _membership_threshold(self, mode_key: str) -> float:
        if mode_key == "loose":
            return 0.30
        if mode_key == "strict":
            return 0.80
        return 0.50

    def _on_membership_ui_changed(self):
        if self.params is not None and hasattr(self.params, "P"):
            self.params.P.cmd_membership_mode = self._membership_mode_key()
            self.params.P.cmd_membership_compare = bool(self.member_compare.isChecked())
            if hasattr(self.params, "save_toml"):
                try:
                    self.params.save_toml()
                except Exception:
                    pass
        self._redraw()

    def _pick_existing_membership_col(self):
        cands = [
            "gaia_pmem",
            "pmem_gaia",
            "membership_prob_gaia",
            "membership_prob",
            "pmem",
        ]
        for c in cands:
            if c not in self.df.columns:
                continue
            v = pd.to_numeric(self.df[c], errors="coerce").to_numpy(float)
            if np.isfinite(v).sum() >= 10:
                return c, v
        return None, None

    def _merge_columns_from_gaia_derived(self, needed_cols):
        self.df = _merge_gaia_columns_from_catalog(self.df, self.result_dir, needed_cols)

    def _ensure_membership_columns_from_master(self):
        needed = [
            "source_id", "gaia_source_id",
            "pmra", "pmdec", "parallax",
            "pmra_error", "pmdec_error", "parallax_error",
            "ruwe", "visibility_periods_used",
            "gaia_pmem", "pmem_gaia", "membership_prob_gaia",
        ]
        missing = [c for c in needed if c not in self.df.columns]
        if not missing:
            return self._merge_columns_from_gaia_derived(needed)

        try:
            master, _, _ = _load_master_table(self.result_dir)
        except Exception:
            return
        if master.empty:
            return

        key = None
        if "ID" in self.df.columns and "ID" in master.columns:
            key = "ID"
        elif "source_id" in self.df.columns and "source_id" in master.columns:
            key = "source_id"
            self.df["source_id"] = parse_int64_series(self.df["source_id"]).astype("Int64")
            master["source_id"] = parse_int64_series(master["source_id"]).astype("Int64")
        elif "gaia_source_id" in self.df.columns and "gaia_source_id" in master.columns:
            key = "gaia_source_id"
            self.df["gaia_source_id"] = parse_int64_series(self.df["gaia_source_id"]).astype("Int64")
            master["gaia_source_id"] = parse_int64_series(master["gaia_source_id"]).astype("Int64")
        if key is None:
            return

        add_cols = [c for c in missing if c in master.columns and c != key]
        if add_cols:
            use = master[[key] + add_cols].copy()
            use = use.drop_duplicates(subset=[key], keep="first")
            try:
                self.df = self.df.merge(use, on=key, how="left")
            except Exception:
                return

        self._merge_columns_from_gaia_derived(needed)

    @staticmethod
    def _logpdf_gauss(x: np.ndarray, mu: np.ndarray, cov: np.ndarray) -> np.ndarray:
        d = x.shape[1]
        cov_r = np.asarray(cov, float) + np.eye(d) * 1e-6
        sign, logdet = np.linalg.slogdet(cov_r)
        if sign <= 0:
            return np.full(x.shape[0], -np.inf, dtype=float)
        try:
            inv = np.linalg.inv(cov_r)
        except Exception:
            return np.full(x.shape[0], -np.inf, dtype=float)
        diff = x - mu[None, :]
        q = np.einsum("ni,ij,nj->n", diff, inv, diff)
        return -0.5 * (d * np.log(2.0 * np.pi) + logdet + q)

    def _fit_two_component_gmm(self, x_fit: np.ndarray):
        n, d = x_fit.shape
        if n < max(30, d * 8):
            return None

        center = np.nanmedian(x_fit, axis=0)
        mad = np.nanmedian(np.abs(x_fit - center), axis=0)
        mad = np.where(np.isfinite(mad) & (mad > 1e-6), mad, 1.0)
        z = (x_fit - center[None, :]) / mad[None, :]
        d2 = np.sum(z * z, axis=1)
        q40 = float(np.nanquantile(d2, 0.40))
        m0 = d2 <= q40
        if m0.sum() < max(12, d * 3) or m0.sum() > (n - max(12, d * 3)):
            order = np.argsort(d2)
            m0 = np.zeros(n, dtype=bool)
            m0[order[: max(n // 2, 1)]] = True
        m1 = ~m0

        def _cov(arr: np.ndarray) -> np.ndarray:
            if arr.shape[0] < (d + 1):
                return np.eye(d, dtype=float)
            c = np.cov(arr, rowvar=False)
            if np.ndim(c) == 0:
                c = np.eye(d, dtype=float) * float(c)
            return np.asarray(c, float) + np.eye(d, dtype=float) * 1e-4

        pi = np.array([max(m0.mean(), 1e-3), max(m1.mean(), 1e-3)], dtype=float)
        pi /= pi.sum()
        mu = np.vstack([
            np.nanmean(x_fit[m0], axis=0),
            np.nanmean(x_fit[m1], axis=0),
        ])
        cov = np.stack([_cov(x_fit[m0]), _cov(x_fit[m1])], axis=0)
        last_ll = -np.inf

        for _ in range(80):
            lp0 = np.log(max(pi[0], 1e-9)) + self._logpdf_gauss(x_fit, mu[0], cov[0])
            lp1 = np.log(max(pi[1], 1e-9)) + self._logpdf_gauss(x_fit, mu[1], cov[1])
            m = np.maximum(lp0, lp1)
            e0 = np.exp(lp0 - m)
            e1 = np.exp(lp1 - m)
            den = e0 + e1 + 1e-12
            r0 = e0 / den
            r1 = e1 / den
            nk = np.array([r0.sum(), r1.sum()], dtype=float)
            if np.any(nk < (d + 2)):
                break
            pi = nk / float(n)
            mu[0] = (r0[:, None] * x_fit).sum(axis=0) / nk[0]
            mu[1] = (r1[:, None] * x_fit).sum(axis=0) / nk[1]
            for k, rk in enumerate((r0, r1)):
                diff = x_fit - mu[k][None, :]
                cov[k] = (diff.T * rk).dot(diff) / max(nk[k], 1.0)
                cov[k] += np.eye(d, dtype=float) * 1e-4
            ll = float(np.sum(m + np.log(den)))
            if np.isfinite(last_ll):
                if abs(ll - last_ll) < 1e-4 * max(1.0, abs(last_ll)):
                    break
            last_ll = ll

        det0 = abs(float(np.linalg.det(cov[0])))
        det1 = abs(float(np.linalg.det(cov[1])))
        cluster_idx = 0 if det0 <= det1 else 1
        return {
            "pi": pi,
            "mu": mu,
            "cov": cov,
            "cluster_idx": int(cluster_idx),
        }

    def _compute_gaia_membership_prob(self):
        req = ("pmra", "pmdec", "parallax")
        if not all(c in self.df.columns for c in req):
            return None, "pm/parallax columns missing"

        pmra = pd.to_numeric(self.df["pmra"], errors="coerce").to_numpy(float)
        pmdec = pd.to_numeric(self.df["pmdec"], errors="coerce").to_numpy(float)
        plx = pd.to_numeric(self.df["parallax"], errors="coerce").to_numpy(float)
        finite = np.isfinite(pmra) & np.isfinite(pmdec) & np.isfinite(plx)
        if int(finite.sum()) < 30:
            return None, "too few stars with finite pm/parallax"

        fit_mask = finite.copy()
        if "ruwe" in self.df.columns:
            ruwe = pd.to_numeric(self.df["ruwe"], errors="coerce").to_numpy(float)
            fit_mask &= (~np.isfinite(ruwe)) | (ruwe <= 2.0)
        if "visibility_periods_used" in self.df.columns:
            vpu = pd.to_numeric(self.df["visibility_periods_used"], errors="coerce").to_numpy(float)
            fit_mask &= (~np.isfinite(vpu)) | (vpu >= 8.0)
        if int(fit_mask.sum()) < 25:
            fit_mask = finite.copy()

        x_fit = np.column_stack([pmra[fit_mask], pmdec[fit_mask], plx[fit_mask]])
        model = self._fit_two_component_gmm(x_fit)
        if model is None:
            return None, "GMM fit failed"

        x_all = np.column_stack([pmra[finite], pmdec[finite], plx[finite]])
        pi = np.asarray(model["pi"], float)
        mu = np.asarray(model["mu"], float)
        cov = np.asarray(model["cov"], float)
        k_cluster = int(model["cluster_idx"])

        lp0 = np.log(max(pi[0], 1e-9)) + self._logpdf_gauss(x_all, mu[0], cov[0])
        lp1 = np.log(max(pi[1], 1e-9)) + self._logpdf_gauss(x_all, mu[1], cov[1])
        m = np.maximum(lp0, lp1)
        e0 = np.exp(lp0 - m)
        e1 = np.exp(lp1 - m)
        den = e0 + e1 + 1e-12
        r0 = e0 / den
        r1 = e1 / den
        p_cluster = r0 if k_cluster == 0 else r1

        pmem = np.full(len(self.df), np.nan, dtype=float)
        pmem[finite] = p_cluster
        note = f"gaia_gmm_3d(valid={int(finite.sum())}, fit={int(fit_mask.sum())})"
        return pmem, note

    def _ensure_membership_prob(self):
        if self._membership_ready:
            return self._membership_prob
        self._membership_ready = True

        c_exist, v_exist = self._pick_existing_membership_col()
        if c_exist is not None:
            self._membership_prob = np.clip(np.asarray(v_exist, float), 0.0, 1.0)
            self._membership_source = c_exist
            self._membership_note = f"existing:{c_exist}"
            return self._membership_prob

        self._ensure_membership_columns_from_master()
        c_exist, v_exist = self._pick_existing_membership_col()
        if c_exist is not None:
            self._membership_prob = np.clip(np.asarray(v_exist, float), 0.0, 1.0)
            self._membership_source = c_exist
            self._membership_note = f"existing:{c_exist}"
            return self._membership_prob

        p_auto, note = self._compute_gaia_membership_prob()
        if p_auto is not None:
            self._membership_prob = np.clip(np.asarray(p_auto, float), 0.0, 1.0)
            self.df["gaia_pmem"] = self._membership_prob
            self._membership_source = "gaia_gmm_3d"
            self._membership_note = note
            return self._membership_prob

        self._membership_prob = None
        self._membership_source = "none"
        self._membership_note = note
        return None

    def _current_membership_mask(self):
        mode = self._membership_mode_key()
        if mode == "off":
            return None, mode, np.nan
        prob = self._ensure_membership_prob()
        if prob is None:
            return None, mode, self._membership_threshold(mode)
        thr = self._membership_threshold(mode)
        mask = np.isfinite(prob) & (prob >= thr)
        return mask, mode, thr

    def _save_membership_csv(self):
        prob = self._ensure_membership_prob()
        if prob is None:
            self.info_text.setPlainText(
                "Membership not available.\n"
                f"Reason: {self._membership_note}"
            )
            return
        out = self.result_dir / "cmd_with_gaia_membership.csv"
        save_cols = ["ID", "source_id", "gaia_source_id", "pmra", "pmdec", "parallax"]
        keep = [c for c in save_cols if c in self.df.columns]
        df_out = self.df[keep].copy() if keep else pd.DataFrame(index=self.df.index)
        df_out["gaia_pmem"] = prob
        df_out.to_csv(out, index=False, na_rep="NaN")
        self.info_text.setPlainText(f"Saved: {out}")

    def _build_figure(self):
        self.figure.clear()
        view_name = self.available_views[self.view_mode]

        # Single view modes
        # Padded right margin + cax pulled inward keeps the colorbar
        # label "Teff (K) + OBAFGKM-like color" fully on-canvas even
        # with long tick labels like "35000 K (O)".  top=0.88 leaves
        # breathing room for the plot title.
        if view_name == "inst":
            self.ax_inst = self.figure.add_subplot(1, 1, 1)
            self.ax_std = None
            self.ax_gaia = None
            self.figure.subplots_adjust(left=0.10, right=0.80, bottom=0.18, top=0.88)
            self.cax = self.figure.add_axes([0.84, 0.18, 0.018, 0.70])
        elif view_name == "std":
            self.ax_inst = None
            self.ax_std = self.figure.add_subplot(1, 1, 1)
            self.ax_gaia = None
            self.figure.subplots_adjust(left=0.10, right=0.80, bottom=0.18, top=0.88)
            self.cax = self.figure.add_axes([0.84, 0.18, 0.018, 0.70])
        elif view_name == "gaia":
            self.ax_inst = None
            self.ax_std = None
            self.ax_gaia = self.figure.add_subplot(1, 1, 1)
            self.figure.subplots_adjust(left=0.10, right=0.80, bottom=0.18, top=0.88)
            self.cax = self.figure.add_axes([0.84, 0.18, 0.018, 0.70])
        elif view_name == "all":
            # Show all available CMDs
            if self.has_std and self.gaia_mode is not None:
                self.ax_inst = self.figure.add_subplot(1, 3, 1)
                self.ax_std = self.figure.add_subplot(1, 3, 2)
                self.ax_gaia = self.figure.add_subplot(1, 3, 3)
                self.figure.subplots_adjust(left=0.055, right=0.85, bottom=0.16, top=0.85, wspace=0.30)
                self.cax = self.figure.add_axes([0.88, 0.16, 0.015, 0.68])
            elif self.has_std:
                self.ax_inst = self.figure.add_subplot(1, 2, 1)
                self.ax_std = self.figure.add_subplot(1, 2, 2)
                self.ax_gaia = None
                self.figure.subplots_adjust(left=0.07, right=0.85, bottom=0.16, top=0.85, wspace=0.30)
                self.cax = self.figure.add_axes([0.88, 0.16, 0.015, 0.68])
            elif self.gaia_mode is not None:
                self.ax_inst = self.figure.add_subplot(1, 2, 1)
                self.ax_std = None
                self.ax_gaia = self.figure.add_subplot(1, 2, 2)
                self.figure.subplots_adjust(left=0.07, right=0.85, bottom=0.16, top=0.85, wspace=0.30)
                self.cax = self.figure.add_axes([0.88, 0.16, 0.015, 0.68])
            else:
                self.ax_inst = self.figure.add_subplot(1, 1, 1)
                self.ax_std = None
                self.ax_gaia = None
                self.figure.subplots_adjust(left=0.13, right=0.80, bottom=0.16, top=0.85)
                self.cax = self.figure.add_axes([0.84, 0.18, 0.018, 0.66])
        else:
            self.ax_inst = self.figure.add_subplot(1, 1, 1)
            self.ax_std = None
            self.ax_gaia = None
            self.figure.subplots_adjust(left=0.13, right=0.80, bottom=0.14, top=0.85)
            self.cax = self.figure.add_axes([0.84, 0.18, 0.018, 0.66])

        self.figure.patch.set_facecolor("black")
        for ax in [self.ax_inst, self.ax_std, self.ax_gaia]:
            if ax is None:
                continue
            ax.set_facecolor("black")
            for sp in ax.spines.values():
                sp.set_color("white")
            ax.tick_params(colors="white")
            ax.xaxis.label.set_color("white")
            ax.yaxis.label.set_color("white")
            ax.title.set_color("white")

        self.cax.set_facecolor("black")
        sm = mpl.cm.ScalarMappable(norm=self.ob_norm, cmap=self.ob_cmap)
        sm.set_array([])
        cbar = self.figure.colorbar(sm, cax=self.cax)
        cbar.set_label("Teff (K) + OBAFGKM-like color", fontsize=9, color="white", labelpad=14)
        ticks = [35000, 20000, 10000, 7500, 6000, 4500, 3200]
        labels = ["35000 K (O)", "20000 K (B)", "10000 K (A)", " 7500 K (F)", " 6000 K (G)", " 4500 K (K)", " 3200 K (M)"]
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(labels)
        cbar.ax.tick_params(colors="white")
        for sp in cbar.ax.spines.values():
            sp.set_color("white")

    def _safe_float(self, series):
        return pd.to_numeric(series, errors="coerce").to_numpy(float)

    def _teff_from_color_index(self, color_x: np.ndarray, mode: str):
        if color_x.size == 0:
            return np.full_like(color_x, np.nan, dtype=float)
        parts = mode.split("-", 1)
        if len(parts) == 2:
            return _teff_from_color(color_x, parts[0], parts[1], self.teff_vmin, self.teff_vmax)
        return np.full_like(color_x, np.nan, dtype=float)

    def _get_y_mode(self, yval):
        if yval in self.y_scalar_opts:
            return ("scalar", yval)
        if isinstance(yval, str) and "-" in yval:
            a, b = yval.split("-", 1)
            if (a, b) in self.y_color_pairs:
                return ("color", (a, b))
        return (None, None)

    def _compute_arrays_and_mask(self, system: str, x_pair, y_choice, snr_cut: float, membership_mask=None):
        a, b = x_pair
        col_ax = f"mag_{system}_{a}"
        col_bx = f"mag_{system}_{b}"
        if (col_ax not in self.df.columns) or (col_bx not in self.df.columns):
            return np.array([]), np.array([]), np.zeros(len(self.df), bool), np.array([])

        Ax = self._safe_float(self.df[col_ax])
        Bx = self._safe_float(self.df[col_bx])
        xcolor = Ax - Bx
        x = xcolor

        y_mode, y_param = self._get_y_mode(y_choice)
        involved = set([a, b])

        if y_mode == "scalar":
            by = y_param
            col_y = f"mag_{system}_{by}"
            if col_y not in self.df.columns:
                return np.array([]), np.array([]), np.zeros(len(self.df), bool), np.array([])
            y = self._safe_float(self.df[col_y])
            # Display-only manual ZP: optional extra shift of the Instrumental
            # magnitude axis. mag_inst already carries the IRAF Z=25 convention
            # (baked in at Step 7), so this defaults to 0; it only adds a manual
            # nudge on top. Only the magnitude (scalar) axis of the Instrumental
            # view is shifted; colors and the Std view are untouched (a constant
            # ZP cancels in any a-b color).
            if system == "inst" and getattr(self, "manual_zp_check", None) is not None \
                    and self.manual_zp_check.isChecked():
                y = y + float(self.manual_zp_spin.value())
            involved.add(by)
        elif y_mode == "color":
            ya, yb = y_param
            col_ya = f"mag_{system}_{ya}"
            col_yb = f"mag_{system}_{yb}"
            if (col_ya not in self.df.columns) or (col_yb not in self.df.columns):
                return np.array([]), np.array([]), np.zeros(len(self.df), bool), np.array([])
            Ya = self._safe_float(self.df[col_ya])
            Yb = self._safe_float(self.df[col_yb])
            y = Ya - Yb
            involved.update([ya, yb])
        else:
            return np.array([]), np.array([]), np.zeros(len(self.df), bool), np.array([])

        mask = np.isfinite(x) & np.isfinite(y)

        if snr_cut > 0 and self.has_snr:
            for band in involved:
                sc = f"snr_{band}"
                if sc in self.df.columns:
                    sv = self._safe_float(self.df[sc])
                    # Only exclude stars with known (finite) SNR below threshold;
                    # NaN SNR means unmeasured → keep (do not reject unknowns)
                    mask &= ~(np.isfinite(sv) & (sv < snr_cut))

        if membership_mask is not None and len(membership_mask) == len(mask):
            mask &= np.asarray(membership_mask, bool)

        return x[mask], y[mask], mask, xcolor[mask]

    def _mag_error_array(self, system: str, band: str) -> np.ndarray:
        candidates = []
        if system == "std":
            candidates.append(f"mag_std_err_{band}")
        candidates.append(f"mag_inst_err_{band}")
        for col in candidates:
            if col in self.df.columns:
                arr = self._safe_float(self.df[col])
                return np.where(np.isfinite(arr) & (arr >= 0), arr, np.nan)
        return np.full(len(self.df), np.nan, dtype=float)

    def _quadrature_error(self, *arrays: np.ndarray) -> np.ndarray:
        if not arrays:
            return np.full(len(self.df), np.nan, dtype=float)
        stack = np.vstack([np.asarray(a, dtype=float) for a in arrays])
        finite = np.isfinite(stack).all(axis=0)
        out = np.full(stack.shape[1], np.nan, dtype=float)
        out[finite] = np.sqrt(np.sum(stack[:, finite] ** 2, axis=0))
        return out

    def _compute_cmd_error_arrays(self, system: str, x_pair, y_choice, mask: np.ndarray):
        if mask is None or len(mask) != len(self.df):
            return np.array([]), np.array([])
        a, b = x_pair
        xerr_full = self._quadrature_error(
            self._mag_error_array(system, a),
            self._mag_error_array(system, b),
        )

        y_mode, y_param = self._get_y_mode(y_choice)
        if y_mode == "scalar":
            yerr_full = self._mag_error_array(system, y_param)
        elif y_mode == "color":
            ya, yb = y_param
            yerr_full = self._quadrature_error(
                self._mag_error_array(system, ya),
                self._mag_error_array(system, yb),
            )
        else:
            yerr_full = np.full(len(self.df), np.nan, dtype=float)

        return xerr_full[mask], yerr_full[mask]

    def _plot_cmd_errorbars(self, ax, x, y, mask: np.ndarray, system: str, x_pair, y_choice):
        show_x = getattr(self, "xerr_check", None) is not None and self.xerr_check.isChecked()
        show_y = getattr(self, "yerr_check", None) is not None and self.yerr_check.isChecked()
        if not (show_x or show_y) or len(x) == 0:
            return

        xerr, yerr = self._compute_cmd_error_arrays(system, x_pair, y_choice, mask)
        if len(xerr) != len(x) or len(yerr) != len(y):
            return

        finite = np.isfinite(x) & np.isfinite(y)
        if show_x:
            finite &= np.isfinite(xerr) & (xerr > 0)
        if show_y:
            finite &= np.isfinite(yerr) & (yerr > 0)
        idx = np.flatnonzero(finite)
        if idx.size == 0:
            return

        max_bars = 400
        if idx.size > max_bars:
            idx = idx[np.linspace(0, idx.size - 1, max_bars).astype(int)]

        ax.errorbar(
            np.asarray(x)[idx],
            np.asarray(y)[idx],
            xerr=np.asarray(xerr)[idx] if show_x else None,
            yerr=np.asarray(yerr)[idx] if show_y else None,
            fmt="none",
            ecolor="#DDE7F0",
            elinewidth=0.55,
            alpha=0.35,
            capsize=0,
            zorder=1,
        )

    def _compute_gaia_arrays_and_mask(self, snr_cut: float, membership_mask=None):
        if self.gaia_mode is None:
            return np.array([]), np.array([]), np.zeros(len(self.df), bool), np.array([])

        if self.gaia_mode == "inst":
            G = self._safe_float(self.df["gaia_G_inst"])
            C = self._safe_float(self.df["gaia_BP_RP_inst"])
        else:
            G = self._safe_float(self.df["gaia_G_syn"])
            C = self._safe_float(self.df["gaia_BP_RP_syn"])

        mask = np.isfinite(G) & np.isfinite(C)

        if snr_cut > 0 and self.has_snr:
            for band in ("g", "r", "i"):
                sc = f"snr_{band}"
                if sc in self.df.columns:
                    sv = self._safe_float(self.df[sc])
                    mask &= ~(np.isfinite(sv) & (sv < snr_cut))

        if membership_mask is not None and len(membership_mask) == len(mask):
            mask &= np.asarray(membership_mask, bool)

        return C[mask], G[mask], mask, C[mask]

    def _style_axis(self, ax):
        ax.set_facecolor("black")
        for sp in ax.spines.values():
            sp.set_color("white")
        ax.tick_params(colors="white")
        ax.xaxis.label.set_color("white")
        ax.yaxis.label.set_color("white")
        ax.title.set_color("white")

    def _apply_y_orientation(self):
        axes = [self.ax_inst, self.ax_std, self.ax_gaia]
        for ax in axes:
            if ax is None:
                continue
            ymin, ymax = ax.get_ylim()
            if self.invert_y.isChecked():
                if ymin < ymax:
                    ax.set_ylim(ymax, ymin)
            else:
                if ymin > ymax:
                    ax.set_ylim(ymax, ymin)

    def _redraw(self):
        ax_primary = self.ax_inst or self.ax_std or self.ax_gaia
        self._plot_cache = {}
        if not self.x_pairs:
            if ax_primary is not None:
                ax_primary.clear()
                ax_primary.set_title("No available X color index", fontsize=10, color="white")
            self.canvas.draw_idle()
            return

        x_text = self.x_combo.currentText()
        if x_text not in [f"{a}-{b}" for (a, b) in self.x_pairs]:
            self.canvas.draw_idle()
            return
        a, b = x_text.split("-", 1)
        x_pair = (a, b)

        yval = self.y_combo.currentText()
        snr_cut = float(self.snr_spin.value())
        member_mask, member_mode, member_thr = self._current_membership_mask()
        member_active = (member_mode != "off") and (member_mask is not None)
        member_compare = bool(self.member_compare.isChecked()) and member_active

        plx_mask = self._parallax_mask()
        plx_active = plx_mask is not None

        roi_mask = self._roi_mask()
        roi_active = roi_mask is not None

        # background (grey dots): spatial/parallax pre-filter (everything inside ROI or parallax range)
        if plx_active and roi_active:
            bg_mask = plx_mask & roi_mask
        elif plx_active:
            bg_mask = plx_mask
        elif roi_active:
            bg_mask = roi_mask
        else:
            bg_mask = None

        # Gaia quality (RUWE / C*) filters the *drawn* sample, background
        # included: the point of enabling them is to see the CMD without the
        # contaminated stars, so leaving them in the grey layer would defeat it.
        qual_mask = self._quality_mask()
        if qual_mask is not None:
            bg_mask = qual_mask if bg_mask is None else (bg_mask & qual_mask)

        # foreground mask: member & spatial filters combined
        spatial_mask = bg_mask  # reuse combined spatial pre-filter
        spatial_active = spatial_mask is not None
        if member_active and spatial_active:
            effective_mask = np.asarray(member_mask, bool) & spatial_mask
        elif member_active:
            effective_mask = member_mask
        elif spatial_active:
            effective_mask = spatial_mask
        else:
            effective_mask = None

        if self.ax_inst is not None:
            self.ax_inst.clear()
        if self.ax_std is not None:
            self.ax_std.clear()
        if self.ax_gaia is not None:
            self.ax_gaia.clear()

        if self.ax_inst is not None:
            x_i_all, y_i_all, mask_i_all, xcol_i_all = self._compute_arrays_and_mask("inst", x_pair, yval, snr_cut, bg_mask)
            if member_active or plx_active or roi_active:
                x_i, y_i, mask_i, xcol_i = self._compute_arrays_and_mask("inst", x_pair, yval, snr_cut, effective_mask)
            else:
                x_i, y_i, mask_i, xcol_i = x_i_all, y_i_all, mask_i_all, xcol_i_all
            teff_i = self._teff_from_color_index(xcol_i, f"{a}-{b}")
        else:
            x_i_all, y_i_all = np.array([]), np.array([])
            x_i, y_i, mask_i, teff_i = np.array([]), np.array([]), np.zeros(len(self.df), bool), np.array([])

        if self.has_std and self.ax_std is not None:
            x_s_all, y_s_all, mask_s_all, xcol_s_all = self._compute_arrays_and_mask("std", x_pair, yval, snr_cut, bg_mask)
            if member_active or plx_active or roi_active:
                x_s, y_s, mask_s, xcol_s = self._compute_arrays_and_mask("std", x_pair, yval, snr_cut, effective_mask)
            else:
                x_s, y_s, mask_s, xcol_s = x_s_all, y_s_all, mask_s_all, xcol_s_all
            teff_s = self._teff_from_color_index(xcol_s, f"{a}-{b}")
        else:
            x_s_all, y_s_all = np.array([]), np.array([])
            x_s, y_s, mask_s, teff_s = np.array([]), np.array([]), np.zeros(len(self.df), bool), np.array([])

        if self.gaia_mode is not None and self.ax_gaia is not None:
            x_g_all, y_g_all, mask_g_all, xcol_g_all = self._compute_gaia_arrays_and_mask(snr_cut, bg_mask)
            if member_active or plx_active or roi_active:
                x_g, y_g, mask_g, xcol_g = self._compute_gaia_arrays_and_mask(snr_cut, effective_mask)
            else:
                x_g, y_g, mask_g, xcol_g = x_g_all, y_g_all, mask_g_all, xcol_g_all
            teff_g = self._teff_from_color_index(xcol_g, "BP-RP")
        else:
            x_g_all, y_g_all = np.array([]), np.array([])
            x_g, y_g, mask_g, teff_g = np.array([]), np.array([]), np.zeros(len(self.df), bool), np.array([])

        color_title = f"{a}-{b}"

        def _count_text(n, n_all, active):
            return f"N={n}/{n_all}" if active else f"N={n}"

        if self.ax_inst is not None:
            self._style_axis(self.ax_inst)
            if member_compare and len(x_i_all) > 0:
                self.ax_inst.scatter(x_i_all, y_i_all, s=10, alpha=0.22, linewidths=0, rasterized=True, c="#9E9E9E")
            if len(x_i) > 0:
                self._plot_cmd_errorbars(self.ax_inst, x_i, y_i, mask_i, "inst", x_pair, yval)
                self.ax_inst.scatter(x_i, y_i, s=12, alpha=0.92, linewidths=0, rasterized=True, c=teff_i, cmap=self.ob_cmap, norm=self.ob_norm)
                self.ax_inst.set_title(
                    format_cmd_title(
                        self.params,
                        yval,
                        color_title,
                        system_label="Instrumental",
                        count_text=_count_text(len(x_i), len(x_i_all), member_active or plx_active or roi_active),
                        result_dir=self.result_dir,
                    ),
                    fontsize=11,
                    color="white",
                )
                self._plot_cache[self.ax_inst] = {
                    "system": "inst",
                    "x": x_i,
                    "y": y_i,
                    "df_index": np.where(mask_i)[0],
                }
            else:
                self.ax_inst.set_title(
                    format_cmd_title(
                        self.params,
                        yval,
                        color_title,
                        system_label="Instrumental",
                        count_text=_count_text(0, len(x_i_all), member_active or plx_active or roi_active),
                        result_dir=self.result_dir,
                    ),
                    fontsize=11,
                    color="white",
                )

        if self.ax_std is not None:
            std_label = photometric_system_label(a, b, yval)
            self._style_axis(self.ax_std)
            if member_compare and len(x_s_all) > 0:
                self.ax_std.scatter(x_s_all, y_s_all, s=10, alpha=0.22, linewidths=0, rasterized=True, c="#9E9E9E")
            if len(x_s) > 0:
                self._plot_cmd_errorbars(self.ax_std, x_s, y_s, mask_s, "std", x_pair, yval)
                self.ax_std.scatter(x_s, y_s, s=12, alpha=0.92, linewidths=0, rasterized=True, c=teff_s, cmap=self.ob_cmap, norm=self.ob_norm)
                self.ax_std.set_title(
                    format_cmd_title(
                        self.params,
                        yval,
                        color_title,
                        system_label=std_label,
                        count_text=_count_text(len(x_s), len(x_s_all), member_active or plx_active or roi_active),
                        result_dir=self.result_dir,
                    ),
                    fontsize=11,
                    color="white",
                )
                self._plot_cache[self.ax_std] = {
                    "system": "std",
                    "x": x_s,
                    "y": y_s,
                    "df_index": np.where(mask_s)[0],
                }
            else:
                self.ax_std.set_title(
                    format_cmd_title(
                        self.params,
                        yval,
                        color_title,
                        system_label=std_label,
                        count_text=_count_text(0, len(x_s_all), member_active or plx_active or roi_active),
                        result_dir=self.result_dir,
                    ),
                    fontsize=11,
                    color="white",
                )

        if self.ax_gaia is not None:
            self._style_axis(self.ax_gaia)
            if member_compare and len(x_g_all) > 0:
                self.ax_gaia.scatter(x_g_all, y_g_all, s=10, alpha=0.22, linewidths=0, rasterized=True, c="#9E9E9E")
            if len(x_g) > 0:
                self.ax_gaia.scatter(x_g, y_g, s=12, alpha=0.92, linewidths=0, rasterized=True, c=teff_g, cmap=self.ob_cmap, norm=self.ob_norm)
                gaia_label = "Gaia instrumental" if self.gaia_mode == "inst" else "Gaia synthetic"
                self.ax_gaia.set_title(
                    format_cmd_title(
                        self.params,
                        "G",
                        "BP-RP",
                        system_label=gaia_label,
                        count_text=_count_text(len(x_g), len(x_g_all), member_active or plx_active or roi_active),
                        result_dir=self.result_dir,
                    ),
                    fontsize=11,
                    color="white",
                )
                self._plot_cache[self.ax_gaia] = {
                    "system": "gaia",
                    "x": x_g,
                    "y": y_g,
                    "df_index": np.where(mask_g)[0],
                }
            else:
                gaia_label = "Gaia instrumental" if self.gaia_mode == "inst" else "Gaia synthetic"
                self.ax_gaia.set_title(
                    format_cmd_title(
                        self.params,
                        "G",
                        "BP-RP",
                        system_label=gaia_label,
                        count_text=_count_text(0, len(x_g_all), member_active or plx_active or roi_active),
                        result_dir=self.result_dir,
                    ),
                    fontsize=11,
                    color="white",
                )

        x_label = f"{a}-{b} (mag)"
        y_label = f"{yval} (mag)" if yval else ""
        if self.ax_inst is not None:
            self.ax_inst.set_xlabel(x_label)
            self.ax_inst.set_ylabel(y_label)
        if self.ax_std is not None:
            self.ax_std.set_xlabel(x_label)
            self.ax_std.set_ylabel(y_label)
        if self.ax_gaia is not None:
            self.ax_gaia.set_xlabel("BP-RP (mag)")
            self.ax_gaia.set_ylabel("G (mag)")

        self._apply_y_orientation()

        def _rng(arr):
            arr = np.asarray(arr, float)
            arr = arr[np.isfinite(arr)]
            if arr.size == 0:
                return "n/a"
            return f"{arr.min():.0f}-{arr.max():.0f} K"

        lines = [
            f"X={a}-{b}, Y={yval}, SNR>={snr_cut:.0f}",
        ]
        if member_mode == "off":
            lines.append("Membership: Off")
        elif member_mask is None:
            lines.append(f"Membership: {member_mode} (P>={member_thr:.2f}) unavailable | {self._membership_note}")
        else:
            n_valid = int(np.isfinite(self._membership_prob).sum()) if self._membership_prob is not None else 0
            n_mem = int(np.asarray(member_mask, bool).sum())
            lines.append(
                f"Membership: {member_mode} (P>={member_thr:.2f}) | source={self._membership_source} | "
                f"selected={n_mem}/{n_valid} | compare={member_compare}"
            )
        if self.plx_check.isChecked():
            plx = self._parallax_values()
            if plx is None:
                lines.append("Parallax: unavailable")
            else:
                finite = np.isfinite(plx)
                n_sel = int(plx_mask.sum()) if plx_mask is not None else 0
                lines.append(
                    f"Parallax: {self.plx_min_spin.value():.3f}..{self.plx_max_spin.value():.3f} mas | "
                    f"selected={n_sel}/{int(finite.sum())}"
                )
        if self.inst_bands:
            lines.append(f"[Inst] N={len(x_i)}{'/' + str(len(x_i_all)) if member_active else ''} | Teff range: {_rng(teff_i)}")
        if self.has_std:
            lines.append(f"[Cal]  N={len(x_s)}{'/' + str(len(x_s_all)) if member_active else ''} | Teff range: {_rng(teff_s)}")
        if self.gaia_mode is not None:
            lines.append(f"[Gaia:{self.gaia_mode}] N={len(x_g)}{'/' + str(len(x_g_all)) if member_active else ''} | Teff range: {_rng(teff_g)}")
        if not self.has_snr:
            lines.append("(snr_* columns missing: SNR cut disabled)")
        if self.last_pick_info:
            lines.append(self.last_pick_info)
        if self.pick_log:
            lines.append("Pick log (latest 5):")
            lines.extend(self.pick_log[-5:])
        self.info_text.setPlainText("\n".join(lines))

        self.canvas.draw_idle()

    def _fmt_val(self, v, nd=3):
        try:
            if v is None or (isinstance(v, float) and not np.isfinite(v)):
                return "NaN"
            return f"{float(v):.{nd}f}"
        except Exception:
            return str(v)

    def _on_plot_click(self, event):
        if event.inaxes is None or event.xdata is None or event.ydata is None:
            return
        if not getattr(event, "dblclick", False):
            return
        if event.inaxes not in self._plot_cache:
            return
        cache = self._plot_cache[event.inaxes]
        x = cache["x"]
        y = cache["y"]
        if x.size == 0:
            return
        xy_disp = event.inaxes.transData.transform(np.column_stack([x, y]))
        click = np.array([event.x, event.y])
        d2 = np.sum((xy_disp - click) ** 2, axis=1)
        idx = int(np.argmin(d2))
        if d2[idx] > (12.0 ** 2):
            return
        df_idx = int(cache["df_index"][idx])
        row = self.df.iloc[df_idx]

        parts = [f"Pick[{cache['system']}] ID={row.get('ID', 'n/a')}", f"source_id={row.get('source_id', 'n/a')}"]
        for band in ("g", "r", "i"):
            c_inst = f"mag_inst_{band}"
            if c_inst in row:
                parts.append(f"{band}_inst={self._fmt_val(row.get(c_inst))}")
        for band in ("g", "r", "i"):
            c_cal = f"mag_cal_{band}"
            if c_cal in row:
                parts.append(f"{band}_cal={self._fmt_val(row.get(c_cal))}")
        for band in ("g", "r", "i"):
            c_std = f"mag_std_{band}"
            if c_std in row:
                parts.append(f"{band}_std={self._fmt_val(row.get(c_std))}")
        if "gaia_G" in row:
            parts.append(f"gaia_G={self._fmt_val(row.get('gaia_G'))}")
        if "gaia_BP" in row:
            parts.append(f"gaia_BP={self._fmt_val(row.get('gaia_BP'))}")
        if "gaia_RP" in row:
            parts.append(f"gaia_RP={self._fmt_val(row.get('gaia_RP'))}")
        for c in ("pmra", "pmdec", "parallax"):
            if c in row:
                parts.append(f"{c}={self._fmt_val(row.get(c))}")
        if self._membership_prob is not None and 0 <= df_idx < len(self._membership_prob):
            pm = self._membership_prob[df_idx]
            if np.isfinite(pm):
                parts.append(f"Pmem={float(pm):.3f}")
        msg = " | ".join(parts)
        self.last_pick_info = msg
        self.pick_log.append(msg)
        self._redraw()

    def _save_png(self):
        if not self.x_pairs:
            self.info_text.setPlainText("No available X color index")
            return
        a, b = self.x_combo.currentText().split("-", 1)
        yv = self.y_combo.currentText().replace(" ", "")

        if self.has_std and self.gaia_mode is not None:
            mode = f"inst_std_gaia{self.gaia_mode}"
        elif self.has_std:
            mode = "inst_std"
        elif self.gaia_mode is not None:
            mode = f"inst_gaia{self.gaia_mode}"
        else:
            mode = "inst_only"

        mem_tag = ""
        mkey = self._membership_mode_key()
        if mkey != "off":
            mem_tag = f"_mem{mkey}"
            if self.member_compare.isChecked():
                mem_tag += "_cmp"

        out = self.result_dir / f"cmd_{mode}_{a}-{b}_vs_{yv}_snr{int(self.snr_spin.value())}{mem_tag}_OBcolor_dark.png"
        self.figure.savefig(out, dpi=170, bbox_inches="tight", facecolor=self.figure.get_facecolor(), edgecolor="none")
        self.info_text.setPlainText(f"Saved: {out}")

    def keyPressEvent(self, event):
        super().keyPressEvent(event)

    def _switch_view(self, delta: int):
        """Switch between views: inst, std, gaia, all"""
        self.view_mode = (self.view_mode + delta) % len(self.available_views)
        self._update_view_label()
        self._build_figure()
        self._redraw()


class ZPFitPlotWidget(QWidget):
    """ZP linear fit (delta vs color) + per-frame ZP timeline plots."""

    FILT_COLORS = {"g": "#1976D2", "r": "#D32F2F", "i": "#388E3C",
                   "u": "#7B1FA2", "z": "#E65100", "y": "#00695C",
                   "ha": "#AD1457", "r_spec": "#B71C1C"}

    def __init__(self, result_dir: Path, parent=None):
        super().__init__(parent)
        self.result_dir = Path(result_dir)
        self._cal_df = None
        self._frame_df = None
        self._coeff_df = None
        self._artist_map = {}   # id(artist) -> list of row dicts
        self._setup_ui()

    def _filter_color(self, filt: str) -> str:
        """Return a consistent color for any filter name."""
        if filt in self.FILT_COLORS:
            return self.FILT_COLORS[filt]
        import hashlib
        h = int(hashlib.md5(filt.encode()).hexdigest()[:6], 16)
        # Keep it reasonably saturated (not too dark/light)
        r = 80 + (h >> 16 & 0xFF) % 140
        g = 80 + (h >> 8 & 0xFF) % 140
        b = 80 + (h & 0xFF) % 140
        return f"#{r:02X}{g:02X}{b:02X}"

    def _detected_filters(self) -> list[str]:
        """Collect all filter names present in loaded data."""
        filts: set[str] = set()
        if self._cal_df is not None:
            filts.update(c[len("delta_"):] for c in self._cal_df.columns if c.startswith("delta_"))
        if self._frame_df is not None and "filter" in self._frame_df.columns:
            filts.update(str(f) for f in self._frame_df["filter"].dropna().unique())
        if self._coeff_df is not None and "filter" in self._coeff_df.columns:
            filts.update(str(f) for f in self._coeff_df["filter"].dropna().unique())
        return sorted(filts)

    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)

        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel("Filter:"))
        self._filt_combo = QComboBox()
        self._filt_combo.addItem("All")
        self._filt_combo.currentTextChanged.connect(self._redraw)
        ctrl.addWidget(self._filt_combo)
        ctrl.addSpacing(12)
        ctrl.addWidget(QLabel("Date:"))
        self._date_combo = QComboBox()
        self._date_combo.addItem("All")
        self._date_combo.currentTextChanged.connect(self._redraw)
        ctrl.addWidget(self._date_combo)
        ctrl.addSpacing(12)
        btn_reload = QPushButton("Reload")
        btn_reload.clicked.connect(lambda: self.reload())
        ctrl.addWidget(btn_reload)
        btn_save = QPushButton("Save PNG")
        btn_save.clicked.connect(self._save_png)
        ctrl.addWidget(btn_save)
        ctrl.addStretch()
        layout.addLayout(ctrl)

        self._fig = Figure(figsize=(12, 8))
        gs = self._fig.add_gridspec(2, 2, hspace=0.35, wspace=0.30)
        self._ax_fit = self._fig.add_subplot(gs[0, 0])
        self._ax_hist = self._fig.add_subplot(gs[0, 1])
        self._ax_frame = self._fig.add_subplot(gs[1, :])
        self._canvas = FigureCanvas(self._fig)
        self._canvas.mpl_connect("pick_event", self._on_pick)
        self._toolbar = NavigationToolbar(self._canvas, self)
        layout.addWidget(self._toolbar)
        layout.addWidget(self._canvas, 1)

        self._info_label = QLabel("Click a data point to see star info.")
        self._info_label.setStyleSheet(mono_note_style())
        layout.addWidget(self._info_label)

    def reload(self, result_dir: Path = None):
        if result_dir is not None:
            self.result_dir = Path(result_dir)
        out_dir = step10_zp_dir(self.result_dir)
        self._cal_df = None
        self._frame_df = None
        self._coeff_df = None

        cal_path = out_dir / "gaia_sdss_calibrator_by_ID.csv"
        if cal_path.exists() and cal_path.stat().st_size > 0:
            try:
                self._cal_df = pd.read_csv(cal_path)
            except Exception:
                self._cal_df = None

        frame_path = out_dir / "frame_zeropoint.csv"
        if frame_path.exists() and frame_path.stat().st_size > 0:
            try:
                self._frame_df = pd.read_csv(frame_path)
            except Exception:
                self._frame_df = None

        coeff_path = out_dir / "zp_fit_coefficients.csv"
        if coeff_path.exists() and coeff_path.stat().st_size > 0:
            try:
                self._coeff_df = pd.read_csv(coeff_path)
            except Exception:
                self._coeff_df = None

        self._update_filter_combo()
        self._update_date_combo()
        self._redraw()

    def _update_filter_combo(self) -> None:
        """Repopulate filter combo from currently loaded data."""
        prev = self._filt_combo.currentText()
        self._filt_combo.blockSignals(True)
        self._filt_combo.clear()
        self._filt_combo.addItem("All")
        for f in self._detected_filters():
            self._filt_combo.addItem(f)
        idx = self._filt_combo.findText(prev)
        self._filt_combo.setCurrentIndex(max(idx, 0))
        self._filt_combo.blockSignals(False)

    def _update_date_combo(self):
        import re
        prev = self._date_combo.currentText()
        self._date_combo.blockSignals(True)
        self._date_combo.clear()
        self._date_combo.addItem("All")
        if self._frame_df is not None and "file" in self._frame_df.columns:
            def _extract_date(fname):
                m = re.search(r"\d{4}-\d{2}-\d{2}", str(fname))
                return m.group(0) if m else None
            dates = sorted(set(d for d in self._frame_df["file"].apply(_extract_date) if d))
            for d in dates:
                self._date_combo.addItem(d)
        idx = self._date_combo.findText(prev)
        self._date_combo.setCurrentIndex(max(idx, 0))
        self._date_combo.blockSignals(False)

    def _redraw(self):
        self._ax_fit.cla()
        self._ax_hist.cla()
        self._ax_frame.cla()
        self._artist_map.clear()
        self._draw_fit_plot()
        self._draw_zp_hist()
        self._draw_frame_zp()
        self._fig.tight_layout(pad=2.0)
        self._canvas.draw_idle()

    def _draw_fit_plot(self):
        ax = self._ax_fit
        filt_sel = self._filt_combo.currentText()

        if self._cal_df is None:
            ax.set_title("ZP Linear Fit — run ZP calibration first")
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes,
                    fontsize=14, color="gray")
            return

        cal = self._cal_df
        filts = self._detected_filters() if filt_sel == "All" else [filt_sel]

        color_labels = set()
        for filt in filts:
            fc = self._filter_color(filt)

            # Determine color index column from coefficients if available, else guess
            color_col = None
            if self._coeff_df is not None:
                row_c = self._coeff_df[self._coeff_df["filter"] == filt]
                if len(row_c) and "color_col" in row_c.columns:
                    cc = str(row_c["color_col"].iloc[0])
                    if cc != "none":
                        # coeff_df stores bare name ("B_V"); cal column has "color_" prefix
                        color_col = cc if cc.startswith("color_") else f"color_{cc}"
            if color_col is None or color_col not in cal.columns:
                # Fallback: scan cal for any color_* column
                for cand in cal.columns:
                    if cand.startswith("color_"):
                        color_col = cand
                        break
            if color_col is None or color_col not in cal.columns:
                continue

            ref_col = f"ref_{filt}"
            inst_col = f"mag_inst_{filt}"
            err_col = f"mag_inst_err_{filt}"
            delta_col = f"delta_{filt}"

            if delta_col in cal.columns:
                delta = cal[delta_col].to_numpy(float)
            elif ref_col in cal.columns and inst_col in cal.columns:
                delta = cal[ref_col].to_numpy(float) - cal[inst_col].to_numpy(float)
            else:
                continue

            color_x = cal[color_col].to_numpy(float)
            mask = np.isfinite(delta) & np.isfinite(color_x)
            if mask.sum() == 0:
                continue

            x_plot = color_x[mask]
            y_plot = delta[mask]
            stars_idx = np.where(mask)[0]

            err_arr = None
            if err_col in cal.columns:
                err_raw = cal[err_col].to_numpy(float)[mask]
                err_arr = np.where(np.isfinite(err_raw), err_raw, np.nan)

            # Sigma-clip outlier mask re-computed from saved fit line
            inlier = np.ones(len(x_plot), dtype=bool)
            zp_val = ct_val = None
            ct2_val = 0.0
            fit_label = f"{filt} (N={mask.sum()})"
            if self._coeff_df is not None:
                row = self._coeff_df[self._coeff_df["filter"] == filt]
                if len(row):
                    zp_val = float(row["zp"].iloc[0])
                    ct_val = float(row["ct"].iloc[0])
                    if "ct2" in row.columns and np.isfinite(row["ct2"].iloc[0]):
                        ct2_val = float(row["ct2"].iloc[0])
                    N_val = int(row["N"].iloc[0]) if "N" in row.columns else mask.sum()
                    sc_val = float(row["scatter_rms"].iloc[0]) if "scatter_rms" in row.columns else np.nan
                    resid = y_plot - (zp_val + ct_val * x_plot + ct2_val * x_plot**2)
                    med_r = np.nanmedian(resid)
                    mad_r = np.nanmedian(np.abs(resid - med_r)) + 1e-12
                    sig_r = MAD_TO_SIGMA * mad_r
                    inlier = np.abs(resid - med_r) <= 3.0 * sig_r
                    sc_str = f"σ={sc_val:.4f}" if np.isfinite(sc_val) else ""
                    ct2_str = f" CT2={ct2_val:+.3f}" if ct2_val else ""
                    fit_label = f"{filt}: ZP={zp_val:.3f} CT={ct_val:+.3f}{ct2_str} {sc_str} (N={N_val})"

            # Outliers first (behind inliers)
            if (~inlier).any():
                ax.scatter(x_plot[~inlier], y_plot[~inlier],
                           marker="x", c="gray", s=25, alpha=0.35, linewidths=1.0, zorder=2)

            sc = ax.scatter(
                x_plot[inlier], y_plot[inlier],
                c=fc, s=12, alpha=0.60,
                label=fit_label,
                picker=True, pickradius=6, zorder=3,
            )

            if err_arr is not None and np.isfinite(err_arr[inlier]).any():
                ax.errorbar(
                    x_plot[inlier], y_plot[inlier], yerr=err_arr[inlier],
                    fmt="none", ecolor=fc, elinewidth=0.7, capsize=2,
                    alpha=0.35, zorder=2,
                )

            row_list = []
            inlier_idx = np.where(inlier)[0]
            for i, si_local in enumerate(inlier_idx):
                si = stars_idx[si_local]
                row_list.append({
                    "filt": filt,
                    "ID": cal["ID"].iloc[si] if "ID" in cal.columns else "?",
                    "color": float(x_plot[si_local]),
                    "delta": float(y_plot[si_local]),
                    "ref_mag": float(cal[ref_col].iloc[si]) if ref_col in cal.columns else np.nan,
                    "inst_mag": float(cal[inst_col].iloc[si]) if inst_col in cal.columns else np.nan,
                    "err": float(err_arr[si_local]) if err_arr is not None and np.isfinite(err_arr[si_local]) else np.nan,
                })
            self._artist_map[id(sc)] = row_list

            if zp_val is not None and ct_val is not None:
                x_fit = np.linspace(np.nanmin(x_plot) - 0.05, np.nanmax(x_plot) + 0.05, 200)
                y_fit = zp_val + ct_val * x_fit + ct2_val * x_fit**2
                ax.plot(x_fit, y_fit, "-", color=fc, linewidth=2.0, zorder=4)

            color_labels.add(color_col)

        # X-axis label based on color columns used
        if len(color_labels) == 1:
            cc = next(iter(color_labels))
            color_label = cc.replace("color_", "").replace("_", " - ") + " (inst)"
        else:
            color_label = "color index (inst)"
        ax.set_xlabel(color_label)
        ax.set_ylabel("gaia_ref − inst (mag)")
        ax.set_title("ZP Linear Fit: gaia_ref − inst vs color  [× = sigma-clipped outliers]")
        if ax.collections or ax.get_lines():
            ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    def _draw_zp_hist(self):
        """ZP distribution histogram per filter (논문 AutoPHOT Fig.8 방식)."""
        from scipy.stats import gaussian_kde
        ax = self._ax_hist
        filt_sel = self._filt_combo.currentText()

        if self._frame_df is None:
            ax.set_title("ZP Distribution")
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    transform=ax.transAxes, fontsize=12, color="gray")
            return

        df = self._frame_df.copy()
        date_sel = self._date_combo.currentText()
        if date_sel != "All" and "file" in df.columns:
            import re
            df = df[df["file"].apply(lambda f: bool(re.search(re.escape(date_sel), str(f))))].copy()

        filts = (sorted(df["filter"].unique()) if "filter" in df.columns else []) if filt_sel == "All" else [filt_sel]
        has_any = False
        for filt in filts:
            sub = df[df["filter"] == filt] if "filter" in df.columns else pd.DataFrame()
            if sub.empty:
                continue
            zp_vals = pd.to_numeric(sub["zp_frame"], errors="coerce").dropna().to_numpy(float)
            zp_vals = zp_vals[np.isfinite(zp_vals)]
            if len(zp_vals) < 3:
                continue
            has_any = True
            fc = self._filter_color(filt)
            med = float(np.median(zp_vals))
            mad = float(np.median(np.abs(zp_vals - med)))
            sigma = MAD_TO_SIGMA * mad
            n_bins = max(8, min(30, len(zp_vals) // 2))
            ax.hist(zp_vals, bins=n_bins, color=fc, alpha=0.5,
                    label=f"{filt}: μ={med:.3f} σ={sigma:.3f} N={len(zp_vals)}")
            # KDE curve
            if len(zp_vals) >= 5:
                try:
                    kde = gaussian_kde(zp_vals, bw_method="scott")
                    xg = np.linspace(zp_vals.min() - 3 * sigma, zp_vals.max() + 3 * sigma, 200)
                    yk = kde(xg) * len(zp_vals) * (zp_vals.max() - zp_vals.min()) / n_bins
                    ax.plot(xg, yk, color=fc, linewidth=1.5)
                except Exception:
                    pass
            ax.axvline(med, color=fc, linestyle="--", linewidth=1.2, alpha=0.9)

        ax.set_xlabel("ZP (mag)")
        ax.set_ylabel("Count")
        ax.set_title("ZP Distribution")
        if has_any:
            ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    def _draw_frame_zp(self):
        import re
        ax = self._ax_frame
        filt_sel = self._filt_combo.currentText()
        date_sel = self._date_combo.currentText()

        if self._frame_df is None:
            ax.set_title("Per-Frame ZP — run ZP calibration first")
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes,
                    fontsize=14, color="gray")
            return

        df = self._frame_df.copy()
        if date_sel != "All" and "file" in df.columns:
            df = df[df["file"].apply(lambda f: bool(re.search(re.escape(date_sel), str(f))))].copy()

        if "filter" not in df.columns:
            ax.set_title("Per-Frame ZP — no filter column")
            return

        filts = (sorted(df["filter"].unique())) if filt_sel == "All" else [filt_sel]

        # Build shared x-tick labels from the first filter's file list
        _label_sub = None
        for filt in filts:
            sub0 = df[df["filter"] == filt].reset_index(drop=True)
            if not sub0.empty:
                _label_sub = sub0
                break

        for filt in filts:
            sub = df[df["filter"] == filt].reset_index(drop=True)
            if sub.empty:
                continue
            fc = self._filter_color(filt)
            x = np.arange(len(sub))
            y = sub["zp_frame"].to_numpy(float)
            yerr = sub["zp_scatter"].to_numpy(float) if "zp_scatter" in sub.columns else None

            # Point size ∝ n_ref
            has_nref = "n_ref" in sub.columns
            s_size = np.clip(sub["n_ref"].to_numpy(float) * 3.0, 12, 90) if has_nref else np.full(len(sub), 25)

            # Connection line
            ax.plot(x, y, "-", color=fc, alpha=0.45, linewidth=1.2, zorder=2)
            # Error bars (scatter)
            if yerr is not None and np.isfinite(yerr).any():
                ax.errorbar(x, y, yerr=yerr, fmt="none", ecolor=fc,
                           elinewidth=0.7, capsize=2, alpha=0.5, zorder=2)
            # Points with variable size
            n_label = int(sub["n_ref"].mean()) if has_nref else len(sub)
            ax.scatter(x, y, s=s_size, c=fc, alpha=0.85, zorder=3,
                      label=f"{filt}  N_frame={len(sub)}  <n_ref>≈{n_label}",
                      picker=True, pickradius=6)

            # Median ± σ band
            yf = y[np.isfinite(y)]
            if len(yf) >= 3:
                med_zp = float(np.median(yf))
                mad_zp = float(np.median(np.abs(yf - med_zp)))
                sigma_zp = MAD_TO_SIGMA * mad_zp
                ax.axhline(med_zp, color=fc, linestyle="--", linewidth=0.9, alpha=0.7)
                ax.axhspan(med_zp - sigma_zp, med_zp + sigma_zp,
                           color=fc, alpha=0.07, zorder=1)

        # Date-based x-tick labels from the first filter
        if _label_sub is not None and "file" in _label_sub.columns:
            n = len(_label_sub)
            tick_every = max(1, n // 18)
            tick_x = list(range(0, n, tick_every))
            tick_labels = []
            for i in tick_x:
                fname = str(_label_sub["file"].iloc[i])
                m = re.search(r"\d{4}-\d{2}-\d{2}", fname)
                if m:
                    tick_labels.append(m.group(0))
                else:
                    m2 = re.search(r"\d{8}", fname)
                    tick_labels.append(m2.group(0) if m2 else str(i))
            ax.set_xticks(tick_x)
            ax.set_xticklabels(tick_labels, rotation=35, ha="right", fontsize=7)
            xlabel = "Frame  (date from filename;  point size ∝ n_ref)"
        else:
            xlabel = "Frame index  (point size ∝ n_ref)"

        ax.set_xlabel(xlabel)
        ax.set_ylabel("ZP (mag)")
        ax.set_title("Per-Frame Zeropoint")
        if ax.get_lines() or ax.collections:
            ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    def _on_pick(self, event):
        artist = event.artist
        aid = id(artist)
        if aid not in self._artist_map:
            return
        rows = self._artist_map[aid]
        ind = event.ind
        if len(ind) == 0:
            return
        info = rows[ind[0]]
        parts = [
            f"Filter={info['filt']}",
            f"ID={info['ID']}",
            f"color={info['color']:.3f}",
            f"delta={info['delta']:.3f}",
            f"ref_mag={info['ref_mag']:.3f}",
            f"inst_mag={info['inst_mag']:.3f}",
        ]
        if np.isfinite(info["err"]):
            parts.append(f"err={info['err']:.4f}")
        self._info_label.setText(" | ".join(parts))

    def _save_png(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Plot", str(self.result_dir), "PNG Images (*.png)"
        )
        if path:
            self._fig.savefig(path, dpi=150, bbox_inches="tight")


class _AnchorDiscoveryWorker(QThread):
    """Background VizieR catalog discovery for the standard-anchor dialog.

    Network-bound (keyword search + per-candidate probe, ~10-60 s) so it must
    never run on the GUI thread; results come back via signals."""

    found = pyqtSignal(list)
    failed = pyqtSignal(str)

    def __init__(self, target_name, ra_deg, dec_deg, bands, parent=None):
        super().__init__(parent)
        self._args = (str(target_name), float(ra_deg), float(dec_deg), list(bands))

    def run(self):
        try:
            from apex.analysis.cmd.standard_anchor import discover_standard_catalogs

            name, ra, dec, bands = self._args
            self.found.emit(discover_standard_catalogs(name, ra, dec, bands))
        except Exception as exc:
            self.failed.emit(f"{type(exc).__name__}: {exc}")


class ZeropointCalibrationWindow(StepWindowBase):
    """Step 10: Zeropoint & Standardization"""

    def __init__(self, params, file_manager, project_state, main_window):
        self.file_manager = file_manager
        self.worker = None
        self.viewer = None
        self._current_zp_signature = None
        self._zp_cache_validation_result = None

        super().__init__(
            step_index=9,
            step_name="Zeropoint Calibration",
            params=params,
            project_state=project_state,
            main_window=main_window
        )

        self.setup_step_ui()
        self.restore_state()

    def setup_step_ui(self):
        # ── Controls ─────────────────────────────────────────────────────────
        info = QLabel("Build per-frame ZP calibration and standardized catalogs.")
        info.setProperty("role", "info")
        self.content_layout.addWidget(info)

        self.photometry_source_label = QLabel()
        self.photometry_source_label.setProperty("role", "caption")
        _f = self.photometry_source_label.font(); _f.setBold(True)
        self.photometry_source_label.setFont(_f)
        self.content_layout.addWidget(self.photometry_source_label)
        self._refresh_photometry_source_label()

        control_layout = QHBoxLayout()
        btn_params = create_parameter_button("Calibration Parameters")
        btn_params.clicked.connect(self.open_parameters_dialog)
        control_layout.addWidget(btn_params)

        self.run_bar = RunControlBar(
            "Run ZP Calibration", "Open Log",
            run_cb=self.run_analysis,
            stop_cb=self.stop_analysis,
            log_cb=self.show_log_window,
        )
        control_layout.addWidget(self.run_bar)
        self.btn_run = self.run_bar.btn_run
        self.btn_stop = self.run_bar.btn_stop
        self.content_layout.addLayout(control_layout)

        self.progress_label = QLabel("Idle")
        self.progress_label.setStyleSheet("QLabel { padding: 4px; }")
        self.content_layout.addWidget(self.progress_label)

        # ── ZP Fit Plot (takes remaining space) ──────────────────────────────
        self.fit_tab = ZPFitPlotWidget(self.params.P.result_dir)
        self.content_layout.addWidget(self.fit_tab, 1)

        # ── Log window with Workers panel ────────────────────────────────────
        _worker_group = QGroupBox("Workers")
        _worker_group.setMinimumWidth(300)
        _wg_layout = QVBoxLayout(_worker_group)
        _wg_layout.setContentsMargins(5, 5, 5, 5)
        self.worker_panel = WorkerStatusPanel(_worker_group)
        _wg_layout.addWidget(self.worker_panel)

        self.log_window = WorkflowLogWindow(
            self, "Calibration Log", width=850, height=420,
            side_widget=_worker_group,
        )
        self.log_text = self.log_window.log_text

    def _refresh_photometry_source_label(self, info: dict | None = None) -> None:
        provenance = info or resolve_cmd_photometry_input(
            self.params.P.result_dir,
            self.project_state,
        )
        self.photometry_source_label.setText(
            format_photometry_provenance(provenance)
        )

    def log(self, message: str):
        append_timestamped_log(self.log_text, message)

    def show_log_window(self):
        show_raised(self.log_window)

    def open_parameters_dialog(self):
        dialog, layout, buttons = build_scroll_param_dialog(
            self, "Calibration Parameters",
            info_text="Adjust zero-point calibration parameters. Changes apply to the next calibration run.",
            size=(500, 620),
        )

        match_group, match_container = create_collapsible_section("Matching", initial_expanded=True)
        match_form = QFormLayout(match_container)
        match_form.setContentsMargins(0, 0, 0, 0)

        self.param_pix = QDoubleSpinBox()
        self.param_pix.setRange(0.0, 50.0)
        self.param_pix.setValue(float(getattr(self.params.P, "pixel_scale_arcsec", 0.0) or 0.0))
        self.param_pix.setEnabled(False)
        match_form.addRow("Pixel scale (arcsec):", self.param_pix)

        self.param_match = QDoubleSpinBox()
        self.param_match.setRange(0.1, 20.0)
        self.param_match.setValue(float(getattr(self.params.P, "match_tol_px", 5.0)))
        match_form.addRow("Match tol (px):", self.param_match)

        self.param_min_match = QSpinBox()
        self.param_min_match.setRange(3, 1000)
        self.param_min_match.setValue(int(getattr(self.params.P, "min_master_gaia_matches", 10)))
        match_form.addRow("Min Gaia matches:", self.param_min_match)

        layout.addWidget(match_group)

        zp_group, zp_container = create_collapsible_section("Zero-point Fit", initial_expanded=True)
        zp_form = QFormLayout(zp_container)
        zp_form.setContentsMargins(0, 0, 0, 0)

        self.param_cmd_snr = QDoubleSpinBox()
        self.param_cmd_snr.setRange(0.0, 200.0)
        self.param_cmd_snr.setValue(float(getattr(self.params.P, "cmd_snr_calib_min", 20.0)))
        zp_form.addRow("CMD calib SNR min:", self.param_cmd_snr)

        self.param_frame_min = QSpinBox()
        self.param_frame_min.setRange(1, 1000)
        self.param_frame_min.setValue(int(getattr(self.params.P, "frame_zp_min_n", 5)))
        zp_form.addRow("Frame ZP min refs:", self.param_frame_min)

        self.param_apply_ext = QCheckBox("Enable")
        self.param_apply_ext.setChecked(bool(getattr(self.params.P, "cmd_apply_extinction", False)))
        zp_form.addRow("Apply extinction (k·X):", self.param_apply_ext)

        self.param_ext_mode = QComboBox()
        self.param_ext_mode.addItems(["absorb", "two_step"])
        self.param_ext_mode.setCurrentText(str(getattr(self.params.P, "cmd_extinction_mode", "absorb")))
        zp_form.addRow("Extinction mode:", self.param_ext_mode)

        self.param_clip = QDoubleSpinBox()
        self.param_clip.setRange(0.5, 10.0)
        self.param_clip.setValue(float(getattr(self.params.P, "zp_clip_sigma", 3.0)))
        zp_form.addRow("ZP clip sigma:", self.param_clip)

        self.param_iters = QSpinBox()
        self.param_iters.setRange(1, 20)
        self.param_iters.setValue(int(getattr(self.params.P, "zp_fit_iters", 5)))
        zp_form.addRow("ZP fit iters:", self.param_iters)

        self.param_slope = QDoubleSpinBox()
        self.param_slope.setRange(0.1, 5.0)
        self.param_slope.setValue(float(getattr(self.params.P, "zp_slope_absmax", 1.0)))
        zp_form.addRow("ZP slope abs max:", self.param_slope)

        layout.addWidget(zp_group)

        gaia_group, gaia_container = create_collapsible_section("Calibration Star Selection")
        gaia_form = QFormLayout(gaia_container)
        gaia_form.setContentsMargins(0, 0, 0, 0)

        self.param_gaia_snr = QDoubleSpinBox()
        self.param_gaia_snr.setRange(0.0, 200.0)
        self.param_gaia_snr.setValue(float(getattr(self.params.P, "gaia_snr_calib_min", 20.0)))
        self.param_gaia_snr.setToolTip("SNR threshold for calibration reference stars (applied in both global fit and per-frame ZP)")
        gaia_form.addRow("Calib star SNR min:", self.param_gaia_snr)

        self.param_gi_min = QDoubleSpinBox()
        self.param_gi_min.setRange(-2.0, 5.0)
        self.param_gi_min.setDecimals(2)
        self.param_gi_min.setValue(float(getattr(self.params.P, "gaia_gi_min", -0.5)))
        self.param_gi_min.setToolTip("Global lower bound for Gaia BP-RP color. Stars outside this range are excluded from all Jordi fits.")
        gaia_form.addRow("BP-RP min:", self.param_gi_min)

        self.param_gi_max = QDoubleSpinBox()
        self.param_gi_max.setRange(0.5, 6.0)
        self.param_gi_max.setDecimals(2)
        self.param_gi_max.setValue(float(getattr(self.params.P, "gaia_gi_max", 3.5)))
        self.param_gi_max.setToolTip("Global upper bound for Gaia BP-RP color. Stars outside this range are excluded from all Jordi fits.")
        gaia_form.addRow("BP-RP max:", self.param_gi_max)

        layout.addWidget(gaia_group)

        anchor_group, anchor_container = create_collapsible_section(
            "External Standard Anchor")
        anchor_form = QFormLayout(anchor_container)
        anchor_form.setContentsMargins(0, 0, 0, 0)

        anchor_info = QLabel(
            "Gaia 변환 참조가 부정확한 밴드(특히 Johnson U, σ≈0.20)를 "
            "VizieR 표준성 카탈로그(Landolt 계열, Gaia 독립)로 재앵커합니다.\n"
            "M67 실측: U 영점이 -0.13 mag 틀어져 이소크론 [M/H]가 "
            "-0.83으로 rail — 재앵커 후 +0.06(문헌 일치). 모든 밴드가 한 "
            "표준계에 앉아야 U-B가 축퇴를 풉니다."
        )
        anchor_info.setWordWrap(True)
        anchor_form.addRow(anchor_info)

        self.param_anchor_enable = QCheckBox("Enable")
        self.param_anchor_enable.setChecked(
            bool(getattr(self.params.P, "std_anchor_enable", False)))
        self.param_anchor_enable.setToolTip(
            "켜면 Step 10이 wide CMD 테이블 저장 직전에 mag_std_* 를 외부 "
            "표준성 오프셋만큼 이동합니다 (계측 측광은 불변). 오프셋 QC는 "
            "cmd_zeropoint/standard_anchor_offsets.csv 에 남습니다.")
        anchor_form.addRow("Anchor mag_std to standards:", self.param_anchor_enable)

        self.param_anchor_catalog = QLineEdit()
        self.param_anchor_catalog.setText(
            str(getattr(self.params.P, "std_anchor_catalog", "") or ""))
        self.param_anchor_catalog.setPlaceholderText(
            "VizieR catalog id — e.g. J/AJ/106/181 (M67 UBVRI, Montgomery+1993)")
        self.param_anchor_catalog.setToolTip(
            "대상 시야를 덮는 표준성 측광 카탈로그의 VizieR ID. 직접 "
            "<band>mag 컬럼 또는 Vmag+색(B-V, U-B, V-R, V-I)이 있어야 하며, "
            "Gaia에서 유도된 카탈로그는 쓰면 안 됩니다(독립성 상실).")
        anchor_form.addRow("VizieR catalog:", self.param_anchor_catalog)

        self.param_anchor_find = QPushButton("Find catalogs for this field (VizieR)")
        self.param_anchor_find.setToolTip(
            "대상 이름·좌표로 VizieR를 검색해 이 시야를 덮는 표준 측광 "
            "카탈로그 후보를 찾습니다 (네트워크, 10~60초). 후보를 고르면 "
            "위 칸이 채워집니다.")
        self.param_anchor_find.clicked.connect(self._discover_anchor_catalogs)
        anchor_form.addRow("", self.param_anchor_find)

        self.param_anchor_candidates = QComboBox()
        self.param_anchor_candidates.setVisible(False)
        self.param_anchor_candidates.activated.connect(
            self._pick_anchor_candidate)
        anchor_form.addRow("", self.param_anchor_candidates)

        self.param_anchor_radius = QDoubleSpinBox()
        self.param_anchor_radius.setRange(0.1, 10.0)
        self.param_anchor_radius.setDecimals(1)
        self.param_anchor_radius.setValue(
            float(getattr(self.params.P, "std_anchor_match_radius", 1.5)))
        anchor_form.addRow("Match radius (arcsec):", self.param_anchor_radius)

        self.param_anchor_min_stars = QSpinBox()
        self.param_anchor_min_stars.setRange(5, 1000)
        self.param_anchor_min_stars.setValue(
            int(getattr(self.params.P, "std_anchor_min_stars", 20)))
        anchor_form.addRow("Min matched stars:", self.param_anchor_min_stars)

        layout.addWidget(anchor_group)
        layout.addStretch(1)
        add_parameter_reset_button(
            buttons,
            [
                (self.param_match, 1.0),
                (self.param_min_match, 10),
                (self.param_cmd_snr, 50.0),
                (self.param_frame_min, 5),
                (self.param_apply_ext, False),
                (self.param_ext_mode, "absorb"),
                (self.param_clip, 3.0),
                (self.param_iters, 5),
                (self.param_slope, 1.0),
                (self.param_gaia_snr, 20.0),
                (self.param_gi_min, -0.5),
                (self.param_gi_max, 4.5),
                (self.param_anchor_enable, False),
                (self.param_anchor_radius, 1.5),
                (self.param_anchor_min_stars, 20),
            ],
        )
        buttons.accepted.connect(lambda: self.save_parameters(dialog))
        buttons.rejected.connect(dialog.reject)
        dialog.exec_()

    def _discover_anchor_catalogs(self):
        P = self.params.P
        ra = getattr(P, "target_ra_deg", None) or getattr(P, "ra_deg", None)
        dec = getattr(P, "target_dec_deg", None) or getattr(P, "dec_deg", None)
        if ra is None or dec is None:
            QMessageBox.warning(self, "Standard anchor",
                                "대상 좌표(target.ra_deg/dec_deg)가 설정에 없어 "
                                "탐색할 수 없습니다.")
            return
        name = str(getattr(P, "target_name", "") or "").strip()
        if not name:
            QMessageBox.warning(self, "Standard anchor",
                                "대상 이름(target.name)이 설정에 없어 키워드 "
                                "탐색을 할 수 없습니다.")
            return
        self.param_anchor_find.setEnabled(False)
        self.param_anchor_find.setText("Searching VizieR…")
        bands = ["U", "B", "V", "R", "I", "g", "r", "i"]
        self._anchor_discovery = _AnchorDiscoveryWorker(name, ra, dec, bands, self)
        self._anchor_discovery.found.connect(self._anchor_candidates_found)
        self._anchor_discovery.failed.connect(self._anchor_discovery_failed)
        self._anchor_discovery.start()

    def _anchor_discovery_reset_button(self):
        self.param_anchor_find.setEnabled(True)
        self.param_anchor_find.setText("Find catalogs for this field (VizieR)")

    def _anchor_candidates_found(self, candidates):
        self._anchor_discovery_reset_button()
        combo = self.param_anchor_candidates
        combo.clear()
        if not candidates:
            combo.setVisible(False)
            QMessageBox.information(
                self, "Standard anchor",
                "이 시야를 덮는 표준 측광 카탈로그를 VizieR에서 찾지 "
                "못했습니다. Gaia 참조만으로 진행하거나, 표준장 전이 "
                "(같은 밤 표준장 관측)나 분광 [M/H] prior를 고려하세요.")
            return
        self._anchor_candidate_ids = [c.catalog_id for c in candidates]
        for c in candidates:
            field = "시야 내" if c.in_field else "시야 밖?"
            combo.addItem(f"{c.catalog_id} — {'/'.join(c.bands)} ({field}) "
                          f"{c.description[:40]}")
        combo.setVisible(True)
        combo.showPopup()

    def _pick_anchor_candidate(self, index):
        ids = getattr(self, "_anchor_candidate_ids", [])
        if 0 <= index < len(ids):
            self.param_anchor_catalog.setText(ids[index])
            self.param_anchor_enable.setChecked(True)

    def _anchor_discovery_failed(self, message):
        self._anchor_discovery_reset_button()
        QMessageBox.warning(self, "Standard anchor",
                            f"VizieR 탐색 실패: {message}")

    def save_parameters(self, dialog):
        self.params.P.match_tol_px = self.param_match.value()
        self.params.P.min_master_gaia_matches = self.param_min_match.value()
        self.params.P.cmd_snr_calib_min = self.param_cmd_snr.value()
        self.params.P.frame_zp_min_n = self.param_frame_min.value()
        self.params.P.cmd_apply_extinction = self.param_apply_ext.isChecked()
        self.params.P.cmd_extinction_mode = self.param_ext_mode.currentText().strip()
        self.params.P.zp_clip_sigma = self.param_clip.value()
        self.params.P.zp_fit_iters = self.param_iters.value()
        self.params.P.zp_slope_absmax = self.param_slope.value()
        self.params.P.gaia_snr_calib_min = self.param_gaia_snr.value()
        self.params.P.gaia_gi_min = self.param_gi_min.value()
        self.params.P.gaia_gi_max = self.param_gi_max.value()
        self.params.P.std_anchor_enable = self.param_anchor_enable.isChecked()
        self.params.P.std_anchor_catalog = self.param_anchor_catalog.text().strip()
        self.params.P.std_anchor_match_radius = self.param_anchor_radius.value()
        self.params.P.std_anchor_min_stars = self.param_anchor_min_stars.value()
        self._zp_cache_validation_result = None
        self.save_state()
        saved = self.persist_params()
        msg = "Parameters saved to TOML." if saved else "Parameters saved (TOML save failed)."
        QMessageBox.information(dialog, "Success", msg)
        dialog.accept()

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
            return [ZeropointCalibrationWindow._signature_value(v) for v in value]
        if isinstance(value, dict):
            return {
                str(k): ZeropointCalibrationWindow._signature_value(v)
                for k, v in sorted(value.items(), key=lambda item: str(item[0]))
            }
        return str(value)

    @staticmethod
    def _file_signature(path: Path | None) -> dict | None:
        if path is None:
            return None
        try:
            path = Path(path)
            if not path.is_file():
                return None
            stat = path.stat()
            return {
                "path": str(path.resolve()),
                "size": int(stat.st_size),
                "mtime_ns": int(stat.st_mtime_ns),
            }
        except OSError:
            return None

    def _build_zp_output_signature(self) -> dict:
        result_dir = Path(self.params.P.result_dir)
        source_info = resolve_cmd_photometry_input(
            result_dir,
            self.project_state,
        )
        active_phot_dir = Path(source_info["directory"])
        upstream_paths: list[Path] = [
            step5_wcs_dir(result_dir) / "wcs_solve_summary.csv",
            Path(source_info["index_path"]),
        ]
        if source_info["source"] == "psf":
            upstream_paths.extend([
                step7_forced_phot_dir(result_dir) / "photometry_index.csv",
                active_phot_dir / "psf_output_signature.json",
            ])
            active_patterns = ("photometry_*.tsv",)
        else:
            active_patterns = ("photometry_*.tsv", "apcorr_summary.csv")
        for directory, patterns in (
            (active_phot_dir, active_patterns),
            (step9_selection_dir(result_dir), ("*.csv", "*.tsv", "*.json")),
            (tool_extinction_dir(result_dir), ("*.csv", "*.json")),
        ):
            if directory.exists():
                for pattern in patterns:
                    upstream_paths.extend(sorted(directory.glob(pattern)))

        frame_paths: list[Path] = []
        try:
            for filename in self.file_manager.get_file_list():
                frame_paths.append(Path(self.file_manager.get_file_path(filename)))
        except Exception:
            frame_paths = []

        def _unique_signatures(paths: list[Path]) -> list[dict]:
            signatures: list[dict] = []
            seen: set[str] = set()
            for path in paths:
                signature = self._file_signature(path)
                if not signature:
                    continue
                key = signature["path"]
                if key in seen:
                    continue
                seen.add(key)
                signatures.append(signature)
            return sorted(signatures, key=lambda item: item["path"])

        payload = {
            "signature_version": _ZP_SIGNATURE_VERSION,
            "step": "cmd_step10_zeropoint_calibration",
            "photometry_input": {
                "source": source_info["source"],
                "mag_column": source_info["mag_column"],
                "mag_error_column": source_info["mag_error_column"],
            },
            "params": {
                key: self._signature_value(getattr(self.params.P, key, None))
                for key in _ZP_SIGNATURE_PARAMS
            },
            "inputs": {
                "upstream": _unique_signatures(upstream_paths),
                "frames": _unique_signatures(frame_paths),
            },
        }
        encoded = json.dumps(payload, sort_keys=True, ensure_ascii=False, allow_nan=False)
        payload["signature_hash"] = hashlib.sha1(encoded.encode("utf-8")).hexdigest()
        return payload

    def _stored_zp_signature(self) -> dict | None:
        path = step10_zp_dir(self.params.P.result_dir) / _ZP_SIGNATURE_FILE
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            return data if isinstance(data, dict) else None
        except Exception:
            return None

    def _write_zp_signature(self, signature: dict) -> None:
        out_dir = step10_zp_dir(self.params.P.result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / _ZP_SIGNATURE_FILE
        temp_path = path.with_name(path.name + ".tmp")
        temp_path.write_text(
            json.dumps(signature, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False),
            encoding="utf-8",
        )
        temp_path.replace(path)

    def _remove_zp_signature(self) -> None:
        path = step10_zp_dir(self.params.P.result_dir) / _ZP_SIGNATURE_FILE
        try:
            path.unlink()
        except FileNotFoundError:
            pass
        except OSError:
            pass

    def _zp_cache_status(self) -> tuple[bool, str, dict | None]:
        if self._zp_cache_validation_result is not None:
            return self._zp_cache_validation_result
        stored = self._stored_zp_signature()
        if not stored:
            result = (False, "missing signature", None)
            self._zp_cache_validation_result = result
            return result
        if stored.get("signature_version") != _ZP_SIGNATURE_VERSION:
            result = (False, "signature version mismatch", None)
            self._zp_cache_validation_result = result
            return result
        current = self._build_zp_output_signature()
        if stored.get("signature_hash") != current.get("signature_hash"):
            result = (False, "signature hash mismatch", None)
            self._zp_cache_validation_result = result
            return result
        summary = self._existing_output_summary()
        if not summary:
            result = (False, "calibration output missing or empty", None)
            self._zp_cache_validation_result = result
            return result
        result = (True, "ok", summary)
        self._zp_cache_validation_result = result
        return result

    def run_analysis(self):
        if self.worker and self.worker.isRunning():
            return
        self._refresh_photometry_source_label()
        self.log_text.clear()
        self.progress_label.setText("Starting...")
        self._zp_cache_validation_result = None
        self._current_zp_signature = self._build_zp_output_signature()
        self._remove_zp_signature()

        self.worker = ZeropointCalibrationWorker(
            self.params,
            self.params.P.data_dir,
            self.params.P.result_dir,
            self.params.P.cache_dir,
        )
        self.worker.progress.connect(self.on_progress)
        self.worker.log.connect(self.log)
        self.worker.finished.connect(self.on_finished)
        self.worker.error.connect(self.on_error)

        self.run_bar.set_running(True)
        self.worker.start()
        self.show_log_window()

    def stop_analysis(self):
        if self.worker and self.worker.isRunning():
            self.worker.stop()
            self.log("Stop requested")

    def on_progress(self, current, total, filename):
        self.progress_label.setText(f"{current}/{total} | {filename}")
        if hasattr(self, "worker_panel"):
            pct = int(100 * current / max(1, total))
            self.worker_panel.update_worker(0, filename, f"{current}/{total}", pct)

    def on_finished(self, summary):
        self.run_bar.set_running(False)
        if summary.get("stopped"):
            self.progress_label.setText("Stopped")
            self.log("Analysis stopped")
        else:
            self.progress_label.setText("Done")
            self.log("ZP calibration complete")
            self._refresh_photometry_source_label(summary)
            if self._current_zp_signature and self._existing_output_summary():
                try:
                    self._write_zp_signature(self._current_zp_signature)
                    self.log("[ZP][CACHE] Output signature saved for future reuse.")
                except Exception as exc:
                    self.log(f"[ZP][CACHE] Signature write failed: {exc}")
            self._zp_cache_validation_result = None
            self.save_state()
            self.update_navigation_buttons()
            self.fit_tab.reload()
        self._current_zp_signature = None
        self._cleanup_worker()

    def on_error(self, message):
        self.run_bar.set_running(False)
        self.progress_label.setText("Error")
        self.log(f"ERROR: {message}")
        self._current_zp_signature = None
        self._cleanup_worker()

    def _cleanup_worker(self, timeout_ms=5000):
        if not self.worker:
            return True

        worker = self.worker
        if worker.isRunning():
            try:
                worker.stop()
            except Exception:
                pass
            worker.quit()
            if not worker.wait(int(timeout_ms)):
                self.log("Calibration worker is still running; close is deferred.")
                return False

        try:
            worker.deleteLater()
        except Exception:
            pass
        self.worker = None
        return True

    def closeEvent(self, event):
        if self.worker and self.worker.isRunning():
            self.stop_analysis()
        if not self._cleanup_worker(timeout_ms=10000):
            QMessageBox.warning(
                self,
                "Background Task Running",
                "Calibration worker is still stopping. Please wait a few seconds and close again.",
            )
            event.ignore()
            return
        super().closeEvent(event)

    def _existing_output_summary(self) -> dict | None:
        out_dir = step10_zp_dir(self.params.P.result_dir)
        wide_cmd = out_dir / "median_by_ID_filter_wide_cmd.csv"
        wide = out_dir / "median_by_ID_filter_wide.csv"
        if not wide_cmd.exists() and not wide.exists():
            legacy_cmd = self.params.P.result_dir / "median_by_ID_filter_wide_cmd.csv"
            legacy = self.params.P.result_dir / "median_by_ID_filter_wide.csv"
            wide_cmd = legacy_cmd
            wide = legacy
        main_path = wide_cmd if wide_cmd.exists() else wide if wide.exists() else None
        if main_path is None:
            return None

        n_sources = 0
        try:
            n_sources = len(pd.read_csv(main_path, nrows=1000000))
        except Exception:
            return None
        if n_sources <= 0:
            return None

        coeff_path = out_dir / "zp_fit_coefficients.csv"
        frame_path = out_dir / "frame_zeropoint.csv"
        cal_path = out_dir / "gaia_sdss_calibrator_by_ID.csv"
        n_coeff = 0
        n_frames = 0
        n_cal = 0
        try:
            if coeff_path.exists():
                n_coeff = len(pd.read_csv(coeff_path))
        except Exception:
            n_coeff = 0
        try:
            if frame_path.exists():
                n_frames = len(pd.read_csv(frame_path))
        except Exception:
            n_frames = 0
        try:
            if cal_path.exists():
                n_cal = len(pd.read_csv(cal_path))
        except Exception:
            n_cal = 0

        return {
            "main_path": str(main_path),
            "n_sources": int(n_sources),
            "n_coeff": int(n_coeff),
            "n_frames": int(n_frames),
            "n_calibrators": int(n_cal),
        }

    def _try_load_existing_results(self) -> bool:
        valid, reason, summary = self._zp_cache_status()
        if not valid or not summary:
            if self._existing_output_summary():
                try:
                    self.log(f"[ZP][CACHE] Previous output not restored ({reason}).")
                except Exception:
                    pass
            return False
        try:
            self.fit_tab.reload(self.params.P.result_dir)
        except Exception:
            pass
        try:
            export_zp_qc_products(step10_zp_dir(self.params.P.result_dir), self.log)
            export_cmd_qc_products(step10_zp_dir(self.params.P.result_dir), self.log)
            export_gaia_cmd_comparison_products(step10_zp_dir(self.params.P.result_dir), self.log)
        except Exception:
            pass
        parts = [f"{summary.get('n_sources', 0)} sources"]
        if summary.get("n_frames", 0):
            parts.append(f"{summary['n_frames']} frame ZPs")
        if summary.get("n_coeff", 0):
            parts.append(f"{summary['n_coeff']} fit coeffs")
        if summary.get("n_calibrators", 0):
            parts.append(f"{summary['n_calibrators']} calibrators")
        self.progress_label.setText("Loaded previous ZP calibration (" + ", ".join(parts) + ")")
        try:
            self.log("[ZP][CACHE] Loaded previous Step 10 ZP calibration from disk.")
        except Exception:
            pass
        self.update_navigation_buttons()
        return True

    def validate_step(self) -> bool:
        valid, _, _ = self._zp_cache_status()
        return valid

    def save_state(self):
        state_data = {
            "match_tol_px": getattr(self.params.P, "match_tol_px", 5.0),
            "min_master_gaia_matches": getattr(self.params.P, "min_master_gaia_matches", 10),
            "cmd_snr_calib_min": getattr(self.params.P, "cmd_snr_calib_min", 20.0),
            "frame_zp_min_n": getattr(self.params.P, "frame_zp_min_n", 5),
            "cmd_apply_extinction": getattr(self.params.P, "cmd_apply_extinction", False),
            "cmd_extinction_mode": getattr(self.params.P, "cmd_extinction_mode", "absorb"),
            "zp_clip_sigma": getattr(self.params.P, "zp_clip_sigma", 3.0),
            "zp_fit_iters": getattr(self.params.P, "zp_fit_iters", 5),
            "zp_slope_absmax": getattr(self.params.P, "zp_slope_absmax", 1.0),
            "gaia_snr_calib_min": getattr(self.params.P, "gaia_snr_calib_min", 20.0),
            "gaia_gi_min": getattr(self.params.P, "gaia_gi_min", -0.5),
            "gaia_gi_max": getattr(self.params.P, "gaia_gi_max", 4.5),
            "gaia_zp_slope_absmax": getattr(self.params.P, "gaia_zp_slope_absmax", 1.0),
            "gaia_color_slope_absmax": getattr(self.params.P, "gaia_color_slope_absmax", 2.0),
        }
        self.project_state.store_step_data("zeropoint_calibration", state_data)

    def restore_state(self):
        state_data = self.project_state.get_step_data("zeropoint_calibration")
        if not state_data:
            state_data = self.project_state.get_step_data("cmd_analysis")
        if state_data:
            for key, val in state_data.items():
                if key == "pixel_scale_arcsec":
                    continue
                if hasattr(self.params.P, key):
                    setattr(self.params.P, key, val)
        self._try_load_existing_results()
