"""
Variable Star Analysis Tool

Interactive variable-star workflow for Step 10/11 or merged-workspace light curves.

Workflow:
  1. Load a light curve workspace
  2. Run a quick LS/PDM/BLS scan
  3. Route to Single or Multi analysis
  4. Single path: refine/bootstrap -> phase/Fourier -> O-C
  5. Multi path: candidate periods -> multi-mode fit -> isolated mode views
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Optional, List

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from matplotlib import rcParams

_CORR_MODE_RE = re.compile(r"lightcurve_.*?_(global|color|offset|raw)\b", re.IGNORECASE)
_TARGET_ID_RE = re.compile(r"lightcurve_(?:combined_)?ID(\d+)_", re.IGNORECASE)
_CORR_MODE_LABELS = {
    "global": "Global ensemble",
    "color": "Color-dependent",
    "offset": "Nightly offset",
    "raw": "Raw",
}


def _detect_corr_mode_from_df(df: pd.DataFrame, filename: str) -> str:
    if "correction_mode" in df.columns:
        vals = df["correction_mode"].dropna().astype(str).str.strip().str.lower()
        if not vals.empty:
            key = vals.iloc[0]
            if key:
                return _CORR_MODE_LABELS.get(key, key)
    m = _CORR_MODE_RE.search(filename)
    if m:
        return _CORR_MODE_LABELS.get(m.group(1).lower(), m.group(1))
    return ""


def _detect_target_id_from_df(df: pd.DataFrame, filename: str) -> int | None:
    m = _TARGET_ID_RE.search(filename)
    if m:
        try:
            return int(m.group(1))
        except Exception:
            return None
    for col in ("target_id", "star_id", "ID"):
        if col not in df.columns:
            continue
        vals = pd.to_numeric(df[col], errors="coerce").dropna().astype(int)
        uniq = sorted(set(vals.tolist()))
        if len(uniq) == 1:
            return int(uniq[0])
    return None


def _is_primary_lightcurve_path(path: Path) -> bool:
    """Check-star curves are overlays, never target-series candidates."""
    return not path.name.lower().startswith("lightcurve_check")


def _collect_mag_options(df: pd.DataFrame, time_mask: np.ndarray, corr_tag: str = "") -> list[tuple[str, str, np.ndarray]]:
    """Prefer canonical raw/corrected differential columns."""
    options: list[tuple[str, str, np.ndarray]] = []

    if corr_tag == "Raw":
        for col in ("diff_mag_raw", "diff_mag"):
            if col in df.columns:
                arr = pd.to_numeric(df[col], errors="coerce").to_numpy(float)[time_mask]
                if np.any(np.isfinite(arr)):
                    return [(f"Raw: {col}", col, arr)]

    if "diff_mag_raw" in df.columns or "diff_mag_corr" in df.columns:
        for col, label in (("diff_mag_raw", "raw"), ("diff_mag_corr", "corrected")):
            if col in df.columns:
                arr = pd.to_numeric(df[col], errors="coerce").to_numpy(float)[time_mask]
                if np.any(np.isfinite(arr)):
                    options.append((label, col, arr))
        if options:
            return options

    fallback_raw = ["mag_raw", "raw_mag", "inst_mag", "mag"]
    fallback_corr = ["mag_corr", "corr_mag", "calibrated_mag", "mag_ensemble_corr"]
    for col in fallback_raw:
        if col in df.columns:
            arr = pd.to_numeric(df[col], errors="coerce").to_numpy(float)[time_mask]
            if np.any(np.isfinite(arr)):
                options.append((f"Raw: {col}", col, arr))
    for col in fallback_corr:
        if col in df.columns:
            arr = pd.to_numeric(df[col], errors="coerce").to_numpy(float)[time_mask]
            if np.any(np.isfinite(arr)):
                options.append((f"Corr: {col}", col, arr))
    return options


def _series_rank(corr_tag: str, mag_col: str, source_name: str) -> tuple[int, int, str]:
    mode_order = {
        "Global ensemble": 0,
        "Color-dependent": 1,
        "Nightly offset": 2,
        "Raw": 3,
    }
    corrected = any(x in mag_col for x in ("corr", "cal"))
    if corr_tag in {"Nightly offset", "Color-dependent"}:
        corr_order = 1 if corrected else 0
    else:
        corr_order = 0 if corrected else 1
    return mode_order.get(corr_tag, 9), corr_order, source_name.lower()


def _source_priority(source_name: str) -> tuple[int, int, str]:
    lower = source_name.lower()
    is_combined = 1 if "_combined_" in lower else 0
    is_current = 0 if "_current" in lower else 1
    return is_combined, is_current, lower


def _describe_series(corr_tag: str, mag_col: str) -> str:
    if corr_tag == "Raw":
        return "Raw"
    if corr_tag == "Nightly offset":
        return "Offset | corrected" if mag_col == "diff_mag_corr" else "Offset | raw"
    if corr_tag == "Color-dependent":
        return "Color | corrected" if mag_col == "diff_mag_corr" else "Color | raw"
    if corr_tag == "Global ensemble":
        return "Global | corrected" if mag_col == "diff_mag_corr" else "Global | raw"
    if "corr" in mag_col or "cal" in mag_col:
        return "Corrected"
    return "Raw"


def _collect_workspace_series(
    paths: list[Path],
    excluded_files: set[str] | None = None,
    should_stop=None,
) -> list[dict]:
    excluded_files = set(excluded_files or set())
    series_items: list[dict] = []
    for path in paths:
        if should_stop is not None and should_stop():
            return []
        if not _is_primary_lightcurve_path(path):
            continue
        df = pd.read_csv(path)
        if excluded_files and "file" in df.columns:
            df = df[~df["file"].astype(str).isin(excluded_files)].reset_index(drop=True)
        time_col = next(
            (
                column
                for column in (
                    "BJD_TDB", "BJD", "bjd", "HJD", "hjd", "JD", "jd", "time"
                )
                if column in df.columns
            ),
            None,
        )
        if time_col is None:
            continue
        time = pd.to_numeric(df[time_col], errors="coerce").to_numpy(float)
        time_mask = np.isfinite(time)
        if not np.any(time_mask):
            continue
        time = time[time_mask]
        error_col = next(
            (
                column
                for column in ("diff_err_corr", "diff_err", "mag_err", "err", "sigma")
                if column in df.columns
            ),
            None,
        )
        error = (
            pd.to_numeric(df[error_col], errors="coerce").to_numpy(float)[time_mask]
            if error_col
            else None
        )
        filter_col = next(
            (
                column
                for column in ("filter", "Filter", "FILTER", "band", "Band")
                if column in df.columns
            ),
            None,
        )
        filters = df[filter_col].astype(str).to_numpy()[time_mask] if filter_col else None
        night_col = next(
            (column for column in ("night_id", "night", "date") if column in df.columns),
            None,
        )
        night_ids = (
            df[night_col].astype(str).to_numpy()[time_mask]
            if night_col
            else infer_night_ids(time)
        )
        correction = _detect_corr_mode_from_df(df, path.name)
        target_id = _detect_target_id_from_df(df, path.name)
        for label, column, values in _collect_mag_options(
            df, time_mask, corr_tag=correction
        ):
            series_items.append(
                {
                    "key": f"{path.name}::{column}",
                    "time": time,
                    "mag": values,
                    "mag_col": column,
                    "mag_err": error,
                    "filters": filters,
                    "night_id": night_ids,
                    "source": path.name,
                    "corr_tag": correction,
                    "correction_preserves_nightly_baseline": (
                        column == "diff_mag_raw"
                        or correction not in {"Nightly offset", "Color-dependent"}
                    ),
                    "series_label": _describe_series(correction, column),
                    "target_id": target_id,
                }
            )

    multi_target = len(
        {
            item["target_id"]
            for item in series_items
            if item.get("target_id") is not None
        }
    ) > 1
    for item in series_items:
        if multi_target:
            target_id = item.get("target_id")
            item["combo_label"] = (
                f"ID{target_id} | {item['series_label']}"
                if target_id is not None
                else f"{item['source']} | {item['series_label']}"
            )
        else:
            item["combo_label"] = item["series_label"]

    unique_series: dict[str, dict] = {}
    for item in series_items:
        combo_label = item["combo_label"]
        previous = unique_series.get(combo_label)
        if previous is None or _source_priority(item["source"]) < _source_priority(
            previous["source"]
        ):
            unique_series[combo_label] = item
    output = list(unique_series.values())
    output.sort(
        key=lambda item: _series_rank(
            item["corr_tag"], item["mag_col"], item["source"]
        )
    )
    return output


def _target_plot_title(params, target_id: int | None = None) -> str:
    target_name = str(getattr(params.P, "target_name", "") or "").strip()
    if target_name:
        return target_name
    if target_id is not None:
        return f"Target ID {int(target_id)}"
    return "Target"


def _safe_tight_layout(fig: Figure) -> None:
    """Avoid hard crashes when hidden or zero-sized canvases break tight_layout."""
    try:
        canvas = getattr(fig, "canvas", None)
        if canvas is not None:
            try:
                width, height = canvas.get_width_height()
            except Exception:
                width, height = None, None
            if width is not None and height is not None and (width <= 2 or height <= 2):
                fig.subplots_adjust(
                    left=0.08,
                    right=0.98,
                    bottom=0.08,
                    top=0.95,
                    hspace=0.35,
                    wspace=0.25,
                )
                return
        fig.tight_layout()
    except Exception:
        try:
            fig.subplots_adjust(
                left=0.08,
                right=0.98,
                bottom=0.08,
                top=0.95,
                hspace=0.35,
                wspace=0.25,
            )
        except Exception:
            pass

from astropy.timeseries import LombScargle

from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
    QPushButton, QLabel, QDoubleSpinBox, QSpinBox,
    QCheckBox, QTabWidget, QTextEdit, QTableWidget, QTableWidgetItem,
    QHeaderView, QFileDialog, QSplitter, QMessageBox, QComboBox, QLineEdit,
    QColorDialog, QDialog, QGridLayout, QScrollArea, QFrame,
    QProgressBar, QSizePolicy,
)
from PyQt5.QtCore import Qt, QThread, QTimer, pyqtSignal
from PyQt5.QtGui import QColor

from apex.gui.layout_rules import FittedDialog, prevent_collapse, scroll_wrap, tame_canvas
from apex.gui.theme import Tokens, style_button


def _repolish(widget) -> None:
    """Re-evaluate the theme QSS after a dynamic property change."""
    style = widget.style()
    style.unpolish(widget)
    style.polish(widget)


def _set_role(widget, prop: str, value) -> None:
    """Set a theme role property (role= / banner= / status=) and repolish."""
    widget.setProperty(prop, value)
    _repolish(widget)


def _mono_note_style() -> str:
    """Result-callout look for monospace summaries (O-C fit, Fourier, MM).

    Built from the live Tokens so every theme keeps its own surface/border;
    no theme property exists for a mono QLabel, hence the one f-string.
    """
    return (
        f"QLabel {{ background: {Tokens.SURFACE_ALT}; color: {Tokens.TEXT}; "
        f"padding: {Tokens.GAP}px; border: 1px solid {Tokens.BORDER}; "
        f"border-radius: {Tokens.RADIUS_SM}px; "
        f"font-family: 'Cascadia Mono', 'Consolas', monospace; "
        f"font-size: {Tokens.FS_CAPTION}px; }}"
    )

rcParams["axes.unicode_minus"] = False

_FILTER_COLORS = [
    "#1E88E5", "#E53935", "#43A047", "#FB8C00",
    "#8E24AA", "#00ACC1", "#F06292", "#795548",
]
_NAMED_FILTER_COLORS = {
    "u": "#1f77b4", "b": "#1f77b4", "B": "#1f77b4",
    "g": "#2ca02c", "v": "#2ca02c", "V": "#2ca02c",
    "r": "#d62728", "R": "#d62728",
    "i": "#9467bd", "I": "#9467bd",
    "z": "#8c564b",
    "H": "#17becf", "J": "#bcbd22",
    "clear": "#7f7f7f", "l": "#7f7f7f", "unknown": "#7f7f7f",
}

def _filt_color(filt: str, idx: int) -> str:
    return _NAMED_FILTER_COLORS.get(filt, _FILTER_COLORS[idx % len(_FILTER_COLORS)])


def _resolve_check_filter(filters, selected_filter: str | None = None) -> str | None:
    if selected_filter and selected_filter != "__all__":
        return selected_filter
    if filters is None:
        return None
    unique_filters = sorted({str(f) for f in filters if str(f).strip() and str(f).lower() != "nan"})
    return unique_filters[0] if len(unique_filters) == 1 else None

from apex.gui.tools.tool_window_base import ToolWindowBase
from apex.gui.workflow.lc.step11_period_analysis import PeriodAnalysisWorker
from apex.analysis.light_curve.period_analysis_service import run_period_analysis
from apex.analysis.light_curve.period_io_service import median_align_by_filter
from apex.analysis.light_curve.variable_analysis_contract import (
    ValidatedLightCurveBundle,
    VariableAnalysisRequest,
    VariableAnalysisResult,
    coerce_validated_bundle,
)
from apex.analysis.light_curve.variable_analysis_service import (
    evaluate_fixed_period_fourier as _evaluate_service_fourier,
)
from apex.gui.workers.variable_analysis_worker import VariableAnalysisWorker
from apex.analysis.light_curve.period_alias_service import (
    classify_frequency_relation,
    compute_spectral_window,
    evaluate_multimode_result as _service_evaluate_multimode_result,
    fit_multimode_model as _service_fit_multimode_model,
    infer_night_ids,
    periods_are_window_aliases,
    search_multimode_alias_solutions,
)


def _load_check_star_for_plot(result_dir: Path, filt: str | None = None):
    """Load check star CSV from step10 output for plotting. Returns (check_id, df_or_None)."""
    try:
        from apex.gui.workflow.lc.step9_lightcurve_builder import _load_check_star_csv
        check_id, df = _load_check_star_csv(result_dir, filt=filt)
        return check_id, (df if not df.empty else None)
    except Exception:
        return None, None


def _pick_check_overlay_cols(df: pd.DataFrame, preferred_mag_col: str | None = None) -> tuple[str | None, str | None]:
    time_col = next((c for c in ["BJD_TDB", "BJD", "bjd", "HJD", "hjd", "JD", "jd", "time"] if c in df.columns), None)
    mag_candidates: list[str] = []
    preferred = str(preferred_mag_col or "").strip()
    if preferred:
        if preferred == "diff_mag_corr":
            mag_candidates.extend(["diff_mag_corr", "diff_mag_raw", "diff_mag", "mag"])
        elif preferred == "diff_mag_raw":
            mag_candidates.extend(["diff_mag_raw", "diff_mag_corr", "diff_mag", "mag"])
        else:
            mag_candidates.append(preferred)
    mag_candidates.extend(["diff_mag_corr", "diff_mag_raw", "diff_mag", "mag_ensemble_corr", "mag"])
    seen: set[str] = set()
    mag_col = None
    for col in mag_candidates:
        if col in seen:
            continue
        seen.add(col)
        if col in df.columns:
            mag_col = col
            break
    return time_col, mag_col


def _parse_period_list(text: str) -> list[float]:
    periods: list[float] = []
    seen: set[float] = set()
    for token in re.split(r"[\s,;/]+", str(text or "").strip()):
        if not token:
            continue
        try:
            value = float(token)
        except Exception:
            continue
        if not np.isfinite(value) or value <= 0:
            continue
        key = round(value, 12)
        if key in seen:
            continue
        seen.add(key)
        periods.append(float(value))
    return periods


def _format_period_list(periods: list[float]) -> str:
    return ", ".join(f"{float(p):.8f}" for p in periods if np.isfinite(p) and float(p) > 0)


def _scan_top_periods(results: dict | None, max_count: int = 3) -> list[float]:
    if not results:
        return []
    for key in ("raw_ls", "raw_pdm", "raw_bls"):
        data = results.get(key)
        if not data or "error" in data:
            continue
        raw_periods = data.get("top_periods") or []
        periods: list[float] = []
        for value in raw_periods:
            try:
                p = float(value)
            except Exception:
                continue
            if not np.isfinite(p) or p <= 0:
                continue
            if any(abs(p - prev) / max(abs(prev), 1e-12) < 1e-5 for prev in periods):
                continue
            periods.append(p)
            if len(periods) >= int(max_count):
                break
        return periods
    return []


def _period_to_frequency(period: float) -> float:
    try:
        value = float(period)
    except Exception:
        return np.nan
    if not np.isfinite(value) or value <= 0:
        return np.nan
    return 1.0 / value


def _is_1day_alias_period(p1: float, p2: float, tol: float = 0.01) -> bool:
    f1 = _period_to_frequency(p1)
    f2 = _period_to_frequency(p2)
    if not (np.isfinite(f1) and np.isfinite(f2)):
        return False
    for order in (1, 2):
        for sign in (-1.0, 1.0):
            alias = f1 + sign * float(order)
            if alias <= 0:
                continue
            if abs(alias - f2) / max(abs(f2), 1e-12) < tol:
                return True
    return False


def _is_near_harmonic_period(p1: float, p2: float, max_order: int = 6, tol: float = 0.02) -> bool:
    f1 = _period_to_frequency(p1)
    f2 = _period_to_frequency(p2)
    if not (np.isfinite(f1) and np.isfinite(f2)):
        return False
    hi = max(f1, f2)
    lo = min(f1, f2)
    if lo <= 0:
        return False
    ratio = hi / lo
    for order in range(2, int(max_order) + 1):
        if abs(ratio - float(order)) / float(order) < tol:
            return True
    return False


def _recommend_analysis_mode_from_scan(
    results: dict | None,
    alias_analysis: dict | None = None,
) -> tuple[str, str]:
    if not results:
        return "unknown", "Run a scan first."

    available = []
    for key in ("raw_ls", "raw_pdm", "raw_bls"):
        data = results.get(key)
        if not data or "error" in data:
            continue
        best = float(data.get("best_period", np.nan))
        if np.isfinite(best) and best > 0:
            available.append((key, data))
    if not available:
        return "unknown", "No usable scan result."

    primary_key, primary = available[0]
    primary_label = primary_key.split("_", 1)[1].upper()
    top_periods = [float(p) for p in primary.get("top_periods", []) if np.isfinite(p) and float(p) > 0]
    top_powers = [float(v) for v in primary.get("top_powers", []) if np.isfinite(v)]
    if not top_periods:
        return "unknown", "No period peaks found."
    best_period = top_periods[0]
    best_power = top_powers[0] if top_powers else float(primary.get("best_power", np.nan))
    window_peaks = (alias_analysis or {}).get("window_peaks", [])
    baseline_days = float((alias_analysis or {}).get("baseline_days", 1.0))

    for idx, cand_period in enumerate(top_periods[1:], start=1):
        cand_power = top_powers[idx] if idx < len(top_powers) else np.nan
        if periods_are_window_aliases(best_period, cand_period, window_peaks, baseline_days):
            continue
        if not window_peaks and _is_1day_alias_period(best_period, cand_period):
            continue
        if _is_near_harmonic_period(best_period, cand_period):
            continue
        if np.isfinite(cand_power) and np.isfinite(best_power) and best_power > 0:
            rel_power = cand_power / best_power
            if rel_power >= 0.45:
                return "multi", (
                    f"{primary_label} scan shows a comparable secondary non-alias peak "
                    f"({cand_period:.6f} d, {rel_power:.0%} of primary power)."
                )
        else:
            return "multi", f"{primary_label} scan shows a secondary non-alias peak at {cand_period:.6f} d."

    best_periods = [float(data["best_period"]) for _, data in available]
    unique_periods: list[float] = []
    for value in best_periods:
        if any(abs(value - prev) / max(abs(prev), 1e-12) < 0.01 for prev in unique_periods):
            continue
        unique_periods.append(value)
    if len(unique_periods) >= 2:
        for idx in range(len(unique_periods)):
            for jdx in range(idx + 1, len(unique_periods)):
                p1 = unique_periods[idx]
                p2 = unique_periods[jdx]
                if periods_are_window_aliases(p1, p2, window_peaks, baseline_days):
                    continue
                if (not window_peaks and _is_1day_alias_period(p1, p2)) or _is_near_harmonic_period(p1, p2):
                    continue
                return "multi", (
                    f"Different methods prefer distinct non-alias periods ({p1:.6f} d vs {p2:.6f} d)."
                )

    return "single", f"Current scan is dominated by one primary period near {best_period:.6f} d."


def _handoff_mag_columns_match(requested: str, actual: str, correction_mode: str = "") -> bool:
    requested_key = str(requested or "").strip().lower()
    actual_key = str(actual or "").strip().lower()
    if not requested_key:
        return True
    if requested_key == actual_key:
        return True
    mode_key = str(correction_mode or "").strip().lower()
    if mode_key in {"raw", "raw (no correction)"}:
        return {requested_key, actual_key} <= {"diff_mag", "diff_mag_raw"}
    return False


def _select_period_handoff_series_key(
    series_options: dict[str, dict],
    handoff: dict,
) -> str | None:
    """Choose the closest loaded series without silently crossing targets."""
    candidates = list(series_options.items())
    if not candidates:
        return None

    requested_target = handoff.get("target_id")
    if requested_target not in (None, 0, ""):
        exact_target = [
            pair for pair in candidates
            if pair[1].get("target_id") is not None
            and int(pair[1]["target_id"]) == int(requested_target)
        ]
        if exact_target:
            candidates = exact_target

    requested_source = Path(str(handoff.get("source_file", "") or "")).name.lower()
    if requested_source:
        exact_source = [
            pair for pair in candidates
            if str(pair[1].get("source", "")).lower() == requested_source
        ]
        if exact_source:
            candidates = exact_source

    requested_mag = str(handoff.get("mag_col", "") or "")
    correction_mode = str(handoff.get("correction_mode", "") or "")
    exact_mag = [
        pair for pair in candidates
        if _handoff_mag_columns_match(
            requested_mag,
            str(pair[1].get("mag_col", "")),
            correction_mode,
        )
    ]
    if exact_mag:
        candidates = exact_mag

    candidates.sort(
        key=lambda pair: _series_rank(
            str(pair[1].get("corr_tag", "")),
            str(pair[1].get("mag_col", "")),
            str(pair[1].get("source", "")),
        )
    )
    return str(candidates[0][0])


def _route_period_handoff(
    scan_results: dict,
    alias_analysis: dict | None,
    multimode_diagnostic: dict | None,
) -> tuple[str, str, str]:
    """Return (recommended mode, workflow step, user-facing reason)."""
    mode, reason = _recommend_analysis_mode_from_scan(scan_results, alias_analysis)
    diagnostic_status = str((multimode_diagnostic or {}).get("status", "")).upper()
    diagnostic_reason = str((multimode_diagnostic or {}).get("reason", "")).strip()
    if diagnostic_status == "MULTIMODE-SUSPECT":
        mode = "multi"
        reason = diagnostic_reason or "Step 12 detected a possible independent residual mode."
    elif diagnostic_status == "SINGLE-COMPATIBLE":
        mode = "single"
        reason = diagnostic_reason or "Step 12 found no significant independent residual mode."

    alias_status = str((alias_analysis or {}).get("status", "UNASSESSED")).upper()
    candidates = list((alias_analysis or {}).get("candidates", []))
    if alias_status != "RESOLVED":
        count_text = f"{len(candidates)} candidate(s)" if candidates else "period candidates"
        alias_reason = str((alias_analysis or {}).get("reason", "")).strip()
        review_reason = f"Alias status is {alias_status}; keep {count_text} under review."
        if alias_reason:
            review_reason += f" {alias_reason}"
        return mode, "scan", review_reason

    workflow_step = mode if mode in {"single", "multi"} else "scan"
    return mode, workflow_step, reason


def _classify_candidate_period(
    period: float,
    adopted_periods: list[float],
    window_peaks: list[dict] | None = None,
    baseline_days: float = 1.0,
) -> tuple[str, str]:
    if any(abs(float(period) - float(prev)) / max(abs(float(prev)), 1e-12) < 1e-5 for prev in adopted_periods):
        return "duplicate", "already adopted"
    return classify_frequency_relation(
        _period_to_frequency(period),
        [_period_to_frequency(prev) for prev in adopted_periods],
        window_peaks=window_peaks,
        baseline_days=baseline_days,
    )


def _find_harmonic_order(target_freq: float, base_freq: float, max_order: int = 6, tol: float = 0.02) -> int | None:
    if not (np.isfinite(target_freq) and np.isfinite(base_freq)):
        return None
    if target_freq <= 0 or base_freq <= 0:
        return None
    for order in range(2, int(max_order) + 1):
        if abs(target_freq - float(order) * base_freq) / max(abs(target_freq), 1e-12) < tol:
            return int(order)
    return None


def _find_combination_relation(
    target_freq: float,
    prev_freqs: list[float],
    max_coeff: int = 2,
    tol: float = 0.01,
) -> str | None:
    if not np.isfinite(target_freq) or target_freq <= 0:
        return None
    best_match: tuple[float, str] | None = None
    for jdx in range(len(prev_freqs)):
        f_j = float(prev_freqs[jdx])
        if not np.isfinite(f_j) or f_j <= 0:
            continue
        for kdx in range(jdx + 1, len(prev_freqs)):
            f_k = float(prev_freqs[kdx])
            if not np.isfinite(f_k) or f_k <= 0:
                continue
            for coeff_j in range(1, int(max_coeff) + 1):
                for coeff_k in range(1, int(max_coeff) + 1):
                    for sign in (-1, 1):
                        combo = coeff_j * f_j + sign * coeff_k * f_k
                        if combo <= 0:
                            continue
                        rel = abs(target_freq - combo) / max(abs(target_freq), 1e-12)
                        if rel >= tol:
                            continue
                        expr = (
                            f"{coeff_j}f(M{jdx + 1}) {'+' if sign > 0 else '-'} "
                            f"{coeff_k}f(M{kdx + 1})"
                        )
                        if best_match is None or rel < best_match[0]:
                            best_match = (rel, expr)
    return None if best_match is None else f"near {best_match[1]}"


def _classify_fitted_periods(
    periods: list[float],
    window_peaks: list[dict] | None = None,
    baseline_days: float = 1.0,
) -> list[dict]:
    labels: list[dict] = []
    freqs = [_period_to_frequency(period) for period in periods]
    for idx, period in enumerate(periods):
        freq = freqs[idx]
        relation = "independent"
        if idx == 0:
            note = "primary fitted mode"
            labels.append({"relation": relation, "note": note})
            continue

        relation, note = classify_frequency_relation(
            freq,
            freqs[:idx],
            window_peaks=window_peaks,
            baseline_days=baseline_days,
        )
        if relation == "new":
            relation = "independent"
        if relation == "independent":
            for prev_idx in range(idx):
                order = _find_harmonic_order(freq, freqs[prev_idx])
                if order is not None:
                    relation = "harmonic"
                    note = f"near {order}f(M{prev_idx + 1})"
                    break
        if relation == "independent" and idx >= 2:
            combo_note = _find_combination_relation(freq, freqs[:idx])
            if combo_note:
                relation = "combination"
                note = combo_note
        labels.append({"relation": relation, "note": note})
    return labels


def _sanitize_filename_token(text: str) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", str(text or "").strip())
    token = token.strip("._")
    return token or "series"


def _build_multimode_design_matrix(time_rel: np.ndarray, periods: list[float], harmonics: int) -> tuple[np.ndarray, list[dict]]:
    cols = [np.ones(len(time_rel), dtype=float)]
    terms: list[dict] = []
    for mode_index, period in enumerate(periods):
        base_freq = 1.0 / float(period)
        for harmonic in range(1, int(harmonics) + 1):
            omega = 2.0 * np.pi * base_freq * harmonic
            cols.append(np.cos(omega * time_rel))
            terms.append({
                "mode_index": int(mode_index),
                "harmonic": int(harmonic),
                "kind": "cos",
                "omega": float(omega),
            })
            cols.append(np.sin(omega * time_rel))
            terms.append({
                "mode_index": int(mode_index),
                "harmonic": int(harmonic),
                "kind": "sin",
                "omega": float(omega),
            })
    return np.column_stack(cols), terms


def _evaluate_multimode_result(result: dict, times: np.ndarray) -> dict:
    return _service_evaluate_multimode_result(result, times)


def _fit_multimode_model(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: np.ndarray | None,
    periods: list[float],
    harmonics: int = 1,
    night_id: np.ndarray | None = None,
    include_night_offsets: bool = False,
) -> dict:
    return _service_fit_multimode_model(
        time,
        mag,
        mag_err,
        periods=periods,
        harmonics=harmonics,
        night_id=night_id,
        include_night_offsets=bool(include_night_offsets),
    )


def _trend_label(dm_dt: float) -> str:
    if not np.isfinite(dm_dt):
        return "unknown"
    if abs(dm_dt) < 1e-6:
        return "turning"
    return "brightening" if dm_dt < 0 else "fading"


def _fit_fixed_period_fourier(
    time: np.ndarray,
    mag: np.ndarray,
    period: float,
    harmonics: int,
    mag_err: np.ndarray | None = None,
) -> dict:
    t = np.asarray(time, dtype=float)
    y = np.asarray(mag, dtype=float)
    min_points = max(8, 2 * int(harmonics) + 3)
    mask = np.isfinite(t) & np.isfinite(y)
    dy = None
    if mag_err is not None:
        dy_all = np.asarray(mag_err, dtype=float)
        if len(dy_all) != len(t):
            raise ValueError("Magnitude-error array length does not match data")
        weighted_mask = mask & np.isfinite(dy_all) & (dy_all > 0)
        if np.count_nonzero(weighted_mask) >= min_points:
            mask = weighted_mask
            dy = dy_all[mask]
    t = t[mask]
    y = y[mask]
    if len(t) < min_points:
        raise ValueError("Not enough points for phase fit")
    time_ref = float(np.nanmin(t))
    tau = t - time_ref
    omega = 2.0 * np.pi / float(period)
    design = np.column_stack(
        [np.ones(len(t))] +
        [f(k * omega * tau) for k in range(1, int(harmonics) + 1) for f in (np.cos, np.sin)]
    )
    if dy is not None:
        coeff, _, _, _ = np.linalg.lstsq(design / dy[:, None], y / dy, rcond=None)
    else:
        coeff, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    return {
        "coeff": np.asarray(coeff, dtype=float),
        "time_ref": time_ref,
        "period": float(period),
        "weighted": dy is not None,
    }


def _evaluate_fixed_period_fourier(
    time: np.ndarray,
    period: float,
    fit_result: dict,
) -> np.ndarray:
    t = np.asarray(time, dtype=float)
    coeff = np.asarray(fit_result["coeff"], dtype=float)
    time_ref = float(fit_result.get("time_ref", 0.0))
    tau = t - time_ref
    harmonics = max(0, (len(coeff) - 1) // 2)
    model = np.full(len(t), float(coeff[0]), dtype=float)
    if harmonics <= 0:
        return model
    omega = 2.0 * np.pi / float(period)
    idx = 1
    for k in range(1, harmonics + 1):
        model += float(coeff[idx]) * np.cos(k * omega * tau)
        model += float(coeff[idx + 1]) * np.sin(k * omega * tau)
        idx += 2
    return model


def _fourier_shape_parameters(coeff: np.ndarray) -> dict:
    """Return amplitudes and phases for A_k cos(k*omega*tau + phi_k)."""
    values = np.asarray(coeff, dtype=float)
    if len(values) < 3 or (len(values) - 1) % 2:
        raise ValueError("Invalid Fourier coefficient vector")
    cos_coeff = values[1::2]
    sin_coeff = values[2::2]
    amplitudes = np.hypot(cos_coeff, sin_coeff)
    # a*cos(theta) + b*sin(theta) = A*cos(theta + phi), phi=atan2(-b, a)
    phases = np.arctan2(-sin_coeff, cos_coeff)

    def _wrap(value: float) -> float:
        return float((value + np.pi) % (2 * np.pi) - np.pi)

    r21 = amplitudes[1] / amplitudes[0] if len(amplitudes) > 1 and amplitudes[0] > 0 else np.nan
    r31 = amplitudes[2] / amplitudes[0] if len(amplitudes) > 2 and amplitudes[0] > 0 else np.nan
    phi21 = _wrap(phases[1] - 2 * phases[0]) if len(phases) > 1 else np.nan
    phi31 = _wrap(phases[2] - 3 * phases[0]) if len(phases) > 2 else np.nan
    return {
        "amplitudes": amplitudes,
        "phases": phases,
        "r21": float(r21),
        "r31": float(r31),
        "phi21": float(phi21),
        "phi31": float(phi31),
    }


def _detect_phase_epoch(
    time: np.ndarray,
    mag: np.ndarray,
    period: float,
    event_kind: str = "max_light",
    n_bins: int = 30,
) -> tuple[float, float]:
    """Estimate a maximum-light or minimum-light epoch in the first cycle."""
    t = np.asarray(time, dtype=float)
    y = np.asarray(mag, dtype=float)
    mask = np.isfinite(t) & np.isfinite(y)
    t = t[mask]
    y = y[mask]
    if len(t) < 5:
        raise ValueError("Not enough valid points to detect an epoch")
    if not np.isfinite(period) or period <= 0:
        raise ValueError("Period must be positive")
    if event_kind not in {"max_light", "min_light"}:
        raise ValueError(f"Unknown epoch event: {event_kind}")

    t_ref = float(np.nanmin(t))
    phase = ((t - t_ref) / float(period)) % 1.0
    edges = np.linspace(0.0, 1.0, int(n_bins) + 1)
    centers = (edges[:-1] + edges[1:]) / 2.0
    bin_mags = np.full(int(n_bins), np.nan)
    for bin_index in range(int(n_bins)):
        selected = (phase >= edges[bin_index]) & (phase < edges[bin_index + 1])
        if np.count_nonzero(selected) >= 2:
            bin_mags[bin_index] = np.nanmedian(y[selected])
    valid_indices = np.flatnonzero(np.isfinite(bin_mags))
    if len(valid_indices) < 3:
        raise ValueError("Phase coverage is insufficient to detect an epoch")

    if event_kind == "max_light":
        event_index = int(valid_indices[np.argmin(bin_mags[valid_indices])])
        expected_curvature = 1.0  # brightest means a minimum in magnitude
    else:
        event_index = int(valid_indices[np.argmax(bin_mags[valid_indices])])
        expected_curvature = -1.0  # faintest means a maximum in magnitude
    event_phase = float(centers[event_index])

    half_window = 4
    fit_indices = np.array(
        [(event_index + delta) % int(n_bins) for delta in range(-half_window, half_window + 1)]
    )
    local_phase = event_phase + ((centers[fit_indices] - event_phase + 0.5) % 1.0 - 0.5)
    local_mag = bin_mags[fit_indices]
    fit_mask = np.isfinite(local_mag)
    if np.count_nonzero(fit_mask) >= 3:
        try:
            parabola = np.polyfit(local_phase[fit_mask], local_mag[fit_mask], 2)
            if parabola[0] * expected_curvature > 0:
                candidate = float(-parabola[1] / (2.0 * parabola[0]))
                phase_distance = abs((candidate - event_phase + 0.5) % 1.0 - 0.5)
                if phase_distance < 0.15:
                    event_phase = candidate % 1.0
        except Exception:
            pass

    return float(t_ref + event_phase * period), float(event_phase)


# ---------------------------------------------------------------------------
# Worker: fine-grid refinement + bootstrap (LS-based)
# ---------------------------------------------------------------------------

class _WorkspaceLoadWorker(QThread):
    result_ready = pyqtSignal(int, object)
    error = pyqtSignal(int, str)

    def __init__(
        self,
        owner,
        token: int,
        paths: list[Path],
        excluded_files: set[str],
    ):
        super().__init__(owner)
        self.token = int(token)
        self.paths = [Path(path) for path in paths]
        self.excluded_files = set(excluded_files)

    def stop(self):
        self.requestInterruption()

    def run(self):
        try:
            items = _collect_workspace_series(
                self.paths,
                self.excluded_files,
                should_stop=self.isInterruptionRequested,
            )
            if not self.isInterruptionRequested():
                self.result_ready.emit(self.token, items)
        except Exception as exc:
            if not self.isInterruptionRequested():
                self.error.emit(self.token, str(exc))


class RefineBootstrapWorker(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)

    def __init__(self, time, mag, mag_err, center_period,
                 n_bootstrap: int = 300, zoom_factor: int = 100,
                 method: str = "ls", pdm_n_bins: int = 10):
        super().__init__()
        self.time = np.asarray(time, dtype=float)
        self.mag = np.asarray(mag, dtype=float)
        self.mag_err = np.asarray(mag_err, dtype=float) if mag_err is not None else None
        self.center_period = float(center_period)
        self.n_bootstrap = int(n_bootstrap)
        self.zoom_factor = int(zoom_factor)
        self.method = method  # "ls" or "pdm"
        self.pdm_n_bins = pdm_n_bins
        self._stop = False

    def stop(self):
        self._stop = True

    # -- PDM helpers ----------------------------------------------------------

    def _pdm_theta_array(self, t, y, trial_periods):
        """Compute PDM theta for an array of trial periods (vectorized).

        Uses bincount to avoid Python bin loops — ~10-50x faster than naive.
        """
        var_total = np.var(y)
        if var_total == 0:
            return np.ones(len(trial_periods))

        n_bins = self.pdm_n_bins
        t_min = t.min()
        dt = t - t_min
        y2 = y * y
        n_periods = len(trial_periods)
        theta = np.ones(n_periods)

        for i in range(n_periods):
            phase = (dt / trial_periods[i]) % 1.0
            bi = np.clip((phase * n_bins).astype(np.int32), 0, n_bins - 1)
            counts = np.bincount(bi, minlength=n_bins)
            sums = np.bincount(bi, weights=y, minlength=n_bins)
            sum_sq = np.bincount(bi, weights=y2, minlength=n_bins)
            # within-bin SS = sum(x²) - sum(x)²/n, dof = n-1
            good = counts >= 2
            if not good.any():
                continue
            c_g = counts[good]
            ss = sum_sq[good] - sums[good] ** 2 / c_g  # sum of squares
            dof = c_g - 1
            theta[i] = ss.sum() / dof.sum() / var_total
        return theta

    # -- Main run -------------------------------------------------------------

    def run(self):
        try:
            t, y, dy = self._filter_valid()
            if len(t) < 10:
                self.error.emit("Not enough valid data points (< 10)")
                return

            baseline = t.max() - t.min()
            f_center = 1.0 / self.center_period
            df_coarse = 1.0 / (10.0 * baseline)
            df_fine = df_coarse / self.zoom_factor
            half_range = 10.0 * df_coarse
            f_fine = np.arange(f_center - half_range, f_center + half_range + df_fine, df_fine)
            f_fine = f_fine[f_fine > 0]
            p_fine = 1.0 / f_fine[::-1]  # period grid (ascending)

            self.progress.emit(f"Fine grid search ({self.method.upper()})…")

            if self.method == "pdm":
                theta_fine = self._pdm_theta_array(t, y, p_fine)
                power_fine = 1.0 - theta_fine  # higher = better
                refined_period = self._parabola_peak_period(p_fine, power_fine)
            else:
                ls = (LombScargle(t, y, dy) if (dy is not None and np.any(dy > 0))
                      else LombScargle(t, y))
                power_freq = ls.power(f_fine)
                power_fine = power_freq[::-1]  # match p_fine order
                refined_period = self._parabola_peak(f_fine, power_freq)

            boot_periods = []
            n_data = len(t)
            template_harmonics = min(4, max(1, (n_data - 3) // 2))
            template_fit = _fit_fixed_period_fourier(
                t,
                y,
                refined_period,
                template_harmonics,
                mag_err=dy,
            )
            template_model = _evaluate_fixed_period_fourier(t, refined_period, template_fit)
            residuals = y - template_model
            residuals = residuals - np.nanmedian(residuals)
            rng = np.random.default_rng()
            for i in range(self.n_bootstrap):
                if self._stop:
                    break
                if i % 20 == 0:
                    self.progress.emit(f"Bootstrap {i}/{self.n_bootstrap}…")
                # Keep the timestamps fixed so every trial retains the real
                # sampling window. This measures precision around one alias.
                yb = template_model + rng.choice(residuals, n_data, replace=True)

                if self.method == "pdm":
                    theta_b = self._pdm_theta_array(t, yb, p_fine)
                    pwr_b = 1.0 - theta_b
                    boot_periods.append(self._parabola_peak_period(p_fine, pwr_b))
                else:
                    ls_b = (LombScargle(t, yb, dy) if (dy is not None and np.any(dy > 0))
                            else LombScargle(t, yb))
                    pwr_b = ls_b.power(f_fine)
                    boot_periods.append(self._parabola_peak(f_fine, pwr_b))

            boot_periods = np.array(boot_periods)
            med = np.median(boot_periods)
            mad = np.median(np.abs(boot_periods - med))
            keep = np.abs(boot_periods - med) < 5 * 1.4826 * mad
            sigma_p = (float(np.std(boot_periods[keep])) if keep.sum() >= 5
                       else float(np.std(boot_periods)))

            self.finished.emit({
                "refined_period": float(refined_period),
                "sigma_p": sigma_p,
                "boot_periods": boot_periods,
                "fine_periods": p_fine,
                "fine_power": power_fine,
                "method": self.method,
                "uncertainty_scope": "local_conditional",
            })
        except Exception as e:
            import traceback
            self.error.emit(f"{e}\n{traceback.format_exc()}")

    def _parabola_peak(self, freq, power):
        """Parabola interpolation in frequency space → return period."""
        idx = int(np.argmax(power))
        if 0 < idx < len(power) - 1:
            f3 = freq[idx - 1: idx + 2]
            p3 = power[idx - 1: idx + 2]
            try:
                coeffs = np.polyfit(f3, p3, 2)
                if coeffs[0] < 0:
                    fv = -coeffs[1] / (2.0 * coeffs[0])
                    if f3[0] < fv < f3[2]:
                        return 1.0 / fv
            except Exception:
                pass
        return 1.0 / freq[idx]

    def _parabola_peak_period(self, periods, power):
        """Parabola interpolation in period space → return period."""
        idx = int(np.argmax(power))
        if 0 < idx < len(power) - 1:
            p3 = periods[idx - 1: idx + 2]
            pw3 = power[idx - 1: idx + 2]
            try:
                coeffs = np.polyfit(p3, pw3, 2)
                if coeffs[0] < 0:
                    pv = -coeffs[1] / (2.0 * coeffs[0])
                    if p3[0] < pv < p3[2]:
                        return float(pv)
            except Exception:
                pass
        return float(periods[idx])

    def _filter_valid(self):
        mask = np.isfinite(self.time) & np.isfinite(self.mag)
        if self.mag_err is not None:
            mask &= np.isfinite(self.mag_err) & (self.mag_err > 0)
        dy = self.mag_err[mask] if self.mag_err is not None else None
        return self.time[mask], self.mag[mask], dy


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------

class VariableStarToolWindow(ToolWindowBase):
    """Variable Star Analysis Tool — period refinement, O-C, Fourier."""

    def __init__(self, params, project_state, parent=None):
        super().__init__(
            "Variable Star Analysis",
            params=params, project_state=project_state,
            parent=parent, min_size=(1000, 700),
        )
        self.lc_data: Optional[dict] = None
        self.series_options: dict[str, dict] = {}
        self.scan_result: Optional[dict] = None
        self.scan_alias_analysis: Optional[dict] = None
        self.handoff_multimode_diagnostic: Optional[dict] = None
        self.handoff_bundle: Optional[ValidatedLightCurveBundle] = None
        self.automated_result: Optional[VariableAnalysisResult] = None
        self._analysis_worker: Optional[VariableAnalysisWorker] = None
        self._review_kind = ""
        self.refined_period: Optional[float] = None
        self.sigma_period: Optional[float] = None
        self.multimode_result: Optional[dict] = None
        self.mm_candidate_scan: Optional[dict] = None
        self.mm_mode_rows: list[dict] = []
        self.mm_history: list[dict] = []
        self._scan_worker: Optional[PeriodAnalysisWorker] = None
        self._refine_worker: Optional[RefineBootstrapWorker] = None
        self._workspace_load_worker: Optional[_WorkspaceLoadWorker] = None
        self._workspace_load_token = 0
        self._pending_workspace_load: Optional[tuple[int, Path, list[Path]]] = None
        self.filter_colors: dict = {}      # user-customized per-filter colors
        self.filter_visibility: dict = {}  # True=visible, False=hidden
        self.workspace_dir = Path(self.params.P.result_dir)
        self.analysis_mode = "auto"
        self.recommended_mode = "unknown"
        self.recommendation_text = "Load a light curve and run a scan."
        self.workflow_step = "load"
        self.alias_override = False

        self.resize(1200, 800)
        self._build_ui()
        self._load_lc_from_workspace(synchronous=False)

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------

    def _build_ui(self):
        root = self.content_layout

        splitter = QSplitter(Qt.Horizontal)
        root.addWidget(splitter, 1)

        # ---- Left panel (controls) ----
        left_scroll = QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_scroll.setFrameShape(QFrame.NoFrame)
        left_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        # 280 squeezed the control rows below their labels — Run/Stop/Save
        # rendered as "Rur"/"Stoc"/"Savt" and the Workspace row's buttons to
        # 35 px of garbage. 330 is the widest that keeps the whole window
        # inside a 1280 px screen (the right pane needs ~948).
        left_scroll.setMinimumWidth(330)
        left_scroll.setMaximumWidth(380)
        left = QWidget()
        left.setMinimumWidth(308)   # 330 minus the vertical scrollbar
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(4, 4, 4, 4)
        left_scroll.setWidget(left)
        splitter.addWidget(left_scroll)

        # Validated, one-click analysis supplied by Main Step 12.
        self.auto_group = QGroupBox("Validated Analysis")
        self.auto_group.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        # No blanket stylesheet: the old fixed #263238 text turned dark-on-dark
        # in every dark theme; the theme QSS already colours all of these.
        auto_layout = QVBoxLayout(self.auto_group)
        auto_layout.setContentsMargins(8, 8, 8, 8)
        auto_layout.setSpacing(6)
        self.auto_release_label = QLabel(
            "No Main-workflow release is loaded. Standalone data is unvalidated."
        )
        self.auto_release_label.setWordWrap(True)
        self.auto_release_label.setMinimumHeight(78)
        self.auto_release_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self._set_release_banner("neutral")
        auto_layout.addWidget(self.auto_release_label)

        mode_row = QHBoxLayout()
        mode_row.addWidget(QLabel("Mode:"))
        self.auto_branch_combo = QComboBox()
        self.auto_branch_combo.addItem("Auto from Step 12", "auto")
        self.auto_branch_combo.addItem("Single mode", "single")
        self.auto_branch_combo.addItem("Multiple modes", "multi")
        mode_row.addWidget(self.auto_branch_combo, 1)
        auto_layout.addLayout(mode_row)

        self.auto_candidate_row = QWidget()
        candidate_layout = QHBoxLayout(self.auto_candidate_row)
        candidate_layout.setContentsMargins(0, 0, 0, 0)
        self.auto_candidate_label = QLabel("Candidate:")
        candidate_layout.addWidget(self.auto_candidate_label)
        self.auto_candidate_combo = QComboBox()
        candidate_layout.addWidget(self.auto_candidate_combo, 1)
        self.auto_candidate_row.hide()
        auto_layout.addWidget(self.auto_candidate_row)

        action_row = QHBoxLayout()
        # Tight column: no icons here — at the 330 px panel floor, icon +
        # theme padding squeezed the labels to "Rur"/"Stoc"/"Savi".
        action_row.setSpacing(4)
        self.btn_auto_run = QPushButton("Run")
        style_button(self.btn_auto_run, "primary", height=Tokens.H_BUTTON)
        self.btn_auto_run.setEnabled(False)
        self.btn_auto_run.clicked.connect(self._run_automated_analysis)
        action_row.addWidget(self.btn_auto_run)
        self.btn_auto_cancel = QPushButton("Stop")
        style_button(self.btn_auto_cancel, "danger", height=Tokens.H_BUTTON)
        self.btn_auto_cancel.setToolTip("Cancel the running analysis")
        self.btn_auto_cancel.setEnabled(False)
        self.btn_auto_cancel.clicked.connect(self._cancel_automated_analysis)
        action_row.addWidget(self.btn_auto_cancel)
        self.btn_auto_export = QPushButton("Save")
        style_button(self.btn_auto_export, height=Tokens.H_BUTTON)
        self.btn_auto_export.setToolTip("Save the complete analysis result")
        self.btn_auto_export.setEnabled(False)
        self.btn_auto_export.clicked.connect(self._export_automated_result)
        action_row.addWidget(self.btn_auto_export)
        auto_layout.addLayout(action_row)

        self.auto_progress = QProgressBar()
        self.auto_progress.setRange(0, 100)
        self.auto_progress.setValue(0)
        self.auto_progress.setTextVisible(True)
        auto_layout.addWidget(self.auto_progress)
        self.auto_status_label = QLabel("Load a validated Step 12 release to run.")
        self.auto_status_label.setWordWrap(True)
        self.auto_status_label.setProperty("role", "caption")
        auto_layout.addWidget(self.auto_status_label)
        self.auto_advanced_chk = QCheckBox("Advanced manual controls")
        self.auto_advanced_chk.setChecked(True)
        self.auto_advanced_chk.toggled.connect(self._refresh_tool_workflow_ui)
        auto_layout.addWidget(self.auto_advanced_chk)
        left_layout.addWidget(self.auto_group)

        # Light curve
        self.lc_group = QGroupBox("Light Curve")
        lc_form = QFormLayout(self.lc_group)
        self.lc_status = QLabel("Not loaded")
        self.lc_status.setWordWrap(True)
        lc_form.addRow("Status:", self.lc_status)
        ws_row = QWidget()
        # Path + Browse + Load on one line left each button ~56 px in this
        # ~230 px field, which rendered them as unreadable fragments
        # ("ult" / "ws" / "oac"). Path on its own line, buttons beneath, so
        # both keep their full labels without widening the column.
        ws_layout = QVBoxLayout(ws_row)
        ws_layout.setContentsMargins(0, 0, 0, 0)
        ws_layout.setSpacing(4)
        self.workspace_edit = QLineEdit(str(self.workspace_dir))
        btn_workspace = QPushButton("Browse…")
        btn_workspace.clicked.connect(self._browse_workspace)
        btn_reload = QPushButton("Load")
        btn_reload.clicked.connect(lambda: self._load_lc_from_workspace())
        # Full-width row (no form label column): the label column ate ~100 px
        # of a ~265 px group, leaving 165 px for a path plus two buttons.
        ws_label = QLabel("Workspace:")
        ws_layout.addWidget(ws_label)
        ws_layout.addWidget(self.workspace_edit)
        ws_btn_row = QHBoxLayout()
        ws_btn_row.setContentsMargins(0, 0, 0, 0)
        ws_btn_row.setSpacing(4)
        ws_btn_row.addStretch(1)
        ws_btn_row.addWidget(btn_workspace)
        ws_btn_row.addWidget(btn_reload)
        ws_layout.addLayout(ws_btn_row)
        lc_form.addRow(ws_row)
        self.mag_col_combo = QComboBox()
        self.mag_col_combo.setEnabled(False)
        self.mag_col_combo.currentIndexChanged.connect(self._on_mag_col_changed)
        lc_form.addRow("Use data:", self.mag_col_combo)
        self.analysis_filter_combo = QComboBox()
        self.analysis_filter_combo.setEnabled(False)
        self.analysis_filter_combo.currentIndexChanged.connect(self._on_analysis_filter_changed)
        lc_form.addRow("Filter:", self.analysis_filter_combo)
        btn_filt_browser = QPushButton("Colors / Visibility…")
        btn_filt_browser.clicked.connect(self.show_filter_color_browser)
        lc_form.addRow("Display:", btn_filt_browser)
        left_layout.addWidget(self.lc_group)

        # Workflow / routing
        self.workflow_group = QGroupBox("Analysis Workflow")
        workflow_layout = QVBoxLayout(self.workflow_group)
        workflow_layout.setContentsMargins(8, 8, 8, 8)
        workflow_layout.setSpacing(6)

        self.workflow_status_label = QLabel("Load a light curve to begin.")
        self.workflow_status_label.setWordWrap(True)
        self.workflow_status_label.setProperty("role", "info")
        workflow_layout.addWidget(self.workflow_status_label)

        workflow_btn_layout = QGridLayout()
        workflow_btn_layout.setHorizontalSpacing(6)
        workflow_btn_layout.setVerticalSpacing(6)
        self.workflow_buttons = {}
        for idx, (step_key, label) in enumerate((
            ("load", "1. Data"),
            ("scan", "2. Period Review"),
            ("analysis", "3. Analysis"),
        )):
            btn = QPushButton(label)
            btn.setCheckable(True)
            btn.clicked.connect(
                lambda _checked=False, stage=step_key: self._navigate_workflow_stage(stage)
            )
            # One row each: side by side in a ~300 px column these got 81 px
            # against the 177 px "2. Period Review" needs, so every label was
            # cut. Stacking keeps the stage names readable.
            workflow_btn_layout.addWidget(btn, idx, 0)
            self.workflow_buttons[step_key] = btn
        workflow_layout.addLayout(workflow_btn_layout)
        left_layout.addWidget(self.workflow_group)

        # Period scan
        self.scan_group = QGroupBox("Period Scan")
        scan_form = QFormLayout(self.scan_group)
        self.min_p = QDoubleSpinBox()
        self.min_p.setRange(0.001, 500); self.min_p.setDecimals(4)
        self.min_p.setValue(0.05); self.min_p.setSuffix(" d")
        scan_form.addRow("P min:", self.min_p)
        self.max_p = QDoubleSpinBox()
        self.max_p.setRange(0.01, 2000); self.max_p.setDecimals(4)
        self.max_p.setValue(100.0); self.max_p.setSuffix(" d")
        scan_form.addRow("P max:", self.max_p)
        self.spp = QSpinBox()
        self.spp.setRange(5, 50); self.spp.setValue(10)
        scan_form.addRow("Samples/peak:", self.spp)
        method_row = QHBoxLayout()
        self.chk_ls = QCheckBox("LS"); self.chk_ls.setChecked(True)
        self.chk_pdm = QCheckBox("PDM"); self.chk_pdm.setChecked(True)
        self.chk_bls = QCheckBox("BLS"); self.chk_bls.setChecked(False)
        method_row.addWidget(self.chk_ls)
        method_row.addWidget(self.chk_pdm)
        method_row.addWidget(self.chk_bls)
        scan_form.addRow("Methods:", method_row)
        self.pdm_bins = QSpinBox()
        self.pdm_bins.setRange(5, 50); self.pdm_bins.setValue(10)
        scan_form.addRow("PDM bins:", self.pdm_bins)
        btn_scan = QPushButton("Scan")
        style_button(btn_scan, "primary", height=Tokens.H_ACTION)
        btn_scan.clicked.connect(self._run_scan)
        scan_form.addRow(btn_scan)
        self.scan_status = QLabel("")
        self.scan_status.setProperty("role", "caption")
        scan_form.addRow(self.scan_status)
        left_layout.addWidget(self.scan_group)

        self.route_group = QGroupBox("Route")
        route_form = QFormLayout(self.route_group)
        self.route_recommend_label = QLabel("Run a scan to get a single/multi recommendation.")
        self.route_recommend_label.setWordWrap(True)
        self.route_recommend_label.setProperty("role", "info")
        route_form.addRow("Recommendation:", self.route_recommend_label)
        route_btn_row = QHBoxLayout()
        self.btn_route_single = QPushButton("Single-period")
        self.btn_route_single.clicked.connect(lambda: self._route_to_analysis_path("single"))
        route_btn_row.addWidget(self.btn_route_single)
        self.btn_route_multi = QPushButton("Multi-mode")
        self.btn_route_multi.clicked.connect(lambda: self._route_to_analysis_path("multi"))
        route_btn_row.addWidget(self.btn_route_multi)
        route_form.addRow("Open:", route_btn_row)
        left_layout.addWidget(self.route_group)

        # Refine
        self.refine_group = QGroupBox("Local Period Precision")
        refine_form = QFormLayout(self.refine_group)
        rp_row = QHBoxLayout()
        self.center_p = QDoubleSpinBox()
        self.center_p.setRange(0.00001, 10000); self.center_p.setDecimals(8)
        self.center_p.setValue(1.0); self.center_p.setSuffix(" d")
        rp_row.addWidget(self.center_p)
        btn_from_scan = QPushButton("← Best")
        btn_from_scan.setMaximumWidth(55)
        btn_from_scan.clicked.connect(self._set_center_from_scan)
        rp_row.addWidget(btn_from_scan)
        refine_form.addRow("Center P:", rp_row)
        self.n_boot = QSpinBox()
        self.n_boot.setRange(50, 2000); self.n_boot.setValue(300); self.n_boot.setSingleStep(50)
        self.n_boot.setToolTip(
            "Residual resamples estimate local precision around the selected alias; "
            "they do not resolve competing aliases."
        )
        refine_form.addRow("Resamples:", self.n_boot)
        self.refine_method_combo = QComboBox()
        self.refine_method_combo.addItem("Lomb-Scargle", "ls")
        self.refine_method_combo.addItem("PDM", "pdm")
        refine_form.addRow("Method:", self.refine_method_combo)
        self.btn_refine = QPushButton("Refine Selected Alias")
        self.btn_refine.setToolTip(
            "Runs a fine search and fixed-timestamp residual bootstrap around the selected period."
        )
        style_button(self.btn_refine, "primary", height=Tokens.H_ACTION)
        self.btn_refine.setEnabled(False)
        self.btn_refine.clicked.connect(self._run_refine)
        refine_form.addRow(self.btn_refine)
        self.refine_status = QLabel("")
        self.refine_status.setProperty("role", "caption")
        refine_form.addRow(self.refine_status)
        left_layout.addWidget(self.refine_group)

        # Multi-mode fit
        self.mm_group = QGroupBox("Multi-Mode Fit")
        mm_form = QFormLayout(self.mm_group)
        self.mm_periods_edit = QLineEdit()
        self.mm_periods_edit.setPlaceholderText("0.086017, 0.066529")
        mm_form.addRow("Periods (d):", self.mm_periods_edit)
        self.mm_n_modes = QSpinBox()
        self.mm_n_modes.setRange(1, 6); self.mm_n_modes.setValue(2)
        mm_form.addRow("Top peaks:", self.mm_n_modes)
        btn_mm_row = QHBoxLayout()
        btn_mm_best = QPushButton("+ Best")
        btn_mm_best.clicked.connect(self._append_current_period_to_multimode)
        btn_mm_row.addWidget(btn_mm_best)
        btn_mm_peaks = QPushButton("Top Peaks")
        btn_mm_peaks.clicked.connect(self._set_multimode_periods_from_scan)
        btn_mm_row.addWidget(btn_mm_peaks)
        btn_mm_remove = QPushButton("- Remove")
        btn_mm_remove.clicked.connect(self._remove_selected_mode_from_multimode)
        btn_mm_row.addWidget(btn_mm_remove)
        mm_form.addRow(btn_mm_row)
        self.mm_harm = QSpinBox()
        self.mm_harm.setRange(1, 4); self.mm_harm.setValue(2)
        mm_form.addRow("Harmonics/mode:", self.mm_harm)
        self.mm_alias_search_chk = QCheckBox("Compare sampling-window aliases")
        self.mm_alias_search_chk.setChecked(True)
        self.mm_alias_search_chk.setToolTip(
            "Disable for fixed literature-frequency validation."
        )
        mm_form.addRow(self.mm_alias_search_chk)
        self.mm_night_offsets_chk = QCheckBox("Fit nightly zero-point terms")
        self.mm_night_offsets_chk.setChecked(False)
        self.mm_night_offsets_chk.setToolTip(
            "Use only for uncorrected nightly baselines; free offsets can absorb long-period modes."
        )
        mm_form.addRow(self.mm_night_offsets_chk)
        focus_row = QHBoxLayout()
        self.mm_focus_combo = QComboBox()
        self.mm_focus_combo.setEnabled(False)
        self.mm_focus_combo.currentIndexChanged.connect(self._on_multimode_focus_changed)
        focus_row.addWidget(self.mm_focus_combo, 1)
        btn_mm_phase = QPushButton("→ Phase")
        btn_mm_phase.clicked.connect(self._apply_multimode_focus_to_phase)
        focus_row.addWidget(btn_mm_phase)
        mm_form.addRow("Send mode:", focus_row)
        state_row = QHBoxLayout()
        self.mm_state_epoch = QDoubleSpinBox()
        self.mm_state_epoch.setRange(2400000, 2600000); self.mm_state_epoch.setDecimals(6)
        self.mm_state_epoch.setValue(2458000.0)
        self.mm_state_epoch.valueChanged.connect(self._draw_multimode)
        state_row.addWidget(self.mm_state_epoch, 1)
        btn_state_t0 = QPushButton("← T₀")
        btn_state_t0.clicked.connect(lambda: self.mm_state_epoch.setValue(self.t0_edit.value()))
        state_row.addWidget(btn_state_t0)
        mm_form.addRow("State epoch:", state_row)
        self.mm_project_chk = QCheckBox("Show projected beat window")
        self.mm_project_chk.setChecked(True)
        self.mm_project_chk.toggled.connect(self._draw_multimode)
        mm_form.addRow(self.mm_project_chk)
        self.mm_overlay_chk = QCheckBox("Overlay isolated mode on Phase Plot")
        self.mm_overlay_chk.setChecked(False)
        self.mm_overlay_chk.toggled.connect(self._update_phase_plot)
        mm_form.addRow(self.mm_overlay_chk)
        self.mm_overlay_chk.hide()
        self.btn_multimode_fit = QPushButton("Fit Multi-Mode")
        style_button(self.btn_multimode_fit, "primary", height=Tokens.H_ACTION)
        self.btn_multimode_fit.clicked.connect(self._run_multimode_fit)
        mm_form.addRow(self.btn_multimode_fit)
        self.mm_status = QLabel("")
        self.mm_status.setProperty("role", "caption")
        self.mm_status.setWordWrap(True)
        mm_form.addRow(self.mm_status)
        left_layout.addWidget(self.mm_group)

        # Phase plot controls
        self.phase_group = QGroupBox("Phase Plot")
        phase_form = QFormLayout(self.phase_group)
        self.phase_p = QDoubleSpinBox()
        self.phase_p.setRange(0.00001, 10000); self.phase_p.setDecimals(8)
        self.phase_p.setValue(1.0); self.phase_p.setSuffix(" d")
        self.phase_p.valueChanged.connect(self._on_phase_period_changed)
        phase_form.addRow("Period:", self.phase_p)
        self.phase_mode_combo = QComboBox()
        self.phase_mode_combo.addItem("Manual / custom period", "manual")
        self.phase_mode_combo.setEnabled(False)
        self.phase_mode_combo.currentIndexChanged.connect(self._on_phase_mode_changed)
        phase_form.addRow("Mode from fit:", self.phase_mode_combo)
        self.t0_edit = QDoubleSpinBox()
        self.t0_edit.setRange(2400000, 2600000); self.t0_edit.setDecimals(6)
        self.t0_edit.setValue(2458000.0)
        self.t0_edit.valueChanged.connect(self._update_phase_plot)
        self.t0_edit.valueChanged.connect(self._draw_multimode)
        phase_form.addRow("T₀ (BJD):", self.t0_edit)
        self.epoch_kind_combo = QComboBox()
        self.epoch_kind_combo.addItem("Maximum light (pulsator)", "max_light")
        self.epoch_kind_combo.addItem("Minimum light (eclipse)", "min_light")
        self.epoch_kind_combo.currentIndexChanged.connect(self._on_epoch_kind_changed)
        phase_form.addRow("Epoch event:", self.epoch_kind_combo)
        self.phase_data_mode = QComboBox()
        self.phase_data_mode.addItem("Raw composite (diagnostic)", "raw")
        self.phase_data_mode.addItem("Focused mode isolation", "focused")
        self.phase_data_mode.currentIndexChanged.connect(self._update_phase_plot)
        phase_form.addRow("Data view:", self.phase_data_mode)
        self.phase_check_overlay_chk = QCheckBox("Show check star")
        self.phase_check_overlay_chk.setChecked(False)
        self.phase_check_overlay_chk.toggled.connect(self._update_phase_plot)
        phase_form.addRow(self.phase_check_overlay_chk)
        self.phase_fit_chk = QCheckBox("Overlay phase fit")
        self.phase_fit_chk.setChecked(True)
        self.phase_fit_chk.toggled.connect(self._update_phase_plot)
        phase_form.addRow(self.phase_fit_chk)
        self.phase_fit_harm = QSpinBox()
        self.phase_fit_harm.setRange(1, 8)
        self.phase_fit_harm.setValue(4)
        self.phase_fit_harm.valueChanged.connect(self._update_phase_plot)
        phase_form.addRow("Fit harmonics:", self.phase_fit_harm)
        self.btn_detect_t0 = QPushButton("Detect Maximum-light Epoch")
        self.btn_detect_t0.setToolTip(
            "For multi-mode data, isolate a fitted mode before measuring an event epoch."
        )
        self.btn_detect_t0.clicked.connect(self._detect_t0)
        phase_form.addRow(self.btn_detect_t0)
        left_layout.addWidget(self.phase_group)

        left_layout.addStretch()

        # Log (collapsible)
        self.btn_log_toggle = QPushButton("Log ▼")
        self.btn_log_toggle.setCheckable(True)
        self.btn_log_toggle.setChecked(False)
        style_button(self.btn_log_toggle, "ghost", height=Tokens.H_COMPACT)
        self.btn_log_toggle.toggled.connect(self._toggle_log)
        left_layout.addWidget(self.btn_log_toggle)
        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        self.log_box.setMaximumHeight(150)
        self.log_box.setObjectName("Log")
        self.log_box.hide()
        left_layout.addWidget(self.log_box)

        # ---- Right panel (tabs) ----
        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(4, 4, 4, 4)
        splitter.addWidget(right)
        splitter.setSizes([300, 900])
        prevent_collapse(splitter)

        self.tabs = QTabWidget()
        right_layout.addWidget(self.tabs)

        # Periodogram tab
        pg_tab = QWidget()
        pg_layout = QVBoxLayout(pg_tab)
        self.pg_canvas = FigureCanvas(Figure(figsize=(8, 4)))
        pg_layout.addWidget(NavigationToolbar(self.pg_canvas, pg_tab))
        pg_layout.addWidget(tame_canvas(self.pg_canvas), 1)
        self.tabs.addTab(pg_tab, "Periodogram")

        # Refine tab
        ref_tab = QWidget()
        ref_layout = QVBoxLayout(ref_tab)
        self.refine_label = QLabel("Select a resolved period candidate, then refine its local precision.")
        self.refine_label.setProperty("role", "info")
        ref_layout.addWidget(self.refine_label)
        self.refine_result_tabs = QTabWidget()
        ref_left = QWidget()
        ref_left_l = QVBoxLayout(ref_left)
        ref_left_l.setContentsMargins(0, 0, 0, 0)
        self.ref_canvas = FigureCanvas(Figure(figsize=(5, 4)))
        ref_left_l.addWidget(NavigationToolbar(self.ref_canvas, ref_left))
        ref_left_l.addWidget(tame_canvas(self.ref_canvas), 1)
        self.refine_result_tabs.addTab(ref_left, "Fine Search")
        ref_right = QWidget()
        ref_right_l = QVBoxLayout(ref_right)
        ref_right_l.setContentsMargins(0, 0, 0, 0)
        self.boot_canvas = FigureCanvas(Figure(figsize=(5, 4)))
        ref_right_l.addWidget(NavigationToolbar(self.boot_canvas, ref_right))
        ref_right_l.addWidget(tame_canvas(self.boot_canvas), 1)
        self.refine_result_tabs.addTab(ref_right, "Local Bootstrap")
        ref_layout.addWidget(self.refine_result_tabs, 1)
        self.tabs.addTab(ref_tab, "Refine")

        # Multi-mode tab
        mm_tab = self._build_multimode_tab()
        # Tallest page (1025 px): scroll it so it does not set the window's
        # minimum height — see layout_rules.scroll_wrap.
        self.tabs.addTab(scroll_wrap(mm_tab), "Multi-Mode")

        # Phase plot tab
        ph_tab = QWidget()
        ph_layout = QVBoxLayout(ph_tab)
        self.ph_canvas = FigureCanvas(Figure(figsize=(8, 5)))
        ph_layout.addWidget(NavigationToolbar(self.ph_canvas, ph_tab))
        ph_layout.addWidget(tame_canvas(self.ph_canvas), 1)
        self.tabs.addTab(ph_tab, "Phase Plot")

        # O-C tab
        oc_tab = self._build_oc_tab()
        self.tabs.addTab(oc_tab, "O-C")

        # Fourier tab
        fo_tab = self._build_fourier_tab()
        self.tabs.addTab(fo_tab, "Fourier")
        self.tabs.currentChanged.connect(self._sync_left_panel_visibility)
        self._refresh_tool_workflow_ui()

    def _effective_analysis_mode(self) -> str:
        if self.analysis_mode in ("single", "multi"):
            return self.analysis_mode
        if self.recommended_mode in ("single", "multi"):
            return self.recommended_mode
        return "single"

    def _set_tool_tab(self, title: str) -> None:
        if not hasattr(self, "tabs"):
            return
        for idx in range(self.tabs.count()):
            if self.tabs.tabText(idx) == title:
                self.tabs.setCurrentIndex(idx)
                return

    def _navigate_workflow_stage(self, stage: str) -> None:
        if stage == "analysis":
            self._set_workflow_step(self._effective_analysis_mode())
        else:
            self._set_workflow_step(stage)

    def _set_workflow_step(self, step: str) -> None:
        if step not in {"load", "scan", "single", "multi"}:
            return
        if step != "load" and self.lc_data is None:
            step = "load"
        if step in {"single", "multi"} and self.scan_result is None:
            step = "scan"
        self.workflow_step = step
        if step == "scan":
            self._set_tool_tab("Periodogram")
        elif step == "single":
            self._set_tool_tab("Refine")
        elif step == "multi":
            self._set_tool_tab("Multi-Mode")
        self._refresh_tool_workflow_ui()

    def _route_to_analysis_path(self, mode: str) -> None:
        if mode not in {"single", "multi"}:
            return
        self.analysis_mode = mode
        alias_status = str((self.scan_alias_analysis or {}).get("status", "UNASSESSED")).upper()
        if mode == "single" and alias_status != "RESOLVED":
            self.alias_override = True
            self.log(
                f"[ALIAS] Manual Single-period route selected while alias status is {alias_status}. "
                "Refinement uncertainty will remain conditional on this alias."
            )
        if mode == "multi" and not _parse_period_list(self.mm_periods_edit.text()) and self.scan_result:
            self._set_multimode_periods_from_scan()
        self._set_workflow_step(mode)
        if mode == "multi" and self.lc_data is not None:
            self._detect_multimode_candidates()

    def _refresh_tool_workflow_ui(self) -> None:
        has_data = self.lc_data is not None
        has_scan = self.scan_result is not None
        effective_mode = self._effective_analysis_mode()
        mode_text = {
            "auto": f"Auto → {effective_mode.title()}" if effective_mode else "Auto",
            "single": "Single",
            "multi": "Multi",
        }.get(self.analysis_mode, str(self.analysis_mode))
        step_text = {
            "load": "Load",
            "scan": "Scan",
            "single": "Single",
            "multi": "Multi",
        }.get(self.workflow_step, self.workflow_step)
        if hasattr(self, "workflow_status_label"):
            self.workflow_status_label.setText(
                f"Current step: {step_text}\nAnalysis path: {mode_text}\n{self.recommendation_text}"
            )

        enabled_map = {
            "load": True,
            "scan": has_data,
            "analysis": has_scan and self.analysis_mode in {"single", "multi"},
        }
        for key, btn in getattr(self, "workflow_buttons", {}).items():
            btn.blockSignals(True)
            btn.setEnabled(enabled_map.get(key, False))
            is_current = key == self.workflow_step or (
                key == "analysis" and self.workflow_step in {"single", "multi"}
            )
            btn.setChecked(is_current)
            # The checked state alone marks the current stage; the theme
            # paints it, so no per-state stylesheet (and dark mode works).
            btn.blockSignals(False)

        if hasattr(self, "route_recommend_label"):
            banner = {"single": "ok", "multi": "warn"}.get(self.recommended_mode)
            self.route_recommend_label.setProperty("role", None if banner else "info")
            self.route_recommend_label.setProperty("banner", banner)
            _repolish(self.route_recommend_label)
            self.route_recommend_label.setText(self.recommendation_text)

        if hasattr(self, "btn_route_single"):
            self.btn_route_single.setEnabled(has_scan)
        if hasattr(self, "btn_route_multi"):
            self.btn_route_multi.setEnabled(has_scan)
        if hasattr(self, "btn_refine"):
            alias_status = str((self.scan_alias_analysis or {}).get("status", "UNASSESSED")).upper()
            refine_running = bool(self._refine_worker and self._refine_worker.isRunning())
            self.btn_refine.setEnabled(
                has_scan
                and not refine_running
                and (alias_status == "RESOLVED" or self.alias_override)
            )

        self.scan_group.setVisible(self.workflow_step == "scan")
        self.route_group.setVisible(self.workflow_step == "scan")
        self.refine_group.setVisible(self.workflow_step == "single")
        self.mm_group.setVisible(self.workflow_step == "multi")
        self.phase_group.setVisible(self.workflow_step in {"single", "multi"})
        compact_handoff = bool(
            self.handoff_bundle is not None
            and hasattr(self, "auto_advanced_chk")
            and not self.auto_advanced_chk.isChecked()
        )
        if compact_handoff:
            for group in (
                self.lc_group,
                self.workflow_group,
                self.scan_group,
                self.route_group,
                self.refine_group,
                self.mm_group,
                self.phase_group,
            ):
                group.hide()
        if hasattr(self, "btn_log_toggle"):
            self.btn_log_toggle.setVisible(not compact_handoff)
            if compact_handoff:
                self.log_box.hide()
        self._refresh_automated_tab_visibility(compact_handoff)

    def _sync_left_panel_visibility(self, *_):
        self._refresh_tool_workflow_ui()

    def _refresh_automated_tab_visibility(self, compact_handoff: bool) -> None:
        if not hasattr(self, "tabs") or not hasattr(self.tabs, "setTabVisible"):
            return
        branch = (
            self.automated_result.branch
            if self.automated_result is not None and self.automated_result.branch
            else self._effective_analysis_mode()
        )
        visible_when_compact = {
            "Periodogram": False,
            "Refine": True,
            "Multi-Mode": branch == "multi",
            "Phase Plot": True,
            "O-C": True,
            "Fourier": branch != "multi",
        }
        for index in range(self.tabs.count()):
            title = self.tabs.tabText(index)
            self.tabs.setTabVisible(
                index,
                True if not compact_handoff else visible_when_compact.get(title, True),
            )

    def _configure_automated_handoff(self, bundle: ValidatedLightCurveBundle) -> None:
        self.handoff_bundle = bundle
        self.automated_result = None
        self._review_kind = ""
        self.auto_candidate_combo.clear()
        self.auto_candidate_row.hide()
        auto_index = self.auto_branch_combo.findData("auto")
        if auto_index >= 0:
            self.auto_branch_combo.setCurrentIndex(auto_index)
        self.auto_progress.setValue(0)
        self.btn_auto_run.setEnabled(bundle.can_launch)
        self.btn_auto_cancel.setEnabled(False)
        self.btn_auto_export.setEnabled(False)
        alias_status = str(bundle.alias_analysis.get("status", "UNASSESSED")).upper()
        mode_status = str(bundle.multimode_diagnostic.get("status", "UNASSESSED")).upper()
        qc_status = str(bundle.main_qc.get("status", bundle.release_status)).upper()
        filter_label = (
            "All filters" if bundle.analysis_filter == "__all__" else bundle.analysis_filter
        )
        self.auto_release_label.setText(
            f"Main QC: {qc_status}\n"
            f"Alias: {alias_status} | Mode: {mode_status}\n"
            f"ID {bundle.target_id} | {filter_label} | P={bundle.adopted_period:.8f} d"
        )
        if bundle.release_status in {"APPROVED", "OVERRIDDEN"}:
            self._set_release_banner("ok")
            self.auto_status_label.setText("Ready to run the validated advanced analysis.")
        else:
            self._set_release_banner("warn")
            self.auto_status_label.setText(
                "This legacy handoff is unverified. Review provenance before using its result."
            )
        self.auto_advanced_chk.setChecked(False)
        self._refresh_tool_workflow_ui()

    def _set_release_banner(self, kind: str) -> None:
        """Release-state banner via theme roles: neutral -> info chip,
        ok -> green banner, warn -> amber banner. Dark-theme safe, unlike the
        old fixed light-pastel stylesheets."""
        label = self.auto_release_label
        label.setProperty("role", "info" if kind == "neutral" else None)
        label.setProperty("banner", kind if kind in ("ok", "warn") else None)
        _repolish(label)

    def _selected_review_periods(self) -> tuple[float | None, float | None]:
        if not self.auto_candidate_row.isVisible():
            return None, None
        data = self.auto_candidate_combo.currentData()
        try:
            period = float(data)
        except (TypeError, ValueError):
            return None, None
        if not np.isfinite(period) or period <= 0:
            return None, None
        if self._review_kind == "alias":
            return period, None
        if self._review_kind == "second_mode":
            return None, period
        return None, None

    def _run_automated_analysis(self) -> None:
        if self.handoff_bundle is None:
            QMessageBox.information(
                self,
                "Validated Analysis",
                "Open this tool from Main Step 12 to run validated automation.",
            )
            return
        if self._analysis_worker is not None and self._analysis_worker.isRunning():
            return
        primary_override, secondary_override = self._selected_review_periods()
        branch = str(self.auto_branch_combo.currentData() or "auto")
        if self._review_kind == "second_mode":
            branch = "multi"
        request = VariableAnalysisRequest(
            bundle=self.handoff_bundle,
            adopted_period_override=primary_override,
            secondary_period_override=secondary_override,
            analysis_branch=branch,
            bootstrap_resamples=int(self.n_boot.value()),
            refinement_harmonics=int(self.phase_fit_harm.value()),
            single_harmonics=int(self.phase_fit_harm.value()),
            multimode_harmonics=int(self.mm_harm.value()),
            include_night_offsets=False,
        )
        self.automated_result = None
        self.btn_auto_run.setEnabled(False)
        self.btn_auto_cancel.setEnabled(True)
        self.btn_auto_export.setEnabled(False)
        self.auto_progress.setValue(0)
        self.auto_status_label.setText("Starting validated advanced analysis...")
        self._analysis_worker = VariableAnalysisWorker(request, parent=self)
        self._analysis_worker.progress.connect(self._on_automated_progress)
        self._analysis_worker.finished.connect(self._on_automated_finished)
        self._analysis_worker.failed.connect(self._on_automated_failed)
        self._analysis_worker.start()

    def _cancel_automated_analysis(self) -> None:
        if self._analysis_worker is not None and self._analysis_worker.isRunning():
            self._analysis_worker.stop()
            self.auto_status_label.setText("Cancelling after the current calculation...")
            self.btn_auto_cancel.setEnabled(False)

    def _on_automated_progress(self, stage: str, percent: int) -> None:
        self.auto_progress.setValue(int(percent))
        self.auto_status_label.setText(str(stage))

    def _on_automated_finished(self, result: VariableAnalysisResult) -> None:
        self.btn_auto_cancel.setEnabled(False)
        self.btn_auto_run.setEnabled(self.handoff_bundle is not None)
        worker = self._analysis_worker
        self._analysis_worker = None
        if worker is not None:
            worker.deleteLater()
        if result.status == "REVIEW_REQUIRED":
            self._show_automated_review(result)
            return
        if result.status == "CANCELLED":
            self.auto_status_label.setText("Analysis cancelled.")
            return
        if result.status != "COMPLETE":
            self.auto_status_label.setText(f"Analysis failed: {result.error}")
            self.auto_progress.setValue(0)
            return
        self.automated_result = result
        self._review_kind = ""
        self.auto_candidate_row.hide()
        self.btn_auto_export.setEnabled(True)
        self.auto_progress.setValue(100)
        sigma = result.local_period_sigma
        sigma_text = f"{sigma:.2e} d" if sigma is not None and np.isfinite(sigma) else "not estimated"
        observed_cycles = float(result.data_summary.get("observed_cycles", np.nan))
        coverage_text = (
            f" | coverage={observed_cycles:.2f} cycles"
            if np.isfinite(observed_cycles)
            else ""
        )
        warning_text = (
            " | short-baseline local result"
            if result.diagnostics.get("limited_cycle_coverage")
            else ""
        )
        self.auto_status_label.setText(
            f"Complete | {result.branch} | P={result.refined_period:.8f} d | "
            f"sigma_local={sigma_text}{coverage_text}{warning_text}"
        )
        self._install_automated_result(result)
        self.log(
            f"[AUTO] Completed {result.branch} analysis at P={result.refined_period:.8f} d; "
            f"Main QC was reused, not recomputed."
        )

    def _on_automated_failed(self, message: str) -> None:
        self.btn_auto_cancel.setEnabled(False)
        self.btn_auto_run.setEnabled(self.handoff_bundle is not None)
        worker = self._analysis_worker
        self._analysis_worker = None
        if worker is not None:
            worker.deleteLater()
        self.auto_status_label.setText(f"Worker failed: {message}")

    def _show_automated_review(self, result: VariableAnalysisResult) -> None:
        review = result.review
        if review is None:
            self.auto_status_label.setText("Observer review is required.")
            return
        self.auto_candidate_combo.clear()
        self._review_kind = (
            "alias" if review.code == "ALIAS_SELECTION_REQUIRED" else "second_mode"
        )
        self.auto_candidate_label.setText(
            "Primary alias:" if self._review_kind == "alias" else "Second mode:"
        )
        for index, candidate in enumerate(review.candidates):
            raw_period = candidate.get("period", candidate.get("candidate_period"))
            try:
                period = float(raw_period)
            except (TypeError, ValueError):
                continue
            if not np.isfinite(period) or period <= 0:
                continue
            relation = str(
                candidate.get("relation_to_best", candidate.get("relation", "candidate"))
            )
            self.auto_candidate_combo.addItem(
                f"{index + 1}. {period:.8f} d | {relation}",
                period,
            )
        self.auto_candidate_row.setVisible(self.auto_candidate_combo.count() > 0)
        self.auto_status_label.setText(review.message)
        self.auto_progress.setValue(0)
        self.btn_auto_run.setEnabled(self.auto_candidate_combo.count() > 0)
        if self._review_kind == "second_mode":
            index = self.auto_branch_combo.findData("multi")
            if index >= 0:
                self.auto_branch_combo.setCurrentIndex(index)

    def _install_automated_result(self, result: VariableAnalysisResult) -> None:
        self.refined_period = float(result.refined_period)
        self.sigma_period = result.local_period_sigma
        self.phase_p.setValue(self.refined_period)
        self.oc_p.setValue(self.refined_period)
        refinement = dict(result.refinement)
        frequency = np.asarray(refinement.get("fine_frequency", []), dtype=float)
        power = np.asarray(refinement.get("fine_power", []), dtype=float)
        plot_result = {
            "fine_periods": 1.0 / frequency[::-1] if len(frequency) else np.array([]),
            "fine_power": power[::-1],
            "refined_period": self.refined_period,
            "sigma_p": result.local_period_sigma,
            "boot_periods": np.asarray(refinement.get("bootstrap_periods", []), dtype=float),
        }
        self._draw_refine_plots(plot_result)
        if result.branch == "multi":
            successful = [
                (name, fit)
                for name, fit in result.per_filter_models.items()
                if "error" not in fit
            ]
            if successful:
                reference_filter, reference_fit = max(
                    successful,
                    key=lambda item: int(item[1].get("n_points", 0)),
                )
                self.multimode_result = dict(reference_fit)
                self.multimode_result["window_peaks"] = list(
                    (self.handoff_bundle.alias_analysis if self.handoff_bundle else {}).get(
                        "window_peaks", []
                    )
                )
                self.mm_periods_edit.setText(_format_period_list(self.multimode_result["periods"]))
                self.mm_focus_combo.blockSignals(True)
                self.mm_focus_combo.clear()
                for index, period in enumerate(self.multimode_result["periods"]):
                    self.mm_focus_combo.addItem(f"M{index + 1} | P={float(period):.8f} d", index)
                self.mm_focus_combo.setCurrentIndex(0)
                self.mm_focus_combo.setEnabled(True)
                self.mm_focus_combo.blockSignals(False)
                self._sync_phase_mode_combo()
                self._update_multimode_mode_views()
                self.mm_status.setText(
                    f"Automated per-filter fit complete; showing reference filter {reference_filter}."
                )
                self._draw_multimode()
            self.analysis_mode = "multi"
            self.workflow_step = "multi"
            self._set_tool_tab("Multi-Mode")
        else:
            self.multimode_result = None
            self.analysis_mode = "single"
            self.workflow_step = "single"
            self._set_tool_tab("Refine")
        self._update_phase_plot()
        self._refresh_tool_workflow_ui()

    def _export_automated_result(self) -> None:
        if self.automated_result is None or self.handoff_bundle is None:
            return
        try:
            output = (
                Path(self.handoff_bundle.workspace_dir)
                / "lc_period"
                / "variable_analysis_result.json"
            )
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(
                json.dumps(
                    self.automated_result.to_dict(json_safe=True),
                    indent=2,
                    ensure_ascii=False,
                    allow_nan=False,
                ),
                encoding="utf-8",
            )
            self.auto_status_label.setText(f"Exported: {output}")
            self.log(f"[AUTO] Exported result to {output}")
        except OSError as exc:
            QMessageBox.warning(self, "Export", str(exc))

    def _build_oc_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)

        hdr = QHBoxLayout()
        hdr.addWidget(QLabel("Event:"))
        self.oc_event_label = QLabel("Maximum light")
        _f = self.oc_event_label.font(); _f.setBold(True)
        self.oc_event_label.setFont(_f)
        hdr.addWidget(self.oc_event_label)
        hdr.addWidget(QLabel("T₀ (BJD):"))
        self.oc_t0 = QDoubleSpinBox()
        self.oc_t0.setRange(2400000, 2600000); self.oc_t0.setDecimals(6)
        self.oc_t0.setValue(2458000.0); self.oc_t0.setMinimumWidth(130)
        self.oc_t0.valueChanged.connect(self._recompute_oc)
        hdr.addWidget(self.oc_t0)
        hdr.addWidget(QLabel("P (d):"))
        self.oc_p = QDoubleSpinBox()
        self.oc_p.setRange(0.00001, 10000); self.oc_p.setDecimals(8)
        self.oc_p.setValue(1.0); self.oc_p.setMinimumWidth(120)
        self.oc_p.valueChanged.connect(self._recompute_oc)
        hdr.addWidget(self.oc_p)
        btn_from_refine = QPushButton("← Refined")
        btn_from_refine.clicked.connect(self._oc_from_refine)
        hdr.addWidget(btn_from_refine)
        hdr.addStretch()
        layout.addLayout(hdr)
        oc_note = QLabel(
            "Cycle counts are inferred from T₀ and P. Verify epoch n before fitting, "
            "especially across long gaps."
        )
        oc_note.setWordWrap(True)
        oc_note.setProperty("banner", "warn")
        layout.addWidget(oc_note)

        splitter = QSplitter(Qt.Horizontal)

        left = QWidget()
        ll = QVBoxLayout(left)
        ll.setContentsMargins(0, 0, 0, 0)
        self.oc_table = QTableWidget()
        self.oc_table.setColumnCount(4)
        self.oc_table.setHorizontalHeaderLabels(["n", "BJD_obs", "O-C (d)", "err (d)"])
        self.oc_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.oc_table.horizontalHeader().setStretchLastSection(True)
        self.oc_table.setMinimumWidth(300)
        ll.addWidget(self.oc_table, 1)
        btn_row = QHBoxLayout()
        for label, slot in [("Add", self._oc_add), ("Del", self._oc_del),
                             ("Import CSV", self._oc_import), ("Export CSV", self._oc_export)]:
            b = QPushButton(label)
            b.clicked.connect(slot)
            btn_row.addWidget(b)
        btn_row.addStretch()
        ll.addLayout(btn_row)
        splitter.addWidget(left)

        right = QWidget()
        rl = QVBoxLayout(right)
        rl.setContentsMargins(0, 0, 0, 0)
        self.oc_canvas = FigureCanvas(Figure(figsize=(6, 4)))
        rl.addWidget(NavigationToolbar(self.oc_canvas, right))
        rl.addWidget(tame_canvas(self.oc_canvas), 1)
        fit_row = QHBoxLayout()
        fit_row.addWidget(QLabel("Fit:"))
        self.oc_fit_combo = QComboBox()
        self.oc_fit_combo.addItem("None", "none")
        self.oc_fit_combo.addItem("Linear (ΔP)", "linear")
        self.oc_fit_combo.addItem("Quadratic (dP/dt)", "quadratic")
        self.oc_fit_combo.addItem("Quadratic + sinusoid (exploratory)", "quadratic_sine")
        fit_row.addWidget(self.oc_fit_combo)
        btn_fit = QPushButton("Fit & Plot")
        btn_fit.clicked.connect(self._oc_fit)
        fit_row.addWidget(btn_fit)
        fit_row.addStretch()
        rl.addLayout(fit_row)
        self.oc_fit_label = QLabel("")
        self.oc_fit_label.setWordWrap(True)
        self.oc_fit_label.setStyleSheet(_mono_note_style())
        rl.addWidget(self.oc_fit_label)
        splitter.addWidget(right)
        splitter.setSizes([330, 670])
        prevent_collapse(splitter)
        layout.addWidget(splitter, 1)

        return tab

    def _build_fourier_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)

        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel("Harmonics:"))
        self.n_harm = QSpinBox()
        self.n_harm.setRange(1, 8); self.n_harm.setValue(4)
        ctrl.addWidget(self.n_harm)
        ctrl.addWidget(QLabel("Filter:"))
        self.fourier_filter_combo = QComboBox()
        self.fourier_filter_combo.setMinimumWidth(80)
        self.fourier_filter_combo.addItem("All data", "__all__")
        ctrl.addWidget(self.fourier_filter_combo)
        btn_fourier = QPushButton("Decompose")
        style_button(btn_fourier, "primary", height=Tokens.H_BUTTON)
        btn_fourier.clicked.connect(self._run_fourier)
        ctrl.addWidget(btn_fourier)
        ctrl.addStretch()
        layout.addLayout(ctrl)

        help_text = (
            "<b>Cosine-series convention</b>: "
            "m(τ) = A₀ + Σ Aₖ cos(2πkτ/P + φₖ), τ = t − t_ref. "
            "R₂₁ and φ₂₁ describe light-curve shape; they are not a standalone variable-type classifier. "
            "Multi-filter data are decomposed one passband at a time."
        )
        help_label = QLabel(help_text)
        help_label.setWordWrap(True)
        help_label.setProperty("role", "info")
        layout.addWidget(help_label)

        self.fourier_label = QLabel("")
        self.fourier_label.setStyleSheet(_mono_note_style())
        self.fourier_label.setWordWrap(True)
        layout.addWidget(self.fourier_label)

        self.fourier_canvas = FigureCanvas(Figure(figsize=(8, 4)))
        layout.addWidget(NavigationToolbar(self.fourier_canvas, tab))
        layout.addWidget(tame_canvas(self.fourier_canvas), 1)

        return tab

    def _build_multimode_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)

        help_label = QLabel(
            "<b>Multi-Mode Fit</b> — 여러 주기를 동시에 넣는 선형 다중주파수 Fourier fit입니다. "
            "문헌 period를 직접 넣거나 스캔 peak를 가져온 뒤, JD time-series 총합 fit과 "
            "각 모드별 분리 phase curve를 확인합니다."
        )
        help_label.setWordWrap(True)
        help_label.setProperty("role", "info")
        layout.addWidget(help_label)

        self.mm_summary_label = QLabel("Set periods and run Fit Multi-Mode.")
        self.mm_summary_label.setWordWrap(True)
        self.mm_summary_label.setStyleSheet(_mono_note_style())
        layout.addWidget(self.mm_summary_label)

        candidate_ctrl = QHBoxLayout()
        self.btn_mm_detect_candidates = QPushButton("Detect Residual Peaks")
        self.btn_mm_detect_candidates.clicked.connect(self._detect_multimode_candidates)
        candidate_ctrl.addWidget(self.btn_mm_detect_candidates)
        self.btn_mm_adopt_candidate = QPushButton("+ Selected Candidate")
        self.btn_mm_adopt_candidate.clicked.connect(self._append_selected_candidate_to_multimode)
        candidate_ctrl.addWidget(self.btn_mm_adopt_candidate)
        self.btn_mm_adopt_all = QPushButton("Adopt New Peaks")
        self.btn_mm_adopt_all.clicked.connect(self._adopt_all_new_candidates)
        candidate_ctrl.addWidget(self.btn_mm_adopt_all)
        candidate_ctrl.addStretch()
        layout.addLayout(candidate_ctrl)

        self.mm_candidate_status = QLabel("Candidate detection not run yet.")
        self.mm_candidate_status.setWordWrap(True)
        self.mm_candidate_status.setProperty("role", "caption")
        layout.addWidget(self.mm_candidate_status)

        self.mm_candidate_table = QTableWidget()
        self.mm_candidate_table.setColumnCount(5)
        self.mm_candidate_table.setHorizontalHeaderLabels(["#", "Period (d)", "Power", "Relation", "Note"])
        self.mm_candidate_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.mm_candidate_table.horizontalHeader().setStretchLastSection(True)
        self.mm_candidate_table.setMaximumHeight(160)
        self.mm_candidate_table.itemDoubleClicked.connect(lambda *_args: self._append_selected_candidate_to_multimode())
        layout.addWidget(self.mm_candidate_table)

        self.mm_mode_status = QLabel("No joint-fit mode classification yet.")
        self.mm_mode_status.setWordWrap(True)
        self.mm_mode_status.setProperty("role", "caption")
        layout.addWidget(self.mm_mode_status)

        self.mm_mode_table = QTableWidget()
        self.mm_mode_table.setColumnCount(8)
        self.mm_mode_table.setHorizontalHeaderLabels(
            ["Mode", "Period (d)", "Freq (d^-1)", "A1", "A2", "Δm", "Class", "Note"]
        )
        self.mm_mode_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.mm_mode_table.horizontalHeader().setStretchLastSection(True)
        self.mm_mode_table.setMaximumHeight(170)
        layout.addWidget(self.mm_mode_table)

        history_ctrl = QHBoxLayout()
        self.btn_mm_restore_history = QPushButton("Restore Selected Step")
        self.btn_mm_restore_history.clicked.connect(self._restore_selected_history_periods)
        history_ctrl.addWidget(self.btn_mm_restore_history)
        self.btn_mm_undo = QPushButton("Undo Last Change")
        self.btn_mm_undo.clicked.connect(self._rollback_last_multimode_step)
        history_ctrl.addWidget(self.btn_mm_undo)
        self.btn_mm_export_report = QPushButton("Export Report")
        self.btn_mm_export_report.clicked.connect(self._export_multimode_report)
        history_ctrl.addWidget(self.btn_mm_export_report)
        history_ctrl.addStretch()
        layout.addLayout(history_ctrl)

        self.mm_history_status = QLabel("Prewhitening history is empty.")
        self.mm_history_status.setWordWrap(True)
        self.mm_history_status.setProperty("role", "caption")
        layout.addWidget(self.mm_history_status)

        self.mm_history_table = QTableWidget()
        self.mm_history_table.setColumnCount(6)
        self.mm_history_table.setHorizontalHeaderLabels(
            ["#", "Action", "Stage", "Modes", "Best / Added", "Note"]
        )
        self.mm_history_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.mm_history_table.horizontalHeader().setStretchLastSection(True)
        self.mm_history_table.setMaximumHeight(180)
        self.mm_history_table.itemDoubleClicked.connect(lambda *_args: self._restore_selected_history_periods())
        layout.addWidget(self.mm_history_table)

        self.mm_pw_canvas = FigureCanvas(Figure(figsize=(8, 3.2)))
        layout.addWidget(NavigationToolbar(self.mm_pw_canvas, tab))
        layout.addWidget(tame_canvas(self.mm_pw_canvas, min_h=150), 1)

        self.mm_canvas = FigureCanvas(Figure(figsize=(8, 6)))
        layout.addWidget(NavigationToolbar(self.mm_canvas, tab))
        layout.addWidget(tame_canvas(self.mm_canvas, min_h=300), 2)

        return tab

    # ------------------------------------------------------------------
    # Multi-mode
    # ------------------------------------------------------------------

    def _clear_multimode_result(self, clear_inputs: bool = False):
        self.multimode_result = None
        self.mm_candidate_scan = None
        self.mm_mode_rows = []
        self.mm_history = []
        if clear_inputs and hasattr(self, "mm_periods_edit"):
            self.mm_periods_edit.clear()
        if hasattr(self, "mm_status"):
            self.mm_status.setText("")
        if hasattr(self, "mm_focus_combo"):
            self.mm_focus_combo.blockSignals(True)
            self.mm_focus_combo.clear()
            self.mm_focus_combo.setEnabled(False)
            self.mm_focus_combo.blockSignals(False)
        if hasattr(self, "mm_summary_label"):
            self.mm_summary_label.setText("Set periods and run Fit Multi-Mode.")
        if hasattr(self, "phase_mode_combo"):
            self.phase_mode_combo.blockSignals(True)
            self.phase_mode_combo.clear()
            self.phase_mode_combo.addItem("Manual / custom period", "manual")
            self.phase_mode_combo.setEnabled(False)
            self.phase_mode_combo.blockSignals(False)
        if hasattr(self, "mm_candidate_status"):
            self.mm_candidate_status.setText("Candidate detection not run yet.")
        if hasattr(self, "mm_candidate_table"):
            self.mm_candidate_table.setRowCount(0)
        if hasattr(self, "mm_mode_status"):
            self.mm_mode_status.setText("No joint-fit mode classification yet.")
        if hasattr(self, "mm_mode_table"):
            self.mm_mode_table.setRowCount(0)
        if hasattr(self, "mm_history_status"):
            self.mm_history_status.setText("Prewhitening history is empty.")
        if hasattr(self, "mm_history_table"):
            self.mm_history_table.setRowCount(0)
        if hasattr(self, "mm_pw_canvas"):
            fig = self.mm_pw_canvas.figure
            fig.clear()
            self.mm_pw_canvas.draw_idle()
        if hasattr(self, "mm_canvas"):
            fig = self.mm_canvas.figure
            fig.clear()
            self.mm_canvas.draw_idle()

    def _invalidate_multimode_runtime(self, status: str = "") -> None:
        self.multimode_result = None
        self.mm_candidate_scan = None
        self.mm_mode_rows = []
        if hasattr(self, "mm_focus_combo"):
            self.mm_focus_combo.blockSignals(True)
            self.mm_focus_combo.clear()
            self.mm_focus_combo.setEnabled(False)
            self.mm_focus_combo.blockSignals(False)
        if hasattr(self, "mm_summary_label"):
            self.mm_summary_label.setText("Period list changed. Re-run Detect/Fit for an updated multi-mode result.")
        if hasattr(self, "phase_mode_combo"):
            self.phase_mode_combo.blockSignals(True)
            self.phase_mode_combo.clear()
            self.phase_mode_combo.addItem("Manual / custom period", "manual")
            self.phase_mode_combo.setEnabled(False)
            self.phase_mode_combo.blockSignals(False)
        if hasattr(self, "mm_candidate_status"):
            self.mm_candidate_status.setText("Candidate detection needs to be re-run for the current period list.")
        if hasattr(self, "mm_candidate_table"):
            self.mm_candidate_table.setRowCount(0)
        if hasattr(self, "mm_mode_status"):
            self.mm_mode_status.setText("No joint-fit mode classification for the current period list yet.")
        if hasattr(self, "mm_mode_table"):
            self.mm_mode_table.setRowCount(0)
        for canvas_name in ("mm_pw_canvas", "mm_canvas"):
            canvas = getattr(self, canvas_name, None)
            if canvas is None:
                continue
            fig = getattr(canvas, "figure", None)
            if fig is None:
                continue
            fig.clear()
            canvas.draw_idle()
        if hasattr(self, "mm_status"):
            self.mm_status.setText(status or "Period list changed. Re-run Detect/Fit.")
        self._update_phase_plot()

    def _selected_history_entry(self) -> tuple[int, dict] | None:
        if not self.mm_history:
            return None
        rows = sorted({idx.row() for idx in self.mm_history_table.selectedIndexes()}) if hasattr(self, "mm_history_table") else []
        if not rows:
            return None
        row = rows[0]
        if 0 <= row < len(self.mm_history):
            return row, self.mm_history[row]
        return None

    def _restore_selected_history_periods(self) -> None:
        selected = self._selected_history_entry()
        if selected is None:
            QMessageBox.information(self, "Multi-Mode", "복원할 history step을 먼저 선택하세요.")
            return
        row_idx, entry = selected
        periods = _parse_period_list(str(entry.get("periods", "")))
        self.mm_periods_edit.setText(_format_period_list(periods))
        self._invalidate_multimode_runtime(
            f"Restored period list from history step {int(entry.get('step', row_idx + 1))}. Re-run Detect/Fit."
        )
        self._record_multimode_history(
            "restore_history",
            stage="history",
            periods=periods,
            best_period=periods[0] if periods else np.nan,
            note=f"restored period list from step {int(entry.get('step', row_idx + 1))}",
        )
        self.log(
            f"[MULTI] Restored history step {int(entry.get('step', row_idx + 1))}: "
            f"{_format_period_list(periods) or '<empty>'}"
        )

    def _rollback_last_multimode_step(self) -> None:
        current_periods = _parse_period_list(self.mm_periods_edit.text())
        current_text = _format_period_list(current_periods)
        for entry in reversed(self.mm_history[:-1]):
            prev_periods = _parse_period_list(str(entry.get("periods", "")))
            if _format_period_list(prev_periods) == current_text:
                continue
            self.mm_periods_edit.setText(_format_period_list(prev_periods))
            self._invalidate_multimode_runtime("Rolled back to the previous period list. Re-run Detect/Fit.")
            self._record_multimode_history(
                "rollback",
                stage="history",
                periods=prev_periods,
                best_period=prev_periods[0] if prev_periods else np.nan,
                note=f"rolled back using step {int(entry.get('step', 0))}",
            )
            self.log(
                f"[MULTI] Rolled back using step {int(entry.get('step', 0))}: "
                f"{_format_period_list(prev_periods) or '<empty>'}"
            )
            return
        if current_periods:
            self.mm_periods_edit.clear()
            self._invalidate_multimode_runtime("Rolled back to an empty period list.")
            self._record_multimode_history(
                "rollback",
                stage="history",
                periods=[],
                note="rolled back to empty period list",
            )
            self.log("[MULTI] Rolled back to empty period list")
            return
        QMessageBox.information(self, "Multi-Mode", "되돌릴 이전 period list가 없습니다.")

    def _selected_mode_table_index(self) -> int | None:
        if not hasattr(self, "mm_mode_table"):
            return None
        rows = sorted({idx.row() for idx in self.mm_mode_table.selectedIndexes()})
        if not rows:
            return None
        return rows[0]

    def _remove_selected_mode_from_multimode(self) -> None:
        periods = _parse_period_list(self.mm_periods_edit.text())
        if not periods:
            QMessageBox.information(self, "Multi-Mode", "제거할 period가 없습니다.")
            return
        remove_idx = self._selected_mode_table_index()
        if remove_idx is None and self.multimode_result:
            remove_idx = self._selected_multimode_focus_index()
        if remove_idx is None or not (0 <= int(remove_idx) < len(periods)):
            remove_idx = len(periods) - 1
        removed_period = float(periods.pop(int(remove_idx)))
        self.mm_periods_edit.setText(_format_period_list(periods))
        self._invalidate_multimode_runtime(
            f"Removed period {removed_period:.8f} d from the list. Re-run Detect/Fit."
        )
        self._record_multimode_history(
            "remove_mode",
            stage="manual",
            periods=periods,
            best_period=removed_period,
            note=f"removed period at index {int(remove_idx) + 1}",
        )
        self.log(f"[MULTI] Removed period {removed_period:.8f} d from the multi-mode list")

    def _record_multimode_history(
        self,
        action: str,
        *,
        stage: str = "",
        periods: list[float] | None = None,
        best_period: float = np.nan,
        best_power: float = np.nan,
        note: str = "",
        fit_result: dict | None = None,
    ) -> None:
        period_list = [float(p) for p in (periods or []) if np.isfinite(p) and float(p) > 0]
        entry = {
            "step": int(len(self.mm_history) + 1),
            "action": str(action),
            "stage": str(stage),
            "n_modes": int(len(period_list)),
            "periods": _format_period_list(period_list),
            "best_or_added": (
                f"{float(best_period):.8f}" if np.isfinite(best_period) and float(best_period) > 0
                else ""
            ),
            "best_power": (
                float(best_power) if np.isfinite(best_power) else np.nan
            ),
            "fit_rmse_mag": np.nan,
            "fit_wrms": np.nan,
            "note": str(note or ""),
        }
        if fit_result:
            entry["fit_rmse_mag"] = float(fit_result.get("rmse", np.nan))
            entry["fit_wrms"] = float(fit_result.get("wrms", np.nan))
            if not entry["note"]:
                entry["note"] = (
                    f"harm={int(fit_result.get('harmonics', 0))}, "
                    f"cond≈{float(fit_result.get('design_condition', np.nan)):.2e}"
                )
        self.mm_history.append(entry)
        self._update_multimode_history_views()

    def _update_multimode_history_views(self) -> None:
        if not hasattr(self, "mm_history_table"):
            return
        self.mm_history_table.setRowCount(0)
        for row_idx, entry in enumerate(self.mm_history):
            self.mm_history_table.insertRow(row_idx)
            values = [
                str(entry.get("step", row_idx + 1)),
                str(entry.get("action", "")),
                str(entry.get("stage", "")),
                str(entry.get("n_modes", "")),
                str(entry.get("best_or_added", "")),
                str(entry.get("note", "")),
            ]
            for col_idx, value in enumerate(values):
                self.mm_history_table.setItem(row_idx, col_idx, QTableWidgetItem(value))
        if hasattr(self, "mm_history_status"):
            if not self.mm_history:
                self.mm_history_status.setText("Prewhitening history is empty.")
            else:
                latest = self.mm_history[-1]
                self.mm_history_status.setText(
                    f"{len(self.mm_history)} action(s) recorded. Latest: "
                    f"{latest.get('action', '')} | modes={latest.get('n_modes', 0)} | "
                    f"{latest.get('note', '')}"
                )

    def _summarize_multimode_modes(self, result: dict) -> list[dict]:
        periods = [float(p) for p in result.get("periods", [])]
        coeff = np.asarray(result.get("coeff", []), dtype=float)
        terms = result.get("terms", [])
        mode_labels = _classify_fitted_periods(
            periods,
            window_peaks=list(result.get("window_peaks", [])),
            baseline_days=float(result.get("baseline", 1.0)),
        )
        harmonic_map: dict[tuple[int, int], dict[str, float]] = {}

        for coeff_idx, term in enumerate(terms, start=1):
            key = (int(term["mode_index"]), int(term["harmonic"]))
            entry = harmonic_map.setdefault(key, {})
            actual_idx = int(term.get("coefficient_index", coeff_idx))
            entry[str(term["kind"])] = float(coeff[actual_idx])

        rows: list[dict] = []
        harmonics = max(1, int(result.get("harmonics", 1)))
        time_ref = float(result.get("time_ref", 0.0))
        for idx, period in enumerate(periods):
            freq = _period_to_frequency(period)
            amp_h1 = np.nan
            amp_h2 = np.nan
            for harmonic in range(1, harmonics + 1):
                pair = harmonic_map.get((idx, harmonic), {})
                amp = float(np.hypot(float(pair.get("cos", 0.0)), float(pair.get("sin", 0.0))))
                if harmonic == 1:
                    amp_h1 = amp
                elif harmonic == 2:
                    amp_h2 = amp
            phase_grid = np.linspace(0.0, 1.0, 600)
            phase_times = time_ref + phase_grid * period
            eval_phase = _evaluate_multimode_result(result, phase_times)
            component_curve = np.asarray(eval_phase["components"][idx], dtype=float)
            ptp = float(np.nanmax(component_curve) - np.nanmin(component_curve)) if component_curve.size else np.nan
            label = mode_labels[idx] if idx < len(mode_labels) else {"relation": "independent", "note": ""}
            rows.append({
                "mode": f"M{idx + 1}",
                "period_d": float(period),
                "frequency_d^-1": float(freq),
                "a1_mag": float(amp_h1) if np.isfinite(amp_h1) else np.nan,
                "a2_mag": float(amp_h2) if np.isfinite(amp_h2) else np.nan,
                "ptp_mag": float(ptp) if np.isfinite(ptp) else np.nan,
                "relation": str(label.get("relation", "independent")),
                "note": str(label.get("note", "")),
            })
        return rows

    def _update_multimode_mode_views(self) -> None:
        if not hasattr(self, "mm_mode_table"):
            return
        self.mm_mode_table.setRowCount(0)
        self.mm_mode_rows = []
        if not self.multimode_result:
            if hasattr(self, "mm_mode_status"):
                self.mm_mode_status.setText("No joint-fit mode classification yet.")
            return

        rows = self._summarize_multimode_modes(self.multimode_result)
        self.mm_mode_rows = rows
        for row_idx, row in enumerate(rows):
            self.mm_mode_table.insertRow(row_idx)
            values = [
                str(row.get("mode", "")),
                f"{float(row.get('period_d', np.nan)):.8f}",
                f"{float(row.get('frequency_d^-1', np.nan)):.5f}",
                f"{float(row.get('a1_mag', np.nan)):.4f}" if np.isfinite(row.get("a1_mag", np.nan)) else "—",
                f"{float(row.get('a2_mag', np.nan)):.4f}" if np.isfinite(row.get("a2_mag", np.nan)) else "—",
                f"{float(row.get('ptp_mag', np.nan)):.4f}" if np.isfinite(row.get("ptp_mag", np.nan)) else "—",
                str(row.get("relation", "")),
                str(row.get("note", "")),
            ]
            relation = str(row.get("relation", ""))
            for col_idx, value in enumerate(values):
                cell = QTableWidgetItem(value)
                if relation == "independent":
                    cell.setBackground(QColor("#E8F5E9"))
                elif relation == "combination":
                    cell.setBackground(QColor("#FFF3E0"))
                elif relation == "harmonic":
                    cell.setBackground(QColor("#FFF8E1"))
                elif relation == "alias":
                    cell.setBackground(QColor("#FBE9E7"))
                self.mm_mode_table.setItem(row_idx, col_idx, cell)

        if hasattr(self, "mm_mode_status"):
            counts: dict[str, int] = {}
            for row in rows:
                relation = str(row.get("relation", "independent"))
                counts[relation] = counts.get(relation, 0) + 1
            status = (
                f"Mode classification: independent={counts.get('independent', 0)}, "
                f"combination={counts.get('combination', 0)}, "
                f"harmonic={counts.get('harmonic', 0)}, "
                f"alias={counts.get('alias', 0)}"
            )
            independent = [float(row["period_d"]) for row in rows if row.get("relation") == "independent"]
            if len(independent) >= 2:
                p_short = min(independent[0], independent[1])
                p_long = max(independent[0], independent[1])
                if p_long > 0:
                    status += f"  |  P_short/P_long={p_short / p_long:.4f}"
            self.mm_mode_status.setText(status)

    def _current_multimode_export_stem(self) -> str:
        parts = []
        if self.lc_data is not None:
            target_id = self.lc_data.get("target_id")
            if target_id is not None:
                parts.append(f"ID{int(target_id)}")
            analysis_filter = str(self.lc_data.get("analysis_filter") or "").strip()
            if analysis_filter and analysis_filter != "__all__":
                parts.append(analysis_filter)
            corr_tag = str(self.lc_data.get("corr_tag") or "").strip()
            if corr_tag:
                parts.append(corr_tag)
            source = Path(str(self.lc_data.get("source") or "series")).stem
            parts.append(source)
        if not parts:
            parts.append("series")
        return _sanitize_filename_token("_".join(parts))

    def _export_multimode_report(self) -> None:
        if not self.mm_history and not self.mm_mode_rows and not self.multimode_result:
            QMessageBox.information(self, "Multi-Mode", "Export할 multi-mode 결과가 없습니다.")
            return

        try:
            out_dir = self._current_workspace_dir() / "variable_star_tool"
            out_dir.mkdir(parents=True, exist_ok=True)
            stem = self._current_multimode_export_stem()
            written: list[Path] = []

            if self.mm_history:
                history_path = out_dir / f"{stem}_multimode_history.csv"
                pd.DataFrame(self.mm_history).to_csv(history_path, index=False)
                written.append(history_path)
            if self.mm_mode_rows:
                mode_path = out_dir / f"{stem}_multimode_modes.csv"
                pd.DataFrame(self.mm_mode_rows).to_csv(mode_path, index=False)
                written.append(mode_path)
            if self.mm_candidate_scan and self.mm_candidate_scan.get("candidates"):
                cand_path = out_dir / f"{stem}_multimode_candidates.csv"
                pd.DataFrame(self.mm_candidate_scan.get("candidates", [])).to_csv(cand_path, index=False)
                written.append(cand_path)
            if self.multimode_result and self.multimode_result.get("alias_solutions"):
                alias_path = out_dir / f"{stem}_multimode_alias_solutions.csv"
                alias_rows = []
                for row in self.multimode_result.get("alias_solutions", []):
                    alias_rows.append({
                        "rank": row.get("rank"),
                        "periods_d": _format_period_list(row.get("periods", [])),
                        "frequencies_cd": ", ".join(
                            f"{float(value):.8f}" for value in row.get("frequencies_cd", [])
                        ),
                        "bic": row.get("bic"),
                        "delta_bic": row.get("delta_bic"),
                        "rmse": row.get("rmse"),
                    })
                pd.DataFrame(alias_rows).to_csv(alias_path, index=False)
                written.append(alias_path)

            summary_path = out_dir / f"{stem}_multimode_summary.txt"
            summary_lines = [
                self.mm_summary_label.text().strip() if hasattr(self, "mm_summary_label") else "",
                "",
                self.mm_mode_status.text().strip() if hasattr(self, "mm_mode_status") else "",
                self.mm_candidate_status.text().strip() if hasattr(self, "mm_candidate_status") else "",
                self.mm_history_status.text().strip() if hasattr(self, "mm_history_status") else "",
            ]
            summary_path.write_text("\n".join(line for line in summary_lines if line), encoding="utf-8")
            written.append(summary_path)
        except Exception as e:
            QMessageBox.warning(self, "Multi-Mode Export", str(e))
            return

        self.log(f"[MULTI] Exported report to {out_dir}")
        QMessageBox.information(
            self,
            "Multi-Mode Export",
            "Saved:\n" + "\n".join(str(path.name) for path in written),
        )

    def _sync_phase_mode_combo(self):
        if not hasattr(self, "phase_mode_combo"):
            return
        current_data = self.phase_mode_combo.currentData()
        self.phase_mode_combo.blockSignals(True)
        self.phase_mode_combo.clear()
        self.phase_mode_combo.addItem("Manual / custom period", "manual")
        periods = []
        if self.multimode_result:
            periods = [float(p) for p in self.multimode_result.get("periods", [])]
        for idx, period in enumerate(periods):
            self.phase_mode_combo.addItem(f"M{idx + 1} | P={period:.8f} d", idx)
        self.phase_mode_combo.setEnabled(bool(periods))
        restore_idx = self.phase_mode_combo.findData(current_data)
        if restore_idx < 0:
            restore_idx = 0
        self.phase_mode_combo.setCurrentIndex(restore_idx)
        self.phase_mode_combo.blockSignals(False)

    def _selected_phase_mode_index(self) -> int | None:
        if not hasattr(self, "phase_mode_combo"):
            return None
        data = self.phase_mode_combo.currentData()
        try:
            idx = int(data)
        except Exception:
            return None
        if not self.multimode_result:
            return None
        periods = [float(p) for p in self.multimode_result.get("periods", [])]
        if 0 <= idx < len(periods):
            return idx
        return None

    def _on_phase_mode_changed(self):
        idx = self._selected_phase_mode_index()
        if idx is None or not self.multimode_result:
            self._update_phase_plot()
            return
        period = float(self.multimode_result["periods"][idx])
        self.phase_p.blockSignals(True)
        self.phase_p.setValue(period)
        self.phase_p.blockSignals(False)
        mode_idx = self.phase_data_mode.findData("focused")
        if mode_idx >= 0:
            self.phase_data_mode.blockSignals(True)
            self.phase_data_mode.setCurrentIndex(mode_idx)
            self.phase_data_mode.blockSignals(False)
        if hasattr(self, "mm_focus_combo"):
            combo_idx = self.mm_focus_combo.findData(idx)
            if combo_idx >= 0:
                self.mm_focus_combo.blockSignals(True)
                self.mm_focus_combo.setCurrentIndex(combo_idx)
                self.mm_focus_combo.blockSignals(False)
        self._draw_multimode()
        self._update_phase_plot()

    def _on_phase_period_changed(self):
        idx = self._selected_phase_mode_index()
        if idx is not None and self.multimode_result:
            try:
                mode_period = float(self.multimode_result["periods"][idx])
            except Exception:
                mode_period = np.nan
            current_period = float(self.phase_p.value())
            if np.isfinite(mode_period) and np.isfinite(current_period):
                rel = abs(current_period - mode_period) / max(abs(mode_period), 1e-12)
                if rel > 1e-5 and hasattr(self, "phase_mode_combo"):
                    manual_idx = self.phase_mode_combo.findData("manual")
                    if manual_idx >= 0:
                        self.phase_mode_combo.blockSignals(True)
                        self.phase_mode_combo.setCurrentIndex(manual_idx)
                        self.phase_mode_combo.blockSignals(False)
        self._update_phase_plot()

    def _phase_multimode_payload(
        self,
        t_v: np.ndarray,
        m_v: np.ndarray,
        current_period: float,
        night_id: np.ndarray | None = None,
    ) -> dict | None:
        if not self.multimode_result or not np.isfinite(current_period) or current_period <= 0:
            return None
        periods = [float(p) for p in self.multimode_result.get("periods", [])]
        if not periods:
            return None
        idx = self._selected_phase_mode_index()
        if idx is None:
            rel = np.asarray([abs(current_period - p) / max(abs(p), 1e-12) for p in periods], dtype=float)
            idx = int(np.argmin(rel))
            if not np.isfinite(rel[idx]) or rel[idx] > 0.03:
                return None
        eval_obs = _evaluate_multimode_result(self.multimode_result, t_v)
        focus_component = np.asarray(eval_obs["components"][idx], dtype=float)
        iso_mag = m_v - (eval_obs["total"] - (float(eval_obs["intercept"]) + focus_component))
        if night_id is not None:
            offsets = self.multimode_result.get("night_offsets", {})
            iso_mag = iso_mag - np.asarray(
                [float(offsets.get(str(value), 0.0)) for value in night_id],
                dtype=float,
            )
        return {
            "index": idx,
            "period": periods[idx],
            "iso_mag": iso_mag,
        }

    def _append_periods_to_multimode(self, periods_to_add: list[float]) -> list[float]:
        current = _parse_period_list(self.mm_periods_edit.text())
        added: list[float] = []
        for period in periods_to_add:
            try:
                value = float(period)
            except Exception:
                continue
            if not np.isfinite(value) or value <= 0:
                continue
            if any(abs(value - prev) / max(abs(prev), 1e-12) < 1e-5 for prev in current):
                continue
            current.append(value)
            added.append(value)
        if added:
            self.mm_periods_edit.setText(_format_period_list(current))
        return added

    def _build_multimode_scan_payload(self, periods: list[float] | None = None) -> dict | None:
        if self.lc_data is None:
            return None
        adopted_periods = periods if periods is not None else _parse_period_list(self.mm_periods_edit.text())
        time = np.asarray(self.lc_data["time"], dtype=float)
        mag = np.asarray(self.lc_data["mag"], dtype=float)
        mag_err = self.lc_data.get("mag_err")
        err = np.asarray(mag_err, dtype=float) if mag_err is not None else None
        night_data = self.lc_data.get("night_id")
        night_ids = np.asarray(night_data, dtype=object) if night_data is not None else None
        mask = np.isfinite(time) & np.isfinite(mag)
        if err is not None:
            mask &= np.isfinite(err)
        t = time[mask]
        y = mag[mask]
        dy = err[mask] if err is not None else None
        nid = night_ids[mask] if night_ids is not None else None
        if len(t) < 10:
            return None

        window = compute_spectral_window(t)
        window_peaks = list(window.get("peaks", []))
        baseline_days = float(window.get("baseline_days", np.ptp(t)))

        fit_result = None
        signal = y.copy()
        stage_label = "raw"
        if adopted_periods:
            fit_result = _fit_multimode_model(
                t,
                y,
                dy,
                periods=adopted_periods,
                harmonics=int(self.mm_harm.value()),
                night_id=nid,
                include_night_offsets=self.mm_night_offsets_chk.isChecked(),
            )
            signal = np.asarray(fit_result["residual"], dtype=float)
            stage_label = "residual"

        results = run_period_analysis(
            time=t,
            mag_raw=signal,
            mag_corr=None,
            mag_err=dy,
            min_period=self.min_p.value(),
            max_period=self.max_p.value(),
            samples_per_peak=self.spp.value(),
            methods=["ls"],
            pdm_n_bins=self.pdm_bins.value(),
        )
        ls_result = results.get("raw_ls")
        if not ls_result or "error" in ls_result:
            return {
                "stage": stage_label,
                "periods": list(adopted_periods),
                "fit_result": fit_result,
                "time": t,
                "signal": signal,
                "ls_result": ls_result,
                "candidates": [],
            }

        top_periods = [float(p) for p in ls_result.get("top_periods", []) if np.isfinite(p) and float(p) > 0]
        top_powers = [float(v) for v in ls_result.get("top_powers", []) if np.isfinite(v)]
        candidates = []
        for idx, period in enumerate(top_periods[: max(3, int(self.mm_n_modes.value()) + 2)]):
            power = top_powers[idx] if idx < len(top_powers) else np.nan
            relation, note = _classify_candidate_period(
                period,
                adopted_periods,
                window_peaks=window_peaks,
                baseline_days=baseline_days,
            )
            candidates.append(
                {
                    "rank": idx + 1,
                    "period": period,
                    "power": power,
                    "relation": relation,
                    "note": note,
                }
            )

        return {
            "stage": stage_label,
            "periods": list(adopted_periods),
            "fit_result": fit_result,
            "time": t,
            "signal": signal,
            "ls_result": ls_result,
            "candidates": candidates,
            "window_peaks": window_peaks,
            "baseline_days": baseline_days,
        }

    def _detect_multimode_candidates(self):
        if self.lc_data is None:
            QMessageBox.warning(self, "Multi-Mode", "Load a light curve first.")
            return
        self.workflow_step = "multi"
        self._refresh_tool_workflow_ui()
        try:
            payload = self._build_multimode_scan_payload()
        except Exception as e:
            QMessageBox.warning(self, "Multi-Mode", str(e))
            self.mm_candidate_status.setText(f"Candidate detection failed: {e}")
            return
        self.mm_candidate_scan = payload
        self._update_multimode_candidate_views()
        ls_result = payload.get("ls_result") if payload else None
        candidates = payload.get("candidates", []) if payload else []
        best_period = float(ls_result.get("best_period", np.nan)) if ls_result else np.nan
        best_power = float(ls_result.get("best_power", np.nan)) if ls_result else np.nan
        self._record_multimode_history(
            "detect_candidates",
            stage=str(payload.get("stage", "")) if payload else "",
            periods=[float(p) for p in payload.get("periods", [])] if payload else [],
            best_period=best_period,
            best_power=best_power,
            note=f"new={sum(1 for item in candidates if item.get('relation') == 'new')}, total={len(candidates)}",
        )
        if np.isfinite(best_period):
            self.log(
                f"[MULTI] Candidate scan ({payload.get('stage', 'raw')}): "
                f"best={best_period:.8f} d, new={sum(1 for item in candidates if item.get('relation') == 'new')}"
            )

    def _selected_candidate_periods(self) -> list[float]:
        if not self.mm_candidate_scan:
            return []
        rows = sorted({idx.row() for idx in self.mm_candidate_table.selectedIndexes()})
        candidates = self.mm_candidate_scan.get("candidates", [])
        periods = []
        for row in rows:
            if 0 <= row < len(candidates):
                periods.append(float(candidates[row]["period"]))
        return periods

    def _append_selected_candidate_to_multimode(self):
        periods = self._selected_candidate_periods()
        if not periods and self.mm_candidate_scan:
            for item in self.mm_candidate_scan.get("candidates", []):
                if item.get("relation") == "new":
                    periods = [float(item["period"])]
                    break
        if not periods:
            QMessageBox.information(self, "Multi-Mode", "Adopt할 candidate가 없습니다.")
            return
        added = self._append_periods_to_multimode(periods)
        if not added:
            QMessageBox.information(self, "Multi-Mode", "선택한 candidate는 이미 period list에 있습니다.")
            return
        self._invalidate_multimode_runtime(f"Added {len(added)} candidate period(s). Re-run Detect/Fit.")
        merged_periods = _parse_period_list(self.mm_periods_edit.text())
        self._record_multimode_history(
            "adopt_selected",
            stage="candidate",
            periods=merged_periods,
            best_period=added[0] if added else np.nan,
            note=f"added {len(added)} selected candidate(s)",
        )

    def _adopt_all_new_candidates(self):
        if not self.mm_candidate_scan:
            QMessageBox.information(self, "Multi-Mode", "먼저 residual candidate detection을 실행하세요.")
            return
        periods = [float(item["period"]) for item in self.mm_candidate_scan.get("candidates", []) if item.get("relation") == "new"]
        if not periods:
            QMessageBox.information(self, "Multi-Mode", "추가할 new candidate가 없습니다.")
            return
        added = self._append_periods_to_multimode(periods)
        if not added:
            QMessageBox.information(self, "Multi-Mode", "추가할 new candidate가 이미 모두 들어 있습니다.")
            return
        self._invalidate_multimode_runtime(f"Added {len(added)} new candidate period(s). Re-run Fit Multi-Mode.")
        merged_periods = _parse_period_list(self.mm_periods_edit.text())
        self._record_multimode_history(
            "adopt_all_new",
            stage="candidate",
            periods=merged_periods,
            best_period=added[0] if added else np.nan,
            note=f"added {len(added)} new candidate(s)",
        )

    def _update_multimode_candidate_views(self):
        if not hasattr(self, "mm_candidate_table") or not hasattr(self, "mm_pw_canvas"):
            return

        payload = self.mm_candidate_scan or {}
        candidates = payload.get("candidates", [])
        adopted_periods = [float(p) for p in payload.get("periods", [])]
        stage = str(payload.get("stage", "raw"))
        ls_result = payload.get("ls_result")

        self.mm_candidate_table.setRowCount(0)
        for row_idx, item in enumerate(candidates):
            self.mm_candidate_table.insertRow(row_idx)
            values = [
                str(item.get("rank", row_idx + 1)),
                f"{float(item.get('period', np.nan)):.8f}",
                f"{float(item.get('power', np.nan)):.4f}" if np.isfinite(item.get("power", np.nan)) else "—",
                str(item.get("relation", "")),
                str(item.get("note", "")),
            ]
            for col_idx, value in enumerate(values):
                cell = QTableWidgetItem(value)
                relation = str(item.get("relation", ""))
                if relation == "new":
                    cell.setBackground(QColor("#E8F5E9"))
                elif relation in {"alias", "harmonic"}:
                    cell.setBackground(QColor("#FFF8E1"))
                elif relation == "duplicate":
                    cell.setBackground(QColor("#ECEFF1"))
                self.mm_candidate_table.setItem(row_idx, col_idx, cell)

        if hasattr(self, "mm_candidate_status"):
            if not ls_result or "error" in ls_result:
                self.mm_candidate_status.setText("Candidate detection failed or no LS result available.")
            else:
                best_period = float(ls_result.get("best_period", np.nan))
                best_power = float(ls_result.get("best_power", np.nan))
                self.mm_candidate_status.setText(
                    f"{stage.title()} scan: best residual peak {best_period:.8f} d  "
                    f"(power={best_power:.4f}), adopted={len(adopted_periods)}"
                )

        fig = self.mm_pw_canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)
        if not ls_result or "error" in ls_result:
            ax.text(0.5, 0.5, "No candidate scan available", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            _safe_tight_layout(fig)
            self.mm_pw_canvas.draw_idle()
            return

        periods = 1.0 / np.asarray(ls_result["frequency"], dtype=float)
        power = np.asarray(ls_result["power"], dtype=float)
        ax.plot(periods, power, color="#455A64", lw=0.9, alpha=0.95)
        ax.set_xscale("log")
        ax.set_xlabel("Period (days)")
        ax.set_ylabel("LS Power")
        ax.set_title(f"Multi-path candidate scan ({stage})")
        ax.grid(True, alpha=0.3)

        best_period = float(ls_result.get("best_period", np.nan))
        best_power = float(ls_result.get("best_power", np.nan))
        if np.isfinite(best_period):
            ax.axvline(best_period, color="#D32F2F", ls="--", lw=1.3, alpha=0.85, label=f"best {best_period:.6f} d")
            if np.isfinite(best_power):
                ax.scatter([best_period], [best_power], color="#D32F2F", s=36, zorder=5)

        for idx, period in enumerate(adopted_periods):
            ax.axvline(float(period), color="#1565C0", ls=":", lw=1.0, alpha=0.7,
                       label="adopted" if idx == 0 else None)
        for idx, item in enumerate(candidates[:4]):
            if item.get("relation") != "new":
                continue
            ax.axvline(float(item["period"]), color="#FB8C00", ls="-.", lw=1.0, alpha=0.7,
                       label="new candidate" if idx == 0 else None)
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(loc="best", fontsize=8)
        _safe_tight_layout(fig)
        self.mm_pw_canvas.draw_idle()

    def _append_current_period_to_multimode(self):
        current = np.nan
        if self.refined_period is not None and np.isfinite(self.refined_period):
            current = float(self.refined_period)
        elif self.scan_result:
            current, _, _ = self._best_from_results(self.scan_result)
        elif np.isfinite(self.phase_p.value()):
            current = float(self.phase_p.value())
        if not np.isfinite(current) or current <= 0:
            QMessageBox.information(self, "Multi-Mode", "추가할 주기가 없습니다.")
            return
        periods = _parse_period_list(self.mm_periods_edit.text())
        added = False
        if not any(abs(current - p) / max(abs(p), 1e-12) < 1e-5 for p in periods):
            periods.append(current)
            added = True
        self.mm_periods_edit.setText(_format_period_list(periods))
        if not added:
            self.mm_status.setText("Current dominant period is already in the multi-mode list.")
            return
        self._invalidate_multimode_runtime("Added current dominant period. Re-run Detect/Fit.")
        self._record_multimode_history(
            "add_best",
            stage="manual",
            periods=periods,
            best_period=current,
            note="added current dominant period to multi-mode list",
        )

    def _set_multimode_periods_from_scan(self):
        periods: list[float] = []
        analysis = self.scan_alias_analysis or {}
        adopted = float(analysis.get("adopted_period", np.nan))
        if np.isfinite(adopted) and adopted > 0:
            periods.append(adopted)
        raw_periods = _scan_top_periods(self.scan_result, max_count=12)
        window_peaks = list(analysis.get("window_peaks", []))
        baseline_days = float(analysis.get("baseline_days", 1.0))
        for period in raw_periods:
            if len(periods) >= int(self.mm_n_modes.value()):
                break
            relation, _ = _classify_candidate_period(
                period,
                periods,
                window_peaks=window_peaks,
                baseline_days=baseline_days,
            )
            if relation == "new":
                periods.append(float(period))
        if not periods:
            QMessageBox.information(self, "Multi-Mode", "먼저 period scan을 실행하세요.")
            return
        self.mm_periods_edit.setText(_format_period_list(periods))
        self._invalidate_multimode_runtime(f"Loaded {len(periods)} peak period(s) from scan.")
        self._record_multimode_history(
            "seed_from_scan",
            stage="scan",
            periods=periods,
            best_period=periods[0] if periods else np.nan,
            note="loaded top scan peaks into multi-mode list",
        )

    def _selected_multimode_focus_index(self) -> int | None:
        if not self.multimode_result:
            return None
        periods = [float(p) for p in self.multimode_result.get("periods", [])]
        if not periods:
            return None
        data = self.mm_focus_combo.currentData() if hasattr(self, "mm_focus_combo") else None
        try:
            idx = int(data)
            if 0 <= idx < len(periods):
                return idx
        except Exception:
            pass
        phase_idx = self._selected_phase_mode_index()
        if phase_idx is not None:
            return phase_idx
        phase_period = float(self.phase_p.value())
        if np.isfinite(phase_period) and phase_period > 0:
            rel = [abs(phase_period - p) / max(abs(p), 1e-12) for p in periods]
            return int(np.argmin(rel))
        return 0

    def _on_multimode_focus_changed(self):
        if self.mm_overlay_chk.isChecked():
            self._update_phase_plot()

    def _apply_multimode_focus_to_phase(self):
        idx = self._selected_multimode_focus_index()
        if idx is None or not self.multimode_result:
            return
        if hasattr(self, "phase_mode_combo"):
            combo_idx = self.phase_mode_combo.findData(idx)
            if combo_idx >= 0:
                self.phase_mode_combo.setCurrentIndex(combo_idx)
                return
        self.phase_p.setValue(float(self.multimode_result["periods"][idx]))

    def _run_multimode_fit(self):
        if self.lc_data is None:
            QMessageBox.warning(self, "No Data", "Load a light curve first.")
            return

        periods = _parse_period_list(self.mm_periods_edit.text())
        if not periods:
            QMessageBox.warning(self, "Multi-Mode", "한 개 이상 period를 입력하세요.")
            return
        if len(periods) > 6:
            QMessageBox.warning(self, "Multi-Mode", "Mode 수는 6개 이하로 제한하세요.")
            return

        try:
            self.workflow_step = "multi"
            self._refresh_tool_workflow_ui()
            self.mm_status.setText("Fitting…")
            window = compute_spectral_window(self.lc_data["time"])
            if self.mm_alias_search_chk.isChecked():
                result = search_multimode_alias_solutions(
                    self.lc_data["time"],
                    self.lc_data["mag"],
                    self.lc_data.get("mag_err"),
                    seed_periods=periods,
                    harmonics=int(self.mm_harm.value()),
                    window_peaks=window.get("peaks", []),
                    night_id=self.lc_data.get("night_id"),
                    include_night_offsets=self.mm_night_offsets_chk.isChecked(),
                    max_alias_offsets=2 if len(periods) <= 3 else 1,
                    max_solutions=1024,
                )
                periods = [float(p) for p in result["periods"]]
                self.mm_periods_edit.setText(_format_period_list(periods))
            else:
                result = _fit_multimode_model(
                    self.lc_data["time"],
                    self.lc_data["mag"],
                    self.lc_data.get("mag_err"),
                    periods=periods,
                    harmonics=int(self.mm_harm.value()),
                    night_id=self.lc_data.get("night_id"),
                    include_night_offsets=self.mm_night_offsets_chk.isChecked(),
                )
                result["window_peaks"] = list(window.get("peaks", []))
        except Exception as e:
            self.mm_status.setText("Fit failed")
            QMessageBox.warning(self, "Multi-Mode Fit", str(e))
            return

        self.multimode_result = result
        self.mm_focus_combo.blockSignals(True)
        self.mm_focus_combo.clear()
        for idx, period in enumerate(result["periods"]):
            self.mm_focus_combo.addItem(f"M{idx + 1} | P={float(period):.8f} d", idx)
        focus_idx = self._selected_multimode_focus_index()
        if focus_idx is None:
            focus_idx = 0
        self.mm_focus_combo.setCurrentIndex(max(0, min(focus_idx, self.mm_focus_combo.count() - 1)))
        self.mm_focus_combo.setEnabled(self.mm_focus_combo.count() > 0)
        self.mm_focus_combo.blockSignals(False)
        self._sync_phase_mode_combo()
        self._update_multimode_mode_views()
        alias_status = str(result.get("alias_status", "FIXED"))
        self.mm_status.setText(
            f"Fitted {len(result['periods'])} mode(s), {int(result['harmonics'])} "
            f"harmonic(s)/mode | alias status: {alias_status}."
        )
        if result.get("alias_solutions"):
            delta = float(result.get("alias_runner_delta_bic", np.nan))
            self.mm_status.setText(
                self.mm_status.text() + f" Runner-up delta BIC={delta:.2f}."
            )
            reason = str(result.get("alias_status_reason", "")).strip()
            if reason:
                self.mm_status.setText(self.mm_status.text() + f" {reason}")
        cond = float(result.get("design_condition", np.nan))
        n_points = int(result.get("n_points", 0))
        n_params = int(result.get("n_params", 0))
        if (np.isfinite(cond) and cond > 1e8) or (n_params > 0 and n_points < n_params * 4):
            self.mm_status.setText(
                self.mm_status.text() + " Warning: fit may be underconstrained/overfit for this baseline."
            )
        self._draw_multimode()
        try:
            self.mm_candidate_scan = self._build_multimode_scan_payload(periods=list(result["periods"]))
        except Exception:
            self.mm_candidate_scan = None
        self._update_multimode_candidate_views()
        self._record_multimode_history(
            "fit",
            stage="joint_fit",
            periods=[float(p) for p in result["periods"]],
            best_period=float(result["periods"][0]) if result.get("periods") else np.nan,
            note=(
                f"rmse={float(result.get('rmse', np.nan)):.5f}, "
                f"independent={sum(1 for row in self.mm_mode_rows if row.get('relation') == 'independent')}"
            ),
            fit_result=result,
        )
        if self.mm_overlay_chk.isChecked():
            self._update_phase_plot()
        self.log(
            f"Multi-mode fit: {len(result['periods'])} modes, "
            f"periods={_format_period_list(result['periods'])}"
        )

    def _draw_multimode(self):
        if not hasattr(self, "mm_canvas"):
            return
        fig = self.mm_canvas.figure
        fig.clear()

        if not self.multimode_result:
            ax = fig.add_subplot(111)
            ax.text(0.5, 0.5, "No multi-mode fit yet", ha="center", va="center",
                    transform=ax.transAxes, fontsize=12, color="gray")
            ax.set_axis_off()
            _safe_tight_layout(fig)
            self.mm_canvas.draw_idle()
            return

        result = self.multimode_result
        periods = [float(p) for p in result["periods"]]
        t = np.asarray(result["time"], dtype=float)
        y = np.asarray(result["mag"], dtype=float)
        model = np.asarray(result["model"], dtype=float)
        residual = np.asarray(result["residual"], dtype=float)

        obs_start = float(np.min(t))
        obs_end = float(np.max(t))
        obs_span = max(obs_end - obs_start, 1e-9)
        baseline = float(result.get("baseline", obs_span))
        longest_period = max(periods)
        proj_span = max(baseline, longest_period * 3.0)
        beat_lines: list[str] = []
        freqs = [1.0 / p for p in periods]
        for i in range(len(freqs)):
            for j in range(i + 1, len(freqs)):
                dfreq = abs(freqs[i] - freqs[j])
                if dfreq <= 1e-9:
                    continue
                beat = 1.0 / dfreq
                if np.isfinite(beat):
                    proj_span = max(proj_span, min(beat, longest_period * 12.0))
                    beat_lines.append(f"M{i + 1}-M{j + 1} beat≈{beat:.5f} d ({beat*24:.2f} h)")

        n_grid_obs = max(600, min(4000, int(obs_span * 4000)))
        t_dense_obs = np.linspace(obs_start, obs_end, n_grid_obs)
        obs_eval = _evaluate_multimode_result(result, t_dense_obs)
        total_dense_obs = obs_eval["total"]

        n_grid_proj = max(1200, min(5000, int(proj_span * 4000)))
        t_dense_proj = np.linspace(obs_start, obs_start + proj_span, n_grid_proj)
        proj_eval = _evaluate_multimode_result(result, t_dense_proj)
        total_dense_proj = proj_eval["total"]

        obs_bright = float(np.nanmin(total_dense_obs))
        obs_faint = float(np.nanmax(total_dense_obs))
        obs_ptp = float(obs_faint - obs_bright)
        proj_bright = float(np.nanmin(total_dense_proj))
        proj_faint = float(np.nanmax(total_dense_proj))
        proj_ptp = float(proj_faint - proj_bright)

        t0 = float(self.t0_edit.value())
        state_epoch = float(self.mm_state_epoch.value())
        state_eval = _evaluate_multimode_result(result, np.array([state_epoch], dtype=float))
        total_state = float(state_eval["total"][0])
        total_state_trend = _trend_label(float(state_eval["total_derivative"][0]))

        component_lines: list[str] = []
        phase_panels: list[dict] = []
        for idx, period in enumerate(periods):
            phase_grid = np.linspace(0.0, 1.0, 500)
            phase_times = t0 + phase_grid * period
            phase_eval = _evaluate_multimode_result(result, phase_times)
            comp_curve = phase_eval["intercept"] + phase_eval["components"][idx]
            comp_ptp = float(np.nanmax(comp_curve) - np.nanmin(comp_curve))
            comp_phase = ((state_epoch - t0) / period) % 1.0 if period > 0 else np.nan
            comp_state = float(state_eval["components"][idx][0])
            comp_trend = _trend_label(float(state_eval["component_derivatives"][idx][0]))
            component_lines.append(
                f"M{idx + 1}: P={period:.8f} d  f={1.0/period:.5f} d^-1  "
                f"Δm≈{comp_ptp:.4f}  φ={comp_phase:.3f}  contrib={comp_state:+.4f}  {comp_trend}"
            )

            component_obs = np.asarray(result["components"][idx], dtype=float)
            other_model = model - (float(result["intercept"]) + component_obs)
            iso_mag = y - other_model
            obs_phase = ((t - t0) / period) % 1.0
            phase_model = np.linspace(0.0, 2.0, 800)
            phase_model_times = t0 + (phase_model % 1.0) * period
            phase_model_eval = _evaluate_multimode_result(result, phase_model_times)
            model_curve = phase_model_eval["intercept"] + phase_model_eval["components"][idx]
            phase_panels.append({
                "index": idx,
                "period": period,
                "phase_ext": np.concatenate([obs_phase, obs_phase + 1.0]),
                "iso_ext": np.concatenate([iso_mag, iso_mag]),
                "phase_model": phase_model,
                "model_curve": model_curve,
            })

        summary_lines = [
            f"[{len(periods)} mode(s) | {int(result['harmonics'])} harmonic(s)/mode | "
            f"n={int(result.get('n_points', len(t)))} pts | p={int(result.get('n_params', 0))} params]",
            f"Observed-window fit range: {obs_bright:.4f} .. {obs_faint:.4f} mag  (Δm={obs_ptp:.4f})",
            f"Projected beat range: {proj_bright:.4f} .. {proj_faint:.4f} mag  (Δm={proj_ptp:.4f})",
            f"Fit RMS={float(result.get('rmse', np.nan)):.5f} mag  "
            f"WRMS={float(result.get('wrms', np.nan)):.3f}  cond≈{float(result.get('design_condition', np.nan)):.2e}",
            f"State @ {state_epoch:.6f} BJD: model={total_state:.4f} mag  ({total_state_trend})",
        ]
        summary_lines.extend(component_lines)
        summary_lines.extend(beat_lines[:4])
        self.mm_summary_label.setText("\n".join(summary_lines))

        show_phase_panels = min(len(phase_panels), 4)
        phase_cols = 2 if show_phase_panels > 1 else 1
        phase_rows = max(1, (show_phase_panels + phase_cols - 1) // phase_cols)
        gs = fig.add_gridspec(2 + phase_rows, phase_cols, height_ratios=[1.25, 0.9] + [1.0] * phase_rows)
        ax_time = fig.add_subplot(gs[0, :])
        ax_resid = fig.add_subplot(gs[1, :])
        phase_axes = [
            fig.add_subplot(gs[2 + (panel_idx // phase_cols), panel_idx % phase_cols])
            for panel_idx in range(show_phase_panels)
        ]

        show_projection = bool(self.mm_project_chk.isChecked()) if hasattr(self, "mm_project_chk") else False

        ax_time.scatter(t, y, s=12, color="#607D8B", alpha=0.7, label="Data")
        ax_time.plot(t_dense_obs, total_dense_obs, color="#D32F2F", lw=1.6, label="Observed-window fit")
        if show_projection and len(t_dense_proj) > len(t_dense_obs):
            proj_mask = t_dense_proj > obs_end
            if np.any(proj_mask):
                ax_time.plot(
                    t_dense_proj[proj_mask],
                    total_dense_proj[proj_mask],
                    color="#D32F2F",
                    lw=1.3,
                    ls="--",
                    alpha=0.85,
                    label="Projected beat window",
                )
                ax_time.axvspan(obs_start, obs_end, color="#ECEFF1", alpha=0.35, zorder=0)
                ax_time.axvline(obs_end, color="#9E9E9E", ls=":", lw=1.0, alpha=0.7)
                ax_time.set_xlim(obs_start, float(t_dense_proj[-1]))
            else:
                ax_time.set_xlim(obs_start, obs_end)
        else:
            ax_time.set_xlim(obs_start, obs_end)
        ax_time.invert_yaxis()
        ax_time.set_xlabel("Time (BJD)")
        ax_time.set_ylabel("Magnitude")
        if show_projection:
            ax_time.set_title("Observed JD + Projected Beat Window")
        else:
            ax_time.set_title("Observed JD Window + Total Multi-Mode Fit")
        ax_time.grid(True, alpha=0.3)
        ax_time.legend(fontsize=8)

        ax_resid.scatter(t, residual, s=12, color="#3949AB", alpha=0.7)
        ax_resid.axhline(0, color="gray", ls="--", lw=1.0, alpha=0.6)
        ax_resid.set_xlabel("Time (BJD)")
        ax_resid.set_ylabel("Residual (mag)")
        ax_resid.set_title("Residual After Total Fit")
        ax_resid.grid(True, alpha=0.3)

        for ax_phase, panel in zip(phase_axes, phase_panels[:show_phase_panels]):
            idx = int(panel["index"])
            ax_phase.scatter(
                panel["phase_ext"],
                panel["iso_ext"],
                s=10,
                color="#1E88E5",
                alpha=0.72,
                label=f"M{idx + 1} isolated",
            )
            ax_phase.plot(
                panel["phase_model"],
                panel["model_curve"],
                color="#C62828",
                lw=1.5,
                label="Fit",
            )
            ax_phase.invert_yaxis()
            ax_phase.set_xlim(0, 2)
            ax_phase.set_xlabel("Phase")
            ax_phase.set_ylabel("Magnitude")
            ax_phase.set_title(f"M{idx + 1} Isolated Phase (P={float(panel['period']):.8f} d)")
            ax_phase.axvline(0, color="gray", ls=":", alpha=0.4)
            ax_phase.axvline(1, color="gray", ls=":", alpha=0.4)
            ax_phase.grid(True, alpha=0.3)
            ax_phase.legend(fontsize=8)

        _safe_tight_layout(fig)
        self.mm_canvas.draw_idle()

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _current_workspace_dir(self) -> Path:
        text = self.workspace_edit.text().strip() if hasattr(self, "workspace_edit") else ""
        path = Path(text) if text else Path(self.workspace_dir)
        self.workspace_dir = path
        return path

    def _browse_workspace(self):
        start_dir = self._current_workspace_dir()
        start = str(start_dir.parent if start_dir.exists() else Path(self.params.P.result_dir).parent)
        path = QFileDialog.getExistingDirectory(self, "Workspace 선택", start)
        if path:
            self.workspace_edit.setText(path)
            self._load_lc_from_workspace()

    def _load_lc_from_workspace(self, *, synchronous: bool = False):
        """Load all available workspace outputs and expose them in the Use data combo."""
        try:
            from apex.utils.step_paths_lc import list_lightcurve_csvs
            rd = self._current_workspace_dir()
            paths = list_lightcurve_csvs(rd)
            if not paths:
                self._clear_loaded_workspace_state(f"No lightcurve_*.csv found in\n{rd}")
                return
            if synchronous:
                self._load_paths(paths)
                return
            self._workspace_load_token += 1
            self._pending_workspace_load = (
                self._workspace_load_token,
                rd,
                paths,
            )
            self.lc_status.setText(f"Loading workspace...\n{rd}")
            worker = self._workspace_load_worker
            if worker is not None and worker.isRunning():
                worker.requestInterruption()
                return
            self._start_pending_workspace_load()
        except Exception as e:
            self._clear_loaded_workspace_state(f"Workspace load failed: {e}")

    def _start_pending_workspace_load(self) -> None:
        pending = self._pending_workspace_load
        if pending is None:
            return
        self._pending_workspace_load = None
        token, workspace, paths = pending
        try:
            from apex.utils.qc_utils import load_frame_excludes as _lfe
            excluded = set(_lfe(workspace).keys())
        except Exception:
            excluded = set()
        worker = _WorkspaceLoadWorker(self, token, paths, excluded)
        self._workspace_load_worker = worker
        worker.result_ready.connect(self._on_workspace_load_ready)
        worker.error.connect(self._on_workspace_load_error)
        worker.finished.connect(
            lambda current=worker: self._on_workspace_load_finished(current)
        )
        worker.start()

    def _on_workspace_load_ready(self, token: int, items: object) -> None:
        if token != self._workspace_load_token or not isinstance(items, list):
            return
        self._apply_loaded_series_items(items, self._current_workspace_dir())

    def _on_workspace_load_error(self, token: int, message: str) -> None:
        if token == self._workspace_load_token:
            self._clear_loaded_workspace_state(f"Workspace load failed: {message}")

    def _on_workspace_load_finished(self, worker: _WorkspaceLoadWorker) -> None:
        if self._workspace_load_worker is worker:
            self._workspace_load_worker = None
        worker.deleteLater()
        if self._pending_workspace_load is not None:
            self._start_pending_workspace_load()

    def apply_period_handoff(
        self,
        handoff: ValidatedLightCurveBundle | dict,
    ) -> None:
        """Load a Step 12 result without treating unresolved aliases as final."""
        bundle = coerce_validated_bundle(handoff)
        if not bundle.can_launch:
            raise ValueError(
                "Main workflow blocked advanced variable analysis: "
                f"{bundle.release_message}"
            )
        handoff = bundle.to_legacy_handoff()
        for worker in (self._scan_worker, self._refine_worker, self._analysis_worker):
            if worker is not None and worker.isRunning():
                raise RuntimeError("Wait for the current Variable Star Tool analysis to finish.")

        workspace = Path(str(handoff.get("workspace_dir", "") or ""))
        if not workspace.exists():
            raise FileNotFoundError(f"Handoff workspace not found: {workspace}")
        self.workspace_dir = workspace
        self.workspace_edit.setText(str(workspace))
        self._load_lc_from_workspace(synchronous=True)
        if not self.series_options:
            raise ValueError("No usable light-curve series was found in the handoff workspace.")
        self.handoff_bundle = bundle

        series_key = _select_period_handoff_series_key(self.series_options, handoff)
        if series_key is None:
            raise ValueError("No light-curve series matches the Step 12 handoff.")
        selected_item = self.series_options[series_key]

        series_index = self.mag_col_combo.findData(series_key)
        if series_index >= 0:
            self.mag_col_combo.blockSignals(True)
            self.mag_col_combo.setCurrentIndex(series_index)
            self.mag_col_combo.blockSignals(False)

        requested_filter = str(handoff.get("analysis_filter", "__all__") or "__all__")
        if requested_filter == "all":
            requested_filter = "__all__"
        self.analysis_filter_combo.blockSignals(True)
        self.analysis_filter_combo.clear()
        self.analysis_filter_combo.addItem(
            "All" if requested_filter == "__all__" else requested_filter,
            requested_filter,
        )
        self.analysis_filter_combo.setCurrentIndex(0)
        self.analysis_filter_combo.blockSignals(False)
        self._apply_series_option(series_key)
        if self.lc_data is None:
            raise ValueError("The selected handoff light curve could not be loaded.")

        search = dict(handoff.get("search") or {})
        try:
            min_period = float(search.get("min_period", np.nan))
        except (TypeError, ValueError):
            min_period = np.nan
        try:
            max_period = float(search.get("max_period", np.nan))
        except (TypeError, ValueError):
            max_period = np.nan
        if np.isfinite(min_period):
            self.min_p.setValue(min_period)
        if np.isfinite(max_period):
            self.max_p.setValue(max_period)
        if search.get("samples_per_peak") is not None:
            self.spp.setValue(int(search["samples_per_peak"]))
        if search.get("pdm_bins") is not None:
            self.pdm_bins.setValue(int(search["pdm_bins"]))
        methods = {str(value).lower() for value in search.get("methods", [])}
        if methods:
            self.chk_ls.setChecked("ls" in methods)
            self.chk_pdm.setChecked("pdm" in methods)
            self.chk_bls.setChecked("bls" in methods)

        requested_source = Path(str(handoff.get("source_file", "") or "")).name.lower()
        source_matches = (
            not requested_source
            or str(selected_item.get("source", "")).lower() == requested_source
        )
        requested_target = handoff.get("target_id")
        selected_target = selected_item.get("target_id")
        target_matches = (
            requested_target in (None, 0, "")
            or (selected_target is not None and int(selected_target) == int(requested_target))
        )
        mag_matches = _handoff_mag_columns_match(
            str(handoff.get("mag_col", "")),
            str(selected_item.get("mag_col", "")),
            str(handoff.get("correction_mode", "")),
        )

        signature = dict(handoff.get("input_signature") or {})
        loaded_time = np.asarray(self.lc_data.get("time", []), dtype=float)
        loaded_mag = np.asarray(self.lc_data.get("mag", []), dtype=float)
        loaded_err_raw = self.lc_data.get("mag_err")
        loaded_valid = np.isfinite(loaded_time) & np.isfinite(loaded_mag)
        if loaded_err_raw is not None:
            loaded_err = np.asarray(loaded_err_raw, dtype=float)
            loaded_valid &= np.isfinite(loaded_err) & (loaded_err > 0)
        signature_matches = True
        if signature:
            signature_matches &= int(signature.get("n_points", -1)) == int(np.count_nonzero(loaded_valid))
            if np.any(loaded_valid):
                for key, actual in (
                    ("time_min", float(np.nanmin(loaded_time[loaded_valid]))),
                    ("time_max", float(np.nanmax(loaded_time[loaded_valid]))),
                ):
                    expected = signature.get(key)
                    if expected is not None:
                        signature_matches &= bool(np.isclose(float(expected), actual, rtol=0.0, atol=1e-9))

        scan_results = {
            str(key): dict(value)
            for key, value in dict(handoff.get("scan_results") or {}).items()
            if isinstance(value, dict)
        }
        reuse_scan = bool(scan_results) and all(
            (source_matches, target_matches, mag_matches, signature_matches)
        )
        adopted_period = float(handoff.get("adopted_period", np.nan))
        if np.isfinite(adopted_period) and adopted_period > 0:
            self.center_p.setValue(adopted_period)
            self.phase_p.setValue(adopted_period)
            self.oc_p.setValue(adopted_period)

        self.handoff_multimode_diagnostic = dict(
            handoff.get("multimode_diagnostic") or {}
        )
        if not reuse_scan:
            self.scan_result = None
            self.scan_alias_analysis = None
            self.recommended_mode = "unknown"
            self.recommendation_text = (
                "The workspace series does not exactly match the Step 12 input. "
                "The period is loaded as a candidate; run Scan before routing."
            )
            self.workflow_step = "scan"
            self.scan_status.setText("Step 12 candidate loaded - rescan required")
            self.btn_refine.setEnabled(False)
            self._set_tool_tab("Periodogram")
            self._refresh_tool_workflow_ui()
            self._update_phase_plot()
            self.log(
                "[HANDOFF] Candidate only; scan reuse rejected "
                f"(source={source_matches}, target={target_matches}, "
                f"mag={mag_matches}, signature={signature_matches})."
            )
            self.log(f"[MAIN QC] release={bundle.release_status}: {bundle.release_message}")
            self._configure_automated_handoff(bundle)
            self.auto_status_label.setText(
                "The displayed workspace series differs from the release. "
                "Automation will validate and read the released source directly."
            )
            return

        self.scan_result = scan_results
        self.scan_alias_analysis = dict(handoff.get("alias_analysis") or {})
        self.alias_override = False
        if not (np.isfinite(adopted_period) and adopted_period > 0):
            adopted_period, _, _ = self._best_from_results(scan_results)

        self.recommended_mode, workflow_step, self.recommendation_text = _route_period_handoff(
            scan_results,
            self.scan_alias_analysis,
            self.handoff_multimode_diagnostic,
        )
        self.analysis_mode = (
            self.recommended_mode if workflow_step in {"single", "multi"} else "auto"
        )
        self.workflow_step = workflow_step

        diagnostic_status = str(
            self.handoff_multimode_diagnostic.get("status", "")
        ).upper()
        if diagnostic_status == "MULTIMODE-SUSPECT":
            candidate_period = float(
                self.handoff_multimode_diagnostic.get("candidate_period", np.nan)
            )
            mode_periods = [adopted_period]
            if (
                np.isfinite(candidate_period)
                and candidate_period > 0
                and abs(candidate_period - adopted_period) / adopted_period > 1e-5
            ):
                mode_periods.append(candidate_period)
            self.mm_periods_edit.setText(_format_period_list(mode_periods))

        alias_status = str(self.scan_alias_analysis.get("status", "UNASSESSED")).upper()
        self.scan_status.setText(
            f"Step 12 P = {adopted_period:.8f} d | alias={alias_status} | "
            f"mode={diagnostic_status or 'UNASSESSED'}"
        )
        # numpy 스칼라가 섞이면 `and` 가 numpy.bool_ 을 그대로 돌려주고 PyQt5 의
        # setEnabled 는 그것을 거부한다(cluster_structure 에서 실제로 창이 죽었다).
        self.btn_refine.setEnabled(bool(np.isfinite(adopted_period) and adopted_period > 0))
        self._draw_periodogram(scan_results)
        self._update_phase_plot()
        if workflow_step == "single":
            self._set_tool_tab("Refine")
        elif workflow_step == "multi":
            self._set_tool_tab("Multi-Mode")
        else:
            self._set_tool_tab("Periodogram")
        self._refresh_tool_workflow_ui()
        self.lc_status.setText(
            self.lc_status.text()
            + f"\nStep 12 handoff: P={adopted_period:.8f} d, {alias_status}"
            + f"\nMain QC release: {bundle.release_status}"
        )
        self.log(
            f"[HANDOFF] Reused Step 12 scan: P={adopted_period:.8f} d, "
            f"alias={alias_status}, route={workflow_step}."
        )
        self.log(f"[MAIN QC] release={bundle.release_status}: {bundle.release_message}")
        self._configure_automated_handoff(bundle)
        if bundle.release_status in {"APPROVED", "OVERRIDDEN"}:
            QTimer.singleShot(0, self._run_automated_analysis)

    def _clear_loaded_workspace_state(self, status: str):
        self.lc_data = None
        self.handoff_bundle = None
        self.automated_result = None
        self.series_options = {}
        self.scan_result = None
        self.scan_alias_analysis = None
        self.alias_override = False
        self.analysis_mode = "auto"
        self.handoff_multimode_diagnostic = None
        self.refined_period = None
        self.sigma_period = None
        self.multimode_result = None
        self.mm_candidate_scan = None
        self.mm_mode_rows = []
        self.mm_history = []
        self.recommended_mode = "unknown"
        self.recommendation_text = str(status)
        self.workflow_step = "load"
        self.mag_col_combo.blockSignals(True)
        self.mag_col_combo.clear()
        self.mag_col_combo.setEnabled(False)
        self.mag_col_combo.blockSignals(False)
        self.analysis_filter_combo.blockSignals(True)
        self.analysis_filter_combo.clear()
        self.analysis_filter_combo.addItem("All", "__all__")
        self.analysis_filter_combo.setEnabled(False)
        self.analysis_filter_combo.blockSignals(False)
        self.lc_status.setText(status)
        _set_role(self.lc_status, "status", "error")
        if hasattr(self, "auto_release_label"):
            self.auto_release_label.setText(
                "No Main-workflow release is loaded. Standalone data is unvalidated."
            )
            self._set_release_banner("neutral")
            self.auto_status_label.setText("Load a validated Step 12 release to run.")
            self.auto_progress.setValue(0)
            self.btn_auto_run.setEnabled(False)
            self.btn_auto_cancel.setEnabled(False)
            self.btn_auto_export.setEnabled(False)
            self.auto_candidate_row.hide()
            self.auto_advanced_chk.setChecked(True)
        self._clear_multimode_result(clear_inputs=False)
        for canvas_name in ("pg_canvas", "ref_canvas", "boot_canvas", "mm_pw_canvas", "mm_canvas", "ph_canvas", "oc_canvas", "fourier_canvas"):
            canvas = getattr(self, canvas_name, None)
            if canvas is None:
                continue
            fig = getattr(canvas, "figure", None)
            if fig is None:
                continue
            fig.clear()
            canvas.draw_idle()
        self._refresh_tool_workflow_ui()

    def closeEvent(self, event) -> None:
        self._pending_workspace_load = None
        for worker in (
            self._scan_worker,
            self._refine_worker,
            self._analysis_worker,
            self._workspace_load_worker,
        ):
            if worker is not None and worker.isRunning() and hasattr(worker, "stop"):
                worker.stop()
        super().closeEvent(event)

    def _apply_loaded_series_items(
        self, series_items: list[dict], workspace: Path
    ) -> None:
        if not series_items:
            self._clear_loaded_workspace_state("No usable light curve series found")
            return
        self.series_options = {item["key"]: item for item in series_items}
        preferred_label = ""
        try:
            from apex.utils.step_paths_lc import load_detrend_preference

            preferred_mode = load_detrend_preference(workspace)
            if preferred_mode:
                preferred_label = _CORR_MODE_LABELS.get(preferred_mode, "")
        except Exception:
            pass
        self.mag_col_combo.blockSignals(True)
        self.mag_col_combo.clear()
        for item in series_items:
            marker = (
                " *"
                if preferred_label and item["corr_tag"] == preferred_label
                else ""
            )
            self.mag_col_combo.addItem(
                item["combo_label"] + marker,
                item["key"],
            )
        default_index = 0
        if preferred_label:
            default_index = next(
                (
                    index
                    for index, item in enumerate(series_items)
                    if item["corr_tag"] == preferred_label
                    and "corr" in item["mag_col"]
                ),
                0,
            )
        if default_index == 0:
            default_index = next(
                (
                    index
                    for index, item in enumerate(series_items)
                    if item["corr_tag"] == "Global ensemble"
                    and "corr" in item["mag_col"]
                ),
                0,
            )
        self.mag_col_combo.setCurrentIndex(default_index)
        self.mag_col_combo.setEnabled(True)
        self.mag_col_combo.blockSignals(False)
        self._apply_series_option(self.mag_col_combo.currentData())

    def _load_paths(self, paths: list[Path]):
        try:
            rd = self._current_workspace_dir()
            try:
                from apex.utils.qc_utils import load_frame_excludes as _lfe
                excl = set(_lfe(rd).keys())
            except Exception:
                excl = set()

            series_items: list[dict] = []
            for path in paths:
                if not _is_primary_lightcurve_path(path):
                    continue
                df = pd.read_csv(path)
                if excl and "file" in df.columns:
                    df = df[~df["file"].astype(str).isin(excl)].reset_index(drop=True)
                time_col = next(
                    (c for c in ["BJD_TDB", "BJD", "bjd", "HJD", "hjd", "JD", "jd", "time"] if c in df.columns),
                    None
                )
                if time_col is None:
                    continue

                t = pd.to_numeric(df[time_col], errors="coerce").to_numpy(float)
                time_mask = np.isfinite(t)
                if not np.any(time_mask):
                    continue
                t = t[time_mask]

                err_col = next(
                    (c for c in ["diff_err_corr", "diff_err", "mag_err", "err", "sigma"] if c in df.columns),
                    None
                )
                e = pd.to_numeric(df[err_col], errors="coerce").to_numpy(float)[time_mask] if err_col else None

                filter_col = next(
                    (c for c in ["filter", "Filter", "FILTER", "band", "Band"] if c in df.columns), None
                )
                filters = df[filter_col].astype(str).to_numpy()[time_mask] if filter_col else None
                night_col = next(
                    (c for c in ["night_id", "night", "date"] if c in df.columns),
                    None,
                )
                if night_col:
                    night_ids = df[night_col].astype(str).to_numpy()[time_mask]
                else:
                    night_ids = infer_night_ids(t)
                corr_tag = _detect_corr_mode_from_df(df, path.name)
                target_id = _detect_target_id_from_df(df, path.name)

                for label, col, arr in _collect_mag_options(df, time_mask, corr_tag=corr_tag):
                    key = f"{path.name}::{col}"
                    series_items.append({
                        "key": key,
                        "time": t,
                        "mag": arr,
                        "mag_col": col,
                        "mag_err": e,
                        "filters": filters,
                        "night_id": night_ids,
                        "source": path.name,
                        "corr_tag": corr_tag,
                        "correction_preserves_nightly_baseline": (
                            col == "diff_mag_raw"
                            or corr_tag not in {"Nightly offset", "Color-dependent"}
                        ),
                        "series_label": _describe_series(corr_tag, col),
                        "target_id": target_id,
                    })

            if not series_items:
                self._clear_loaded_workspace_state("No usable light curve series found")
                return

            multi_target = len({item["target_id"] for item in series_items if item.get("target_id") is not None}) > 1
            for item in series_items:
                if multi_target:
                    tid = item.get("target_id")
                    item["combo_label"] = f"ID{tid} | {item['series_label']}" if tid is not None else f"{item['source']} | {item['series_label']}"
                else:
                    item["combo_label"] = item["series_label"]

            unique_series: dict[str, dict] = {}
            for item in series_items:
                label = item["combo_label"]
                prev = unique_series.get(label)
                if prev is None or _source_priority(item["source"]) < _source_priority(prev["source"]):
                    unique_series[label] = item
            series_items = list(unique_series.values())
            series_items.sort(key=lambda item: _series_rank(item["corr_tag"], item["mag_col"], item["source"]))
            self.series_options = {item["key"]: item for item in series_items}
            # Read step11 preference
            pref_label = ""
            try:
                from apex.utils.step_paths_lc import load_detrend_preference
                pref_mode = load_detrend_preference(rd)
                if pref_mode:
                    pref_label = _CORR_MODE_LABELS.get(pref_mode, "")
            except Exception:
                pass
            self.mag_col_combo.blockSignals(True)
            self.mag_col_combo.clear()
            for item in series_items:
                star = " *" if (pref_label and item["corr_tag"] == pref_label) else ""
                self.mag_col_combo.addItem(item["combo_label"] + star, item["key"])
            default_idx = 0
            if pref_label:
                default_idx = next(
                    (i for i, item in enumerate(series_items)
                     if item["corr_tag"] == pref_label and "corr" in item["mag_col"]), 0
                )
            if default_idx == 0:
                default_idx = next(
                    (i for i, item in enumerate(series_items)
                     if item["corr_tag"] == "Global ensemble" and "corr" in item["mag_col"]), 0
                )
            self.mag_col_combo.setCurrentIndex(default_idx)
            self.mag_col_combo.setEnabled(True)
            self.mag_col_combo.blockSignals(False)
            self._apply_series_option(self.mag_col_combo.currentData())
        except Exception as e:
            self._clear_loaded_workspace_state(f"Error: {e}")
            self.log(f"[ERROR] {e}")

    def _apply_series_option(self, key: str | None):
        if not key or key not in self.series_options:
            return
        item = self.series_options[key]
        selected_filter = self._refresh_analysis_filter_combo(item.get("filters"))
        t = item["time"]
        mag = item["mag"]
        mag_err = item.get("mag_err")
        filters = item.get("filters")
        night_id = item.get("night_id")
        if selected_filter and selected_filter != "__all__" and filters is not None:
            mask = (filters == selected_filter)
            t = t[mask]
            mag = mag[mask]
            mag_err = mag_err[mask] if mag_err is not None else None
            filters = filters[mask]
            night_id = night_id[mask] if night_id is not None else None
        elif selected_filter == "__all__" and filters is not None:
            mag = median_align_by_filter(mag, filters)
        self.lc_data = {
            "time": t,
            "mag": mag,
            "mag_col": item["mag_col"],
            "mag_err": mag_err,
            "filters": filters,
            "night_id": night_id,
            "source": item["source"],
            "corr_tag": item.get("corr_tag", ""),
            "correction_preserves_nightly_baseline": item.get(
                "correction_preserves_nightly_baseline", True
            ),
            "analysis_filter": selected_filter,
            "series_label": item.get("series_label", item["mag_col"]),
            "target_id": item.get("target_id"),
        }
        n = int(np.sum(np.isfinite(self.lc_data["time"]) & np.isfinite(self.lc_data["mag"])))
        corr_line = f"  [{self.lc_data['corr_tag']}]" if self.lc_data.get("corr_tag") else ""
        workspace_dir = self._current_workspace_dir()
        workspace_name = workspace_dir.name
        workspace_type = ""
        try:
            from apex.utils.run_workspace import load_run_manifest
            run_meta = load_run_manifest(workspace_dir)
            run_type = str(run_meta.get("run_type") or "").strip().lower()
            if run_type:
                workspace_type = f" [{run_type}]"
        except Exception:
            pass
        self.lc_status.setText(
            f"{workspace_name}{workspace_type}\n{self.lc_data['source']}\n{n} pts{corr_line}\n{self.mag_col_combo.currentText()}"
        )
        _set_role(self.lc_status, "status", "ok")
        filt_label = self.lc_data.get("analysis_filter", "__all__")
        filt_info = f", filter={filt_label}" if filt_label and filt_label != "__all__" else ""
        self.log(
            f"Loaded: {self.lc_data['source']}  ({n} pts, {self.lc_data.get('series_label', self.lc_data['mag_col'])}"
            f"{filt_info}, detrend={self.lc_data.get('corr_tag') or 'N/A'})"
        )
        if not self.lc_data.get("correction_preserves_nightly_baseline", True):
            self.log(
                "[WARNING] This corrected series removes nightly target baselines; "
                "use the matching raw series for long-period or multi-night mode searches."
            )
        t0_guess = float(np.nanmin(self.lc_data["time"]))
        self.t0_edit.setValue(t0_guess)
        self.oc_t0.setValue(t0_guess)
        if hasattr(self, "mm_state_epoch"):
            self.mm_state_epoch.setValue(t0_guess)
        self._update_fourier_filter_combo()
        self._clear_multimode_result(clear_inputs=False)
        self.scan_result = None
        self.scan_alias_analysis = None
        self.alias_override = False
        self.analysis_mode = "auto"
        self.refined_period = None
        self.sigma_period = None
        self.recommended_mode = "unknown"
        self.recommendation_text = "Run a period scan, then choose the Single or Multi path."
        self.workflow_step = "scan"
        self._refresh_tool_workflow_ui()

    def _refresh_analysis_filter_combo(self, filters) -> str:
        current = self.analysis_filter_combo.currentData()
        filter_values = filters.tolist() if filters is not None else []
        unique_filters = sorted({str(f) for f in filter_values if str(f).strip() and str(f).lower() != "nan"})
        self.analysis_filter_combo.blockSignals(True)
        self.analysis_filter_combo.clear()
        if unique_filters:
            for f in unique_filters:
                self.analysis_filter_combo.addItem(f, f)
            self.analysis_filter_combo.addItem("All", "__all__")
            target = current if current in unique_filters or current == "__all__" else unique_filters[0]
        else:
            self.analysis_filter_combo.addItem("All", "__all__")
            target = "__all__"
        idx = max(self.analysis_filter_combo.findData(target), 0)
        self.analysis_filter_combo.setCurrentIndex(idx)
        self.analysis_filter_combo.setEnabled(self.analysis_filter_combo.count() > 0)
        self.analysis_filter_combo.blockSignals(False)
        return self.analysis_filter_combo.currentData()

    def _update_fourier_filter_combo(self):
        current = self.fourier_filter_combo.currentData()
        self.fourier_filter_combo.blockSignals(True)
        self.fourier_filter_combo.clear()
        preferred = None
        if self.lc_data is not None:
            filters = self.lc_data.get("filters")
            if filters is not None:
                filter_values = np.asarray(filters, dtype=str)
                unique_filters = sorted(
                    {f for f in filter_values if f.strip() and f.lower() != "nan"}
                )
                for filt in unique_filters:
                    self.fourier_filter_combo.addItem(filt, filt)
                analysis_filter = self.lc_data.get("analysis_filter")
                if analysis_filter in unique_filters:
                    preferred = analysis_filter
                elif unique_filters:
                    preferred = max(
                        unique_filters,
                        key=lambda filt: int(np.count_nonzero(filter_values == filt)),
                    )
        if self.fourier_filter_combo.count() == 0:
            self.fourier_filter_combo.addItem("All data", "__all__")
            preferred = "__all__"
        target = current if self.fourier_filter_combo.findData(current) >= 0 else preferred
        index = self.fourier_filter_combo.findData(target)
        self.fourier_filter_combo.setCurrentIndex(max(index, 0))
        self.fourier_filter_combo.blockSignals(False)

    def _on_mag_col_changed(self):
        key = self.mag_col_combo.currentData()
        if not key:
            return
        self._apply_series_option(key)
        self._update_phase_plot()

    def _on_analysis_filter_changed(self):
        key = self.mag_col_combo.currentData()
        if not key:
            return
        self._apply_series_option(key)
        self._update_phase_plot()

    # ------------------------------------------------------------------
    # Scan
    # ------------------------------------------------------------------

    def _run_scan(self):
        if self.lc_data is None:
            QMessageBox.warning(self, "No Data", "Load a light curve first.")
            return
        if self._scan_worker and self._scan_worker.isRunning():
            return
        methods = []
        if self.chk_ls.isChecked(): methods.append("ls")
        if self.chk_pdm.isChecked(): methods.append("pdm")
        if self.chk_bls.isChecked(): methods.append("bls")
        if not methods:
            QMessageBox.warning(self, "No Method", "메서드를 하나 이상 선택하세요.")
            return

        self.scan_status.setText("Running…")
        self.workflow_step = "scan"
        self._refresh_tool_workflow_ui()
        mag = self.lc_data["mag"]
        self._scan_worker = PeriodAnalysisWorker(
            time=self.lc_data["time"],
            mag_raw=mag,
            mag_corr=None,
            mag_err=self.lc_data.get("mag_err"),
            min_period=self.min_p.value(),
            max_period=self.max_p.value(),
            samples_per_peak=self.spp.value(),
            methods=methods,
            pdm_n_bins=self.pdm_bins.value(),
            night_id=self.lc_data.get("night_id"),
            include_alias_diagnostics=True,
        )
        self._scan_worker.progress.connect(self.scan_status.setText)
        self._scan_worker.finished.connect(self._on_scan_done)
        self._scan_worker.error.connect(lambda e: (self.scan_status.setText("Error"), self.log(e)))
        self._scan_worker.start()

    def _on_scan_done(self, results: dict):
        # results keyed as "raw_ls", "raw_pdm", "raw_bls"
        self.scan_alias_analysis = results.pop("alias_analysis", None)
        self.alias_override = False
        self.analysis_mode = "auto"
        results.pop("multinight", None)
        self.scan_result = results
        best_period, best_power, fap = self._best_from_results(results)
        adopted = float((self.scan_alias_analysis or {}).get("adopted_period", np.nan))
        if np.isfinite(adopted) and adopted > 0:
            best_period = adopted
        self.recommended_mode, self.recommendation_text = _recommend_analysis_mode_from_scan(
            results,
            self.scan_alias_analysis,
        )
        fap_str = f"{fap:.2e}" if np.isfinite(fap) else "—"
        methods_run = [k.split("_", 1)[1].upper() for k in results]
        alias_status = str((self.scan_alias_analysis or {}).get("status", "UNASSESSED"))
        self.scan_status.setText(
            f"Candidate P = {best_period:.6f} d  [{'/'.join(methods_run)}] | {alias_status}"
        )
        self.center_p.setValue(best_period)
        self.phase_p.setValue(best_period)
        self.oc_p.setValue(best_period)
        self.btn_refine.setEnabled(True)
        self._draw_periodogram(results)
        self._update_phase_plot()
        self.log(f"Scan done: P={best_period:.6f} d, power={best_power:.4f}, FAP={fap_str}")
        self.log(f"[ROUTE] Recommended path: {self.recommended_mode} — {self.recommendation_text}")
        self._refresh_tool_workflow_ui()

    def _best_from_results(self, results: dict):
        """Pick best period across all methods (prefer ls > pdm > bls)."""
        for key in ("raw_ls", "raw_pdm", "raw_bls"):
            d = results.get(key)
            if d and "error" not in d and np.isfinite(d.get("best_period", np.nan)):
                return float(d["best_period"]), float(d["best_power"]), float(d.get("fap", np.nan))
        return np.nan, np.nan, np.nan

    def _draw_periodogram(self, results: dict):
        method_labels = {"ls": "Lomb-Scargle", "pdm": "PDM (1-θ)", "bls": "BLS"}
        method_colors = {"ls": "#1E88E5", "pdm": "#E53935", "bls": "#FF9800"}
        y_labels = {"ls": "LS Power", "pdm": "1 - θ", "bls": "BLS Power"}

        n = len(results)
        fig = self.pg_canvas.figure
        fig.clear()
        if n == 0:
            self.pg_canvas.draw_idle()
            return
        axes = fig.subplots(1, n, squeeze=False)[0]

        for ax, (key, data) in zip(axes, results.items()):
            method = key.split("_", 1)[1]
            if "error" in data:
                ax.text(0.5, 0.5, data["error"], ha="center", va="center",
                        transform=ax.transAxes, fontsize=9)
                ax.set_title(method_labels.get(method, method))
                continue
            if "frequency" in data:
                periods = 1.0 / data["frequency"]
            else:
                periods = data["trial_periods"]
            power = data["power"]
            best = data["best_period"]
            color = method_colors.get(method, "#666")
            ax.plot(periods, power, color=color, lw=0.8, alpha=0.9)
            ax.axvline(best, color="red", ls="--", lw=1.5, label=f"P={best:.6f} d")
            ax.scatter([best], [data["best_power"]], color="red", s=50, zorder=5)
            for p in data.get("top_periods", [])[1:4]:
                ax.axvline(p, color="orange", ls=":", lw=0.8, alpha=0.6)
            for candidate in (self.scan_alias_analysis or {}).get("candidates", [])[1:6]:
                candidate_period = float(candidate.get("period", np.nan))
                if np.isfinite(candidate_period) and candidate_period > 0:
                    ax.axvline(
                        candidate_period,
                        color="#C77A12",
                        ls="-.",
                        lw=0.9,
                        alpha=0.75,
                    )
            ax.set_xscale("log")
            ax.set_xlabel("Period (days)")
            ax.set_ylabel(y_labels.get(method, "Power"))
            ax.set_title(f"{method_labels.get(method, method)}\nP={best:.6f} d")
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)

        _safe_tight_layout(fig)
        self.pg_canvas.draw_idle()

    # ------------------------------------------------------------------
    # Refine
    # ------------------------------------------------------------------

    def _set_center_from_scan(self):
        adopted = float((self.scan_alias_analysis or {}).get("adopted_period", np.nan))
        if np.isfinite(adopted) and adopted > 0:
            self.center_p.setValue(adopted)
            return
        if self.scan_result:
            best, _, _ = self._best_from_results(self.scan_result)
            if np.isfinite(best):
                self.center_p.setValue(best)

    def _run_refine(self):
        if self.lc_data is None:
            return
        if self._refine_worker and self._refine_worker.isRunning():
            return
        alias_status = str((self.scan_alias_analysis or {}).get("status", "UNASSESSED")).upper()
        if alias_status != "RESOLVED" and not self.alias_override:
            QMessageBox.warning(
                self,
                "Alias Review Required",
                "The period alias is not resolved. Return to Period Review and explicitly "
                "choose Single-period to refine this candidate as a manual assumption.",
            )
            return
        self.workflow_step = "single"
        self._refresh_tool_workflow_ui()
        self.btn_refine.setEnabled(False)
        self.refine_status.setText("Running…")
        self.tabs.setCurrentIndex(1)  # Refine tab
        refine_method = self.refine_method_combo.currentData() or "ls"
        self._refine_worker = RefineBootstrapWorker(
            self.lc_data["time"], self.lc_data["mag"], self.lc_data.get("mag_err"),
            center_period=self.center_p.value(),
            n_bootstrap=self.n_boot.value(),
            method=refine_method,
        )
        self._refine_worker.progress.connect(self.refine_status.setText)
        self._refine_worker.finished.connect(self._on_refine_done)
        self._refine_worker.error.connect(self._on_refine_error)
        self._refine_worker.start()

    def _on_refine_done(self, result: dict):
        self.btn_refine.setEnabled(True)
        p = result["refined_period"]
        sig = result["sigma_p"]
        method_tag = result.get("method", "ls").upper()
        self.refined_period = p
        self.sigma_period = sig
        self.refine_status.setText(f"[{method_tag}] P = {p:.8f} | σ_local = {sig:.2e} d")
        self.refine_label.setText(
            f"P = {p:.8f} d | σ_local = {sig:.2e} d | "
            f"fixed-window residual bootstrap, N={self.n_boot.value()} | selected alias only"
        )
        self.phase_p.setValue(p)
        self.oc_p.setValue(p)
        self._draw_refine_plots(result)
        self._update_phase_plot()
        self.log(
            f"Refined selected alias: P={p:.8f}, σ_local={sig:.2e} d "
            "(does not include alias-family uncertainty)"
        )

    def _on_refine_error(self, msg: str):
        self.btn_refine.setEnabled(True)
        self.refine_status.setText("Error")
        QMessageBox.warning(self, "Refine Error", msg)

    def _draw_refine_plots(self, result: dict):
        fp = np.asarray(result["fine_periods"], dtype=float)
        pw = np.asarray(result["fine_power"], dtype=float)
        p_best = result["refined_period"]
        sig = result["sigma_p"]

        fig1 = self.ref_canvas.figure
        fig1.clear()
        ax1 = fig1.add_subplot(111)
        if len(fp) and len(pw):
            ax1.plot(fp, pw, color="#7B1FA2", lw=0.8)
        ax1.axvline(p_best, color="red", ls="--", lw=1.5, label=f"P={p_best:.8f} d")
        ax1.set_xlabel("Period (days)")
        ax1.set_ylabel("Relative fit power")
        ax1.set_title("Per-filter Multi-harmonic Period Grid")
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3)
        _safe_tight_layout(fig1)
        self.ref_canvas.draw_idle()

        bp = np.asarray(result["boot_periods"], dtype=float)
        fig2 = self.boot_canvas.figure
        fig2.clear()
        ax2 = fig2.add_subplot(111)
        if len(bp):
            ax2.hist(
                bp,
                bins=min(40, len(bp) // 3 + 1),
                color="#7B1FA2",
                alpha=0.7,
                edgecolor="white",
            )
        else:
            ax2.text(
                0.5,
                0.5,
                "Bootstrap disabled",
                ha="center",
                va="center",
                transform=ax2.transAxes,
            )
        ax2.axvline(p_best, color="red", ls="--", lw=2, label=f"P={p_best:.8f}")
        ax2.axvline(p_best - sig, color="orange", ls=":", lw=1.5)
        ax2.axvline(p_best + sig, color="orange", ls=":", lw=1.5, label=f"σ_local={sig:.2e}")
        ax2.set_xlabel("Bootstrap Period (d)")
        ax2.set_ylabel("Count")
        ax2.set_title(f"Local Residual Bootstrap (N={len(bp)}; selected alias)")
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)
        _safe_tight_layout(fig2)
        self.boot_canvas.draw_idle()

    # ------------------------------------------------------------------
    # Phase plot
    # ------------------------------------------------------------------

    def _on_epoch_kind_changed(self):
        event_kind = self.epoch_kind_combo.currentData() or "max_light"
        if event_kind == "min_light":
            button_text = "Detect Minimum-light Epoch"
            event_text = "Minimum light"
        else:
            button_text = "Detect Maximum-light Epoch"
            event_text = "Maximum light"
        self.btn_detect_t0.setText(button_text)
        self.oc_event_label.setText(event_text)

    def _detect_t0(self):
        if self.lc_data is None:
            return
        t = self.lc_data["time"]
        mag = self.lc_data["mag"]
        period = self.phase_p.value()
        mask = np.isfinite(t) & np.isfinite(mag)
        t_v, m_v = t[mask], mag[mask]
        phase_view_mode = self.phase_data_mode.currentData()
        if self._effective_analysis_mode() == "multi":
            if phase_view_mode != "focused":
                QMessageBox.warning(
                    self,
                    "Epoch Detection",
                    "A composite multi-mode light curve has no stable single event epoch. "
                    "Fit the modes and select Focused mode isolation first.",
                )
                return
            night_ids = self.lc_data.get("night_id")
            night_v = night_ids[mask] if night_ids is not None else None
            focus_payload = self._phase_multimode_payload(
                t_v,
                m_v,
                float(period),
                night_id=night_v,
            )
            if focus_payload is None:
                QMessageBox.warning(
                    self,
                    "Epoch Detection",
                    "Select a fitted mode before detecting its event epoch.",
                )
                return
            m_v = np.asarray(focus_payload["iso_mag"], dtype=float)

        event_kind = self.epoch_kind_combo.currentData() or "max_light"
        try:
            t0_epoch, event_phase = _detect_phase_epoch(
                t_v,
                m_v,
                float(period),
                event_kind=str(event_kind),
            )
        except ValueError as exc:
            QMessageBox.warning(self, "Epoch Detection", str(exc))
            return
        self.t0_edit.setValue(t0_epoch)
        self.oc_t0.setValue(t0_epoch)
        event_label = "maximum light" if event_kind == "max_light" else "minimum light"
        self.log(
            f"T₀ detected ({event_label}): {t0_epoch:.6f} BJD  "
            f"(phase={event_phase:.4f})"
        )

    def _overlay_phase_fourier_fit(
        self,
        ax,
        t_v: np.ndarray,
        m_v: np.ndarray,
        filt_v,
        dy_v: np.ndarray | None = None,
    ):
        if not self.phase_fit_chk.isChecked():
            return
        period = float(self.phase_p.value())
        t0 = float(self.t0_edit.value())
        if not np.isfinite(period) or period <= 0:
            return
        harmonics = int(self.phase_fit_harm.value())
        phase_model = np.linspace(0.0, 2.0, 1000)
        time_model = t0 + (phase_model % 1.0) * period

        if filt_v is not None:
            unique_filters = sorted(set(filt_v))
        else:
            unique_filters = [""]

        for fi, filt in enumerate(unique_filters):
            if filt and not self.filter_visibility.get(filt, True):
                continue
            sel = (filt_v == filt) if filt_v is not None else np.ones(len(t_v), dtype=bool)
            if np.sum(sel) < max(8, 2 * harmonics + 3):
                continue
            service_fit = None
            if self.automated_result is not None and self.automated_result.branch == "single":
                candidate = self.automated_result.per_filter_models.get(str(filt))
                if candidate and "error" not in candidate:
                    fitted_period = float(candidate.get("period", np.nan))
                    if (
                        np.isfinite(fitted_period)
                        and abs(fitted_period - period) / period <= 1e-5
                    ):
                        service_fit = candidate
            if service_fit is not None:
                model_mag = _evaluate_service_fourier(time_model, service_fit)
                color = self.filter_colors.get(filt) or _filt_color(filt, fi)
                label = f"{filt} validated fit" if filt else "Validated fit"
                ax.plot(
                    phase_model,
                    model_mag,
                    color=color,
                    lw=1.6,
                    alpha=0.95,
                    zorder=4,
                    label=label,
                )
                continue
            try:
                fit_result = _fit_fixed_period_fourier(
                    t_v[sel],
                    m_v[sel],
                    period,
                    harmonics,
                    mag_err=dy_v[sel] if dy_v is not None else None,
                )
            except Exception:
                continue
            model_mag = _evaluate_fixed_period_fourier(time_model, period, fit_result)
            color = self.filter_colors.get(filt) or _filt_color(filt, fi)
            label = f"{filt} fit" if filt else f"Fit (N={harmonics})"
            ax.plot(phase_model, model_mag, color=color, lw=1.6, alpha=0.95, zorder=4, label=label)

    def _overlay_multimode_phase_plot(self, ax, t_v: np.ndarray, m_v: np.ndarray, phase: np.ndarray):
        if not self.multimode_result or not self.mm_overlay_chk.isChecked():
            return
        focus_idx = self._selected_multimode_focus_index()
        if focus_idx is None:
            return
        focus_period = float(self.multimode_result["periods"][focus_idx])
        current_period = float(self.phase_p.value())
        if not np.isfinite(current_period) or current_period <= 0:
            return
        if abs(current_period - focus_period) / max(abs(focus_period), 1e-12) > 0.02:
            return

        eval_obs = _evaluate_multimode_result(self.multimode_result, t_v)
        focus_component = np.asarray(eval_obs["components"][focus_idx], dtype=float)
        iso_mag = m_v - (eval_obs["total"] - (float(eval_obs["intercept"]) + focus_component))
        phase_ext = np.concatenate([phase, phase + 1.0])
        iso_ext = np.concatenate([iso_mag, iso_mag])
        ax.scatter(
            phase_ext, iso_ext, s=14, facecolors="none", edgecolors="#263238",
            linewidths=0.8, alpha=0.55, zorder=4, label=f"M{focus_idx + 1} isolated"
        )

        t0 = float(self.t0_edit.value())
        phase_model = np.linspace(0.0, 2.0, 800)
        phase_times = t0 + (phase_model % 1.0) * focus_period
        eval_model = _evaluate_multimode_result(self.multimode_result, phase_times)
        focus_model = eval_model["intercept"] + eval_model["components"][focus_idx]
        ax.plot(phase_model, focus_model, color="#212121", lw=1.4, zorder=5, label=f"M{focus_idx + 1} fit")

    def _update_phase_plot(self):
        if self.lc_data is None:
            return
        period = self.phase_p.value()
        t0 = self.t0_edit.value()
        t = self.lc_data["time"]
        mag = self.lc_data["mag"]
        mag_err = self.lc_data.get("mag_err")
        filters = self.lc_data.get("filters")
        night_ids = self.lc_data.get("night_id")
        mask = np.isfinite(t) & np.isfinite(mag)
        t_v, m_v = t[mask], mag[mask]
        dy_v = np.asarray(mag_err, dtype=float)[mask] if mag_err is not None else None
        filt_v = filters[mask] if filters is not None else None
        night_v = night_ids[mask] if night_ids is not None else None

        phase = ((t_v - t0) / period) % 1.0
        plot_mag = m_v
        title_suffix = ""
        phase_view_mode = self.phase_data_mode.currentData() if hasattr(self, "phase_data_mode") else "raw"
        focus_payload = None
        if phase_view_mode == "focused":
            focus_payload = self._phase_multimode_payload(
                t_v,
                m_v,
                float(period),
                night_id=night_v,
            )
            if focus_payload is not None:
                plot_mag = np.asarray(focus_payload["iso_mag"], dtype=float)
                title_suffix = f"  [M{int(focus_payload['index']) + 1} isolated]"
            else:
                title_suffix = "  [select fitted mode]"
        elif self.multimode_result:
            title_suffix = "  [raw composite]"

        fig = self.ph_canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)

        unique_filters = sorted(set(filt_v)) if filt_v is not None else [""]
        for fi, filt in enumerate(unique_filters):
            if not self.filter_visibility.get(filt, True):
                continue
            sel = (filt_v == filt) if filt_v is not None else np.ones(len(phase), bool)
            ph_sel = np.concatenate([phase[sel], phase[sel] + 1.0])
            m_sel = np.concatenate([plot_mag[sel], plot_mag[sel]])
            color = self.filter_colors.get(filt) or _filt_color(filt, fi)
            label = filt if filt else "data"
            if dy_v is not None:
                dy_sel = np.concatenate([dy_v[sel], dy_v[sel]])
                ax.errorbar(ph_sel, m_sel, yerr=dy_sel,
                            fmt="o", color=color, markersize=3,
                            elinewidth=0.5, capsize=0, alpha=0.7, label=label)
            else:
                ax.scatter(ph_sel, m_sel, color=color, s=12, alpha=0.7, label=label)

        target_title = _target_plot_title(self.params, self.lc_data.get("target_id"))
        ax.invert_yaxis()
        ax.set_xlabel("Phase")
        ax.set_ylabel("Magnitude")
        ax.set_title(f"{target_title}  |  P={period:.6f} d   T₀={t0:.4f}{title_suffix}")
        ax.set_xlim(0, 2)
        ax.axvline(0, color="gray", ls=":", alpha=0.4)
        ax.axvline(1, color="gray", ls=":", alpha=0.4)
        ax.grid(True, alpha=0.3)
        legend_title = "Filter" if (len(unique_filters) > 1 or (unique_filters and unique_filters[0] != "")) else None

        # Check star phase-folded overlay
        if self.phase_check_overlay_chk.isChecked():
            try:
                _rd = self._current_workspace_dir()
                _check_filter = _resolve_check_filter(
                    self.lc_data.get("filters"),
                    self.lc_data.get("analysis_filter"),
                )
                _ck_id, _ck_df = _load_check_star_for_plot(_rd, filt=_check_filter)
                if _ck_df is not None and not _ck_df.empty:
                    _t_col, _y_col = _pick_check_overlay_cols(_ck_df, self.lc_data.get("mag_col"))
                    if _t_col and _y_col:
                        _ct = pd.to_numeric(_ck_df[_t_col], errors="coerce").to_numpy(float)
                        _cm = pd.to_numeric(_ck_df[_y_col], errors="coerce").to_numpy(float)
                        _mask = np.isfinite(_ct) & np.isfinite(_cm)
                        if _mask.any():
                            _ck_label = f"Check ID {_ck_id}" if _ck_id is not None else "Check"
                            _ck_label += " (median shifted)"
                            _phase = ((_ct[_mask] - t0) / period) % 1.0
                            _phase_ext = np.concatenate([_phase, _phase + 1.0])
                            _cm_shifted = _cm[_mask] - np.nanmedian(_cm[_mask]) + np.nanmedian(plot_mag)
                            _mag_ext = np.concatenate([_cm_shifted, _cm_shifted])
                            ax.scatter(
                                _phase_ext,
                                _mag_ext,
                                s=8,
                                color="#FFD700",
                                alpha=0.4,
                                zorder=2,
                                label=_ck_label,
                                marker="^",
                            )
            except Exception:
                pass

        self._overlay_phase_fourier_fit(ax, t_v, plot_mag, filt_v, dy_v)
        if phase_view_mode != "focused":
            self._overlay_multimode_phase_plot(ax, t_v, m_v, phase)
        elif focus_payload is None and self.multimode_result:
            ax.text(
                0.5,
                0.03,
                "Select M1/M2 from 'Mode from fit' for a meaningful single-mode phase curve.",
                transform=ax.transAxes,
                ha="center",
                va="bottom",
                fontsize=9,
                color="#6D4C41",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="#FFF3E0", edgecolor="#FFCC80", alpha=0.9),
            )
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(fontsize=8, title=legend_title)

        _safe_tight_layout(fig)
        self.ph_canvas.draw_idle()

    # ------------------------------------------------------------------
    # O-C
    # ------------------------------------------------------------------

    def _oc_from_refine(self):
        if self.refined_period is not None:
            self.oc_p.setValue(self.refined_period)
        if self.t0_edit.value():
            self.oc_t0.setValue(self.t0_edit.value())

    def _oc_add(self):
        r = self.oc_table.rowCount()
        self.oc_table.insertRow(r)
        for c in range(4):
            self.oc_table.setItem(r, c, QTableWidgetItem(""))

    def _oc_del(self):
        rows = sorted({i.row() for i in self.oc_table.selectedIndexes()}, reverse=True)
        for r in rows:
            self.oc_table.removeRow(r)
        self._draw_oc()

    def _oc_import(self):
        path, _ = QFileDialog.getOpenFileName(self, "Import O-C CSV", "", "CSV (*.csv)")
        if not path:
            return
        try:
            df = pd.read_csv(path)
            bjd_col = next(
                (c for c in df.columns if c.lower() in ("bjd_obs", "bjd", "hjd", "jd")), None
            )
            if bjd_col is None:
                QMessageBox.warning(self, "Error", "BJD 컬럼을 찾을 수 없습니다.")
                return
            err_col = next(
                (c for c in df.columns if c.lower() in ("err", "error", "sigma")), None
            )
            epoch_col = next(
                (c for c in df.columns if c.lower() in ("n", "epoch", "cycle", "e")), None
            )
            t0 = self.oc_t0.value()
            p = self.oc_p.value()
            self.oc_table.setRowCount(0)
            for _, row_data in df.iterrows():
                bjd = float(row_data[bjd_col])
                epoch_value = pd.to_numeric(row_data[epoch_col], errors="coerce") if epoch_col else np.nan
                n = (
                    int(round(float(epoch_value)))
                    if np.isfinite(epoch_value)
                    else (int(round((bjd - t0) / p)) if p > 0 else 0)
                )
                oc = bjd - (t0 + n * p)
                err = f"{float(row_data[err_col]):.6f}" if err_col else ""
                r = self.oc_table.rowCount()
                self.oc_table.insertRow(r)
                self.oc_table.setItem(r, 0, QTableWidgetItem(str(n)))
                self.oc_table.setItem(r, 1, QTableWidgetItem(f"{bjd:.6f}"))
                self.oc_table.setItem(r, 2, QTableWidgetItem(f"{oc:.6f}"))
                self.oc_table.setItem(r, 3, QTableWidgetItem(err))
            self._draw_oc()
        except Exception as e:
            QMessageBox.warning(self, "Import Error", str(e))

    def _oc_export(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export O-C CSV", "", "CSV (*.csv)")
        if not path:
            return
        rows = []
        for r in range(self.oc_table.rowCount()):
            rows.append([
                self.oc_table.item(r, c).text() if self.oc_table.item(r, c) else ""
                for c in range(4)
            ])
        pd.DataFrame(rows, columns=["n", "BJD_obs", "O-C (d)", "err (d)"]).to_csv(path, index=False)
        self.log(f"Exported O-C to {path}")

    def _recompute_oc(self):
        t0 = self.oc_t0.value()
        p = self.oc_p.value()
        if p <= 0:
            return
        for r in range(self.oc_table.rowCount()):
            bjd_item = self.oc_table.item(r, 1)
            if not bjd_item or not bjd_item.text():
                continue
            try:
                bjd = float(bjd_item.text())
                n = int(round((bjd - t0) / p))
                oc = bjd - (t0 + n * p)
                self.oc_table.setItem(r, 0, QTableWidgetItem(str(n)))
                self.oc_table.setItem(r, 2, QTableWidgetItem(f"{oc:.6f}"))
            except ValueError:
                pass
        self._draw_oc()

    def _get_oc_arrays(self):
        ns, ocs, errs = [], [], []
        for r in range(self.oc_table.rowCount()):
            try:
                n_it = self.oc_table.item(r, 0)
                oc_it = self.oc_table.item(r, 2)
                er_it = self.oc_table.item(r, 3)
                if not n_it or not oc_it or not n_it.text() or not oc_it.text():
                    continue
                ns.append(int(n_it.text()))
                ocs.append(float(oc_it.text()))
                errs.append(float(er_it.text()) if (er_it and er_it.text()) else np.nan)
            except (ValueError, TypeError):
                pass
        return np.array(ns), np.array(ocs), np.array(errs)

    def _draw_oc(self, fit_result=None):
        fig = self.oc_canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)
        ns, ocs, errs = self._get_oc_arrays()
        if len(ns) == 0:
            ax.text(0.5, 0.5, "O-C 데이터 없음", ha="center", va="center",
                    transform=ax.transAxes, fontsize=12, color="gray")
            _safe_tight_layout(fig)
            self.oc_canvas.draw_idle()
            return
        oc_min = ocs * 1440
        has_err = np.any(np.isfinite(errs))
        if has_err:
            err_pl = np.where(np.isfinite(errs), errs * 1440, 0.0)
            ax.errorbar(ns, oc_min, yerr=err_pl, fmt="o", color="#1565C0",
                        markersize=5, elinewidth=1, capsize=3, label="O-C")
        else:
            ax.scatter(ns, oc_min, color="#1565C0", s=40, zorder=5, label="O-C")
        ax.axhline(0, color="gray", ls="--", lw=1, alpha=0.6)
        if fit_result is not None:
            n_fit = np.linspace(ns.min(), ns.max(), 500)
            ax.plot(n_fit, fit_result["func"](n_fit) * 1440, color="red",
                    lw=1.5, label=fit_result["label"])
        ax.set_xlabel("Epoch (n)")
        ax.set_ylabel("O-C (minutes)")
        ax.set_title("O-C Diagram")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        _safe_tight_layout(fig)
        self.oc_canvas.draw_idle()

    def _oc_fit(self):
        ns, ocs, errs = self._get_oc_arrays()
        mode = self.oc_fit_combo.currentData() or "none"
        if mode == "none":
            self.oc_fit_label.clear()
            self._draw_oc()
            return
        minimum_points = {
            "linear": 3,
            "quadratic": 5,
            "quadratic_sine": 12,
        }.get(str(mode), 3)
        if len(ns) < minimum_points or len(np.unique(ns)) < minimum_points:
            QMessageBox.warning(
                self,
                "O-C Fit",
                f"This model requires at least {minimum_points} distinct timing epochs.",
            )
            return
        p = self.oc_p.value()
        sigma = None
        if np.any(np.isfinite(errs)):
            fill = np.nanmedian(errs[np.isfinite(errs)])
            sigma = np.where(np.isfinite(errs), errs, fill)
            sigma = np.clip(sigma, 1e-9, None)
        fit_kwargs = {
            "sigma": sigma,
            "absolute_sigma": sigma is not None,
        }

        fit_result = None
        text = ""
        try:
            if mode == "linear":
                def model(n, a, b): return a + b * n
                popt, pcov = curve_fit(model, ns, ocs, **fit_kwargs)
                perr = np.sqrt(np.diag(pcov))
                dp = popt[1]
                text = (
                    f"O-C = {popt[0]*1440:.3f} + {dp*1440:.4f}·n  (min)\n"
                    f"ΔP = {dp:.3e} ± {perr[1]:.1e} d/cycle\n"
                    f"Corrected P = {p + dp:.8f} d"
                )
                fit_result = {"func": lambda n: model(n, *popt), "label": f"Linear (ΔP={dp:.2e})"}

            elif mode == "quadratic":
                def model(n, a, b, c): return a + b * n + c * n**2
                popt, pcov = curve_fit(model, ns, ocs, maxfev=10000, **fit_kwargs)
                perr = np.sqrt(np.diag(pcov))
                c = popt[2]
                dP_dt = 2 * c / p
                dP_yr = dP_dt * 365.25
                arrow = "▲ 주기 증가" if dP_dt > 0 else "▼ 주기 감소"
                text = (
                    f"c = {c:.3e} ± {perr[2]:.1e} d/cycle²\n"
                    f"dP/dt = {dP_dt:.3e} d/d  ({dP_yr:.3e} d/yr)\n"
                    f"{arrow}"
                )
                fit_result = {"func": lambda n: model(n, *popt), "label": f"Quadratic (dP/dt={dP_dt:.2e})"}

            elif mode == "quadratic_sine":
                n_span = ns.max() - ns.min()
                if n_span <= 10:
                    QMessageBox.warning(
                        self,
                        "O-C Fit",
                        "The timing baseline is too short for an exploratory sinusoidal term.",
                    )
                    return
                A_g = (ocs.max() - ocs.min()) / 2
                p3_g = max(10.0, n_span / 2.0)

                def model(n, a, b, c, A, p3_n, phi):
                    return a + b*n + c*n**2 + A * np.sin(2*np.pi*n/p3_n + phi)

                bounds = (
                    [-np.inf, -np.inf, -np.inf, 0, 10, -np.pi],
                    [np.inf, np.inf, np.inf, np.inf, n_span * 10, np.pi],
                )
                popt, pcov = curve_fit(
                    model, ns, ocs, p0=[0, 0, 0, A_g, p3_g, 0],
                    bounds=bounds, maxfev=30000, **fit_kwargs,
                )
                perr = np.sqrt(np.diag(pcov))
                c, A, p3_n = popt[2], abs(popt[3]), abs(popt[4])
                dP_dt = 2 * c / p
                p3_days = p3_n * p
                p3_err_days = float(perr[4] * p) if len(perr) > 4 else np.nan
                p3_yr = p3_days / 365.25
                text = (
                    f"dP/dt = {dP_dt:.3e} d/d\n"
                    f"Candidate modulation period = {p3_days:.1f} ± {p3_err_days:.1f} d "
                    f"({p3_yr:.2f} yr)\n"
                    f"O-C semi-amplitude = {A*1440:.2f} min\n"
                    "Exploratory model only; a physical third-body interpretation requires "
                    "independent validation and adequate cycle coverage."
                )
                fit_result = {
                    "func": lambda n: model(n, *popt),
                    "label": f"Exploratory sinusoid ({p3_yr:.2f} yr)",
                }

        except Exception as e:
            text = f"피팅 실패: {e}"

        self.oc_fit_label.setText(text)
        self._draw_oc(fit_result)
        mode_label = self.oc_fit_combo.currentText()
        self.log(f"[O-C Fit] {mode_label}: {text.splitlines()[0]}")

    # ------------------------------------------------------------------
    # Fourier decomposition
    # ------------------------------------------------------------------

    def _run_fourier(self):
        if self.lc_data is None:
            QMessageBox.warning(self, "No Data", "Load a light curve first.")
            return
        period = self.phase_p.value()
        if period <= 0:
            return

        t = self.lc_data["time"]
        mag = self.lc_data["mag"]
        filters = self.lc_data.get("filters")
        mag_err = self.lc_data.get("mag_err")
        t0 = self.t0_edit.value()

        sel_filter = self.fourier_filter_combo.currentData() or "__all__"
        if sel_filter != "__all__" and filters is not None:
            fmask = np.asarray(filters, dtype=str) == str(sel_filter)
        else:
            fmask = np.ones(len(t), bool)

        mask = np.isfinite(t) & np.isfinite(mag) & fmask
        t_v, m_v = t[mask], mag[mask]
        dy_v = np.asarray(mag_err, dtype=float)[mask] if mag_err is not None else None

        if len(t_v) < 20:
            QMessageBox.warning(self, "Fourier", "데이터가 너무 적습니다 (< 20).")
            return

        nh = self.n_harm.value()
        phase = ((t_v - t0) / period) % 1.0
        try:
            fit_result = _fit_fixed_period_fourier(
                t_v,
                m_v,
                float(period),
                int(nh),
                mag_err=dy_v,
            )
        except Exception as e:
            QMessageBox.warning(self, "Fourier", f"LSQ failed: {e}")
            return

        coeff = np.asarray(fit_result["coeff"], dtype=float)
        a0 = coeff[0]
        shape = _fourier_shape_parameters(coeff)
        amp_k = shape["amplitudes"]
        phi_k = shape["phases"]
        R21 = shape["r21"]
        phi21 = shape["phi21"]
        R31 = shape["r31"]
        phi31 = shape["phi31"]

        ph_model = np.linspace(0, 1, 1000)
        t_model = t0 + ph_model * period
        mag_model = _evaluate_fixed_period_fourier(t_model, period, fit_result)
        model_at_data = _evaluate_fixed_period_fourier(t_v, period, fit_result)

        amp_ptp = float(np.max(mag_model) - np.min(mag_model))
        obs_range = float(np.percentile(m_v, 95) - np.percentile(m_v, 5))
        residuals = m_v - model_at_data
        rms_res = float(np.std(residuals))

        filt_label = str(sel_filter) if sel_filter != "__all__" else "All data"
        fit_kind = "weighted" if fit_result.get("weighted") else "unweighted"
        text = f"[{filt_label} | {nh} harmonics | n={len(t_v)} | {fit_kind}]\n"
        text += f"t_ref = {fit_result['time_ref']:.6f} BJD\n"
        text += f"Model peak-to-peak : {amp_ptp:.4f} mag\n"
        text += f"Observed 5-95%     : {obs_range:.4f} mag\n"
        text += f"Residual RMS       : {rms_res:.4f} mag\n"
        text += "\nHarmonic coefficients (cosine-series convention)\n"
        text += f"A₀ = {a0:.4f} mag\n"
        for k in range(nh):
            text += f"A{k+1} = {amp_k[k]:.4f}   φ{k+1} = {phi_k[k]:.4f} rad\n"
        text += "\nShape parameters\n"
        text += f"R₂₁ = A₂/A₁     = {R21:.4f}\n"
        text += f"φ₂₁ = φ₂−2φ₁   = {phi21:.4f} rad\n"
        if np.isfinite(R31):
            text += f"R₃₁ = A₃/A₁     = {R31:.4f}\n"
            text += f"φ₃₁ = φ₃−3φ₁   = {phi31:.4f} rad\n"
        text += "\nR₂₁ and φ₂₁ describe morphology; do not assign a variable class from them alone."
        self.fourier_label.setText(text)

        # Plot
        fig = self.fourier_canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)

        color = (
            self.filter_colors.get(str(sel_filter)) or _filt_color(str(sel_filter), 0)
            if sel_filter != "__all__" else "#1E88E5"
        )
        ax.scatter(phase, m_v, color=color, s=8, alpha=0.6, label=filt_label)

        ax.plot(ph_model, mag_model, color="red", lw=1.5, label=f"Fourier (N={nh})")
        # Mark max/min of model
        ax.axhline(np.max(mag_model), color="gray", ls=":", lw=0.8, alpha=0.7)
        ax.axhline(np.min(mag_model), color="gray", ls=":", lw=0.8, alpha=0.7)
        ax.invert_yaxis()
        ax.set_xlabel("Phase")
        ax.set_ylabel("Magnitude")
        r21_str = f"{R21:.3f}" if np.isfinite(R21) else "—"
        ax.set_title(f"Fourier [{filt_label}]   Δm={amp_ptp:.3f}   R₂₁={r21_str}")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        _safe_tight_layout(fig)
        self.fourier_canvas.draw_idle()

        self.log(f"Fourier [{filt_label}] N={nh}: Δm={amp_ptp:.3f}, R21={R21:.4f}, φ21={phi21:.4f}")

    # ------------------------------------------------------------------
    # Filter color/visibility browser
    # ------------------------------------------------------------------

    def _apply_filter_swatch_style(self, button, color: str):
        # The background IS the data (the filter's plot colour) — the one
        # legitimate per-button background. Borders come from the theme.
        button.setStyleSheet(
            f"QPushButton {{ background-color: {color}; "
            f"border: 1px solid {Tokens.BORDER_STRONG}; "
            f"border-radius: 3px; min-width: 28px; min-height: 20px; }}"
            f"QPushButton:hover {{ border: 2px solid {Tokens.ACCENT}; }}"
        )

    def show_filter_color_browser(self):
        if self.lc_data is None:
            QMessageBox.information(self, "Filter Browser", "먼저 라이트커브를 로드하세요.")
            return
        filters = self.lc_data.get("filters")
        if filters is None:
            QMessageBox.information(self, "Filter Browser", "이 CSV에는 filter 컬럼이 없습니다.")
            return
        keys = sorted(set(filters))

        # Initialize missing entries
        for i, k in enumerate(keys):
            if k not in self.filter_colors:
                self.filter_colors[k] = _filt_color(k, i)
            if k not in self.filter_visibility:
                self.filter_visibility[k] = True

        dialog = FittedDialog(self)
        dialog.setWindowTitle("필터 색상 / 표시 설정")
        dialog.setMinimumWidth(360)
        layout = QVBoxLayout(dialog)

        info = QLabel("필터별 색상을 변경하거나 표시/숨김을 설정합니다.\n변경 즉시 위상 그래프에 반영됩니다.")
        info.setProperty("role", "caption")
        info.setWordWrap(True)
        layout.addWidget(info)

        grid = QGridLayout()
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(6)
        swatch_btns: dict = {}
        vis_chks: dict = {}

        for row, key in enumerate(keys):
            # Visibility checkbox
            chk = QCheckBox(key)
            chk.setChecked(self.filter_visibility.get(key, True))
            _f = chk.font(); _f.setBold(True); chk.setFont(_f)
            def _on_vis(state, k=key):
                self.filter_visibility[k] = bool(state)
                self._update_phase_plot()
            chk.stateChanged.connect(_on_vis)
            grid.addWidget(chk, row, 0)
            vis_chks[key] = chk

            # Color swatch
            swatch = QPushButton("")
            swatch.setFixedSize(28, 20)
            swatch.setFocusPolicy(Qt.NoFocus)
            self._apply_filter_swatch_style(swatch, self.filter_colors[key])
            grid.addWidget(swatch, row, 1)
            swatch_btns[key] = swatch

            # Browse button
            browse = QPushButton("Browse…")
            def _on_browse(_checked=False, k=key, sw=None, dlg=dialog):
                current = QColor(self.filter_colors.get(k, "#888888"))
                picked = QColorDialog.getColor(current, dlg, f"{k} 색상 선택")
                if picked.isValid():
                    self.filter_colors[k] = picked.name()
                    self._apply_filter_swatch_style(sw, picked.name())
                    self._update_phase_plot()
            # bind swatch via default arg
            browse.clicked.connect(lambda _c=False, k=key, sw=swatch: _on_browse(_c, k, sw, dialog))
            grid.addWidget(browse, row, 2)

        layout.addLayout(grid)

        btn_row = QHBoxLayout()
        btn_reset = QPushButton("색상 초기화")
        def _reset():
            for i, k in enumerate(keys):
                self.filter_colors[k] = _filt_color(k, i)
                self._apply_filter_swatch_style(swatch_btns[k], self.filter_colors[k])
            self._update_phase_plot()
        btn_reset.clicked.connect(_reset)
        btn_row.addWidget(btn_reset)

        btn_all = QPushButton("모두 표시")
        def _show_all():
            for k in keys:
                self.filter_visibility[k] = True
                vis_chks[k].setChecked(True)
            self._update_phase_plot()
        btn_all.clicked.connect(_show_all)
        btn_row.addWidget(btn_all)

        btn_row.addStretch()
        btn_close = QPushButton("닫기")
        btn_close.clicked.connect(dialog.accept)
        btn_row.addWidget(btn_close)
        layout.addLayout(btn_row)

        dialog.exec_()

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _toggle_log(self, checked: bool):
        self.log_box.setVisible(checked)
        self.btn_log_toggle.setText("Log ▲" if checked else "Log ▼")

    def log(self, msg: str):
        self.log_box.append(msg)
