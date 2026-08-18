"""Step 6 reference-build QC, drawn once for both the window and the batch.

The window plotted `ref_frame_stats.csv` — a file the batch already writes —
into a canvas and never saved it, so a headless run produced the numbers and no
picture of them. Fourth instance of the same shape this week, after the Step 8
comparison figure, the Step 12 isochrone plots and the Step 4 QC overview.

The window keeps the date selector; the drawing moved here and takes the chosen
date as an argument.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from apex.utils.step_paths import step6_refbuild_dir


def export_refbuild_qc(params, *, date_key: str = "All") -> list[Path]:
    """Write `step6_refbuild_overview.png`. Returns what was written."""
    out_dir = Path(step6_refbuild_dir(params.P.result_dir))
    if not (out_dir / "ref_frame_stats.csv").exists():
        return []
    out_dir.mkdir(parents=True, exist_ok=True)
    fig = Figure(figsize=(11.0, 4.6), dpi=120)
    if not draw_refbuild_overview(fig, params, date_key=date_key):
        return []
    path = out_dir / "step6_refbuild_overview.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    return [path]


def draw_refbuild_overview(fig, params, *, date_key='All', summary=None) -> bool:
    fig.clear()
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)

    stats_path = step6_refbuild_dir(params.P.result_dir) / "ref_frame_stats.csv"
    if not stats_path.exists():
        ax1.text(0.5, 0.5, "No ref stats available", ha="center", va="center")
        ax2.axis("off")
        return False
    try:
        df = pd.read_csv(stats_path)
    except Exception:
        ax1.text(0.5, 0.5, "Failed to read ref stats", ha="center", va="center")
        ax2.axis("off")
        return False
    if df.empty:
        ax1.text(0.5, 0.5, "No ref stats available", ha="center", va="center")
        ax2.axis("off")
        return False

    # 창은 여기서 날짜 콤보를 다시 채웠다. 배치에는 콤보가 없으므로
    # 고른 날짜를 인자로 받는다 — 'All' 이면 전부.
    selected_date = None
    if date_key and str(date_key).strip().lower() != "all":
        selected_date = str(date_key).strip()

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
        meta_path = step6_refbuild_dir(params.P.result_dir) / "ref_build_meta.json"
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
    return True
