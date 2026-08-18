"""Step 7 QC figures, drawn once for both the window and the batch.

Two panels the Step 7 window has always drawn and never saved: how far the
forced positions had to be recentred, and the growth curve the aperture
correction is read off. A headless run produced `centering_stats.csv` and
`apcorr_summary.csv` and no picture of either.

The growth curve is the awkward one. Unlike every other Step 7 product it never
reaches disk — `run_forced_photometry` hands it to an `apcorr_cb` callback and
the window plots it live. So the batch has to draw during the run, through that
same callback, rather than rebuilding afterwards.

Bodies transcribed from `step7_forced_aperture_phot.py` so both routes draw the
same thing; the window keeps its canvases, its table selection and its
formatting helpers.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from apex.utils.step_paths import step7_forced_phot_dir


def _numeric(df: pd.DataFrame, name: str) -> pd.Series:
    if name not in df.columns:
        return pd.Series(np.nan, index=df.index, dtype=float)
    return pd.to_numeric(df[name], errors="coerce")


def _boolean(df: pd.DataFrame, name: str) -> pd.Series:
    if name not in df.columns:
        return pd.Series(False, index=df.index, dtype=bool)
    return df[name].astype(str).str.lower().isin(("1", "true", "yes"))


def _finite(value, default=np.nan) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if np.isfinite(out) else default


def init_center_axes(ax_hist, ax_scatter) -> None:
    ax_hist.set_title("Center error", fontsize=9)
    ax_hist.set_xlabel("center_error_px", fontsize=9)
    ax_hist.set_ylabel("N", fontsize=9)
    ax_hist.grid(True, alpha=0.3)
    ax_scatter.set_title("Error vs mag_err", fontsize=9)
    ax_scatter.set_xlabel("center_error_px", fontsize=9)
    ax_scatter.set_ylabel("mag_err", fontsize=9)
    ax_scatter.grid(True, alpha=0.3)


def draw_center_shift(ax_hist, ax_scatter, result_dir, fname: str,
                      row: Optional[dict] = None,
                      outlier_default: float = 1.0) -> bool:
    """Recentring distance: its distribution, and what it cost in magnitude error."""
    ax_hist.cla()
    ax_scatter.cla()
    init_center_axes(ax_hist, ax_scatter)
    if not fname:
        return False

    path = Path(step7_forced_phot_dir(result_dir)) / f"photometry_{fname}.tsv"
    if not path.exists():
        ax_hist.set_title(f"{fname} not found", fontsize=9)
        return False
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as exc:                        # noqa: BLE001
        ax_hist.set_title(f"Read failed: {exc}", fontsize=9)
        return False

    center_error = _numeric(df, "center_error_px")
    mag_err = _numeric(df, "mag_err")
    outlier = _boolean(df, "centroid_outlier")

    finite_center = center_error[np.isfinite(center_error)]
    if len(finite_center) > 0:
        ax_hist.hist(finite_center.to_numpy(float), bins=30,
                     color="#1565C0", alpha=0.78)
    outlier_px = _finite((row or {}).get("centroid_outlier_px"), outlier_default)
    if np.isfinite(outlier_px):
        ax_hist.axvline(outlier_px, color="#E53935", ls="--", lw=1.2,
                        label=f"outlier={outlier_px:.2f}px")
        ax_hist.legend(fontsize=7, frameon=False)

    both = center_error.notna() & mag_err.notna()
    if both.any():
        idx = np.where(both.to_numpy())[0]
        if len(idx) > 5000:                         # the window's own sample cap
            idx = idx[np.linspace(0, len(idx) - 1, 5000).astype(int)]
        colors = np.where(outlier.iloc[idx].to_numpy(bool), "#E53935", "#2E7D32")
        ax_scatter.scatter(center_error.iloc[idx].to_numpy(float),
                           mag_err.iloc[idx].to_numpy(float),
                           s=8, c=colors, alpha=0.55, linewidths=0)

    title = str(fname)
    if row:
        p90 = _finite(row.get("center_error_p90_px"))
        rate = _finite(row.get("centroid_outlier_rate"))
        if np.isfinite(p90):
            title += f" | p90={p90:.3f}px"
        if np.isfinite(rate):
            title += f" | outliers={100.0 * rate:.2f}%"
    ax_hist.set_title(title, fontsize=9)
    return True


def init_gc_axes(ax_mag, ax_err) -> None:
    ax_mag.set_ylabel("Inst Magnitude", fontsize=9)
    ax_mag.set_title("Growth Curve", fontsize=9)
    ax_mag.grid(True, alpha=0.3)
    ax_err.set_xlabel("Aperture radius (px)", fontsize=9)
    ax_err.set_ylabel("Median mag_err", fontsize=9)
    ax_err.set_title("Error vs Aperture (U-shape)", fontsize=9)
    ax_err.grid(True, alpha=0.3)


def draw_growth_curve(ax_mag, ax_err, gc: dict) -> bool:
    """The curve the aperture correction is read off, with the radii marked."""
    ax_mag.cla()
    ax_err.cla()
    init_gc_axes(ax_mag, ax_err)
    if not gc:
        return False

    fwhm = _finite(gc.get("fwhm_px"), 0.0)
    r_ap = _finite(gc.get("r_ap_px"), 0.0)
    r_ref = _finite(gc.get("r_ref_px"), 0.0)
    r_opt = _finite(gc.get("r_opt_px"), 0.0)
    fname = str(gc.get("fname", "") or gc.get("file", ""))
    apcorr_val = _finite(gc.get("apcorr"))

    radii = gc.get("radii_px", [])
    encs = gc.get("enclosed_frac", [])
    errs = gc.get("mag_err", [])

    if len(radii) and len(encs):
        arr = np.asarray(encs, dtype=float)
        arr = np.where((arr > 0) & np.isfinite(arr), arr, np.nan)
        ax_mag.plot(radii, -2.5 * np.log10(arr), "-o", color="#1565C0",
                    lw=1.5, markersize=5, markeredgecolor="white",
                    markeredgewidth=0.5)
    if len(radii) and len(errs):
        ax_err.plot(radii, errs, "-s", color="#E53935", lw=1.8, markersize=5,
                    markeredgecolor="white", markeredgewidth=0.5)

    for ax in (ax_mag, ax_err):
        if r_opt > 0:
            ax.axvline(r_opt, color="#7B1FA2", lw=1.6, ls="-", alpha=0.85,
                       label=f"r_opt={r_opt:.1f}px")
        if r_ap > 0:
            ax.axvline(r_ap, color="#E53935", lw=1.2, ls="--", alpha=0.8,
                       label=f"r_ap={r_ap:.1f}px")
        if r_ref > 0:
            ax.axvline(r_ref, color="#43A047", lw=1.2, ls="--", alpha=0.8,
                       label=f"r_ref={r_ref:.1f}px")
        if fwhm > 0:
            ax.axvline(fwhm, color="#6D4C41", lw=1.0, ls=":", alpha=0.8,
                       label=f"FWHM={fwhm:.2f}px")

    ax_mag.invert_yaxis()
    title = fname
    if np.isfinite(apcorr_val):
        title += f"  |  apcorr={apcorr_val:.4f}"
    ax_mag.set_title(title, fontsize=9)
    ax_mag.legend(fontsize=7, frameon=False, loc="best")
    ax_err.legend(fontsize=7, frameon=False, loc="best")
    return True


def export_forcedphot_qc(result_dir, params=None,
                         growth_curves=None) -> list[Path]:
    """Write both figures. `growth_curves` maps frame name to the run's payload.

    The centre-shift panel is rebuilt from disk. The growth curve cannot be —
    it is never written there — so the caller has to have kept what the run
    handed to `apcorr_cb`. One figure per filter is enough to show the shape;
    fifteen near-identical curves would only bury it.
    """
    from matplotlib.figure import Figure

    out_dir = Path(step7_forced_phot_dir(result_dir))
    if not out_dir.exists():
        return []
    saved: list[Path] = []

    stats_path = out_dir / "centering_stats.csv"
    if stats_path.exists():
        try:
            stats = pd.read_csv(stats_path)
        except Exception:                           # noqa: BLE001
            stats = pd.DataFrame()
        if not stats.empty:
            row = stats.iloc[0].to_dict()
            fig = Figure(figsize=(10.0, 3.8), dpi=120)
            ax_hist = fig.add_subplot(121)
            ax_scatter = fig.add_subplot(122)
            outlier_default = _finite(
                getattr(getattr(params, "P", None), "centroid_outlier_px", 1.0), 1.0)
            if draw_center_shift(ax_hist, ax_scatter, result_dir,
                                 str(row.get("file", "")), row, outlier_default):
                fig.tight_layout()
                path = out_dir / "step7_center_shift.png"
                fig.savefig(path, dpi=160, bbox_inches="tight")
                saved.append(path)

    for fname, gc in sorted((growth_curves or {}).items()):
        payload = dict(gc)
        payload.setdefault("fname", fname)
        fig = Figure(figsize=(6.5, 6.0), dpi=120)
        ax_mag = fig.add_subplot(211)
        ax_err = fig.add_subplot(212)
        if draw_growth_curve(ax_mag, ax_err, payload):
            fig.tight_layout()
            path = out_dir / f"step7_growth_curve_{Path(fname).stem}.png"
            fig.savefig(path, dpi=160, bbox_inches="tight")
            saved.append(path)
    return saved


def one_per_filter(growth_curves: dict) -> dict:
    """Keep the first curve of each filter. The rest look the same."""
    picked: dict = {}
    seen: set = set()
    for fname, gc in sorted((growth_curves or {}).items()):
        key = str((gc or {}).get("filter", "")) or "?"
        if key in seen:
            continue
        seen.add(key)
        picked[fname] = gc
    return picked
