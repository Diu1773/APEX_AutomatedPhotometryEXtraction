"""Step 5 astrometry QC — the one step that had no figure at all.

Steps 4, 6, 7, 8, 10 and 12 all drew something the window kept to itself, so
fixing them meant moving code. Step 5 is different: the window never plotted the
solve either. The numbers were there — `frame_wcs_qc.csv` carries 36 columns per
frame — and nothing looked at them together.

Which four panels was decided by measurement rather than taste. Across 137 real
frames from five clusters, the coefficient of variation of each candidate:

    resid_vs_radius_slope   1.75     <- moves most, and changes a decision
    match_rate_cat          1.12
    center_offset_arcsec    0.73
    edge_resid_ratio        0.68
    resid_mad_px            0.47
    rms_px                  0.35
    n_match                 0.30
    scale_delta_pct         0.05
    inlier_rate             0.02
    match_rate              0.01     <- a flat line; a wasted panel

So `match_rate`, `inlier_rate` and `scale_delta_pct` are reported as numbers in
the title and not given axes. A panel that is flat on every real dataset says
nothing and takes the space of one that would.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from apex.utils.step_paths import step5_wcs_dir

PASS_COLOR, FAIL_COLOR = "#212121", "#D32F2F"


def _num(df: pd.DataFrame, name: str) -> np.ndarray:
    if name not in df.columns:
        return np.full(len(df), np.nan)
    return pd.to_numeric(df[name], errors="coerce").to_numpy(float)


def load_wcs_qc(result_dir) -> pd.DataFrame:
    """The per-frame astrometry QC table Step 5 writes."""
    path = Path(step5_wcs_dir(result_dir)) / "frame_wcs_qc.csv"
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except Exception:                               # noqa: BLE001
        return pd.DataFrame()


def draw_wcs_qc_overview(fig, df: pd.DataFrame) -> bool:
    """Four panels: residual size, distortion left over, yield, pointing."""
    fig.clear()
    if df is None or df.empty:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No WCS QC table", transform=ax.transAxes,
                ha="center", va="center", color="gray")
        ax.set_xticks([]); ax.set_yticks([])
        return False

    passed = df.get("wcs_qc_pass")
    ok = (passed.astype(str).str.lower().isin(("true", "1", "yes")).to_numpy()
          if passed is not None else np.ones(len(df), dtype=bool))
    colors = np.where(ok, PASS_COLOR, FAIL_COLOR)
    x = np.arange(len(df), dtype=float)

    ax_res = fig.add_subplot(2, 2, 1)
    ax_dist = fig.add_subplot(2, 2, 2)
    ax_yield = fig.add_subplot(2, 2, 3)
    ax_point = fig.add_subplot(2, 2, 4)

    # ① Residual size. The headline number, and its tail — a good median with a
    #    bad p99 means a few frames carry the error, not the solution.
    rms = _num(df, "rms_px")
    p99 = _num(df, "resid_p99_px")
    ax_res.scatter(x, rms, s=26, c=colors, marker="o", label="rms")
    ax_res.scatter(x, p99, s=20, c=colors, marker="^", alpha=0.55, label="p99")
    ax_res.set_title("Residual size", fontsize=9)
    ax_res.set_xlabel("Frame index", fontsize=9)
    ax_res.set_ylabel("px", fontsize=9)
    ax_res.legend(fontsize=7, frameon=False)
    ax_res.grid(True, alpha=0.3)

    # ② Distortion left on the table. A positive slope, or an edge ratio above
    #    one, means the residual grows outward — the signature of a SIP order
    #    too low for this optic. Negative means the *centre* is worse, which is
    #    not distortion and usually means crowding.
    slope = _num(df, "resid_vs_radius_slope")
    edge = _num(df, "edge_resid_ratio")
    ax_dist.scatter(slope, edge, s=30, c=colors)
    ax_dist.axvline(0.0, color="#1565C0", lw=1.0, ls="--", alpha=0.7,
                    label="no radial trend")
    ax_dist.axhline(1.0, color="#43A047", lw=1.0, ls=":", alpha=0.7,
                    label="edge = centre")
    ax_dist.set_title("Distortion left over", fontsize=9)
    ax_dist.set_xlabel("resid vs radius slope", fontsize=9)
    ax_dist.set_ylabel("edge / centre residual", fontsize=9)
    ax_dist.legend(fontsize=7, frameon=False)
    ax_dist.grid(True, alpha=0.3)

    # ③ Yield. Detected against matched, with the 1:1 line — a frame far below
    #    it detected sources the solver could not identify.
    n_det = _num(df, "n_detect")
    n_match = _num(df, "n_match")
    ax_yield.scatter(n_det, n_match, s=30, c=colors)
    finite = np.isfinite(n_det) & np.isfinite(n_match)
    if finite.any():
        hi = float(np.nanmax([n_det[finite].max(), n_match[finite].max()]))
        ax_yield.plot([0, hi], [0, hi], "k--", lw=0.8, alpha=0.5, zorder=0)
    ax_yield.set_title("Matched vs detected", fontsize=9)
    ax_yield.set_xlabel("n_detect", fontsize=9)
    ax_yield.set_ylabel("n_match", fontsize=9)
    ax_yield.grid(True, alpha=0.3)

    # ④ Pointing. How far the solution moved from the hint — a mount property,
    #    but a large scatter here is why a solve occasionally needs a blind run.
    offset = _num(df, "center_offset_arcsec")
    ax_point.scatter(x, offset, s=26, c=colors)
    if np.isfinite(offset).any():
        ax_point.axhline(float(np.nanmedian(offset)), color="#1565C0", lw=1.2,
                         ls="-", alpha=0.8,
                         label=f"median {np.nanmedian(offset):.0f}\"")
        ax_point.legend(fontsize=7, frameon=False)
    ax_point.set_title("Pointing offset", fontsize=9)
    ax_point.set_xlabel("Frame index", fontsize=9)
    ax_point.set_ylabel("arcsec", fontsize=9)
    ax_point.grid(True, alpha=0.3)

    fig.suptitle(_title(df, ok), fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return True


def _title(df: pd.DataFrame, ok: np.ndarray) -> str:
    """The three flat metrics belong here, as numbers, not as panels."""
    parts = [f"Step 5 WCS QC | PASS {int(ok.sum())}/{len(df)}"]
    for col, label, fmt in (("match_rate", "match", "{:.4f}"),
                            ("inlier_rate", "inlier", "{:.3f}"),
                            ("scale_delta_pct", "scale Δ", "{:+.2f}%")):
        values = _num(df, col)
        if np.isfinite(values).any():
            parts.append(f"{label} {fmt.format(float(np.nanmedian(values)))}")
    solver = df.get("solver")
    if solver is not None and len(solver):
        parts.append(str(solver.dropna().iloc[0]) if solver.notna().any() else "")
    return "  |  ".join(p for p in parts if p)


def export_wcs_qc(result_dir, params=None) -> list[Path]:
    """Write `step5_wcs_qc_overview.png`. Returns what was written."""
    from matplotlib.figure import Figure

    df = load_wcs_qc(result_dir)
    if df.empty:
        return []
    out_dir = Path(step5_wcs_dir(result_dir))
    out_dir.mkdir(parents=True, exist_ok=True)
    fig = Figure(figsize=(11.0, 7.4), dpi=120)
    if not draw_wcs_qc_overview(fig, df):
        return []
    path = out_dir / "step5_wcs_qc_overview.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    return [path]
