"""Figure (all-APEX) — NGC 6811 CMD from the fully-APEX reduction.

Same as the paper's CMD figure but from the raw->science all-APEX pipeline
(Step 0 detector calibration -> Steps 1-7 -> Step 10), not AIPPI-preprocessed
frames. Panel (a): the Johnson CMD (V vs B-V), APEX ground-based aperture
photometry vs the Gaia-transformed space-based reference for the SAME stars,
with sigma-clipped main-sequence ridgelines. Panel (b): the ridgeline residual
(APEX - Gaia) vs V, quantifying the agreement.

    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig_apex_cmd.py
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

ZP = Path(r"E:/APEX_validation/reprocess/NGC6811/result/cmd_zeropoint")
CMD_CSV = ZP / "median_by_ID_filter_wide_cmd.csv"
CAL_CSV = ZP / "gaia_sdss_calibrator_by_ID.csv"
OUTDIR = REPO / "validation" / "paper" / "figures"


def ridgeline(color, mag, edges):
    """Sigma-clipped median colour per mag bin (rejects field/binary outliers)."""
    xs, ys = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = np.isfinite(color) & np.isfinite(mag) & (mag >= lo) & (mag < hi)
        if k.sum() < 10:
            continue
        c = color[k]
        for _ in range(2):
            med, sig = np.median(c), 1.4826 * np.median(np.abs(c - np.median(c)))
            keep = np.abs(c - med) < 2.5 * max(sig, 1e-3)
            if keep.sum() < 8 or keep.all():
                break
            c = c[keep]
        xs.append(0.5 * (lo + hi)); ys.append(float(np.median(c)))
    return np.array(ys), np.array(xs)


def main() -> int:
    cmd = pd.read_csv(CMD_CSV)
    cal = pd.read_csv(CAL_CSV)[["ID", "ref_B", "ref_V"]]
    df = cmd.merge(cal, on="ID", how="left")

    apex_V = pd.to_numeric(df["mag_std_V"], errors="coerce").to_numpy(float)
    apex_BV = (pd.to_numeric(df["mag_std_B"], errors="coerce")
               - pd.to_numeric(df["mag_std_V"], errors="coerce")).to_numpy(float)
    gaia_V = pd.to_numeric(df["ref_V"], errors="coerce").to_numpy(float)
    gaia_BV = (pd.to_numeric(df["ref_B"], errors="coerce")
               - pd.to_numeric(df["ref_V"], errors="coerce")).to_numpy(float)

    edges = np.arange(12.8, 17.8, 0.4)
    apex_rl, apex_rlm = ridgeline(apex_BV, apex_V, edges)
    gaia_rl, gaia_rlm = ridgeline(gaia_BV, gaia_V, edges)
    common = np.intersect1d(apex_rlm, gaia_rlm)
    ai, gi = np.isin(apex_rlm, common), np.isin(gaia_rlm, common)
    resid = apex_rl[ai] - gaia_rl[gi]
    ridge_rms = float(np.sqrt(np.nanmean(resid ** 2)))
    n_ok = int(np.isfinite(apex_V + apex_BV).sum())

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 4.0),
                                     gridspec_kw={"width_ratios": [1.4, 1.0]})

    ax_a.scatter(apex_BV, apex_V, s=6, alpha=0.35, color=C["data"],
                 edgecolors="none", label="APEX (all-APEX, ground)", zorder=2)
    ax_a.scatter(gaia_BV, gaia_V, s=6, alpha=0.35, color=C["accent"],
                 edgecolors="none", label="Gaia-transformed (space)", zorder=2)
    ax_a.plot(apex_rl, apex_rlm, "-", color=C["data"], lw=2.0, zorder=4)
    ax_a.plot(gaia_rl, gaia_rlm, "--", color=C["model"], lw=2.0, zorder=5,
              label="MS ridgelines")
    ax_a.set_xlim(-0.1, 1.8); ax_a.set_ylim(18.4, 11.0)
    ax_a.set_xlabel(r"$B-V$"); ax_a.set_ylabel(r"$V$")
    ax_a.text(0.05, 0.05, f"ridgeline RMS = {ridge_rms*1000:.0f} mmag\nN = {n_ok}",
              transform=ax_a.transAxes, va="bottom", ha="left", fontsize=8.0,
              bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
                    "alpha": 0.85, "edgecolor": PALETTE["grey"]})
    ax_a.legend(loc="upper right", fontsize=7.2)
    ax_a.set_title("(a) Johnson CMD — all-APEX vs Gaia", loc="left")

    ax_b.axhline(0.0, color=PALETTE["grey"], lw=0.8, ls="--")
    ax_b.plot(apex_rlm[ai], resid * 1000, "o-", color=C["model"], lw=1.4, ms=4)
    ax_b.set_xlabel(r"$V$"); ax_b.set_ylabel(r"ridgeline $\Delta(B-V)$ (mmag)")
    ax_b.set_ylim(-60, 60)
    ax_b.set_title("(b) MS ridgeline residual", loc="left")

    fig.suptitle("NGC 6811 CMD — fully-APEX (raw$\\to$science) reduction",
                 fontsize=8.0, y=1.01, color="#333333")
    fig.tight_layout()
    paths = save_fig(fig, "fig_apex_cmd", OUTDIR)
    plt.close(fig)
    print(f"N={n_ok}  ridgeline RMS={ridge_rms*1000:.1f} mmag")
    print(f"saved: {paths['png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
