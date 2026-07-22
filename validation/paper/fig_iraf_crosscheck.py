"""Figure — independent-engine cross-check on real data (APEX vs IRAF/DAOPHOT).

Reproduces the new-code-vs-established-code test of Schechter+1993 (DoPHOT vs
DAOPHOT, PASP 105,1342) and AutoPhOT Fig.14: measure the same stars, at the same
fixed sky coordinates, with two independent flux integrators and check that the
residual is flat with magnitude (no systematic). Data: fully-APEX-reduced NGC 6811
V frame, 499 stars, APEX forced photometry vs IRAF `phot` (DAOPHOT via PyRAF).

    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig_iraf_crosscheck.py
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
from matplotlib.gridspec import GridSpec

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

CSV = REPO / "benchmark" / "runs" / "ngc6811_iraf_allapex_v1" / "phot_fixed_coords" / "fixed_comparison.csv"
OUTDIR = REPO / "validation" / "paper" / "figures"


def main() -> int:
    d = pd.read_csv(CSV)
    a = pd.to_numeric(d["apex_mag_iraf_units"], errors="coerce").to_numpy(float)
    i = pd.to_numeric(d["iraf_mag"], errors="coerce").to_numpy(float)
    dc = pd.to_numeric(d["delta_mag_units_centered"], errors="coerce").to_numpy(float)
    k = np.isfinite(a) & np.isfinite(i) & np.isfinite(dc)
    a, i, dc = a[k], i[k], dc[k]
    mad = 1.4826 * np.median(np.abs(dc - np.median(dc)))
    rms = float(np.sqrt(np.mean(dc ** 2)))
    r = float(np.corrcoef(a, i)[0, 1])

    fig = plt.figure(figsize=(DOUBLE_COL * 0.66, 3.5))
    gs = GridSpec(1, 2, width_ratios=[3.1, 1.0], wspace=0.05, figure=fig)
    ax = fig.add_subplot(gs[0, 0])
    axh = fig.add_subplot(gs[0, 1], sharey=ax)

    ax.axhline(0.0, color=PALETTE["grey"], lw=0.9, ls="--", zorder=1)
    ax.plot(i, dc, "o", color=C["data"], ms=3.2, alpha=0.5, mew=0, zorder=3)
    # binned running median — Schechter's systematic check
    bw = 0.6
    edges = np.arange(np.floor(i.min()), np.ceil(i.max()) + bw, bw)
    bx, bmed = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (i >= lo) & (i < hi)
        if m.sum() >= 5:
            bx.append(0.5 * (lo + hi)); bmed.append(np.median(dc[m]))
    ax.plot(bx, bmed, "-", color=C["accent"], lw=1.8, zorder=5,
            label="binned median")
    ax.set_xlabel("IRAF magnitude")
    ax.set_ylabel(r"$\Delta$mag  (APEX $-$ IRAF)")
    ax.set_ylim(-0.16, 0.16)
    ax.legend(loc="upper left", fontsize=7.2)
    ax.set_title("real NGC 6811 V — residual flat to the faint end", loc="left")
    ax.text(0.03, 0.06,
            f"N = {len(a)}\nMAD = {1000*mad:.1f} mmag\nRMS = {1000*rms:.1f} mmag\n"
            rf"$r$ = {r:.5f}",
            transform=ax.transAxes, va="bottom", fontsize=7.4,
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
                  "alpha": .85, "edgecolor": PALETTE["grey"]})

    # marginal histogram
    axh.hist(dc, bins=30, orientation="horizontal", color=C["data"], alpha=0.7)
    axh.axhline(0.0, color=PALETTE["grey"], lw=0.9, ls="--")
    axh.set_xticks([])
    axh.tick_params(labelleft=False)

    paths = save_fig(fig, "fig_iraf_crosscheck", OUTDIR)
    plt.close(fig)
    print(f"[iraf] N={len(a)} MAD={1000*mad:.1f}mmag RMS={1000*rms:.1f}mmag r={r:.5f}")
    print(f"saved: {paths['png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
