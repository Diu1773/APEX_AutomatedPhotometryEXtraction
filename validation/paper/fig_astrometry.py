"""Figure — Astrometric solution validated against Gaia (real data).

Reproduces the validation style of Ofek 2019 (PASP 131,054504) and Masci 2019
(ZTF, PASP 131,018003): for every science frame, match detected sources to Gaia
DR3 and report the residual RMS and the solve reliability. No fitting, no
synthetic data — this is the WCS the pipeline actually delivered on 66 real
frames of three clusters (Moravian C3-61000).

    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig_astrometry.py
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

REPRO = Path(r"E:\APEX_validation\reprocess")
OUTDIR = REPO / "validation" / "paper" / "figures"
TARGETS = {"M13": "M13", "M67": "M67", "NGC6811": "NGC 6811"}
COLORS = {"M13": C["data"], "M67": C["reference"], "NGC6811": C["accent"]}


def load():
    rows = []
    for key in TARGETS:
        p = REPRO / key / "result" / "step5_wcs" / "frame_wcs_qc.csv"
        if p.exists():
            d = pd.read_csv(p)
            d["target"] = key
            rows.append(d)
    d = pd.concat(rows, ignore_index=True)
    d["rms_as"] = d["rms_px"] * d["pixscale"]      # residual RMS in arcsec
    return d


def main() -> int:
    d = load()
    n = len(d)
    n_solved = int(d["wcs_ok"].sum())
    med = float(d["rms_as"].median())

    fig = plt.figure(figsize=(DOUBLE_COL, 3.3))
    gs = GridSpec(1, 2, width_ratios=[1.0, 1.15], wspace=0.30, figure=fig)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])

    # ── (a) per-frame residual RMS by target (strip) ──
    rng = np.random.default_rng(3)
    for i, (key, label) in enumerate(TARGETS.items()):
        sub = d[d["target"] == key]["rms_as"].to_numpy()
        x = i + rng.uniform(-0.16, 0.16, len(sub))
        axa.plot(x, sub, "o", color=COLORS[key], ms=4, alpha=0.7, mew=0, zorder=3)
        axa.plot(i, np.median(sub), "_", color="k", ms=22, mew=2.2, zorder=4)
    axa.axhline(med, color=PALETTE["grey"], lw=0.9, ls="--", zorder=1)
    axa.axhline(0.30, color=C["reference"], lw=0.9, ls=":", zorder=1, alpha=0.7)
    axa.text(2.42, 0.305, "PP (Mommert 2017)", fontsize=6.5, color=C["reference"],
             va="bottom", ha="right")
    axa.set_xticks(range(len(TARGETS)))
    axa.set_xticklabels(list(TARGETS.values()))
    axa.set_ylabel("Gaia residual RMS (arcsec)")
    axa.set_xlim(-0.5, len(TARGETS) - 0.5)
    axa.set_ylim(0, max(0.42, d["rms_as"].max() * 1.08))
    axa.set_title(f"(a) residual vs Gaia (median {med:.2f}″)", loc="left")
    axa.annotate(f"N = {n} frames, {n_solved}/{n} solved",
                 xy=(0.04, 0.06), xycoords="axes fraction", fontsize=7.4,
                 bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
                       "alpha": .85, "edgecolor": PALETTE["grey"]})

    # ── (b) residual vs number of matched Gaia stars (stability) ──
    for key, label in TARGETS.items():
        sub = d[d["target"] == key]
        axb.plot(sub["n_match"], sub["rms_as"], "o", color=COLORS[key], ms=4,
                 alpha=0.75, mew=0, label=label, zorder=3)
    axb.axhline(med, color=PALETTE["grey"], lw=0.9, ls="--", zorder=1)
    axb.set_xlabel("Gaia stars matched per frame")
    axb.set_ylabel("Gaia residual RMS (arcsec)")
    axb.set_ylim(0, max(0.42, d["rms_as"].max() * 1.08))
    axb.legend(loc="upper right", fontsize=7.0)
    axb.set_title("(b) stable across star count", loc="left")

    paths = save_fig(fig, "fig_astrometry", OUTDIR)
    plt.close(fig)
    print(f"[astrom] N={n} solved={n_solved} median_rms={med:.3f}arcsec "
          f"match_rate_med={d['match_rate'].median():.3f} "
          f"inlier_med={d['inlier_rate'].median():.3f} nmatch_med={int(d['n_match'].median())}")
    print(f"saved: {paths['png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
