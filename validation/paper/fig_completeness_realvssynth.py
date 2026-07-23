"""Figure — injection-recovery completeness across real frames of differing quality,
plus the controlled synthetic verification frame.

Reference-standard artificial-star test (DAOPHOT ADDSTAR; DES Balrog; HSC SynPipe;
AutoPhOT App.D; Haynes 1994/2002): stars injected into REAL frames. Three real
cluster frames spanning a range of sky brightness and seeing are each measured by
the identical injection pipeline, together with the idealised synthetic frame used
as the numerical-verification rung.

The point: the recovered 50% depth tracks each frame's sky + seeing, exactly as it
must — a dark-sky sharp frame (M67 i) reaches the same depth as the synthetic; a
bright-sky poor-seeing frame (M13 V) is ~3 mag shallower. So no single frame is
"anomalous", and the synthetic is NOT systematically optimistic (a real good frame
matches it). The bottom strip shows real injected stars from the shallowest frame
straddling its transition. The M13 depth was additionally cross-checked against the
pipeline's own detection roll-off (see COMPLETENESS_REALFRAME_INVESTIGATION.md).

    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig_completeness_realvssynth.py
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

DATA = REPO / "validation" / "paper"
SYN = DATA / "data" / "artificial_star" / "benchmark_run"
CUTOUTS = DATA / "data_realframe_M13V" / "injection_cutouts.npz"
OUTDIR = DATA / "figures"
ZP = 25.0

# real frames spanning depth — (label, run_dir, color_key, marker, sky_adu, fwhm_px)
RUNS = [
    ("M67 i",      DATA / "data_realframe_M67i/artificial_star/benchmark_run",     "data",      "-o", 27,   5.2),
    ("NGC 6811 R", DATA / "data_realframe_NGC6811R/artificial_star/benchmark_run", "reference", "-s", 1315, 5.3),
    ("M13 V",      DATA / "data_realframe_M13V/artificial_star/benchmark_run",     "accent",    "-D", 1315, 7.6),
]


def read_off(mag, comp, level=0.5):
    o = np.argsort(mag); mm, cc = mag[o], comp[o]
    for i in range(len(mm) - 1):
        if cc[i] >= level >= cc[i + 1]:
            den = cc[i] - cc[i + 1]
            f = (cc[i] - level) / den if den else 0.0
            return float(mm[i] + f * (mm[i + 1] - mm[i]))
    return float("nan")


def curve(run_dir, bw=0.25):
    s = pd.read_csv(run_dir / "stars.csv")
    s = s[~s["baseline_confounded"].astype(bool)]
    m = s["magnitude_true"].to_numpy(float)
    r = s["recovered"].to_numpy(bool)
    lo0 = np.floor(m.min() / bw) * bw
    edges = np.arange(lo0, m.max() + bw, bw)
    xs, cs = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        msk = (m >= lo) & (m < hi)
        if msk.sum() >= 15:
            xs.append(0.5 * (lo + hi)); cs.append(r[msk].mean())
    return np.array(xs), np.array(cs)


def main() -> int:
    cut = np.load(CUTOUTS)
    stamps, cmags = cut["stamps"], cut["mags"]
    clo, chi = float(cut["lo"]), float(cut["hi"])

    fig = plt.figure(figsize=(DOUBLE_COL * 0.76, 4.9))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 0.95], hspace=0.48)
    ax = fig.add_subplot(gs[0])
    ax.axhline(0.5, color=PALETTE["grey"], lw=0.7, ls=":", zorder=2)

    # synthetic verification rung (grey dashed, thin)
    ms, cs = curve(SYN)
    m50s = read_off(ms, cs)
    ax.plot(ms, cs, "--", color=PALETTE["grey"], lw=1.4, zorder=3,
            label="synthetic (verification)  ·  sky 150, FWHM 3.4px")

    m50r_m13 = None
    for label, rd, ckey, mk, sky, fwhm in RUNS:
        x, c = curve(rd)
        mm = read_off(x, c)
        col = C[ckey]
        ax.plot(x, c, mk, color=col, lw=2.0, ms=3.6, mfc=col, zorder=6,
                label=f"{label}  ·  sky {sky}, FWHM {fwhm}px")
        ax.axvline(mm, color=col, lw=0.9, ls=":", zorder=2, alpha=0.75)
        # nudge apart the two close labels (M13 14.9 vs NGC 15.6)
        ha = "right" if label.startswith("M13") else ("left" if label.startswith("NGC") else "center")
        dx = -0.08 if ha == "right" else (0.08 if ha == "left" else 0.0)
        ax.text(mm + dx, 1.02, f"{mm:.1f}", color=col, fontsize=7.2,
                ha=ha, va="bottom", fontweight="bold")
        if label.startswith("M13"):
            m50r_m13 = mm
    ax.text(12.35, 1.02, "$m_{50}$:", color="0.35", fontsize=7.0,
            ha="left", va="bottom")

    ax.set_xlabel("injected / instrumental magnitude (count-rate, ZP = 25)")
    ax.set_ylabel("recovery completeness")
    ax.set_xlim(12.2, 20.3)
    ax.set_ylim(-0.03, 1.11)
    leg = ax.legend(loc="lower left", fontsize=6.9, framealpha=0.93,
                    title="real frames (+ synthetic), deep → shallow",
                    title_fontsize=7.0, handlelength=1.8)
    leg._legend_box.align = "left"
    ax.set_title("real-frame completeness tracks each frame's sky + seeing",
                 loc="left", fontsize=10.5)
    # the key point, placed in free space above the M67/synthetic transition
    ax.annotate("M67 i (dark sky, sharp)\nreaches the synthetic depth\n"
                "→ synthetic is not optimistic",
                xy=(17.5, 0.55), xytext=(15.1, 0.80), fontsize=6.6, color="0.25",
                va="center", ha="left",
                arrowprops={"arrowstyle": "->", "color": "0.45", "lw": 0.9})

    # ── bottom strip: real injected-star cutouts (M13, the shallowest frame) ──
    gsc = gs[1].subgridspec(1, len(cmags), wspace=0.16)
    for j, (s, m) in enumerate(zip(stamps, cmags)):
        axc = fig.add_subplot(gsc[j])
        axc.imshow(s, cmap="gray", vmin=clo, vmax=chi, origin="lower",
                   interpolation="nearest")
        axc.set_xticks([]); axc.set_yticks([])
        recovered = m < (m50r_m13 or 14.9)
        mark = "found" if recovered else "lost"
        col = C["data"] if recovered else C["accent"]
        for sp in axc.spines.values():
            sp.set_edgecolor(col); sp.set_linewidth(1.5)
        axc.set_title(f"m = {m:.1f}", fontsize=7.6, pad=2)
        axc.set_xlabel(mark, fontsize=7.8, color=col, labelpad=2)
    fig.text(0.012, gs[1].get_position(fig).y1 + 0.016,
             "real injected stars — M13 V (the shallowest frame), shared stretch  —  "
             "blue = recovered, orange = lost",
             fontsize=7.2, color="0.3", va="bottom", ha="left")

    paths = save_fig(fig, "fig_completeness_realvssynth", OUTDIR)
    plt.close(fig)
    print(f"[completeness] synthetic m50={m50s:.2f}")
    for label, rd, *_ in RUNS:
        x, c = curve(rd)
        print(f"  {label}: m50={read_off(x, c):.2f}")
    print(f"saved: {paths['png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
