"""Figure — APEX workflow and architecture (Section 2).

What the reader should get in one look:
  * the pipeline is a fixed sequence of steps, each writing a checkable product;
  * Steps 0-7 are shared and the two science modes branch only after photometry;
  * every step lives in a Qt-free core that the GUI merely drives, which is why the
    headless validation of Section 3 exercises exactly the code the interface runs;
  * the small label under each step is the subsection that validates it — the
    map between what the tool does and where this paper checks it.

Monochrome, per the paper figure convention.

    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig_architecture.py
"""
from __future__ import annotations
import sys
from pathlib import Path
REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from apex_paper_style import apply_paper_style, save_fig, DOUBLE_COL

apply_paper_style()
OUTDIR = REPO / "validation" / "paper" / "figures"

INK, MUT, LINE = "#1a1a1a", "#6b6b6b", "#9a9a9a"
F_SHARED, F_MODE = "#f0f0f0", "#e3e3e3"
ACCENT = "#3f3f3f"

SHARED = [("0  Calibration", "§3.2"), ("1  File selection", "—"),
          ("2  Crop", "—"), ("3  Sky QC", "§3.3"),
          ("4  Detection", "§3.4–3.5"), ("5  WCS solve", "§3.6"),
          ("6  Master cat.", "§5.2"), ("7  Photometry", "§3.7–3.9")]
CMD = [("8  PSF phot.", "§3.10–3.11"), ("9  Master IDs", "—"),
       ("10  Zeropoint", "§3.12"), ("11  CMD", "§3.12"),
       ("12  Isochrone", "§5.2")]
LC = [("8  Target sel.", "—"), ("9  Light curve", "§4"),
      ("10  Detrend", "§3.13"), ("11  Period", "§3.13, §4")]


def box(ax, x, y, w, h, label, sec, fill, fs=7.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.004,rounding_size=0.010",
                                lw=0.8, edgecolor=LINE, facecolor=fill, zorder=2))
    ax.text(x + w / 2, y + h * 0.55, label, ha="center", va="center",
            fontsize=fs, color=INK, zorder=3)
    if sec:
        ax.text(x + w / 2, y - 0.028, sec, ha="center", va="center", fontsize=6.0,
                color=ACCENT if sec != "—" else MUT,
                fontweight="semibold" if sec != "—" else "normal", zorder=3)


def arrow(ax, x0, y0, x1, y1, lw=0.8):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", color=LINE,
                                 lw=lw, mutation_scale=7, shrinkA=1, shrinkB=1, zorder=1))


def main() -> int:
    fig, ax = plt.subplots(1, 1, figsize=(DOUBLE_COL, 3.9))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    w, h = 0.205, 0.105
    xs = [0.045 + i * (w + 0.028) for i in range(4)]

    # header ------------------------------------------------------------
    ax.text(0.045, 0.955, "graphical layer (PyQt5)   drives  →   ", fontsize=7.2,
            color=MUT, style="italic", va="center")
    ax.text(0.36, 0.955, "Qt-free core (apex.analysis / apex.pipeline)",
            fontsize=7.4, color=INK, fontweight="bold", va="center")
    ax.text(0.845, 0.912, "same core used in validation", fontsize=7.0,
            color=ACCENT, va="center")

    # shared chain, two rows of four ------------------------------------
    ax.text(0.045, 0.885, "shared measurement chain", fontsize=7.6,
            color=INK, fontweight="bold", va="center")
    y1, y2 = 0.735, 0.545
    for i, (lab, sec) in enumerate(SHARED[:4]):
        box(ax, xs[i], y1, w, h, lab, sec, F_SHARED)
        if i:
            arrow(ax, xs[i] - 0.028, y1 + h / 2, xs[i], y1 + h / 2)
    for i, (lab, sec) in enumerate(SHARED[4:]):
        box(ax, xs[i], y2, w, h, lab, sec, F_SHARED)
        if i:
            arrow(ax, xs[i] - 0.028, y2 + h / 2, xs[i], y2 + h / 2)
    # row 1 -> row 2 wrap
    ax.plot([xs[3] + w + 0.012, xs[3] + w + 0.012], [y1 + h / 2, y2 + h / 2],
            color=LINE, lw=0.8)
    ax.plot([xs[3] + w, xs[3] + w + 0.012], [y1 + h / 2, y1 + h / 2], color=LINE, lw=0.8)
    arrow(ax, xs[3] + w + 0.012, y2 + h / 2, xs[0] - 0.014, y2 + h / 2, lw=0.0)
    ax.plot([xs[0] - 0.014, xs[3] + w + 0.012], [y2 + h / 2 + 0.075, y2 + h / 2 + 0.075],
            color=LINE, lw=0.8, ls=":")
    ax.plot([xs[0] - 0.014, xs[0] - 0.014], [y2 + h / 2 + 0.075, y2 + h / 2],
            color=LINE, lw=0.8, ls=":")
    arrow(ax, xs[0] - 0.014, y2 + h / 2, xs[0], y2 + h / 2)

    # branch -------------------------------------------------------------
    ybr = y2 - 0.075
    ax.plot([0.5, 0.5], [y2 - 0.045, ybr], color=LINE, lw=0.8)
    ax.plot([0.018, 0.5], [ybr, ybr], color=LINE, lw=0.8)
    ax.text(0.5, y2 - 0.062, "photometry complete — modes branch", fontsize=6.4,
            color=MUT, ha="center", style="italic")

    # mode rows ----------------------------------------------------------
    wm, hm = 0.163, 0.095
    yc, yl = 0.315, 0.098
    ax.text(0.045, yc + hm + 0.032, "CMD mode", fontsize=7.4, color=INK, fontweight="bold")
    for i, (lab, sec) in enumerate(CMD):
        x = 0.045 + i * (wm + 0.022)
        box(ax, x, yc, wm, hm, lab, sec, F_MODE, fs=6.6)
        if i:
            arrow(ax, x - 0.022, yc + hm / 2, x, yc + hm / 2)
    arrow(ax, 0.127, ybr, 0.127, yc + hm + 0.020)

    ax.text(0.045, yl + hm + 0.032, "LC mode", fontsize=7.4, color=INK, fontweight="bold")
    for i, (lab, sec) in enumerate(LC):
        x = 0.045 + i * (wm + 0.022)
        box(ax, x, yl, wm, hm, lab, sec, F_MODE, fs=6.6)
        if i:
            arrow(ax, x - 0.022, yl + hm / 2, x, yl + hm / 2)
    ax.plot([0.018, 0.018], [ybr, yl + hm / 2], color=LINE, lw=0.8)
    ax.plot([0.018, 0.127], [ybr, ybr], color=LINE, lw=0.8)
    arrow(ax, 0.018, yl + hm / 2, 0.045, yl + hm / 2)

    # products: stated in the caption, not crowded into the panel

    ax.text(0.045, 0.022,
            "label below each step = validation or application section",
            fontsize=6.2, color=ACCENT, style="italic")

    paths = save_fig(fig, "fig_architecture", OUTDIR)
    plt.close(fig)
    print(f"saved: {paths['png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
