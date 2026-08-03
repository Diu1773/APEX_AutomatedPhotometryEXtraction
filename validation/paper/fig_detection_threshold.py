"""Figure — detection threshold: spurious contamination measured from the frame
itself, and the safe threshold floor moving frame to frame.

Method (Serra, Jurek & Floer 2012; applied to optical images as in Molino+2014):
the noise of a background-subtracted frame is symmetric about zero while real
sources are positive, so re-running the *same* detector on the sign-flipped image
counts the noise-origin spurious detections directly. No WCS, no external
catalogue, no injection — it runs inside step 4, before plate solving.

    purity_est = (N_pos - N_neg) / N_pos

Panel (a): purity versus threshold for five real frames (two globular, two open
clusters; B, R and g'). The floor at which purity crosses 0.95 differs frame to
frame — including between two filters of the *same* cluster on the same night —
while the collapse point is common at sigma 1.5 -> 1.2.

Panel (b): validation against Gaia DR3. Detections matched to a Gaia source
within 2.0 arcsec are counted real; the rest are the realized spurious count.
The estimate needs no catalogue and reproduces it.

Inputs: data_detection_threshold/*.json, written by
    validation/gui_tools/sigma_qc_scan.py <fits> [gaia.ecsv]

    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig_detection_threshold.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import matplotlib.pyplot as plt

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

SRC = REPO / "validation" / "paper" / "data_detection_threshold"
OUTDIR = REPO / "validation" / "paper" / "figures"
PUR_MIN = 0.95
NOMINAL = 3.2

# key, label, colour, linestyle, marker.  M13 appears twice (R and B) on purpose:
# it is the same cluster and the same night, so the two curves isolate the filter.
SERIES = [
    ("M13_R",   "M13 $R$",       C["accent"],    "-",  "o"),
    ("M13",     "M13 $B$",       C["accent"],    "--", "s"),
    ("M3",      "M3 $B$",        C["bad"],       "-",  "^"),
    ("M67",     "M67 $g'$",      C["data"],      "-",  "D"),
    ("NGC6811", "NGC 6811 $B$",  C["reference"], "-",  "v"),
]


def load() -> dict:
    out = {}
    for key, *_ in SERIES:
        path = SRC / f"{key}.json"
        if path.exists() and path.stat().st_size > 2:
            out[key] = json.loads(path.read_text(encoding="utf-8"))
    if not out:
        raise SystemExit(f"no scan JSON under {SRC}")
    return out


def main() -> int:
    data = load()
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COL, 3.0))

    # ── (a) contamination vs threshold ───────────────────────────────────────
    # Plotted as contamination N_-/N_+ on a log axis rather than purity: it is
    # the quantity ALHAMBRA thresholds on (3 per cent; Molino+2014) and it keeps
    # the well-behaved high-sigma end readable instead of piled against 1.
    ax = axes[0]
    ax.axvspan(0.95, 1.35, color=PALETTE["grey"], alpha=0.16, lw=0)
    ax.axhline(1 - PUR_MIN, color=PALETTE["grey"], ls=":", lw=0.9)
    ax.axvline(NOMINAL, color=PALETTE["black"], ls="-", lw=0.7, alpha=0.55)

    floors = {}
    for key, label, colour, ls, mk in SERIES:
        d = data.get(key)
        if not d:
            continue
        sig = np.array([r["sigma"] for r in d["rows"]], float)
        npos = np.array([r["n_pos"] for r in d["rows"]], float)
        nneg = np.array([r["n_neg"] for r in d["rows"]], float)
        cont = nneg / npos
        # A frame with zero negative detections is an upper limit, not a zero.
        # Upper limits are drawn at 1/N_+ as open carets and kept OUT of the
        # connecting line — joining them would draw a slope that is an artifact
        # of N_+ changing, not of contamination changing.
        is_ul = nneg == 0
        plot_y = np.where(is_ul, 1.0 / npos, cont)
        ax.plot(sig, np.where(is_ul, np.nan, cont), ls=ls, marker=mk,
                color=colour, label=label, ms=3.2, lw=1.1)
        if is_ul.any():
            ax.plot(sig[is_ul], plot_y[is_ul], ls="none", marker="v",
                    color=colour, ms=3.6, mfc="none", mew=0.9, zorder=4)
        ok = sig[cont <= 1 - PUR_MIN]
        if ok.size:
            floors[key] = float(ok.min())
            j = int(np.argmin(np.abs(sig - floors[key])))
            ax.plot(sig[j], plot_y[j], marker=mk, color=colour, ms=8.0,
                    mfc="none", mew=1.3, zorder=5)

    ax.set_yscale("log")
    ax.invert_xaxis()
    ax.set_xlabel(r"detection threshold  $\sigma$")
    ax.set_ylabel(r"estimated contamination  $N_-/N_+$")
    # Room below the lowest upper limit so the legend never overlaps data.
    ax.set_ylim(1.1e-5, 2.2)
    ax.set_xlim(5.3, 0.88)
    ax.text(3.14, 0.42, "pipeline\ndefault", fontsize=6.5, ha="right", va="top",
            color=PALETTE["black"], alpha=0.75)
    ax.text(1.15, 1.9e-5, "collapse", fontsize=6.5, ha="center",
            color=PALETTE["grey"])
    ax.text(5.15, 0.062, "5 per cent gate", fontsize=6.5, color=PALETTE["grey"])
    ax.text(5.15, 2.4e-4, r"open $\triangledown$: $N_-\!=\!0$ (upper limit)",
            fontsize=6.2, color=PALETTE["grey"])
    ax.legend(loc="lower left", fontsize=6.8, ncol=2, handlelength=1.7,
              borderaxespad=0.4, columnspacing=0.8)
    ax.set_title(r"(a) safe floor differs frame to frame", loc="left", fontsize=9)

    # ── (b) estimate vs Gaia-realized ────────────────────────────────────────
    ax = axes[1]
    lim = (0.6, 2.2e4)
    grid = np.array(lim)
    ax.plot(grid, grid, "-", color=PALETTE["black"], lw=0.8, zorder=1)
    ax.fill_between(grid, grid * 0.5, grid * 2.0, color=PALETTE["grey"],
                    alpha=0.15, lw=0, zorder=0)

    n_pt = 0
    for key, label, colour, ls, mk in SERIES:
        d = data.get(key)
        if not d:
            continue
        xs = [r["fp_gaia"] for r in d["rows"] if "fp_gaia" in r]
        ys = [r["n_neg"] for r in d["rows"] if "fp_gaia" in r]
        if not xs:
            continue
        n_pt += len(xs)
        ax.plot(xs, ys, ls="none", marker=mk, color=colour, ms=3.6,
                label=label, alpha=0.9, zorder=3)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("realized spurious, Gaia-unmatched")
    ax.set_ylabel("estimated spurious, sign-flipped")
    ax.text(1.1, 6.0e3, "factor 2", fontsize=7, color=PALETTE["grey"], rotation=41)
    ax.legend(loc="lower right", fontsize=6.8, ncol=1, borderaxespad=0.4)
    ax.set_title(f"(b) estimate validated, {n_pt} points", loc="left",
                 fontsize=9)

    fig.tight_layout(pad=0.5, w_pad=1.8)
    paths = save_fig(fig, "fig_detection_threshold", OUTDIR)
    for ext, p in paths.items():
        print(f"[saved] {p}")
    print("purity>=0.95 floors:", {k: v for k, v in floors.items()})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
