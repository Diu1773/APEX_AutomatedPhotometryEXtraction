"""Figure — the built-in quad plate solver against ASTAP and astrometry.net.

Same eight M13 exposures, one copy tree per engine, the same step-4 detection
list shared by all three (the solution is written into the FITS header of the
copy, so nothing overwrites anything). Two comparisons:

(a) agreement between the solutions themselves: a 5x5 pixel grid per frame is
    projected to sky through each engine's WCS and the pairwise angular
    separation is taken (median over the grid).

(b) each solution's Gaia residual, measured independently of the pipeline's QC
    code: Gaia DR3 stars are projected to pixels through the header WCS and
    matched to the step-4 detections within 2 arcsec; the RMS distance is the
    residual. Computed identically for all three engines, no clipping, no
    refinement — so the numbers are comparable across engines.

Inputs: data_wcs_engines/engine_cross.json, independent_residuals.json
(produced by run_wcs_engine_cross.py and the independent-residual pass).

    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig_wcs_engines.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import matplotlib.pyplot as plt

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

DATA = REPO / "validation" / "paper" / "data_wcs_engines"
OUTDIR = REPO / "validation" / "paper" / "figures"
PIX = 0.395   # arcsec / px

PAIRS = [("internal|astap", "built-in vs ASTAP", C["data"], "o"),
         ("internal|astnet", "built-in vs astrometry.net", C["model"], "s"),
         ("astap|astnet", "ASTAP vs astrometry.net", C["reference"], "^")]
ENG = [("internal", "built-in", C["data"], "o"),
       ("astnet", "astrometry.net", C["model"], "s"),
       ("astap", "ASTAP", C["reference"], "^")]


def short(fname: str) -> str:
    # pp_messier13-0001-B.fit -> 1B
    core = fname.split("-")
    return core[1].lstrip("0") + core[2][0]


def main() -> int:
    cross = json.loads((DATA / "engine_cross.json").read_text(encoding="utf-8"))
    resid = json.loads((DATA / "independent_residuals.json").read_text(encoding="utf-8"))

    frames = [r["file"] for r in cross["compare"]["frames"] if "pairs" in r]
    labels = [short(f) for f in frames]
    x = np.arange(len(frames))

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COL, 2.6))

    # (a) pairwise solution separations
    ax = axes[0]
    for key, lab, colour, mk in PAIRS:
        y = [next(r["pairs"][key]["median"] for r in cross["compare"]["frames"]
                  if r["file"] == f) for f in frames]
        ax.plot(x, y, marker=mk, ls="-", ms=3.4, lw=1.0, color=colour, label=lab)
    ax.axhline(PIX, color=PALETTE["grey"], ls=":", lw=0.9)
    ax.text(0.02, PIX + 0.03, "1 px", fontsize=7, color=PALETTE["grey"], ha="left")
    ax.set_xticks(x, labels)
    ax.set_xlabel("frame")
    ax.set_ylabel("solution separation (arcsec)")
    ax.set_ylim(0, 1.25)
    ax.legend(fontsize=6.4, loc="upper right", handlelength=1.6)
    ax.set_title("(a) agreement between the three solutions", loc="left", fontsize=9)

    # (b) independent Gaia residuals
    ax = axes[1]
    for key, lab, colour, mk in ENG:
        rows = {r["file"]: r["rms_px"] for r in resid[key]}
        y = [rows[f] for f in frames]
        ax.plot(x, y, marker=mk, ls="-", ms=3.4, lw=1.0, color=colour, label=lab)
    ax.set_xticks(x, labels)
    ax.set_xlabel("frame")
    ax.set_ylabel("Gaia residual RMS (px)")
    ax.set_ylim(0, 2.3)
    ax.legend(fontsize=6.4, loc="upper right", handlelength=1.6)
    ax.set_title("(b) Gaia residual, one procedure for all", loc="left", fontsize=9)

    fig.tight_layout(pad=0.5, w_pad=1.6)
    for ext, p in save_fig(fig, "fig_wcs_engines", OUTDIR).items():
        print(f"[saved] {p}")

    for key, lab, *_ in PAIRS:
        v = [r["pairs"][key]["median"] for r in cross["compare"]["frames"] if "pairs" in r]
        print(f"{lab:<28} 중앙값 {np.median(v):.3f}\"  범위 {min(v):.3f}-{max(v):.3f}\"")
    for key, lab, *_ in ENG:
        v = [r["rms_px"] for r in resid[key]]
        print(f"{lab:<28} rms 중앙값 {np.median(v):.3f}px = {np.median(v)*PIX:.3f}\"")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
