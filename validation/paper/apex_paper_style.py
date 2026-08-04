"""Shared publication-quality matplotlib style for APEX validation figures.

Every APEX paper figure imports this module so the whole figure set is visually
consistent: one serif font family, one colour-blind-safe palette, one
line/grid weight scale, and one save path that emits BOTH a
vector PDF (for LaTeX) and a 300-dpi PNG (for quick viewing / the docs site).

Usage
-----
    from apex_paper_style import apply_paper_style, save_fig, PALETTE, C
    apply_paper_style()
    fig, ax = plt.subplots(figsize=(3.4, 2.8))   # single-column width
    ...
    save_fig(fig, "fig1_completeness", outdir)     # -> .pdf + .png

Design notes
------------
* Figure widths follow journal columns: SINGLE_COL = 3.4 in (~86 mm, one
  column of a two-column article); DOUBLE_COL = 7.0 in (full text width).
* Semantic colours use the Okabe–Ito colour-blind-safe set.  Markers, line
  styles, and hatches still distinguish categories in grayscale printouts.
* 300 dpi PNG + true vector PDF is the standard "camera-ready" combination.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: never require a display
import matplotlib.pyplot as plt


# ── Journal monochrome palette ─────────────────────────────────────────────
PALETTE: dict[str, str] = {
    "black": "#000000",
    "orange": "#E69F00",
    "skyblue": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "grey": "#777777",
}

# Semantic aliases (use these in figures so intent is explicit).
C = {
    "data": PALETTE["blue"],        # measured / recovered points (dark)
    "model": PALETTE["vermillion"], # fitted model / theory line (mid)
    "reference": PALETTE["green"],  # external tool / truth reference
    "accent": PALETTE["purple"],    # secondary series (light)
    "floor": PALETTE["grey"],       # noise floor / guide lines
    "bad": PALETTE["vermillion"],
    "good": PALETTE["green"],
}

# Journal column widths in inches.
SINGLE_COL = 3.4
DOUBLE_COL = 7.0


def apply_paper_style() -> None:
    """Install the shared rcParams. Call once at the top of every figure script."""
    plt.rcParams.update(
        {
            # Fonts — serif body, Computer-Modern-style math to match LaTeX.
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman", "Times", "serif"],
            "mathtext.fontset": "cm",
            "font.size": 9.0,
            "axes.titlesize": 9.5,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "legend.fontsize": 8.0,
            "figure.titlesize": 10.0,
            # Lines & markers.
            "lines.linewidth": 1.4,
            "lines.markersize": 4.0,
            "lines.markeredgewidth": 0.6,
            # Axes & spines.
            "axes.linewidth": 0.8,
            "axes.grid": True,
            "axes.axisbelow": True,
            "axes.prop_cycle": plt.cycler(
                color=[
                    PALETTE["blue"],
                    PALETTE["vermillion"],
                    PALETTE["green"],
                    PALETTE["purple"],
                    PALETTE["grey"],
                    PALETTE["skyblue"],
                    PALETTE["black"],
                ]
            ),
            # Grid.
            "grid.color": "#B0B0B0",
            "grid.linewidth": 0.5,
            "grid.alpha": 0.35,
            # Ticks — inward, on all four sides (journal convention).
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 3.5,
            "ytick.major.size": 3.5,
            "xtick.minor.size": 2.0,
            "ytick.minor.size": 2.0,
            "xtick.minor.visible": True,
            "ytick.minor.visible": True,
            # Legend.
            "legend.frameon": False,
            "legend.handlelength": 1.6,
            "legend.columnspacing": 1.0,
            # Figure / save.
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
            "pdf.fonttype": 42,   # embed TrueType (editable text in the PDF)
            "ps.fonttype": 42,
        }
    )


def save_fig(fig, stem: str, outdir: str | Path) -> dict[str, Path]:
    """Save ``fig`` as both a vector PDF and a 300-dpi PNG under ``outdir``.

    Returns the two written paths keyed by extension.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for ext in ("pdf", "png"):
        path = outdir / f"{stem}.{ext}"
        fig.savefig(path)
        paths[ext] = path
    return paths
