# -*- coding: utf-8 -*-
"""Figure 4 — absolute detector-preprocessing products.

The user-facing comparison is intentionally a table, not an APEX-minus-reference
plot.  It lists the measured product level and robust spread for APEX, the
independent Python ``ccdproc`` reduction, and the independent IRAF ``ccdproc``
task on the same NGC 6811 calibration set.

Run::

    .venv-deploy\\Scripts\\python -X utf8 validation\\paper\\fig12_preproc_crosscheck.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO / "validation" / "paper"))

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

DATA = REPO / "validation" / "paper" / "data"
OUTDIR = REPO / "validation" / "paper" / "figures"
CAPDIR = REPO / "validation" / "paper" / "captions"

STAGES = [
    ("master_bias", "master bias"),
    ("master_dark", "master dark (60 s)"),
    ("master_flat", "master flat"),
    ("full_pipeline", "science after full chain"),
]
METHODS = [
    ("apex", "APEX"),
    ("python_ccdproc", "Python\nccdproc"),
    ("iraf_ccdproc", "IRAF\nccdproc"),
]


def _value(row: dict) -> str:
    return f"{float(row['median']):.4g} ± {float(row['robust_sigma']):.4g}"


def main() -> int:
    summary = json.loads(
        (DATA / "preproc_absolute_summary.json").read_text(encoding="utf-8")
    )

    fig, ax = plt.subplots(figsize=(DOUBLE_COL, 3.25))
    ax.axis("off")
    ax.set_title(
        "Figure 4 — absolute calibration products (median ± robust σ)",
        loc="left", fontsize=9.3, pad=5,
    )

    headers = ["stage", "APEX", "Python\nccdproc", "IRAF\nccdproc"]
    rows = []
    for key, label in STAGES:
        rows.append([label] + [_value(summary[method][key]) for method, _ in METHODS])

    table = ax.table(
        cellText=rows,
        colLabels=headers,
        cellLoc="center",
        colLoc="center",
        colWidths=[0.31, 0.22, 0.235, 0.235],
        bbox=[0.02, 0.25, 0.96, 0.63],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.0)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor(PALETTE["grey"])
        cell.set_linewidth(0.45)
        cell.PAD = 0.12
        if row == 0:
            cell.set_facecolor("#EAF2F8" if col == 1 else
                               "#EEF7EE" if col == 2 else
                               "#FFF2E5" if col == 3 else "#F2F3F5")
            cell.get_text().set_weight("bold")
        elif col == 1:
            cell.get_text().set_color(C["data"])
        elif col == 2:
            cell.get_text().set_color("#2A7F3F")
        elif col == 3:
            cell.get_text().set_color(PALETTE["orange"])
        if col == 0:
            cell.get_text().set_ha("left")

    ax.text(
        0.02, 0.16,
        "All entries are absolute output values; no APEX-reference subtraction is plotted.",
        transform=ax.transAxes, va="bottom", fontsize=7.2, color=PALETTE["grey"],
    )
    ax.text(
        0.02, 0.105,
        "Master flat is unit-median normalized. Robust σ = 1.4826 × MAD over the full frame.",
        transform=ax.transAxes, va="bottom", fontsize=7.2, color=PALETTE["grey"],
    )
    ax.text(
        0.02, 0.05,
        f"NGC 6811 B 60 s · Moravian C3-61000 (2×2) · 8 bias / 8 dark / 5 flat · {summary['input'].split('(')[-1].rstrip(')')}",
        transform=ax.transAxes, va="bottom", fontsize=6.7, color=PALETTE["grey"],
    )

    fig.tight_layout(rect=(0, 0.02, 1, 0.98))
    paths = save_fig(fig, "fig12_preproc_crosscheck", OUTDIR)
    plt.close(fig)

    CAPDIR.mkdir(parents=True, exist_ok=True)
    (CAPDIR / "fig12_preproc_crosscheck.md").write_text(
        """# Figure — absolute detector-preprocessing products

The table lists the absolute output level of each calibration product as the
full-frame median ± robust σ (1.4826 × MAD), rather than subtracting APEX from
either reference.  Columns are the production APEX reduction, an independent
Python `ccdproc` reduction, and an independent IRAF `ccdproc` task run through
PyRAF.  Rows are the master bias, the 60-s master dark, the unit-median
normalised master flat, and the science frame after the complete bias–dark–flat
chain.  The Python column agrees with APEX at the displayed precision; the
separate pixel-residual audit remains in `data/iraf_preproc_stats.json` and is
not used as the plotted quantity.  Inputs are eight bias frames, eight 60-s
darks, five B flats, and one 60-s NGC 6811 B science frame from 2026-06-11
(Moravian C3-61000, 2×2).  Cosmetic repair was disabled because this table
isolates detector calibration arithmetic.
""",
        encoding="utf-8",
    )

    print("=== absolute detector preprocessing products ===")
    for key, label in STAGES:
        values = " | ".join(
            f"{name}={_value(summary[method][key])}"
            for method, name in METHODS
        )
        print(f"  {label:25s} {values}")
    for ext, path in paths.items():
        print(f"wrote {ext}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
