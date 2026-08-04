# -*- coding: utf-8 -*-
r"""Figure 4 — detector preprocessing cross-checks.

The comparison is deliberately made at one level: pixel values after the
detector-calibration stages.  APEX is the production reduction; the two
independent reductions are Python ``ccdproc`` and the IRAF ``ccdproc`` task.
The table reports measured differences, while the bar panel gives the same
maximum absolute differences on a log scale.  Exact zeros are printed as
zeros in the table and placed at a labelled display floor in the plot; they
are not replaced by an arbitrary numerical value in the reported result.

Run::

    .venv-deploy\Scripts\python -X utf8 validation\paper\fig12_preproc_crosscheck.py

The small IRAF summary is generated in WSL/PyRAF from the same raw frames and
is stored in ``data/iraf_preproc_stats.json`` (the large FITS intermediates
remain ignored by the repository).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO / "validation" / "paper"))

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

DATA = REPO / "validation" / "paper" / "data"
OUTDIR = REPO / "validation" / "paper" / "figures"
CAPDIR = REPO / "validation" / "paper" / "captions"

STAGES = [
    ("master_bias", "master bias"),
    ("master_dark", "master dark"),
    ("master_flat", "master flat"),
    ("full_pipeline", "full pipeline"),
]


def _fmt(value: float) -> str:
    return "0" if value == 0 else f"{value:.1e}"


def main() -> int:
    py = json.loads((DATA / "calib_crosscheck_ngc6811.json").read_text(encoding="utf-8"))
    iraf = json.loads((DATA / "iraf_preproc_stats.json").read_text(encoding="utf-8"))
    py_steps = py["steps"]
    # The figure compares each independent implementation to the same APEX
    # output, so the values are directly comparable even though the two
    # references are not compared to one another.
    py_vals = [float(py_steps[key]["max_abs"]) for key, _ in STAGES]
    iraf_vals = [float(iraf[key]["max_abs"]) for key, _ in STAGES]
    py_sig = [float(py_steps[key]["robust_sigma"]) for key, _ in STAGES]
    iraf_sig = [float(iraf[key]["robust_sigma"]) for key, _ in STAGES]

    fig, (axa, axb) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL, 3.05),
        gridspec_kw={"width_ratios": [1.42, 1.08]},
    )

    # ── (a) compact numeric table ────────────────────────────────────────
    axa.axis("off")
    axa.set_title("(a) pixel-level comparison", loc="left", fontsize=8.6, pad=3)
    headers = ["stage", "Python\nmax |Δ| / σ", "IRAF\nmax |Δ| / σ", "interpretation"]
    rows = []
    for (key, label), p, ps, i, is_ in zip(STAGES, py_vals, py_sig, iraf_vals, iraf_sig):
        ptxt = f"{_fmt(p)}\n{_fmt(ps)}"
        itxt = f"{_fmt(i)}\n{_fmt(is_)}"
        if p == 0 and i == 0:
            note = "bit exact"
        elif key == "master_flat":
            note = "flat norm."
        else:
            note = "chain arithmetic"
        rows.append([label, ptxt, itxt, note])
    table = axa.table(
        cellText=rows, colLabels=headers, cellLoc="center", colLoc="center",
        colWidths=[0.25, 0.23, 0.23, 0.29], bbox=[0.0, 0.18, 1.0, 0.70],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(6.35)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor(PALETTE["grey"])
        cell.set_linewidth(0.35)
        if row == 0:
            cell.set_facecolor("#EAF2F8") if col == 1 else None
            cell.set_facecolor("#FFF2E5") if col == 2 else cell.get_facecolor()
            cell.get_text().set_weight("bold")
        elif col == 1:
            cell.get_text().set_color(C["data"])
        elif col == 2:
            cell.get_text().set_color(PALETTE["orange"])
        if col == 0:
            cell.get_text().set_ha("left")
    axa.text(
        0.0, 0.08,
        "entries are max |APEX − reference| / robust σ; 0 means every pixel matched",
        transform=axa.transAxes, va="bottom", fontsize=5.8, color=PALETTE["grey"],
    )
    axa.text(
        0.0, 0.02,
        "IRAF flat was divided by its median before comparison (same convention as APEX)",
        transform=axa.transAxes, va="bottom", fontsize=5.8, color=PALETTE["grey"],
    )

    # ── (b) grouped bars ─────────────────────────────────────────────────
    axb.set_title("(b) maximum pixel difference", loc="left", fontsize=8.6, pad=3)
    x = np.arange(len(STAGES), dtype=float)
    width = 0.34
    floor = 1e-6
    py_plot = np.maximum(py_vals, floor)
    iraf_plot = np.maximum(iraf_vals, floor)
    bars_py = axb.bar(x - width / 2, py_plot, width, label="APEX − Python ccdproc",
                      color=C["data"], edgecolor=C["data"], zorder=3)
    bars_iraf = axb.bar(x + width / 2, iraf_plot, width, label="APEX − IRAF ccdproc",
                        color=PALETTE["orange"], edgecolor=PALETTE["orange"],
                        hatch="///", zorder=3)
    axb.set_yscale("log")
    axb.set_ylim(floor, 100)
    axb.set_ylabel("max |Δ| (DN)")
    axb.set_xticks(x)
    axb.set_xticklabels([label.replace(" ", "\n", 1) for _, label in STAGES], fontsize=6.5)
    axb.legend(loc="upper left", fontsize=5.9, handlelength=1.3,
               borderaxespad=0.2, ncol=1)
    axb.axhline(3.5, color=PALETTE["grey"], ls=":", lw=0.8, zorder=2)
    axb.axhline(41.0, color=PALETTE["grey"], ls="--", lw=0.8, zorder=2)
    for bars, vals in ((bars_py, py_vals), (bars_iraf, iraf_vals)):
        for bar, value in zip(bars, vals):
            if value == 0:
                continue
            axb.annotate(f"{value:.1e}",
                         (bar.get_x() + bar.get_width() / 2, value),
                         xytext=(0, 3), textcoords="offset points", ha="center",
                         va="bottom", fontsize=5.8)

    prov1 = (
        f"NGC 6811 B 60 s · Moravian C3-61000 (2×2) · {py['night']} · "
        f"8 bias / 8 dark (60 s) / 5 flat"
    )
    prov2 = (
        "same raw frames and median-combine convention · IRAF ccdred/PyRAF · "
        "cosmetic repair disabled"
    )
    fig.tight_layout(rect=(0, 0.17, 1, 1), w_pad=1.25)
    fig.text(0.995, 0.077, prov1, ha="right", va="bottom", fontsize=5.7,
             color=PALETTE["grey"])
    fig.text(0.995, 0.020, prov2, ha="right", va="bottom", fontsize=5.7,
             color=PALETTE["grey"])

    paths = save_fig(fig, "fig12_preproc_crosscheck", OUTDIR)
    plt.close(fig)

    CAPDIR.mkdir(parents=True, exist_ok=True)
    (CAPDIR / "fig12_preproc_crosscheck.md").write_text(
        f"""# Figure — detector preprocessing comparison

**(a)** Pixel-level comparison of the APEX reduction with two independent
implementations applied to the same NGC 6811 raw set: the Python `ccdproc`
package and the IRAF `ccdproc` task (run through PyRAF). Each entry is
`max |APEX − reference| / robust σ` in DN. The bias and dark masters are bit
identical for both references. The flat comparison uses the same unit-median
normalisation; the remaining {iraf['master_flat']['max_abs']:.2e} DN maximum
is below the displayed detector scales. **(b)** The maximum differences are
shown as grouped bars; exact zeros are placed at a display floor of
{floor:.0e} DN, while the table retains the measured zero. Dotted and dashed
lines mark the 3.5 DN read noise and 41 DN sky-shot-noise scales. The full
pipeline difference is {py_vals[-1]:.2e} DN for Python `ccdproc` and
{iraf_vals[-1]:.2f} DN for IRAF `ccdproc`; the latter reflects the independent
IRAF combination/flat-correction path, not a photometry comparison. Inputs:
8 bias, 8 darks (60 s), 5 flats and one 60 s B-band light, Moravian C3-61000,
night {py['night']}. Cosmetic repair was disabled because it intentionally
changes pixels and is validated separately.
""",
        encoding="utf-8",
    )

    print("=== detector preprocessing comparison ===")
    for (_, label), p, i in zip(STAGES, py_vals, iraf_vals):
        print(f"  {label:14s} Python={p:.3e}  IRAF={i:.3e}")
    for ext, path in paths.items():
        print(f"wrote {ext}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
