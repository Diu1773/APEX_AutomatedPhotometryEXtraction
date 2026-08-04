"""Figure 3 — detector constants reproduced by three reduction paths.

The figure is deliberately a cross-tool comparison, not a comparison with a
camera catalogue value.  APEX, Python ``ccdproc`` and the IRAF ``ccdproc``
task use the same 2026-06-11 calibration set and the same estimators: a
flat-pair photon-transfer fit, a bias-pair read-noise estimate, and a linear
dark-current ladder.  The pixel-level residuals of the reductions are shown
separately in Figure 4; here the reader can see the actual constants that are
fed to the photometric error model.

Run::

    .venv-deploy\\Scripts\\python -X utf8 validation\\paper\\fig11_detector.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO / "validation" / "paper"))

from apex_paper_style import C, DOUBLE_COL, PALETTE, apply_paper_style, save_fig

apply_paper_style()

DATA = REPO / "validation" / "paper" / "data"
OUTDIR = REPO / "validation" / "paper" / "figures"
CAPDIR = REPO / "validation" / "paper" / "captions"

METHODS = ["APEX", "Python ccdproc", "IRAF ccdproc"]
METHOD_COLORS = [C["data"], PALETTE["orange"], PALETTE["green"]]


def _fmt(value: float, key: str) -> str:
    if key == "gain":
        return f"{value:.3f}"
    if key == "read_noise":
        return f"{value:.2f}"
    return f"{value:.4f}"


def main() -> int:
    payload = json.loads((DATA / "detector_crosscheck.json").read_text(encoding="utf-8"))
    metrics = payload["metrics"]

    fig = plt.figure(figsize=(DOUBLE_COL, 3.35))
    grid = fig.add_gridspec(
        1, 2, width_ratios=[1.28, 1.0], wspace=0.34,
        left=0.06, right=0.98, top=0.85, bottom=0.22,
    )

    # ── (a) the values themselves ────────────────────────────────────────
    ax_table = fig.add_subplot(grid[0, 0])
    ax_table.axis("off")
    ax_table.set_title("(a) detector constants", loc="left", fontsize=8.8, pad=3)
    headers = ["quantity", "unit", "APEX", "Python\nccdproc", "IRAF\nccdproc"]
    rows = []
    for metric in metrics:
        vals = metric["values"]
        errs = metric["errors"]
        cells = [metric["label"], metric["unit"]]
        for method in METHODS:
            err = float(errs[method])
            val = float(vals[method])
            cells.append(f"{_fmt(val, metric['key'])}"
                         + (f" ± {_fmt(err, metric['key'])}" if err > 0 else ""))
        rows.append(cells)
    table = ax_table.table(
        cellText=rows, colLabels=headers, cellLoc="center", colLoc="center",
        colWidths=[0.24, 0.13, 0.17, 0.23, 0.23], bbox=[0.0, 0.23, 1.0, 0.62],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(5.9)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor(PALETTE["grey"])
        cell.set_linewidth(0.4)
        if row == 0:
            cell.get_text().set_weight("bold")
            if col >= 2:
                cell.get_text().set_color(METHOD_COLORS[col - 2])
        elif col >= 2:
            cell.get_text().set_color(METHOD_COLORS[col - 2])
        if col == 0:
            cell.get_text().set_ha("left")
    ax_table.text(
        0.0, 0.13,
        "same raw calibration set · same PTC and dark-ladder estimators",
        transform=ax_table.transAxes, fontsize=6.4, color=PALETTE["grey"],
    )
    ax_table.text(
        0.0, 0.05,
        "pixel-level residuals are given in Figure 4; no header/vendor value is used",
        transform=ax_table.transAxes, fontsize=6.4, color=PALETTE["grey"],
    )

    # ── (b) grouped bars, one scale per physical quantity ────────────────
    sub = grid[0, 1].subgridspec(3, 1, hspace=0.78)
    for idx, metric in enumerate(metrics):
        ax = fig.add_subplot(sub[idx, 0])
        vals = np.array([float(metric["values"][m]) for m in METHODS])
        errs = np.array([float(metric["errors"][m]) for m in METHODS])
        x = np.arange(len(METHODS))
        ax.bar(x, vals, yerr=errs if np.any(errs > 0) else None,
               color=METHOD_COLORS, width=0.68, capsize=2.5,
               edgecolor="white", linewidth=0.35, zorder=3)
        ax.set_ylabel(metric["unit"], fontsize=6.9, labelpad=2)
        ax.set_title(metric["label"], loc="left", fontsize=7.5, pad=1)
        ax.grid(axis="y", alpha=0.25)
        ax.set_xticks(x)
        if idx == len(metrics) - 1:
            ax.set_xticklabels(["APEX", "Python", "IRAF"], fontsize=6.4)
        else:
            ax.set_xticklabels([])
        lo = max(0.0, float(vals.min()) * 0.88)
        hi = float(vals.max() + max(errs.max(), vals.max() * 0.04)) * 1.12
        ax.set_ylim(lo, hi)
        ax.tick_params(axis="y", labelsize=6.4, pad=1)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    fig.suptitle(
        "Moravian C3-61000 (Sony IMX455, 2×2): constants reproduced by three reductions",
        fontsize=8.3, y=0.96, color="#333333",
    )
    fig.text(
        0.98, 0.055,
        "8 bias · 5 B flats · 10–480 s dark ladder · 2026-06-11",
        ha="right", va="bottom", fontsize=6.1, color=PALETTE["grey"],
    )

    paths = save_fig(fig, "fig11_detector", OUTDIR)
    plt.close(fig)

    CAPDIR.mkdir(parents=True, exist_ok=True)
    gain = metrics[0]["values"]["APEX"]
    gerr = metrics[0]["errors"]["APEX"]
    rn = metrics[1]["values"]["APEX"]
    dark = metrics[2]["values"]["APEX"]
    (CAPDIR / "fig11_detector.md").write_text(
        f"""# Figure 3 — detector constants from three reductions

**Figure 3.** Detector constants obtained from the same Moravian C3-61000
(Sony IMX455, 2×2) calibration set by three reduction paths: APEX, the Python
`ccdproc` package, and the IRAF `ccdproc` task. **(a)** The table reports the
values used by the error model: flat-pair photon-transfer gain
({gain:.3f} ± {gerr:.3f} e⁻/ADU), bias-pair read noise ({rn:.2f} e⁻), and the
linear dark-ladder slope ({dark:.4f} e⁻/s). **(b)** The grouped bars show one
physical scale per quantity, so the agreement is visible without putting
incommensurate units on a single axis. The three values are identical at the
shown precision because bias/dark subtraction is additive and cancels in the
flat-pair variance, while all reductions use the same median dark ladder; this
is an agreement result, not a claim that the pixel arrays are byte-identical.
The pixel-level residuals, including the independent IRAF flat-normalisation
and full-chain differences, are reported in Figure 4. No FITS-header, vendor,
or laboratory gain is used in this comparison. Inputs: 8 bias frames, 5 B
flats, and the 10–480 s dark ladder from 2026-06-11.
""",
        encoding="utf-8",
    )

    print("=== fig11 detector cross-check ===")
    for metric in metrics:
        values = ", ".join(f"{m}={metric['values'][m]:.6g}" for m in METHODS)
        print(f"{metric['label']:12s}: {values}")
    for ext, path in paths.items():
        print(f"wrote {ext}: {path}  exists={path.exists()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
