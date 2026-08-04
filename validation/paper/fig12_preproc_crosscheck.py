# -*- coding: utf-8 -*-
r"""Per-stage calibration cross-check against the Python ccdproc package.

Every calibration stage remains an explicit datum. The left panel is a compact
numeric audit table (maximum absolute difference, robust scatter, and status),
while the right panel compares the only non-zero end-to-end difference with
detector and sky noise. This is a pixel-level arithmetic test, not only a chart
comparison.

Run: .venv-deploy\Scripts\python -X utf8 validation\paper\fig12_preproc_crosscheck.py
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

STEPS = [
    ("master_bias", "master bias"),
    ("master_dark", "master dark"),
    ("master_flat", "master flat"),
    ("bias_subtract", "bias subtract"),
    ("dark_subtract", "dark subtract"),
    ("flat_correct", "flat field"),
    ("full_pipeline", "full pipeline"),
]


def main() -> int:
    res = json.loads((DATA / "calib_crosscheck_ngc6811.json").read_text(encoding="utf-8"))
    ptc = json.loads((DATA / "detector_ptc.json").read_text(encoding="utf-8"))
    gain = ptc["gain_e_per_adu"]
    rn_dn = ptc["read_noise_adu"]
    sky_shot_dn = float((max(res["light_median"], 1.0) / gain) ** 0.5)
    steps = res["steps"]
    worst = max(steps[key]["max_abs"] for key, _ in STEPS)

    fig, (axa, axb) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL, 2.95),
        gridspec_kw={"width_ratios": [1.55, 1.0]},
    )

    # A table is more honest than plotting six exact zeros at an arbitrary
    # logarithmic display floor. It retains both measured metrics and makes
    # the full-chain rounding residual explicit.
    axa.axis("off")
    axa.set_title("(a) stage-wise pixel audit", loc="left", pad=4)
    headers = ["stage", "max |Δ|\n(DN)", "robust σ\n(DN)", "result"]
    rows = []
    for key, label in STEPS:
        mx = float(steps[key]["max_abs"])
        rs = float(steps[key]["robust_sigma"])
        if mx == 0.0:
            mx_s, rs_s, status = "0", "0", "bit exact"
        else:
            mx_s, rs_s, status = f"{mx:.2e}", f"{rs:.2e}", "float32 round"
        rows.append([label, mx_s, rs_s, status])
    table = axa.table(
        cellText=rows, colLabels=headers, cellLoc="center", colLoc="center",
        colWidths=[0.35, 0.22, 0.22, 0.21], bbox=[0.0, 0.08, 1.0, 0.78],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(6.6)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#777777")
        cell.set_linewidth(0.35)
        if row == 0:
            cell.set_text_props(weight="bold")
            cell.set_facecolor("#eeeeee")
        elif col == 0:
            cell.get_text().set_ha("left")
    axa.text(
        0.0, 0.01,
        "master construction + bias/dark/flat application; zero means every pixel matched",
        transform=axa.transAxes, va="bottom", fontsize=5.8, color=PALETTE["grey"],
    )

    names = ["APEX vs Python\nccdproc", "read\nnoise", "sky shot\nnoise"]
    vals = [max(worst, 1e-5), rn_dn, sky_shot_dn]
    xb = list(range(3))
    bars = axb.bar(
        xb, vals, color=["#252525", "#777777", "#c7c7c7"],
        edgecolor="#252525", width=0.62, zorder=3,
    )
    for bar, hatch in zip(bars, ["", "///", "..."]):
        bar.set_hatch(hatch)
    axb.set_yscale("log")
    axb.set_ylim(1e-5, sky_shot_dn * 30)
    axb.set_xticks(xb)
    axb.set_xticklabels(names, fontsize=6.6)
    for xi, value in zip(xb, vals):
        axb.annotate(
            f"{value:.0e}" if value < 1e-2 else f"{value:.1f}",
            (xi, value), xytext=(0, 3), textcoords="offset points",
            ha="center", va="bottom", fontsize=6.6,
        )
    axb.set_ylabel("scale (DN)")
    axb.set_title("(b) below detector and sky noise", loc="left")

    prov1 = (
        f"NGC 6811 $B$ single 60 s · Moravian C3-61000 (2$\\times$2) · "
        f"night {res['night']} · {res['n_bias']} bias / {res['n_dark']} dark "
        f"(60 s) / {res['n_flat']} flat"
    )
    prov2 = (
        f"APEX vs {res['reference']} · cosmetic step off (validated by injection) "
        "· generator: calib_crosscheck_ngc6811.py"
    )
    fig.tight_layout(rect=(0, 0.17, 1, 1), w_pad=1.2)
    fig.text(0.995, 0.078, prov1, ha="right", va="bottom", fontsize=5.7,
             color=PALETTE["grey"])
    fig.text(0.995, 0.020, prov2, ha="right", va="bottom", fontsize=5.7,
             color=PALETTE["grey"])

    paths = save_fig(fig, "fig12_preproc_crosscheck", OUTDIR)
    plt.close(fig)

    CAPDIR.mkdir(parents=True, exist_ok=True)
    (CAPDIR / "fig12_preproc_crosscheck.md").write_text(
        f"""# Figure — per-step preprocessing cross-check vs Python ccdproc

**(a)** Numeric audit against the independent Python `ccdproc` package
(`ccdproc` is not the IRAF task of the same name). The table reports the
maximum absolute difference and robust scatter for each stage; six rows are
bit-identical at every pixel, while the full chain leaves only
{steps['full_pipeline']['max_abs']:.1e} DN (robust σ =
{steps['full_pipeline']['robust_sigma']:.1e} DN) from float32 rounding.
**(b)** The end-to-end difference beside the read-noise ({rn_dn:.2f} DN) and
sky-shot-noise ({sky_shot_dn:.0f} DN) scales. Inputs: {res['n_bias']} bias,
{res['n_dark']} darks (60 s), {res['n_flat']} flats, one 60 s NGC 6811 $B$
light (Moravian C3-61000, 2×2, night {res['night']}); reference
{res['reference']}. The cosmetic (L.A.Cosmic + hot-pixel) stage is disabled
here because it repairs pixels by design and is validated separately by
injection. Generator: `calib_crosscheck_ngc6811.py`.
""",
        encoding="utf-8",
    )

    print("=== preprocessing cross-check ===")
    for key, _ in STEPS:
        print(f"  {key:14s} max|delta|={steps[key]['max_abs']:.2e}")
    for ext, path in paths.items():
        print(f"wrote {ext}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
