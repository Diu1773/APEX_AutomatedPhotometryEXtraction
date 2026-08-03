# -*- coding: utf-8 -*-
r"""Per-stage calibration cross-check against astropy ccdproc.

Every calibration stage remains an explicit datum. Bit-identical stages are
shown as labelled zero markers rather than empty image panels, and the sole
non-zero end-to-end difference is compared with detector and sky noise.

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
        1, 2, figsize=(DOUBLE_COL, 2.75),
        gridspec_kw={"width_ratios": [1.55, 1.0]},
    )

    # Exact zeros are placed at a display floor. Their labels retain the
    # measured value so the visual position cannot be read as a measurement.
    x = list(range(len(STEPS)))
    display_floor = 1e-6
    y = [max(steps[key]["max_abs"], display_floor) for key, _ in STEPS]
    axa.vlines(x, display_floor, y, color=PALETTE["grey"], lw=0.9, zorder=2)
    for xi, ((key, _label), yi) in enumerate(zip(STEPS, y)):
        exact_zero = steps[key]["max_abs"] == 0.0
        axa.plot(
            xi, yi, marker="o" if exact_zero else "s", ms=5.0,
            mfc="white" if exact_zero else C["data"], mec=C["data"],
            ls="none", zorder=3,
        )
        note = "0" if exact_zero else f"{steps[key]['max_abs']:.1e}"
        axa.annotate(
            note, (xi, yi), xytext=(0, 6), textcoords="offset points",
            ha="center", va="bottom", fontsize=6.3,
        )
    axa.set_yscale("log")
    axa.set_ylim(5e-7, 4e-3)
    axa.set_xticks(x)
    axa.set_xticklabels(
        [label.replace(" ", "\n", 1) for _, label in STEPS], fontsize=6.6,
    )
    axa.set_ylabel(r"max $|\mathrm{APEX}-\mathrm{ccdproc}|$ (DN)")
    axa.set_title("(a) six stages are bit-identical", loc="left")
    axa.text(
        0.01, 0.95,
        "six open circles = measured zero (shown at a display floor)",
        transform=axa.transAxes, va="top", fontsize=5.8, color=PALETTE["grey"],
    )

    names = ["APEX vs\nccdproc", "read\nnoise", "sky shot\nnoise"]
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
        f"""# Figure — per-step preprocessing cross-check vs ccdproc

**(a)** Maximum absolute APEX−ccdproc difference at every calibration stage.
Master bias/dark/flat construction and the three individual corrections are
bit-identical; open circles are labelled measured zeros and placed on a display
floor only so that all stages remain visible. The full chain differs by at most
{steps['full_pipeline']['max_abs']:.1e} DN (robust σ =
{steps['full_pipeline']['robust_sigma']:.1e} DN) because of float32 rounding.
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
