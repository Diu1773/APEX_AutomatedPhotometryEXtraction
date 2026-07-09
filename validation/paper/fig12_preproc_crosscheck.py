"""Figure 12 — Per-step preprocessing cross-check against astropy ccdproc.

Every APEX calibration stage (master bias/dark/flat construction, then bias/dark/
flat application, then the full pipeline) is run alongside the equivalent
operation in ccdproc — the community-standard, independently-maintained astropy
CCD reduction package — on the SAME real Moravian C3-61000 frames (NGC 6811,
B band, 2x2). This is the calibration analogue of the sep (Fig 4) and IRAF
(Fig 5) photometry cross-checks: an independent implementation, so agreement is
evidence the reduction arithmetic is standard and correct.

Panel (a): the per-step APEX-vs-ccdproc pixel disagreement (max |delta|), against
the detector's own read noise and the sky shot noise. Panel (b): the same three
numbers as a noise budget. Every stage agrees to <= 5e-4 DN — four to nine orders
of magnitude below any real noise term — so the two independent pipelines are
numerically identical.

(This quantifies cross-implementation agreement, NOT ground truth; the
ground-truth gate is the synthetic inject->recover test, Fig 10 / calibration
tests. Cosmetic correction uses astroscrappy = the L.A.Cosmic reference itself.)

Run: .venv-deploy\\Scripts\\python validation\\paper\\fig12_preproc_crosscheck.py
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

DATA = REPO / "validation" / "paper" / "data"
OUTDIR = REPO / "validation" / "paper" / "figures"
CAPDIR = REPO / "validation" / "paper" / "captions"

STEP_LABELS = {
    "master_bias": "master bias\n(combine)",
    "master_dark": "master dark\n(combine)",
    "master_flat": "master flat\n(combine+norm)",
    "bias_subtract": "bias subtract",
    "dark_subtract": "dark subtract\n(exp-scaled)",
    "flat_correct": "flat field",
    "full_pipeline": "full pipeline",
}
ORDER = list(STEP_LABELS)
FLOOR = 1e-9   # plot exact-zero (bit-identical) steps at the axis floor


def main() -> int:
    res = json.loads((DATA / "calib_crosscheck_ngc6811.json").read_text())
    ptc = json.loads((DATA / "detector_ptc.json").read_text())
    gain = ptc["gain_e_per_adu"]
    rn_dn = ptc["read_noise_adu"]
    sky_dn = res["light_median"]
    sky_shot_dn = float(np.sqrt(max(sky_dn, 1.0) / gain))   # sqrt(sky_e)/gain

    steps = res["steps"]
    maxabs = np.array([steps[s]["max_abs"] for s in ORDER])
    plotx = np.where(maxabs > 0, maxabs, FLOOR)
    exact = maxabs == 0

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 3.3),
                                   gridspec_kw={"width_ratios": [2.1, 1.0]})

    # (a) per-step lollipop
    y = np.arange(len(ORDER))[::-1]
    axa.hlines(y, FLOOR, plotx, color=PALETTE["grey"], lw=1.2, zorder=1)
    axa.scatter(plotx, y, s=42, color=C["data"], zorder=3)
    for xi, yi, ex in zip(plotx, y, exact):
        if ex:
            axa.annotate("= 0 (bit-identical)", (FLOOR, yi), xytext=(4, 0),
                         textcoords="offset points", va="center", ha="left",
                         fontsize=6.6, color=C["model"])
    axa.axvline(rn_dn, color=C["reference"], lw=1.4, ls="--",
                label=f"read noise ({rn_dn:.2f} DN)")
    axa.axvline(sky_shot_dn, color=C["bad"], lw=1.4, ls=":",
                label=f"sky shot noise ({sky_shot_dn:.0f} DN)")
    axa.set_xscale("log")
    axa.set_xlim(FLOOR, sky_shot_dn * 4)
    axa.set_yticks(y)
    axa.set_yticklabels([STEP_LABELS[s] for s in ORDER], fontsize=7)
    axa.set_xlabel(r"APEX $-$ ccdproc  max$|\Delta|$  (DN)")
    axa.legend(loc="lower right", fontsize=6.8)
    axa.set_title("(a) Every preprocessing step vs ccdproc", loc="left")

    # (b) noise budget
    worst = float(np.max(maxabs))
    names = ["APEX vs\nccdproc", "read\nnoise", "sky shot\nnoise"]
    vals = [max(worst, FLOOR), rn_dn, sky_shot_dn]
    cols = [C["data"], C["reference"], C["bad"]]
    xb = np.arange(3)
    axb.bar(xb, vals, color=cols, width=0.62, zorder=3)
    axb.set_yscale("log")
    axb.set_ylim(1e-5, sky_shot_dn * 4)
    axb.set_xticks(xb)
    axb.set_xticklabels(names, fontsize=7)
    axb.set_ylabel("DN")
    for xi, v in zip(xb, vals):
        axb.annotate(f"{v:.0e}" if v < 1e-2 else f"{v:.1f}", (xi, v),
                     xytext=(0, 3), textcoords="offset points",
                     ha="center", va="bottom", fontsize=6.8)
    axb.set_title("(b) Noise budget", loc="left")

    fig.suptitle(
        "APEX reduction is numerically identical to the community-standard "
        "ccdproc at every stage (NGC 6811, B, C3-61000).",
        fontsize=7.4, y=1.02, color="#333333")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    paths = save_fig(fig, "fig12_preproc_crosscheck", OUTDIR)
    plt.close(fig)

    CAPDIR.mkdir(parents=True, exist_ok=True)
    (CAPDIR / "fig12_preproc_crosscheck.md").write_text(
        f"""# Figure 12 — Per-step preprocessing cross-check vs ccdproc

**Figure 12.** Each APEX detector-calibration stage compared, pixel-for-pixel,
against the equivalent operation in **astropy ccdproc** — the community-standard
Python CCD-reduction package — on the same real Moravian C3-61000 frames
(NGC 6811, B, 2x2; {res['n_bias']} bias, {res['n_dark']} darks, {res['n_flat']}
flats, one science frame). **(a)** The maximum per-pixel disagreement for master
bias/dark/flat construction, bias/dark/flat application, and the full pipeline.
Master bias, master dark, and all three application steps are **bit-identical**
(delta = 0); master-flat construction and the full pipeline agree to
{steps['master_flat']['max_abs']:.0e} and {steps['full_pipeline']['max_abs']:.0e}
DN (float32 rounding). All stages sit four to nine orders of magnitude below the
detector read noise ({rn_dn:.2f} DN) and the sky shot noise
({sky_shot_dn:.0f} DN). **(b)** The same three quantities as a noise budget.
The two independently-written pipelines are numerically identical, so APEX's
reduction implements the standard bias/dark/flat arithmetic correctly. This is a
cross-implementation check (analogous to the sep and IRAF photometry
cross-checks), not a ground-truth validation — the latter is the synthetic
inject->recover test. Cosmetic correction uses astroscrappy, the L.A.Cosmic
reference implementation.
""", encoding="utf-8")

    print("=== fig12 preprocessing cross-check ===")
    for s in ORDER:
        st = steps[s]
        print(f"  {s:16s} max|Δ|={st['max_abs']:.2e} DN  σ={st['robust_sigma']:.2e}")
    print(f"read noise {rn_dn:.2f} DN, sky shot {sky_shot_dn:.1f} DN")
    for ext, p in paths.items():
        print(f"wrote {ext}: {p}  exists={p.exists()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
