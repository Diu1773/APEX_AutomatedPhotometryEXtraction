# -*- coding: utf-8 -*-
"""Figure — per-step preprocessing cross-check against astropy ccdproc.

2026-08-03 rework (user: the old panel (a) log-lollipop put five bit-identical
steps on a fake 1e-9 floor and read as an empty plot — "각 처리단계마다 plot을
넣던가").  Panel (a) is now one difference *map* per calibration step: the
max-pooled |APEX - ccdproc| image itself.  Six steps are exactly zero
everywhere — a blank map is the honest picture of bit-identity — and the full
end-to-end chain shows only float32 rounding dust, orders of magnitude below
read noise.  Panel (b) keeps the noise budget.

Data: data/calib_crosscheck_ngc6811.json + data/calib_crosscheck_maps.npz,
both written by calib_crosscheck_ngc6811.py (the committed generator; it also
records every input frame name).

Run: .venv-deploy\\Scripts\\python -X utf8 validation\\paper\\fig12_preproc_crosscheck.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

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
    maps = np.load(DATA / "calib_crosscheck_maps.npz")
    gain = ptc["gain_e_per_adu"]
    rn_dn = ptc["read_noise_adu"]
    sky_shot_dn = float(np.sqrt(max(res["light_median"], 1.0) / gain))
    steps = res["steps"]
    worst = max(steps[k]["max_abs"] for k, _ in STEPS)

    fig = plt.figure(figsize=(DOUBLE_COL, 4.15))
    gs = fig.add_gridspec(2, 4, hspace=0.42, wspace=0.14,
                          left=0.035, right=0.985, top=0.90, bottom=0.145)

    vmax = 1e-3
    norm = colors.Normalize(vmin=0.0, vmax=vmax)
    last_im = None
    for i, (key, label) in enumerate(STEPS):
        ax = fig.add_subplot(gs[i // 4, i % 4])
        m = maps[key]
        last_im = ax.imshow(m, cmap="Greys", norm=norm, origin="lower",
                            aspect="auto", interpolation="nearest")
        mx = steps[key]["max_abs"]
        ax.set_title(f"({chr(97 + i)}) {label}", loc="left", fontsize=7.6)
        if mx == 0.0:
            note, colr = "$\\Delta = 0$  (bit-identical)", C["model"]
        else:
            note, colr = f"max$|\\Delta|$ = {mx:.1e} DN", C["data"]
        ax.text(0.5, 0.5, note, transform=ax.transAxes, ha="center", va="center",
                fontsize=7.2, color=colr,
                bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.4))
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_linewidth(0.5)

    # (h) noise budget — the same disagreement against the real noise terms
    axb = fig.add_subplot(gs[1, 3])
    names = ["APEX vs\nccdproc", "read\nnoise", "sky shot\nnoise"]
    vals = [max(worst, 1e-5), rn_dn, sky_shot_dn]
    cols = [C["data"], C["reference"], C["bad"]]
    xb = np.arange(3)
    axb.bar(xb, vals, color=cols, width=0.62, zorder=3)
    axb.set_yscale("log")
    axb.set_ylim(1e-5, sky_shot_dn * 30)   # 막대 위 주석이 잘리지 않게
    axb.set_xticks(xb)
    axb.set_xticklabels(names, fontsize=6.6)
    axb.tick_params(labelsize=6.4)
    for xi, v in zip(xb, vals):
        axb.annotate(f"{v:.0e}" if v < 1e-2 else f"{v:.1f}", (xi, v),
                     xytext=(0, 2.5), textcoords="offset points",
                     ha="center", va="bottom", fontsize=6.6)
    axb.set_title("(h) noise budget (DN)", loc="left", fontsize=7.6)

    # 공유 색막대 — 지도 (a)-(g) 의 |Δ| 눈금
    cax = fig.add_axes([0.035, 0.078, 0.30, 0.02])
    cb = fig.colorbar(last_im, cax=cax, orientation="horizontal")
    cb.set_label(r"$|\mathrm{APEX} - \mathrm{ccdproc}|$  (DN)", fontsize=6.4)
    cb.ax.tick_params(labelsize=6.0)

    # 데이터 명세 (그림 안에 출처를 박는다). 한 줄로 쓰면 그림 폭을 넘어
    # tight-bbox 가 캔버스를 왼쪽으로 늘린다(2026-08-03 실측) — 두 줄로 나눈다.
    prov1 = (f"NGC 6811 $B$ single 60 s · Moravian C3-61000 (2$\\times$2) · night {res['night']} · "
             f"{res['n_bias']} bias / {res['n_dark']} dark (60 s) / {res['n_flat']} flat")
    prov2 = (f"APEX vs {res['reference']} · cosmetic step off (validated separately by "
             f"injection) · gen: calib_crosscheck_ngc6811.py")
    fig.text(0.985, 0.052, prov1, ha="right", va="bottom", fontsize=5.9,
             color=PALETTE["grey"])
    fig.text(0.985, 0.012, prov2, ha="right", va="bottom", fontsize=5.9,
             color=PALETTE["grey"])

    paths = save_fig(fig, "fig12_preproc_crosscheck", OUTDIR)
    plt.close(fig)

    CAPDIR.mkdir(parents=True, exist_ok=True)
    (CAPDIR / "fig12_preproc_crosscheck.md").write_text(
        f"""# Figure — per-step preprocessing cross-check vs ccdproc

**(a)–(g)** |APEX − ccdproc| difference map for every calibration stage
(8×8 max-pooled over the full {maps['shape'][0]}×{maps['shape'][1]} frame):
master bias/dark/flat construction, bias/dark/flat application, and the full
end-to-end pipeline. Six stages are **bit-identical** (Δ = 0 at every pixel);
only the full chain shows float32 rounding at max|Δ| =
{steps['full_pipeline']['max_abs']:.1e} DN (robust σ =
{steps['full_pipeline']['robust_sigma']:.1e} DN). **(h)** That worst
disagreement against the detector read noise ({rn_dn:.2f} DN) and sky shot
noise ({sky_shot_dn:.0f} DN) — more than three orders of magnitude below any
real noise term. Inputs: {res['n_bias']} bias, {res['n_dark']} darks (60 s),
{res['n_flat']} flats, one 60 s NGC 6811 $B$ light (Moravian C3-61000, 2×2,
night {res['night']}); reference {res['reference']}. The cosmetic
(L.A.Cosmic + hot-pixel) stage is disabled here — it repairs ~1% of pixels by
design and is validated separately by injection; this figure isolates the
bias/dark/flat arithmetic. Generator: `calib_crosscheck_ngc6811.py` (input
file names recorded in the JSON).
""", encoding="utf-8")

    print("=== fig12 rework ===")
    for k, _ in STEPS:
        print(f"  {k:14s} max|Δ|={steps[k]['max_abs']:.2e}")
    for ext, p in paths.items():
        print(f"wrote {ext}: {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
