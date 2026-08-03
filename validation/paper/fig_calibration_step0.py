# -*- coding: utf-8 -*-
"""Figure — detector calibration (Step 0): the frames applied, and their effect.

2026-08-03 second design. The first attempt put three abstract bar charts next
to one before/after pair; the user's objections were all correct:
  * the before/after cutouts looked identical, so the panel proved nothing
  * "recovered - true" was an opaque axis label
  * a "star cores touched 0.000%" bar sat on the same axis as two ~100% bars
    and read as noise
  * what a reader actually wants is the calibration frames themselves, and the
    light frame before and after each operation

So this figure shows exactly that: the three masters as images (top row), the
same science frame after each operation (bottom row), and the sky profile that
makes the vignette removal visible. The quantitative validation lives in the
text and in Figs 3-5 (detector constants, ccdproc, cross-instrument).

Data: data/calib_stages.npz, written by calib_crosscheck_ngc6811.py from the
real 2026-06-11 night (8 bias, 8x60 s dark, 5 B flats, one 60 s NGC 6811 B).

Run: .venv-deploy\\Scripts\\python -X utf8 validation\\paper\\fig_calibration_step0.py
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

DATA = REPO / "validation" / "paper" / "data"
OUTDIR = REPO / "validation" / "paper" / "figures"
CAPDIR = REPO / "validation" / "paper" / "captions"


def show(ax, img, title, lo=1.0, hi=99.0, vmin=None, vmax=None, cmap="gray"):
    v = img[np.isfinite(img)]
    if vmin is None:
        vmin, vmax = np.percentile(v, lo), np.percentile(v, hi)
    ax.imshow(img, cmap=cmap, origin="lower", aspect="equal",
              vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, loc="left", fontsize=7.6, pad=2.5)
    return vmin, vmax


def note(ax, text, color=None):
    """이미지 패널 아래 한 줄. 눈금이 없는 패널 전용."""
    ax.text(0.5, -0.045, text, transform=ax.transAxes, fontsize=6.6,
            ha="center", va="top", color=color or PALETTE["grey"])


def inset_note(ax, text, x=0.5, y=0.055, ha="center", color=None):
    """눈금이 있는 선 그래프는 축 안에 적는다 — 축 밖에 두면 눈금과 겹친다."""
    ax.text(x, y, text, transform=ax.transAxes, fontsize=6.6,
            ha=ha, va="bottom", color=color or PALETTE["grey"],
            bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.0))


def main() -> int:
    d = np.load(DATA / "calib_stages.npz")
    res = json.loads((DATA / "calib_crosscheck_ngc6811.json").read_text(encoding="utf-8"))

    mb, md, mf = d["master_bias"], d["master_dark"], d["master_flat"]
    raw, bd, cal = d["light_raw"], d["light_bd"], d["light_cal"]
    p_bd, p_cal, p_flat = d["prof_bd"], d["prof_cal"], d["prof_flat"]

    fig = plt.figure(figsize=(DOUBLE_COL, 4.35))
    gs = fig.add_gridspec(2, 4, wspace=0.16, hspace=0.40,
                          left=0.035, right=0.985, top=0.925, bottom=0.145)

    # ── 윗줄: 적용되는 세 마스터 ──
    ax = fig.add_subplot(gs[0, 0])
    show(ax, mb, "(a) master bias")
    note(ax, f"median {float(d['bias_med']):.0f} DN · 8 frames")

    ax = fig.add_subplot(gs[0, 1])
    show(ax, md, "(b) master dark", lo=1.0, hi=99.5)
    note(ax, f"median {float(d['dark_med']):.1f} DN in {float(d['dark_exp']):.0f} s · 8 frames")

    ax = fig.add_subplot(gs[0, 2])
    show(ax, mf, "(c) master flat", lo=0.5, hi=99.5)
    note(ax, f"median {float(d['flat_med']):.3f} · 5 frames")

    # (d) flat 의 가로 프로파일 — 무엇이 나눠지는지
    axd = fig.add_subplot(gs[0, 3])
    x = np.linspace(0, 1, p_flat.size)
    axd.plot(x, p_flat, color=C["reference"], lw=1.4)
    axd.axhline(1.0, color=PALETTE["grey"], lw=0.8, ls=":")
    lo, hi = float(np.min(p_flat)), float(np.max(p_flat))
    axd.set_ylim(lo - 0.015, hi + 0.02)
    axd.set_xlim(0, 1)
    axd.set_xticks([0, 0.5, 1.0]); axd.set_xticklabels(["left", "centre", "right"],
                                                       fontsize=6.6)
    axd.set_ylabel("flat", fontsize=7.0, labelpad=1.5)
    axd.tick_params(labelsize=6.4)
    axd.set_title("(d) what the flat divides out", loc="left", fontsize=7.6, pad=2.5)
    inset_note(axd, f"{(hi - lo) * 100:.0f}% across the field")

    # ── 아랫줄: 같은 과학 프레임, 연산마다 ──
    ax = fig.add_subplot(gs[1, 0])
    show(ax, raw, "(e) raw science frame")
    note(ax, f"NGC 6811 $B$, {float(d['light_exp']):.0f} s")

    # (f)(g) 는 같은 눈금이라야 flat 효과가 보인다
    v = bd[np.isfinite(bd)]
    vmin, vmax = np.percentile(v, 1.0), np.percentile(v, 99.0)
    ax = fig.add_subplot(gs[1, 1])
    show(ax, bd, "(f) $-$ bias $-$ dark", vmin=vmin, vmax=vmax)
    note(ax, "offset and dark current gone")
    ax = fig.add_subplot(gs[1, 2])
    show(ax, cal, "(g) $\\div$ flat  = calibrated", vmin=vmin, vmax=vmax)
    note(ax, "same greyscale as (f)")

    # (h) 하늘 프로파일 — 비네팅이 실제로 평평해진다
    axh = fig.add_subplot(gs[1, 3])
    xb = np.linspace(0, 1, p_bd.size)
    n_bd = p_bd / float(np.median(p_bd))
    n_cal = p_cal / float(np.median(p_cal))
    axh.plot(xb, n_bd, color=C["bad"], lw=1.3, label="before flat (f)")
    axh.plot(xb, n_cal, color=C["data"], lw=1.3, label="after flat (g)")
    axh.axhline(1.0, color=PALETTE["grey"], lw=0.8, ls=":")
    axh.set_xlim(0, 1)
    axh.set_xticks([0, 0.5, 1.0]); axh.set_xticklabels(["left", "centre", "right"],
                                                       fontsize=6.6)
    axh.set_ylabel("sky / median", fontsize=7.0, labelpad=1.5)
    axh.tick_params(labelsize=6.4)
    axh.legend(fontsize=6.2, loc="lower left", handlelength=1.4,
               framealpha=0.9, borderpad=0.25)
    axh.set_title("(h) sky flattens", loc="left", fontsize=7.6, pad=2.5)
    span_before = (float(np.max(n_bd)) - float(np.min(n_bd))) * 100
    span_after = (float(np.max(n_cal)) - float(np.min(n_cal))) * 100
    # 곡선이 가운데서 솟으므로 왼쪽 위가 유일하게 빈 자리다
    inset_note(axh, f"peak-to-peak\n{span_before:.1f}% $\\to$ {span_after:.1f}%",
               x=0.035, y=0.80, ha="left")

    # 한 줄로 쓰면 그림 폭을 넘어 tight-bbox 가 캔버스를 왼쪽으로 늘린다
    # (2026-08-03 실측: 그림 왼쪽에 빈 여백 2인치). 두 줄로 나눈다.
    prov1 = (f"All frames real: Moravian C3-61000 (2$\\times$2), night {res['night']} · "
             f"{res['n_bias']} bias, {res['n_dark']} dark ({float(d['dark_exp']):.0f} s), "
             f"{res['n_flat']} flat, one {float(d['light_exp']):.0f} s NGC 6811 $B$ light")
    prov2 = (f"displayed at 1/{int(d['shrink'])} scale (block-mean) · "
             "gen: calib_crosscheck_ngc6811.py $\\to$ fig_calibration_step0.py")
    fig.text(0.985, 0.048, prov1, ha="right", va="bottom", fontsize=5.9,
             color=PALETTE["grey"])
    fig.text(0.985, 0.010, prov2, ha="right", va="bottom", fontsize=5.9,
             color=PALETTE["grey"])

    paths = save_fig(fig, "fig_calibration_step0", OUTDIR)
    plt.close(fig)

    CAPDIR.mkdir(parents=True, exist_ok=True)
    (CAPDIR / "fig_calibration_step0.md").write_text(
        f"""# Figure — detector calibration (Step 0)

**(a)-(c)** The three master frames APEX builds from the night's own
calibration exposures: bias (median {float(d['bias_med']):.0f} DN),
dark ({float(d['dark_med']):.1f} DN in {float(d['dark_exp']):.0f} s) and flat
(unit median). **(d)** The flat's horizontal profile — a
{(float(np.max(p_flat)) - float(np.min(p_flat))) * 100:.0f} per cent
response gradient across the field, which is what the division removes.
**(e)-(g)** The same NGC 6811 $B$ science frame after each operation: raw,
then bias- and dark-subtracted, then flat-fielded. Panels (f) and (g) share a
greyscale so the change is the flat-fielding alone. **(h)** The sky profile
before and after flat-fielding, each normalised to its own median: the
peak-to-peak gradient falls from {span_before:.1f} to {span_after:.1f} per
cent. All frames are real (Moravian C3-61000, 2x2, night {res['night']};
{res['n_bias']} bias, {res['n_dark']} darks, {res['n_flat']} flats), shown at
1/{int(d['shrink'])} scale. Numerical validation of these operations is
separate: recovery of injected truth (text), the detector constants (Fig. 3),
pixel-for-pixel agreement with ccdproc (Fig. 4) and reproduction on two more
cameras (Fig. 5).
""", encoding="utf-8")

    print("=== fig_calibration_step0 (v2) ===")
    print(f"  flat gradient {(float(np.max(p_flat)) - float(np.min(p_flat))) * 100:.1f}%")
    print(f"  sky peak-to-peak {span_before:.2f}% -> {span_after:.2f}%")
    print(f"  bias {float(d['bias_med']):.1f} DN · dark {float(d['dark_med']):.2f} DN"
          f" · flat {float(d['flat_med']):.4f}")
    for ext, p in paths.items():
        print(f"wrote {ext}: {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
