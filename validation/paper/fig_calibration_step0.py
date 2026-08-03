# -*- coding: utf-8 -*-
"""Figure — detector calibration (Step 0): what it does, and that it is exact.

Replaces the legacy 4-panel fig10_calibration dropped on 2026-08-03 (user:
"그림2는 좀 버려라"). That figure mixed a synthetic-injection claim into a
caption whose panels did not show it, and carried stale numbers. This one shows
the four things §3.2 actually asserts, each from its own committed artefact:

  (a) real M13 V frame before/after — vignette and fixed pattern removed
  (b) controlled truth recovery — inject known bias/dark/flat, invert, and the
      residual sits at the injected read-noise floor with no systematic offset
  (c) cosmetic stage — injected cosmic rays and hot pixels removed with the
      stars untouched
  (d) cross-implementation — the same raw reduced by the AstralImage/AIPPI
      engine is bit-identical across 19 datasets

Every number is recomputed here, not copied: (b) reruns the calibration on the
same fixture the regression test uses (seed 1234), (c) reruns
apex.benchmark.cosmetic_validate (seed fixed in that module), (d) reads
calibration_equivalence_multi.json.

Run: .venv-deploy\\Scripts\\python -X utf8 validation\\paper\\fig_calibration_step0.py
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))
sys.path.insert(0, str(REPO / "tests"))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

DATA = REPO / "validation" / "paper" / "data"
OUTDIR = REPO / "validation" / "paper" / "figures"
CAPDIR = REPO / "validation" / "paper" / "captions"


def stretch(img, lo=25.0, hi=99.5):
    v = img[np.isfinite(img)]
    a, b = np.percentile(v, lo), np.percentile(v, hi)
    x = np.clip((img - a) / max(b - a, 1e-6), 0, 1)
    return np.arcsinh(8 * x) / np.arcsinh(8.0)


def truth_recovery() -> dict:
    """(b) 알려진 bias/dark/flat 을 주입하고 0단계로 되돌린다.

    tests/test_calibration.py 의 fixture(시드 1234)를 그대로 쓴다 — 회귀 테스트가
    지키는 것과 같은 구성이라야 그림의 수치가 테스트와 어긋나지 않는다.
    """
    import test_calibration as T
    from apex.analysis import calibration as cal

    opts = cal.CalibrationOptions(combine_method="median", pedestal_mode="none",
                                  cosmetic_enable=False)
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        bias_p, dark_p, flat_p, light_p, (B, Dr, F, S) = T._make_dataset(tmp)
        mbias, _ = cal.build_master_bias(bias_p, opts)
        mdark, dexp, _ = cal.build_master_dark(dark_p, opts, master_bias=mbias)
        mflat, _ = cal.build_master_flat(flat_p, opts, master_bias=mbias,
                                         master_dark=mdark, dark_exp=dexp)
        calibrated, _h, _q = cal.calibrate_light_file(
            light_p, opts, master_bias=mbias, master_dark=mdark,
            dark_exp=dexp, master_flat=mflat)
        resid = calibrated.astype(np.float64) - S
        med = float(np.nanmedian(resid))
        mad = float(np.nanmedian(np.abs(resid - med)) * 1.4826)
        no_flat, _, _ = cal.calibrate_light_file(
            light_p, opts, master_bias=mbias, master_dark=mdark, dark_exp=dexp)
        cy = cx = 5
        return {"resid": resid[np.isfinite(resid)].ravel(),
                "med": med, "mad": mad,
                "rn": 3.0,                                   # 주입 읽기잡음 (T._write)
                "mbias_rms": float(np.sqrt(np.nanmean((mbias - B) ** 2))),
                "vign_before": abs(float(no_flat[cy, cx]) - S[cy, cx]),
                "vign_after": abs(float(calibrated[cy, cx]) - S[cy, cx])}


def main() -> int:
    tr = truth_recovery()
    from apex.benchmark import cosmetic_validate
    cm = cosmetic_validate.run()
    equiv = [r for r in json.loads(
        (REPO / "validation" / "paper" / "calibration_equivalence_multi.json")
        .read_text(encoding="utf-8")) if not r.get("error")]

    raw = np.load(DATA / "calib_m13_raw.npy").astype(float)
    cal_img = np.load(DATA / "calib_m13_cal.npy").astype(float)

    # 1x4 는 DOUBLE_COL 에서 패널당 1.75in 뿐이라 제목·주석이 서로 겹친다
    # (2026-08-03 렌더 검수). 2x2 로 짜면 패널당 3.4in 라 여유가 생긴다.
    # 높이는 조판 예산에 걸린다. 판면 전폭 690px 에서 세로가 454px 를 넘으면
    # 지면의 그림 자리 하나를 통째로 먹어 다음 그림이 다음 지면으로 밀린다
    # (2026-08-03: 4.9in 판에서 참조–그림 거리가 2.2 -> 4.3쪽). 4.1in = 405px.
    fig = plt.figure(figsize=(DOUBLE_COL, 4.1))
    gs = fig.add_gridspec(2, 2, wspace=0.26, hspace=0.52,
                          left=0.078, right=0.985, top=0.905, bottom=0.155)

    # ── (a) 실측 프레임 보정 전/후 + 반경 프로파일 ──
    # 컷아웃만 나란히 두면 비네팅 제거가 눈에 안 보인다(2026-08-03 렌더 검수).
    # flat 자체의 반경 프로파일을 함께 보여 무엇이 나눠졌는지 드러낸다.
    gsa = gs[0, 0].subgridspec(2, 1, height_ratios=[2.1, 1.0], hspace=0.30)
    gsi = gsa[0].subgridspec(1, 2, wspace=0.04)
    for k, (img, lab) in enumerate(((raw, "raw"), (cal_img, "calibrated"))):
        ax = fig.add_subplot(gsi[k])
        ax.imshow(stretch(img), cmap="gray", origin="lower", aspect="equal",
                  interpolation="nearest")
        ax.text(0.05, 0.955, lab, transform=ax.transAxes, fontsize=7,
                color="white", va="top",
                bbox=dict(fc="black", ec="none", alpha=0.5, pad=1.3))
        ax.set_xticks([]); ax.set_yticks([])
        if k == 0:
            ax.set_title("(a) vignette + fixed pattern removed",
                         loc="left", fontsize=8.2)
    prof = np.load(DATA / "calib_m13_profile.npy").astype(float)
    axp = fig.add_subplot(gsa[1])
    x = np.arange(prof.shape[1]) / (prof.shape[1] - 1)
    axp.plot(x, prof[0], color=C["reference"], lw=1.3, label="master flat")
    axp.axhline(1.0, color=PALETTE["grey"], lw=0.8, ls=":")
    lo, hi = float(np.min(prof[0])), float(np.max(prof[0]))
    axp.set_ylim(lo - 0.02, hi + 0.03)
    axp.set_xlim(0, 1)
    axp.set_xticks([0, 0.5, 1.0]); axp.set_xticklabels(["edge", "centre", "edge"],
                                                       fontsize=6.8)
    axp.set_ylabel("flat", fontsize=7.4)
    axp.tick_params(labelsize=6.6)
    axp.text(0.5, 0.06, f"divides out a {(1 - lo) * 100:.0f}% edge-to-centre gradient",
             transform=axp.transAxes, fontsize=6.8, ha="center", va="bottom",
             color=PALETTE["grey"])

    # ── (b) 참값 복원 잔차 ──
    axb = fig.add_subplot(gs[0, 1])
    r = tr["resid"]
    lim = 4.0 * tr["rn"]
    axb.hist(r, bins=np.linspace(-lim, lim, 90), color=C["data"],
             density=True, zorder=2)
    xs = np.linspace(-lim, lim, 400)
    axb.plot(xs, np.exp(-0.5 * (xs / tr["rn"]) ** 2) / (tr["rn"] * np.sqrt(2 * np.pi)),
             color=C["reference"], lw=1.5, zorder=3,
             label=f"injected floor $\\sigma={tr['rn']:.1f}$ DN")
    axb.axvline(0.0, color=PALETTE["grey"], lw=0.8, ls=":", zorder=1)
    axb.set_xlim(-lim, lim)
    axb.set_ylim(0, 0.185)                       # 주석 자리를 남긴다
    axb.set_xlabel("recovered $-$ true  (DN)", fontsize=7.6)
    axb.set_ylabel("density", fontsize=7.6)
    axb.tick_params(labelsize=6.8)
    axb.legend(fontsize=6.6, loc="upper right", handlelength=1.5,
               framealpha=0.9, borderpad=0.3)
    axb.text(0.035, 0.94,
             f"offset {tr['med']:+.3f} DN\nMAD {tr['mad']:.2f} DN",
             transform=axb.transAxes, fontsize=7.2, va="top", color=C["data"])
    axb.set_title("(b) truth recovery: offset-free, at the noise floor",
                  loc="left", fontsize=8.2)

    # ── (c) cosmetic ──
    axc = fig.add_subplot(gs[1, 0])
    names = ["cosmic rays\nremoved", "hot pixels\nremoved", "star cores\ntouched"]
    vals = [cm["cr_completeness"] * 100, cm["hot_completeness"] * 100,
            cm["star_core_false_positive"] * 100]
    cols = [C["data"], C["data"], C["bad"]]
    xb = np.arange(3)
    axc.bar(xb, vals, color=cols, width=0.58, zorder=3)
    axc.axhline(100, color=PALETTE["grey"], lw=0.8, ls=":", zorder=1)
    axc.set_ylim(0, 132)                          # 주석 자리
    axc.set_yticks([0, 25, 50, 75, 100])
    axc.set_xticks(xb); axc.set_xticklabels(names, fontsize=7.0)
    axc.set_ylabel("injected pixels (%)", fontsize=7.6)
    axc.tick_params(labelsize=6.8)
    for xi, v in zip(xb, vals):
        axc.annotate(f"{v:.0f}%" if v >= 1 else f"{v:.3f}%", (xi, v),
                     xytext=(0, 3), textcoords="offset points",
                     ha="center", va="bottom", fontsize=7.4,
                     color=C["bad"] if xi == 2 else "black")
    axc.text(0.5, 0.97,
             f"aperture flux change {cm['flux_preservation_med_mmag']:.2f} mmag (median)",
             transform=axc.transAxes, fontsize=7.0, va="top", ha="center",
             color=PALETTE["grey"])
    axc.set_title("(c) cosmetics removed, stars untouched",
                  loc="left", fontsize=8.2)

    # ── (d) AIPPI 교차구현 ──
    axd = fig.add_subplot(gs[1, 1])
    keys = ["calibrated_rms", "master_bias_rms", "master_dark_rms", "master_flat_rms"]
    labs = ["calibrated\nframe", "master\nbias", "master\ndark", "master\nflat"]
    FLOOR = 1e-4
    worst = [max(float(np.max([r[k] for r in equiv])), FLOOR) for k in keys]
    exact = [all(r[k] == 0.0 for r in equiv) for k in keys]
    xb = np.arange(4)
    axd.bar(xb, worst, width=0.58, zorder=3,
            color=[C["model"] if e else C["data"] for e in exact])
    axd.set_yscale("log")
    axd.set_ylim(FLOOR * 0.5, 12.0)               # 주석 자리
    axd.set_xticks(xb); axd.set_xticklabels(labs, fontsize=7.0)
    axd.set_ylabel("worst difference RMS (DN)", fontsize=7.6)
    axd.tick_params(labelsize=6.8)
    for xi, w, e in zip(xb, worst, exact):
        axd.annotate("= 0" if e else f"{w:.2g}", (xi, w),
                     xytext=(0, 3), textcoords="offset points",
                     ha="center", va="bottom", fontsize=7.4,
                     color=C["model"] if e else "black")
    axd.text(0.5, 0.97, f"{len(equiv)} datasets · APEX vs AstralImage/AIPPI",
             transform=axd.transAxes, fontsize=7.0, va="top", ha="center",
             color=PALETTE["grey"])
    axd.set_title("(d) independent engine: calibrated frames bit-identical",
                  loc="left", fontsize=8.2)

    prov1 = ("(a) M13 $V$ 60 s, Moravian C3-61000 (2$\\times$2) · "
             "(b) synthetic bias/dark/flat injected and inverted through APEX Step 0, "
             "seed 1234 (same fixture as the regression test)")
    prov2 = (f"(c) apex.benchmark.cosmetic_validate: {cm['n_stars']} stars, "
             f"{cm['n_cr']} CR pixels, {cm['n_hot']} hot pixels injected · "
             f"(d) {len(equiv)} cluster/variable datasets, 9 bands · "
             "gen: fig_calibration_step0.py")
    fig.text(0.985, 0.042, prov1, ha="right", va="bottom", fontsize=5.9,
             color=PALETTE["grey"])
    fig.text(0.985, 0.008, prov2, ha="right", va="bottom", fontsize=5.9,
             color=PALETTE["grey"])

    paths = save_fig(fig, "fig_calibration_step0", OUTDIR)
    plt.close(fig)

    CAPDIR.mkdir(parents=True, exist_ok=True)
    (CAPDIR / "fig_calibration_step0.md").write_text(
        f"""# Figure — detector calibration (Step 0)

**(a)** A real M13 $V$ frame before and after APEX Step 0: the optical vignette
and the fixed-pattern structure are removed. **(b)** Controlled truth recovery —
known bias, dark and flat are injected into a synthetic frame and inverted
through Step 0; the residual against the true science frame has systematic
offset {tr['med']:+.3f} DN and scatter MAD {tr['mad']:.2f} DN, the latter equal
to the injected read-noise floor ({tr['rn']:.1f} DN). The master bias is
recovered to RMS {tr['mbias_rms']:.2f} DN and the unit-median flat reduces the
corner vignette residual from {tr['vign_before']:.1f} to {tr['vign_after']:.1f}
DN. **(c)** The optional cosmetic stage (L.A.Cosmic via astroscrappy) removes
{cm['cr_completeness']*100:.0f} per cent of injected cosmic-ray pixels and
{cm['hot_completeness']*100:.0f} per cent of hot pixels while touching
{cm['star_core_false_positive']*100:.3f} per cent of star-core pixels; aperture
fluxes shift by {cm['flux_preservation_med_mmag']:.2f} mmag (median).
**(d)** The same raw frames reduced by the independent AstralImage/AIPPI engine:
across {len(equiv)} datasets in 9 bands the calibrated frames are bit-identical
(difference RMS and maximum both 0 DN), and the master frames agree to
$\\leq${max(float(np.max([r['master_bias_rms'] for r in equiv])), 0):.2f} DN RMS.
Generator: `fig_calibration_step0.py`; every number is recomputed at figure time.
""", encoding="utf-8")

    print("=== fig_calibration_step0 ===")
    print(f"  (b) offset {tr['med']:+.4f} DN  MAD {tr['mad']:.3f} DN  "
          f"bias RMS {tr['mbias_rms']:.3f}  vign {tr['vign_before']:.1f}->{tr['vign_after']:.1f}")
    print(f"  (c) CR {cm['cr_completeness']*100:.1f}%  hot {cm['hot_completeness']*100:.1f}%  "
          f"starFP {cm['star_core_false_positive']*100:.4f}%  "
          f"dmag {cm['flux_preservation_med_mmag']:.3f} mmag")
    print(f"  (d) {len(equiv)} datasets; worst per key: " +
          ", ".join(f"{k}={max(float(np.max([r[k] for r in equiv])),0):.3g}" for k in keys))
    for ext, p in paths.items():
        print(f"wrote {ext}: {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
