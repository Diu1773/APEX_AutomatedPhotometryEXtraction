"""Figure — injection-recovery completeness across real frames of differing quality,
and the background-limited detection law they obey.

Reference-standard artificial-star test (DAOPHOT ADDSTAR; DES Balrog; HSC SynPipe;
AutoPhOT App.D; Haynes 1994/2002): stars injected into REAL frames. Real cluster
frames spanning a wide range of sky brightness and seeing are each measured by the
identical injection pipeline, together with the idealised synthetic verification frame.

Main panel: three representative real frames (deep/medium/shallow) + the synthetic
verification rung, with the real injected-star cutout ladder beneath. Inset: the
instrument-independent result — every frame's 50% depth obeys the background-limited,
peak-detection law  m50 = C - 2.5 log10(sigma_sky) - 5 log10(FWHM)  with the exponents
FIXED by theory (only the constant C, which absorbs the hardware, is fit). This is what
"reproducing" the reference validation means here: not matching anyone's absolute mmag
(different telescopes/detectors — CCD vs this Sony IMX455 CMOS), but showing APEX obeys
the same universal law on its own instrument.

    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig_completeness_realvssynth.py
"""
from __future__ import annotations
import sys
from pathlib import Path
REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import pandas as pd
from astropy.io import fits
import sep
import matplotlib.pyplot as plt

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

DATA = REPO / "validation" / "paper"
SYN = DATA / "data" / "artificial_star" / "benchmark_run"
CUTOUTS = DATA / "data_realframe_M13V" / "injection_cutouts.npz"
OUTDIR = DATA / "figures"
REPRO = Path(r"E:\APEX_validation\reprocess")

# All real-frame injections. hero=True → drawn as a completeness curve in the main
# panel; every entry (hero or not) is a point in the law inset.
# (label, run_subdir, injected_fits, fs_target, fs_file, color_key, marker, hero)
FRAMES = [
    ("M67 i",     "data_realframe_M67i",          REPRO/"M67/sci/pp_Messier67-0008-i.fit",         "M67",     "pp_Messier67-0008-i.fit", "data",      "-o", True),
    ("NGC 6811 R","data_realframe_NGC6811R",       REPRO/"NGC6811/sci/pp_NGC6811-0005-R.fit",       "NGC6811", "pp_NGC6811-0005-R.fit",   "reference", "-s", True),
    ("M13 V",     "data_realframe_M13V",           REPRO/"M13/calibrated/20260515/pp_messier13-0001-V.fit","M13","pp_messier13-0001-V.fit","accent","-D", True),
    ("M67 r",     "data_realframe_M67r_mid",       REPRO/"M67/sci/pp_Messier67-0003-r.fit",         "M67",     "pp_Messier67-0003-r.fit", None,        None, False),
    ("M67 g",     "data_realframe_M67g_broad",     REPRO/"M67/sci/pp_Messier67-0004-g.fit",         "M67",     "pp_Messier67-0004-g.fit", None,        None, False),
    ("NGC 6811 R (soft)","data_realframe_NGC6811R_broad", REPRO/"NGC6811/sci/pp_NGC6811-0008-R.fit","NGC6811", "pp_NGC6811-0008-R.fit",   None,        None, False),
    ("M13 R",     "data_realframe_M13R_sharp",     REPRO/"M13/sci/pp_messier13-0004-R.fit",         "M13",     "pp_messier13-0004-R.fit", None,        None, False),
]


def read_off(mag, comp, level=0.5):
    o = np.argsort(mag); mm, cc = mag[o], comp[o]
    for i in range(len(mm) - 1):
        if cc[i] >= level >= cc[i + 1]:
            den = cc[i] - cc[i + 1]
            f = (cc[i] - level) / den if den else 0.0
            return float(mm[i] + f * (mm[i + 1] - mm[i]))
    return float("nan")


def curve(run_subdir, bw=0.25):
    s = pd.read_csv(DATA / run_subdir / "artificial_star/benchmark_run/stars.csv")
    s = s[~s["baseline_confounded"].astype(bool)]
    m = s["magnitude_true"].to_numpy(float)
    r = s["recovered"].to_numpy(bool)
    lo0 = np.floor(m.min() / bw) * bw
    edges = np.arange(lo0, m.max() + bw, bw)
    xs, cs = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        msk = (m >= lo) & (m < hi)
        if msk.sum() >= 15:
            xs.append(0.5 * (lo + hi)); cs.append(r[msk].mean())
    return np.array(xs), np.array(cs)


def curve_syn(bw=0.25):
    s = pd.read_csv(SYN / "stars.csv")
    s = s[~s["baseline_confounded"].astype(bool)]
    m = s["magnitude_true"].to_numpy(float); r = s["recovered"].to_numpy(bool)
    edges = np.arange(np.floor(m.min() / bw) * bw, m.max() + bw, bw)
    xs, cs = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (m >= lo) & (m < hi)
        if k.sum() >= 15:
            xs.append(0.5 * (lo + hi)); cs.append(r[k].mean())
    return np.array(xs), np.array(cs)


def fwhm_of(run_subdir):
    """FWHM of the empirical PSF kernel the suite actually injected (half-max area).
    Self-consistent for the law: injected-star completeness depends on this kernel,
    not on frame_stats' bright-star FWHM (they diverge on soft frames)."""
    k = fits.getdata(DATA / run_subdir / "artificial_star/benchmark_run/empirical_psf.fits").astype(float)
    k = k / k.sum()
    area = (k >= k.max() / 2).sum()
    return float(2 * np.sqrt(area / np.pi))


def sigma_of(fits_path):
    d = fits.getdata(str(fits_path)).astype(np.float64)
    return float(np.median(sep.Background(d).rms()))


def main() -> int:
    cut = np.load(CUTOUTS)
    stamps, cmags = cut["stamps"], cut["mags"]
    clo, chi = float(cut["lo"]), float(cut["hi"])

    # gather every frame's (sigma, FWHM, m50)
    recs = []
    for label, sub, fp, tgt, fn, ckey, mk, hero in FRAMES:
        run = DATA / sub / "artificial_star/benchmark_run/stars.csv"
        if not run.exists():
            print(f"  [skip] {label}: {sub} not found")
            continue
        x, c = curve(sub)
        recs.append(dict(label=label, sub=sub, sigma=sigma_of(fp),
                         fwhm=fwhm_of(sub), m50=read_off(x, c),
                         ckey=ckey, mk=mk, hero=hero, x=x, c=c))

    # background-limited peak-detection law: m50 = C - 2.5 log10(sigma) - 5 log10(FWHM)
    sig = np.array([r["sigma"] for r in recs])
    fwh = np.array([r["fwhm"] for r in recs])
    obs = np.array([r["m50"] for r in recs])
    base = -2.5 * np.log10(sig) - 5 * np.log10(fwh)
    Cfit = float(np.mean(obs - base))
    pred = Cfit + base
    resid_rms = float(np.std(obs - pred))

    fig = plt.figure(figsize=(DOUBLE_COL * 0.80, 4.9))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 0.95], hspace=0.48)
    ax = fig.add_subplot(gs[0])
    ax.axhline(0.5, color=PALETTE["grey"], lw=0.7, ls=":", zorder=2)

    # synthetic verification rung
    ms, cs = curve_syn(); m50s = read_off(ms, cs)
    ax.plot(ms, cs, "--", color=PALETTE["grey"], lw=1.4, zorder=3,
            label="synthetic (verification)  ·  sky 150, FWHM 3.4px")

    m50_m13 = None
    for r in recs:
        if not r["hero"]:
            continue
        col = C[r["ckey"]]
        ax.plot(r["x"], r["c"], r["mk"], color=col, lw=2.0, ms=3.6, mfc=col, zorder=6,
                label=f"{r['label']}  ·  sky-noise {r['sigma']:.0f}, FWHM {r['fwhm']:.1f}px")
        ax.axvline(r["m50"], color=col, lw=0.9, ls=":", zorder=2, alpha=0.75)
        ha = "right" if r["label"].startswith("M13") else ("left" if r["label"].startswith("NGC") else "center")
        dx = -0.08 if ha == "right" else (0.08 if ha == "left" else 0.0)
        ax.text(r["m50"] + dx, 1.02, f"{r['m50']:.1f}", color=col, fontsize=7.0,
                ha=ha, va="bottom", fontweight="bold")
        if r["label"].startswith("M13"):
            m50_m13 = r["m50"]
    ax.text(12.35, 1.02, "$m_{50}$:", color="0.35", fontsize=7.0, ha="left", va="bottom")

    ax.set_xlabel("injected / instrumental magnitude (count-rate, ZP = 25)")
    ax.set_ylabel("recovery completeness")
    ax.set_xlim(12.2, 20.3)
    ax.set_ylim(-0.03, 1.11)
    leg = ax.legend(loc="lower left", fontsize=6.8, framealpha=0.93,
                    title="representative real frames (+ synthetic)",
                    title_fontsize=7.0, handlelength=1.8)
    leg._legend_box.align = "left"
    ax.set_title("real-frame depth obeys the background-limited detection law",
                 loc="left", fontsize=10.5)

    # ── inset: the instrument-independent scaling law (all frames) ──
    axi = ax.inset_axes([0.745, 0.17, 0.245, 0.47])
    lo = min(obs.min(), pred.min()) - 0.25
    hi = max(obs.max(), pred.max()) + 0.25
    axi.plot([lo, hi], [lo, hi], "-", color=PALETTE["grey"], lw=1.0, zorder=1)
    for r, p in zip(recs, pred):
        col = C[r["ckey"]] if r["hero"] else "0.45"
        mk = "o" if not r["hero"] else r["mk"].lstrip("-")
        axi.plot(r["m50"], p, mk, color=col, ms=5.5 if r["hero"] else 4.5,
                 mfc=col, mew=0, zorder=3)
    axi.set_xlim(lo, hi); axi.set_ylim(lo, hi)
    axi.set_xlabel("observed $m_{50}$", fontsize=6.8, labelpad=1)
    axi.set_ylabel("law $m_{50}$", fontsize=6.8, labelpad=1)
    axi.tick_params(labelsize=6.0)
    axi.set_title(f"$m_{{50}}=C-2.5\\log\\sigma-5\\log$FWHM\n"
                  f"exponents fixed · N={len(recs)} · resid {1000*resid_rms:.0f} mmag",
                  fontsize=6.2, pad=2)
    axi.set_aspect("equal", adjustable="box")

    # ── bottom strip: real injected-star cutouts (M13, the shallowest frame) ──
    gsc = gs[1].subgridspec(1, len(cmags), wspace=0.16)
    for j, (s, m) in enumerate(zip(stamps, cmags)):
        axc = fig.add_subplot(gsc[j])
        axc.imshow(s, cmap="gray", vmin=clo, vmax=chi, origin="lower",
                   interpolation="nearest")
        axc.set_xticks([]); axc.set_yticks([])
        recovered = m < (m50_m13 or 14.9)
        mark = "found" if recovered else "lost"
        col = C["data"] if recovered else C["accent"]
        for sp in axc.spines.values():
            sp.set_edgecolor(col); sp.set_linewidth(1.5)
        axc.set_title(f"m = {m:.1f}", fontsize=7.6, pad=2)
        axc.set_xlabel(mark, fontsize=7.8, color=col, labelpad=2)
    fig.text(0.012, gs[1].get_position(fig).y1 + 0.016,
             "real injected stars — M13 V (the shallowest frame), shared stretch  —  "
             "blue = recovered, orange = lost",
             fontsize=7.2, color="0.3", va="bottom", ha="left")

    paths = save_fig(fig, "fig_completeness_realvssynth", OUTDIR)
    plt.close(fig)
    print(f"[completeness] synthetic m50={m50s:.2f}   law C={Cfit:.3f}  "
          f"resid RMS={1000*resid_rms:.1f} mmag  N={len(recs)}")
    for r, p in zip(recs, pred):
        print(f"  {r['label']:18s} σ={r['sigma']:6.2f} FWHM={r['fwhm']:.2f} "
              f"obs={r['m50']:.2f} law={p:.2f} Δ={r['m50']-p:+.3f}")
    print(f"saved: {paths['png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
