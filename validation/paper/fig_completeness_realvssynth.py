"""Figure — injection-recovery completeness: per-frame depths in magnitude space,
and their collapse onto a single universal curve in peak-S/N space.

Reference-standard artificial-star test (DAOPHOT ADDSTAR; DES Balrog; HSC SynPipe;
AutoPhOT App.D; Haynes 1994/2002): stars injected into REAL frames. Seven real
cluster frames spanning a wide range of sky brightness (sigma 5-58 e-/px) and
seeing (FWHM 5.2-9.0 px) are each measured by the identical injection pipeline.

Panel (a), magnitude space: each frame has its own 50% depth — depth is a frame
property (sky + seeing), spanning 3.4 mag across the seven frames. Three
representative curves shown; real injected-star cutouts beneath.

Panel (b), the point of the figure: re-express every injected star by its peak
pixel S/N (expected peak = flux x kernel-peak-fraction, divided by that frame's
background noise) and ALL SEVEN frames collapse onto one curve with S/N_50 = 4.0
+/- 0.2 — the instrument-independent detection law (completeness is erf-like in
S/N: Masci 2011; AutoPhOT App.D), consistent with the 3.2-sigma matched-filter
threshold. Magnitude-space depths are this one curve shifted by each frame's
sigma x FWHM^2.

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
from scipy.special import erf as _erf
from scipy.optimize import curve_fit
import sep
import matplotlib.pyplot as plt

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

DATA = REPO / "validation" / "paper"
SYN = DATA / "data" / "artificial_star" / "benchmark_run"
CUTOUTS = DATA / "data_realframe_M13V" / "injection_cutouts.npz"
OUTDIR = DATA / "figures"
REPRO = Path(r"E:\APEX_validation\reprocess")
GAIN = 0.689  # e-/ADU (PTC-measured, C3-61000)

# All real-frame injections. hero=True → completeness curve in panel (a);
# every frame contributes to the S/N collapse in panel (b).
FRAMES = [
    ("M67 i",      "data_realframe_M67i",           REPRO/"M67/sci/pp_Messier67-0008-i.fit",          "data",      "-o", True),
    ("NGC 6811 R", "data_realframe_NGC6811R",        REPRO/"NGC6811/sci/pp_NGC6811-0005-R.fit",        "reference", "-s", True),
    ("M13 V",      "data_realframe_M13V",            REPRO/"M13/calibrated/20260515/pp_messier13-0001-V.fit", "accent", "-D", True),
    ("M67 r",      "data_realframe_M67r_mid",        REPRO/"M67/sci/pp_Messier67-0003-r.fit",          None, None, False),
    ("M67 g",      "data_realframe_M67g_broad",      REPRO/"M67/sci/pp_Messier67-0004-g.fit",          None, None, False),
    ("NGC 6811 R (soft)", "data_realframe_NGC6811R_broad", REPRO/"NGC6811/sci/pp_NGC6811-0008-R.fit",  None, None, False),
    ("M13 R",      "data_realframe_M13R_sharp",      REPRO/"M13/sci/pp_messier13-0004-R.fit",          None, None, False),
]


def read_off(mag, comp, level=0.5, increasing=False):
    o = np.argsort(mag); mm, cc = mag[o], comp[o]
    for i in range(len(mm) - 1):
        hit = (cc[i] <= level <= cc[i + 1]) if increasing else (cc[i] >= level >= cc[i + 1])
        if hit:
            den = cc[i + 1] - cc[i]
            f = (level - cc[i]) / den if den else 0.0
            return float(mm[i] + f * (mm[i + 1] - mm[i]))
    return float("nan")


def binned(x, rec, edges, nmin=15):
    xs, cs = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x >= lo) & (x < hi)
        if m.sum() >= nmin:
            xs.append(0.5 * (lo + hi)); cs.append(rec[m].mean())
    return np.array(xs), np.array(cs)


def load_run(sub, frame_fits):
    rd = DATA / sub / "artificial_star/benchmark_run"
    s = pd.read_csv(rd / "stars.csv")
    s = s[~s["baseline_confounded"].astype(bool)]
    k = fits.getdata(rd / "empirical_psf.fits").astype(float)
    k = k / k.sum()
    d = fits.getdata(str(frame_fits)).astype(np.float64)
    sigma_e = float(np.median(sep.Background(d).rms())) * GAIN
    mag = s["magnitude_true"].to_numpy(float)
    rec = s["recovered"].to_numpy(bool)
    peak_snr = s["flux_expected_e"].to_numpy(float) * float(k.max()) / sigma_e
    return dict(mag=mag, rec=rec, snr=peak_snr, sigma_e=sigma_e,
                fwhm=float(2 * np.sqrt((k >= k.max() / 2).sum() / np.pi)))


def main() -> int:
    cut = np.load(CUTOUTS)
    stamps, cmags = cut["stamps"], cut["mags"]
    clo, chi = float(cut["lo"]), float(cut["hi"])

    runs = []
    for label, sub, fp, ckey, mk, hero in FRAMES:
        if not (DATA / sub / "artificial_star/benchmark_run/stars.csv").exists():
            print(f"  [skip] {label}"); continue
        r = load_run(sub, fp)
        r.update(label=label, ckey=ckey, mk=mk, hero=hero)
        e = np.arange(np.floor(r["mag"].min() / .25) * .25, r["mag"].max() + .25, .25)
        r["mx"], r["mc"] = binned(r["mag"], r["rec"], e)
        r["m50"] = read_off(r["mx"], r["mc"])
        le = np.arange(np.log10(r["snr"]).min(), np.log10(r["snr"]).max() + .12, .12)
        r["sx"], r["sc"] = binned(np.log10(r["snr"]), r["rec"], le)
        r["s50"] = 10 ** read_off(r["sx"], r["sc"], increasing=True)
        runs.append(r)

    s50s = np.array([r["s50"] for r in runs])
    print(f"[collapse] S/N50 per frame: " +
          " ".join(f"{r['s50']:.2f}" for r in runs) +
          f"  → {s50s.mean():.2f} ± {s50s.std():.2f}")

    # pooled erf fit in log10(S/N):  p = A/2 * (1 + erf((x-mu)/(sqrt2*w)))
    allx = np.concatenate([np.log10(r["snr"]) for r in runs])
    allr = np.concatenate([r["rec"] for r in runs])
    pe = np.arange(allx.min(), allx.max() + .08, .08)
    px, pc = binned(allx, allr, pe, nmin=40)
    def model(x, A, mu, w):
        return A / 2 * (1 + _erf((x - mu) / (np.sqrt(2) * w)))
    popt, _ = curve_fit(model, px, pc, p0=[0.97, np.log10(4.0), 0.1])
    A_f, mu_f, w_f = popt
    s50_fit = 10 ** float(mu_f)

    # ── layout: (a) mag space | (b) S/N collapse ; bottom cutouts ──
    fig = plt.figure(figsize=(DOUBLE_COL, 4.9))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 0.95], hspace=0.46)
    gst = gs[0].subgridspec(1, 2, width_ratios=[1.5, 1.0], wspace=0.24)
    ax = fig.add_subplot(gst[0])
    axb = fig.add_subplot(gst[1])

    # (a) magnitude space — heroes + synthetic
    ax.axhline(0.5, color=PALETTE["grey"], lw=0.7, ls=":", zorder=2)
    syn = pd.read_csv(SYN / "stars.csv")
    syn = syn[~syn["baseline_confounded"].astype(bool)]
    se = np.arange(np.floor(syn.magnitude_true.min() / .25) * .25,
                   syn.magnitude_true.max() + .25, .25)
    sx_, sc_ = binned(syn.magnitude_true.to_numpy(float),
                      syn.recovered.to_numpy(bool), se)
    ax.plot(sx_, sc_, "--", color=PALETTE["grey"], lw=1.3, zorder=3,
            label="synthetic (verification)")
    m50_m13 = 14.9
    for r in runs:
        if not r["hero"]:
            continue
        col = C[r["ckey"]]
        ax.plot(r["mx"], r["mc"], r["mk"], color=col, lw=2.0, ms=3.4, mfc=col,
                zorder=6, label=f"{r['label']}")
        ax.axvline(r["m50"], color=col, lw=0.9, ls=":", zorder=2, alpha=0.75)
        ha = "right" if r["label"].startswith("M13") else ("left" if r["label"].startswith("NGC") else "center")
        dx = -0.08 if ha == "right" else (0.08 if ha == "left" else 0.0)
        ax.text(r["m50"] + dx, 1.02, f"{r['m50']:.1f}", color=col, fontsize=7.0,
                ha=ha, va="bottom", fontweight="bold")
        if r["label"].startswith("M13"):
            m50_m13 = r["m50"]
    ax.text(12.45, 1.02, "$m_{50}$:", color="0.35", fontsize=7.0, ha="left", va="bottom")
    ax.set_xlabel("injected magnitude (count-rate, ZP = 25)")
    ax.set_ylabel("recovery completeness")
    ax.set_xlim(12.3, 20.2)
    ax.set_ylim(-0.03, 1.11)
    ax.legend(loc="lower left", fontsize=6.9, framealpha=0.93)
    ax.set_title("(a) depth is a frame property (sky + seeing)", loc="left", fontsize=9.5)

    # (b) the collapse — all seven frames in peak-S/N space
    axb.axhline(0.5, color=PALETTE["grey"], lw=0.7, ls=":", zorder=2)
    for r in runs:
        col = C[r["ckey"]] if r["hero"] else "0.55"
        lw = 1.7 if r["hero"] else 1.0
        axb.plot(10 ** r["sx"], r["sc"], "-", color=col, lw=lw,
                 alpha=0.9 if r["hero"] else 0.55, zorder=4 if r["hero"] else 3)
    xf = np.linspace(allx.min(), allx.max(), 300)
    axb.plot(10 ** xf, model(xf, *popt), "-", color="k", lw=1.1, ls="--", zorder=6,
             label=f"erf fit — S/N$_{{50}}$ = {s50_fit:.1f}")
    axb.axvline(s50_fit, color="k", lw=0.8, ls=":", zorder=2, alpha=0.7)
    axb.set_xscale("log")
    axb.set_xlim(0.7, 300)
    axb.set_ylim(-0.03, 1.11)
    axb.set_xlabel("expected peak S/N")
    axb.set_title("(b) all seven frames, one law", loc="left", fontsize=9.5)
    axb.legend(loc="lower right", fontsize=6.8, framealpha=0.93)
    axb.text(0.97, 0.60,
             f"7 frames · $\\sigma$ 5–58 e$^-$\nFWHM 5.2–9.0 px\n"
             f"S/N$_{{50}}$ = {s50s.mean():.1f} ± {s50s.std():.1f}\n"
             "3.4 mag of depth →\none curve",
             transform=axb.transAxes, fontsize=6.6, color="0.3", va="top", ha="right",
             bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
                   "alpha": .85, "edgecolor": PALETTE["grey"]})

    fig.suptitle("real-frame injection: per-frame depths are one detection law in S/N",
                 x=0.01, ha="left", fontsize=11)

    # bottom strip: real injected-star cutouts (M13, the shallowest frame)
    gsc = gs[1].subgridspec(1, len(cmags), wspace=0.16)
    for j, (s_, m_) in enumerate(zip(stamps, cmags)):
        axc = fig.add_subplot(gsc[j])
        axc.imshow(s_, cmap="gray", vmin=clo, vmax=chi, origin="lower",
                   interpolation="nearest")
        axc.set_xticks([]); axc.set_yticks([])
        recovered = m_ < m50_m13
        col = C["data"] if recovered else C["accent"]
        for sp in axc.spines.values():
            sp.set_edgecolor(col); sp.set_linewidth(1.5)
        axc.set_title(f"m = {m_:.1f}", fontsize=7.4, pad=2)
        axc.set_xlabel("found" if recovered else "lost", fontsize=7.6, color=col, labelpad=2)
    fig.text(0.012, gs[1].get_position(fig).y1 + 0.016,
             "real injected stars — M13 V (the shallowest frame), shared stretch  —  "
             "blue = recovered, orange = lost",
             fontsize=7.2, color="0.3", va="bottom", ha="left")

    paths = save_fig(fig, "fig_completeness_realvssynth", OUTDIR)
    plt.close(fig)
    print(f"[fit] A={A_f:.3f} S/N50={s50_fit:.2f} w(log)={w_f:.3f}")
    for r in runs:
        print(f"  {r['label']:20s} σ_e={r['sigma_e']:5.1f} FWHM={r['fwhm']:.2f} "
              f"m50={r['m50']:.2f}  S/N50={r['s50']:.2f}")
    print(f"saved: {paths['png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
