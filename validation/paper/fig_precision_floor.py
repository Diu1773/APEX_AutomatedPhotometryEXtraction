"""Figure — photometric precision floor (per-star RMS vs magnitude). PRELIMINARY:
single night, N=10 epochs.

Reproduces the ensemble-photometry method validation of Honeycutt 1992 (PASP
104,435): the per-star scatter across frames, after a differential-ensemble
zeropoint correction, plotted against magnitude. (Collier Cameron+2006 / Kovacs+
2005 are NOT claimed here — their point is the before/after-detrending floor drop
on thousands of epochs; that reproduction is deferred to the LC re-run with SYSREM.) Two things are checked
at once — (i) the empirical scatter follows the pipeline's OWN reported photon-noise
error bar (error model validated), and (ii) at the bright end the scatter hits a
systematic floor set by flat-field / PSF / scintillation residuals.

Data: fully-APEX-reduced M67 (Moravian C3-61000), r band, 10 frames, 1073 stars
matched by the forced-photometry master catalogue. Ensemble ZP per frame = median
residual of the brightest 40% of stars (Honeycutt differential ensemble). No
external catalogue, no code shared with any reference tool.

    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig_precision_floor.py
"""
from __future__ import annotations
import sys, glob, os
from pathlib import Path
REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

BASE = Path(r"E:\APEX_validation\reprocess\M67\result\step7_forced_phot")
OUTDIR = REPO / "validation" / "paper" / "figures"
FILT = "r"


def load_filter(filt):
    files = sorted(glob.glob(str(BASE / f"photometry_pp_Messier67-*-{filt}.fit.tsv")))
    frames = []
    for f in files:
        d = pd.read_csv(f, sep="\t")
        d = d[(d.bad_phot_flag == False) & (d.is_saturated == False)
              & np.isfinite(d.mag_inst) & np.isfinite(d.mag_err)]
        frames.append(d[["master_id", "mag_inst", "mag_err"]].copy())
    magw = pd.concat([fr.set_index("master_id")["mag_inst"] for fr in frames], axis=1)
    errw = pd.concat([fr.set_index("master_id")["mag_err"] for fr in frames], axis=1)
    nfr = magw.shape[1]
    # Honeycutt ensemble ZP: per-frame median residual of the bright 40%
    star_med = magw.median(axis=1)
    resid = magw.sub(star_med, axis=0)
    bright = star_med < star_med.quantile(0.40)
    zp = resid[bright].median(axis=0)
    magc = magw.sub(zp, axis=1)
    nobs = magc.notna().sum(axis=1)
    ok = nobs >= max(5, nfr - 3)
    ref = magc[ok].median(axis=1)
    rms = magc[ok].std(axis=1, ddof=1)
    err = errw[ok].median(axis=1)         # pipeline's own reported error
    return nfr, ref.to_numpy(), rms.to_numpy(), err.to_numpy()


def binned(x, y, edges, stat=np.median, nmin=8):
    bx, by = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x >= lo) & (x < hi) & np.isfinite(y)
        if m.sum() >= nmin:
            bx.append(0.5 * (lo + hi)); by.append(stat(y[m]))
    return np.array(bx), np.array(by)


def main() -> int:
    nfr, ref, rms, err = load_filter(FILT)
    k = np.isfinite(ref) & np.isfinite(rms) & (rms > 0)
    ref, rms, err = ref[k], rms[k], err[k]
    edges = np.arange(np.floor(ref.min()), np.ceil(ref.max()) + 0.5, 0.5)
    bx, brms = binned(ref, rms, edges)
    ex, berr = binned(ref, err, edges)
    floor = float(np.median(rms[ref < np.quantile(ref, 0.25)]))

    # extra floors for the caption/text
    floors = {}
    for f in ("g", "r", "i"):
        try:
            _, rf, rr, _ = load_filter(f)
            floors[f] = float(np.median(rr[rf < np.quantile(rf, 0.25)]))
        except Exception:
            floors[f] = float("nan")

    fig, ax = plt.subplots(1, 1, figsize=(DOUBLE_COL * 0.66, 3.6))
    ax.plot(ref, 1000 * rms, "o", color=C["data"], ms=2.6, alpha=0.35, mew=0,
            zorder=2, label="per-star RMS (10 frames)")
    ax.plot(bx, 1000 * brms, "-", color=C["accent"], lw=2.0, zorder=5,
            label="binned RMS (empirical)")
    ax.plot(ex, 1000 * berr, "--", color=C["reference"], lw=1.8, zorder=4,
            label="reported error (photon-noise model)")
    ax.axhline(1000 * floor, color=PALETTE["grey"], lw=1.0, ls=":", zorder=3)
    ax.text(ref.min() + 0.1, 1000 * floor * 1.12,
            f"systematic floor ≈ {1000*floor:.0f} mmag",
            fontsize=7.2, color="0.3", va="bottom")

    ax.set_yscale("log")
    ax.set_xlabel("instrumental magnitude (count-rate, ZP = 25)")
    ax.set_ylabel("photometric scatter (mmag)")
    ax.set_ylim(2, 800)
    ax.legend(loc="upper left", fontsize=7.0, framealpha=0.9)
    ax.set_title("scatter tracks the photon-noise model down to a mmag floor",
                 loc="left")
    ax.text(0.97, 0.05,
            f"M67 {FILT} · {nfr} frames · {k.sum()} stars\n"
            f"floors: g {1000*floors['g']:.0f} · r {1000*floors['r']:.0f} · "
            f"i {1000*floors['i']:.0f} mmag (preliminary)",
            transform=ax.transAxes, va="bottom", ha="right", fontsize=6.6,
            color="0.3",
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
                  "alpha": .85, "edgecolor": PALETTE["grey"]})

    paths = save_fig(fig, "fig_precision_floor", OUTDIR)
    plt.close(fig)
    print(f"[precision] M67 {FILT}: floor={1000*floor:.1f}mmag  "
          f"floors g/r/i={1000*floors['g']:.1f}/{1000*floors['r']:.1f}/"
          f"{1000*floors['i']:.1f}  nstars={k.sum()}  frames={nfr}")
    print(f"saved: {paths['png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
