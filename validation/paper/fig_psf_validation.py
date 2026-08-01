"""Figure — PSF photometry against forced aperture photometry, across cameras.

Method identical to validation/psf_crossinstrument/report_psf_cross.py: join the
step-8 PSF table to the step-7 forced-aperture table on det_uid, keep stars with
SNR > 20 in both and no crowding/saturation flag, take D = mag_psf - mag_inst,
remove each frame's median (the EPSF-normalisation offset that the frame
zeropoint absorbs at calibration), and judge by the MAD of what remains.

(a) star-level distributions of (D - frame median) for six datasets on three
    cameras, crowded globular to 2.0-degree wide field.
(b) the M45 full-frame set: radial profile of the same quantity from field
    centre to corner — the single-EPSF-per-frame adequacy check.

Distilled arrays are written to data_psf_validation/ so the figure can be
rebuilt without the E: trees.

    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig_psf_validation.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

DATADIR = REPO / "validation" / "paper" / "data_psf_validation"
OUTDIR = REPO / "validation" / "paper" / "figures"
RNG = np.random.default_rng(20260802)
CAP = 30000          # per-dataset sample cap stored/plotted

REPRO = Path(r"E:\APEX_validation\reprocess")
XI = Path(r"E:\APEX_validation\psf_crossinstrument")

# key, result tree, label (two lines: field / camera), pixel scale "/px
SETS = [
    ("M13",     REPRO / "M13" / "result",             "M13\nglobular",        "C3-61000", 0.395),
    ("NGC6811", REPRO / "NGC6811" / "result",         "NGC 6811\nopen",       "C3-61000", 0.395),
    ("M67",     XI / "m67_lco" / "result",            "M67\nopen",            "QHY600",   0.74),
    ("Proxima", XI / "qhy600" / "result",             "Proxima\nfield",       "QHY600",   0.74),
    ("M45",     XI / "m45_wide" / "result",           "M45\n2.0$^\\circ$ wide", "QHY600", 0.74),
    ("N5985",   XI / "sinistro" / "result",           "NGC 5985\nfield",      "Sinistro", 0.389),
]

PSF_COLS = ["det_uid", "mag_psf", "snr_psf", "crowding_unreliable_psf",
            "saturated_psf", "x_fit", "y_fit"]
AP_COLS = ["det_uid", "mag_inst", "snr"]


def load_dataset(rdir: Path):
    """Pooled (D - frame median) plus positions; per-frame medians kept."""
    frames = sorted((rdir / "cmd_psf").glob("photometry_*.tsv"))
    pooled, meds, xy = [], [], []
    for f in frames:
        fname = f.name[len("photometry_"):-len(".tsv")]
        ap = rdir / "step7_forced_phot" / f"photometry_{fname}.tsv"
        if not ap.exists():
            continue
        p = pd.read_csv(f, sep="\t", usecols=PSF_COLS)
        a = pd.read_csv(ap, sep="\t", usecols=AP_COLS)
        m = p.dropna(subset=["det_uid"]).merge(a.dropna(subset=["det_uid"]), on="det_uid")
        ok = m[np.isfinite(m.mag_psf) & np.isfinite(m.mag_inst)
               & (m.snr_psf > 20) & (m.snr > 20)
               & (m.crowding_unreliable_psf == 0) & (m.saturated_psf == 0)]
        if len(ok) < 5:
            continue
        d = (ok.mag_psf - ok.mag_inst).to_numpy(float)
        med = float(np.median(d))
        meds.append(med)
        pooled.append(d - med)
        xy.append(np.c_[ok.x_fit.to_numpy(float), ok.y_fit.to_numpy(float)])
    dd = np.concatenate(pooled)
    xy = np.vstack(xy)
    return dd, xy, meds, len(pooled)


def main() -> int:
    DATADIR.mkdir(parents=True, exist_ok=True)
    stats, store = [], {}
    for key, rdir, label, cam, scale in SETS:
        dd, xy, meds, nfr = load_dataset(rdir)
        mad = float(1.4826 * np.median(np.abs(dd)))
        stats.append(dict(key=key, camera=cam, n_frames=nfr, n_stars=int(dd.size),
                          mad_mag=mad, offset_med=float(np.median(meds)),
                          offset_spread=float(np.std(meds))))
        idx = RNG.permutation(dd.size)[:CAP]
        store[key] = dd[idx]
        if key == "M45":
            cx = (xy[:, 0].min() + xy[:, 0].max()) / 2
            cy = (xy[:, 1].min() + xy[:, 1].max()) / 2
            r_deg = np.hypot(xy[:, 0] - cx, xy[:, 1] - cy) * scale / 3600.0
            store["M45_r"] = r_deg[idx]
            m45_full = (dd, r_deg)
        print(f"{key:<8} frames={nfr:>2} stars={dd.size:>7,} MAD={mad*1000:6.1f} mmag "
              f"offset={np.median(meds):+.3f}±{np.std(meds):.3f}", flush=True)

    np.savez_compressed(DATADIR / "psf_datasets.npz", **store)
    (DATADIR / "summary.json").write_text(json.dumps(stats, indent=1), encoding="utf-8")

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COL, 2.7),
                             gridspec_kw={"width_ratios": [1.5, 1.0]})

    # (a) per-dataset distributions.  Stars beyond +-0.2 mag are dropped, not
    # clipped — clipping piles them at the boundary and draws a false bulge.
    ax = axes[0]
    pos = np.arange(len(SETS))
    data, trimmed = [], []
    for k, *_ in SETS:
        v = store[k]
        keep = np.abs(v) <= 0.2
        data.append(v[keep])
        trimmed.append(1 - keep.mean())
    vp = ax.violinplot(data, positions=pos, widths=0.75, showextrema=False)
    cam_col = {"C3-61000": C["data"], "QHY600": C["reference"], "Sinistro": C["accent"]}
    for body, (key, _r, _l, cam, _s) in zip(vp["bodies"], SETS):
        body.set_facecolor(cam_col[cam]); body.set_alpha(0.55)
        body.set_edgecolor("none")
    for i, s in enumerate(stats):
        ax.text(pos[i], 0.175, f"{s['mad_mag']*1000:.0f}", ha="center", fontsize=7.2)
        ax.text(pos[i], 0.142, f"{s['n_frames']} fr", ha="center", fontsize=5.4,
                color=PALETTE["grey"])
        ax.text(pos[i], 0.118, f"{s['n_stars']:,}", ha="center", fontsize=5.4,
                color=PALETTE["grey"])
    ax.text(-0.62, 0.175, "MAD\n(mmag)", fontsize=5.4, color=PALETTE["grey"],
            ha="center", va="center")
    ax.axhline(0, color=PALETTE["black"], lw=0.6)
    ax.set_xticks(pos, [l for _k, _r, l, _c, _s in SETS], fontsize=6.6)
    ax.set_xlim(-0.95, len(SETS) - 0.4)
    ax.set_ylabel(r"$m_{\rm PSF}-m_{\rm ap}$ (mag)")
    ax.set_ylim(-0.21, 0.21)
    ax.text(0.015, 0.03, "frame median removed", transform=ax.transAxes,
            fontsize=5.8, color=PALETTE["grey"])
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c, alpha=0.55)
               for c in cam_col.values()]
    ax.legend(handles, list(cam_col), fontsize=6.2, loc="lower right", ncol=1,
              handlelength=1.1, labelspacing=0.25, frameon=False,
              borderaxespad=0.2)
    ax.set_title("(a) star-level agreement, three cameras", loc="left", fontsize=9)
    print("trim >0.2 mag:", " ".join(f"{k}:{t*100:.1f}%"
                                     for (k, *_), t in zip(SETS, trimmed)))

    # (b) M45 radial profile
    ax = axes[1]
    dd, r_deg = m45_full
    bins = np.linspace(0, r_deg.max(), 9)
    mid = 0.5 * (bins[1:] + bins[:-1])
    med = [float(np.median(dd[(r_deg >= a) & (r_deg < b)])) for a, b in zip(bins, bins[1:])]
    mad = [float(1.4826 * np.median(np.abs(dd[(r_deg >= a) & (r_deg < b)]
                                            - m))) for (a, b), m in zip(zip(bins, bins[1:]), med)]
    sub = RNG.permutation(dd.size)[:12000]
    ax.plot(r_deg[sub], dd[sub], ".", ms=1.0, color=PALETTE["grey"], alpha=0.25,
            rasterized=True)
    ax.errorbar(mid, med, yerr=mad, fmt="o-", ms=3.2, lw=1.1, capsize=2,
                color=C["model"], label="binned median $\\pm$ MAD")
    inner_span = (max(med[:-1]) - min(med[:-1])) * 1000
    ax.annotate(f"corner bin\n{med[-1]*1000:+.0f} mmag", (mid[-1], med[-1]),
                textcoords="offset points", xytext=(-52, 26), fontsize=6.2,
                color=C["model"],
                arrowprops=dict(arrowstyle="-", lw=0.7, color=C["model"]))
    ax.axhline(0, color=PALETTE["black"], lw=0.6)
    ax.set_xlabel("distance from field centre (deg)")
    ax.set_ylabel(r"$m_{\rm PSF}-m_{\rm ap}$ (mag)")
    ax.set_ylim(-0.14, 0.14)
    ax.legend(fontsize=6.4, loc="upper left")
    ax.set_title("(b) M45 $2.0^\\circ$ field, one EPSF per frame",
                 loc="left", fontsize=9)
    print(f"\nM45 반경: 모서리 제외 span {inner_span:.1f} mmag, "
          f"모서리 bin {med[-1]*1000:+.1f} mmag")

    fig.tight_layout(pad=0.5, w_pad=1.5)
    for ext, p in save_fig(fig, "fig_psf_validation", OUTDIR).items():
        print(f"[saved] {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
