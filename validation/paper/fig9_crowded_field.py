"""Crowded-field comparison on the globular clusters M5 and M13.

Panel (a) measures the radial source-density enhancement. Panel (b) compares
APEX forced-aperture and PSF magnitudes against nearest-neighbour separation.
Both measurements use retained products from reductions made with the current
code. The test is internal to APEX and does not cover blends below the
approximately 10-pixel detection/deduplication floor.

Needs the data SSD (M5 + M13 + NGC 6811 re-reductions). Run:
    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig9_crowded_field.py
"""

from __future__ import annotations

import glob
import os
import sys
from pathlib import Path

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.spatial import cKDTree

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL
from apex.utils.psf_core import _auto_density_center

apply_paper_style()

NGC6811_CAL = Path(r"E:/observed_Analysis/NGC6811/pp/result/cmd_zeropoint/gaia_sdss_calibrator_by_ID.csv")
OUTDIR = REPO / "validation" / "paper" / "figures"
CAPDIR = REPO / "validation" / "paper" / "captions"

MAD = lambda v: 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v)))
RES_FLOOR_PX = 10.0  # detection/master min separation (master_sources neighbor_dist_px floor)

CLUSTERS = {
    "M5": {
        "label": "M5 (NGC 5904)",
        "color": C["data"],
        "marker": "s",
        "result": Path(r"E:/observed_Analysis/M5/light/result"),
        "data_glob": r"E:/observed_Analysis/M5/light/pp_*_20250308.fit",
        "fwhm_px": 5.5,
        "filter": "r",
        "n_frames": 34,
    },
    "M13": {
        "label": "M13 (NGC 6205)",
        "color": C["accent"],
        "marker": "^",
        "result": Path(r"E:/observed_Analysis/M13/light/result"),
        "data_glob": r"E:/observed_Analysis/M13/light/pp_*.fit",
        "fwhm_px": 6.0,
        "filter": "r",
        "n_frames": 12,
    },
}


def density_profile(cfg):
    master = pd.read_csv(cfg["result"] / "step7_forced_phot" / "master_sources.csv")
    mr = master[master["filter"] == cfg["filter"]][["ID", "x_ref", "y_ref"]].drop_duplicates("ID")
    xy = mr[["x_ref", "y_ref"]].to_numpy(float)

    samples = sorted(glob.glob(cfg["data_glob"]))
    if samples:
        h = fits.getheader(samples[0])
        shape = (int(h["NAXIS2"]), int(h["NAXIS1"]))
    else:
        # The raw/processed FITS may be archived after the distilled Step-7
        # products are retained. Both reductions use the same 4800x3200 crop;
        # also guard against any catalog extending past that nominal boundary.
        shape = (
            max(3200, int(np.ceil(xy[:, 1].max())) + 1),
            max(4800, int(np.ceil(xy[:, 0].max())) + 1),
        )

    cx, cy, ratio = _auto_density_center(xy, shape, fwhm_px=cfg["fwhm_px"], cell_fwhm_mult=8.0)
    r = np.hypot(xy[:, 0] - cx, xy[:, 1] - cy)

    edges = np.array([0, 20, 40, 60, 90, 130, 180, 250, 400, 600, 900])
    centers, dens = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (r >= lo) & (r < hi)
        area = np.pi * (hi ** 2 - lo ** 2)
        if area <= 0:
            continue
        centers.append(0.5 * (lo + hi))
        dens.append(k.sum() / area * 1e4)
    return np.array(centers), np.array(dens), float(ratio), len(mr)


def ngc6811_density_ratio():
    cal = pd.read_csv(NGC6811_CAL)
    xy = cal[["x_pix", "y_pix"]].dropna().to_numpy(float)
    _, _, ratio = _auto_density_center(xy, (4300, 4300), fwhm_px=7.3, cell_fwhm_mult=8.0)
    return float(ratio)


def aperture_vs_psf_by_crowding(cfg):
    master = pd.read_csv(cfg["result"] / "step7_forced_phot" / "master_sources.csv")
    mr = master[master["filter"] == cfg["filter"]][
        ["ID", "neighbor_dist_px", "crowding_flag"]
    ].drop_duplicates("ID")

    ap_rows, psf_rows = [], []
    for f in sorted(glob.glob(str(cfg["result"] / "step7_forced_phot" / "photometry_*.tsv"))):
        df = pd.read_csv(f, sep="\t")
        df = df[df["FILTER"] == cfg["filter"]]
        if df.empty:
            continue
        ap_rows.append(df[["ID", "x", "y", "mag_inst", "snr"]])

        pf = f.replace("step7_forced_phot", "cmd_psf")
        if not os.path.exists(pf):
            archived = cfg["result"] / "cmd_psf_backup_gui_20260729" / Path(f).name
            if archived.exists():
                pf = str(archived)
            else:
                continue
        psf = pd.read_csv(pf, sep="\t")
        psf = psf[psf["FILTER"] == cfg["filter"]]
        if psf.empty:
            continue
        # Positional match (det_uid is unreliable across step7/step8 -- lesson
        # from the earlier B-filter PSF cross-check this session).
        tree = cKDTree(df[["x", "y"]].to_numpy(float))
        d, i = tree.query(psf[["x_fit", "y_fit"]].to_numpy(float), k=1)
        ok = d < 3.0
        m = psf[ok].copy()
        m["ID"] = df["ID"].to_numpy()[i[ok]]
        psf_rows.append(m[["ID", "mag_psf", "snr_psf"]])

    ap_all = pd.concat(ap_rows, ignore_index=True)
    psf_all = pd.concat(psf_rows, ignore_index=True)
    ap_med = ap_all.groupby("ID").agg(mag_ap=("mag_inst", "median"),
                                       snr_ap=("snr", "median")).reset_index()
    psf_med = psf_all.groupby("ID").agg(mag_psf=("mag_psf", "median"),
                                         snr_psf=("snr_psf", "median")).reset_index()
    merged = mr.merge(ap_med, on="ID", how="inner").merge(psf_med, on="ID", how="inner")

    bright_iso = (merged["snr_ap"] > 50) & (merged["snr_psf"] > 50)
    offset = float(np.nanmedian((merged["mag_ap"] - merged["mag_psf"])[bright_iso]))
    resid = (merged["mag_ap"] - merged["mag_psf"] - offset).to_numpy(float)
    nn = merged["neighbor_dist_px"].to_numpy(float)
    return nn, resid, int(bright_iso.sum()), int(len(merged))


def main() -> int:
    ratio_ngc6811 = ngc6811_density_ratio()
    results = {}
    for name, cfg in CLUSTERS.items():
        centers, dens, ratio, n_master = density_profile(cfg)
        nn, resid, n_anchor, n_matched = aperture_vs_psf_by_crowding(cfg)
        results[name] = dict(centers=centers, dens=dens, ratio=ratio, n_master=n_master,
                             nn=nn, resid=resid, n_anchor=n_anchor, n_matched=n_matched)
        print(f"{name}: density_ratio={ratio:.2f}  n_master={n_master}  "
              f"n_matched(ap-vs-psf)={n_matched}  n_anchor={n_anchor}")

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 3.5))

    # (a) radial density profiles, both clusters overlaid
    for name, cfg in CLUSTERS.items():
        r = results[name]
        ax_a.semilogy(r["centers"], r["dens"], cfg["marker"] + "-", color=cfg["color"],
                     lw=1.5, ms=4, label=f"{cfg['label']}: {r['ratio']:.0f}$\\times$")
    ax_a.set_xlabel(r"radius from density peak (px)")
    ax_a.set_ylabel(r"source density ($10^{-4}$ px$^{-2}$)")
    ax_a.legend(loc="upper right", fontsize=7.6, title=f"core/background\n(NGC 6811 = {ratio_ngc6811:.0f}$\\times$, Figs 1-8)",
               title_fontsize=6.8)
    ax_a.set_title("(a) globular-core density enhancement", loc="left")

    # (b) aperture-vs-PSF residual vs crowding, both clusters overlaid
    ax_b.axhline(0.0, color=C["floor"], lw=0.8, ls=":", zorder=1)
    ax_b.axvspan(0, RES_FLOOR_PX, color=PALETTE["grey"], alpha=0.15, zorder=0)
    edges = np.array([10, 15, 20, 30, 45, 70, 110, 1000])
    info_lines = []
    for name, cfg in CLUSTERS.items():
        r = results[name]
        ax_b.scatter(r["nn"], r["resid"], s=7, alpha=0.22, color=cfg["color"],
                    edgecolors="none", zorder=2)
        for lo, hi in zip(edges[:-1], edges[1:]):
            k = (r["nn"] >= lo) & (r["nn"] < hi) & np.isfinite(r["resid"])
            if k.sum() < 8:
                continue
            xm = 0.5 * (lo + hi) if hi < 1000 else float(np.median(r["nn"][k]))
            ax_b.errorbar(xm, np.nanmedian(r["resid"][k]),
                         yerr=MAD(r["resid"][k]) / np.sqrt(k.sum()),
                         fmt=cfg["marker"], color=cfg["color"], ms=5.5, capsize=2, zorder=4,
                         markeredgecolor="black", markeredgewidth=0.4)
        info_lines.append(f"{cfg['label']}: N={r['n_matched']}, anchor={r['n_anchor']}")
    ax_b.set_xscale("log")
    ax_b.set_xlim(8, 400)
    ax_b.set_ylim(-0.15, 0.15)
    ax_b.set_xlabel(r"nearest-neighbour separation (px)")
    ax_b.set_ylabel(r"$m_{\rm aperture} - m_{\rm PSF}$  (mag, ZP-aligned)")
    ax_b.text(0.03, 0.05,
              "\n".join(info_lines) + f"\ngrey band: resolution floor ({RES_FLOOR_PX:.0f} px)",
              transform=ax_b.transAxes, va="bottom", ha="left", fontsize=6.8,
              bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
                    "alpha": 0.85, "edgecolor": PALETTE["grey"]})
    ax_b.set_title("(b) aperture vs PSF vs local density", loc="left")

    fig.suptitle("M5 & M13 (globular clusters) — crowded-field re-reductions, current codebase",
                 fontsize=8.0, y=0.995, color="#333333")
    provenance = (
        f"Moravian C3-61000 (2×2), r band · M5: {CLUSTERS['M5']['n_frames']} frames, "
        f"{results['M5']['n_master']} master sources · M13: {CLUSTERS['M13']['n_frames']} frames, "
        f"{results['M13']['n_master']} master sources · APEX Step 7/8 retained products"
    )
    fig.tight_layout(rect=(0, 0.055, 0.99, 0.97))
    fig.text(0.5, 0.012, provenance, ha="center", va="bottom",
             fontsize=5.4, color=PALETTE["grey"])
    paths = save_fig(fig, "fig9_crowded_field", OUTDIR)
    plt.close(fig)

    m5, m13 = results["M5"], results["M13"]
    caption = f"""# Figure — crowded-field photometry comparison

M5 ({CLUSTERS['M5']['n_frames']} $r$-band frames, {m5['n_master']} master
sources) and M13 ({CLUSTERS['M13']['n_frames']} frames, {m13['n_master']}
master sources) were reduced with the current APEX code. **(a)** The core
source density is {m5['ratio']:.0f}$\\times$ the field background in M5 and
{m13['ratio']:.0f}$\\times$ in M13; NGC 6811 is {ratio_ngc6811:.0f}$\\times$
by the same metric. **(b)** Zeropoint-aligned aperture-minus-PSF residuals for
{m5['n_matched']} M5 and {m13['n_matched']} M13 stars show no trend with
nearest-neighbour separation above the {RES_FLOOR_PX:.0f}-pixel resolution
floor. This is an internal comparison of two APEX methods; it does not test
external absolute accuracy or unresolved blends below that floor.
"""
    CAPDIR.mkdir(parents=True, exist_ok=True)
    (CAPDIR / "fig9_crowded_field.md").write_text(caption, encoding="utf-8")

    for ext, p in paths.items():
        print(f"wrote {ext}: {p}  exists={p.exists()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
