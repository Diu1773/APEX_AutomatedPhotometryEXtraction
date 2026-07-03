"""Crowding-vs-drift discrimination for the B-filter faint bias (NGC 6811).

Answers: is the residual faint-end drift (B +0.03 / V -0.02 vs Gaia-transformed
reference) concentrated in CROWDED stars? If yes, the driver is neighbour
contamination — Gaia BP/RP prism windows (fix: C* cut, apex/utils/gaia_quality)
and/or APEX aperture blending. If the drift is the same for isolated stars,
crowding is exonerated and the remainder points at the transformation itself.

Requires the external data SSD (default path on E:). Run when mounted:
    .venv-deploy/Scripts/python.exe validation/b_faint_crowding_check.py
    # or with an explicit calibrator CSV:
    ... b_faint_crowding_check.py --csv <path to gaia_sdss_calibrator_by_ID.csv>

Self-contained: numpy/pandas/scipy only; re-fits the linear ZP+CT on the
reference color internally (does not depend on the committed coefficients).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

DEFAULT_CSV = (
    r"E:\observed_Analysis\NGC6811\pp\result\cmd_zeropoint\gaia_sdss_calibrator_by_ID.csv"
)

MAD = lambda v: 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v)))


def robust_linfit(x, y, w, clip=3.0, iters=5):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    coef = (np.nan, np.nan)
    for _ in range(iters):
        if m.sum() < 10:
            break
        coef = np.polyfit(x[m], y[m], 1, w=np.sqrt(w[m]))
        r = y - np.polyval(coef, x)
        s = MAD(r[m])
        m2 = m & (np.abs(r - np.nanmedian(r[m])) < clip * max(s, 1e-6))
        if m2.sum() == m.sum():
            break
        m = m2
    return coef


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", default=DEFAULT_CSV)
    ap.add_argument("--snr-cut", type=float, default=20.0, help="fit-side SNR cut")
    args = ap.parse_args()

    path = Path(args.csv)
    if not path.exists():
        print(f"[ABORT] calibrator CSV not found: {path}")
        print("        Mount the data SSD (or pass --csv) and rerun.")
        return 2

    df = pd.read_csv(path)
    x = pd.to_numeric(df["x_pix"], errors="coerce").to_numpy(float)
    y = pd.to_numeric(df["y_pix"], errors="coerce").to_numpy(float)
    ok_xy = np.isfinite(x) & np.isfinite(y)
    nn_dist = np.full(len(df), np.nan)
    if ok_xy.sum() > 10:
        tree = cKDTree(np.column_stack([x[ok_xy], y[ok_xy]]))
        d, _ = tree.query(np.column_stack([x[ok_xy], y[ok_xy]]), k=2)
        nn_dist[ok_xy] = d[:, 1]  # distance to nearest OTHER star
    med_nn = float(np.nanmedian(nn_dist))
    crowded = nn_dist < med_nn  # bottom half = crowded, top half = isolated
    print(f"N={len(df)}  median NN distance = {med_nn:.1f} px "
          f"(crowded = closer than this)")

    bands = [b for b in ("B", "V", "R") if f"ref_{b}" in df.columns]
    pair = {"B": ("B", "V"), "V": ("B", "V"), "R": ("V", "R")}
    for band in bands:
        fa, fb = pair[band]
        if f"ref_{fa}" not in df.columns or f"ref_{fb}" not in df.columns:
            continue
        cref = df[f"ref_{fa}"].to_numpy(float) - df[f"ref_{fb}"].to_numpy(float)
        delta = df[f"ref_{band}"].to_numpy(float) - df[f"mag_inst_{band}"].to_numpy(float)
        err = pd.to_numeric(df[f"mag_inst_err_{band}"], errors="coerce").to_numpy(float)
        w = np.where(np.isfinite(err) & (err > 0), 1.0 / err**2, np.nan)
        snr = pd.to_numeric(df[f"snr_{band}"], errors="coerce").to_numpy(float)
        coef = robust_linfit(cref[snr >= args.snr_cut], delta[snr >= args.snr_cut],
                             w[snr >= args.snr_cut])
        resid = delta - np.polyval(coef, cref)
        ref = df[f"ref_{band}"].to_numpy(float)

        print(f"\n===== {band} (zp={coef[1]:+.4f}, ct={coef[0]:+.4f}) =====")
        print(f"  {'mag bin':>12} | {'isolated med (N)':>22} | {'crowded med (N)':>22} | diff")
        for lo, hi in [(13, 16), (16, 17), (17, 17.5), (17.5, 18), (18, 19)]:
            k = np.isfinite(resid) & (ref >= lo) & (ref < hi)
            iso, crw = k & ~crowded, k & crowded
            if iso.sum() < 10 or crw.sum() < 10:
                continue
            mi, mc = np.nanmedian(resid[iso]), np.nanmedian(resid[crw])
            print(f"  {lo:5.1f}-{hi:4.1f}   | {mi:+.4f}  (N={iso.sum():4d})      "
                  f"| {mc:+.4f}  (N={crw.sum():4d})      | {mc - mi:+.4f}")
        bright = np.isfinite(resid) & (ref < 16)
        faint = np.isfinite(resid) & (ref >= 17.5)
        for tag, sel in (("isolated", ~crowded), ("crowded ", crowded)):
            b_med = np.nanmedian(resid[bright & sel])
            f_med = np.nanmedian(resid[faint & sel])
            print(f"  bright->faint drift [{tag}] = {f_med - b_med:+.4f} mag")

    # If ruwe / excess-factor columns exist (post step-6 requery), also show
    # the drift after the Gaia quality cut for comparison.
    try:
        from apex.utils.gaia_quality import gaia_quality_mask  # repo import if run from root

        if "phot_bp_rp_excess_factor" in df.columns or "ruwe" in df.columns:
            mq = gaia_quality_mask(df)
            print(f"\nGaia quality cut available: keeps {int(mq.sum())}/{len(df)} "
                  f"— rerun the tables above on df[mq] to quantify the C* effect.")
    except ImportError:
        pass

    print("\nVerdict guide: |crowded - isolated| >~ 0.01-0.02 mag at faint bins "
          "=> contamination-driven (C*/blending); ~0 => transformation-side.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
