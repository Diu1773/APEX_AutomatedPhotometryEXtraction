"""Aperture/annulus sweep on a real frame — the experiment that identified the
B-band faint bias driver (sky overestimation by a PSF-wing/neighbour-halo
contaminated annulus at 4-6xFWHM; NGC 6811, 2026-07).

Discrimination logic:
  * sky-overestimation  -> faint drift GROWS with r_ap (n_pix grows) and
    SHRINKS when the annulus moves outward (cleaner sky).
  * decentering loss    -> faint drift SHRINKS with larger r_ap.

Measured (isolated stars, drift = median resid faint>=17.5 minus bright<15):
  pp_...0008-B (480 s): annulus 4-6 +0.031  ->  6-9 +0.005
  pp_...0002-B (240 s): annulus 4-6 +0.087  ->  6-9 +0.044  (shorter exposure,
  larger relative bias — exactly the delta_sky*n_pix/F prediction)

Run (data SSD mounted):
  .venv-deploy/Scripts/python.exe validation/b_faint_annulus_sweep.py
  ... --frame pp_NGC6811-0002-B_20260611_2.fit --fwhm 7.229
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO))

from astropy.io import fits
from scipy.spatial import cKDTree

from apex.utils.photometry_utils import phot_vectorized

MAD = lambda v: 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v)))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--frame", default="pp_NGC6811-0008-B_20260611_2.fit")
    ap.add_argument("--fwhm", type=float, default=7.303948,
                    help="frame FWHM px (see step7 centering_stats.csv)")
    ap.add_argument("--data-dir", default=r"E:/observed_Analysis/NGC6811/pp")
    ap.add_argument("--band", default="B", choices=["B", "V", "R"])
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    result = data_dir / "result"
    fits_path = data_dir / args.frame
    tsv = result / "step7_forced_phot" / f"photometry_{args.frame}.tsv"
    cal_path = result / "cmd_zeropoint" / "gaia_sdss_calibrator_by_ID.csv"
    for p in (fits_path, tsv, cal_path):
        if not p.exists():
            print(f"[ABORT] missing: {p}")
            return 2

    band = args.band
    other = "V" if band != "V" else "B"
    df = pd.read_csv(tsv, sep="\t").merge(
        pd.read_csv(cal_path)[["ID", f"ref_{band}", f"ref_{other}", "x_pix", "y_pix"]],
        on="ID", how="inner", suffixes=("", "_cal"),
    )
    img = fits.getdata(fits_path).astype(float)
    gain = float(pd.to_numeric(df["gain_e_per_adu"], errors="coerce").median())
    rn = float(pd.to_numeric(df["rdnoise_e"], errors="coerce").median())
    expt = float(pd.to_numeric(df["exptime"], errors="coerce").median())
    fwhm = float(args.fwhm)

    xy = df[["x_pix", "y_pix"]].to_numpy(float)
    d2, _ = cKDTree(xy).query(xy, k=2)
    iso_mask = d2[:, 1] > float(np.nanmedian(d2[:, 1]))
    pos = df[["x", "y"]].to_numpy(float)
    ref = df[f"ref_{band}"].to_numpy(float)
    cref = (df[f"ref_{'B' if band != 'V' else 'V'}"] - df[f"ref_{other}"]).to_numpy(float) \
        if band == "B" else (df[f"ref_{other}"] - df[f"ref_{band}"]).to_numpy(float)
    print(f"{args.frame}: N={len(df)} gain={gain:.4f} RN={rn:.2f} exptime={expt:.0f}s FWHM={fwhm:.2f}px")

    def run(r_ap_f, r_in_f, r_out_f):
        flux_e, _, snr, sky, _, _ = phot_vectorized(
            img, pos, r_ap=r_ap_f * fwhm, r_in=r_in_f * fwhm, r_out=r_out_f * fwhm,
            gain=gain, rn_param_e=rn, sky_sigma_mode="local",
            sat_adu=65000.0, datamax_adu=60000.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            mag = -2.5 * np.log10(np.where(flux_e > 0, flux_e, np.nan) / expt)
        delta = ref - mag
        good = np.isfinite(delta) & np.isfinite(cref) & np.isfinite(snr)
        m = good & (snr > 50)
        c = (np.nan, np.nan)
        for _ in range(4):
            if m.sum() < 20:
                return None
            c = np.polyfit(cref[m], delta[m], 1)
            r = delta - np.polyval(c, cref)
            s = MAD(r[m])
            m2 = m & (np.abs(r - np.nanmedian(r[m])) < 3 * max(s, 1e-6))
            if m2.sum() == m.sum():
                break
            m = m2
        r = delta - np.polyval(c, cref)
        out = {}
        for tag, sel in (("all", good), ("iso", good & iso_mask)):
            b, f = sel & (ref < 15), sel & (ref >= 17.5)
            out[tag] = (float(np.nanmedian(r[f]) - np.nanmedian(r[b]))
                        if b.sum() >= 10 and f.sum() >= 10 else np.nan)
        out["sky"] = float(np.nanmedian(sky))
        return out

    print("\n--- r_ap sweep (annulus fixed 4-6 FWHM) ---")
    for rr in (0.5, 0.65, 0.8, 1.0, 1.25):
        o = run(rr, 4.0, 6.0)
        if o:
            print(f"  r_ap={rr:4.2f}xFWHM  drift all={o['all']:+.4f}  iso={o['iso']:+.4f}  sky={o['sky']:.2f}")
    print("--- annulus sweep (r_ap fixed 0.8 FWHM) ---")
    for ri, ro in ((2.5, 4.0), (4.0, 6.0), (6.0, 9.0), (9.0, 13.0)):
        o = run(0.8, ri, ro)
        if o:
            print(f"  ann {ri:4.1f}-{ro:<4.1f}xFWHM  drift all={o['all']:+.4f}  iso={o['iso']:+.4f}  sky={o['sky']:.2f}")
    print("\nInterpretation: drift grows with r_ap AND shrinks with an outward "
          "annulus => sky overestimation; drift shrinking with r_ap => decentering.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
