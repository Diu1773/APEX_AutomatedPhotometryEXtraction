"""Cross-match two CMD tables (median_by_ID_filter_wide_cmd.csv) by sky position
and report per-band calibrated-magnitude agreement (median offset + robust sigma).

Used to validate the full-APEX (Step 0) reprocessing against a reference
result (e.g. AIPPI-preprocessed). Independent runs assign different integer
IDs, so matching is positional (RA/Dec), not by ID.

    python scripts/compare_cmd.py <new_cmd.csv> <ref_cmd.csv> [--tol-arcsec 1.0] [--json out.json]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u


def robust(d: np.ndarray):
    d = d[np.isfinite(d)]
    if d.size == 0:
        return 0, float("nan"), float("nan")
    med = float(np.median(d))
    sig = float(np.median(np.abs(d - med)) * 1.4826)
    return int(d.size), med, sig


def compare(new_csv: str, ref_csv: str, tol_arcsec: float = 1.0) -> dict:
    new = pd.read_csv(new_csv)
    ref = pd.read_csv(ref_csv)
    cn = SkyCoord(new.ra_deg.values * u.deg, new.dec_deg.values * u.deg)
    cr = SkyCoord(ref.ra_deg.values * u.deg, ref.dec_deg.values * u.deg)
    idx, sep, _ = cn.match_to_catalog_sky(cr)
    m = sep.arcsec < tol_arcsec
    N = new[m].reset_index(drop=True)
    R = ref.iloc[idx[m]].reset_index(drop=True)
    out = {"n_new": int(len(new)), "n_ref": int(len(ref)),
           "n_matched": int(m.sum()), "tol_arcsec": tol_arcsec, "bands": {}}
    bands = [b for b in ("B", "V", "R", "g", "r", "i", "z", "I")
             if f"mag_cal_{b}" in N.columns and f"mag_cal_{b}" in R.columns]
    for b in bands:
        n, med, sig = robust((N[f"mag_cal_{b}"] - R[f"mag_cal_{b}"]).values)
        out["bands"][b] = {"n": n, "median_dmag": med, "sigma_mad": sig}
    # a representative color if both members present
    if {"mag_cal_B", "mag_cal_V"} <= set(N.columns):
        n, med, sig = robust(((N.mag_cal_B - N.mag_cal_V) - (R.mag_cal_B - R.mag_cal_V)).values)
        out["color_BV"] = {"n": n, "median_dcolor": med, "sigma_mad": sig}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("new_csv")
    ap.add_argument("ref_csv")
    ap.add_argument("--tol-arcsec", type=float, default=1.0)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    res = compare(a.new_csv, a.ref_csv, a.tol_arcsec)
    print(f"NEW={res['n_new']} REF={res['n_ref']} matched(<{a.tol_arcsec}\")={res['n_matched']}")
    print("band  N      median_dMag  robust_sigma(MAD)")
    for b, s in res["bands"].items():
        print(f"  {b}   {s['n']:5d}   {s['median_dmag']:+.4f}     {s['sigma_mad']:.4f}")
    if "color_BV" in res:
        c = res["color_BV"]
        print(f"  (B-V) color: N={c['n']}  median={c['median_dcolor']:+.4f}  sigma={c['sigma_mad']:.4f}")
    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=2), encoding="utf-8")
        print(f"[json] {a.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
