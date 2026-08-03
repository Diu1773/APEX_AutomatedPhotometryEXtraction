"""Re-anchor the M67 U-band zero-point on external Johnson U standards.

Why: the U-B degeneracy test (REPORT_UB_DEGENERACY.md) failed because the U
reference came from the Gaia G-U approx polynomial (sigma = 0.200 mag).  M67 is
a classic photometric standard field — Montgomery, Marschall & Janes 1993
(AJ 106, 181; VizieR J/AJ/106/181, table3) provides UBVRI for 1,456 stars.

What this does:
  1. fetch MMJ93 table3 (cached to CSV next to the result dir),
  2. cross-match with APEX Step-10 wide CMD table (ra/dec, 1.5"),
  3. measure dU = U_MMJ - mag_std_U (median offset + colour-term check),
     with dB / dV as controls (their ZPs were already good, ~0.025/0.030),
  4. write a corrected copy of the wide table under result_ustd/cmd_zeropoint/
     with mag_std_U shifted by the measured offset (APEX instrumental
     photometry untouched — only the zero-point anchor changes).

Run:
    .venv-deploy/Scripts/python.exe -X utf8 \
        validation/psf_crossinstrument/recalibrate_u_mmj93.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(r"E:/APEX_validation/psf_crossinstrument/m67_ubv")
EXT_DIR = BASE / "external_ref"
MATCH_ARCSEC = 1.5

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--source-result", default="result",
                 help="result dir (under m67_ubv/) whose wide table to anchor, "
                      "e.g. result (aperture) or result_psf (PSF)")
_ap.add_argument("--suffix", default="",
                 help="suffix for output dirs: result{suffix}_ustd / _uonly")
ARGS = _ap.parse_args()

WIDE = (BASE / ARGS.source_result / "cmd_zeropoint"
        / "median_by_ID_filter_wide_cmd.csv")
OUT_RESULT = BASE / f"result{ARGS.suffix}_ustd"
OUT_UONLY = BASE / f"result{ARGS.suffix}_uonly"
RES_CSV = EXT_DIR / f"mmj93_match_residuals{ARGS.suffix or ''}.csv"


def fetch_mmj93() -> pd.DataFrame:
    cache = EXT_DIR / "mmj93_table3.csv"
    if cache.exists():
        return pd.read_csv(cache)
    from astroquery.vizier import Vizier

    v = Vizier(columns=["**"], row_limit=-1)
    tab = v.get_catalogs("J/AJ/106/181")[0]  # table3
    df = tab.to_pandas()
    EXT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache, index=False, encoding="utf-8")
    return df


def main() -> int:
    mmj = fetch_mmj93()
    # Johnson magnitudes from V + colours; require U-B present with error
    for c in ("Vmag", "B-V", "U-B"):
        mmj[c] = pd.to_numeric(mmj[c], errors="coerce")
    mmj["Bmag"] = mmj["Vmag"] + mmj["B-V"]
    mmj["Umag"] = mmj["Bmag"] + mmj["U-B"]
    mmj = mmj.dropna(subset=["Umag", "_RA.icrs", "_DE.icrs"]).copy()
    print(f"MMJ93 table3: {len(mmj)} stars with U")

    wide = pd.read_csv(WIDE)
    ok = wide["mag_std_U"].notna() & wide["ra_deg"].notna()
    w = wide[ok].copy()
    print(f"APEX wide table: {len(wide)} rows, {len(w)} with U + position")

    from astropy.coordinates import SkyCoord
    import astropy.units as u

    c_apex = SkyCoord(w["ra_deg"].values * u.deg, w["dec_deg"].values * u.deg)
    ra_num = pd.to_numeric(mmj["_RA.icrs"], errors="coerce")
    if ra_num.notna().all():
        c_mmj = SkyCoord(ra_num.values * u.deg,
                         pd.to_numeric(mmj["_DE.icrs"]).values * u.deg)
    else:  # sexagesimal strings from VizieR
        c_mmj = SkyCoord(mmj["_RA.icrs"].astype(str).values,
                         mmj["_DE.icrs"].astype(str).values,
                         unit=(u.hourangle, u.deg))
    idx, sep, _ = c_apex.match_to_catalog_sky(c_mmj)
    m = sep.arcsec < MATCH_ARCSEC
    print(f"matches < {MATCH_ARCSEC}\": {m.sum()}  "
          f"(sep median {np.median(sep.arcsec[m]):.2f}\")")

    mm = mmj.iloc[idx[m]].reset_index(drop=True)
    ww = w[m].reset_index(drop=True)

    res = pd.DataFrame({
        "ID": ww["ID"],
        "sep_arcsec": sep.arcsec[m],
        "BV_mmj": mm["B-V"],
        "dU": mm["Umag"].values - ww["mag_std_U"].values,
        "dB": mm["Bmag"].values - ww["mag_std_B"].values,
        "dV": mm["Vmag"].values - ww["mag_std_V"].values,
        "U_mmj": mm["Umag"].values,
        "eUB_mmj": pd.to_numeric(mm["e_U-B"], errors="coerce").values,
        "snr_U": ww["snr_U"].values,
    })

    def stats(col: str) -> tuple[float, float, int]:
        v = res[col].dropna().values
        med = float(np.median(v))
        mad = float(1.4826 * np.median(np.abs(v - med)))
        return med, mad, len(v)

    print("\n  band   median offset   robust sigma   N   (MMJ - APEX_std)")
    for col in ("dU", "dB", "dV"):
        med, mad, n = stats(col)
        print(f"  {col}    {med:+.4f}         {mad:.4f}       {n}")

    # colour dependence of dU (bright, well-measured stars only)
    good = res.dropna(subset=["dU", "BV_mmj"])
    good = good[good["snr_U"] > 20]
    if len(good) >= 10:
        coef = np.polyfit(good["BV_mmj"], good["dU"], 1)
        pred = np.polyval(coef, good["BV_mmj"])
        rms = float(np.std(good["dU"] - pred))
        print(f"\ndU vs (B-V): slope {coef[0]:+.4f} mag/mag, "
              f"intercept {coef[1]:+.4f}, rms {rms:.4f}  (N={len(good)}, snr_U>20)")

    dU_med, dU_mad, nU = stats("dU")
    dB_med, _, _ = stats("dB")
    dV_med, _, _ = stats("dV")

    # Two corrected copies (APEX instrumental photometry untouched — only the
    # zero-point anchor changes):
    #   result_ustd  — U, B, V all re-anchored to MMJ93 (one homogeneous
    #                  Johnson system; Sandquist 2004 confirms the B direction)
    #   result_uonly — only U shifted (ablation: is B's 0.03-0.05 also needed?)
    variants = {
        OUT_RESULT: {"U": dU_med, "B": dB_med, "V": dV_med},
        OUT_UONLY: {"U": dU_med},
    }
    for out_dir, shifts in variants.items():
        out_zp = out_dir / "cmd_zeropoint"
        out_zp.mkdir(parents=True, exist_ok=True)
        wide_c = wide.copy()
        for band, dz in shifts.items():
            wide_c[f"mag_std_{band}"] = wide_c[f"mag_std_{band}"] + dz
        out_csv = out_zp / "median_by_ID_filter_wide_cmd.csv"
        wide_c.to_csv(out_csv, index=False, encoding="utf-8")
        print(f"{out_dir.name}: shifts {shifts} -> {out_csv}")

    res.to_csv(RES_CSV, index=False, encoding="utf-8")
    print(f"\ndU = {dU_med:+.4f} (robust sigma {dU_mad:.4f}, N={nU})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
