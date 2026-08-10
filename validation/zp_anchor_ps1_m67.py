"""Where is APEX's standard system anchored, and by how much does it differ
from independent photometry?

Step 10 calibrates M67 g/r/i against Jordi+2010 — a *transformation* from Gaia
photometry, not measured SDSS/PS1 standards. A transformation carries its own
zeropoint and colour-term systematics. This script measures that anchor
difference directly: match the calibrated catalogue to PS1 DR2 (independent,
measured photometry in nearly the same bands) and fit

    mag_std_X - ps1_X = zp + ct * (g-r)_PS1        (iterative 3-sigma clip)

The fitted colour term absorbs the small PS1<->SDSS system difference without
importing literature transform coefficients (whose own errors would then be
part of the audit). The number that matters is the offset at the sample's
median colour: that is the anchor systematic a user of APEX magnitudes
inherits.

Context: the distance-modulus rail investigation asked whether the zeropoints
are ~0.06 mag bright. The turn-off comparison said no; this is the direct
photometric test of the same question.

Fit range 14.0 < PS1 < 16.5: PS1 saturates around 13.5 and the faint end
mixes in the B-filter-style reference drifts already documented elsewhere.
"""

import json

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.vizier import Vizier
from pathlib import Path

WIDE = Path(r"E:\APEX_validation\reprocess\M67\result\cmd_zeropoint"
            r"\median_by_ID_filter_wide_cmd.csv")
CAL = Path(r"E:\APEX_validation\reprocess\M67\result\cmd_zeropoint"
           r"\gaia_sdss_calibrator_by_ID.csv")
OUT = Path("validation/paper/data/zp_anchor_ps1_m67.json")

mad = lambda v: 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v)))

wide = pd.read_csv(WIDE)
cal = pd.read_csv(CAL)

ra0 = float(np.nanmedian(wide["ra_deg"]))
de0 = float(np.nanmedian(wide["dec_deg"]))
coords = SkyCoord(wide["ra_deg"].to_numpy(float),
                  wide["dec_deg"].to_numpy(float), unit="deg")
radius = float(np.nanmax(coords.separation(SkyCoord(ra0, de0, unit="deg")).deg))
print(f"M67 field {ra0:.4f} {de0:+.4f}  radius {radius:.3f} deg  "
      f"catalog rows {len(wide)}")

v = Vizier(columns=["RAJ2000", "DEJ2000", "gmag", "e_gmag", "rmag", "e_rmag",
                    "imag", "e_imag", "Qual"],
           column_filters={"gmag": "<20.5"}, row_limit=-1)
ps1 = v.query_region(SkyCoord(ra0, de0, unit="deg"),
                     radius=(radius + 0.02) * u.deg,
                     catalog="II/349/ps1")[0].to_pandas()
ps1 = ps1.dropna(subset=["gmag", "rmag", "imag"])
print(f"PS1 rows: {len(ps1)}")

c_ps1 = SkyCoord(ps1["RAJ2000"].to_numpy(float),
                 ps1["DEJ2000"].to_numpy(float), unit="deg")
idx, sep, _ = coords.match_to_catalog_sky(c_ps1)
ok = sep.arcsec < 1.0
df = wide[ok].reset_index(drop=True)
m = ps1.iloc[idx[ok]].reset_index(drop=True)
for band in "gri":
    df[f"ps1_{band}"] = m[f"{band}mag"].to_numpy(float)
df["ps1_gr"] = df["ps1_g"] - df["ps1_r"]
print(f"matched < 1 arcsec: {len(df)}")

# The calibrator table carries the Jordi-transformed reference (ref_g/r/i) for
# a subset of the same IDs — merge it in so the two references can be compared
# head-to-head against PS1 on identical stars.
df = df.merge(cal[["ID", "ref_g", "ref_r", "ref_i"]], on="ID", how="left")


def anchor_offset(label, series, band):
    """Fitted offset at the median colour, colour term absorbed."""
    y = np.asarray(series, float) - df[f"ps1_{band}"].to_numpy(float)
    x = df["ps1_gr"].to_numpy(float)
    r_ps1 = df[f"ps1_{band}"].to_numpy(float)
    good = np.isfinite(y) & np.isfinite(x) & (r_ps1 > 14.0) & (r_ps1 < 16.5)
    if good.sum() < 30:
        print(f"{label:26s} N<{30}")
        return None
    mm = good.copy()
    for _ in range(4):
        c = np.polyfit(x[mm], y[mm], 1)
        resid = y - np.polyval(c, x)
        s = mad(resid[mm])
        nxt = mm & (np.abs(resid - np.nanmedian(resid[mm])) < 3 * max(s, 1e-6))
        if nxt.sum() == mm.sum():
            break
        mm = nxt
    x0 = float(np.nanmedian(x[mm]))
    offset = float(np.polyval(c, x0))
    scatter = float(mad((y - np.polyval(c, x))[mm]))
    print(f"{label:26s} offset@gr={x0:.2f}: {offset*1000:+7.1f} mmag   "
          f"ct={c[0]:+.4f}   MAD={scatter*1000:5.1f} mmag   N={int(mm.sum())}")
    return {"label": label, "band": band, "offset_mmag": offset * 1000,
            "color_term": float(c[0]), "mad_mmag": scatter * 1000,
            "n": int(mm.sum()), "median_gr": x0}


print("\n=== APEX calibrated magnitudes vs PS1 (the anchor test) ===")
results = []
for band in "gri":
    r = anchor_offset(f"mag_std_{band} - ps1_{band}", df[f"mag_std_{band}"], band)
    if r:
        results.append(r)

print("\n=== Jordi+2010 reference vs PS1 (the two references head-to-head) ===")
for band in "gri":
    r = anchor_offset(f"ref_{band}   - ps1_{band}", df[f"ref_{band}"], band)
    if r:
        results.append(r)

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps({
    "field": {"ra": ra0, "dec": de0, "radius_deg": radius},
    "n_matched": int(len(df)),
    "fit_range_ps1": [14.0, 16.5],
    "results": results,
}, indent=1), encoding="utf-8")
print(f"\nsaved -> {OUT}")
