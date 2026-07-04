"""PS1 cross-check: is the B faint drift a Gaia-BP reference artifact?

Tests (all with the same zp+ct-on-color fit machinery used for Gaia):
  A. APEX B/V instrumental vs PS1 (g, g-r)  -> drift vs magnitude
  B. Gaia-transformed ref_B/ref_V vs PS1    -> the two REFERENCES head-to-head
If A shows no drift while ref-vs-PS1 (B) drifts +0.03-0.04 faint-ward,
Gaia BP is convicted directly.
"""
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.vizier import Vizier
from scipy.spatial import cKDTree

mad = lambda v: 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v)))

cal = pd.read_csv(r"E:/observed_Analysis/NGC6811/pp/result/cmd_zeropoint/gaia_sdss_calibrator_by_ID.csv")
ra0, de0 = float(np.nanmedian(cal["ra_deg"])), float(np.nanmedian(cal["dec_deg"]))
sep_max = float(np.nanmax(SkyCoord(cal["ra_deg"], cal["dec_deg"], unit="deg").separation(
    SkyCoord(ra0, de0, unit="deg")).deg))
print(f"field center {ra0:.4f} {de0:+.4f}  radius {sep_max:.3f} deg  N_apex={len(cal)}")

v = Vizier(columns=["RAJ2000", "DEJ2000", "gmag", "e_gmag", "rmag", "e_rmag", "Qual"],
           column_filters={"gmag": "<20.5"}, row_limit=-1)
tbl = v.query_region(SkyCoord(ra0, de0, unit="deg"), radius=(sep_max + 0.02) * u.deg,
                     catalog="II/349/ps1")[0].to_pandas()
tbl = tbl.dropna(subset=["gmag", "rmag"])
print(f"PS1 rows: {len(tbl)}")

c_apex = SkyCoord(cal["ra_deg"].to_numpy(float), cal["dec_deg"].to_numpy(float), unit="deg")
c_ps1 = SkyCoord(tbl["RAJ2000"].to_numpy(float), tbl["DEJ2000"].to_numpy(float), unit="deg")
idx, sep, _ = c_apex.match_to_catalog_sky(c_ps1)
ok = sep.arcsec < 1.0
df = cal[ok].copy().reset_index(drop=True)
m = tbl.iloc[idx[ok]].reset_index(drop=True)
df["g"], df["r"] = m["gmag"].to_numpy(float), m["rmag"].to_numpy(float)
df["gr"] = df["g"] - df["r"]
print(f"matched <1 arcsec: {len(df)}")

xy = df[["x_pix", "y_pix"]].to_numpy(float)
d2, _ = cKDTree(xy).query(xy, k=2)
iso = d2[:, 1] > np.nanmedian(d2[:, 1])

# PS1 saturates ~g<13.5: bright anchor = 14.5-15.5 in ref mag, faint >= 17.5
def drift_table(label, apex_mag, refmag, snr, sel_extra=None):
    delta = apex_mag  # fit absorbs sign/zeropoint; use delta = X - model(color)
    good = np.isfinite(delta) & np.isfinite(df["gr"]) & np.isfinite(refmag) & (df["g"] > 13.5)
    if sel_extra is not None:
        good &= sel_extra
    fitm = good & (snr > 20) if snr is not None else good & refmag.between(14.0, 16.5)
    x = df["gr"].to_numpy(float); y = np.asarray(delta, float)
    mm = fitm.to_numpy() if hasattr(fitm, "to_numpy") else fitm
    for _ in range(4):
        c = np.polyfit(x[mm], y[mm], 1)
        r = y - np.polyval(c, x)
        s = mad(r[mm]); mm2 = mm & (np.abs(r - np.nanmedian(r[mm])) < 3 * max(s, 1e-6))
        if mm2.sum() == mm.sum():
            break
        mm = mm2
    r = y - np.polyval(c, x)
    rm = np.asarray(refmag, float)
    g = good.to_numpy() if hasattr(good, "to_numpy") else good
    out = []
    for lo, hi in ((14.5, 15.5), (16.0, 17.0), (17.0, 17.5), (17.5, 18.5)):
        k = g & (rm >= lo) & (rm < hi) & np.isfinite(r)
        if k.sum() < 12:
            out.append(f"{lo}-{hi}: N<12")
            continue
        out.append(f"{lo}-{hi}: {np.nanmedian(r[k]):+.4f} (N={k.sum()})")
    b = g & (rm >= 14.5) & (rm < 15.5) & np.isfinite(r)
    f = g & (rm >= 17.5) & np.isfinite(r)
    dr = np.nanmedian(r[f]) - np.nanmedian(r[b]) if (b.sum() >= 12 and f.sum() >= 12) else np.nan
    print(f"{label:34s} drift(14.5-15.5 vs >=17.5) = {dr:+.4f}")
    print("    " + " | ".join(out))
    return dr

print("\n=== A. APEX instrumental vs PS1 (measurement test) ===")
for band in ("B", "V"):
    am = -(df[f"mag_inst_{band}"].to_numpy(float))  # sign so brighter=larger irrelevant; fit absorbs
    drift_table(f"APEX {band} vs PS1(g,g-r) [all]", df[f"mag_inst_{band}"], df[f"ref_{band}"], df[f"snr_{band}"])
    drift_table(f"APEX {band} vs PS1(g,g-r) [iso]", df[f"mag_inst_{band}"], df[f"ref_{band}"], df[f"snr_{band}"], pd.Series(iso))

print("\n=== B. Gaia-transformed refs vs PS1 (reference-vs-reference) ===")
for band in ("B", "V"):
    drift_table(f"Gaia ref_{band} vs PS1(g,g-r)", df[f"ref_{band}"], df[f"ref_{band}"], None)
