"""APEX PDM vs PyAstronomy pyPDM on the identical light curve and trial grid.

The engine scorecard's A axis: APEX carries its own Stellingwerf (1978) PDM
because astropy.timeseries has none. Absence justifies existence, but not
correctness — that needs an independent implementation on the same input.

Input: the committed YZ Boo run (E:\APEX_validation\reprocess\YZBoo_2n),
series `diff_mag_raw` over BJD_TDB (364 points), and the trial periods of the
committed APEX periodogram. APEX's theta is read from that CSV — nothing is
re-run on the APEX side, so this compares the *shipped* output.

pyPDM caveat: its Scanner only generates uniform grids, but PyPDM.pdmEquiBinCover
accepts arbitrary period arrays via the frequency interface? No — simplest
robust route: call PyPDM with a dense scanner over the same range, then
interpolate both theta curves onto the APEX grid for correlation, and compare
the argmin directly. Bin count matched to APEX (10); pyPDM's covers=3 variant
is also reported since it is that package's default.
"""

import json

import numpy as np
from astropy.table import Table
from PyAstronomy.pyTiming import pyPDM

LC = (r"E:\APEX_validation\reprocess\YZBoo_2n\result\lc_lightcurve"
      r"\lightcurve_ID153_raw.csv")
PG = (r"E:\APEX_validation\reprocess\YZBoo_2n\result\lc_period"
      r"\periodogram_all_raw_pdm_ID153.csv")
OUT = r"C:\ast_v\pdm_compare.json"

lc = Table.read(LC, format="ascii.csv")
t = np.asarray(lc["BJD_TDB"], float)
y = np.asarray(lc["diff_mag_raw"], float)
filt = np.asarray(lc["filter"]).astype(str)
# The "all" series mixes g/r/i differential magnitudes whose per-filter
# offsets (~0.8 mag) dwarf the 0.4 mag pulsation. Fed raw, any PDM sees no
# signal (theta ~ 0.92 everywhere — measured on the first attempt). APEX's
# combined periodogram removes the per-filter level first; do the same here.
for b in set(filt):
    m = filt == b
    y[m] = y[m] - np.nanmedian(y[m])
ok = np.isfinite(t) & np.isfinite(y)
t, y = t[ok], y[ok]
print(f"light curve points: {len(t)}  baseline {t.max()-t.min():.2f} d")

pg = Table.read(PG, format="ascii.csv")
periods = np.asarray(pg["period"], float)
theta_apex = np.asarray(pg["theta"], float)
print(f"APEX periodogram: {len(periods)} trials, "
      f"P {periods.min():.4f}-{periods.max():.4f} d")

# pyPDM is a pure-Python loop over trials, so scanning APEX's full 0.01-10 d
# grid at matching resolution would take hours (about 1e6 trials). The audit
# question — do two implementations agree — is answered in the window that
# contains the science: around the YZ Boo period and its aliases.
WIN = (0.05, 0.30)
in_win = (periods >= WIN[0]) & (periods <= WIN[1])
periods, theta_apex = periods[in_win], theta_apex[in_win]
print(f"comparison window {WIN[0]}-{WIN[1]} d: {len(periods)} APEX trials")
scanner = pyPDM.Scanner(minVal=WIN[0], maxVal=WIN[1], dVal=2e-5, mode="period")
pdm = pyPDM.PyPDM(t, y)

results = {}
for label, method in (("equi_10bins", lambda s: pdm.pdmEquiBin(10, s)),
                      ("cover_10x3", lambda s: pdm.pdmEquiBinCover(10, 3, s))):
    p2, theta2 = method(scanner)
    theta_on_apex = np.interp(periods, p2, theta2)
    best_py = float(p2[np.argmin(theta2)])
    best_apex = float(periods[np.argmin(theta_apex)])
    corr = float(np.corrcoef(theta_apex, theta_on_apex)[0, 1])
    results[label] = {
        "best_period_pypdm": best_py,
        "best_period_apex": best_apex,
        "delta_period_d": best_py - best_apex,
        "theta_correlation_on_apex_grid": corr,
        "theta_min_pypdm": float(np.min(theta2)),
        "theta_min_apex": float(np.min(theta_apex)),
    }
    print(f"[{label}] best P: pyPDM {best_py:.6f}  APEX {best_apex:.6f}  "
          f"delta {best_py-best_apex:+.6f} d   theta corr {corr:.4f}   "
          f"theta_min {np.min(theta2):.4f} vs {np.min(theta_apex):.4f}")

results["literature_period_d"] = 0.104092
json.dump(results, open(OUT, "w"), indent=1)
print(f"\nliterature YZ Boo P = 0.104092 d")
print(f"saved -> {OUT}")
