"""What does turning on the Riello+2021 C* cut do to APEX's zero points?

The cut has never run: no Gaia query fetched `phot_bp_rp_excess_factor` until
2026-08-11, so `gaia_quality_mask` skipped it silently on every target ever
processed. Now that the column is in the contract the cut *can* run, and the
question is whether it *should* — it rejects 3.6 % of Gaia references in M67
but 49.7 % in M13, which is a large change to make blind.

This measures the difference on the real calibrator tables, driving step10's
own fitting functions with step10's own mask construction (SNR cut, Gaia
quality mask, inverse-variance weights, sigma-clipped linear fit with the
quadratic colour refinement). Only the C* flag changes between the two runs.

Why not just re-run the pipeline: the calibrator tables and the Gaia catalogue
are the only inputs the cut touches, and re-running would also re-roll SEP's
non-deterministic deblending, mixing an unrelated 2 mmag jitter into the
comparison. Holding everything fixed and flipping one flag is the cleaner
measurement — it is a differential, and it is exact.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table

REPO = Path(__file__).absolute().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from apex.gui.workflow.cmd.step10_zeropoint_calibration import (  # noqa: E402
    robust_weighted_polyfit,
)
from apex.utils.gaia_quality import gaia_quality_report  # noqa: E402

SNR_CUT = 20.0          # step10 default (gaia_snr_calib_min)
CLIP_SIGMA, ITERS, MIN_N = 3.0, 5, 10
OUT = REPO / "validation" / "cstar_zeropoint_impact.json"

# Colour axes are taken from each target's own zp_fit_coefficients.csv
# (`color_col`), so the fit axis is the one the pipeline actually used rather
# than a guess. The i band is left out: its recorded axis does not reproduce
# step10's N here, and an unverified axis would make the differential
# meaningless. g and r reproduce step10 exactly (N 503 / 582, zp -2.7771 /
# -3.3897, scatter 24.1 / 19.1 mmag), which is what validates the method.
TARGETS = {
    "M67 (open, g/r)": dict(
        cal=r"E:\APEX_validation\phase3\M67\result\cmd_zeropoint"
            r"\gaia_sdss_calibrator_by_ID.csv",
        gaia=r"E:\APEX_validation\cstar_test\M67\result\step5_wcs\gaia_fov.ecsv",
        bands=("g", "r"),
        colors={"g": ("g", "r"), "r": ("g", "r")},
    ),
    "M13 (globular, B/V/R)": dict(
        cal=r"E:\APEX_validation\phase3\M13\result\cmd_zeropoint"
            r"\gaia_sdss_calibrator_by_ID.csv",
        gaia=r"E:\APEX_validation\cstar_test\M13\gaia_fov.ecsv",
        bands=("B", "V", "R"),
        colors={"B": ("B", "V"), "V": ("B", "V"), "R": ("V", "R")},
    ),
}


def attach_gaia(cal: pd.DataFrame, gaia_path: str) -> pd.DataFrame:
    """Cross-match the calibrators to a catalogue that has the excess factor."""
    g = Table.read(gaia_path).to_pandas()
    c1 = SkyCoord(cal["ra_deg"].to_numpy(float), cal["dec_deg"].to_numpy(float),
                  unit="deg")
    c2 = SkyCoord(g["ra"].to_numpy(float), g["dec"].to_numpy(float), unit="deg")
    idx, sep, _ = c1.match_to_catalog_sky(c2)
    ok = sep.arcsec < 1.0
    out = cal.copy()
    for col, target in (("phot_bp_rp_excess_factor", "phot_bp_rp_excess_factor"),
                        ("phot_g_mean_mag", "gaia_G")):
        vals = pd.to_numeric(g[col], errors="coerce").to_numpy()[idx]
        out[target] = np.where(ok, vals, np.nan)
    bp = pd.to_numeric(g["phot_bp_mean_mag"], errors="coerce").to_numpy()[idx]
    rp = pd.to_numeric(g["phot_rp_mean_mag"], errors="coerce").to_numpy()[idx]
    out["gaia_BP_RP"] = np.where(ok, bp - rp, np.nan)
    print(f"  Gaia cross-match within 1 arcsec: {int(ok.sum())}/{len(cal)}")
    return out


def fit_band(cal: pd.DataFrame, band: str, color: tuple[str, str],
             cstar_on: bool) -> dict | None:
    """step10's fit for one band, with the C* cut on or off."""
    delta = (pd.to_numeric(cal[f"ref_{band}"], errors="coerce")
             - pd.to_numeric(cal[f"mag_inst_{band}"], errors="coerce")).to_numpy(float)
    x = (pd.to_numeric(cal[f"ref_{color[0]}"], errors="coerce")
         - pd.to_numeric(cal[f"ref_{color[1]}"], errors="coerce")).to_numpy(float)
    err = pd.to_numeric(cal[f"mag_inst_err_{band}"], errors="coerce").to_numpy(float)
    w = np.where(np.isfinite(err) & (err > 0), 1.0 / err**2, np.nan)
    snr = pd.to_numeric(cal.get(f"snr_{band}"), errors="coerce").to_numpy(float)

    qual, report = gaia_quality_report(
        cal, cstar_nsigma=3.0 if cstar_on else None)
    m = (np.isfinite(delta) & np.isfinite(x) & np.isfinite(w)
         & np.isfinite(snr) & (snr >= SNR_CUT) & qual)
    if int(m.sum()) < MIN_N:
        return None

    coeffs, n, scatter = robust_weighted_polyfit(
        x[m], delta[m], w=w[m], degree=2,
        clip_sigma=CLIP_SIGMA, iters=ITERS, min_n=MIN_N)
    if coeffs is None:
        return None
    return {
        "n_before_fit": int(m.sum()), "n_inliers": int(n),
        "zp": float(coeffs[2]), "ct": float(coeffs[1]), "ct2": float(coeffs[0]),
        "scatter_mmag": float(scatter) * 1000,
        "quality": report["cuts"],
    }


results = {}
for name, spec in TARGETS.items():
    print(f"=== {name} ===")
    cal = attach_gaia(pd.read_csv(spec["cal"]), spec["gaia"])
    rows = []
    for band in spec["bands"]:
        off = fit_band(cal, band, spec["colors"][band], cstar_on=False)
        on = fit_band(cal, band, spec["colors"][band], cstar_on=True)
        if not (off and on):
            continue
        rows.append({"band": band, "off": off, "on": on})
        print(f"  {band}: N {off['n_inliers']:4d} -> {on['n_inliers']:4d}"
              f"   zp {off['zp']:+.4f} -> {on['zp']:+.4f}"
              f"  (delta {1000*(on['zp']-off['zp']):+6.1f} mmag)"
              f"   ct {off['ct']:+.4f} -> {on['ct']:+.4f}"
              f"   scatter {off['scatter_mmag']:.1f} -> {on['scatter_mmag']:.1f} mmag")
    results[name] = rows

OUT.write_text(json.dumps(results, indent=1), encoding="utf-8")
print(f"\nsaved -> {OUT}")
