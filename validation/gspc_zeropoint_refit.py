"""Refit the zero points against Gaia DR3 synthetic photometry, and see what moves.

APEX anchors its standard magnitudes on colour-colour transformations of Gaia
broadband photometry (Jordi+2010 for g/r/i, Pancino+2022 and Riello+2021 for
B/V/R). Those predict a star's magnitude in another system from three numbers,
which assumes stars of equal BP-RP have equal spectra. Metal-poor stars break
that assumption, and the break is measurable: comparing APEX's reference to the
Gaia Synthetic Photometry Catalogue (GSPC, Gaia Collaboration 2023 — synthetic
magnitudes integrated from the XP spectra and standardised against SDSS
Stripe 82 / PS1) gives a star-to-star scatter of 13 mmag in M67's g but
**52 mmag in M13's B**, the metal-poor globular in the blue band where line
blanketing bites hardest.

A constant offset would be harmless — the zero point absorbs it. Star-to-star
scatter is not absorbed; it lands in the CMD. So the question this answers is
whether swapping the reference for GSPC tightens the zero-point fit, and by how
much.

Method: step10's own fit path, with only the reference magnitudes swapped.
The mask (SNR >= 20, Gaia quality, finite colour/weights) and the fit
(inverse-variance weights, sigma-clipped, quadratic colour refinement) are the
production ones. Validated earlier: this reproduction returns exactly step10's
published N, zp and scatter for M67 g/r and M13 B when fed the current
reference.

GSPC covers G <= 17.65. That sounds restrictive but is not, for these data:
98-99.6 % of the stars that survive the SNR cut are already brighter than the
limit (M13's faintest calibrator is G = 16.3). Stars without GSPC are dropped,
never back-filled from the transformation — mixing two reference systems would
put a magnitude-dependent step at G = 17.65, which is the one systematic shape
photometry can least afford.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

REPO = Path(__file__).absolute().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from apex.gui.workflow.cmd.step10_zeropoint_calibration import (  # noqa: E402
    robust_weighted_polyfit,
)
from apex.utils.gaia_quality import gaia_quality_report  # noqa: E402

SNR_CUT, CLIP, ITERS, MIN_N = 20.0, 3.0, 5, 10
OUT = REPO / "validation" / "gspc_zeropoint_refit.json"

TARGETS = {
    "M67": dict(
        ra=132.825, dec=11.80,
        bands={"g": "g_sdss_mag", "r": "r_sdss_mag", "i": "i_sdss_mag"},
        colors={"g": ("g", "r"), "r": ("g", "r"), "i": ("r", "i")},
    ),
    "M13": dict(
        ra=250.421, dec=36.460,
        bands={"B": "b_jkc_mag", "V": "v_jkc_mag", "R": "r_jkc_mag"},
        colors={"B": ("B", "V"), "V": ("B", "V"), "R": ("V", "R")},
    ),
}


def fetch_gspc(ra: float, dec: float, cols: list[str]) -> pd.DataFrame:
    Gaia.ROW_LIMIT = -1
    sel = ", ".join(f"p.{c}" for c in cols)
    q = f"""SELECT s.ra, s.dec, {sel}
    FROM gaiadr3.gaia_source AS s
    LEFT JOIN gaiadr3.synthetic_photometry_gspc AS p ON s.source_id = p.source_id
    WHERE 1=CONTAINS(POINT('ICRS', s.ra, s.dec),
                     CIRCLE('ICRS', {ra}, {dec}, 0.4241))
      AND s.phot_g_mean_mag <= 18"""
    return Gaia.launch_job_async(q, dump_to_file=False).get_results().to_pandas()


def zp_for(cal: pd.DataFrame, band: str, color: tuple[str, str],
           ref_of: dict[str, str]) -> dict | None:
    ref = pd.to_numeric(cal[ref_of[band]], errors="coerce").to_numpy(float)
    inst = pd.to_numeric(cal[f"mag_inst_{band}"], errors="coerce").to_numpy(float)
    err = pd.to_numeric(cal[f"mag_inst_err_{band}"], errors="coerce").to_numpy(float)
    snr = pd.to_numeric(cal.get(f"snr_{band}"), errors="coerce").to_numpy(float)
    ca = pd.to_numeric(cal[ref_of[color[0]]], errors="coerce").to_numpy(float)
    cb = pd.to_numeric(cal[ref_of[color[1]]], errors="coerce").to_numpy(float)

    delta, x = ref - inst, ca - cb
    w = np.where(np.isfinite(err) & (err > 0), 1.0 / err**2, np.nan)
    qual, _ = gaia_quality_report(cal, cstar_nsigma=None)
    m = (np.isfinite(delta) & np.isfinite(x) & np.isfinite(w)
         & np.isfinite(snr) & (snr >= SNR_CUT) & qual)
    if int(m.sum()) < MIN_N:
        return None
    c, n, s = robust_weighted_polyfit(x[m], delta[m], w=w[m], degree=2,
                                      clip_sigma=CLIP, iters=ITERS, min_n=MIN_N)
    if c is None:
        return None
    return {"n_in": int(m.sum()), "n_inliers": int(n), "zp": float(c[2]),
            "ct": float(c[1]), "ct2": float(c[0]),
            "scatter_mmag": float(s) * 1000}


results = {}
for name, spec in TARGETS.items():
    print(f"=== {name} ===")
    cal = pd.read_csv(Path(rf"E:\APEX_validation\phase3\{name}\result"
                           r"\cmd_zeropoint\gaia_sdss_calibrator_by_ID.csv"))
    g = fetch_gspc(spec["ra"], spec["dec"], list(spec["bands"].values()))
    idx, sep, _ = SkyCoord(cal["ra_deg"], cal["dec_deg"], unit="deg") \
        .match_to_catalog_sky(SkyCoord(g["ra"], g["dec"], unit="deg"))
    ok = sep.arcsec < 1.0
    for band, gcol in spec["bands"].items():
        vals = g[gcol].to_numpy()[idx]
        cal[f"gspc_{band}"] = np.where(ok, vals, np.nan)

    ref_current = {b: f"ref_{b}" for b in spec["bands"]}
    ref_gspc = {b: f"gspc_{b}" for b in spec["bands"]}

    rows = []
    for band, color in spec["colors"].items():
        cur = zp_for(cal, band, color, ref_current)
        new = zp_for(cal, band, color, ref_gspc)
        if not (cur and new):
            continue
        rows.append({"band": band, "current": cur, "gspc": new})
        print(f"  {band}: N {cur['n_inliers']:4d} -> {new['n_inliers']:4d}"
              f"   zp {cur['zp']:+.4f} -> {new['zp']:+.4f}"
              f"   ct {cur['ct']:+.4f} -> {new['ct']:+.4f}"
              f"   scatter {cur['scatter_mmag']:5.1f} -> {new['scatter_mmag']:5.1f} mmag"
              f"  ({100*(new['scatter_mmag']/cur['scatter_mmag']-1):+.0f} %)")
    results[name] = rows

OUT.write_text(json.dumps(results, indent=1), encoding="utf-8")
print(f"\nsaved -> {OUT}")
