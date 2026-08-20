"""Can a single field tell a real curvature from a fitted one?

The transfer test needs five fields. A run has one. So the question is whether
anything measurable *inside* one fit separates the B-band curvature (which
transfers between clusters) from the R-band one (which makes other clusters
worse by up to 17 mmag).

Candidates, all computable from the fit itself:
  - |ct2| / sigma(ct2)      how well the curvature is constrained
  - colour baseline         a quadratic over a narrow span chases noise
  - leverage balance        how many calibrators sit at the ends of the span
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, os.getcwd())

import numpy as np
import pandas as pd

from apex.utils.gaia_quality import gaia_quality_report

SNR_CUT = 20.0
ROOT = Path(r"E:\APEX_validation\reprocess")
CLUSTERS = ["M13", "M3", "M5", "M67", "NGC6811"]
MAD = 1.4826

# From d012_transfer.py: does this band's curvature help other fields?
VERDICT = {"B": "진짜", "G": "진짜(작음)", "V": "가짜", "R": "가짜", "I": "가짜"}


def fits():
    for cluster in CLUSTERS:
        base = ROOT / cluster / "result" / "cmd_zeropoint"
        cal, coeff = (base / "gaia_sdss_calibrator_by_ID.csv",
                      base / "zp_fit_coefficients.csv")
        if not cal.exists() or not coeff.exists():
            continue
        table = pd.read_csv(cal)
        m_qual = np.asarray(gaia_quality_report(table, cstar_nsigma=None)[0], bool)
        for _, row in pd.read_csv(coeff).iterrows():
            filt, pair = str(row["filter"]), str(row["color_col"])
            dcol, ccol = f"delta_{filt}", f"color_{pair}"
            if dcol not in table.columns or ccol not in table.columns:
                continue
            x = pd.to_numeric(table[ccol], errors="coerce").to_numpy(float)
            y = pd.to_numeric(table[dcol], errors="coerce").to_numpy(float)
            err = pd.to_numeric(table.get(f"mag_inst_err_{filt}", 0.01),
                                errors="coerce").to_numpy(float)
            with np.errstate(divide="ignore", invalid="ignore"):
                w = 1.0 / np.square(err)
            snr = pd.to_numeric(table.get(f"snr_{filt}", np.inf),
                                errors="coerce").to_numpy(float)
            keep = (np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
                    & np.isfinite(snr) & (snr >= SNR_CUT) & m_qual)
            if keep.sum() < 50:
                continue
            yield cluster, filt.upper(), pair, x[keep], y[keep], w[keep]


print(f"{'성단':<9} {'필터':<3} {'색':<5} {'ct2':>9} {'σ(ct2)':>8} {'ct2/σ':>7} "
      f"{'색범위':>7} {'양끝%':>7}  {'교차검증':<10}")
print("-" * 84)

rows = []
for cluster, band, pair, x, y, w in fits():
    c = np.polyfit(x, y, 2, w=np.sqrt(w))
    r = y - np.polyval(c, x)
    med = np.median(r)
    keep = np.abs(r - med) <= 3.0 * max(MAD * np.median(np.abs(r - med)), 1e-6)
    xc, yc, wc = x[keep], y[keep], w[keep]

    # Weighted least squares covariance, scaled by the fit's own chi2/dof so
    # sigma reflects the actual residual spread rather than the quoted errors.
    A = np.column_stack([xc * xc, xc, np.ones_like(xc)]) * np.sqrt(wc)[:, None]
    coef, *_ = np.linalg.lstsq(A, yc * np.sqrt(wc), rcond=None)
    resid = yc * np.sqrt(wc) - A @ coef
    dof = max(len(xc) - 3, 1)
    cov = np.linalg.inv(A.T @ A) * float(resid @ resid) / dof
    ct2, sig_ct2 = float(coef[0]), float(np.sqrt(cov[0, 0]))

    span = float(np.percentile(xc, 97.5) - np.percentile(xc, 2.5))
    lo, hi = np.percentile(xc, [20, 80])
    ends = float(((xc <= lo) | (xc >= hi)).mean() * 100.0)

    print(f"{cluster:<9} {band:<3} {pair.replace('_','-'):<5} {ct2:+9.4f} "
          f"{sig_ct2:8.4f} {abs(ct2)/max(sig_ct2,1e-9):7.1f} {span:7.3f} "
          f"{ends:6.1f}%  {VERDICT.get(band,'?'):<10}")
    rows.append({"cluster": cluster, "band": band, "pair": pair, "ct2": ct2,
                 "sigma": sig_ct2, "snr": abs(ct2) / max(sig_ct2, 1e-9),
                 "span": span, "ends": ends, "verdict": VERDICT.get(band, "?")})

df = pd.DataFrame(rows)
print()
for name, col, unit in (("ct2/σ (구속 정도)", "snr", ""),
                        ("색 범위", "span", " mag"),
                        ("|ct2|", "ct2", "")):
    real = df[df["verdict"].str.startswith("진짜")][col].abs()
    fake = df[df["verdict"] == "가짜"][col].abs()
    print(f"  {name:<18} 진짜 {real.min():.3f}~{real.max():.3f}{unit}"
          f"   ·   가짜 {fake.min():.3f}~{fake.max():.3f}{unit}"
          f"   {'← 갈린다' if real.min() > fake.max() or real.max() < fake.min() else '← 안 갈린다'}")
df.to_csv("validation/color_term_quadratic_criterion.csv", index=False)
