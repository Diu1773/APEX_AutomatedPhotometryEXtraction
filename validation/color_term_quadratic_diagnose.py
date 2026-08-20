"""D-012: why does the quadratic colour term switch on?

TRACK framed this as "only in R band". Across five clusters it is on in 11 of
15 fits and in every band, so the question is not which band — it is what the
acceptance rule can actually distinguish.

This re-runs the code's own two fits on the code's own saved calibrator tables
and reports, for each, whether the quadratic buys anything a free parameter
would not buy by chance.
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, os.getcwd())

import numpy as np
import pandas as pd
from scipy import stats

from apex.analysis.cmd.zeropoint_runner import (
    ZeropointCalibrationRunner, robust_weighted_polyfit,
)

# `_robust_linfit` is a method but touches no state, so call it unbound rather
# than reimplementing it — a reimplementation would be a second copy to drift.
_robust_linfit = ZeropointCalibrationRunner._robust_linfit.__get__(object())

from apex.utils.gaia_quality import gaia_quality_report

SNR_CUT = 20.0          # gaia_snr_calib_min default
ROOT = Path(r"E:\APEX_validation\reprocess")
CLUSTERS = ["M13", "M3", "M5", "M67", "NGC6811"]

def _recorded_ct2(row) -> float:
    """What that run fitted, whichever column the run happened to write.

    Runs from 2026-08-21 record the fit in `ct2_fitted` and apply it only when
    configured, so `ct2` is 0 on a default run. Reading `ct2` alone would make
    this script's own evidence vanish as workspaces are reprocessed — the
    spread across fields is the finding, and it has to stay readable.
    """
    for column in ("ct2_fitted", "ct2"):
        if column in row.index:
            value = float(row[column])
            if value == value:               # not NaN
                return value
    return 0.0


print(f"{'성단':<9} {'필터':<3} {'N(1차)':>7} {'N(2차)':>7} "
      f"{'산포1차':>9} {'산포2차':>9} {'ct1차':>8} {'ct2차':>8} {'ct2':>8} "
      f"{'F':>7} {'p':>9}  채택?  (기록 N 대조)")
print("-" * 108)

rows = []
for cluster in CLUSTERS:
    cal = ROOT / cluster / "result" / "cmd_zeropoint" / "gaia_sdss_calibrator_by_ID.csv"
    coeff_path = ROOT / cluster / "result" / "cmd_zeropoint" / "zp_fit_coefficients.csv"
    if not cal.exists() or not coeff_path.exists():
        continue
    table = pd.read_csv(cal)
    coeffs = pd.read_csv(coeff_path)

    for _, crow in coeffs.iterrows():
        filt = str(crow["filter"])
        pair = str(crow["color_col"])
        dcol, ccol = f"delta_{filt}", f"color_{pair}"
        if dcol not in table.columns or ccol not in table.columns:
            continue
        ecol = f"mag_inst_err_{filt}"
        x = pd.to_numeric(table[ccol], errors="coerce").to_numpy(float)
        y = pd.to_numeric(table[dcol], errors="coerce").to_numpy(float)
        if ecol in table.columns:
            err = pd.to_numeric(table[ecol], errors="coerce").to_numpy(float)
            with np.errstate(divide="ignore", invalid="ignore"):
                w = 1.0 / np.square(err)
            w[~np.isfinite(w)] = np.nan
        else:
            w = np.ones_like(y)
        # The code's own fit sample, not all rows: SNR cut + Gaia quality cut.
        # Fitting every row gave N=1162 where the run recorded 433 — a different
        # sample answers a different question.
        snr_col = f"snr_{filt}"
        if snr_col in table.columns:
            sv = pd.to_numeric(table[snr_col], errors="coerce").to_numpy(float)
            m_snr = np.isfinite(sv) & (sv >= SNR_CUT)
        else:
            m_snr = np.ones(len(table), dtype=bool)
        m_qual, _report = gaia_quality_report(table, cstar_nsigma=None)
        m_qual = np.asarray(m_qual, dtype=bool)

        keep = (np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
                & m_snr & m_qual)
        x, y, w = x[keep], y[keep], w[keep]
        if x.size < 30:
            continue

        zp1, ct1, n1, s1 = _robust_linfit(x, y, w=w, slope_absmax=0.8)
        q, n2, s2 = robust_weighted_polyfit(x, y, w=w, degree=2)
        if q is None:
            continue
        ct2, ctq, zpq = float(q[0]), float(q[1]), float(q[2])

        # The code's rule, verbatim.
        adopted = bool(
            abs(ct2) <= 0.25 and abs(ctq) <= 0.8
            and np.isfinite(s2) and (not np.isfinite(s1) or s2 <= s1 + 1e-6)
        )

        # An F-test on the SAME inlier set — the only comparison that means
        # anything, since the two robust fits clip independently.
        lin = np.polyfit(x, y, 1, w=np.sqrt(w))
        r1 = y - np.polyval(lin, x)
        med = np.median(r1)
        sig = 1.4826 * np.median(np.abs(r1 - med))
        common = np.abs(r1 - med) <= 3.0 * max(sig, 1e-6)
        xc, yc, wc = x[common], y[common], w[common]
        c1 = np.polyfit(xc, yc, 1, w=np.sqrt(wc))
        c2 = np.polyfit(xc, yc, 2, w=np.sqrt(wc))
        rss1 = float(np.sum(wc * (yc - np.polyval(c1, xc)) ** 2))
        rss2 = float(np.sum(wc * (yc - np.polyval(c2, xc)) ** 2))
        n = xc.size
        f_stat = (rss1 - rss2) / (rss2 / (n - 3)) if rss2 > 0 and n > 3 else np.nan
        p_val = float(1.0 - stats.f.cdf(f_stat, 1, n - 3)) if np.isfinite(f_stat) else np.nan

        n_saved = int(crow["N"])
        flag = "" if abs(n1 - n_saved) <= max(5, 0.02 * n_saved) else "  ← N 불일치"
        print(f"{cluster:<9} {filt:<3} {n1:7d} {n2:7d} {s1:9.5f} {s2:9.5f} "
              f"{ct1:+8.4f} {ctq:+8.4f} {ct2:+8.4f} {f_stat:7.1f} {p_val:9.2e}"
              f"  {'예' if adopted else '아니오'}{flag}")
        rows.append({"cluster": cluster, "filter": filt, "n1": n1, "n2": n2,
                     "s1": s1, "s2": s2, "ct1": ct1, "ctq": ctq, "ct2": ct2,
                     "F": f_stat, "p": p_val, "adopted": adopted,
                     "saved_ct2": _recorded_ct2(crow)})

df = pd.DataFrame(rows)
print()
print(f"  적합 {len(df)} 개 중 2 차가 산포를 나쁘게 만든 경우: "
      f"{int((df['s2'] > df['s1'] + 1e-6).sum())} 개")
print(f"  「산포 안 나빠짐」 규칙이 채택: {int(df['adopted'].sum())} 개")
print(f"  F 검정 p<0.05 로 정당화되는 것: {int((df['p'] < 0.05).sum())} 개")
print(f"  둘이 엇갈리는 경우: "
      f"{int((df['adopted'] != (df['p'] < 0.05)).sum())} 개")

print("\n  같은 밴드가 성단마다 다른 곡률을 낸다 (같은 카메라·같은 필터):")
for band, group in df.groupby(df["filter"].str.upper()):
    saved = group["saved_ct2"].to_numpy(float)
    print(f"    {band:<2} ct2 = " + ", ".join(f"{v:+.4f}" for v in saved)
          + f"   (폭 {saved.max() - saved.min():.4f})")

df.to_csv("validation/color_term_quadratic_table.csv", index=False)
