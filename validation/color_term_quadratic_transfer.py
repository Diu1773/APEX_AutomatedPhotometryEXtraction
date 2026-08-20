"""Is the quadratic colour term a property of the instrument, or of the field?

A colour term describes how one filter+detector differs from the standard
passband. It belongs to the instrument. The zero point moves with night and
airmass, but the *shape* of the transformation should not move with which
cluster is in the frame.

So: fit the curvature on one field, apply it to another, and see whether it
helps. If the curvature is real, a coefficient measured on M13 should reduce
scatter on M3. If it is the fit absorbing that field's colour distribution,
transferring it will make things worse.

Same camera throughout (Moravian C3-61000), five fields, 2025 season.
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


def load(cluster):
    base = ROOT / cluster / "result" / "cmd_zeropoint"
    cal = base / "gaia_sdss_calibrator_by_ID.csv"
    coeff = base / "zp_fit_coefficients.csv"
    if not cal.exists() or not coeff.exists():
        return {}
    table = pd.read_csv(cal)
    out = {}
    m_qual = np.asarray(gaia_quality_report(table, cstar_nsigma=None)[0], bool)
    for _, row in pd.read_csv(coeff).iterrows():
        filt, pair = str(row["filter"]), str(row["color_col"])
        dcol, ccol, ecol = f"delta_{filt}", f"color_{pair}", f"mag_inst_err_{filt}"
        if dcol not in table.columns or ccol not in table.columns:
            continue
        x = pd.to_numeric(table[ccol], errors="coerce").to_numpy(float)
        y = pd.to_numeric(table[dcol], errors="coerce").to_numpy(float)
        err = pd.to_numeric(table.get(ecol, 0.01), errors="coerce").to_numpy(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            w = 1.0 / np.square(err)
        snr = pd.to_numeric(table.get(f"snr_{filt}", np.inf), errors="coerce").to_numpy(float)
        keep = (np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
                & np.isfinite(snr) & (snr >= SNR_CUT) & m_qual)
        if keep.sum() < 50:
            continue
        out[(filt.upper(), pair)] = (x[keep], y[keep], w[keep])
    return out


def clipped(x, y, w, degree):
    """Fit, sigma-clip once, refit — and report robust scatter on the inliers."""
    c = np.polyfit(x, y, degree, w=np.sqrt(w))
    r = y - np.polyval(c, x)
    med = np.median(r)
    keep = np.abs(r - med) <= 3.0 * max(MAD * np.median(np.abs(r - med)), 1e-6)
    c = np.polyfit(x[keep], y[keep], degree, w=np.sqrt(w[keep]))
    return c, keep


def scatter_with(x, y, w, ct2, keep):
    """Robust scatter after removing a FIXED curvature, refitting only zp+ct."""
    y_flat = y - ct2 * x * x
    c = np.polyfit(x[keep], y_flat[keep], 1, w=np.sqrt(w[keep]))
    r = y_flat[keep] - np.polyval(c, x[keep])
    return float(MAD * np.median(np.abs(r - np.median(r))))


data = {c: load(c) for c in CLUSTERS}
bands = {}
for cluster, fits in data.items():
    for (band, pair) in fits:
        bands.setdefault((band, pair), []).append(cluster)

print("곡률을 한 성단에서 재고 다른 성단에 적용한다 (같은 카메라·같은 필터)\n")
rows = []
for (band, pair), members in sorted(bands.items()):
    if len(members) < 2:
        continue
    print(f"  ── {band} 밴드 ({pair.replace('_', '-')} 색) · 성단 {len(members)} 개")
    own = {}
    for cluster in members:
        x, y, w = data[cluster][(band, pair)]
        c2, keep = clipped(x, y, w, 2)
        own[cluster] = (float(c2[0]), keep, x, y, w)
    header = "      적용대상    자기 곡률   자기산포   곡률0 산포"
    header += "".join(f"  ←{m[:4]:>6}" for m in members)
    print(header)
    for target in members:
        ct2_own, keep, x, y, w = own[target]
        s_own = scatter_with(x, y, w, ct2_own, keep)
        s_lin = scatter_with(x, y, w, 0.0, keep)
        line = (f"      {target:<10} {ct2_own:+9.4f} {s_own:10.5f} {s_lin:11.5f}")
        for donor in members:
            if donor == target:
                line += f"  {'—':>6}"
                continue
            s_cross = scatter_with(x, y, w, own[donor][0], keep)
            # Positive = donor's curvature made this field worse than no
            # curvature at all, in millimagnitudes of robust scatter.
            delta_mmag = (s_cross - s_lin) * 1000.0
            line += f"  {delta_mmag:+6.1f}"
            rows.append({"band": band, "target": target, "donor": donor,
                         "ct2_donor": own[donor][0], "ct2_own": ct2_own,
                         "d_mmag": delta_mmag})
        print(line)
    print()

df = pd.DataFrame(rows)
if not df.empty:
    worse = int((df["d_mmag"] > 0).sum())
    print(f"  남의 곡률을 빌려 쓴 {len(df)} 경우 중 **{worse} 개가 더 나빠졌다** "
          f"(중앙값 {df['d_mmag'].median():+.1f} mmag, 최악 {df['d_mmag'].max():+.1f} mmag)")
    print("  숫자는 「곡률 0 대비」 robust 산포 변화. 양수 = 안 쓰느니만 못하다.")
    df.to_csv("validation/color_term_quadratic_transfer.csv", index=False)
