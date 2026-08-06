"""CR on/off photometric impact — NGC 6811.

Compares APEX forced photometry between the cosmic-ray-corrected run
(result/, cosmetic_enable=True) and the uncorrected run (result_nocr/),
matching star-by-star on source_id within each frame.

Read-only: touches nothing but the two existing result trees.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(r"E:/APEX_validation/reprocess/NGC6811")
ON = BASE / "result" / "step7_forced_phot"
OFF = BASE / "result_nocr" / "step7_forced_phot"


def mad(x: np.ndarray) -> float:
    return float(np.median(np.abs(x - np.median(x))))


def pick_mag_column(cols) -> str:
    for cand in ("mag_inst", "mag", "mag_aper", "magnitude"):
        if cand in cols:
            return cand
    mags = [c for c in cols if c.lower().startswith("mag") and "err" not in c.lower()]
    if not mags:
        raise SystemExit(f"no magnitude column found; columns = {list(cols)[:60]}")
    return mags[0]


def main() -> None:
    frames = sorted(p.name for p in ON.glob("photometry_*.tsv"))
    if not frames:
        raise SystemExit(f"no photometry TSVs under {ON}")

    rows = []
    per_frame = []
    magcol = None

    for name in frames:
        f_on, f_off = ON / name, OFF / name
        if not f_off.exists():
            print(f"  [skip] {name}: no counterpart in result_nocr")
            continue

        a = pd.read_csv(f_on, sep="\t", low_memory=False)
        b = pd.read_csv(f_off, sep="\t", low_memory=False)
        if magcol is None:
            magcol = pick_mag_column(a.columns)
            print(f"magnitude column: {magcol}\n")

        keep = ["source_id", magcol]
        if "snr" in a.columns:
            keep.append("snr")
        a, b = a[keep].copy(), b[keep].copy()

        m = a.merge(b, on="source_id", suffixes=("_on", "_off"))
        m = m.dropna(subset=[f"{magcol}_on", f"{magcol}_off"])
        if m.empty:
            continue

        d = (m[f"{magcol}_on"] - m[f"{magcol}_off"]).to_numpy(float)
        finite = np.isfinite(d)
        d, m = d[finite], m.loc[finite]
        if d.size == 0:
            continue

        changed = int(np.sum(np.abs(d) > 1e-6))
        per_frame.append(
            (name.replace("photometry_pp_NGC6811-", "").replace(".fit.tsv", ""),
             len(d), changed, float(np.median(d)) * 1e3, mad(d) * 1e3,
             float(np.max(np.abs(d))) * 1e3)
        )
        rows.append(pd.DataFrame({"d": d, "mag": m[f"{magcol}_off"].to_numpy(float)}))

    if not rows:
        raise SystemExit("no overlapping stars found")

    print(f"{'frame':>12} {'N':>6} {'changed':>8} {'median':>9} {'MAD':>8} {'|max|':>9}")
    print(f"{'':>12} {'':>6} {'':>8} {'(mmag)':>9} {'(mmag)':>8} {'(mmag)':>9}")
    for r in per_frame:
        print(f"{r[0]:>12} {r[1]:>6} {r[2]:>8} {r[3]:>9.3f} {r[4]:>8.3f} {r[5]:>9.1f}")

    allr = pd.concat(rows, ignore_index=True)
    d = allr["d"].to_numpy()
    n_changed = int(np.sum(np.abs(d) > 1e-6))

    print("\n=== ALL FRAMES POOLED ===")
    print(f"  matched measurements : {len(d)}")
    print(f"  changed (>1e-6 mag)  : {n_changed}  ({100*n_changed/len(d):.2f}%)")
    print(f"  median delta         : {np.median(d)*1e3:+.4f} mmag")
    print(f"  MAD                  : {mad(d)*1e3:.4f} mmag")
    print(f"  RMS                  : {np.sqrt(np.mean(d**2))*1e3:.4f} mmag")
    print(f"  max |delta|          : {np.max(np.abs(d))*1e3:.2f} mmag")
    for q in (0.90, 0.99, 0.999):
        print(f"  |delta| {q*100:5.1f} pct    : {np.quantile(np.abs(d), q)*1e3:.4f} mmag")

    print("\n=== magnitude dependence (CR-off mag bins) ===")
    finite = np.isfinite(allr["mag"])
    sub = allr[finite]
    if len(sub) > 100:
        edges = np.nanquantile(sub["mag"], np.linspace(0, 1, 9))
        print(f"{'mag range':>16} {'N':>7} {'median':>10} {'MAD':>9}")
        for lo, hi in zip(edges[:-1], edges[1:]):
            s = sub[(sub["mag"] >= lo) & (sub["mag"] < hi)]
            if len(s) < 10:
                continue
            dd = s["d"].to_numpy()
            print(f"{lo:7.2f}–{hi:6.2f} {len(dd):7d} "
                  f"{np.median(dd)*1e3:+10.3f} {mad(dd)*1e3:9.3f}")


if __name__ == "__main__":
    sys.exit(main())
