"""Compare two reprocess trees star by star.

Written for the 2026-08-07 reprocess, which changed two things at once:

* cosmetic (cosmic-ray / hot-pixel) correction turned on for M67, M5 and M3 —
  M13 and NGC 6811 already had it;
* gain 0.689 -> 0.68 e-/ADU and read noise 2.5 -> 2.35 e- everywhere.

Gain enters the CCD noise equation, not the flux sum, so it must move the
reported errors and leave the magnitudes alone.  That gives a sharp test: on
the two targets whose calibration did not change, the magnitudes must be
essentially identical, and any drift there means something else moved.

Read-only.  Matches on ``source_id`` (Gaia-derived, so stable across runs).

    .venv-deploy/Scripts/python.exe -X utf8 validation/compare_reprocess_runs.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

OLD = Path(r"E:/APEX_validation/reprocess")
NEW = Path(r"E:/APEX_validation/reprocess_cr")
TARGETS = ("NGC6811", "M67", "M13", "M3", "M5")
# Targets that already had cosmetic correction: only the constants changed, so
# their magnitudes are the control sample.
UNCHANGED_CALIBRATION = {"NGC6811", "M13"}
MAG = "mag_inst"
ERR_CANDIDATES = ("mag_err", "mag_inst_err", "magerr")


def mad(x: np.ndarray) -> float:
    return float(np.median(np.abs(x - np.median(x)))) if x.size else float("nan")


def load_frame(path: Path, columns):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    keep = [c for c in columns if c in df.columns]
    return df[keep]


def compare_target(target: str) -> dict | None:
    old_dir = OLD / target / "result" / "step7_forced_phot"
    new_dir = NEW / target / "result" / "step7_forced_phot"
    if not old_dir.is_dir() or not new_dir.is_dir():
        print(f"  [skip] {target}: step7 output missing on one side")
        return None

    err_col = None
    dmag, derr = [], []
    matched = frames = 0

    for new_file in sorted(new_dir.glob("photometry_*.tsv")):
        old_file = old_dir / new_file.name
        if not old_file.exists():
            continue
        if err_col is None:
            head = pd.read_csv(new_file, sep="\t", nrows=1)
            err_col = next((c for c in ERR_CANDIDATES if c in head.columns), None)
        columns = ["source_id", MAG] + ([err_col] if err_col else [])
        a, b = load_frame(old_file, columns), load_frame(new_file, columns)
        merged = a.merge(b, on="source_id", suffixes=("_old", "_new")).dropna(
            subset=[f"{MAG}_old", f"{MAG}_new"])
        if merged.empty:
            continue
        frames += 1
        matched += len(merged)
        dmag.append((merged[f"{MAG}_new"] - merged[f"{MAG}_old"]).to_numpy(float))
        if err_col:
            ratio = merged[f"{err_col}_new"] / merged[f"{err_col}_old"]
            derr.append(ratio.replace([np.inf, -np.inf], np.nan).dropna().to_numpy(float))

    if not dmag:
        print(f"  [skip] {target}: no overlapping measurements")
        return None

    d = np.concatenate(dmag)
    d = d[np.isfinite(d)]
    result = {
        "target": target,
        "frames": frames,
        "matched": int(d.size),
        "dmag_median_mmag": float(np.median(d)) * 1e3,
        "dmag_mad_mmag": mad(d) * 1e3,
        "dmag_p99_mmag": float(np.quantile(np.abs(d), 0.99)) * 1e3,
        "calibration_changed": target not in UNCHANGED_CALIBRATION,
    }
    if derr:
        r = np.concatenate(derr)
        r = r[np.isfinite(r) & (r > 0)]
        result["err_ratio_median"] = float(np.median(r)) if r.size else float("nan")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--old", type=Path, default=OLD)
    parser.add_argument("--new", type=Path, default=NEW)
    parser.add_argument("--out", type=Path, default=None, help="write JSON here")
    args = parser.parse_args()

    globals()["OLD"], globals()["NEW"] = args.old, args.new
    print(f"old: {args.old}\nnew: {args.new}\n")

    rows = []
    for target in TARGETS:
        print(f"{target} …", flush=True)
        row = compare_target(target)
        if row:
            rows.append(row)

    if not rows:
        print("nothing compared")
        return 1

    # gain 0.68 / 0.689: sigma scales as 1/sqrt(N_e), so the reported error
    # should grow by sqrt(0.689/0.68) - 1 = +0.66 %.
    expected_ratio = float(np.sqrt(0.689 / 0.68))

    print(f"\n{'대상':>9} {'프레임':>5} {'대조수':>8} {'Δmag 중앙':>11} {'MAD':>9} "
          f"{'99%':>9} {'오차비':>8}  보정변경")
    print(f"{'':>9} {'':>5} {'':>8} {'(mmag)':>11} {'(mmag)':>9} {'(mmag)':>9} {'':>8}")
    for r in rows:
        print(f"{r['target']:>9} {r['frames']:>5} {r['matched']:>8} "
              f"{r['dmag_median_mmag']:>11.3f} {r['dmag_mad_mmag']:>9.3f} "
              f"{r['dmag_p99_mmag']:>9.1f} "
              f"{r.get('err_ratio_median', float('nan')):>8.4f}"
              f"  {'예' if r['calibration_changed'] else '아니오 (대조군)'}")

    print(f"\n예상 오차비 (gain 0.689→0.68): {expected_ratio:.4f}")

    control = [r for r in rows if not r["calibration_changed"]]
    if control:
        worst = max(abs(r["dmag_median_mmag"]) for r in control)
        print(f"대조군 최대 |Δmag 중앙값|: {worst:.3f} mmag", end="  ")
        print("— 통과 (보정이 안 바뀐 대상은 등급도 안 바뀌어야 한다)"
              if worst < 1.0 else "— 확인 필요")

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(
            {"old": str(args.old), "new": str(args.new),
             "expected_err_ratio": expected_ratio, "targets": rows},
            ensure_ascii=False, indent=1), encoding="utf-8")
        print(f"\nsaved: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
