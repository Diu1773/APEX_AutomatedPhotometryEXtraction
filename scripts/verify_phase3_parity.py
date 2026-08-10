"""Phase 3 clean run vs the preserved baseline: does adopting the optimisations
change the science?

The clean run is the first full reprocess with all three adopted changes live —
streaming master combine (O1), per-stage worker caps (O2), and process-pool
forced photometry. Each was parity-checked in isolation; this checks them
together, end to end, on every target.

The gates are not the same at every stage, and that is deliberate:

* Step 0 calibrated frames -> BIT-IDENTICAL. Master combine is single-threaded
  and streaming was proven byte-for-byte equal (30/30 frames, max |delta| = 0).
  Any difference here is a defect, not noise.

* Step 7 photometry -> STATISTICAL. Bit equality is unavailable downstream of
  source detection: SEP's deblending is not deterministic for blended sources,
  which moves ~2.7 % of measurements with a median |delta| of 0.016 mmag and a
  15.9 mmag tail (benchmark/perf/20260807/RESULTS.md, B4). The thresholds below
  are set from that measured noise floor, not chosen for convenience.

Reports every target; exits non-zero if any gate fails, so it can gate a commit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

# Measured noise floor from B4 (same input, three repeats). A run that stays
# inside these has not changed the science; one that leaves them has.
GATE_CHANGED_FRACTION = 0.05      # B4 measured 2.57 % between repeats
GATE_MEDIAN_MAD_MMAG = 0.1        # B4 measured 0.0000
GATE_MAX_ABS_MMAG = 30.0          # B4 measured 15.9
GATE_NDETECT_TOL = 5
GATE_CATALOG_TOL = 2


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def compare_calibrated(base: Path, new: Path) -> dict:
    """Step 0 products must be byte-identical."""
    base_files = {p.name: p for p in sorted(base.glob("*.fit*"))}
    new_files = {p.name: p for p in sorted(new.glob("*.fit*"))}
    common = sorted(set(base_files) & set(new_files))
    mismatched = [n for n in common
                  if sha256(base_files[n]) != sha256(new_files[n])]
    return {
        "n_baseline": len(base_files),
        "n_new": len(new_files),
        "n_common": len(common),
        "n_mismatched": len(mismatched),
        "mismatched": mismatched[:10],
        "pass": len(mismatched) == 0 and len(common) > 0,
    }


def compare_photometry(base_dir: Path, new_dir: Path) -> dict:
    """Step 7 products are matched by source_id and judged statistically."""
    def load(d: Path) -> pd.DataFrame | None:
        frames = []
        for tsv in sorted(d.glob("*.tsv")):
            try:
                t = pd.read_csv(tsv, sep="\t")
            except Exception:
                continue
            if {"source_id", "mag_inst"} <= set(t.columns):
                t["__frame"] = tsv.name
                frames.append(t[["__frame", "source_id", "mag_inst"]])
        return pd.concat(frames, ignore_index=True) if frames else None

    b, n = load(base_dir), load(new_dir)
    if b is None or n is None:
        return {"pass": False, "reason": "photometry tables not found",
                "baseline_found": b is not None, "new_found": n is not None}

    merged = b.merge(n, on=["__frame", "source_id"], suffixes=("_b", "_n"))
    if merged.empty:
        return {"pass": False, "reason": "no source_id overlap",
                "n_baseline": len(b), "n_new": len(n)}

    delta = (merged["mag_inst_n"] - merged["mag_inst_b"]).to_numpy(float)
    delta = delta[np.isfinite(delta)]
    changed = delta != 0
    frac = float(changed.mean())
    med_mad = float(1.4826 * np.median(np.abs(delta - np.median(delta)))) * 1000
    max_abs = float(np.max(np.abs(delta))) * 1000 if delta.size else 0.0
    med_changed = (float(np.median(np.abs(delta[changed]))) * 1000
                   if changed.any() else 0.0)

    return {
        "n_matched": int(len(merged)),
        "changed_fraction": round(frac, 4),
        "median_mad_mmag": round(med_mad, 4),
        "median_abs_of_changed_mmag": round(med_changed, 4),
        "max_abs_mmag": round(max_abs, 3),
        "pass": (frac <= GATE_CHANGED_FRACTION
                 and med_mad <= GATE_MEDIAN_MAD_MMAG
                 and max_abs <= GATE_MAX_ABS_MMAG),
    }


def compare_target(baseline_root: Path, new_root: Path, target: str) -> dict:
    b = baseline_root / target
    n = new_root / target
    out: dict = {"target": target}
    if not n.exists():
        return {**out, "pass": False, "reason": "clean-run target missing"}
    if not b.exists():
        return {**out, "pass": None, "reason": "no baseline to compare"}

    out["calibrated"] = compare_calibrated(b / "sci", n / "sci")
    out["photometry"] = compare_photometry(
        b / "result" / "step7_forced_phot", n / "result" / "step7_forced_phot")
    checks = [v.get("pass") for v in (out["calibrated"], out["photometry"])]
    out["pass"] = all(c is True for c in checks)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", type=Path,
                    default=Path(r"E:\APEX_validation\reprocess"))
    ap.add_argument("--new", type=Path,
                    default=Path(r"E:\APEX_validation\phase3"))
    ap.add_argument("--targets", default="NGC6811,M67,M13,M3,M5")
    ap.add_argument("--output", type=Path,
                    default=REPO / "benchmark" / "perf" / "20260811"
                    / "phase3_parity.json")
    a = ap.parse_args()

    results = [compare_target(a.baseline, a.new, t)
               for t in a.targets.split(",") if t]

    for r in results:
        verdict = {True: "PASS", False: "FAIL", None: "n/a"}[r.get("pass")]
        print(f"\n=== {r['target']}: {verdict} ===")
        if "reason" in r:
            print(f"  {r['reason']}")
        cal = r.get("calibrated")
        if cal:
            print(f"  calibrated: {cal['n_common']} common, "
                  f"{cal['n_mismatched']} byte-mismatched "
                  f"({'PASS' if cal['pass'] else 'FAIL'})")
            if cal["mismatched"]:
                print(f"    first mismatches: {cal['mismatched']}")
        ph = r.get("photometry")
        if ph and "n_matched" in ph:
            print(f"  photometry: {ph['n_matched']:,} matched | "
                  f"changed {ph['changed_fraction']*100:.2f}% | "
                  f"median MAD {ph['median_mad_mmag']:.4f} mmag | "
                  f"max |delta| {ph['max_abs_mmag']:.2f} mmag "
                  f"({'PASS' if ph['pass'] else 'FAIL'})")
        elif ph:
            print(f"  photometry: {ph.get('reason')}")

    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(
        {"gates": {"changed_fraction": GATE_CHANGED_FRACTION,
                   "median_mad_mmag": GATE_MEDIAN_MAD_MMAG,
                   "max_abs_mmag": GATE_MAX_ABS_MMAG},
         "results": results}, indent=1), encoding="utf-8")
    print(f"\nsaved -> {a.output}")

    failed = [r["target"] for r in results if r.get("pass") is False]
    if failed:
        print(f"FAILED: {', '.join(failed)}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
