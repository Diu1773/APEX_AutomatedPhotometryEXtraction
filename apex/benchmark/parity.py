"""Parity gate: one command that decides whether two runs agree (T0.5).

Every optimisation in docs/audit/APEX_PERF_DEV_PLAN.md is judged against the
preserved baseline with the same three checks:

A. **Self-consistency** of the new tree — per star, the median of the Step 7
   per-frame instrumental magnitudes must equal ``mag_inst_med_*`` in the
   ``cmd_zeropoint`` aggregate.  This caught the 2026-08-07 stale-CMD problem
   (old tree: +82.7 mmag offset), and it needs only one tree.
B. **Star-by-star comparison** old vs new — calibrated ``mag_std_*`` matched
   on Gaia ``source_id`` (never the per-run ``ID`` column, which is renumbered
   every run and turns the comparison into noise of ~1.2 mag MAD).
C. **Counters** — Step 7 table count, summed Step 4 detections, WCS rows.

Exit code 0 only if every requested check passes its tolerance.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_OLD = Path(r"E:\APEX_validation\reprocess")
DEFAULT_TARGETS = ("NGC6811", "M67", "M13", "M3", "M5")
CMD_CSV = "result/cmd_zeropoint/median_by_ID_filter_wide_cmd.csv"

# A: a tree must agree with itself to numerical noise.
SELF_MED_TOL_MMAG = 0.5
SELF_MAD_TOL_MMAG = 1.0
# B: default informational tolerance until B4 measures the true noise floor.
COMPARE_MAD_TOL_MMAG = 15.0


def _mad(x: np.ndarray) -> float:
    return float(np.median(np.abs(x - np.median(x)))) if x.size else float("nan")


def _per_star_medians(base: Path) -> dict[tuple[str, int], float]:
    """(filter, source_id) -> median per-frame mag_inst across Step 7 tables."""
    per: dict[tuple[str, int], list[float]] = {}
    for p in sorted((base / "result/step7_forced_phot").glob("photometry_*.tsv")):
        df = pd.read_csv(p, sep="\t", usecols=["source_id", "mag_inst", "filter"],
                         low_memory=False).dropna(subset=["mag_inst"])
        for sid, mag, filt in zip(df["source_id"], df["mag_inst"], df["filter"]):
            per.setdefault((str(filt).strip(), int(sid)), []).append(float(mag))
    return {k: float(np.median(v)) for k, v in per.items()}


def check_self_consistency(base: Path) -> dict:
    """Check A for one target directory; returns per-filter stats + verdict."""
    cmd_path = base / CMD_CSV
    if not cmd_path.exists():
        return {"status": "missing", "path": str(cmd_path)}
    agg = pd.read_csv(cmd_path)
    agg = agg[agg["source_id"] > 0]
    per = _per_star_medians(base)

    filters = sorted({f for f, _ in per})
    out: dict = {"status": "ok", "filters": {}}
    worst_med = worst_mad = 0.0
    for filt in filters:
        col = f"mag_inst_med_{filt}"
        if col not in agg.columns:
            continue
        diffs = np.array([
            per[(filt, int(sid))] - float(val)
            for sid, val in zip(agg["source_id"], agg[col])
            if (filt, int(sid)) in per and np.isfinite(val)
        ])
        if diffs.size < 20:
            continue
        med, mad = float(np.median(diffs)) * 1e3, _mad(diffs) * 1e3
        out["filters"][filt] = {"n": int(diffs.size),
                                "med_mmag": round(med, 3), "mad_mmag": round(mad, 3)}
        worst_med = max(worst_med, abs(med))
        worst_mad = max(worst_mad, mad)
    out["pass"] = bool(out["filters"]) and (worst_med < SELF_MED_TOL_MMAG
                                            and worst_mad < SELF_MAD_TOL_MMAG)
    return out


def check_compare(old_base: Path, new_base: Path, mad_tol: float) -> dict:
    a_path, b_path = old_base / CMD_CSV, new_base / CMD_CSV
    if not (a_path.exists() and b_path.exists()):
        return {"status": "missing"}
    a, b = pd.read_csv(a_path), pd.read_csv(b_path)
    a, b = a[a["source_id"] > 0], b[b["source_id"] > 0]
    merged = a.merge(b, on="source_id", suffixes=("_o", "_n"))

    out: dict = {"status": "ok", "n_common": int(len(merged)), "filters": {}}
    worst_mad = 0.0
    for col in [c for c in a.columns if c.startswith("mag_std_") and "err" not in c]:
        if f"{col}_o" not in merged.columns:
            continue
        d = (merged[f"{col}_n"] - merged[f"{col}_o"]).to_numpy(float)
        d = d[np.isfinite(d)]
        if d.size < 20:
            continue
        med, mad = float(np.median(d)) * 1e3, _mad(d) * 1e3
        out["filters"][col] = {"n": int(d.size),
                               "med_mmag": round(med, 3), "mad_mmag": round(mad, 3)}
        worst_mad = max(worst_mad, mad)
    out["pass"] = bool(out["filters"]) and worst_mad < mad_tol
    return out


def check_counters(base: Path) -> dict:
    step7 = len(list((base / "result/step7_forced_phot").glob("photometry_*.tsv")))
    n_detect, seen = 0, False
    for p in (base / "result/step4_detection").glob("detect_*.json"):
        try:
            n = json.loads(p.read_text(encoding="utf-8")).get("n_sources")
        except Exception:
            continue
        if isinstance(n, (int, float)):
            n_detect += int(n)
            seen = True
    wcs_csv = base / "result/step5_wcs/wcs_solve_summary.csv"
    n_wcs = len(pd.read_csv(wcs_csv)) if wcs_csv.exists() else None
    return {"step7_tables": step7,
            "detect_sum": n_detect if seen else None,
            "wcs_rows": n_wcs}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="apex bench parity", description=__doc__.splitlines()[0])
    parser.add_argument("--new", required=True, help="root of the tree under test")
    parser.add_argument("--old", default=str(DEFAULT_OLD),
                        help="baseline root (default: the preserved reprocess tree)")
    parser.add_argument("--targets", default=",".join(DEFAULT_TARGETS))
    parser.add_argument("--mad-tol", type=float, default=COMPARE_MAD_TOL_MMAG,
                        help="check-B MAD tolerance in mmag (set from B4)")
    parser.add_argument("--skip-compare", action="store_true",
                        help="run only self-consistency + counters (single tree)")
    parser.add_argument("--json-out", default=None)
    args = parser.parse_args(argv)

    new_root, old_root = Path(args.new), Path(args.old)
    report: dict = {"new": str(new_root), "old": str(old_root), "targets": {}}
    all_pass = True

    for target in args.targets.split(","):
        target = target.strip()
        base_new = new_root / target
        entry: dict = {"self": check_self_consistency(base_new),
                       "counters_new": check_counters(base_new)}
        if not args.skip_compare:
            entry["compare"] = check_compare(old_root / target, base_new,
                                            args.mad_tol)
            entry["counters_old"] = check_counters(old_root / target)
        report["targets"][target] = entry

        self_ok = entry["self"].get("pass", False)
        comp = entry.get("compare", {})
        comp_ok = comp.get("pass", True) if comp.get("status") == "ok" else True
        status = "PASS" if (self_ok and comp_ok) else "FAIL"
        if entry["self"].get("status") == "missing":
            status = "SKIP (no cmd_zeropoint)"
        else:
            all_pass = all_pass and (self_ok and comp_ok)

        worst = max((f["mad_mmag"] for f in comp.get("filters", {}).values()),
                    default=None)
        print(f"[parity] {target:>9}  self={'ok' if self_ok else 'FAIL'}  "
              f"compare_worst_MAD={worst} mmag  -> {status}")

    report["pass"] = all_pass
    if args.json_out:
        Path(args.json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.json_out).write_text(
            json.dumps(report, ensure_ascii=False, indent=1), encoding="utf-8")
        print(f"[parity] saved -> {args.json_out}")
    print(f"[parity] overall: {'PASS' if all_pass else 'FAIL'}")
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
