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
# The CMD table is gated against the zero-point fit's own scatter, which is
# 20-24 mmag RMS on these targets. A fifth of that is comfortably "does not
# change the science" while still catching a real regression.
GATE_CMD_MMAG = 5.0
GATE_NDETECT_TOL = 5
GATE_CATALOG_TOL = 2


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def is_wcs_card(key: str) -> bool:
    """True for any card the astrometric step owns.

    Three families, all written by Step 5: the SIP distortion terms (A_i_j,
    B_i_j and their inverses), the core WCS keywords, and APEX's own record of
    how the solve went (WCSNST = stars used, WCSRMD = median residual, …).
    Measured on NGC 6811, the two frames whose headers differed between runs
    differed in 26 and 29 cards and every one of them was in this set.
    """
    return (key[:2] in ("A_", "B_")
            or key.startswith(("AP_", "BP_", "CD", "CRPIX", "CRVAL", "PC",
                               "CDELT", "CTYPE", "CUNIT", "WCS", "LONPOLE",
                               "LATPOLE", "EQUINOX", "RADESYS")))


def compare_calibrated(base: Path, new: Path) -> dict:
    """Step 0 products must be byte-identical — in the PIXELS.

    Hashing the whole file is the wrong test and the first run proved it: two
    of NGC 6811's twenty-one frames differed while every one of the 42,000
    photometric measurements was identical to the last bit. The differing cards
    were the SIP distortion coefficients (A_i_j, B_i_j) — the WCS solution that
    Step 5 writes back into the calibrated frame's header. A whole-file hash
    therefore does not test Step 0 at all; it silently tests the astrometric
    solver's run-to-run determinism, which is downstream of SEP's documented
    non-deterministic deblending (B4) and is not what this gate is for.

    So the gate compares pixel data, and WCS/header differences are counted and
    reported as an observation rather than a failure.
    """
    from astropy.io import fits

    base_files = {p.name: p for p in sorted(base.glob("*.fit*"))}
    new_files = {p.name: p for p in sorted(new.glob("*.fit*"))}
    common = sorted(set(base_files) & set(new_files))

    pixel_mismatch, header_only, wcs_only = [], [], []
    for name in common:
        with fits.open(base_files[name]) as hb, fits.open(new_files[name]) as hn:
            if not np.array_equal(hb[0].data, hn[0].data):
                pixel_mismatch.append(name)
                continue
            kb, kn = hb[0].header, hn[0].header
            diff = {k for k in set(kb) | set(kn)
                    if str(kb.get(k)) != str(kn.get(k))}
        if diff:
            header_only.append(name)
            if all(is_wcs_card(k) for k in diff):
                wcs_only.append(name)

    return {
        "n_baseline": len(base_files),
        "n_new": len(new_files),
        "n_common": len(common),
        "n_pixel_mismatched": len(pixel_mismatch),
        "pixel_mismatched": pixel_mismatch[:10],
        "n_header_differs": len(header_only),
        "n_header_differs_wcs_only": len(wcs_only),
        "header_differs": header_only[:10],
        "pass": len(pixel_mismatch) == 0 and len(common) > 0,
    }


def compare_photometry(base_dir: Path, new_dir: Path) -> dict:
    """Step 7 products are matched row-for-row and judged statistically.

    The key is `master_id`, not `source_id`, and that distinction produced a
    false failure before it was fixed. In M13 exactly one Gaia `source_id` per
    frame is carried by two master entries 4.7 px apart — two detections in the
    crowded core cross-matched to the same Gaia star. Merging on `source_id`
    then joins baseline row 1 against new row 2 and reports their difference as
    a photometric change: M13 came out at max |delta| = 755 mmag when the two
    values were simply each other, swapped. `master_id` is unique per frame
    (1501 of 1501 rows checked), so it is the honest key.

    The duplicate cross-matches are still counted and reported, because a Gaia
    source claimed by two master entries is worth knowing about even though it
    is not a parity question.
    """
    def load(d: Path) -> pd.DataFrame | None:
        frames = []
        for tsv in sorted(d.glob("*.tsv")):
            try:
                t = pd.read_csv(tsv, sep="\t")
            except Exception:
                continue
            if {"master_id", "mag_inst"} <= set(t.columns):
                t["__frame"] = tsv.name
                cols = ["__frame", "master_id", "mag_inst"]
                if "source_id" in t.columns:
                    cols.append("source_id")
                frames.append(t[cols])
        return pd.concat(frames, ignore_index=True) if frames else None

    b, n = load(base_dir), load(new_dir)
    if b is None or n is None:
        return {"pass": False, "reason": "photometry tables not found",
                "baseline_found": b is not None, "new_found": n is not None}

    dup_srcid = (int(b.duplicated(["__frame", "source_id"]).sum())
                 if "source_id" in b.columns else 0)
    merged = b.merge(n, on=["__frame", "master_id"], suffixes=("_b", "_n"))
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
        "n_duplicate_source_id_rows": dup_srcid,
        "changed_fraction": round(frac, 4),
        "median_mad_mmag": round(med_mad, 4),
        "median_abs_of_changed_mmag": round(med_changed, 4),
        "max_abs_mmag": round(max_abs, 3),
        "pass": (frac <= GATE_CHANGED_FRACTION
                 and med_mad <= GATE_MEDIAN_MAD_MMAG
                 and max_abs <= GATE_MAX_ABS_MMAG),
    }


def compare_cmd_table(base: Path, new: Path) -> dict:
    """The CMD table is where the chain actually ends for the paper.

    Step 0 and Step 7 parity do not by themselves guarantee it: the zero-point
    and colour-term fit sits between them and the table, and it is refitted per
    run. Measured on M67, a Step 7 perturbation of 3.18 % of measurements with
    a 2.04 mmag maximum becomes a 0.70 mmag shift of the whole g sequence —
    every star moves, because they all share one zero point.

    So "fraction changed" is the wrong gate here: it reads 100 % for a zero-point
    shift of less than a millimagnitude. The gate is instead the SIZE of the
    difference measured against the zero-point fit's own scatter (M67: 24 mmag
    RMS). A run-to-run difference far below the fit's intrinsic scatter has not
    changed the science. Single-star outliers are reported but not gated — they
    are the SEP deblending tail already characterised in B4.
    """
    if not base.exists() or not new.exists():
        return {"pass": None, "reason": "CMD table missing",
                "baseline": base.exists(), "new": new.exists()}
    b, n = pd.read_csv(base), pd.read_csv(new)
    key = "ID" if "ID" in b.columns and "ID" in n.columns else None
    if key is None:
        return {"pass": None, "reason": "no ID column to match on"}
    bands = sorted({c for c in b.columns if c.startswith("mag_std_")
                    and not c.endswith("_err")} & set(n.columns))
    if not bands:
        return {"pass": None, "reason": "no mag_std_* columns"}

    merged = b[[key] + bands].merge(n[[key] + bands], on=key,
                                    suffixes=("_b", "_n"))
    per_band, worst_max, worst_frac = {}, 0.0, 0.0
    worst_mad = worst_shift = 0.0
    for band in bands:
        d = (merged[f"{band}_n"] - merged[f"{band}_b"]).to_numpy(float)
        d = d[np.isfinite(d)]
        if not d.size:
            continue
        mx = float(np.max(np.abs(d))) * 1000
        fr = float((d != 0).mean())
        med = float(np.median(d)) * 1000
        mad = float(1.4826 * np.median(np.abs(d - np.median(d)))) * 1000
        per_band[band] = {"n": int(d.size), "changed_fraction": round(fr, 4),
                          "max_abs_mmag": round(mx, 3),
                          "median_shift_mmag": round(med, 4),
                          "median_mad_mmag": round(mad, 4)}
        worst_max, worst_frac = max(worst_max, mx), max(worst_frac, fr)
        worst_mad = max(worst_mad, mad)
        worst_shift = max(worst_shift, abs(med))

    return {
        "n_stars_baseline": int(len(b)), "n_stars_new": int(len(n)),
        "n_matched": int(len(merged)), "bands": per_band,
        "worst_max_abs_mmag": round(worst_max, 3),
        "worst_changed_fraction": round(worst_frac, 4),
        "worst_median_mad_mmag": round(worst_mad, 4),
        "worst_median_shift_mmag": round(worst_shift, 4),
        "gate_mmag": GATE_CMD_MMAG,
        "pass": worst_mad <= GATE_CMD_MMAG and worst_shift <= GATE_CMD_MMAG,
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
    cmd_rel = Path("result") / "cmd_zeropoint" / "median_by_ID_filter_wide_cmd.csv"
    out["cmd_table"] = compare_cmd_table(b / cmd_rel, n / cmd_rel)
    checks = [v.get("pass")
              for v in (out["calibrated"], out["photometry"], out["cmd_table"])]
    # A stage with no baseline to compare (pass is None) must not veto the rest.
    out["pass"] = (all(c is not False for c in checks)
                   and any(c is True for c in checks))
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
                  f"{cal['n_pixel_mismatched']} pixel-mismatched "
                  f"({'PASS' if cal['pass'] else 'FAIL'})")
            if cal["pixel_mismatched"]:
                print(f"    pixel mismatches: {cal['pixel_mismatched']}")
            if cal["n_header_differs"]:
                print(f"    header differs on {cal['n_header_differs']} frame(s), "
                      f"{cal['n_header_differs_wcs_only']} of them WCS/SIP only "
                      f"(observation, not a failure)")
        ph = r.get("photometry")
        if ph and "n_matched" in ph:
            print(f"  photometry: {ph['n_matched']:,} matched | "
                  f"changed {ph['changed_fraction']*100:.2f}% | "
                  f"median MAD {ph['median_mad_mmag']:.4f} mmag | "
                  f"max |delta| {ph['max_abs_mmag']:.2f} mmag "
                  f"({'PASS' if ph['pass'] else 'FAIL'})")
            if ph.get("n_duplicate_source_id_rows"):
                print(f"    note: {ph['n_duplicate_source_id_rows']} row(s) share "
                      f"a Gaia source_id with another master entry "
                      f"(matched on master_id, which is unique)")
        elif ph:
            print(f"  photometry: {ph.get('reason')}")
        cm = r.get("cmd_table")
        if cm and "n_matched" in cm:
            print(f"  CMD table:  {cm['n_matched']:,} stars matched "
                  f"({cm['n_stars_baseline']} vs {cm['n_stars_new']}) | "
                  f"shift {cm['worst_median_shift_mmag']:.2f} | "
                  f"MAD {cm['worst_median_mad_mmag']:.2f} mmag "
                  f"(gate {cm['gate_mmag']}) | outlier max "
                  f"{cm['worst_max_abs_mmag']:.1f} mmag "
                  f"({'PASS' if cm['pass'] else 'FAIL'})")
        elif cm:
            print(f"  CMD table:  {cm.get('reason')}")

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
