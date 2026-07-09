"""Per-step preprocessing cross-check: APEX vs astropy ``ccdproc``.

``ccdproc`` (astropy-affiliated) is the community-standard Python CCD reduction
package.  This module runs every APEX calibration stage — master bias/dark/flat
construction and bias/dark/flat application — and the equivalent ``ccdproc``
operation on the SAME real frames, then compares them pixel-for-pixel.  It is the
calibration analogue of APEX's photometry cross-checks against ``sep`` (Fig 4)
and IRAF/DAOPHOT (Fig 5): an *independent* implementation, not the author's own
imaging tool, so agreement is evidence the arithmetic is standard and correct.

It is NOT a ground-truth validation on its own (two codes could share a bug);
the ground-truth gate is the synthetic inject->recover test in
``calibration_validate.py``.  This quantifies cross-implementation agreement.

Run (needs ``ccdproc`` in the environment)::

    .venv-deploy\\Scripts\\python -m apex.benchmark.calibration_crosscheck \\
        --bias E:\\bias --dark <dark_dir> --flat <flat_dir> --light <light.fit> \\
        --dark-exp 60 --filter B --out crosscheck.json
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from pathlib import Path
from typing import Callable, Dict, List, Optional

import numpy as np
from astropy.io import fits

from apex.analysis import calibration as cal
from apex.analysis.calibration import CalibrationOptions


# --- steps we cross-check, in pipeline order ------------------------------
STEP_ORDER = ("master_bias", "master_dark", "master_flat",
              "bias_subtract", "dark_subtract", "flat_correct", "full_pipeline")


def _pick(pattern: str, n: int, pred: Callable) -> List[str]:
    out: List[str] = []
    for f in sorted(glob.glob(pattern, recursive=True)):
        try:
            h = fits.getheader(f)
        except Exception:
            continue
        if pred(h):
            out.append(f)
            if len(out) >= n:
                break
    return out


def _compare(a, b) -> Dict[str, float]:
    a = np.asarray(a, np.float64)
    b = np.asarray(b, np.float64)
    d = a - b
    fin = np.isfinite(d)
    med = float(np.median(d[fin]))
    sigma = float(1.4826 * np.median(np.abs(d[fin] - med)))
    return {
        "delta_median": med,
        "robust_sigma": sigma,
        "max_abs": float(np.nanmax(np.abs(d[fin]))),
        "frac_gt_1e3_pct": float(100.0 * np.mean(np.abs(d[fin]) > 1e-3)),
    }


def run_crosscheck(bias: List[str], dark: List[str], flat: List[str],
                   light: str, opts: Optional[CalibrationOptions] = None) -> Dict:
    """Cross-check every calibration stage against ccdproc. Returns a result
    dict keyed by step. Requires ``ccdproc`` (imported lazily)."""
    import ccdproc
    from astropy.nddata import CCDData
    import astropy.units as u

    opts = opts or CalibrationOptions(combine_method="median",
                                      pedestal_mode="none", cosmetic_enable=False)
    CC = lambda arr: CCDData(np.asarray(arr, np.float64), unit="adu")
    results: Dict[str, Dict] = {}

    # 1) master bias — APEX median combine vs ccdproc.combine(median)
    apex_mb, _ = cal.build_master_bias(bias, opts)
    cc_mb = ccdproc.combine(bias, method="median", unit="adu").data
    results["master_bias"] = _compare(apex_mb, cc_mb)

    # 2) master dark — both bias-subtracted with the same master bias
    apex_md, dexp, _ = cal.build_master_dark(dark, opts, master_bias=apex_mb)
    cc_darks = [ccdproc.subtract_bias(CCDData.read(p, unit="adu"), CC(apex_mb))
                for p in dark]
    cc_md = ccdproc.combine(cc_darks, method="median").data
    results["master_dark"] = _compare(apex_md, cc_md)

    # 3) master flat — per-frame median scale, combine, renormalise to unit median
    apex_mf, _ = cal.build_master_flat(flat, opts, master_bias=apex_mb)
    cc_flats = [ccdproc.subtract_bias(CCDData.read(p, unit="adu"), CC(apex_mb))
                for p in flat]
    cc_mf = ccdproc.combine(cc_flats, method="median",
                            scale=lambda a: 1.0 / np.ma.median(a)).data
    cc_mf = cc_mf / np.nanmedian(cc_mf)
    apex_mf_n = apex_mf / np.nanmedian(apex_mf)
    results["master_flat"] = _compare(apex_mf_n, cc_mf)

    # shared masters for the application steps
    mb, md, mf = apex_mb, apex_md, apex_mf_n
    raw = fits.getdata(light).astype(np.float64)
    lexp = float(fits.getheader(light).get("EXPTIME", dexp))

    apex_b = raw - mb
    cc_b = ccdproc.subtract_bias(CC(raw), CC(mb)).data
    results["bias_subtract"] = _compare(apex_b, cc_b)

    ratio = lexp / dexp if dexp else 1.0
    apex_d = apex_b - md * ratio
    cc_d = ccdproc.subtract_dark(CC(cc_b), CC(md), scale=True,
                                 dark_exposure=dexp * u.s,
                                 data_exposure=lexp * u.s).data
    results["dark_subtract"] = _compare(apex_d, cc_d)

    apex_f = apex_d / mf
    cc_f = ccdproc.flat_correct(CC(cc_d), CC(mf), norm_value=1.0).data
    results["flat_correct"] = _compare(apex_f, cc_f)

    apex_full, _, _ = cal.calibrate_light_file(light, opts, master_bias=mb,
                                               master_dark=md, dark_exp=dexp,
                                               master_flat=mf)
    results["full_pipeline"] = _compare(apex_full, cc_f)

    return {
        "reference": "astropy ccdproc",
        "n_bias": len(bias), "n_dark": len(dark), "n_flat": len(flat),
        "dark_exp": dexp, "light": os.path.basename(light),
        "light_median": float(np.nanmedian(np.asarray(apex_full, np.float64))),
        "steps": results,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="APEX vs ccdproc per-step cross-check")
    ap.add_argument("--bias", required=True)
    ap.add_argument("--dark", required=True)
    ap.add_argument("--flat", required=True)
    ap.add_argument("--light", required=True)
    ap.add_argument("--filter", default="B")
    ap.add_argument("--dark-exp", type=float, default=60.0)
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)

    bias = _pick(os.path.join(a.bias, "*.fit*"), a.n, lambda h: True)
    dark = _pick(os.path.join(a.dark, "**", "*.fit*"), a.n,
                 lambda h: float(h.get("EXPTIME", 0)) == a.dark_exp)
    flat = _pick(os.path.join(a.flat, "**", "*.fit*"), a.n,
                 lambda h: str(h.get("FILTER")) == a.filter)
    res = run_crosscheck(bias, dark, flat, a.light)

    print(f"APEX vs {res['reference']}  (bias {res['n_bias']} / dark {res['n_dark']} "
          f"/ flat {res['n_flat']} / light {res['light']})")
    for step in STEP_ORDER:
        s = res["steps"][step]
        print(f"  {step:16s} Δmed={s['delta_median']:+.2e}  σ={s['robust_sigma']:.2e}"
              f"  max|Δ|={s['max_abs']:.2e}  >1e-3DN={s['frac_gt_1e3_pct']:.4f}%")
    worst = max(s["frac_gt_1e3_pct"] for s in res["steps"].values())
    print("ALL STEPS agree to <1e-3 DN" if worst < 0.01 else "SOME steps exceed 1e-3 DN")

    if a.out:
        Path(a.out).write_text(json.dumps(res, indent=2), encoding="utf-8")
        print("wrote", a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
