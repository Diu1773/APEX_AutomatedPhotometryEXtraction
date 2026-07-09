"""Multi-dataset APEX-vs-AIPPI calibration equivalence (real data).

Extends the single-night M13 check (calibration_equivalence.py) to several real
targets / nights / filters / exposures.  For each dataset it auto-selects a
light frame whose filter has flats and whose (exposure, temperature) has darks,
builds master bias/dark/flat with BOTH APEX (apex.analysis.calibration) and the
AstralImage/AIPPI engine (core.engine.CalibrationEngine), calibrates that light
with both, and reports per-pixel agreement.

Needs the external raw data on the E: drive and the AstralImage repo checked out
beside this one.  Run with the project venv:

    .venv-deploy/Scripts/python -m validation.paper.calibration_equivalence_multi
"""

from __future__ import annotations

import glob
import json
import sys
import warnings
from pathlib import Path

import numpy as np

_ASTRAL = r"C:\Users\bmffr\Desktop\Result\Astralimage"
if _ASTRAL not in sys.path:
    sys.path.insert(0, _ASTRAL)

from apex.analysis import calibration as cal                     # noqa: E402
from apex.analysis import calibration_scan as scan               # noqa: E402
from apex.analysis.calibration import CalibrationOptions         # noqa: E402

BIAS_DIR = r"E:\bias"
RAW = r"E:\observe_raw_Analysis"

DATASETS = [
    ("M13", rf"{RAW}\M13_20260515"),
    ("M67", rf"{RAW}\M67_20260208"),
    ("M3", rf"{RAW}\M3"),
    ("AE_UMa", rf"{RAW}\AE UMa"),
    ("20260611", rf"{RAW}\20260611"),
    ("YZ_Boo", rf"{RAW}\YZbootis"),
]

N_BIAS, N_DARK, N_FLAT = 12, 6, 5


def _rms(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    d = a - b
    fin = np.isfinite(d)
    if not fin.any():
        return float("nan"), float("nan")
    return float(np.sqrt(np.mean(d[fin] ** 2))), float(np.max(np.abs(d[fin])))


def _pick_lights(g, max_filters):
    """One light per distinct filter that has both flats and a matching dark."""
    seen, chosen = set(), []
    for l in g["light"]:
        if l.filt in seen:
            continue
        if scan.match_flat(g["flat"], l.filt) and scan.match_dark(g["dark"], l.exp, l.temp):
            seen.add(l.filt)
            chosen.append(l)
        if len(chosen) >= max_filters:
            break
    return chosen


def _meta(paths, exp, filt=""):
    return [{"path": p, "exp": exp, "filter": filt, "filter_raw": filt} for p in paths]


def equiv_rows(name, root, eng, LightCtx, max_filters=4):
    """One equivalence row per matchable filter (first matchable night)."""
    from astropy.io import fits
    frames = scan.scan_folder(root) + scan.scan_folder(BIAS_DIR)
    opts = CalibrationOptions(combine_method="median", pedestal_mode="adaptive")
    out = Path.cwd()
    for night in scan.nights(frames):
        g = scan.group_for_night(frames, night)
        if not g["bias"]:
            continue
        lights = _pick_lights(g, max_filters)
        if not lights:
            continue
        bias_p = [f.path for f in g["bias"]][:N_BIAS]
        A_b, _ = cal.build_master_bias(bias_p, opts)
        I_b = eng.stack_images(_meta(bias_p, 0.0), combine_method="median", step_name="bias")
        bias_rms, _ = _rms(A_b, I_b)

        dark_cache = {}
        rows = []
        for l in lights:
            dk = scan.match_dark(g["dark"], l.exp, l.temp)
            ff = scan.match_flat(g["flat"], l.filt)
            if dk not in dark_cache:
                dark_p = [f.path for f in g["dark"][dk]][:N_DARK]
                A_d, dexp, _ = cal.build_master_dark(dark_p, opts, master_bias=A_b)
                I_d = eng.stack_images(_meta(dark_p, l.exp), master_bias=I_b,
                                       combine_method="median", step_name="dark")
                dark_cache[dk] = (A_d, dexp, I_d)
            A_d, dexp, I_d = dark_cache[dk]
            flat_p = [f.path for f in g["flat"][ff]][:N_FLAT]
            A_f, _ = cal.build_master_flat(flat_p, opts, master_bias=A_b,
                                           master_dark=A_d, dark_exp=dexp)
            I_f = eng.stack_images(_meta(flat_p, dk[0], ff), master_bias=I_b, master_dark=I_d,
                                   dark_exp=l.exp, is_flat=True, combine_method="median",
                                   step_name="flat")
            A_cal, _, _ = cal.calibrate_light_file(
                l.path, opts, master_bias=A_b, master_dark=A_d, dark_exp=dexp, master_flat=A_f)
            eng.output_prefix = "ai"
            ctx = LightCtx(
                meta={"path": l.path, "exp": l.exp, "filter": l.filt,
                      "filter_raw": l.filt, "date": night},
                master_bias=A_b, master_dark=A_d, master_dark_exp=dexp, master_flat=A_f,
                pedestal_mode="adaptive", output_dir=str(out), session_date=night,
                defect_map=None, cc_apply_to="lights", mask_source="none", append_exp=False,
            )
            _ok, oname = eng.process_single_light(ctx)
            I_cal = fits.getdata(out / oname).astype(float)
            dark_rms, _ = _rms(A_d, I_d)
            flat_rms, _ = _rms(A_f / np.nanmedian(A_f), I_f / np.nanmedian(I_f))
            cal_rms, cal_max = _rms(A_cal, I_cal)
            rows.append({
                "dataset": name, "night": night, "filter": l.filt, "exp": l.exp,
                "n_bias": len(bias_p), "n_dark": min(N_DARK, len(g["dark"][dk])),
                "n_flat": len(flat_p),
                "master_bias_rms": bias_rms, "master_dark_rms": dark_rms,
                "master_flat_rms": flat_rms, "calibrated_rms": cal_rms,
                "calibrated_max": cal_max,
            })
        return rows
    return [{"dataset": name, "error": "no matchable night/light"}]


def main() -> int:
    warnings.filterwarnings("ignore")
    from core.engine import CalibrationEngine
    from core.defs import LightProcessingContext
    eng = CalibrationEngine(logger_func=lambda m: None)
    eng.cc_enable = False
    eng.overscan_enable = False
    eng.dark_optimize_enable = False

    rows = []
    print(f"{'dataset':10s} {'night':9s} {'filt':5s} {'exp':>6s} "
          f"{'biasRMS':>9s} {'darkRMS':>9s} {'flatRMS':>9s} {'calRMS':>9s} {'calMax':>9s}")
    print("-" * 92, flush=True)
    for name, root in DATASETS:
        try:
            drows = equiv_rows(name, root, eng, LightProcessingContext)
        except Exception as exc:
            drows = [{"dataset": name, "error": f"{type(exc).__name__}: {exc}"}]
        for r in drows:
            rows.append(r)
            if r.get("error"):
                print(f"{name:10s} ERROR: {r['error']}", flush=True)
            else:
                print(f"{r['dataset']:10s} {r['night']:9s} {r['filter']:5s} {r['exp']:6.0f} "
                      f"{r['master_bias_rms']:9.3f} {r['master_dark_rms']:9.3f} "
                      f"{r['master_flat_rms']:9.2e} {r['calibrated_rms']:9.3f} "
                      f"{r['calibrated_max']:9.3f}", flush=True)
    Path(__file__).with_name("calibration_equivalence_multi.json").write_text(json.dumps(rows, indent=2))
    print("-" * 92)
    ok = [r for r in rows if not r.get("error")]
    filts = sorted({r["filter"] for r in ok})
    print(f"{len(ok)} (dataset,filter) comparisons across filters {filts}. "
          f"max calibrated max|diff| = {max((r['calibrated_max'] for r in ok), default=float('nan')):.3f} DN")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
