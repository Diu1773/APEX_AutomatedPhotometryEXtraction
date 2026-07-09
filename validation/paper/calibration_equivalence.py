"""Real-data equivalence check: APEX Step 0 vs the AstralImage/AIPPI engine.

Evidence for the manuscript's raw->science claim. On a real M13 session
(Moravian C3-61000, 4788x3194) it builds master bias/dark/flat and calibrates a
science frame with BOTH APEX (`apex.analysis.calibration`) and the AstralImage
AIPPI engine (`core.engine.CalibrationEngine`, documented to match PixInsight
WBPP within tolerance), then compares per pixel.

This is NOT a committed unit test: it needs the external raw data on E:\ and the
AstralImage repo checked out next to this one. The self-contained honesty gate
is `apex/benchmark/calibration_validate.py` (synthetic injection->recovery).

Result recorded 2026-07-06 (12 bias / 6 dark / 5 flat_V, shared masters):
    master bias   RMS diff = 0.30 DN   (median 512.0 == 512.0)
    master dark   RMS diff = 0.30 DN   (median   1.0 ==   1.0)
    master flat   RMS diff = 1.5e-5    (normalised, median 1.0)
    calibrated light frame: per-pixel diff = 0.000 DN  (bit-identical, 100%)

Usage (paths hard-coded to this machine's layout; edit as needed)::

    python validation/paper/calibration_equivalence.py
"""

from __future__ import annotations

import glob
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits

# The AstralImage repo (AIPPI engine) must be importable.
_ASTRAL = r"C:\Users\bmffr\Desktop\Result\Astralimage"
if _ASTRAL not in sys.path:
    sys.path.insert(0, _ASTRAL)

from apex.analysis import calibration as cal          # noqa: E402
from apex.analysis.calibration import CalibrationOptions  # noqa: E402

RAW = r"E:\observe_raw_Analysis\M13_20260515"
BIAS_DIR = r"E:\bias"
DARK_EXP = 60.0
FLAT_EXP = 3.0


def _cmp(name, A, I):
    A = np.asarray(A, float)
    I = np.asarray(I, float)
    d = A - I
    fin = np.isfinite(d)
    rms = float(np.sqrt(np.nanmean(d[fin] ** 2)))
    mx = float(np.nanmax(np.abs(d[fin])))
    print(f"  {name:16s} APEXmed={np.nanmedian(A):.4f} AIPPImed={np.nanmedian(I):.4f} "
          f"RMSdiff={rms:.3e} max|d|={mx:.3e}")
    return rms, mx


def main() -> int:
    warnings.filterwarnings("ignore")
    from core.engine import CalibrationEngine
    from core.defs import LightProcessingContext

    bias = sorted(glob.glob(BIAS_DIR + r"\bias-*.fit"))[:12]
    dark = sorted(glob.glob(RAW + r"\Dark\*.fit"))[:6]
    flat = sorted(glob.glob(RAW + r"\flats\flat_V_*.fit"))[:5]
    light = sorted(glob.glob(RAW + r"\M13\messier13-*-V.fit"))[0]
    if not (bias and dark and flat):
        print("Raw data not found — edit RAW/BIAS_DIR for your layout.")
        return 2

    opts = CalibrationOptions(combine_method="median", pedestal_mode="adaptive")
    A_b, _ = cal.build_master_bias(bias, opts)
    A_d, dexp, _ = cal.build_master_dark(dark, opts, master_bias=A_b)
    A_f, _ = cal.build_master_flat(flat, opts, master_bias=A_b, master_dark=A_d, dark_exp=dexp)

    eng = CalibrationEngine(logger_func=lambda m: None)
    eng.cc_enable = False
    eng.overscan_enable = False
    eng.dark_optimize_enable = False

    def meta(paths, exp, filt=""):
        return [{"path": p, "exp": exp, "filter": filt, "filter_raw": filt} for p in paths]

    I_b = eng.stack_images(meta(bias, 0.0), combine_method="median", step_name="bias")
    I_d = eng.stack_images(meta(dark, DARK_EXP), master_bias=I_b,
                           combine_method="median", step_name="dark")
    I_f = eng.stack_images(meta(flat, FLAT_EXP, "V"), master_bias=I_b, master_dark=I_d,
                           dark_exp=DARK_EXP, is_flat=True, combine_method="median",
                           step_name="flat")

    print("Master frames (APEX vs AIPPI, same inputs):")
    _cmp("bias", A_b, I_b)
    _cmp("dark", A_d, I_d)
    _cmp("flat(norm)", A_f / np.nanmedian(A_f), I_f / np.nanmedian(I_f))

    # Calibrated light frame with SHARED (APEX) masters — isolates the per-frame
    # arithmetic. Both should be bit-identical.
    out = Path.cwd()
    A_cal, _, _ = cal.calibrate_light_file(light, opts, master_bias=A_b,
                                           master_dark=A_d, dark_exp=dexp, master_flat=A_f)
    eng.output_prefix = "ai"
    ctx = LightProcessingContext(
        meta={"path": light, "exp": DARK_EXP, "filter": "V", "filter_raw": "V", "date": "20260515"},
        master_bias=A_b, master_dark=A_d, master_dark_exp=dexp, master_flat=A_f,
        pedestal_mode="adaptive", output_dir=str(out), session_date="20260515",
        defect_map=None, cc_apply_to="lights", mask_source="none", append_exp=False,
    )
    ok, name = eng.process_single_light(ctx)
    I_cal = fits.getdata(out / name).astype(float)
    print("Calibrated light frame (shared masters):")
    rms, mx = _cmp("light", A_cal, I_cal)
    print(f"  -> {'BIT-IDENTICAL' if mx == 0.0 else 'agree within %.3g DN' % mx}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
