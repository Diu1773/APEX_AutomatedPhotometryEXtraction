"""Generate the real M13 before/after cutout that fig10_calibration.py consumes.

Builds master bias/dark/flat for the M13 V session, calibrates one V light with
APEX, and saves a centre cutout of the raw and calibrated frames (small .npy) to
validation/paper/data/ so the figure does not need the E:\\ raw data at plot time.

Run with the project venv:
    .venv-deploy/Scripts/python -m validation.paper._make_calib_figdata
"""

from __future__ import annotations

import glob
from pathlib import Path

import numpy as np
from astropy.io import fits

from apex.analysis import calibration as cal
from apex.analysis.calibration import CalibrationOptions

RAW = r"E:\observe_raw_Analysis\M13_20260515"
BIAS = r"E:\bias"
OUT = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction\validation\paper\data")
CROP = 640


def _center_crop(a, n):
    h, w = a.shape
    y0, x0 = (h - n) // 2, (w - n) // 2
    return a[y0:y0 + n, x0:x0 + n]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    bias = sorted(glob.glob(BIAS + r"\bias-*.fit"))[:12]
    dark = sorted(glob.glob(RAW + r"\Dark\*.fit"))[:8]
    flat = sorted(glob.glob(RAW + r"\flats\flat_V_*.fit"))[:5]
    light = sorted(glob.glob(RAW + r"\M13\messier13-0001-V.fit"))[0]

    opts = CalibrationOptions(combine_method="median", pedestal_mode="adaptive")
    mb, _ = cal.build_master_bias(bias, opts)
    md, dexp, _ = cal.build_master_dark(dark, opts, master_bias=mb)
    mf, _ = cal.build_master_flat(flat, opts, master_bias=mb, master_dark=md, dark_exp=dexp)

    raw = fits.getdata(light).astype(np.float32)
    cald, _hdr, _qc = cal.calibrate_light_file(
        light, opts, master_bias=mb, master_dark=md, dark_exp=dexp, master_flat=mf)

    np.save(OUT / "calib_m13_raw.npy", _center_crop(raw, CROP).astype(np.float32))
    np.save(OUT / "calib_m13_cal.npy", _center_crop(np.asarray(cald, np.float32), CROP))

    # Full-width column-median profiles (star-robust): the master flat traces
    # the instrument's optical vignette, and the calibrated sky is flat after
    # dividing by it — the direct before/after of flat-fielding.
    cal_arr = np.asarray(cald, np.float64)
    flat_prof = np.nanmedian(mf.astype(np.float64), axis=0)   # vignette shape
    cal_prof = np.nanmedian(cal_arr, axis=0)                  # calibrated sky
    np.save(OUT / "calib_m13_profile.npy",
            np.vstack([flat_prof, cal_prof]).astype(np.float32))
    print(f"wrote {OUT/'calib_m13_raw.npy'}, calib_m13_cal.npy ({CROP}x{CROP}), "
          f"calib_m13_profile.npy ({flat_prof.size} cols)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
