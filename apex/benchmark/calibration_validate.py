"""Synthetic detector-calibration validation (Step 0 honesty gate).

Injects a KNOWN bias / dark / flat signature into synthetic frames, builds
masters with :mod:`apex.analysis.calibration`, calibrates a science frame, and
quantifies how well the injected detector signature is recovered:

* master-bias / master-dark / master-flat error vs. the injected truth,
* calibrated-frame residual vs. the true science frame (systematic + scatter),
* hot-pixel removal, and the vignette that a *missing* flat would leave behind.

Runs with ZERO external data and no Qt imports.  Prints a short report and
writes a JSON manifest, mirroring the other ``apex/benchmark`` harnesses so the
result can be cited in the paper's validation section.

Usage::

    python -m apex.benchmark.calibration_validate [--size 256] [--n 15] [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path
from typing import Dict

import numpy as np
from astropy.io import fits

from apex.analysis import calibration as cal
from apex.analysis.calibration import CalibrationOptions


DARK_EXP = 10.0
LIGHT_EXP = 10.0
FLAT_EXP = 10.0
BIAS_BASE = 300.0
DARK_RATE = 0.5
SKY = 120.0
READ_NOISE = 3.0


# ── truth model: raw = bias + dark_rate*exptime + flat * science ──────────────

def _truth(size: int):
    col = np.linspace(0.0, 20.0, size)[np.newaxis, :]
    bias = (BIAS_BASE + col) * np.ones((size, size))

    rate = DARK_RATE * np.ones((size, size))
    hot = [(size // 8, size // 8, 50.0), (size - 10, size // 2, 90.0)]
    for (y, x, v) in hot:
        rate[y, x] = v

    yy, xx = np.mgrid[0:size, 0:size]
    cy = cx = (size - 1) / 2.0
    r2 = (xx - cx) ** 2 + (yy - cy) ** 2
    flat = 1.0 - 0.25 * (r2 / r2.max())
    flat = flat / np.median(flat)

    science = SKY * np.ones((size, size))
    for (y, x, flux) in [(size // 3, size // 2, 6000.0), (2 * size // 3, size // 4, 3000.0)]:
        science += flux * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * 2.0 ** 2))
    return bias, rate, flat, science, hot


def _write(path: Path, data: np.ndarray, exptime: float, rng) -> str:
    noisy = data + rng.normal(0.0, READ_NOISE, size=data.shape)
    hdu = fits.PrimaryHDU(noisy.astype(np.float32))
    hdu.header["EXPTIME"] = float(exptime)
    hdu.header["FILTER"] = "V"
    hdu.writeto(path, overwrite=True)
    return str(path)


def run(size: int = 256, n: int = 15, out_dir: Path | None = None,
        seed: int = 20260706) -> Dict:
    rng = np.random.default_rng(seed)
    bias, rate, flat, science, hot = _truth(size)

    tmp = out_dir or Path(tempfile.mkdtemp(prefix="apex_calib_val_"))
    tmp.mkdir(parents=True, exist_ok=True)

    bias_p, dark_p, flat_p = [], [], []
    for i in range(n):
        bias_p.append(_write(tmp / f"bias_{i}.fits", bias, 0.0, rng))
        dark_p.append(_write(tmp / f"dark_{i}.fits", bias + rate * DARK_EXP, DARK_EXP, rng))
        flat_p.append(_write(tmp / f"flat_{i}.fits", bias + rate * FLAT_EXP + flat * 20000.0, FLAT_EXP, rng))
    light_p = _write(tmp / "light_0.fits", bias + rate * LIGHT_EXP + flat * science, LIGHT_EXP, rng)

    opts = CalibrationOptions(combine_method="median", pedestal_mode="none")
    mbias, _ = cal.build_master_bias(bias_p, opts)
    mdark, dexp, _ = cal.build_master_dark(dark_p, opts, master_bias=mbias)
    mflat, _ = cal.build_master_flat(flat_p, opts, master_bias=mbias, master_dark=mdark, dark_exp=dexp)
    calibrated, _, qc = cal.calibrate_light_file(
        light_p, opts, master_bias=mbias, master_dark=mdark, dark_exp=dexp, master_flat=mflat)
    no_flat, _, _ = cal.calibrate_light_file(
        light_p, opts, master_bias=mbias, master_dark=mdark, dark_exp=dexp)

    def _rms(a):
        v = a[np.isfinite(a)]
        return float(np.sqrt(np.mean(v * v))) if v.size else float("nan")

    resid = calibrated.astype(np.float64) - science
    med_off = float(np.nanmedian(resid))
    scatter = float(np.nanmedian(np.abs(resid - med_off)) * 1.4826)
    corner = (5, 5)
    vign_left = abs(float(no_flat[corner]) - science[corner])
    vign_removed = abs(float(calibrated[corner]) - science[corner])
    hot_ok = all(abs(float(calibrated[y, x]) - science[y, x]) < 50.0 for (y, x, _) in hot)

    report = {
        "size": size, "n_cal_frames": n, "dark_exptime_recovered": dexp,
        "master_bias_err_rms": _rms(mbias - bias),
        "master_dark_err_rms": _rms(mdark - rate * DARK_EXP),
        "master_flat_median": float(np.median(mflat)),
        "master_flat_err_rms": _rms(mflat - flat),
        "calibrated_median_offset": med_off,
        "calibrated_scatter_mad": scatter,
        "single_frame_read_noise": READ_NOISE,
        "vignette_left_if_no_flat_dn": vign_left,
        "vignette_after_flat_dn": vign_removed,
        "hot_pixels_removed": bool(hot_ok),
    }
    # PASS criteria: no systematic offset, scatter at the single-frame floor,
    # flat removes the vignette, hot pixels gone.
    report["PASS"] = bool(
        abs(med_off) < 1.0
        and scatter < 2.0 * READ_NOISE
        and vign_left > 10.0
        and vign_removed < 2.5 * READ_NOISE
        and hot_ok
    )

    (tmp / "calibration_validation.json").write_text(json.dumps(report, indent=2))
    return report


def _print(r: Dict) -> None:
    print("APEX detector-calibration synthetic validation")
    print("-" * 52)
    print(f"  frames/master           : {r['n_cal_frames']}   size {r['size']}x{r['size']}")
    print(f"  dark exptime recovered  : {r['dark_exptime_recovered']:.3f} s")
    print(f"  master bias  err (RMS)  : {r['master_bias_err_rms']:.3f} DN")
    print(f"  master dark  err (RMS)  : {r['master_dark_err_rms']:.3f} DN")
    print(f"  master flat  median     : {r['master_flat_median']:.4f}  (target 1.0)")
    print(f"  master flat  err (RMS)  : {r['master_flat_err_rms']:.4f}")
    print(f"  calibrated offset       : {r['calibrated_median_offset']:+.3f} DN  (target 0)")
    print(f"  calibrated scatter (MAD): {r['calibrated_scatter_mad']:.3f} DN  (floor ~{r['single_frame_read_noise']:.1f})")
    print(f"  vignette if no flat     : {r['vignette_left_if_no_flat_dn']:.2f} DN -> after flat {r['vignette_after_flat_dn']:.2f} DN")
    print(f"  hot pixels removed      : {r['hot_pixels_removed']}")
    print("-" * 52)
    print(f"  RESULT: {'PASS' if r['PASS'] else 'FAIL'}")


def main() -> int:
    ap = argparse.ArgumentParser(description="Synthetic detector-calibration validation")
    ap.add_argument("--size", type=int, default=256)
    ap.add_argument("--n", type=int, default=15, help="calibration frames per master")
    ap.add_argument("--out", type=Path, default=None, help="output dir (default: temp)")
    args = ap.parse_args()
    report = run(size=args.size, n=args.n, out_dir=args.out)
    _print(report)
    return 0 if report["PASS"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
