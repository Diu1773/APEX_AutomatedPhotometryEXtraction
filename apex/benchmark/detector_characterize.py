"""First-principles detector characterisation from the data itself.

Measures the Moravian C3-61000 gain, read noise and dark current directly from
APEX's own bias/dark/flat frames — the FITS ``EGAIN`` keyword on these frames is
unreliable (nominal ~0.05 e-/ADU implies an absurd <4 ke- full well), so the
real values are measured, not trusted.

Gain + read noise use the mature, binning-aware Janesick photon-transfer
implementation from AstralImage (``core/camera_calib.py::measure_ptc``): read
noise from a bias-pair difference, gain from the slope of var(flat-pair diff)/2
vs signal (the flat DIFFERENCE cancels PRNU and vignette).  Our frames are 2x2
binned (XPIXSZ 7.52um = 2x IMX455 3.76um), so per-unbinned-pixel values differ
from the stored-pixel values by the binning factor — both are reported.

Dark current is the slope of source-free background vs exposure across the
20260611 10-480 s ladder.

Qt-free; run with the project venv:
    .venv-deploy/Scripts/python -m apex.benchmark.detector_characterize
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits

_ASTRAL = r"C:\Users\bmffr\Desktop\Result\Astralimage"
if _ASTRAL not in sys.path:
    sys.path.insert(0, _ASTRAL)

BIAS = r"E:\bias"
DARK_LADDER = r"E:\observe_raw_Analysis\20260611"
FLAT_DIRS = [r"E:\observe_raw_Analysis\M13_20260515\flats",
             r"E:\observe_raw_Analysis\20260530"]
BOX = 1200


def _loader(p):
    return np.asarray(fits.getdata(p), dtype=np.float64)


def _center(a, half=BOX):
    h, w = a.shape
    cy, cx = h // 2, w // 2
    return a[cy - half:cy + half, cx - half:cx + half]


def _flat_pairs(bias_level):
    """Same-level flat pairs (consecutive after sorting by level, per filter)."""
    by_filt = defaultdict(list)
    for d in FLAT_DIRS:
        for p in glob.glob(d + r"\**\*.fit", recursive=True):
            try:
                h = fits.getheader(p)
            except Exception:
                continue
            if "FLAT" in str(h.get("IMAGETYP", "")).upper():
                by_filt[str(h.get("FILTER", "?"))].append(p)
    pairs = []
    for ps in by_filt.values():
        lev = sorted((float(np.median(_center(_loader(p)))), p) for p in ps)
        for i in range(0, len(lev) - 1, 2):
            pairs.append((lev[i][1], lev[i + 1][1]))
    return pairs


def dark_current(ladder_dir, bias_level, gain):
    by_exp = defaultdict(list)
    for p in glob.glob(ladder_dir + r"\**\*.fit", recursive=True):
        try:
            h = fits.getheader(p)
        except Exception:
            continue
        if "DARK" in str(h.get("IMAGETYP", "")).upper():
            by_exp[float(h.get("EXPTIME", 0.0))].append(p)
    exps, levels = [], []
    for t in sorted(by_exp):
        stack = [_center(_loader(p)) for p in by_exp[t][:6]]
        m = np.median(np.stack(stack, 0), axis=0)
        exps.append(t)
        levels.append(float(np.median(m)) - bias_level)
    exps, levels = np.array(exps), np.array(levels)
    A = np.vstack([exps, np.ones_like(exps)]).T
    slope, intercept = np.linalg.lstsq(A, levels, rcond=None)[0]
    resid = levels - (A @ [slope, intercept])
    ss = 1 - np.sum(resid ** 2) / np.sum((levels - levels.mean()) ** 2)
    return {"dn_per_s": float(slope), "e_per_s": float(slope * gain),
            "intercept_dn": float(intercept), "linearity_r2": float(ss),
            "exptime": exps.tolist(), "level_dn": levels.tolist()}


def run() -> dict:
    from core.camera_calib import measure_ptc

    bias = sorted(glob.glob(BIAS + r"\bias-*.fit"))
    bias_level = float(np.median(_center(_loader(bias[0]))))
    pairs = _flat_pairs(bias_level)
    r = measure_ptc(
        pairs, bias_pair=(bias[0], bias[1]), bias_level=bias_level,
        box=None, loader=_loader,
    )
    gain_eff = float(r.gain_eff)
    dark = dark_current(DARK_LADDER, bias_level, gain_eff)

    return {
        "camera": "Moravian C3-61000 (Sony IMX455)",
        "binning": int(getattr(r, "binning", 2) or 2),
        "bias_level_dn": bias_level,
        "gain_eff_e_per_adu": gain_eff,
        "read_noise_eff_e": float(r.read_noise_eff),
        "gain_pixel_e_per_adu": float(r.gain_pixel),
        "read_noise_pixel_e": float(r.read_noise_pixel),
        "header_egain": (float(r.header_egain) if r.header_egain is not None else None),
        "header_ratio_measured_over_egain": (float(r.header_ratio) if r.header_ratio is not None else None),
        "n_ptc_points": int(np.asarray(r.signal_adu).size),
        "dark_current_e_per_s": dark["e_per_s"],
        "dark_current_dn_per_s": dark["dn_per_s"],
        "dark_linearity_r2": dark["linearity_r2"],
        "ptc": {"signal_adu": np.asarray(r.signal_adu).tolist(),
                "variance_adu2": np.asarray(r.variance_adu2).tolist(),
                "slope": float(r.fit_slope), "intercept": float(r.fit_intercept)},
        "dark": dark,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="Detector characterisation (gain/RN/dark) from data")
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).resolve().parents[2] / "validation" / "paper" / "data")
    args = ap.parse_args()
    r = run()
    print("Detector characterisation — Moravian C3-61000 (measured from data)")
    print("-" * 62)
    print(f"  binning              : {r['binning']}x{r['binning']}")
    print(f"  bias level           : {r['bias_level_dn']:.1f} DN")
    print(f"  gain (stored pixel)  : {r['gain_eff_e_per_adu']:.4f} e-/ADU")
    print(f"  read noise (stored)  : {r['read_noise_eff_e']:.2f} e-")
    print(f"  gain (unbinned pixel): {r['gain_pixel_e_per_adu']:.4f} e-/ADU")
    print(f"  read noise (unbinned): {r['read_noise_pixel_e']:.2f} e-")
    print(f"  header EGAIN         : {r['header_egain']}  (measured/header = {r['header_ratio_measured_over_egain']})")
    print(f"  dark current         : {r['dark_current_dn_per_s']:.4f} DN/s = {r['dark_current_e_per_s']:.4f} e-/s")
    print(f"  dark linearity R^2   : {r['dark_linearity_r2']:.5f}")
    print(f"  PTC points           : {r['n_ptc_points']}")
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "detector_characterization.json").write_text(json.dumps(r, indent=2))
    print(f"  wrote {args.out / 'detector_characterization.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
