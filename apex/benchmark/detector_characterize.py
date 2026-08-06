"""Detector characterisation from the calibration frames themselves.

Measures the gain, read noise and dark current of the camera used for the
paper's photometry directly from its own bias, dark and flat exposures.  The
FITS ``EGAIN`` keyword on this camera (Moravian C3-61000, Sony IMX455) is the
sensor's nominal register value rather than the realised conversion gain — it
is wrong by more than a factor of ten — and the manufacturer publishes only a
read-noise ceiling and a full-well figure, never a gain.  Measuring is
therefore the only option, not a preference.

Gain and read noise come from :mod:`apex.analysis.detector_ptc`, the same code
the GUI tool and the pipeline's error model use.  Dark current is the slope of
the source-free background against exposure time across a dark ladder.

Qt-free.  Run with the project venv::

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

from apex.analysis.detector_ptc import characterize_detector, read_box, robust_center

# Defaults describe the paper's KNUEMAO data set; override on the command line
# for any other camera.
DEFAULT_BIAS = r"E:\bias"
DEFAULT_DARK_LADDER = r"E:\observe_raw_Analysis\20260611"
DEFAULT_FLAT_DIRS = (
    r"E:\observe_raw_Analysis\M13_20260515\flats",
    r"E:\observe_raw_Analysis\20260530",
)
BOX_RADIUS = 300
SIGNAL_FLOOR = 20000.0


def _find(root: str, pattern: str = "*.fit") -> list[str]:
    return sorted(glob.glob(str(Path(root) / "**" / pattern), recursive=True))


def _of_type(paths, keyword: str) -> list[str]:
    out = []
    for p in paths:
        try:
            imagetyp = str(fits.getheader(p).get("IMAGETYP", "")).upper()
        except Exception:  # noqa: BLE001 — an unreadable frame is simply skipped
            continue
        if keyword in imagetyp:
            out.append(p)
    return out


def dark_current(ladder_dir: str, bias_level: float, gain: float,
                 box_radius: int = BOX_RADIUS) -> dict:
    """Dark current from the slope of background level against exposure time."""
    by_exposure: dict[float, list[str]] = defaultdict(list)
    for p in _find(ladder_dir):
        try:
            header = fits.getheader(p)
        except Exception:  # noqa: BLE001
            continue
        if "DARK" in str(header.get("IMAGETYP", "")).upper():
            by_exposure[float(header.get("EXPTIME", 0.0))].append(p)

    exposures, levels = [], []
    for t in sorted(by_exposure):
        stack = [read_box(p, box_radius) for p in by_exposure[t][:6]]
        median_frame = np.median(np.stack(stack, axis=0), axis=0)
        exposures.append(t)
        levels.append(robust_center(median_frame) - bias_level)

    exposures, levels = np.array(exposures), np.array(levels)
    design = np.vstack([exposures, np.ones_like(exposures)]).T
    slope, intercept = np.linalg.lstsq(design, levels, rcond=None)[0]
    residual = levels - design @ np.array([slope, intercept])
    ss_tot = float(np.sum((levels - levels.mean()) ** 2))
    r_squared = 1.0 - float(residual @ residual) / ss_tot if ss_tot > 0 else float("nan")

    return {
        "dn_per_s": float(slope),
        "e_per_s": float(slope * gain),
        "intercept_dn": float(intercept),
        "linearity_r2": float(r_squared),
        "exptime": exposures.tolist(),
        "level_dn": levels.tolist(),
    }


def run(bias_dir: str = DEFAULT_BIAS,
        flat_dirs: tuple[str, ...] = DEFAULT_FLAT_DIRS,
        dark_ladder: str = DEFAULT_DARK_LADDER,
        signal_floor: float = SIGNAL_FLOOR,
        verbose: bool = True) -> dict:
    log = (lambda m: print(m, flush=True)) if verbose else (lambda m: None)

    bias_paths = _find(bias_dir, "bias-*.fit")
    flat_paths: list[str] = []
    for d in flat_dirs:
        flat_paths.extend(_of_type(_find(d), "FLAT"))
    log(f"bias {len(bias_paths)} frames, flat {len(flat_paths)} frames")

    ptc = characterize_detector(
        bias_paths, flat_paths,
        box_radius=BOX_RADIUS, signal_floor=signal_floor,
        fix_intercept=True, log_fn=log,
    )
    log("\n" + ptc.summary())

    bias_level = robust_center(read_box(bias_paths[0], BOX_RADIUS))
    dark = dark_current(dark_ladder, bias_level, ptc.gain_eff)
    log(f"dark current    : {dark['e_per_s']:.5f} e-/s "
        f"(R^2 = {dark['linearity_r2']:.4f})")

    return {
        "camera": "Moravian C3-61000 (Sony IMX455)",
        "binning": ptc.binning,
        "bias_level_dn": bias_level,
        "signal_floor_adu": signal_floor,
        "gain_eff_e_per_adu": ptc.gain_eff,
        "gain_eff_err": ptc.gain_eff_err,
        "read_noise_eff_e": ptc.read_noise_eff,
        "read_noise_adu": ptc.read_noise_adu,
        "gain_pixel_e_per_adu": ptc.gain_pixel,
        "read_noise_pixel_e": ptc.read_noise_pixel,
        "header_egain": ptc.header_egain,
        "header_ratio_measured_over_egain": ptc.header_ratio,
        "n_ptc_points": ptc.n_pairs,
        "signal_min_adu": ptc.signal_min,
        "signal_max_adu": ptc.signal_max,
        "fit_r_squared": ptc.r_squared,
        "fit_intercept_adu2": ptc.fit_intercept,
        "max_residual_frac": ptc.max_residual_frac,
        "dark_current_e_per_s": dark["e_per_s"],
        "dark_current_dn_per_s": dark["dn_per_s"],
        "dark_linearity_r2": dark["linearity_r2"],
        "ptc": {
            "signal_adu": [p.signal_adu for p in ptc.points],
            "variance_adu2": [p.variance_adu2 for p in ptc.points],
            "filter": [p.filter_name for p in ptc.points],
            "slope": ptc.fit_slope,
            "intercept": ptc.fit_intercept,
        },
        "dark": dark,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Detector characterisation (gain / read noise / dark) from data")
    parser.add_argument("--bias", default=DEFAULT_BIAS)
    parser.add_argument("--flats", nargs="+", default=list(DEFAULT_FLAT_DIRS))
    parser.add_argument("--darks", default=DEFAULT_DARK_LADDER)
    parser.add_argument("--signal-floor", type=float, default=SIGNAL_FLOOR)
    parser.add_argument(
        "--out", type=Path,
        default=Path(__file__).resolve().parents[2] / "validation" / "paper" / "data")
    args = parser.parse_args(argv)

    r = run(args.bias, tuple(args.flats), args.darks, args.signal_floor)

    print("\nDetector characterisation — Moravian C3-61000 (measured from data)")
    print("-" * 66)
    print(f"  binning              : {r['binning']}x{r['binning']}")
    print(f"  bias level           : {r['bias_level_dn']:.1f} DN")
    print(f"  gain (stored pixel)  : {r['gain_eff_e_per_adu']:.4f} "
          f"+- {r['gain_eff_err']:.4f} e-/ADU")
    print(f"  read noise (stored)  : {r['read_noise_eff_e']:.3f} e-")
    print(f"  gain (unbinned pixel): {r['gain_pixel_e_per_adu']:.4f} e-/ADU")
    print(f"  read noise (unbinned): {r['read_noise_pixel_e']:.3f} e-")
    print(f"  header EGAIN         : {r['header_egain']}  "
          f"(measured/header = {r['header_ratio_measured_over_egain']:.2f}x)")
    print(f"  dark current         : {r['dark_current_dn_per_s']:.5f} DN/s "
          f"= {r['dark_current_e_per_s']:.5f} e-/s")
    print(f"  dark linearity R^2   : {r['dark_linearity_r2']:.5f}")
    print(f"  PTC points           : {r['n_ptc_points']}  "
          f"({r['signal_min_adu']:,.0f}–{r['signal_max_adu']:,.0f} ADU)")

    args.out.mkdir(parents=True, exist_ok=True)
    out_file = args.out / "detector_characterization.json"
    out_file.write_text(json.dumps(r, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"  wrote {out_file}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
