"""Synthetic truth validation for CMD Step12 isochrone fitting.

This runner builds a small deterministic isochrone grid, generates observed CMD
points from known parameters, runs the Step12 fitter, and checks recovery
tolerances. It is intended as a science-regression guard, not as an astrophysical
calibration dataset.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).absolute().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from apex.analysis.cmd.isochrone_fitter_v2 import FitBounds, IsochroneFitterV2


TRUTH = {
    "log_age": 8.50,
    "metallicity": 0.00,
    "distance_mod": 10.00,
    "extinction_color": 0.24,
}

TOLERANCES = {
    "log_age": 0.04,
    "metallicity": 0.06,
    "distance_mod": 0.10,
    "extinction_color": 0.04,
}


def build_synthetic_iso_data() -> np.ndarray:
    """Return PARSEC-like columns: [unused, mh, log_age, band1, band2, mass]."""
    rows: list[np.ndarray] = []
    masses = np.linspace(0.28, 1.85, 120)
    for log_age in (8.00, 8.25, 8.50, 8.75, 9.00):
        for mh in (-0.30, 0.00, 0.20):
            turnoff = 1.45 - 0.48 * (log_age - 8.5) - 0.10 * mh
            evolved = np.clip(masses - turnoff, 0.0, None)

            color = (
                0.22
                + 0.42 * (1.85 - masses)
                + 0.32 * (log_age - 8.5)
                + 0.10 * mh
                + 0.18 * mh * (1.85 - masses)
                + 0.28 * evolved**2
            )
            mag = (
                2.20
                + 4.30 * (1.85 - masses)
                + 1.15 * (log_age - 8.5)
                + 0.25 * mh
                + 0.45 * mh * (1.85 - masses) ** 2
                - 1.35 * evolved
                + 0.70 * evolved**2
            )

            band1 = mag
            band2 = mag - color
            for m, b1, b2 in zip(masses, band1, band2):
                row = np.zeros(6, dtype=float)
                row[1] = mh
                row[2] = log_age
                row[3] = b1
                row[4] = b2
                row[5] = m
                rows.append(row)
    return np.asarray(rows, dtype=float)


def make_fitter(iso_data: np.ndarray) -> IsochroneFitterV2:
    return IsochroneFitterV2(
        ROOT / "validation" / "synthetic_isochrone.dat",
        col_mh=1,
        col_age=2,
        col_g=3,
        col_r=4,
        col_mag=3,
        col_mass=5,
        fit_fraction=1.00,
        R_band1=3.303,
        R_band2=2.285,
        R_mag=3.303,
        iso_data=iso_data,
        use_bilinear=False,
    )


def build_observed_cmd(seed: int = 42) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    iso_data = build_synthetic_iso_data()
    truth_fitter = make_fitter(iso_data)
    color, mag, _ = truth_fitter._interpolate_isochrone(
        TRUTH["log_age"],
        TRUTH["metallicity"],
        TRUTH["distance_mod"],
        TRUTH["extinction_color"],
    )

    rng = np.random.default_rng(seed)
    keep = np.linspace(8, len(color) - 12, 85, dtype=int)
    obs_color = color[keep] + rng.normal(0.0, 0.007, size=len(keep))
    obs_mag = mag[keep] + rng.normal(0.0, 0.010, size=len(keep))

    color_err = np.full(len(obs_color), 0.025, dtype=float)
    mag_err = np.full(len(obs_mag), 0.030, dtype=float)
    return obs_color, obs_mag, color_err, mag_err


def run_validation(seed: int = 42) -> dict[str, Any]:
    iso_data = build_synthetic_iso_data()
    fitter = make_fitter(iso_data)
    obs_color, obs_mag, color_err, mag_err = build_observed_cmd(seed=seed)

    bounds = FitBounds(
        log_age=(8.25, 8.75),
        metallicity=(-0.30, 0.20),
        distance_mod=(9.70, 10.30),
        extinction_gr=(0.10, 0.38),
    )

    start = time.perf_counter()
    grid = fitter.fit_grid_scan(
        obs_color,
        obs_mag,
        color_err,
        mag_err,
        bounds=bounds,
        snr_min=5.0,
        initial_dm=9.95,
        initial_egr=0.22,
        local_maxiter=80,
    )
    elapsed = time.perf_counter() - start
    result = grid.best_fit

    recovered = {
        "log_age": float(result.log_age),
        "metallicity": float(result.metallicity),
        "distance_mod": float(result.distance_mod),
        "extinction_color": float(result.extinction_gr),
    }
    deltas = {key: abs(recovered[key] - TRUTH[key]) for key in TRUTH}
    checks = {
        key: {
            "truth": TRUTH[key],
            "recovered": recovered[key],
            "abs_delta": deltas[key],
            "tolerance": TOLERANCES[key],
            "passed": deltas[key] <= TOLERANCES[key],
        }
        for key in TRUTH
    }
    passed = all(item["passed"] for item in checks.values()) and bool(result.converged)

    return {
        "scenario": "cmd_step12_synthetic_truth",
        "passed": bool(passed),
        "seed": seed,
        "n_stars": int(result.n_stars),
        "n_grid_cells_evaluated": int(grid.n_evaluated),
        "elapsed_sec": float(elapsed),
        "fit": {
            "converged": bool(result.converged),
            "chi2": float(result.chi2),
            "reduced_chi2": float(result.reduced_chi2),
            "fit_mode": result.fit_mode,
        },
        "checks": checks,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", type=Path, help="Optional JSON report path.")
    args = parser.parse_args(argv)

    report = run_validation(seed=args.seed)
    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")

    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
