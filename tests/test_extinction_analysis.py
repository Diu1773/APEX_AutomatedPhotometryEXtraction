from __future__ import annotations

import numpy as np
import pandas as pd

from apex.analysis.extinction import (
    fit_ensemble_extinction_once,
    robust_weighted_linfit,
)


def _orthogonal_frame_offsets(airmass: np.ndarray, raw: np.ndarray) -> np.ndarray:
    dx = airmass - float(np.mean(airmass))
    offsets = raw - float(np.mean(raw))
    denom = float(np.dot(dx, dx))
    if denom > 0:
        offsets = offsets - dx * (float(np.dot(offsets, dx)) / denom)
    return offsets


def test_ensemble_extinction_recovers_synthetic_k1_with_frame_offsets():
    true_k1 = 0.22
    airmass = np.array([1.02, 1.18, 1.34, 1.57, 1.78, 1.96], dtype=float)
    raw_offsets = np.array([0.012, -0.018, 0.006, 0.020, -0.010, -0.004], dtype=float)
    frame_offsets = _orthogonal_frame_offsets(airmass, raw_offsets)
    star_mags = 12.0 + 0.17 * np.arange(10, dtype=float)

    rows = []
    for frame_idx, x in enumerate(airmass):
        for star_idx, star_mag in enumerate(star_mags):
            rows.append({
                "source_id": 1000 + star_idx,
                "file": f"frame_{frame_idx:02d}.fits",
                "airmass": x,
                "mag_inst": star_mag + true_k1 * x + frame_offsets[frame_idx],
                "snr": 250.0,
            })

    fit = fit_ensemble_extinction_once(
        pd.DataFrame(rows),
        clip_sigma=3.0,
        iters=5,
        min_n=20,
    )

    assert fit is not None
    assert abs(fit["k1"] - true_k1) < 1e-8
    assert fit["n_used"] == len(rows)
    assert fit["n_star"] == len(star_mags)
    assert fit["n_frame"] == len(airmass)
    assert float(fit["scatter"]) < 1e-10


def test_robust_weighted_linfit_rejects_large_outlier():
    x = np.linspace(1.0, 2.0, 16)
    y = 0.31 * x + 1.7
    y[-1] += 1.2
    w = np.ones_like(x)

    k, zp, scatter, n_used, outlier_frac = robust_weighted_linfit(
        x, y, w=w, clip_sigma=3.0, iters=5, min_n=8
    )

    assert abs(k - 0.31) < 1e-8
    assert abs(zp - 1.7) < 1e-8
    assert scatter < 1e-10
    assert n_used == len(x) - 1
    assert outlier_frac > 0
