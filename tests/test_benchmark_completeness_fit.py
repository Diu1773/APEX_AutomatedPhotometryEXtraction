import numpy as np
import pandas as pd

from apex.benchmark.metrics import (
    fit_completeness_logistic,
    magnitude_point_summary,
)


def test_logistic_completeness_fit_recovers_known_m50():
    rng = np.random.default_rng(123)
    magnitudes = np.repeat(np.arange(14.0, 15.61, 0.2), 300)
    m50_true = 14.8
    width_true = 0.18
    probability = 1.0 / (1.0 + np.exp((magnitudes - m50_true) / width_true))
    recovered = rng.random(len(magnitudes)) < probability
    stars = pd.DataFrame(
        {
            "magnitude_true": magnitudes,
            "recovered": recovered,
            "baseline_confounded": False,
            "trial": np.tile(np.arange(30), len(magnitudes) // 30),
            "forced_mag_error": 0.0,
        }
    )

    fit = fit_completeness_logistic(stars, bootstrap_samples=100, seed=7)

    assert abs(fit["m50"] - m50_true) < 0.05
    assert abs(fit["width_mag"] - width_true) < 0.05
    assert fit["m50_ci95_low"] < fit["m50"] < fit["m50_ci95_high"]


def test_magnitude_point_summary_has_wilson_interval():
    stars = pd.DataFrame(
        {
            "magnitude_true": [15.0] * 10,
            "recovered": [True] * 7 + [False] * 3,
            "baseline_confounded": False,
            "forced_mag_error": np.linspace(-0.02, 0.02, 10),
        }
    )

    summary = magnitude_point_summary(stars).iloc[0]

    assert summary["completeness"] == 0.7
    assert summary["ci95_low"] < 0.7 < summary["ci95_high"]
