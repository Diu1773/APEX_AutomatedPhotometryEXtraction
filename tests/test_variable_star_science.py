from __future__ import annotations

import numpy as np

from apex.gui.tools.variable_star import (
    _detect_phase_epoch,
    _fit_fixed_period_fourier,
    _fourier_shape_parameters,
)


def _phase_distance(value: float, expected: float) -> float:
    return abs((value - expected + 0.5) % 1.0 - 0.5)


def _angle_distance(value: float, expected: float) -> float:
    return abs((value - expected + np.pi) % (2.0 * np.pi) - np.pi)


def test_weighted_fourier_fit_uses_documented_cosine_phase_convention():
    period = 0.104
    time = 2459000.0 + np.linspace(0.0, 2.1, 500)
    tau = time - np.min(time)
    omega = 2.0 * np.pi / period
    phi1 = 0.43
    phi2 = -0.71
    mag = (
        10.2
        + 0.18 * np.cos(omega * tau + phi1)
        + 0.05 * np.cos(2.0 * omega * tau + phi2)
    )
    mag_err = np.linspace(0.008, 0.025, len(time))

    fit = _fit_fixed_period_fourier(
        time,
        mag,
        period,
        harmonics=2,
        mag_err=mag_err,
    )
    shape = _fourier_shape_parameters(fit["coeff"])

    assert fit["weighted"] is True
    assert np.isclose(shape["amplitudes"][0], 0.18, atol=1e-8)
    assert np.isclose(shape["amplitudes"][1], 0.05, atol=1e-8)
    assert _angle_distance(shape["phases"][0], phi1) < 1e-8
    assert _angle_distance(shape["phases"][1], phi2) < 1e-8


def test_epoch_detection_distinguishes_maximum_and_minimum_light():
    period = 0.104
    reference_maximum = 2459000.031
    time = 2458999.97 + np.linspace(0.0, 1.1, 1200)
    phase = (time - reference_maximum) / period
    mag = 10.0 - 0.22 * np.cos(2.0 * np.pi * phase)

    maximum_epoch, _ = _detect_phase_epoch(
        time,
        mag,
        period,
        event_kind="max_light",
    )
    minimum_epoch, _ = _detect_phase_epoch(
        time,
        mag,
        period,
        event_kind="min_light",
    )

    maximum_phase = ((maximum_epoch - reference_maximum) / period) % 1.0
    minimum_phase = ((minimum_epoch - reference_maximum) / period) % 1.0
    assert _phase_distance(maximum_phase, 0.0) < 0.02
    assert _phase_distance(minimum_phase, 0.5) < 0.02
