from __future__ import annotations

import numpy as np
import pytest

from apex.utils.gaia_transforms import GAIA_TO_BAND, get_gaia_to_band


def _poly_eval(x: float, coeffs) -> float:
    return float(sum(value * x**power for power, value in enumerate(coeffs)))


def test_default_b_transform_uses_pancino_dwarf_relation():
    coeffs, color_min, color_max, source, sigma = GAIA_TO_BAND["B"]

    assert source == "Pancino+2022 dwarf"
    assert (color_min, color_max) == (-0.4, 3.5)
    assert sigma == pytest.approx(0.0248)
    assert _poly_eval(1.0, coeffs) == pytest.approx(-0.9727869442)


def test_pancino_b_restores_main_sequence_color_scale():
    b_coeffs = GAIA_TO_BAND["B"][0]
    v_coeffs = GAIA_TO_BAND["V"][0]

    bp_rp = 1.5
    g_minus_b = _poly_eval(bp_rp, b_coeffs)
    g_minus_v = _poly_eval(bp_rp, v_coeffs)
    synthetic_b_minus_v = g_minus_v - g_minus_b

    assert synthetic_b_minus_v > 1.0
    assert np.isfinite(synthetic_b_minus_v)


def test_pancino_source_can_be_selected_explicitly():
    transforms = get_gaia_to_band("pancino2022")

    assert set(transforms) == {"B"}
