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


def test_choices_lead_with_auto_and_name_the_bands_each_source_covers():
    """The dropdown is built from this, so coverage has to be in the entry.

    Choosing one source is not only a change of coefficients — it drops every
    band that source has no relation for. Jordi+2010 on a Johnson B/V night
    leaves no reference at all, and the only other sign of it is one skip line
    per band in the log.
    """
    from apex.utils.gaia_transforms import (
        GAIA_TRANSFORM_AUTO, GAIA_TRANSFORM_TABLES, gaia_transform_choices,
    )

    choices = gaia_transform_choices()
    assert choices[0][0] == GAIA_TRANSFORM_AUTO
    assert set(choices[0][2]) == set(GAIA_TO_BAND)

    named = {key: bands for key, _label, bands in choices[1:]}
    assert set(named) == set(GAIA_TRANSFORM_TABLES)
    for key, bands in named.items():
        assert set(bands) == set(GAIA_TRANSFORM_TABLES[key]), key
        assert bands, f"{key} covers nothing"

    assert set(named["jordi2010"]) == {"g", "r", "i", "z"}
    assert "V" not in named["jordi2010"], (
        "a Johnson night on this source has no reference — the label must say so")


def test_the_configured_transform_source_reaches_the_parameters(tmp_path):
    """`gaia.transform_source` has to survive the two-stage loader.

    It did not at first: the key was read in the namespace constructor but had
    no row in `CMD_TOML_KEY_MAP`, so nothing ever put it into `raw` and every
    value in the file came back as "auto". `gaia.cstar_cut` was in exactly that
    state already, with a code comment telling the user to set it.
    """
    import json

    from apex.config.parameters_cmd import read_params

    for value, expected_bands in (
        ("auto", set(GAIA_TO_BAND)),
        ("jordi2010", {"g", "r", "i", "z"}),
        ("riello2021", {"B", "V", "R", "I", "U"}),
    ):
        cfg = tmp_path / f"apex_config_{value}.json"
        cfg.write_text(json.dumps({
            "io": {"data_dir": str(tmp_path), "result_dir": str(tmp_path)},
            "gaia": {"transform_source": value, "cstar_cut": True},
        }), encoding="utf-8")

        P = read_params(cfg).P
        assert getattr(P, "gaia_transform_source", None) == value
        assert set(get_gaia_to_band(P.gaia_transform_source)) == expected_bands
        assert getattr(P, "gaia_cstar_cut", None) is True, (
            "gaia.cstar_cut is documented as the user's switch and must arrive")
