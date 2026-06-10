from types import SimpleNamespace

import pytest

from apex.utils.noise_params import infer_binning, resolve_effective_noise_params


def test_infer_binning_from_common_header_keys():
    assert infer_binning({"XBINNING": 2, "YBINNING": 3}) == (2, 3)
    assert infer_binning({"CCDSUM": "4 4"}) == (4, 4)
    assert infer_binning({"BINNING": "2x2"}) == (2, 2)


def test_header_egain_wins_and_manual_rdnoise_scales_from_reference_binning():
    params = SimpleNamespace(
        gain_e_per_adu=0.012380952,
        rdnoise_e=1.28,
        binning_default=1,
        noise_use_fits_header=True,
        noise_reference_binning=1,
        noise_scale_by_binning=True,
    )
    header = {"EGAIN": 0.049523808, "XBINNING": 2, "YBINNING": 2}

    noise = resolve_effective_noise_params(params, header)

    assert noise.gain_e_per_adu == pytest.approx(0.049523808)
    assert noise.rdnoise_e == pytest.approx(2.56)
    assert noise.gain_source == "header:EGAIN"
    assert noise.rdnoise_source == "manual*sqrtbin(4/1)"


def test_measured_values_win_by_default_over_bad_fits_noise_headers():
    params = SimpleNamespace(
        gain_e_per_adu=0.689,
        rdnoise_e=2.5,
        binning_default=1,
        noise_reference_binning=1,
    )
    header = {"EGAIN": 0.049523808, "RDNOISE": 1.39, "XBINNING": 1, "YBINNING": 1}

    noise = resolve_effective_noise_params(params, header)

    assert noise.gain_e_per_adu == pytest.approx(0.689)
    assert noise.rdnoise_e == pytest.approx(2.5)
    assert noise.gain_source == "manual"
    assert noise.rdnoise_source == "manual"


def test_missing_manual_values_fall_back_to_fits_noise_headers():
    params = SimpleNamespace(
        gain_e_per_adu=None,
        rdnoise_e=None,
        binning_default=2,
        noise_use_fits_header=False,
        noise_reference_binning=1,
        noise_scale_by_binning=False,
    )
    header = {"EGAIN": 0.049523808, "RDNOISE": 1.39, "XBINNING": 2, "YBINNING": 2}

    noise = resolve_effective_noise_params(params, header)

    assert noise.gain_e_per_adu == pytest.approx(0.049523808)
    assert noise.rdnoise_e == pytest.approx(1.39)
    assert noise.gain_source == "header:EGAIN"
    assert noise.rdnoise_source == "header:RDNOISE"


def test_manual_values_are_legacy_effective_without_reference_binning():
    params = SimpleNamespace(
        gain_e_per_adu=0.1,
        rdnoise_e=1.39,
        binning_default=2,
        noise_use_fits_header=False,
        noise_reference_binning=None,
        noise_scale_by_binning=True,
    )
    header = {"XBINNING": 2, "YBINNING": 2}

    noise = resolve_effective_noise_params(params, header)

    assert noise.gain_e_per_adu == pytest.approx(0.1)
    assert noise.rdnoise_e == pytest.approx(1.39)
    assert noise.gain_source == "manual"
    assert noise.rdnoise_source == "manual"


def test_header_rdnoise_wins_without_extra_binning_scaling():
    params = SimpleNamespace(
        gain_e_per_adu=0.01,
        rdnoise_e=1.28,
        binning_default=1,
        noise_use_fits_header=True,
        noise_reference_binning=1,
        noise_scale_by_binning=True,
    )
    header = {"EGAIN": 0.04, "RDNOISE": 2.7, "XBINNING": 2, "YBINNING": 2}

    noise = resolve_effective_noise_params(params, header)

    assert noise.gain_e_per_adu == pytest.approx(0.04)
    assert noise.rdnoise_e == pytest.approx(2.7)
    assert noise.rdnoise_source == "header:RDNOISE"
