"""The profile-error term stays off by default, and the reason is measured.

DAOPHOT adds `proferr` (5 % of the counts) to the fit variance, so a bright
core stops dominating the solution. Measured on M13 alone it looked like a
clear win and nearly became the default. Widening to three instruments on
2026-08-15 reversed the sign on one of them: on the QHY600 field every star
fainter than the bright anchor came back ~20 mmag too bright, because a core
that no longer holds its own flux lets the fainter neighbours absorb it.

So the value is an option for crowded fields, not a default — and the thing
worth pinning is that turning it off reproduces the old arithmetic exactly.
A default that quietly moved would invalidate every product measured before it.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from apex.config.parameter_map import PSF_TOML_KEY_MAP
from apex.config.parameters_cmd import read_params
from apex.gui.workflow.cmd.step8_psf_photometry import _fit_variance


def test_the_key_reaches_the_runtime_parameters(tmp_path):
    # A map row is (path, attr) or (path, attr, kind, default) since
    # 2026-08-16, so it no longer converts to a dict.
    assert any(row[0] == ("psf", "profile_error_frac") for row in PSF_TOML_KEY_MAP)
    config = tmp_path / "apex_config.json"
    config.write_text(json.dumps({"psf": {"profile_error_frac": 0.05}}),
                      encoding="utf-8")
    assert read_params(config).P.psf_profile_error_frac == pytest.approx(0.05)


def test_the_default_is_off(tmp_path):
    config = tmp_path / "apex_config.json"
    config.write_text("{}", encoding="utf-8")
    assert read_params(config).P.psf_profile_error_frac == 0.0


def test_off_is_photon_and_background_noise_and_nothing_else():
    counts = np.array([0.0, 50.0, 500.0, 20000.0])
    variance = _fit_variance(counts, background_rms=18.08, gain=1.0,
                             profile_error_frac=0.0)
    assert variance == pytest.approx(18.08 ** 2 + counts)


def test_on_adds_a_term_that_grows_with_the_counts_themselves():
    """Quadratic in counts, so it only bites where the star is bright."""
    counts = np.array([50.0, 20000.0])   # sky-subtracted ADU
    off = _fit_variance(counts, background_rms=18.08, gain=1.0,
                        profile_error_frac=0.0)
    on = _fit_variance(counts, background_rms=18.08, gain=1.0,
                       profile_error_frac=0.05)
    assert on == pytest.approx(off + (0.05 * counts) ** 2)
    faint, bright = (on - off) / on
    assert faint < 0.02 < 0.9 < bright


def test_negative_counts_do_not_reduce_the_variance():
    """A sky-subtracted pixel can go negative; the term must not cancel noise."""
    variance = _fit_variance(np.array([-400.0]), background_rms=10.0, gain=1.0,
                             profile_error_frac=0.05)
    assert variance == pytest.approx([100.0])
