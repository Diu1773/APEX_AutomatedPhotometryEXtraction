"""The grouper's knobs must actually reach the grouper.

Two ways to lose a setting were both hit for real on 2026-08-14:

* registering a key in `parameters_cmd.TOML_KEY_MAP`, which is reassigned to
  `CMD_TOML_KEY_MAP` right after it is defined — the real map is
  `parameter_map.PSF_TOML_KEY_MAP`;
* letting the value through config but clamping it at the call site, where
  `hard_max_size=25` and `max_fraction=0.10` were written as literals. A
  configured group size of 60 became 25 and only a tenth of each frame was ever
  solved jointly, with nothing in the log saying so.

Neither shows up as an error. The run succeeds and quietly measures something
else, which is how a grouping experiment ran twice before anyone noticed the
setting had not applied. These tests pin the path from JSON to the value the
grouper receives, and pin the defaults so the fix cannot change existing runs.
"""

from __future__ import annotations

import json

import pytest

from apex.analysis.psf_policy import local_group_policy
from apex.config.parameter_map import PSF_TOML_KEY_MAP
from apex.config.parameters_cmd import read_params

GROUPER_KEYS = {
    "psf_final_pass_max_iter",
    "psf_use_grouper",
    "psf_grouper_max_size",
    "psf_grouper_radius_fwhm",
    "psf_grouper_budget_frac",
    "psf_grouper_budget_cap",
}


def test_grouper_keys_are_in_the_live_map():
    """All five knobs sit in the map step 8 actually reads."""
    mapped = {flat for _, flat in PSF_TOML_KEY_MAP}
    assert GROUPER_KEYS <= mapped, GROUPER_KEYS - mapped


def _workspace(tmp_path, psf: dict) -> object:
    cfg = tmp_path / "apex_config.json"
    cfg.write_text(json.dumps({"psf": psf}), encoding="utf-8")
    return read_params(cfg).P


@pytest.mark.parametrize("key,value", [
    ("grouper_max_size", 60),
    ("grouper_radius_fwhm", 2.5),
    ("grouper_budget_frac", 1.0),
    ("grouper_budget_cap", 0),
    # The pass that sets every published flux was pinned at two Newton steps
    # by a literal, which solves an isolated star and not a blended group.
    ("final_pass_max_iter", 30),
])
def test_configured_value_survives_the_round_trip(tmp_path, key, value):
    params = _workspace(tmp_path, {key: value})
    assert getattr(params, f"psf_{key}") == value


def test_defaults_preserve_the_previous_behaviour(tmp_path):
    """An untouched workspace must group exactly as it did before the fix."""
    params = _workspace(tmp_path, {})
    assert params.psf_grouper_budget_frac == pytest.approx(0.10)
    assert params.psf_grouper_budget_cap == 200
    assert params.psf_final_pass_max_iter == 2


def test_policy_honours_a_raised_ceiling():
    """The clamp that silently cut 60 down to 25 is gone."""
    max_size, budget = local_group_policy(
        1652, enabled=True, requested_max_size=60, hard_max_size=60,
        max_fraction=1.0, absolute_cap=0,
    )
    assert max_size == 60
    assert budget == 1652, "a fraction of 1.0 must make every source eligible"


def test_policy_still_bounds_by_default():
    """Raising the ceiling must not remove the bound for default workspaces."""
    max_size, budget = local_group_policy(
        1652, enabled=True, requested_max_size=60, hard_max_size=25,
        max_fraction=0.10, absolute_cap=200,
    )
    assert max_size == 25
    assert budget == 166
