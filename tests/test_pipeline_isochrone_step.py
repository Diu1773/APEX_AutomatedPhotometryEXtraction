"""Step 12 runs headless, and refuses to run on defaults.

The isochrone service never needed Qt. What kept the step deferred was that the
settings deciding the answer — which colours, how wide the age window, whether
there is a reddening prior — could not be written in a config file, and library
defaults do not fail loudly: the default 0.2-6 Gyr window cannot reach a
globular, and without an E(B-V) prior an open cluster rails at the floor (both
measured, `validation/psf_crossinstrument/REPORT_UB_DEGENERACY.md`).

So the settings became config rows and the step blocks until the decisive ones
are present. These tests pin both halves: the translation from config to
`IsochroneFitConfig`, and the refusal.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from apex.config.parameters_cmd import read_params
from apex.pipeline.base import DeferredStep, StepStatus
from apex.pipeline.registry import get_steps
from apex.pipeline.steps.isochrone import (
    IsochroneStep,
    _colors,
    _prior,
    build_fit_config,
    missing_decisive_settings,
)


def _params(tmp_path, isochrone: dict | None = None):
    config = tmp_path / "apex_config.json"
    body = {"io": {"result_dir": str(tmp_path / "result"), "data_dir": str(tmp_path)}}
    if isochrone:
        body["isochrone"] = isochrone
    config.write_text(json.dumps(body), encoding="utf-8")
    return read_params(config)


def test_step_12_is_no_longer_deferred():
    step = next(s for s in get_steps("cmd") if s.index == 12)
    assert not isinstance(step, DeferredStep)
    assert isinstance(step, IsochroneStep)


def test_colours_are_read_as_pairs(tmp_path):
    params = _params(tmp_path, {"colors": "B-V, V-R"})
    assert _colors(params) == [("B", "V"), ("V", "R")]
    assert _colors(_params(tmp_path, {"colors": ""})) == []
    assert _colors(_params(tmp_path, {"colors": "nonsense"})) == []


def test_a_prior_is_value_and_sigma_or_nothing(tmp_path):
    params = _params(tmp_path, {"ecolor_prior": "0.02, 0.01"})
    assert _prior(params, "iso_ecolor_prior") == (0.02, 0.01)
    assert _prior(_params(tmp_path), "iso_ecolor_prior") is None
    assert _prior(_params(tmp_path, {"ecolor_prior": "0.02"}), "iso_ecolor_prior") is None


def test_the_config_reaches_the_fit_settings(tmp_path):
    params = _params(tmp_path, {
        "colors": "B-V",
        "mag_band": "V",
        "age_min": 9.0, "age_max": 10.2,
        "mh_prior": "-1.5,0.2",
        "n_walkers": 16, "n_steps": 200,
        "file_path": "grid.dat",
    })
    config = build_fit_config(params)
    assert config.colors == [("B", "V")]
    assert config.mag_band == "V"
    assert config.age_bounds == (9.0, 10.2)
    assert config.mh_prior == (-1.5, 0.2)
    assert (config.n_walkers, config.n_steps) == (16, 200)
    assert config.iso_file == "grid.dat"


def test_it_refuses_rather_than_guessing(tmp_path):
    """A batch that produces nothing beats one that produces a confident wrong age."""
    params = _params(tmp_path)
    missing = missing_decisive_settings(params)
    assert any("colors" in m for m in missing)
    assert any("file_path" in m for m in missing)

    zp_dir = tmp_path / "result" / "cmd_zeropoint"
    zp_dir.mkdir(parents=True)
    (zp_dir / "median_by_ID_filter_wide_cmd.csv").write_text("ID\n1\n", encoding="utf-8")

    import logging

    from apex.pipeline.context import RunContext

    ctx = RunContext(mode="cmd", params=params,
                     result_dir=Path(params.P.result_dir),
                     data_dir=Path(params.P.data_dir),
                     logger=logging.getLogger("test"))
    result = IsochroneStep().run(ctx)
    assert result.status == StepStatus.BLOCKED
    assert "colors" in result.message


def test_settings_present_means_it_would_run(tmp_path):
    params = _params(tmp_path, {"colors": "B-V", "file_path": "grid.dat"})
    assert missing_decisive_settings(params) == []


def test_the_step_module_carries_no_qt():
    import io
    import tokenize

    source = (Path(__file__).absolute().parents[1]
              / "apex/pipeline/steps/isochrone.py").read_text(encoding="utf-8")
    code = " ".join(
        tok.string
        for tok in tokenize.generate_tokens(io.StringIO(source).readline)
        if tok.type not in (tokenize.COMMENT, tokenize.STRING)
    )
    for word in ("PyQt5", "QThread", "apex.gui"):
        assert word not in code
