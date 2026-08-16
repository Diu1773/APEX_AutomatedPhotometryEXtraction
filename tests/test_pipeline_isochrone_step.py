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


def test_it_reads_the_key_the_service_actually_writes():
    """`convergence_ok`, not `converged`.

    The first version of this step asked the summary for `converged`. The
    service writes `convergence_ok`, so the missing key came back None and
    every fit — including the ones that converged — was reported as a failure.
    A batch that mislabels its own success is worse than one that is silent.
    """
    import inspect

    from apex.pipeline.steps import isochrone as mod

    source = inspect.getsource(mod.IsochroneStep.run)
    assert 'summary.get("convergence_ok")' in source
    assert 'summary.get("converged")' not in source

    service = inspect.getsource(
        __import__("apex.analysis.cmd.isochrone_fit_service", fromlist=["x"]))
    assert '"convergence_ok"' in service, "서비스가 쓰는 키 이름이 바뀌었다"


def test_the_record_says_what_the_fit_was(tmp_path, monkeypatch):
    """A posterior without its bounds and priors cannot be judged."""
    import json
    import logging
    from types import SimpleNamespace

    from apex.pipeline.context import RunContext

    params = _params(tmp_path, {"colors": "B-V", "file_path": "grid.dat",
                                "n_walkers": 8, "n_steps": 10})
    zp_dir = tmp_path / "result" / "cmd_zeropoint"
    zp_dir.mkdir(parents=True)
    (zp_dir / "median_by_ID_filter_wide_cmd.csv").write_text("ID\n1\n", encoding="utf-8")

    def fake_fit(df, config, make_figures=True, progress_cb=None):
        return SimpleNamespace(
            summary={"convergence_ok": True, "age_gyr": [1, 2, 3]},
            n_stars=42, member_meta={"applied": False}, warnings=["조심"],
        )

    import apex.analysis.cmd.isochrone_fit_service as service
    monkeypatch.setattr(service, "fit_cluster_isochrone", fake_fit)

    ctx = RunContext(mode="cmd", params=params,
                     result_dir=Path(params.P.result_dir),
                     data_dir=Path(params.P.data_dir),
                     logger=logging.getLogger("test"))
    result = IsochroneStep().run(ctx)
    assert result.status == StepStatus.OK
    assert "did NOT" not in result.message

    written = json.loads(
        (Path(params.P.result_dir) / "cmd_isochrone"
         / "isochrone_fit_summary.json").read_text(encoding="utf-8"))
    assert written["n_stars"] == 42
    assert written["warnings"] == ["조심"]
    assert written["settings"]["n_walkers"] == 8
    assert written["settings"]["seed"] is not None
    assert written["settings"]["age_bounds"]
    assert written["wide_table"].endswith("median_by_ID_filter_wide_cmd.csv")


# ---------------------------------------------------------------------------
# The search box: where the walls are decides what can be found
# ---------------------------------------------------------------------------

def test_the_metallicity_box_can_hold_a_globular(tmp_path):
    """The default `[M/H]` floor used to be -1.0. M13 is near -1.5.

    Worse than "cannot reach it": the fit narrows this box around an [M/H]
    prior with max(lo,·)/min(hi,·) and never widens it, so a metal-poor prior
    inside a (-1.0, +0.5) box produced lo = -1.000, hi = -1.360 — inverted. The
    grid mask `(mh >= lo) & (mh <= hi)` then selects nothing, and the failure
    reads as bad data rather than as walls in the wrong place.
    """
    config = build_fit_config(_params(tmp_path, {"colors": "B-V", "file_path": "g.dat"}))
    lo, hi = config.mh_bounds
    assert lo <= -1.5, f"[M/H] 하한 {lo} 는 구상성단을 못 담는다"
    assert hi >= 0.0
    assert config.ecolor_bounds[1] >= 1.0, "E(colour) 상한이 데스크톱 판보다 좁다"


def test_an_out_of_box_prior_is_named_before_the_mcmc(tmp_path):
    from apex.pipeline.steps.isochrone import check_bounds

    ok = build_fit_config(_params(tmp_path, {
        "colors": "B-V", "file_path": "g.dat", "mh_prior": "-1.56, 0.10"}))
    assert check_bounds(ok) == [], "이 사전값은 이제 상자 안에 있다"

    narrow = build_fit_config(_params(tmp_path, {
        "colors": "B-V", "file_path": "g.dat",
        "mh_min": -1.0, "mh_max": 0.5, "mh_prior": "-1.56, 0.10"}))
    problems = check_bounds(narrow)
    assert problems, "뒤집힌 상자를 잡아내지 못한다"
    assert "[M/H]" in problems[0]


def test_it_blocks_instead_of_spending_an_mcmc_on_an_empty_grid(tmp_path):
    import logging

    from apex.pipeline.context import RunContext

    params = _params(tmp_path, {
        "colors": "B-V", "file_path": "g.dat",
        "mh_min": -1.0, "mh_max": 0.5, "mh_prior": "-1.56, 0.10"})
    zp_dir = tmp_path / "result" / "cmd_zeropoint"
    zp_dir.mkdir(parents=True)
    (zp_dir / "median_by_ID_filter_wide_cmd.csv").write_text("ID\n1\n", encoding="utf-8")

    ctx = RunContext(mode="cmd", params=params,
                     result_dir=Path(params.P.result_dir),
                     data_dir=Path(params.P.data_dir),
                     logger=logging.getLogger("test"))
    result = IsochroneStep().run(ctx)
    assert result.status == StepStatus.BLOCKED
    assert "inverted" in result.message


def test_the_desktop_dialog_no_longer_decides_the_boxes_itself():
    """Same setting, one place — the point of the whole 2026-08-17 sweep.

    The dialog hardcoded `mh_bounds=(-1.0, 0.5)`, `dm_bounds` and
    `ecolor_bounds` while the headless step read them from the configuration.
    Two authorities for one number is how they came to disagree.
    """
    import inspect
    import re

    source = (Path(__file__).absolute().parents[1]
              / "apex/gui/workflow/cmd/step12_isochrone_model.py").read_text(encoding="utf-8")
    body = re.sub(r"#[^\n]*", "", source)
    for name in ("mh_bounds=", "dm_bounds=", "ecolor_bounds="):
        assert name not in body, f"창이 아직 {name} 를 스스로 정한다"
    assert "build_fit_config(self.params" in body, "창이 공용 빌더를 안 쓴다"
    assert "check_bounds(cfg)" in body, "창이 뒤집힌 상자를 안 막는다"

    # And both callers must be translating with the same function.
    from apex.analysis.cmd import isochrone_config
    from apex.pipeline.steps import isochrone as step
    assert step.build_fit_config is isochrone_config.build_fit_config
    assert inspect.getmodule(step.build_fit_config) is isochrone_config


def test_the_result_survives_its_own_output_directory_vanishing(tmp_path, monkeypatch):
    """A 40-minute posterior must not be lost to a missing folder.

    Measured 2026-08-17: a 32 x 6000 M13 chain ran to completion and then died
    with FileNotFoundError writing the summary. The directory was created before
    the fit and was gone by the end — an empty folder under a temp root does not
    reliably survive three quarters of an hour. Create it again at write time.
    """
    import logging
    import shutil
    from types import SimpleNamespace

    from apex.pipeline.context import RunContext

    params = _params(tmp_path, {"colors": "B-V", "file_path": "grid.dat"})
    zp_dir = tmp_path / "result" / "cmd_zeropoint"
    zp_dir.mkdir(parents=True)
    (zp_dir / "median_by_ID_filter_wide_cmd.csv").write_text("ID\n1\n", encoding="utf-8")

    out_dir = Path(params.P.result_dir) / "cmd_isochrone"

    def fit_then_lose_the_folder(df, config, make_figures=True, progress_cb=None):
        # Whatever the fit was doing, the folder is gone when it returns.
        shutil.rmtree(out_dir, ignore_errors=True)
        return SimpleNamespace(summary={"convergence_ok": True}, n_stars=7,
                               member_meta={}, warnings=[])

    import apex.analysis.cmd.isochrone_fit_service as service
    monkeypatch.setattr(service, "fit_cluster_isochrone", fit_then_lose_the_folder)

    ctx = RunContext(mode="cmd", params=params,
                     result_dir=Path(params.P.result_dir),
                     data_dir=Path(params.P.data_dir),
                     logger=logging.getLogger("test"))
    result = IsochroneStep().run(ctx)
    assert result.status == StepStatus.OK
    assert (out_dir / "isochrone_fit_summary.json").exists(), "결과를 잃었다"
