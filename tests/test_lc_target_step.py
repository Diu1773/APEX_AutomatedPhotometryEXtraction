"""The LC branch stopped at Step 7, and not because the science needed a window.

Fourteen Qt-free services in `apex.analysis.light_curve` — some seven thousand
lines — already build, detrend and period-analyse. What stopped the branch is
smaller and harder than a port: nothing in the configuration could say *which
star* the light curve is of, and a batch run has nobody to click it.

Same gate the isochrone fit sat behind, same answer: config rows and a refusal
rather than a guess. A light curve of the wrong star does not look wrong — it
just belongs to something else — so picking the brightest, or the one nearest
the field centre, would be worse than producing nothing.

This suite also stands in for the guard that caught the first attempt: six
config rows were added with nothing reading them, and
`test_settings_nobody_reads` refused them. The ceiling was not raised.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from apex.analysis.light_curve.target_config import (
    missing_target_settings, read_selection, read_target, resolve_comparisons,
    write_selection,
)
from apex.pipeline.base import StepStatus
from apex.pipeline.context import RunContext
from apex.pipeline.steps.lc_target import LcTargetStep
from apex.utils.step_paths import master_catalog_path


def _params(tmp_path, **over):
    body = dict(result_dir=str(tmp_path), data_dir=str(tmp_path),
                lc_target_id=-1, lc_target_name="", lc_comparison_ids="",
                lc_comparison_mode="auto", lc_comparison_count=10, lc_filter="")
    body.update(over)
    return SimpleNamespace(P=SimpleNamespace(**body))


def _master(tmp_path, n=25):
    """Write the master catalog where Step 6 writes it, in Step 6's format.

    This helper used to invent `master_sources.csv`, which is also the name the
    step guessed — so the suite passed while the step blocked on every real
    workspace. Both now ask `master_catalog_path()`, and the binding test below
    holds that demand against what Step 6 itself calls done.
    """
    path = master_catalog_path(tmp_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"ID": range(1, n + 1)}).to_csv(path, sep="\t", index=False)


def _ctx(tmp_path, params):
    return RunContext(mode="lc", params=params, result_dir=tmp_path,
                      data_dir=tmp_path, logger=logging.getLogger("test"))


def test_it_refuses_rather_than_picking_a_star(tmp_path):
    _master(tmp_path)
    result = LcTargetStep().run(_ctx(tmp_path, _params(tmp_path)))
    assert result.status == StepStatus.BLOCKED
    assert "target_id" in result.message


def test_manual_mode_needs_the_list_it_promises(tmp_path):
    params = _params(tmp_path, lc_target_id=7, lc_comparison_mode="manual")
    missing = missing_target_settings(params)
    assert any("comparison_ids" in m for m in missing)


def test_manual_comparisons_are_taken_verbatim(tmp_path):
    """The user chose those stars; a batch must not quietly substitute others."""
    params = _params(tmp_path, lc_target_id=7, lc_comparison_mode="manual",
                     lc_comparison_ids="3, 11; 19")
    target = read_target(params)
    assert target.comparison_ids == [3, 11, 19]
    catalog = pd.DataFrame({"ID": range(1, 26)})
    assert resolve_comparisons(target, catalog) == [3, 11, 19]


def test_auto_mode_excludes_the_target_and_respects_the_count(tmp_path):
    params = _params(tmp_path, lc_target_id=7, lc_comparison_count=5)
    picked = resolve_comparisons(read_target(params), pd.DataFrame({"ID": range(1, 26)}))
    assert len(picked) == 5
    assert 7 not in picked


def test_auto_mode_prefers_the_stability_ranking_when_there_is_one(tmp_path):
    params = _params(tmp_path, lc_target_id=7, lc_comparison_count=3)
    catalog = pd.DataFrame({"ID": range(1, 26)})
    stability = pd.DataFrame({"ID": [20, 4, 15, 7], "rms": [0.01, 0.02, 0.03, 0.001]})
    picked = resolve_comparisons(read_target(params), catalog, stability)
    assert picked == [20, 4, 15], "가장 안정한 것부터, 대상은 빼고"


def test_a_target_outside_the_catalog_blocks(tmp_path):
    _master(tmp_path, n=10)
    result = LcTargetStep().run(_ctx(tmp_path, _params(tmp_path, lc_target_id=999)))
    assert result.status == StepStatus.BLOCKED
    assert "not in the master catalog" in result.message


def test_a_resolved_target_is_written_where_the_next_steps_look(tmp_path):
    _master(tmp_path)
    params = _params(tmp_path, lc_target_id=7, lc_target_name="AE UMa",
                     lc_comparison_count=4, lc_filter="V")
    result = LcTargetStep().run(_ctx(tmp_path, params))
    assert result.status == StepStatus.OK
    assert "ID 7" in result.message and "AE UMa" in result.message

    body = read_selection(tmp_path)
    assert body["target_id"] == 7
    assert body["target_name"] == "AE UMa"
    assert body["filter"] == "V"
    assert len(body["comparison_ids"]) == 4
    assert body["source"] == "config"


def test_the_gate_is_step_8_of_the_lc_pipeline():
    """This file is about the gate, not the pipeline's length.

    It used to assert the whole LC step list, and so did two other files — so
    every new LC step broke three tests that had nothing to do with it. The
    full shape now lives in `test_pipeline_runner.test_registry_shared_steps_shape`
    alone.
    """
    from apex.pipeline.registry import get_steps

    by_index = {s.index: s for s in get_steps("lc")}
    assert 8 in by_index
    assert by_index[8].key == "lctarget"
    assert isinstance(by_index[8], LcTargetStep)


def test_the_new_settings_are_actually_read():
    """The guard that caught the first attempt — six rows nothing read."""
    from apex.config.config_audit import unread_settings

    dead = set(unread_settings(mode="lc")["dead"])
    for name in ("lc_target_id", "lc_target_name", "lc_comparison_ids",
                 "lc_comparison_mode", "lc_comparison_count", "lc_filter"):
        assert name not in dead, f"{name} 을 읽는 코드가 없다"


def test_the_gate_asks_for_the_file_step_6_writes(tmp_path):
    """The input this step demands must be the file its producer calls done.

    This is the test that was missing. The step named `master_sources.csv`; the
    refbuild step writes `ref_catalog.tsv`; the fixture above invented the
    former, so nothing disagreed and `LcTargetStep` returned BLOCKED on every
    workspace in existence for as long as it existed. Comparing against a
    literal would have reproduced the same mistake, so this compares against
    what `RefBuildStep.is_complete` actually looks at.
    """
    from apex.pipeline.steps.refbuild import RefBuildStep

    required = LcTargetStep().inputs(_ctx(tmp_path, _params(tmp_path)))
    assert len(required) == 1
    demanded = required[0]

    step6 = RefBuildStep()
    ctx = _ctx(tmp_path, _params(tmp_path))
    assert not step6.is_complete(ctx)

    demanded.parent.mkdir(parents=True, exist_ok=True)
    demanded.write_text("ID\n1\n", encoding="utf-8")

    # Writing exactly what step 8 demands must be enough to make step 6 complete.
    assert step6.is_complete(ctx), (
        f"LcTargetStep demands {demanded.name}, but RefBuildStep does not "
        "consider that file its output — the two disagree about the master catalog"
    )


def test_it_reads_a_tab_separated_catalog(tmp_path):
    """Step 6 writes TSV. Read as CSV it is one column whose name is the header line."""
    path = master_catalog_path(tmp_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"ID": [1, 2, 3], "source_id": [10, 20, 30]}).to_csv(
        path, sep="\t", index=False)

    params = _params(tmp_path, lc_target_id=2, lc_comparison_mode="auto")
    result = LcTargetStep().run(_ctx(tmp_path, params))
    assert result.status == StepStatus.OK, result.message

    selection = read_selection(tmp_path)
    assert selection["target_id"] == 2
    assert selection["comparison_ids"], "auto mode found no comparisons in a TSV catalog"
