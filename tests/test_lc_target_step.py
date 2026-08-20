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


# ── the screen the batch run was missing ───────────────────────────────────
#
# `auto` mode ranked by *catalogue order* for its first eight days, because the
# stability screen that ranks by measured steadiness lived inside the window. On
# 364 YZ Boo frames the two ensembles' averages came out 0.68 mag apart and
# nothing in either output said the runs had used different stars.

def _photometry(tmp_path, filters=("g",), n_frames=40, n_stars=14):
    """Forced-photometry output shaped the way the loader actually reads it.

    The first version of this fixture invented `<stem>_phot.csv`; the loader
    looks for `photometry_<frame>.tsv`, tab-separated. That is the same class of
    mistake as the `master_sources.csv` one — a fixture agreeing with a guess
    instead of with the producer — so the filename comes from the loader's own
    resolver rather than from memory.
    """
    import numpy as np

    from apex.utils.step_paths import step7_forced_phot_dir

    out = Path(step7_forced_phot_dir(tmp_path))
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(11)
    index_rows = []
    for flt in filters:
        for frame in range(n_frames):
            name = f"{flt}_{frame:04d}.fit"
            # Scatter grows with star index, so the stability ranking has a real
            # order to find rather than ties broken by row position.
            table = pd.DataFrame([{
                "source_id": 1000 + star,
                "mag": 12.0 + 0.2 * star + rng.normal(0, 0.002 * star),
                "mag_err": 0.004,
            } for star in range(1, n_stars + 1)])
            table.to_csv(out / f"photometry_{name}.tsv", sep="	", index=False)
            index_rows.append({"file": name, "filter": flt, "JD": 2460000.0 + frame})
    pd.DataFrame(index_rows).to_csv(out / "photometry_index.csv", index=False)


def _master_with_source_ids(tmp_path, n=14):
    path = master_catalog_path(tmp_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "ID": range(1, n + 1),
        "source_id": [1000 + i for i in range(1, n + 1)],
        "phot_bp_mean_mag": [15.0 + 0.1 * i for i in range(1, n + 1)],
        "phot_rp_mean_mag": [14.0 + 0.1 * i for i in range(1, n + 1)],
    }).to_csv(path, sep="\t", index=False)


def test_filter_all_is_a_sentinel_not_a_filter_name(tmp_path):
    """`lightcurve.filter=all` must not be handed to the photometry loader.

    The window writes `all` to mean "these roles apply to every filter". Asking
    for a filter literally named `all` matched nothing, and the step reported
    "no usable photometry" about a workspace holding hundreds of frames.
    """
    from apex.pipeline.steps.lc_target import _screen
    from apex.analysis.light_curve.target_config import LcTarget

    _master_with_source_ids(tmp_path)
    _photometry(tmp_path)
    catalog = pd.read_csv(master_catalog_path(tmp_path), sep="\t")
    ctx = _ctx(tmp_path, _params(tmp_path, lc_target_id=1, lc_filter="all"))

    result, note = _screen(ctx, LcTarget(target_id=1, filter_key="all"),
                           catalog, 1001)
    assert result is not None, note
    assert "most frames" in note
    assert "all" not in note.split("(")[0].split()


def test_a_named_filter_with_no_photometry_says_so(tmp_path):
    from apex.pipeline.steps.lc_target import _screen
    from apex.analysis.light_curve.target_config import LcTarget

    _master_with_source_ids(tmp_path)
    _photometry(tmp_path, filters=("g",))
    catalog = pd.read_csv(master_catalog_path(tmp_path), sep="\t")
    ctx = _ctx(tmp_path, _params(tmp_path, lc_target_id=1, lc_filter="V"))

    result, note = _screen(ctx, LcTarget(target_id=1, filter_key="V"),
                           catalog, 1001)
    assert result is None
    assert "no photometry" in note and "g" in note


def test_auto_mode_ranks_by_stability_and_says_which_screen_ran(tmp_path):
    _master_with_source_ids(tmp_path)
    _photometry(tmp_path)
    ctx = _ctx(tmp_path, _params(tmp_path, lc_target_id=1,
                                 lc_comparison_mode="auto", lc_filter=""))

    out = LcTargetStep().run(ctx)
    assert out.status is StepStatus.OK, out.message
    assert "measured" in out.message and "adopted" in out.message

    written = json.loads(
        read_selection_path(tmp_path).read_text(encoding="utf-8"))
    assert written["selected_by"] == "stability"
    assert written["check_id"] not in (None, written["target_id"])
    assert written["target_id"] not in written["comparison_ids"]
    # The per-star verdict lands beside the selection, not only in a log line.
    report = next(p for p in out.outputs if "comparison_screening" in p)
    table = pd.read_csv(report, sep="\t")
    assert {"star_id", "coverage", "eligible", "basic_reason", "role"} <= set(table.columns)
    assert (table["role"] == "target").sum() == 1


def test_catalog_order_is_labelled_as_a_fallback_not_a_selection(tmp_path):
    """With no photometry there is no ranking, and the file must admit it."""
    _master_with_source_ids(tmp_path)
    ctx = _ctx(tmp_path, _params(tmp_path, lc_target_id=1,
                                 lc_comparison_mode="auto"))

    out = LcTargetStep().run(ctx)
    assert out.status is StepStatus.OK, out.message
    assert "fell back to catalogue order" in out.message
    written = json.loads(read_selection_path(tmp_path).read_text(encoding="utf-8"))
    assert written["selected_by"] == "catalog_order"
    assert "check_id" not in written


def test_manual_mode_is_never_second_guessed_by_the_screen(tmp_path):
    """The user named those stars; a batch run must not substitute others."""
    _master_with_source_ids(tmp_path)
    _photometry(tmp_path)
    ctx = _ctx(tmp_path, _params(tmp_path, lc_target_id=1,
                                 lc_comparison_mode="manual",
                                 lc_comparison_ids="4,5,6"))

    out = LcTargetStep().run(ctx)
    assert out.status is StepStatus.OK, out.message
    written = json.loads(read_selection_path(tmp_path).read_text(encoding="utf-8"))
    assert written["comparison_ids"] == [4, 5, 6]
    assert written["selected_by"] == "manual"


def test_the_window_screens_through_the_shared_function(tmp_path, monkeypatch):
    """Structural parity: the window must not keep its own copy of the screen.

    Asserting the two agree on one dataset is weaker than asserting there is
    only one implementation — an assertion that keeps holding as either side
    changes.
    """
    import apex.gui.workflow.lc.step8_target_selection as window_module
    from apex.analysis.light_curve import comparison_screening

    assert window_module.build_candidate_pool is comparison_screening.build_candidate_pool
    assert window_module.screen_measurements is comparison_screening.screen_measurements


def read_selection_path(result_dir) -> Path:
    from apex.utils.step_paths_lc import step8_selection_dir

    return Path(step8_selection_dir(result_dir)) / "lc_target_selection.json"
