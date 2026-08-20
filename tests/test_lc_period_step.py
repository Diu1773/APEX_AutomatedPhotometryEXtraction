"""LC step 11: the period, headless.

Nothing had to move for this one. `PeriodAnalysisWorker` is a pass-through — it
takes arrays, calls `run_period_analysis()`, emits the dict — and loading,
file choice and saving were already services. What was missing was a caller
that does not need someone to type a range into a spin box.

Measured against the analysis the window saved (YZ Boo, two nights, 364 points,
`E:/APEX_validation/reprocess/YZBoo_2n`), with the same filter grouping and the
same input file: all four methods agree to 0.0e+00, and PDM lands 0.19 % from
the literature period (0.104092 d). That needs the workspace, so it lives in
`docs/audit/LC_HEADLESS_STEPS_8_TO_11.md` rather than here; what is testable here is
the wiring and the refusals.
"""

from __future__ import annotations

import json
import logging
from types import SimpleNamespace

import pandas as pd
import pytest

from apex.analysis.light_curve.target_config import LcTarget, write_selection
from apex.pipeline.base import StepStatus
from apex.pipeline.context import RunContext
from apex.pipeline.steps.lc_period import LcPeriodStep, _methods
from apex.utils.step_paths_lc import step9_lc_dir, step11_period_dir


def _params(tmp_path, **over):
    body = dict(
        result_dir=str(tmp_path), data_dir=str(tmp_path),
        lc_target_id=-1, lc_filter="",
        lc_period_min_days=0.01, lc_period_max_days=10.0,
        lc_period_samples_per_peak=10, lc_period_methods="ls,pdm",
        lc_period_pdm_bins=10,
    )
    body.update(over)
    return SimpleNamespace(P=SimpleNamespace(**body))


def _ctx(tmp_path, params=None):
    return RunContext(mode="lc", params=params or _params(tmp_path),
                      result_dir=tmp_path, data_dir=tmp_path,
                      logger=logging.getLogger("test"))


def _curve(tmp_path, target_id=5, period=0.1, n=200):
    """A sinusoid at a known period, so the answer is checkable."""
    import numpy as np

    out = step9_lc_dir(tmp_path)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20260819)
    t = np.sort(rng.uniform(0, 3.0, n))
    mag = 12.0 + 0.30 * np.sin(2 * np.pi * t / period) + rng.normal(0, 0.004, n)
    pd.DataFrame({
        "file": [f"f{i:04d}.fit" for i in range(n)],
        "filter": ["V"] * n,
        "BJD_TDB": t + 2460000.0,
        "mag": mag,
        "mag_err": np.full(n, 0.004),
        "diff_mag_raw": mag - 12.0,
        "diff_err": np.full(n, 0.006),
    }).to_csv(out / f"lightcurve_ID{target_id}_raw.csv", index=False)


def test_it_blocks_when_no_target_was_chosen(tmp_path):
    result = LcPeriodStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.BLOCKED
    assert "Step 8" in result.message or "target_id" in result.message


def test_it_blocks_when_step_9_has_produced_no_curve(tmp_path):
    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    result = LcPeriodStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.BLOCKED
    assert "Step 9" in result.message


def test_it_blocks_on_an_empty_search_window(tmp_path):
    """A period outside the range is simply not found, so an inverted range
    would return a confident wrong answer rather than nothing."""
    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    _curve(tmp_path)
    params = _params(tmp_path, lc_period_min_days=5.0, lc_period_max_days=0.5)
    result = LcPeriodStep().run(_ctx(tmp_path, params))
    assert result.status == StepStatus.BLOCKED
    assert "search window" in result.message


def test_it_recovers_a_known_period(tmp_path):
    """0.1 d sinusoid in, 0.1 d out — the wiring carries real numbers."""
    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    _curve(tmp_path, period=0.1)
    result = LcPeriodStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.OK, result.message

    files = sorted(step11_period_dir(tmp_path).glob("period_analysis_*.json"))
    assert files, "no analysis written"
    body = json.loads(files[0].read_text(encoding="utf-8"))
    periods = [v.get("best_period") for v in (body.get("results") or {}).values()
               if isinstance(v, dict) and v.get("best_period")]
    assert periods, f"no method reported a period: {body.get('results', {}).keys()}"
    assert any(abs(p - 0.1) < 0.002 for p in periods), (
        f"injected 0.1 d, recovered {periods}")


def test_it_analyses_the_combined_curve_and_each_band(tmp_path):
    """The combined curve is the answer; the per-band runs are the check.

    Splitting by band only was this step's behaviour until 2026-08-21, and it
    cost the answer. On YZ Boo each band holds 124 / 119 / 21 points against 364
    combined, and that density is what separates the true period from its
    one-cycle alias: per band the step returned g = 0.105272 d with the alias
    unresolved (+1.13 % from the literature), while the combined curve resolved
    to 0.104209 d (+0.11 %). Disagreement between bands is still a real
    diagnostic, so both are produced.
    """
    import numpy as np

    write_selection(tmp_path, LcTarget(target_id=5), [1, 2])
    out = step9_lc_dir(tmp_path)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(7)
    rows = []
    for flt in ("g", "r"):
        t = np.sort(rng.uniform(0, 3.0, 150))
        mag = 12.0 + 0.3 * np.sin(2 * np.pi * t / 0.1) + rng.normal(0, 0.004, 150)
        rows.append(pd.DataFrame({
            "file": [f"{flt}{i:04d}.fit" for i in range(150)], "filter": flt,
            "BJD_TDB": t + 2460000.0, "mag": mag, "mag_err": 0.004,
            "diff_mag_raw": mag - 12.0, "diff_err": 0.006,
        }))
    pd.concat(rows, ignore_index=True).to_csv(out / "lightcurve_ID5_raw.csv", index=False)

    result = LcPeriodStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.OK, result.message
    names = {p.name for p in step11_period_dir(tmp_path).glob("period_analysis_*.json")}
    assert names == {"period_analysis_all_ID5.json",
                     "period_analysis_g_ID5.json",
                     "period_analysis_r_ID5.json"}, names


def test_a_single_band_curve_is_not_analysed_twice(tmp_path):
    """With one band, `all` and that band are the same rows — one run, not two."""
    write_selection(tmp_path, LcTarget(target_id=5), [1, 2])
    _curve(tmp_path)

    result = LcPeriodStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.OK, result.message
    names = {p.name for p in step11_period_dir(tmp_path).glob("period_analysis_*.json")}
    assert names == {"period_analysis_V_ID5.json"}, names


def test_the_configured_filter_wins(tmp_path):
    write_selection(tmp_path, LcTarget(target_id=5), [1, 2])
    _curve(tmp_path)
    result = LcPeriodStep().run(_ctx(tmp_path, _params(tmp_path, lc_filter="all")))
    assert result.status == StepStatus.OK, result.message
    names = {p.name for p in step11_period_dir(tmp_path).glob("period_analysis_*.json")}
    assert names == {"period_analysis_all_ID5.json"}, names


def test_methods_parse_and_never_come_back_empty(tmp_path):
    """Two different empties, deliberately two different answers.

    An unset value means "I did not choose", so it gets the default pair.
    A value that is present but yields no method means something was meant and
    could not be read — that falls back to the one method that always works,
    rather than silently running the pair the user did not ask for.
    """
    assert _methods(_params(tmp_path, lc_period_methods="LS, PDM")) == ["ls", "pdm"]
    assert _methods(_params(tmp_path, lc_period_methods="ls;bls")) == ["ls", "bls"]
    assert _methods(_params(tmp_path, lc_period_methods="")) == ["ls", "pdm"]
    assert _methods(_params(tmp_path, lc_period_methods="  ,  ")) == ["ls"]
    assert _methods(_params(tmp_path, lc_period_methods=None)) == ["ls", "pdm"]


def test_the_new_settings_are_actually_read():
    """The guard that rejected six unread rows when step 8 was built."""
    from apex.config.config_audit import unread_settings

    dead = set(unread_settings(mode="lc")["dead"])
    for name in ("lc_period_min_days", "lc_period_max_days",
                 "lc_period_samples_per_peak", "lc_period_methods",
                 "lc_period_pdm_bins"):
        assert name not in dead, f"{name} 을 읽는 코드가 없다"


def test_the_lc_pipeline_now_reaches_step_11():
    """Step 10 was deferred for one day, on the reading that its calculation
    took its inputs from widgets. It did not — `_sync_state_from_controls` had
    been copying them into plain attributes all along — so it runs headless
    too, and the LC branch has no deferred step left."""
    from apex.pipeline.registry import get_steps

    steps = get_steps("lc")
    assert [s.index for s in steps] == list(range(1, 12))
    assert [s.key for s in steps][-4:] == [
        "lctarget", "lclightcurve", "lcdetrend", "lcperiod"]


def test_the_detrend_recognises_a_window_produced_result(tmp_path):
    """A batch run must not redo work the window already did."""
    from apex.pipeline.registry import get_steps
    from apex.utils.step_paths_lc import step10_detrend_dir

    step = {s.index: s for s in get_steps("lc")}[10]
    ctx = _ctx(tmp_path)
    assert not step.is_complete(ctx)
    out = step10_detrend_dir(tmp_path)
    out.mkdir(parents=True, exist_ok=True)
    # `_current` is what both paths write last; a bare `_offset` is one mode's
    # file and does not mean the step finished.
    (out / "lightcurve_ID5_current.csv").write_text("a\n1\n")
    assert step.is_complete(ctx)


# ── the summary figure ──────────────────────────────────────────────────────
#
# The window drew it and nobody else could. It is the picture that says whether
# a period is real — a clean fold at the adopted period, or a smear that says
# the peak was an alias of the sampling window.

PLOT_SHARED = [
    "_update_summary_plot",
    "_summary_data_type",
    "_summary_period",
    "_summary_method_result",
    "_median_center_for_plot",
    "_check_star_plot_data",
    "_check_star_ls_result",
    "_two_harmonic_phase_model",
    "_phase_plot_check_time_column",
    "_summary_uses_compact_layout",
    "_summary_uses_stacked_layout",
]

# What the window supplies from widgets, and the batch path is handed. An
# accidental addition here means the batch figure silently uses a default the
# user never chose; an accidental removal means the window stops reading its
# own controls — which is exactly what happened while this was being built.
PLOT_WINDOW_ONLY = [
    "_summary_figure",
    "_summary_canvas_width",
    "_summary_redraw",
    "_search_window",
    "_show_alias_marks",
    "_load_check_star",
]


def test_the_figure_code_does_not_reach_for_qt():
    import subprocess
    import sys

    code = ("import sys;"
            "import apex.analysis.light_curve.period_plot as m;"
            "print('PyQt5' in sys.modules)")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "False"


def test_the_figure_does_not_import_the_gui_layer_at_module_scope():
    """`analysis` must not depend on `gui`. The palette is looked up lazily so
    a themed window still recolours its plot, and a batch run needs nothing."""
    import ast
    import pathlib

    src = pathlib.Path("apex/analysis/light_curve/period_plot.py").read_text(encoding="utf-8")
    for node in ast.parse(src).body:
        mod = getattr(node, "module", "") or ""
        assert not mod.startswith("apex.gui"), (
            f"module-scope import of {mod} — analysis may not depend on gui")


def test_the_window_draws_the_analysis_figure():
    pytest.importorskip("PyQt5")
    from apex.analysis.light_curve.period_plot import PeriodSummaryPlotter
    from apex.gui.workflow.lc.step11_period_analysis import PeriodAnalysisWindow

    drifted = [n for n in PLOT_SHARED
               if getattr(PeriodAnalysisWindow, n, None) is not getattr(PeriodSummaryPlotter, n)]
    assert not drifted, (
        f"the window no longer inherits {drifted} — a second copy of the "
        "summary figure has appeared")


def test_the_window_still_reads_its_own_controls():
    """The regression this test exists for: moving the figure out left the
    window inheriting the batch defaults, so its spin boxes stopped mattering."""
    pytest.importorskip("PyQt5")
    from apex.analysis.light_curve.period_plot import PeriodSummaryPlotter
    from apex.gui.workflow.lc.step11_period_analysis import PeriodAnalysisWindow

    for name in PLOT_WINDOW_ONLY:
        assert getattr(PeriodAnalysisWindow, name) is not getattr(PeriodSummaryPlotter, name), (
            f"{name} must be the window's own — otherwise the window uses the "
            "batch default instead of what the user set")


def test_the_plotter_needs_nothing_the_class_does_not_declare():
    """Same guard as the light-curve builder, for the same reason: an attribute
    the window happened to set in its constructor fails deep inside a draw."""
    import ast
    import inspect

    from apex.analysis.light_curve.period_plot import PeriodSummaryPlotter

    cls = next(n for n in ast.parse(inspect.getsource(PeriodSummaryPlotter)).body
               if isinstance(n, ast.ClassDef))
    read, written = set(), set()
    for node in ast.walk(cls):
        if (isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name)
                and node.value.id == "self"):
            (written if isinstance(node.ctx, ast.Store) else read).add(node.attr)

    provided = set(dir(PeriodSummaryPlotter)) | written | {"__dict__"}
    missing = sorted(read - provided)
    assert not missing, f"the figure reads {missing} and the class never provides them"


def test_it_draws_a_figure_for_a_real_period(tmp_path):
    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    _curve(tmp_path, period=0.1)
    result = LcPeriodStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.OK, result.message

    figs = list(step11_period_dir(tmp_path).glob("period_summary_*.png"))
    assert figs, "no summary figure written"
    assert figs[0].stat().st_size > 20_000, "figure is suspiciously small"


def test_alias_resolution_runs_by_default(tmp_path):
    """Off by default it would adopt whichever peak Lomb-Scargle liked, and on
    a multi-night run that is routinely an alias of the sampling window."""
    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    _curve(tmp_path, period=0.1)
    LcPeriodStep().run(_ctx(tmp_path))

    body = json.loads(sorted(step11_period_dir(tmp_path)
                             .glob("period_analysis_*.json"))[0].read_text(encoding="utf-8"))
    assert body.get("alias_analysis"), "no alias analysis in the saved result"
    assert body["alias_analysis"].get("adopted_period")


def test_alias_resolution_can_be_turned_off(tmp_path):
    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    _curve(tmp_path, period=0.1)
    params = _params(tmp_path, lc_period_resolve_aliases=False)
    LcPeriodStep().run(_ctx(tmp_path, params))

    body = json.loads(sorted(step11_period_dir(tmp_path)
                             .glob("period_analysis_*.json"))[0].read_text(encoding="utf-8"))
    assert not body.get("alias_analysis")


def test_nightly_offset_detrending_forces_the_raw_series(tmp_path):
    """Nightly-offset detrending removes the between-night baseline a
    multi-night period search needs. Using the corrected series after it moved
    YZ Boo's adopted period from 0.11 % off the literature value to 1.22 %.
    """
    import numpy as np

    from apex.pipeline.steps.lc_period import _resolve_aliases

    n = 200
    rng = np.random.default_rng(11)
    t = np.sort(rng.uniform(0, 3.0, n))
    raw = 12.0 + 0.3 * np.sin(2 * np.pi * t / 0.1) + rng.normal(0, 0.004, n)
    lc_data = {
        "time": t, "mag_raw": raw, "mag_corr": raw + 5.0, "mag_err": np.full(n, 0.004),
        "night_id": np.zeros(n, dtype=int),
        "correction_preserves_nightly_baseline": False,
    }
    out = _resolve_aliases(lc_data, {}, 0.01, 10.0, 10)
    assert out is not None
    assert out["input_series"] == "raw", (
        "the corrected series was used after a correction that does not "
        "preserve the nightly baseline")
    assert out["nightly_baseline_preserved"] is False


def test_an_unresolved_alias_analysis_does_not_override_the_peak():
    """The service reports RESOLVED / AMBIGUOUS / INSUFFICIENT. Taking its top
    candidate regardless throws away the one thing it was asked to determine.

    AE UMa measured the cost. Two nights, 2.13 d baseline; the service returned
    AMBIGUOUS ("leave-one-night-out agreement is 50%") with rank 1 = 0.082514 d
    and rank 2 = 0.086079 d. The plain Lomb-Scargle peak was 0.086011 d and the
    literature period is 0.086017 d — so obeying rank 1 turned a 0.01 % answer
    into a 4.07 % one, and folded the figure at the wrong period.
    """
    from apex.analysis.light_curve.period_plot import PeriodSummaryPlotter

    plotter = PeriodSummaryPlotter()
    plotter.lc_data = {"mag_corr": None}
    plotter.results = {"raw_ls": {"best_period": 0.086011}}

    plotter.alias_analysis = {"status": "AMBIGUOUS", "adopted_period": 0.082514}
    assert plotter._summary_period() == pytest.approx(0.086011), (
        "an AMBIGUOUS alias analysis must not override the periodogram peak")

    plotter.alias_analysis = {"status": "INSUFFICIENT", "adopted_period": 0.082514}
    assert plotter._summary_period() == pytest.approx(0.086011)

    # ...and when it did resolve, it wins — that is what takes YZ Boo from
    # 9.2 % off the literature period to 0.11 %.
    plotter.alias_analysis = {"status": "RESOLVED", "adopted_period": 0.104209}
    assert plotter._summary_period() == pytest.approx(0.104209)

    plotter.alias_analysis = None
    assert plotter._summary_period() == pytest.approx(0.086011)


def test_the_step_reports_the_same_period_the_figure_folds_at(tmp_path):
    """One rule, two consumers. They disagreed until 2026-08-19: the message
    printed the alias candidate while the figure folded at the peak."""
    import numpy as np

    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    _curve(tmp_path, period=0.1)
    result = LcPeriodStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.OK, result.message

    body = json.loads(sorted(step11_period_dir(tmp_path)
                             .glob("period_analysis_*.json"))[0].read_text(encoding="utf-8"))
    status = str((body.get("alias_analysis") or {}).get("status", "")).upper()
    if status == "RESOLVED":
        assert "[alias" not in result.message
    else:
        assert "[alias" in result.message, (
            "an unresolved alias analysis must be visible in the step message, "
            f"not silently dropped (status={status})")


def test_the_filter_count_in_the_message_counts_filters(tmp_path):
    """It counted `written`, which also holds the figure and four periodogram
    CSVs — so a single-filter run reported "2 filter(s)"."""
    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    _curve(tmp_path)
    result = LcPeriodStep().run(_ctx(tmp_path, _params(tmp_path, lc_filter="all")))
    assert result.status == StepStatus.OK, result.message
    assert "1 filter(s)" in result.message, result.message


def test_it_writes_the_candidates_rather_than_only_a_pick(tmp_path):
    """Two nights do not determine one period, so the run lists what it found.

    Measured both ways on real objects: YZ Boo's Lomb-Scargle peak is 9.2 % from
    the literature value and its PDM peak 0.19 %; on AE UMa it is 0.01 % against
    0.04 %, the other way round. The software has no basis to choose between
    them, and choosing anyway publishes one of those errors as a number. So the
    step lays the candidates out — every periodogram peak and every ranked
    alias — and the answer is picked against what the object is known to do.
    """
    import csv

    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    _curve(tmp_path, period=0.1)
    result = LcPeriodStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.OK, result.message

    tables = sorted(step11_period_dir(tmp_path).glob("period_candidates_*.csv"))
    assert tables, "no candidate table written"
    rows = list(csv.DictReader(tables[0].open(encoding="utf-8")))
    assert rows, "candidate table is empty"

    columns = set(rows[0])
    assert {"method", "period_days", "period_hours", "source"} <= columns, columns
    assert any(r["source"] == "periodogram" for r in rows), (
        "the periodogram peaks themselves must be in the table")
    # hours as well as days: a period is compared against a literature value,
    # and those are quoted both ways.
    for r in rows:
        assert abs(float(r["period_hours"]) - float(r["period_days"]) * 24.0) < 1e-9
    assert any(abs(float(r["period_days"]) - 0.1) < 0.002 for r in rows), (
        f"injected 0.1 d is not among the candidates: "
        f"{[r['period_days'] for r in rows]}")


# ── the candidate table is the answer, so it has to be readable ────────────
#
# D-014 closed on "do not pick; lay the candidates out and let a person compare
# with the literature". That makes this table the deliverable, and its first
# version worked against itself: it sorted by period, so the resolver's rank-1
# candidate sat sixth of twelve; and `strength` held periodogram power (0-1,
# higher is stronger) beside raw BIC (about -2400, lower is better), so the
# strongest candidate carried the smallest-looking number.

def _candidates(tmp_path, target_id=5, flt="all"):
    from apex.pipeline.steps.lc_period import _write_candidates

    results = {
        "raw_ls": {"best_period": 0.094468, "best_power": 0.878},
        "corr_pdm": {"best_period": 0.105297, "best_power": 0.874},
    }
    alias = {
        "status": "AMBIGUOUS",
        "input_series": "raw",
        "candidates": [
            {"rank": 3, "period": 0.094554, "bic": -2370.86},
            {"rank": 1, "period": 0.104151, "bic": -2430.67},
            {"rank": 2, "period": 0.115729, "bic": -2379.85},
        ],
    }
    _write_candidates(tmp_path, target_id, flt, results, alias)
    return pd.read_csv(
        step11_period_dir(tmp_path) / f"period_candidates_{flt}_ID{target_id}.csv")


def test_the_ranked_candidate_leads_the_table(tmp_path):
    table = _candidates(tmp_path)
    assert table.iloc[0]["rank"] == 1
    assert table.iloc[0]["period_days"] == pytest.approx(0.104151)
    ranked = table[table["source"] == "alias candidate"]
    assert list(ranked["rank"]) == [1, 2, 3]
    # Periodogram peaks come after the ranked candidates, strongest first.
    peaks = table[table["source"] == "periodogram"]
    assert list(peaks["period_days"]) == pytest.approx([0.094468, 0.105297])
    assert table.index.get_loc(peaks.index[0]) > table.index.get_loc(ranked.index[-1])


def test_strength_says_what_kind_of_number_it_is(tmp_path):
    table = _candidates(tmp_path)
    kinds = set(table["strength_kind"])
    assert any("delta BIC" in k for k in kinds)
    assert any("periodogram power" in k for k in kinds)
    # Never one unlabelled column holding both.
    for _, row in table.iterrows():
        assert isinstance(row["strength_kind"], str) and row["strength_kind"]


def test_delta_bic_is_zero_for_the_best_and_positive_for_the_rest(tmp_path):
    """Raw BIC of -2430 in a column called `strength` reads as the weakest."""
    table = _candidates(tmp_path)
    ranked = table[table["source"] == "alias candidate"]
    assert ranked.iloc[0]["strength"] == pytest.approx(0.0)
    assert (ranked.iloc[1:]["strength"] > 0).all()
    assert ranked.iloc[1]["strength"] == pytest.approx(2430.67 - 2379.85, abs=1e-6)


def test_an_unresolved_run_says_so_next_to_its_front_runner(tmp_path):
    """The top row must not read as an answer when the run could not confirm it."""
    table = _candidates(tmp_path)
    note = str(table.iloc[0]["note"])
    assert "rank 1 of 3" in note
    assert "ambiguous" in note.lower()
    assert "adopted" not in note
