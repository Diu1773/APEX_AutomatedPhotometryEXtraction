"""Step 11 was deferred with the reason "nothing to port". That was half right.

The window loads the Step 10 table and opens an interactive viewer — region
selection, parallax and proper-motion sliders, quality masks, a Teff colour bar.
There was no calculation trapped behind a `QThread`, so nothing to free, and the
registry said so.

The reason held for the viewer and not for the figure. A run that measures a
cluster should leave the picture of it, and the picture needs none of that
machinery. So this is a new export rather than a move — which is why it needed a
decision and not a refactor (2026-08-19).
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from matplotlib.figure import Figure

from apex.analysis.cmd.cmd_plot import choose_axes, draw_cmd, export_cmd_plot


def _table(tmp_path, bands=("B", "V", "R"), n=300):
    out = tmp_path / "cmd_zeropoint"
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20260819)
    body = {"ID": np.arange(n)}
    base = rng.uniform(12.0, 18.0, n)
    for i, b in enumerate(bands):
        body[f"mag_std_{b}"] = base + 0.4 * i + rng.normal(0, 0.05, n)
        body[f"mag_std_err_{b}"] = np.abs(rng.normal(0.02, 0.01, n))
    pd.DataFrame(body).to_csv(out / "median_by_ID_filter_wide_cmd.csv", index=False)
    return SimpleNamespace(P=SimpleNamespace(result_dir=str(tmp_path)))


def test_the_config_decides_the_axes_when_it_says(tmp_path):
    """So this figure and the isochrone fit are drawn on the same axes."""
    params = _table(tmp_path)
    params.P.iso_colors = "B-V"
    params.P.iso_mag_band = "V"
    df = pd.read_csv(tmp_path / "cmd_zeropoint" / "median_by_ID_filter_wide_cmd.csv")
    blue, red, mag, why = choose_axes(df, params)
    assert (blue, red, mag) == ("B", "V", "V")
    assert "isochrone.colors" in why


def test_it_falls_back_to_the_widest_baseline_and_says_so(tmp_path):
    """An unexplained axis choice is how two figures look like two datasets."""
    _table(tmp_path)
    df = pd.read_csv(tmp_path / "cmd_zeropoint" / "median_by_ID_filter_wide_cmd.csv")
    blue, red, _, why = choose_axes(df, None)
    assert (blue, red) == ("B", "R")
    assert "widest" in why


def test_the_magnitude_axis_is_inverted(tmp_path):
    """A CMD with brighter at the bottom is not a CMD."""
    params = _table(tmp_path)
    df = pd.read_csv(tmp_path / "cmd_zeropoint" / "median_by_ID_filter_wide_cmd.csv")
    fig = Figure(figsize=(7, 7.6))
    assert draw_cmd(fig, df, params)
    lo, hi = fig.axes[0].get_ylim()
    assert lo > hi, "밝은 별이 위로 가야 한다"


def test_one_band_cannot_make_a_colour(tmp_path):
    _table(tmp_path, bands=("V",))
    df = pd.read_csv(tmp_path / "cmd_zeropoint" / "median_by_ID_filter_wide_cmd.csv")
    fig = Figure(figsize=(7, 7.6))
    assert not draw_cmd(fig, df, None)
    assert "two calibrated bands" in fig.axes[0].texts[0].get_text()


def test_the_export_writes_the_figure(tmp_path):
    params = _table(tmp_path)
    written = export_cmd_plot(tmp_path, params)
    assert [p.name for p in written] == ["step11_cmd.png"]
    assert written[0].stat().st_size > 0


def test_a_missing_table_blocks_rather_than_crashes(tmp_path):
    from apex.pipeline.steps.cmdplot import CmdPlotStep
    from apex.pipeline.base import StepStatus
    from apex.pipeline.context import RunContext
    import logging

    params = SimpleNamespace(P=SimpleNamespace(result_dir=str(tmp_path)))
    ctx = RunContext(mode="cmd", params=params, result_dir=tmp_path,
                     data_dir=tmp_path, logger=logging.getLogger("test"))
    assert CmdPlotStep().run(ctx).status == StepStatus.BLOCKED


def test_step_11_is_no_longer_deferred():
    from apex.pipeline.base import DeferredStep
    from apex.pipeline.registry import get_steps
    from apex.pipeline.steps.cmdplot import CmdPlotStep

    step = next(s for s in get_steps("cmd") if s.index == 11)
    assert isinstance(step, CmdPlotStep)
    assert not isinstance(step, DeferredStep)

    # Only the interactive editor should be left.
    deferred = {s.key for s in get_steps("cmd") if isinstance(s, DeferredStep)}
    assert deferred == {"masterid"}
