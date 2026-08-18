"""Step 5 had no figure — not even in the window.

Every other QC figure this week was a move: the window drew it and kept it to
itself. Step 5 was the exception. `frame_wcs_qc.csv` carries 36 columns per
frame and nothing ever looked at them together, in either route.

Which four panels was decided by measuring the coefficient of variation of every
candidate across 137 real frames from five clusters, not by taste. `match_rate`
came out at 0.01 — a flat line on every dataset — so it is a number in the title
and not an axis. `resid_vs_radius_slope` came out at 1.75, the most variable, and
it is the only one that changes a decision.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from matplotlib.figure import Figure

from apex.analysis.wcs_qc_plots import (
    draw_wcs_qc_overview, export_wcs_qc, load_wcs_qc,
)


def _table(tmp_path, n=12, slope=-0.4):
    out = tmp_path / "step5_wcs"
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20260818)
    pd.DataFrame({
        "file": [f"f{i:03d}.fit" for i in range(n)],
        "solver": ["apex_internal"] * n,
        "rms_px": np.abs(rng.normal(0.42, 0.06, n)),
        "resid_p99_px": np.abs(rng.normal(2.8, 0.4, n)),
        "resid_vs_radius_slope": rng.normal(slope, 0.1, n),
        "edge_resid_ratio": np.abs(rng.normal(1.2, 0.2, n)),
        "n_detect": rng.integers(600, 1300, n),
        "n_match": rng.integers(600, 1300, n),
        "center_offset_arcsec": np.abs(rng.normal(39.0, 0.6, n)),
        "match_rate": np.full(n, 0.9989),
        "inlier_rate": np.full(n, 0.90),
        "scale_delta_pct": np.full(n, 0.51),
        "wcs_qc_pass": [True] * (n - 1) + [False],
    }).to_csv(out / "frame_wcs_qc.csv", index=False)
    return out


def test_the_table_is_read_from_where_step5_writes_it(tmp_path):
    _table(tmp_path)
    df = load_wcs_qc(tmp_path)
    assert len(df) == 12
    assert "resid_vs_radius_slope" in df.columns


def test_a_missing_table_is_not_an_error(tmp_path):
    assert load_wcs_qc(tmp_path).empty
    assert export_wcs_qc(tmp_path) == []


def test_the_four_panels_are_drawn(tmp_path):
    _table(tmp_path)
    fig = Figure(figsize=(11, 7.4))
    assert draw_wcs_qc_overview(fig, load_wcs_qc(tmp_path))
    titles = {ax.get_title() for ax in fig.axes}
    assert titles == {"Residual size", "Distortion left over",
                      "Matched vs detected", "Pointing offset"}


def test_the_flat_metrics_are_in_the_title_not_on_an_axis(tmp_path):
    """A panel that is flat on every dataset takes the space of one that is not."""
    _table(tmp_path)
    fig = Figure(figsize=(11, 7.4))
    draw_wcs_qc_overview(fig, load_wcs_qc(tmp_path))
    title = fig._suptitle.get_text()
    assert "match 0.9989" in title
    assert "inlier 0.900" in title
    assert "scale Δ +0.51%" in title
    assert "PASS 11/12" in title


def test_an_empty_table_says_so_rather_than_drawing_empty_axes(tmp_path):
    fig = Figure(figsize=(11, 7.4))
    assert not draw_wcs_qc_overview(fig, pd.DataFrame())
    assert "No WCS QC table" in fig.axes[0].texts[0].get_text()


def test_the_export_writes_the_overview(tmp_path):
    _table(tmp_path)
    written = export_wcs_qc(tmp_path)
    assert [p.name for p in written] == ["step5_wcs_qc_overview.png"]
    assert written[0].stat().st_size > 0


def test_the_step_writes_it_too():
    import io
    import tokenize

    source = (Path(__file__).absolute().parents[1]
              / "apex/pipeline/steps/wcs.py").read_text(encoding="utf-8")
    code = " ".join(
        tok.string
        for tok in tokenize.generate_tokens(io.StringIO(source).readline)
        if tok.type not in (tokenize.COMMENT, tokenize.STRING)
    )
    assert "export_wcs_qc" in code
