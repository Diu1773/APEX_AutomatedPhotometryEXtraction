"""The comparison screen: what it drops, what it says, and who calls it.

Two of these exist because the code failed in a way that produced a *plausible*
answer rather than an error. `_as_id_set` raised on a numpy array and the caller
caught it, so a pipeline run silently ranked comparisons by catalogue order and
reported success. And `lightcurve.filter` may hold the sentinel `all`, which is
not a filter name — asking the loader for it returned nothing, and the step said
"no usable photometry" about a workspace with 364 frames.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from apex.analysis.light_curve.comparison_screening import (
    MIN_COVERAGE,
    ScreeningFunnel,
    _as_id_set,
    build_candidate_pool,
    colors_from_catalog,
    screen_measurements,
)


def _measurements(n_frames: int = 40, n_stars: int = 12, *,
                  partial: dict | None = None, seed: int = 7) -> pd.DataFrame:
    """A long frame/star table shaped like the photometry loader's output."""
    rng = np.random.default_rng(seed)
    partial = partial or {}
    rows = []
    for star in range(1, n_stars + 1):
        keep = partial.get(star, n_frames)
        base = 12.0 + 0.25 * star
        for frame in range(keep):
            rows.append({
                "star_id": star,
                "frame": f"f{frame:04d}.fit",
                "mag": base + rng.normal(0.0, 0.004),
                "mag_err": 0.004,
            })
    return pd.DataFrame(rows)


def test_as_id_set_accepts_a_numpy_array():
    """`if not values` on an ndarray raises; the caller then swallows it."""
    values = np.array([3, 1, 2], dtype="int64")
    assert _as_id_set(values) == {1, 2, 3}
    assert _as_id_set(np.array([], dtype="int64")) == set()
    assert _as_id_set(None) == set()


def test_pool_drops_low_coverage_and_records_the_reason():
    frame = _measurements(n_frames=40, n_stars=8, partial={5: 10, 6: 12})
    pool, report, target_mag = build_candidate_pool(frame, target_id=1, target_count=4)

    assert 5 not in pool and 6 not in pool
    reasons = dict(zip(report["star_id"], report["basic_reason"]))
    assert reasons[5] == "low_coverage"
    assert reasons[6] == "low_coverage"
    assert reasons[1] == "target"
    assert 1 not in pool
    assert np.isfinite(target_mag)
    # Every measured star keeps a row, dropped ones included — that is what
    # makes "3 out of 4000" distinguishable from "3 out of 5".
    assert len(report) == 8


def test_pool_reason_records_what_removed_a_star_first():
    """A star already out on coverage is not relabelled by a later screen."""
    frame = _measurements(n_frames=40, n_stars=6, partial={4: 5})
    _pool, report, _mag = build_candidate_pool(
        frame, target_id=1, target_count=3, variable_ids=[4], manual_rejects=[3]
    )
    reasons = dict(zip(report["star_id"], report["basic_reason"]))
    assert reasons[4] == "low_coverage"
    assert reasons[3] == "manual_reject"


def test_pool_applies_manual_gaia_and_simbad_rejects():
    frame = _measurements(n_frames=30, n_stars=10)
    pool, report, _mag = build_candidate_pool(
        frame, target_id=1, target_count=5,
        manual_rejects=[2], variable_ids=[3], external_rejects=[4],
    )
    reasons = dict(zip(report["star_id"], report["basic_reason"]))
    assert reasons[2] == "manual_reject"
    assert reasons[3] == "gaia_variable"
    assert reasons[4] == "simbad_variable"
    assert {2, 3, 4}.isdisjoint(pool)


def test_pool_prefers_the_ids_it_is_asked_to_prefer():
    frame = _measurements(n_frames=30, n_stars=12)
    pool, _report, _mag = build_candidate_pool(
        frame, target_id=1, target_count=4, prefer_ids=[11, 12], pool_cap=6
    )
    assert pool[:2] == [11, 12]


def test_is_variable_callable_matches_the_id_list():
    """The window passes a predicate, a batch run passes a set — same answer."""
    frame = _measurements(n_frames=30, n_stars=8)
    by_ids, report_ids, _ = build_candidate_pool(
        frame, target_id=1, target_count=4, variable_ids=[5, 6])
    by_call, report_call, _ = build_candidate_pool(
        frame, target_id=1, target_count=4, is_variable=lambda sid: sid in (5, 6))
    assert by_ids == by_call
    assert list(report_ids["basic_reason"]) == list(report_call["basic_reason"])


def test_screen_refuses_rather_than_returning_a_short_ensemble():
    """Three comparisons plus a check star, or no result at all."""
    frame = _measurements(n_frames=30, n_stars=3)
    with pytest.raises(ValueError, match="coverage"):
        screen_measurements(frame, target_id=1, filter_key="g")


def test_screen_refuses_when_the_target_was_never_measured():
    frame = _measurements(n_frames=30, n_stars=8)
    with pytest.raises(ValueError, match="no usable measurements"):
        screen_measurements(frame, target_id=999, filter_key="g")


def test_screen_reports_a_funnel_that_adds_up():
    frame = _measurements(n_frames=40, n_stars=14, partial={13: 8, 14: 9})
    result = screen_measurements(frame, target_id=1, filter_key="g", desired_count=5)

    f = result.funnel
    assert f.measured == 14
    assert f.coverage == 12                 # two below the coverage floor
    assert f.eligible == 11                 # minus the target
    assert f.pool == len(result.candidate_ids)
    assert f.adopted == len(result.selected_ids)
    assert f.measured >= f.coverage >= f.eligible >= f.pool >= f.adopted
    assert "measured" in f.as_text() and "adopted" in f.as_text()
    assert f.as_dict()["measured"] == 14


def test_coverage_floor_is_the_documented_one():
    """A star at exactly the floor stays in; the screen is `>=`, not `>`."""
    frame = _measurements(n_frames=100, n_stars=6, partial={4: int(100 * MIN_COVERAGE)})
    _pool, report, _mag = build_candidate_pool(frame, target_id=1, target_count=3)
    reasons = dict(zip(report["star_id"], report["basic_reason"]))
    assert reasons[4] == ""


def test_colors_never_invent_a_zero():
    """A missing colour must be NaN — 0.0 reads as 'same colour as the target'."""
    catalog = pd.DataFrame({
        "source_id": [10, 11, 12],
        "phot_bp_mean_mag": [15.0, np.nan, 14.0],
        "phot_rp_mean_mag": [14.0, 13.0, np.nan],
    })
    colors = colors_from_catalog(catalog, np.array([10, 11, 12, 99]))
    assert colors[10] == pytest.approx(1.0)
    assert np.isnan(colors[11])
    assert np.isnan(colors[12])
    assert 99 not in colors


def test_colors_prefer_a_precomputed_column():
    catalog = pd.DataFrame({
        "source_id": [10],
        "color_gr": [0.42],
        "phot_bp_mean_mag": [15.0],
        "phot_rp_mean_mag": [14.0],
    })
    assert colors_from_catalog(catalog, [10])[10] == pytest.approx(0.42)


def test_funnel_default_is_all_zero():
    assert ScreeningFunnel().as_dict() == {
        "measured": 0, "coverage": 0, "eligible": 0, "pool": 0, "adopted": 0,
    }


def test_the_stability_search_return_is_carried_whole():
    """Listing the keys a caller needs makes the rest vanish without an error.

    The window's stability report reads `removed_ids` and its preview reads
    `active_ids`. Both come from `select_stable_comparisons`, and a first cut of
    this extraction copied out `metrics` alone — emptying both silently, on a
    path no test touched.
    """
    frame = _measurements(n_frames=40, n_stars=14)
    result = screen_measurements(frame, target_id=1, filter_key="g", desired_count=5)

    assert {"selected_ids", "active_ids", "metrics", "residuals", "removed_ids"} <= set(
        result.stability
    )
    assert result.stability["metrics"] is result.metrics
