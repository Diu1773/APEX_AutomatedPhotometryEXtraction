"""A posterior sitting on its own bound must say so.

Twice now a bound has been read as a result. M67's [M/H] came back at −0.83
and looked like a metal-poor cluster until the U zero-point was re-anchored;
NGC 6811's age sat at the sampler's floor without an E(B−V) prior. Both were
the edge of the box, not a measurement.

The trap is still live in the GUI: its default age window is 0.2–6 Gyr, which
cannot reach a globular cluster at all. Fitting M13 with the defaults returns
6 Gyr — a number that means "at least 6", printed exactly like a number that
means "6". These pin the warning that tells the difference.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from apex.analysis.cmd.isochrone_fit_service import _railed_parameters

BOUNDS = SimpleNamespace(age_bounds=(8.5, 10.2), mh_bounds=(-1.0, 0.5))


def _summary(log_age=9.5, mh=-0.2):
    return {"log_age": [log_age - 0.05, log_age, log_age + 0.05],
            "metallicity": [mh - 0.05, mh, mh + 0.05]}


def test_an_age_on_the_ceiling_is_flagged():
    """M13 fitted inside the GUI's 0.2-6 Gyr default."""
    messages = _railed_parameters(_summary(log_age=10.2), BOUNDS)
    assert any("상한" in m and "나이" in m for m in messages)


def test_an_age_on_the_floor_is_flagged():
    messages = _railed_parameters(_summary(log_age=8.5), BOUNDS)
    assert any("하한" in m and "나이" in m for m in messages)


def test_a_metallicity_on_the_floor_is_flagged():
    """The M67 [M/H] rail, which a bad U zero-point produced."""
    messages = _railed_parameters(_summary(mh=-1.0), BOUNDS)
    assert any("하한" in m and "[M/H]" in m for m in messages)


def test_a_posterior_in_the_middle_says_nothing():
    assert _railed_parameters(_summary(), BOUNDS) == []


def test_the_edge_scales_with_the_box_not_the_units():
    """2 % of a wide box is a wider absolute margin than 2 % of a narrow one."""
    narrow = SimpleNamespace(age_bounds=(9.0, 9.1), mh_bounds=(-1.0, 0.5))
    # 9.09 is 10 % from the top of a 0.1-wide box: not railed.
    assert not any("나이" in m for m in _railed_parameters(_summary(log_age=9.09), narrow))
    # The same 0.01 gap in the default 1.7-wide box is well inside 2 %.
    assert any("나이" in m for m in _railed_parameters(_summary(log_age=10.19), BOUNDS))


@pytest.mark.parametrize("value", [float("nan"), float("inf")])
def test_a_broken_median_is_not_a_rail(value):
    assert _railed_parameters(_summary(log_age=value), BOUNDS) == []


def test_a_degenerate_box_is_skipped():
    flat = SimpleNamespace(age_bounds=(9.0, 9.0), mh_bounds=(-1.0, 0.5))
    assert _railed_parameters(_summary(log_age=9.0), flat) == []
