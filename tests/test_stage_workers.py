"""Per-stage worker ceilings (O2).

The ceilings are measured, not assumed: on the benchmark fixture forced
photometry runs ~4x SLOWER at every parallel count than serially, while detect
and wcs saturate at 2 workers (benchmark/perf/20260807/RESULTS.md).  These
tests pin the resulting policy so a later refactor cannot silently restore the
slow default.
"""

from __future__ import annotations

import os

import pytest

from apex.utils.constants import STAGE_WORKER_CAPS, get_parallel_workers


class _Params:
    """Minimal stand-in for the parameter object (`params.P.max_workers`)."""

    def __init__(self, max_workers):
        self.P = type("P", (), {"max_workers": max_workers})()


@pytest.fixture(autouse=True)
def _no_env_override(monkeypatch):
    monkeypatch.delenv("APEX_MAX_WORKERS", raising=False)


def test_forcedphot_runs_serially_by_default():
    assert get_parallel_workers(stage="forcedphot") == 1


def test_detect_and_wcs_are_capped_but_parallel():
    for stage in ("detect", "wcs"):
        workers = get_parallel_workers(stage=stage)
        assert 1 < workers <= STAGE_WORKER_CAPS[stage]


def test_unknown_stage_keeps_the_general_default():
    general = get_parallel_workers()
    assert get_parallel_workers(stage="isochrone") == general
    assert get_parallel_workers(stage=None) == general


def test_stage_cap_bounds_a_larger_configured_maximum():
    """`max_workers` is a maximum, so the stage ceiling may lower it."""
    params = _Params(16)
    assert get_parallel_workers(params, stage="forcedphot") == 1
    assert get_parallel_workers(params, stage="detect") == STAGE_WORKER_CAPS["detect"]
    # Without a stage the configured value still stands.
    assert get_parallel_workers(params) == 16


def test_configured_value_below_the_cap_is_respected():
    params = _Params(2)
    assert get_parallel_workers(params, stage="detect") == 2


def test_env_override_bypasses_the_stage_cap(monkeypatch):
    """The sweep must be able to probe above the ceilings it derives."""
    monkeypatch.setenv("APEX_MAX_WORKERS", "12")
    assert get_parallel_workers(stage="forcedphot") == 12
    assert get_parallel_workers(_Params(2), stage="detect") == 12


def test_result_is_always_at_least_one():
    for stage in list(STAGE_WORKER_CAPS) + [None, "unknown"]:
        assert get_parallel_workers(_Params(1), stage=stage) >= 1
