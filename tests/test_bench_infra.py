"""Benchmark infrastructure: measurement envelope, worker override, counters."""

from __future__ import annotations

import time

import numpy as np
import pytest

from apex.benchmark import resources
from apex.utils import photometry_loader
from apex.utils.constants import get_parallel_workers


def test_measure_records_wall_env_and_rss():
    with resources.measure("unit-test", extra={"tag": 1}) as metrics:
        time.sleep(0.05)
        _ = np.zeros((256, 256))
    assert metrics["wall_s"] >= 0.04
    assert metrics["tag"] == 1
    assert metrics["env"]["python"]
    assert "packages" in metrics["env"]
    if resources.HAS_PSUTIL:
        assert metrics["peak_rss_mb"] > 0
        assert metrics["peak_rss_total_mb"] >= metrics["peak_rss_mb"]
    else:  # pragma: no cover - psutil is an optional extra
        assert metrics["peak_rss_mb"] is None


def test_measure_command_captures_child(tmp_path):
    if not resources.HAS_PSUTIL:
        pytest.skip("psutil not installed (bench extra)")
    import sys

    rc, metrics = resources.measure_command(
        [sys.executable, "-c", "import time; x = bytearray(50_000_000); time.sleep(0.6)"],
        label="child-alloc")
    assert rc == 0
    # The 50 MB allocation lives in the child; the poller must see it.
    assert metrics["peak_rss_total_mb"] > 40


def test_env_var_pins_worker_count(monkeypatch):
    monkeypatch.setenv("APEX_MAX_WORKERS", "3")
    assert get_parallel_workers() == 3
    monkeypatch.setenv("APEX_MAX_WORKERS", "1")
    assert get_parallel_workers() == 1


def test_env_var_invalid_values_fall_through(monkeypatch):
    monkeypatch.setenv("APEX_MAX_WORKERS", "0")
    assert get_parallel_workers() >= 2
    monkeypatch.setenv("APEX_MAX_WORKERS", "banana")
    assert get_parallel_workers() >= 2
    monkeypatch.delenv("APEX_MAX_WORKERS")
    assert get_parallel_workers() >= 2


def test_env_var_is_bounded_against_typos(monkeypatch):
    import os

    monkeypatch.setenv("APEX_MAX_WORKERS", "99999")
    assert get_parallel_workers() <= 4 * (os.cpu_count() or 4)


def test_load_counters_reset_and_shape():
    photometry_loader.reset_load_counters()
    counters = photometry_loader.get_load_counters()
    assert set(counters) == {"frames_loaded", "rows_loaded",
                             "sidmap_built", "sidmap_reused"}
    assert all(v == 0 for v in counters.values())
    # get_load_counters must be a snapshot, not the live dict
    counters["frames_loaded"] = 99
    assert photometry_loader.get_load_counters()["frames_loaded"] == 0
