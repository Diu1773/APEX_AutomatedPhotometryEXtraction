"""RAM admission control for the worker pool (O2).

A count picked from the core count alone is a promise the machine may not be
able to keep: on a 16 GB laptop with 1.9 GB free, four detection workers do
not fit, and what they cost is not just their own memory — they evict the OS
page cache holding the frames, so the next stage re-reads all of them from
disk (measured: forced photometry 222.6 s -> 387.4 s on identical input).

The per-worker figure is measured, not assumed, and is carried as a multiple
of the frame so a different camera changes the frame size rather than the
constant.
"""

from __future__ import annotations

import pytest

from apex.utils import constants as C
from apex.utils.constants import (
    WORKER_MEM_FRAME_MULTIPLE,
    get_parallel_workers,
    get_worker_decisions,
    reset_worker_decisions,
    workers_for_memory,
)

FRAME_61MP = 3194 * 4788 * 4          # the NGC 6811 fixture: 61.2 MB float32


@pytest.fixture(autouse=True)
def _clean(monkeypatch):
    monkeypatch.delenv("APEX_MAX_WORKERS", raising=False)
    reset_worker_decisions()


class _Params:
    def __init__(self, max_workers):
        self.P = type("P", (), {"max_workers": max_workers})()


# ── the memory model ───────────────────────────────────────────────────────

def test_worker_budget_scales_with_available_ram():
    plenty = workers_for_memory(FRAME_61MP, available_mb=16_000)
    tight = workers_for_memory(FRAME_61MP, available_mb=1_900)
    assert plenty > tight >= 1


def test_a_bigger_frame_buys_fewer_workers():
    """The cost is per frame, so a 4x larger sensor gets ~4x fewer workers."""
    small = workers_for_memory(FRAME_61MP, available_mb=16_000)
    large = workers_for_memory(FRAME_61MP * 4, available_mb=16_000)
    assert large < small
    assert large == pytest.approx(small / 4, abs=1)


def test_the_measured_case_reproduces():
    """1.9 GB free and a 61 MP frame is where this laptop actually sat."""
    # budget = 1900 * 0.6 - 400 = 740 MB;  per worker = 6.1 * 61.2 = 373 MB
    assert workers_for_memory(FRAME_61MP, available_mb=1_900) == 1
    # With 8 GB free the same frame supports several.
    assert workers_for_memory(FRAME_61MP, available_mb=8_000) >= 4


def test_never_returns_zero_workers():
    assert workers_for_memory(FRAME_61MP, available_mb=1) == 1


def test_unknown_inputs_yield_no_opinion():
    """None means 'cannot decide', so the caller keeps the CPU-based count."""
    assert workers_for_memory(None, available_mb=16_000) is None
    assert workers_for_memory(0, available_mb=16_000) is None


def test_frame_multiple_matches_the_measurement():
    """Guards the documented fit: 374 MB per worker on a 61.2 MB frame."""
    per_worker_mb = WORKER_MEM_FRAME_MULTIPLE * FRAME_61MP / 1e6
    assert per_worker_mb == pytest.approx(374, abs=10)


# ── how it composes with the other limits ──────────────────────────────────

def test_memory_can_lower_the_stage_cap(monkeypatch):
    monkeypatch.setattr(C, "available_ram_mb", lambda: 1_900)
    assert get_parallel_workers(_Params(16), stage="detect",
                                frame_bytes=FRAME_61MP) == 1


def test_memory_never_raises_the_stage_cap(monkeypatch):
    monkeypatch.setattr(C, "available_ram_mb", lambda: 64_000)
    workers = get_parallel_workers(_Params(16), stage="detect",
                                   frame_bytes=FRAME_61MP)
    assert workers == C.STAGE_WORKER_CAPS["detect"]


def test_without_a_frame_size_behaviour_is_unchanged(monkeypatch):
    monkeypatch.setattr(C, "available_ram_mb", lambda: 500)
    assert get_parallel_workers(_Params(16), stage="detect") == \
        C.STAGE_WORKER_CAPS["detect"]


def test_env_override_still_wins(monkeypatch):
    """The sweep must be able to probe above every ceiling, RAM included."""
    monkeypatch.setenv("APEX_MAX_WORKERS", "12")
    monkeypatch.setattr(C, "available_ram_mb", lambda: 500)
    assert get_parallel_workers(_Params(1), stage="detect",
                                frame_bytes=FRAME_61MP) == 12


# ── the decision ledger ────────────────────────────────────────────────────

def test_every_decision_records_why(monkeypatch):
    monkeypatch.setattr(C, "available_ram_mb", lambda: 1_900)
    get_parallel_workers(_Params(16), stage="detect", frame_bytes=FRAME_61MP)

    entry = get_worker_decisions()[-1]
    assert entry["stage"] == "detect"
    assert entry["workers"] == 1
    assert entry["source"] == "config"
    assert entry["configured"] == 16
    assert entry["stage_cap"] == 4
    assert entry["memory_cap"] == 1
    assert entry["available_ram_mb"] == 1_900


def test_env_decisions_are_recorded_too(monkeypatch):
    monkeypatch.setenv("APEX_MAX_WORKERS", "8")
    get_parallel_workers(stage="wcs")
    entry = get_worker_decisions()[-1]
    assert entry["source"] == "env" and entry["workers"] == 8


def test_the_ledger_is_a_snapshot_not_the_live_list():
    get_parallel_workers(stage="detect")
    snapshot = get_worker_decisions()
    snapshot.clear()
    assert get_worker_decisions()


def test_run_manifest_carries_the_decisions(monkeypatch):
    """`auto` is not reproducible unless the manifest says what it chose."""
    from apex.pipeline.runner import RunReport

    monkeypatch.setattr(C, "available_ram_mb", lambda: 1_900)
    get_parallel_workers(_Params(16), stage="detect", frame_bytes=FRAME_61MP)

    manifest = RunReport(mode="cmd").to_dict()
    assert manifest["worker_decisions"][-1]["stage"] == "detect"
    assert manifest["worker_decisions"][-1]["workers"] == 1
