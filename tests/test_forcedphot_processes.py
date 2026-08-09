"""Forced photometry may run in processes instead of threads (O2 follow-up).

Threads cannot help this stage: measured on NGC 6811 it stays on one core
however many it gets (83 % CPU at 1 worker, 109 % at 12) because photutils'
per-aperture statistics hold the interpreter lock, so 12 threads cost 4.2x wall
time and buy nothing. Separate processes each own an interpreter lock.

The parity that matters is checked end-to-end elsewhere — the same 21 frames
produce byte-identical `photometry_*.tsv` either way, and Step 7 alone drops
from 244.2 s to 69.6 s. These tests pin the plumbing that makes that possible:
the payload must survive pickling, the split must not lose or duplicate a
frame, and a machine that cannot spawn processes must still get its photometry.
"""

from __future__ import annotations

import pickle

import pandas as pd
import pytest

from apex.analysis import forced_photometry as FP


def _tasks(n_per_filter: int, filters=("B", "V", "R")) -> list:
    master = pd.DataFrame({"source_id": [1, 2, 3]})
    return [(filt, f"{filt}_{i:03d}.fits", master)
            for filt in filters for i in range(n_per_filter)]


# ── splitting ──────────────────────────────────────────────────────────────

@pytest.mark.parametrize("n_chunks", [2, 3, 4, 7, 12])
def test_every_frame_appears_exactly_once(n_chunks):
    tasks = _tasks(7)
    chunks = FP._chunk_frames(tasks, n_chunks)
    flat = [t[1] for chunk in chunks for t in chunk]
    assert sorted(flat) == sorted(t[1] for t in tasks)
    assert len(flat) == len(set(flat)), "a frame was duplicated across chunks"


def test_a_chunk_holds_one_filter():
    """Each child re-loads the master catalogue, so mixing filters costs extra."""
    for chunk in FP._chunk_frames(_tasks(5), 4):
        assert len({t[0] for t in chunk}) == 1


def test_a_single_frame_is_not_split():
    tasks = _tasks(1, filters=("V",))
    assert FP._chunk_frames(tasks, 8) == [tasks]


def test_one_chunk_requested_returns_everything():
    tasks = _tasks(3)
    assert FP._chunk_frames(tasks, 1) == [tasks]


def test_no_empty_chunks():
    for n in range(1, 13):
        assert all(chunk for chunk in FP._chunk_frames(_tasks(4), n))


# ── the process boundary ───────────────────────────────────────────────────

def test_the_payload_survives_pickling():
    """Everything handed to a worker process must cross the boundary."""
    payload = {
        "files": ["a.fits", "b.fits"],
        "params": None,                    # real params pickle at ~14 kB
        "data_dir": "E:/data",
        "cache_dir": "E:/result/cache",
        "result_dir": "E:/result",
        "output_dir": "E:/result/step7_forced_phot",
    }
    assert pickle.loads(pickle.dumps(payload)) == payload


def test_the_worker_entry_point_is_importable_by_name():
    """A closure cannot be pickled; the entry point has to be module-level."""
    from apex.analysis.forced_photometry import _forced_phot_chunk

    assert _forced_phot_chunk.__module__ == "apex.analysis.forced_photometry"
    assert pickle.loads(pickle.dumps(_forced_phot_chunk)) is _forced_phot_chunk


def test_a_single_worker_declines_the_process_path(monkeypatch):
    """One process is strictly worse than one thread — it pays spawn for nothing."""
    monkeypatch.setattr(FP, "get_parallel_workers", lambda *a, **k: 1)
    assert FP._run_tasks_in_processes(
        _tasks(2), None, "d", "c", "r", "o") is None


def test_a_pool_failure_falls_back_instead_of_losing_data(monkeypatch):
    """Frozen builds and restricted environments must still get photometry."""
    monkeypatch.setattr(FP, "get_parallel_workers", lambda *a, **k: 4)

    said = []
    result = FP._run_tasks_in_processes(
        _tasks(2), object(), "d", "c", "r", "o", log=said.append)

    assert result is None, "an unusable pool must return None, not raise"
    assert any("thread" in m for m in said), "the fallback should be logged"
