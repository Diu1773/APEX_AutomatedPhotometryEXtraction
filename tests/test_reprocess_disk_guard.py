"""Disk guards in scripts/reprocess_batch.py.

The old guard only looked at E: and only between targets, so a single target
larger than the floor (YZBoo_2n is 42 GB) could fill the drive mid-run, and C:
— routinely the tighter drive — was never checked at all.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).absolute().parents[1]


def _load():
    spec = importlib.util.spec_from_file_location(
        "reprocess_batch", REPO / "scripts" / "reprocess_batch.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["reprocess_batch"] = module
    spec.loader.exec_module(module)
    return module


rb = _load()


def test_both_drives_are_guarded(monkeypatch):
    seen = []

    def fake_free(drive):
        seen.append(drive)
        return 999.0

    monkeypatch.setattr(rb, "free_gb", fake_free)
    rb.check_disk("unit test")
    assert rb.OUT_DRIVE in seen, "the output drive must be checked"
    assert rb.REPO_DRIVE in seen, "the repo/venv drive must be checked too"


def test_low_output_drive_raises(monkeypatch):
    monkeypatch.setattr(rb, "free_gb",
                        lambda d: 1.0 if d == rb.OUT_DRIVE else 999.0)
    with pytest.raises(rb.DiskFull, match=r"free 1 GB"):
        rb.check_disk("after Step 0")


def test_low_repo_drive_raises(monkeypatch):
    """C: used to be invisible to the guard."""
    monkeypatch.setattr(rb, "free_gb",
                        lambda d: 999.0 if d == rb.OUT_DRIVE else 1.0)
    with pytest.raises(rb.DiskFull, match=r"free 1 GB"):
        rb.check_disk("before Steps 1-7")


def test_guard_message_names_the_step(monkeypatch):
    monkeypatch.setattr(rb, "free_gb", lambda d: 0.0)
    with pytest.raises(rb.DiskFull, match="before Step 10"):
        rb.check_disk("before Step 10")


def test_step0_refuses_a_target_that_cannot_fit(monkeypatch, tmp_path):
    """Estimate + floor must exceed free space -> refuse before the 3-hour run."""
    monkeypatch.setattr(rb, "REPROCESS", tmp_path)
    monkeypatch.setattr(rb, "estimate_step0_gb", lambda raw: 60.0)
    monkeypatch.setattr(rb, "free_gb", lambda d: 30.0)
    monkeypatch.setattr(rb, "log", lambda msg: None)
    with pytest.raises(rb.DiskFull, match="needs ~60 GB"):
        rb.run_step0("FAKE", str(tmp_path / "raw"))


def test_step0_proceeds_when_it_fits(monkeypatch, tmp_path):
    monkeypatch.setattr(rb, "REPROCESS", tmp_path)
    monkeypatch.setattr(rb, "estimate_step0_gb", lambda raw: 20.0)
    monkeypatch.setattr(rb, "free_gb", lambda d: 200.0)
    monkeypatch.setattr(rb, "log", lambda msg: None)
    called = {}

    class _Ok:
        returncode = 0

    def _run(*_args, **_kwargs):
        called["ran"] = True
        return _Ok()

    monkeypatch.setattr(rb.subprocess, "run", _run)
    rb.run_step0("FAKE", str(tmp_path / "raw"))
    assert called.get("ran")


def test_estimate_reads_the_raw_tree(tmp_path):
    (tmp_path / "n1").mkdir()
    (tmp_path / "n1" / "a.fit").write_bytes(b"x" * 2_000_000)
    (tmp_path / "n1" / "b.fits").write_bytes(b"x" * 3_000_000)
    (tmp_path / "n1" / "notes.txt").write_bytes(b"x" * 9_000_000)   # ignored
    assert rb.estimate_step0_gb(str(tmp_path)) == pytest.approx(0.005, rel=0.01)
    assert rb.estimate_step0_gb(str(tmp_path / "missing")) == 0.0
