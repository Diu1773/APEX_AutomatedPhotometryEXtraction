from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("PyQt5")
from PyQt5.QtCore import QCoreApplication

from apex.gui.workflow.cmd.step12_isochrone_model import IsoCacheWorker


def test_iso_cache_worker_writes_complete_atomic_sidecar(tmp_path):
    app = QCoreApplication.instance() or QCoreApplication([])
    source = tmp_path / "iso.dat"
    source.write_text("# test\n1 2 3\n4 5 6\n", encoding="utf-8")

    worker = IsoCacheWorker(source)
    worker.start()

    assert worker.wait(10000)
    cache_path = source.with_suffix(source.suffix + ".apex.npy")
    assert np.array_equal(np.load(cache_path), np.array([[1, 2, 3], [4, 5, 6]]))
    assert list(tmp_path.glob("*.tmp")) == []
    assert app is not None


def test_iso_cache_worker_cancelled_before_start_does_not_publish(tmp_path):
    app = QCoreApplication.instance() or QCoreApplication([])
    source = tmp_path / "iso.dat"
    source.write_text("# test\n1 2 3\n", encoding="utf-8")

    worker = IsoCacheWorker(source)
    worker.cancel()
    worker.start()

    assert worker.wait(10000)
    assert not source.with_suffix(source.suffix + ".apex.npy").exists()
    assert list(tmp_path.glob("*.tmp")) == []
    assert app is not None
