from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("astropy")
pytest.importorskip("pandas")

from apex.core.file_manager import FileManager


def _params(tmp_path: Path, prefix: str):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    return SimpleNamespace(
        P=SimpleNamespace(
            data_dir=data_dir,
            filename_prefix=prefix,
            result_dir=tmp_path / "result",
        )
    )


def test_scan_files_empty_prefix_matches_all_fits(tmp_path):
    params = _params(tmp_path, "")
    data_dir = Path(params.P.data_dir)
    for name in ("raw_b.fits", "pp_raw_a.fit", "calibrated.fit.fz", "notes.txt"):
        (data_dir / name).touch()

    files = FileManager(params).scan_files()

    assert files == ["calibrated.fit.fz", "pp_raw_a.fit", "raw_b.fits"]


def test_scan_files_nonempty_prefix_filters_fits(tmp_path):
    params = _params(tmp_path, "pp_")
    data_dir = Path(params.P.data_dir)
    for name in ("raw_b.fits", "pp_raw_a.fit", "pp_notes.txt"):
        (data_dir / name).touch()

    files = FileManager(params).scan_files()

    assert files == ["pp_raw_a.fit"]
