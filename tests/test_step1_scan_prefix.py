from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("astropy")
pytest.importorskip("pandas")

import pandas as pd

import apex.core.file_manager as file_manager_module
from apex.core.file_manager import FileManager
from apex.utils.step_paths import step1_dir


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


def test_read_headers_reuses_manifest_validated_cache(tmp_path, monkeypatch):
    params = _params(tmp_path, "")
    source = Path(params.P.data_dir) / "raw.fits"
    source.write_bytes(b"not opened when cache is valid")
    manager = FileManager(params)
    manager.scan_files()

    output_dir = step1_dir(params.P.result_dir)
    pd.DataFrame(
        {
            "Filename": [source.name],
            "DATE-OBS": ["2026-01-01T00:00:00"],
            "FILTER": ["V"],
        }
    ).to_csv(output_dir / "headers.csv", index=False)
    (output_dir / "headers_cache_manifest.json").write_text(
        json.dumps(manager._header_cache_manifest()),
        encoding="utf-8",
    )

    def fail_open(*_args, **_kwargs):
        raise AssertionError("valid header cache should avoid FITS I/O")

    monkeypatch.setattr(file_manager_module.fits, "open", fail_open)

    cached = manager.read_headers()

    assert cached["Filename"].tolist() == [source.name]
    assert manager.df_headers is cached
