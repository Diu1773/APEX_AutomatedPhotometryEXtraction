"""Headless Step 2 (crop) parity + skip-marker tests.

Covers:
  * ``apex.analysis.crop.run_crop`` produces the exact cropped pixel data and
    FITS header keywords the GUI worker produced (the GUI now delegates to this
    same function, so this also pins the shared contract).
  * The cropped slice equals the reference slice computed independently.
  * ``CropStep`` writes a ``skipped: true`` marker when no crop is configured
    and returns OK; and crops + writes a rectangle marker when configured.
"""

from __future__ import annotations

import json
import types
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from apex.analysis.crop import run_crop, _is_crop_cache_valid
from apex.utils.step_paths import step2_cropped_dir, crop_rect_path, step1_dir


def _make_params(tmp_path: Path):
    data_dir = tmp_path / "data"
    result_dir = tmp_path / "result"
    data_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)
    P = types.SimpleNamespace(
        data_dir=data_dir,
        result_dir=result_dir,
        cache_dir=result_dir / "cache",
        max_workers=1,
        parallel_max_workers=1,
        parallel_mode="thread",
    )
    P.cache_dir.mkdir(parents=True, exist_ok=True)

    class _Params:
        def __init__(self, P):
            self.P = P

        def get_file_path(self, filename):
            return Path(self.P.data_dir) / filename

    return _Params(P), data_dir, result_dir


def _write_fits(path: Path, data: np.ndarray, **hdr_kw):
    hdr = fits.Header()
    hdr["OBJECT"] = "TEST"
    for k, v in hdr_kw.items():
        hdr[k] = v
    fits.PrimaryHDU(data=data, header=hdr).writeto(path, overwrite=True)


def test_run_crop_pixels_and_headers(tmp_path):
    params, data_dir, result_dir = _make_params(tmp_path)
    rng = np.random.default_rng(0)
    files = []
    originals = {}
    for i in range(3):
        d = rng.normal(1000, 50, (300, 400)).astype(np.float32)
        fn = f"frame_{i:03d}.fits"
        _write_fits(data_dir / fn, d)
        files.append(fn)
        originals[fn] = d

    x0, y0, x1, y1 = 50, 40, 350, 260
    n_total, errors = run_crop(files, params, x0, y0, x1, y1)
    assert n_total == 3
    assert errors == []

    cropped_dir = step2_cropped_dir(result_dir)
    for fn in files:
        cpath = cropped_dir / fn
        assert cpath.exists()
        with fits.open(cpath) as hdul:
            chdr = hdul[0].header
            cdata = hdul[0].data
        # Pixel parity: exactly the reference slice.
        ref = originals[fn][y0:y1, x0:x1]
        np.testing.assert_array_equal(cdata, ref)
        # Header keywords match the GUI's documented schema.
        assert bool(chdr["CROPPED"]) is True
        assert int(chdr["CROP_X0"]) == x0
        assert int(chdr["CROP_Y0"]) == y0
        assert int(chdr["CROP_X1"]) == x1
        assert int(chdr["CROP_Y1"]) == y1
        assert int(chdr["NAXIS1"]) == (x1 - x0)
        assert int(chdr["NAXIS2"]) == (y1 - y0)


def test_run_crop_cache_reuse(tmp_path):
    params, data_dir, result_dir = _make_params(tmp_path)
    d = np.arange(300 * 400, dtype=np.float32).reshape(300, 400)
    _write_fits(data_dir / "a.fits", d)
    rect = (50, 40, 350, 260)
    run_crop(["a.fits"], params, *rect)
    cpath = step2_cropped_dir(result_dir) / "a.fits"
    orig = data_dir / "a.fits"
    assert _is_crop_cache_valid(orig, cpath, *rect)
    # A different rectangle invalidates the cache.
    assert not _is_crop_cache_valid(orig, cpath, 0, 0, 100, 100)


def _run_crop_step(params, result_dir, data_dir, tmp_path):
    from apex.pipeline.context import RunContext
    from apex.pipeline.steps.crop import CropStep
    import logging

    ctx = RunContext(
        mode="cmd",
        params=params,
        result_dir=Path(result_dir),
        data_dir=Path(data_dir),
        logger=logging.getLogger("test"),
    )
    return CropStep().run(ctx), ctx


def _write_selection(result_dir, filenames):
    s1 = step1_dir(result_dir)
    s1.mkdir(parents=True, exist_ok=True)
    (s1 / "selection.json").write_text(
        json.dumps({"filenames": list(filenames)}), encoding="utf-8"
    )


def test_cropstep_skip_marker_when_unconfigured(tmp_path):
    params, data_dir, result_dir = _make_params(tmp_path)
    # crop disabled (default)
    params.P.crop_enable = False
    params.P.crop_x0 = params.P.crop_y0 = params.P.crop_x1 = params.P.crop_y1 = None
    _write_fits(data_dir / "a.fits", np.zeros((100, 100), np.float32))
    _write_selection(result_dir, ["a.fits"])

    result, _ = _run_crop_step(params, result_dir, data_dir, tmp_path)
    assert result.status == "ok"
    marker = crop_rect_path(result_dir)
    assert marker.exists()
    data = json.loads(marker.read_text())
    assert data.get("skipped") is True
    # No cropped FITS produced when skipped.
    assert not (step2_cropped_dir(result_dir) / "a.fits").exists()


def test_cropstep_config_driven_crop(tmp_path):
    params, data_dir, result_dir = _make_params(tmp_path)
    rng = np.random.default_rng(1)
    d = rng.normal(500, 20, (300, 400)).astype(np.float32)
    _write_fits(data_dir / "a.fits", d)
    _write_selection(result_dir, ["a.fits"])

    params.P.crop_enable = True
    params.P.crop_x0 = 60
    params.P.crop_y0 = 30
    params.P.crop_x1 = 360
    params.P.crop_y1 = 250

    result, _ = _run_crop_step(params, result_dir, data_dir, tmp_path)
    assert result.status == "ok"
    marker = crop_rect_path(result_dir)
    md = json.loads(marker.read_text())
    assert md["skipped"] is False
    assert (md["x0"], md["y0"], md["x1"], md["y1"]) == (60, 30, 360, 250)

    cpath = step2_cropped_dir(result_dir) / "a.fits"
    assert cpath.exists()
    with fits.open(cpath) as hdul:
        np.testing.assert_array_equal(hdul[0].data, d[30:250, 60:360])
