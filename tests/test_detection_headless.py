"""Headless source-detection (Step 4) tests.

Exercises the Qt-free ``apex.analysis.detection.run_detection`` and the
``DetectStep`` pipeline wrapper on a small synthetic frame. No Qt, no real data.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("astropy")
pytest.importorskip("photutils")
pytest.importorskip("pandas")

import pandas as pd
from astropy.io import fits

from apex.analysis.detection import _bounded_frame_fwhm, run_detection
from apex.config.parameters_cmd import read_params
from apex.pipeline.context import RunContext
from apex.pipeline.steps.detect import DetectStep
from apex.utils.step_paths import step1_dir, step4_dir

_EXAMPLE_TOML = Path(__file__).resolve().parents[1] / "parameters.example.toml"

# Columns the headless CSV is documented to contain (pre-quality + quality).
_EXPECTED_CSV_COLUMNS = {
    "det_uid", "id", "x", "y", "sharpness",
    "roundness1", "roundness2", "dao_peak", "dao_flux",
    "peak_adu", "fwhm_px", "fwhm_status", "source_type",
    "quality_score", "anchor_candidate", "apcorr_candidate",
    "epsf_candidate", "psf_seed_candidate",
}


def _make_synthetic_frame(n_stars=20, fwhm_px=6.0, background=1000.0,
                          noise_level=10.0, seed=7, size=256):
    rng = np.random.RandomState(seed)
    h = w = size
    img = np.full((h, w), background, dtype=float)
    y, x = np.mgrid[0:h, 0:w]
    sigma = fwhm_px / 2.355
    for _ in range(n_stars):
        cx = rng.uniform(20, w - 20)
        cy = rng.uniform(20, h - 20)
        flux = rng.uniform(2000, 12000)
        img += flux * np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * sigma ** 2))
    img += rng.normal(0, noise_level, (h, w))
    return np.clip(img, 0, 65535).astype(np.float32)


def _write_frame(path: Path, filt="R") -> None:
    hdu = fits.PrimaryHDU(_make_synthetic_frame())
    hdu.header["FILTER"] = filt
    hdu.header["EXPTIME"] = 60.0
    hdu.header["DATE-OBS"] = "2026-06-20T00:00:00"
    hdu.writeto(path, overwrite=True)


def _build_params(tmp_path: Path, name: str):
    data_dir = tmp_path / "data"
    result_dir = tmp_path / "result"
    cache_dir = result_dir / "cache"
    data_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    fits_path = data_dir / name
    _write_frame(fits_path)

    params = read_params(_EXAMPLE_TOML)
    params.P.data_dir = data_dir
    params.P.result_dir = result_dir
    params.P.cache_dir = cache_dir
    params.P.file_path_map = {name: str(fits_path)}
    # segmentation engine is always available (no optional 'sep' dependency)
    params.P.detect_engine = "segm"
    params.P.max_workers = 1
    return params, data_dir, result_dir, cache_dir, fits_path


def test_frame_fwhm_rejects_compact_artifacts_and_broad_blends():
    value, count = _bounded_frame_fwhm(
        [8.68, 8.67, 21.68, 35.44, 8.58, 16.24, 8.43, 9.19, 8.94,
         7.71, 16.92, 2.05, 2.26, 1.00, 0.99, 1.99, 1.00, 1.09],
        minimum_px=3.5,
        maximum_px=12.0,
        fallback_px=6.0,
    )

    assert count == 7
    assert value == pytest.approx(8.67, abs=0.15)


def test_frame_fwhm_keeps_dominant_bad_seeing_mode_above_qc_maximum():
    value, count = _bounded_frame_fwhm(
        [
            1.00, 1.01, 1.99, 2.00,
            3.998, 4.992, 6.998,
            10.128, 10.134, 10.200, 10.259, 10.298, 10.321,
            10.537, 10.726, 10.867, 10.868, 11.090,
            12.940, 13.044, 14.849, 18.141, 18.662,
        ],
        minimum_px=3.0,
        maximum_px=10.0,
        fallback_px=6.0,
    )

    assert value > 10.0
    assert value == pytest.approx(10.4, abs=0.2)
    assert count >= 10


def test_run_detection_writes_valid_outputs(tmp_path):
    name = "synth_R.fits"
    params, data_dir, result_dir, cache_dir, _ = _build_params(tmp_path, name)

    summary = run_detection(
        [name], params, data_dir, cache_dir,
        use_cropped=False, filter_sigma_map={},
    )

    # Returned summary
    assert isinstance(summary, dict)
    assert summary.get("total_files") == 1
    assert summary.get("total_sources", 0) > 0

    s4 = step4_dir(result_dir)
    json_path = s4 / f"detect_{name}.json"
    csv_path = s4 / f"detect_{name}.csv"
    assert json_path.exists(), "detection JSON not written"
    assert csv_path.exists(), "detection CSV not written"

    meta = json.loads(json_path.read_text())
    assert meta["n_sources"] > 0
    assert meta["filter"] == "R"
    assert meta["detect_engine"] == "segm"
    assert "detect_method" in meta
    assert np.isfinite(meta["fwhm_px"])
    assert meta["fwhm_estimator"] == "lower_bounded_mad_median"
    assert "fwhm_within_config" in meta
    assert meta["fwhm_qc_min_px"] == pytest.approx(params.P.fwhm_px_min)
    assert meta["fwhm_qc_max_px"] == pytest.approx(params.P.fwhm_px_max)

    df = pd.read_csv(csv_path)
    assert len(df) == meta["n_sources"]
    missing = _EXPECTED_CSV_COLUMNS - set(df.columns)
    assert not missing, f"CSV missing documented columns: {missing}"
    # id is 1-indexed and contiguous; det_uid is 0-indexed.
    assert list(df["id"]) == list(range(1, len(df) + 1))
    assert list(df["det_uid"]) == list(range(len(df)))


def test_filter_sigma_override_is_applied(tmp_path):
    """The per-filter detection sigma override must reach sigma_used.

    Guards the filter_sigma_map lookup path (which the empty-map smoke test does
    not exercise): with no override, sigma_used == base detect_sigma; with an
    override keyed by the frame's own FILTER value, sigma_used == the override.
    """
    name = "synth_R.fits"
    params, data_dir, result_dir, cache_dir, _ = _build_params(tmp_path, name)
    base_sigma = float(params.P.detect_sigma)
    override = base_sigma + 1.0  # distinct from the base so the assertion is real

    # Baseline: empty map -> base sigma.
    run_detection([name], params, data_dir, cache_dir, use_cropped=False,
                  filter_sigma_map={})
    meta_base = json.loads((step4_dir(result_dir) / f"detect_{name}.json").read_text())
    assert meta_base["sigma_used"] == pytest.approx(base_sigma)

    # Override keyed by the frame's FILTER ("R") -> override applied. Force a
    # rerun by clearing the step4 output first.
    for p in step4_dir(result_dir).glob("detect_*"):
        p.unlink()
    for p in cache_dir.glob("detect_*"):
        p.unlink()
    run_detection([name], params, data_dir, cache_dir, use_cropped=False,
                  filter_sigma_map={"R": override})
    meta_ovr = json.loads((step4_dir(result_dir) / f"detect_{name}.json").read_text())
    assert meta_ovr["sigma_used"] == pytest.approx(override)


def test_detect_step_is_complete_after_run(tmp_path):
    name = "synth_R.fits"
    params, data_dir, result_dir, cache_dir, _ = _build_params(tmp_path, name)

    ctx = RunContext(
        mode="cmd",
        params=params,
        result_dir=result_dir,
        data_dir=data_dir,
        logger=logging.getLogger("test.detect"),
    )

    step = DetectStep()
    assert step.index == 4 and step.key == "detect"

    # Not complete before running.
    assert step.is_complete(ctx) is False

    # Required input (Step 1 selection.json) must be present for run().
    sel_dir = step1_dir(result_dir)
    sel_dir.mkdir(parents=True, exist_ok=True)
    (sel_dir / "selection.json").write_text(
        json.dumps({"filenames": [name], "target": None}),
        encoding="utf-8",
    )

    result = step.run(ctx)
    assert result.ok, result.message
    assert result.status == "ok"

    # Now is_complete() detects the produced detect_*.json.
    assert step.is_complete(ctx) is True
