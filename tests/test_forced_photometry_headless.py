"""Headless forced-photometry (Step 7) tests.

Exercises the Qt-free ``apex.analysis.forced_photometry.run_forced_photometry``
and the ``ForcedPhotStep`` pipeline wrapper on a small synthetic field. No Qt,
no real data.
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
from astropy.wcs import WCS

from apex.analysis.forced_photometry import run_forced_photometry
from apex.config.parameters_cmd import read_params
from apex.pipeline.context import RunContext
from apex.pipeline.steps.forcedphot import ForcedPhotStep
from apex.utils.step_paths import (
    step1_dir,
    step4_dir,
    step5_wcs_dir,
    step6_refbuild_dir,
    step7_forced_phot_dir,
)

_EXAMPLE_TOML = Path(__file__).resolve().parents[1] / "parameters.example.toml"

_N_STARS = 40
_SIZE = 256
_FWHM_PX = 5.0


def _make_wcs(size, crval=(150.0, 2.0), scale_arcsec=0.4):
    w = WCS(naxis=2)
    w.wcs.crpix = [size / 2.0, size / 2.0]
    w.wcs.cdelt = [-scale_arcsec / 3600.0, scale_arcsec / 3600.0]
    w.wcs.crval = list(crval)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return w


def _make_frame(size, xs, ys, fluxes, fwhm_px, seed):
    rng = np.random.RandomState(seed)
    img = np.full((size, size), 800.0, dtype=float)
    yy, xx = np.mgrid[0:size, 0:size]
    sigma = fwhm_px / 2.355
    for cx, cy, f in zip(xs, ys, fluxes):
        img += f * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * sigma ** 2))
    img += rng.normal(0, 8.0, (size, size))
    return np.clip(img, 0, 65535).astype(np.float32)


def _build_dataset(tmp_path: Path):
    data_dir = tmp_path / "data"
    result_dir = tmp_path / "result"
    cache_dir = result_dir / "cache"
    for d in (data_dir, result_dir, cache_dir):
        d.mkdir(parents=True, exist_ok=True)

    rng = np.random.RandomState(3)
    xs = rng.uniform(30, _SIZE - 30, _N_STARS)
    ys = rng.uniform(30, _SIZE - 30, _N_STARS)
    fluxes = rng.uniform(3000, 40000, _N_STARS)
    wcs = _make_wcs(_SIZE)
    ra, dec = wcs.all_pix2world(xs, ys, 0)

    names = ["frame_R_001.fits", "frame_R_002.fits"]
    file_path_map = {}
    s4 = step4_dir(result_dir)
    s4.mkdir(parents=True, exist_ok=True)
    for fi, name in enumerate(names):
        img = _make_frame(_SIZE, xs, ys, fluxes, _FWHM_PX, seed=11 + fi)
        hdr = wcs.to_header()
        hdr["FILTER"] = "R"
        hdr["EXPTIME"] = 60.0
        hdr["GAIN"] = 1.5
        hdr["RDNOISE"] = 6.0
        path = data_dir / name
        fits.PrimaryHDU(img, header=hdr).writeto(path, overwrite=True)
        file_path_map[name] = str(path)

        det = pd.DataFrame({
            "det_uid": range(_N_STARS),
            "id": range(1, _N_STARS + 1),
            "x": xs,
            "y": ys,
            "elongation": np.full(_N_STARS, 1.05),
            "roundness1": np.full(_N_STARS, 0.02),
            "roundness2": np.full(_N_STARS, -0.03),
            "sharpness": np.full(_N_STARS, 0.6),
            "peak_adu": fluxes + 800.0,
            "fwhm_px": np.full(_N_STARS, _FWHM_PX),
            "fwhm_status": ["ok"] * _N_STARS,
            "anchor_candidate": [True] * _N_STARS,
            "apcorr_candidate": [True] * _N_STARS,
            "quality_score": np.full(_N_STARS, 0.9),
            "quality_flags": [""] * _N_STARS,
        })
        det.to_csv(s4 / f"detect_{name}.csv", index=False)
        (s4 / f"detect_{name}.json").write_text(
            json.dumps({"fwhm_px": _FWHM_PX, "sky_med": 800.0, "sky_sigma": 8.0}),
            encoding="utf-8",
        )

    step5_wcs_dir(result_dir).mkdir(parents=True, exist_ok=True)
    (step5_wcs_dir(result_dir) / "wcs_solve_summary.csv").write_text(
        "file,solved\n" + "\n".join(f"{n},True" for n in names) + "\n", encoding="utf-8"
    )

    s6 = step6_refbuild_dir(result_dir)
    s6.mkdir(parents=True, exist_ok=True)
    master = pd.DataFrame({
        "master_id": range(1, _N_STARS + 1),
        "source_id": range(1, _N_STARS + 1),
        "ID": range(1, _N_STARS + 1),
        "ra_deg": ra,
        "dec_deg": dec,
        "x_ref": xs,
        "y_ref": ys,
        "crowding_flag": [False] * _N_STARS,
        "neighbor_dist_px": np.full(_N_STARS, 50.0),
    })
    master.to_csv(s6 / "ref_catalog_R.tsv", sep="\t", index=False)

    params = read_params(_EXAMPLE_TOML)
    params.P.data_dir = data_dir
    params.P.result_dir = result_dir
    params.P.cache_dir = cache_dir
    params.P.file_path_map = file_path_map
    params.P.max_workers = 1
    return params, data_dir, result_dir, cache_dir, names


def test_run_forced_photometry_writes_finite_photometry(tmp_path):
    params, data_dir, result_dir, cache_dir, names = _build_dataset(tmp_path)

    summary = run_forced_photometry(
        names, params, data_dir, cache_dir, result_dir=result_dir,
    )

    assert isinstance(summary, dict)
    index_rows = summary["index_rows"]
    assert len(index_rows) == len(names)
    assert all(r["status"] == "ok" for r in index_rows)

    out = step7_forced_phot_dir(result_dir)
    assert (out / "photometry_index.csv").exists()
    assert (out / "apcorr_summary.csv").exists()
    assert (out / "centering_stats.csv").exists()

    for name in names:
        tsv = out / f"photometry_{name}.tsv"
        assert tsv.exists(), f"missing per-frame TSV for {name}"
        df = pd.read_csv(tsv, sep="\t")
        assert len(df) == _N_STARS
        # Non-trivial behavior: magnitudes and fluxes must be present and finite.
        for col in ("flux", "flux_err", "mag_inst", "mag_err", "snr"):
            assert col in df.columns
            assert np.isfinite(df[col]).all(), f"{col} has non-finite values"
        assert (df["flux"] > 0).all()
        assert (df["snr"] > 0).all()
        # All sources detected (master positions match injected stars exactly).
        assert int(df["detected_flag"].sum()) == _N_STARS
        # Aperture correction was actually derived from the growth curve.
        assert (df["apcorr"] > 1.0).all()
        assert df["apcorr_applied"].all()


def test_forced_phot_step_is_complete_after_run(tmp_path):
    params, data_dir, result_dir, cache_dir, names = _build_dataset(tmp_path)

    ctx = RunContext(
        mode="cmd",
        params=params,
        result_dir=result_dir,
        data_dir=data_dir,
        logger=logging.getLogger("test.forcedphot"),
    )

    step = ForcedPhotStep()
    assert step.index == 7 and step.key == "forcedphot"

    # Not complete before running.
    assert step.is_complete(ctx) is False

    # Required input (Step 1 selection.json) must be present for run().
    sel_dir = step1_dir(result_dir)
    sel_dir.mkdir(parents=True, exist_ok=True)
    (sel_dir / "selection.json").write_text(
        json.dumps({"filenames": names, "target": None}),
        encoding="utf-8",
    )

    result = step.run(ctx)
    assert result.ok, result.message
    assert result.status == "ok"

    # Now is_complete() detects photometry_index.csv.
    assert step.is_complete(ctx) is True


def test_forced_phot_step_blocked_without_selection(tmp_path):
    params, data_dir, result_dir, cache_dir, names = _build_dataset(tmp_path)
    ctx = RunContext(
        mode="cmd",
        params=params,
        result_dir=result_dir,
        data_dir=data_dir,
        logger=logging.getLogger("test.forcedphot"),
    )
    result = ForcedPhotStep().run(ctx)
    assert result.status == "blocked"
