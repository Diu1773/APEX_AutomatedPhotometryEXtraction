"""Headless master-catalog-build (Step 6) tests.

Exercises the Qt-free ``apex.analysis.refbuild.run_refbuild`` and the
``RefBuildStep`` pipeline wrapper on a small synthetic field. No Qt, no real data.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("astropy")
pytest.importorskip("pandas")

import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from apex.analysis.refbuild import run_refbuild
from apex.config.parameters_cmd import read_params
from apex.pipeline.context import RunContext
from apex.pipeline.steps.refbuild import RefBuildStep
from apex.utils.step_paths import (
    step1_dir,
    step4_dir,
    step5_wcs_dir,
    step6_refbuild_dir,
)

_EXAMPLE_TOML = Path(__file__).resolve().parents[1] / "parameters.example.toml"

_RA0, _DEC0 = 200.0, -30.0
_PIX = 0.0005  # deg/pix


def _make_wcs(nx=512, ny=512):
    w = WCS(naxis=2)
    w.wcs.crpix = [nx / 2.0, ny / 2.0]
    w.wcs.cdelt = [-_PIX, _PIX]
    w.wcs.crval = [_RA0, _DEC0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return w


def _build_dataset(tmp_path: Path):
    data_dir = tmp_path / "data"
    result_dir = tmp_path / "result"
    cache_dir = result_dir / "cache"
    s4 = step4_dir(result_dir)
    s5 = step5_wcs_dir(result_dir)
    for d in (data_dir, result_dir, cache_dir, s4, s5):
        d.mkdir(parents=True, exist_ok=True)

    nx = ny = 512
    wcs = _make_wcs(nx, ny)
    hdr = wcs.to_header()
    hdr["EXPTIME"] = 30.0
    hdr["FILTER"] = "V"
    hdr["GAIN"] = 1.5
    hdr["RDNOISE"] = 8.0

    # Gaia catalog within the FOV.
    n = 120
    rng = np.random.default_rng(7)
    half = 0.5 * nx * _PIX
    ra = _RA0 + rng.uniform(-half, half, n) / np.cos(np.deg2rad(_DEC0))
    dec = _DEC0 + rng.uniform(-half, half, n)
    gaia = Table()
    gaia["source_id"] = (np.arange(n) + 4000000000000000000).astype(np.int64)
    gaia["ra"] = ra
    gaia["dec"] = dec
    gaia["phot_g_mean_mag"] = rng.uniform(12.0, 19.0, n)
    gaia["phot_bp_mean_mag"] = gaia["phot_g_mean_mag"] + 0.4
    gaia["phot_rp_mean_mag"] = gaia["phot_g_mean_mag"] - 0.4
    gaia.write(str(s5 / "gaia_fov.ecsv"), format="ascii.ecsv", overwrite=True)

    sky = np.column_stack([np.asarray(gaia["ra"]), np.asarray(gaia["dec"])])
    xy = wcs.all_world2pix(sky, 0)
    on = (xy[:, 0] >= 5) & (xy[:, 0] < nx - 5) & (xy[:, 1] >= 5) & (xy[:, 1] < ny - 5)
    xy = xy[on]

    names = ["cluster_V_20250115_001.fits", "cluster_V_20250115_002.fits"]
    file_path_map = {}
    summary_rows = []
    jit = np.random.default_rng(11)
    for fi, name in enumerate(names):
        img = (jit.standard_normal((ny, nx)) * 5.0 + 100.0).astype(np.float32)
        path = data_dir / name
        fits.PrimaryHDU(img, header=hdr).writeto(path, overwrite=True)
        file_path_map[name] = str(path)

        det_xy = xy + jit.normal(0.0, 0.05, xy.shape)
        keep = np.ones(len(det_xy), dtype=bool)
        if fi == 0:
            keep[:3] = False
        else:
            keep[-3:] = False
        dxy = det_xy[keep]
        ndet = len(dxy)

        det = pd.DataFrame({
            "det_uid": np.arange(ndet, dtype=int),
            "x": dxy[:, 0],
            "y": dxy[:, 1],
            "elongation": np.full(ndet, 1.1),
            "roundness": np.full(ndet, 0.05),
            "sharpness": np.full(ndet, 0.5),
            "dao_flux": jit.uniform(1000.0, 50000.0, ndet),
            "dao_peak": jit.uniform(500.0, 30000.0, ndet),
            "peak_adu": jit.uniform(500.0, 30000.0, ndet),
            "fwhm_px": np.full(ndet, 3.2),
            "anchor_candidate": np.ones(ndet, dtype=bool),
            "apcorr_candidate": np.ones(ndet, dtype=bool),
        })
        det.to_csv(s4 / f"detect_{name}.csv", index=False)
        (s4 / f"detect_{name}.json").write_text(
            json.dumps({
                "filter": "V", "fwhm_px": 3.2, "n_sources": int(ndet),
                "sat_star_count": int(fi), "median_elongation": 1.1,
                "median_roundness": 0.05, "sky_med": 100.0, "sky_sigma": 5.0,
            }),
            encoding="utf-8",
        )
        summary_rows.append({
            "file": name, "wcs_ok": True, "match_n": ndet, "n_match": ndet,
            "n_catalog_in_fov": n, "match_rate": 0.95, "match_rate_cat": 0.9,
            "match_rate_eff": 0.95, "resid_med": 0.3, "resid_max": 0.8,
            "rms_px": 0.2, "wcs_qc_pass": True, "wcs_qc_reason": "",
            "gaia_source": "fixture",
        })

    pd.DataFrame(summary_rows).to_csv(s5 / "wcs_solve_summary.csv", index=False)

    params = read_params(_EXAMPLE_TOML)
    params.P.data_dir = data_dir
    params.P.result_dir = result_dir
    params.P.cache_dir = cache_dir
    params.P.file_path_map = file_path_map
    return params, data_dir, result_dir, cache_dir, names


def _run_kwargs(params, names):
    return dict(
        params=params,
        data_dir=params.P.data_dir,
        result_dir=params.P.result_dir,
        cache_dir=params.P.cache_dir,
        file_list=names,
        ref_filter="V",
        sat_drop_pct=20.0,
        elong_drop_pct=20.0,
        ref_cat_max_sources=0,
        ref_cat_min_sources=50,
        ref_cat_max_elong=1.5,
        ref_cat_max_abs_round=0.4,
        ref_cat_sharp_min=0.2,
        ref_cat_sharp_max=1.0,
        ref_cat_min_peak_adu=0.0,
        wcs_match_radius_arcsec=2.0,
        wcs_min_match_rate=0.2,
        wcs_min_match_n=10,
        wcs_max_sep_med_arcsec=1.5,
        wcs_max_sep_p90_arcsec=2.5,
        wcs_max_dup_rate=0.5,
        ref_per_date=True,
        ref_master_union=True,
        ref_union_min_frames=1,
        ref_build_mode="hybrid",
        gaia_mag_limit=18.0,
    )


def test_run_refbuild_builds_master_with_finite_coords(tmp_path):
    params, data_dir, result_dir, cache_dir, names = _build_dataset(tmp_path)

    summary = run_refbuild(**_run_kwargs(params, names))

    assert isinstance(summary, dict)
    assert summary["ref_frame"] in names
    assert summary["ref_filter"] == "V"
    assert summary["n_sources"] > 0

    out = step6_refbuild_dir(result_dir)
    cat = out / "ref_catalog.tsv"
    assert cat.exists() and cat.stat().st_size > 0
    assert (out / "ref_frame_stats.csv").exists()
    assert (out / "ref_build_meta.json").exists()

    df = pd.read_csv(cat, sep="\t")
    assert len(df) == summary["n_sources"]
    for col in ("ID", "source_id", "ra_deg", "dec_deg", "x_ref", "y_ref"):
        assert col in df.columns
    # Coordinates must be finite.
    assert np.isfinite(pd.to_numeric(df["ra_deg"], errors="coerce")).all()
    assert np.isfinite(pd.to_numeric(df["dec_deg"], errors="coerce")).all()
    assert np.isfinite(pd.to_numeric(df["x_ref"], errors="coerce")).all()
    assert np.isfinite(pd.to_numeric(df["y_ref"], errors="coerce")).all()
    # Union build merged the two frames into a single, deduped, dense ID set.
    assert df["ID"].is_unique
    # Some sources matched Gaia (positive hybrid IDs).
    sid = pd.to_numeric(df["source_id"], errors="coerce")
    assert (sid > 0).any()


def test_refbuild_step_is_complete_after_run(tmp_path):
    params, data_dir, result_dir, cache_dir, names = _build_dataset(tmp_path)

    ctx = RunContext(
        mode="cmd",
        params=params,
        result_dir=result_dir,
        data_dir=data_dir,
        logger=logging.getLogger("test.refbuild"),
    )

    step = RefBuildStep()
    assert step.index == 6 and step.key == "refbuild"

    # Not complete before running.
    assert step.is_complete(ctx) is False

    sel_dir = step1_dir(result_dir)
    sel_dir.mkdir(parents=True, exist_ok=True)
    (sel_dir / "selection.json").write_text(
        json.dumps({"filenames": names, "target": None}),
        encoding="utf-8",
    )

    result = step.run(ctx)
    assert result.ok, result.message
    assert result.status == "ok"

    # Now is_complete() detects ref_catalog.tsv.
    assert step.is_complete(ctx) is True


def test_refbuild_step_blocked_without_selection(tmp_path):
    params, data_dir, result_dir, cache_dir, names = _build_dataset(tmp_path)
    ctx = RunContext(
        mode="cmd",
        params=params,
        result_dir=result_dir,
        data_dir=data_dir,
        logger=logging.getLogger("test.refbuild"),
    )
    result = RefBuildStep().run(ctx)
    assert result.status == "blocked"
