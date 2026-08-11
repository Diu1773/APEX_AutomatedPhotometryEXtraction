"""Headless Step 5 (WCS) parity/smoke tests — INTERNAL engine only.

External solvers (ASTAP, astrometry.net/solve-field) are not available in the
test environment, so only the internal Python engine path is parity-tested here.

Parity strategy (internal engine):
  * Build a synthetic field with a known WCS (same fixture style as
    tests/test_wcs_solver_internal.py): inject Gaia stars at true sky positions,
    write a synthetic FITS + a Step 4 detection CSV, and pre-seed
    step5_wcs/gaia_fov.ecsv so no network query is needed.
  * BASELINE: drive the CURRENT GUI worker (apex.gui ... InternalWcsWorker) under
    a QCoreApplication with the internal engine; capture the FITS WCS headers and
    the wcs_solve_summary.csv it produces.
  * After extraction: re-run apex.analysis.wcs_solve.run_wcs_solve(engine=
    "internal") on an identical copy of the inputs; diff CRVAL/CRPIX/CD within
    rtol=1e-9 and the summary CSV equal modulo paths/timestamps.

If PyQt5 is unavailable the GUI-baseline half is skipped and the test still
exercises the headless orchestrator end-to-end (so the internal path is always
covered at least at the orchestration level).
"""

from __future__ import annotations

import json
import types
from pathlib import Path

import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

pytest.importorskip("scipy")

from apex.utils.step_paths import step1_dir, step4_dir, step5_wcs_dir


# ── synthetic field ──────────────────────────────────────────────────────────

APPROX_RA = 250.4218
APPROX_DEC = 36.4613
SCALE = 0.42
SHAPE = (1024, 1024)


def _synthetic_wcs(rotation_deg=12.0, crpix_offset=(0.0, 0.0)) -> WCS:
    ny, nx = SHAPE
    scale_deg = SCALE / 3600.0
    theta = np.deg2rad(rotation_deg)
    c, s = np.cos(theta), np.sin(theta)
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crpix = [nx / 2.0 + crpix_offset[0], ny / 2.0 + crpix_offset[1]]
    w.wcs.crval = [APPROX_RA, APPROX_DEC]
    w.wcs.cd = np.array([[-scale_deg * c, scale_deg * s],
                         [scale_deg * s, scale_deg * c]])
    return w


def _build_field(seed=42, n_true=40):
    rng = np.random.default_rng(seed)
    true_wcs = _synthetic_wcs(rotation_deg=12.0, crpix_offset=(-40.0, 30.0))
    true_xy = np.column_stack([
        rng.uniform(120.0, SHAPE[1] - 120.0, n_true),
        rng.uniform(120.0, SHAPE[0] - 120.0, n_true),
    ])
    sky = SkyCoord(true_wcs.pixel_to_world(true_xy[:, 0], true_xy[:, 1]))
    detected_xy = true_xy + rng.normal(0.0, 0.03, true_xy.shape)
    gmag = rng.uniform(12.0, 16.5, n_true)
    return true_wcs, detected_xy, sky, gmag


def _make_params(tmp_path: Path):
    data_dir = tmp_path / "data"
    result_dir = tmp_path / "result"
    cache_dir = result_dir / "cache"
    for d in (data_dir, result_dir, cache_dir):
        d.mkdir(parents=True, exist_ok=True)
    P = types.SimpleNamespace(
        data_dir=data_dir,
        result_dir=result_dir,
        cache_dir=cache_dir,
        pixel_scale_arcsec=SCALE,
        target_ra_deg=APPROX_RA,
        target_dec_deg=APPROX_DEC,
        fov_w_deg=0.0,
        fov_h_deg=0.0,
        gaia_radius_fudge=1.35,
        gaia_mag_max=20.0,
        gaia_wcs_mag_max=20.0,
        gaia_allow_no_cache=True,
        max_workers=1,
        parallel_max_workers=1,
        parallel_mode="thread",
        wcs_internal_min_matches=8,
        wcs_internal_rms_max_arcsec=2.0,
        # QC params (defaults)
        wcs_qc_require_wcs_ok=True,
        wcs_qc_min_match_n=20,
        wcs_qc_min_match_rate=0.20,
        wcs_qc_max_rms_px=2.5,
        wcs_qc_max_p99_px=5.0,
        wcs_qc_min_inlier_rate=0.50,
        wcs_qc_match_radius_arcsec=2.0,
        wcs_qc_clip_sigma=3.0,
        fwhm_seed_px=6.0,
    )

    class _Params:
        def __init__(self, P):
            self.P = P

        def get_file_path(self, filename):
            return Path(self.P.data_dir) / filename

    return _Params(P), data_dir, result_dir, cache_dir


def _seed_inputs(tmp_path: Path, detected_xy, sky, gmag, fnames=("f1.fits", "f2.fits")):
    params, data_dir, result_dir, cache_dir = _make_params(tmp_path)

    # Synthetic image data (flat noise — the solver uses the detection CSV, not
    # pixel peaks, when Step 4 detections are present).
    rng = np.random.default_rng(7)
    for fn in fnames:
        data = rng.normal(100.0, 5.0, SHAPE).astype(np.float32)
        hdr = fits.Header()
        hdr["DATE-OBS"] = "2023-06-15T04:00:00"
        hdr["OBJCTRA"] = "16 41 41.2"
        hdr["OBJCTDEC"] = "+36 27 40"
        fits.PrimaryHDU(data=data, header=hdr).writeto(data_dir / fn, overwrite=True)

        # Step 4 detection CSV (x,y + flux).
        s4 = step4_dir(result_dir)
        s4.mkdir(parents=True, exist_ok=True)
        df_rows = {
            "x": detected_xy[:, 0],
            "y": detected_xy[:, 1],
            "flux": np.linspace(5000.0, 1000.0, len(detected_xy)),
        }
        import pandas as pd
        pd.DataFrame(df_rows).to_csv(s4 / f"detect_{fn}.csv", index=False)

    # Pre-seed the Gaia cache so no network query is performed.
    s5 = step5_wcs_dir(result_dir)
    s5.mkdir(parents=True, exist_ok=True)
    n = len(sky)
    tab = Table({
        "source_id": np.arange(1, n + 1, dtype=np.int64),
        "ra": sky.ra.deg,
        "dec": sky.dec.deg,
        "phot_g_mean_mag": gmag,
        "phot_variable_flag": ["NOT_AVAILABLE"] * n,
        "ruwe": np.ones(n),
        "pmra": np.zeros(n),
        "pmdec": np.zeros(n),
    })
    tab.write(s5 / "gaia_fov.ecsv", format="ascii.ecsv", overwrite=True)
    (s5 / "gaia_fov_meta.json").write_text(json.dumps({
        "center_ra_deg": APPROX_RA,
        "center_dec_deg": APPROX_DEC,
        "mag_max": 99.0,
        "radius_deg": 1.0,
    }), encoding="utf-8")

    # Step 1 selection.json + targets_simbad.tsv (fallback target source).
    s1 = step1_dir(result_dir)
    s1.mkdir(parents=True, exist_ok=True)
    (s1 / "selection.json").write_text(json.dumps({"filenames": list(fnames)}), encoding="utf-8")
    return params, data_dir, result_dir, cache_dir, list(fnames)


def _wcs_keys_from_fits(path: Path) -> dict | None:
    with fits.open(path) as hdul:
        hdr = hdul[0].header
        if not hdr.get("WCS_OK", False):
            return None
        w = WCS(hdr, relax=True)
        if not w.has_celestial:
            return None
        cd = w.wcs.cd if (hasattr(w.wcs, "cd") and w.wcs.cd is not None) else None
        return {
            "crval": np.asarray(w.wcs.crval, float),
            "crpix": np.asarray(w.wcs.crpix, float),
            "cd": np.asarray(cd, float) if cd is not None else None,
            "wcssrc": str(hdr.get("WCSSRC", "")),
        }


def _read_summary(result_dir: Path):
    import pandas as pd
    p = step5_wcs_dir(result_dir) / "wcs_solve_summary.csv"
    if not p.exists():
        return None
    return pd.read_csv(p)


# ── tests ────────────────────────────────────────────────────────────────────

def test_headless_internal_solve_smoke(tmp_path):
    """run_wcs_solve(internal) solves the synthetic field and writes Step 5 outputs."""
    _, detected_xy, sky, gmag = _build_field()
    params, data_dir, result_dir, cache_dir, fnames = _seed_inputs(tmp_path / "h", detected_xy, sky, gmag)

    from apex.analysis.wcs_solve import run_wcs_solve

    summary = run_wcs_solve(
        fnames, params, data_dir, result_dir, cache_dir,
        engine="internal", use_cropped=False,
    )
    assert int(summary.get("ok", 0)) == len(fnames), summary

    # Summary CSV written with one row per frame.
    df = _read_summary(result_dir)
    assert df is not None and len(df) == len(fnames)

    # Each FITS got a valid celestial WCS tagged APEX_INTERNAL near the truth.
    for fn in fnames:
        keys = _wcs_keys_from_fits(data_dir / fn)
        assert keys is not None, f"{fn} has no WCS"
        assert keys["wcssrc"] == "APEX_INTERNAL"
        assert abs(keys["crval"][0] - APPROX_RA) < 0.05
        assert abs(keys["crval"][1] - APPROX_DEC) < 0.05


def test_gui_vs_headless_internal_parity(tmp_path):
    """GUI InternalWcsWorker (internal engine) ≡ headless run_wcs_solve."""
    pytest.importorskip("PyQt5")
    from PyQt5.QtCore import QCoreApplication
    import sys as _sys

    app = QCoreApplication.instance() or QCoreApplication(_sys.argv[:1])

    _, detected_xy, sky, gmag = _build_field()

    # ── BASELINE: current GUI worker ──
    pb, data_b, result_b, cache_b, fnames = _seed_inputs(tmp_path / "base", detected_xy, sky, gmag)
    from apex.gui.workflow.step5_wcs_plate_solving import InternalWcsWorker
    wkr = InternalWcsWorker(
        fnames, pb, pb.P.data_dir, pb.P.result_dir,
        use_cropped=False, min_matches=8, rms_max_arcsec=2.0,
    )
    done = {}
    wkr.frame_done.connect(lambda fn, info: done.__setitem__(fn, info))
    fin = {}
    wkr.finished.connect(lambda s: fin.update(s))
    wkr.run()  # synchronous (no thread start) — runs the relocated body inline
    # The GUI window normally writes the summary CSV; replicate that here so the
    # baseline has the same artifact set as the headless path.
    from apex.analysis.wcs_solve import _write_internal_summary_csv
    _write_internal_summary_csv(result_b, {k: v for k, v in done.items() if v.get("ok")})

    base_ok = int(fin.get("ok", 0))
    assert base_ok == len(fnames), fin

    # ── HEADLESS: orchestrator ──
    ph, data_h, result_h, cache_h, _ = _seed_inputs(tmp_path / "head", detected_xy, sky, gmag)
    from apex.analysis.wcs_solve import run_wcs_solve
    summary = run_wcs_solve(
        fnames, ph, data_h, result_h, cache_h,
        engine="internal", use_cropped=False,
    )
    assert int(summary.get("ok", 0)) == len(fnames), summary

    # ── DIFF FITS WCS (CRVAL/CRPIX/CD within rtol=1e-9) ──
    for fn in fnames:
        kb = _wcs_keys_from_fits(data_b / fn)
        kh = _wcs_keys_from_fits(data_h / fn)
        assert kb is not None and kh is not None
        np.testing.assert_allclose(kb["crval"], kh["crval"], rtol=1e-9, atol=1e-12)
        np.testing.assert_allclose(kb["crpix"], kh["crpix"], rtol=1e-9, atol=1e-12)
        if kb["cd"] is not None and kh["cd"] is not None:
            np.testing.assert_allclose(kb["cd"], kh["cd"], rtol=1e-9, atol=1e-15)
        assert kb["wcssrc"] == kh["wcssrc"] == "APEX_INTERNAL"

    # ── DIFF summary CSV (numeric cols equal modulo paths/timestamps) ──
    db = _read_summary(result_b).set_index("fname").sort_index()
    dh = _read_summary(result_h).set_index("fname").sort_index()
    assert list(db.index) == list(dh.index)
    for col in ("n_matches", "rms_arcsec", "rotation_deg",
                "scale_arcsec_per_px", "center_ra_deg", "center_dec_deg",
                "wcs_qc_pass", "n_match"):
        if col in db.columns and col in dh.columns:
            vb = db[col].to_numpy()
            vh = dh[col].to_numpy()
            if vb.dtype.kind in "fi" and vh.dtype.kind in "fi":
                np.testing.assert_allclose(vb, vh, rtol=1e-9, atol=1e-12,
                                           err_msg=f"summary col {col} differs")
            else:
                assert list(map(str, vb)) == list(map(str, vh)), f"summary col {col} differs"
