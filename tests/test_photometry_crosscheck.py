"""Fast, Qt-free tests for the photometric cross-check harness.

Build a synthetic frame with known Gaussian stars, measure the "truth" with
photutils (so the APEX side is well-defined on the SAME apertures), write a
synthetic APEX-style Step 7 TSV + index, then run ``crosscheck_sep`` and assert
that sep agrees with photutils to << 0.05 mag on identical apertures. Also
assert ``crosscheck_iraf`` skips cleanly when IRAF is absent.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from apex.benchmark import photometry_crosscheck as pc


# ── synthetic frame fixture ─────────────────────────────────────────────────────

FWHM_PX = 4.0
GAIN = 1.5
SKY = 200.0
ZP = 25.0


def _gaussian_sigma(fwhm: float) -> float:
    return fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))


def _make_synthetic_frame(seed: int = 7):
    """A 512x512 frame with bright, isolated Gaussian stars on a flat sky."""
    rng = np.random.default_rng(seed)
    h = w = 512
    sigma = _gaussian_sigma(FWHM_PX)

    n_stars = 60
    margin = 40
    xs = rng.uniform(margin, w - margin, n_stars)
    ys = rng.uniform(margin, h - margin, n_stars)
    # Keep stars isolated (annulus radius ~ 4*FWHM = 16 px); enforce min spacing.
    keep_x, keep_y = [], []
    for x, y in zip(xs, ys):
        if all(math.hypot(x - kx, y - ky) > 60 for kx, ky in zip(keep_x, keep_y)):
            keep_x.append(x)
            keep_y.append(y)
    xs = np.asarray(keep_x)
    ys = np.asarray(keep_y)
    n_stars = xs.size

    # Bright, well-detected fluxes (total ADU).
    fluxes = rng.uniform(2.0e4, 5.0e5, n_stars)

    yy, xx = np.mgrid[0:h, 0:w]
    img = np.full((h, w), SKY, dtype=np.float64)
    norm = 1.0 / (2.0 * math.pi * sigma * sigma)
    for x, y, f in zip(xs, ys, fluxes):
        img += f * norm * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2.0 * sigma * sigma))
    # Poisson + small read noise so SNR is high but finite.
    img = rng.poisson(np.clip(img, 0, None)).astype(np.float64)
    img += rng.normal(0.0, 3.0, img.shape)
    return img, xs, ys


def _measure_truth_photutils(img, xs, ys):
    """Measure the truth with photutils on APEX's aperture geometry."""
    from photutils.aperture import (
        CircularAperture,
        CircularAnnulus,
        ApertureStats,
        aperture_photometry,
    )
    from astropy.stats import SigmaClip

    scales = {
        "aperture_scale": pc.DEFAULT_APERTURE_SCALE,
        "ref_ap_scale": pc.DEFAULT_REF_AP_SCALE,
        "ann_in_scale": pc.DEFAULT_ANN_IN_SCALE,
        "dann_scale": pc.DEFAULT_DANN_SCALE,
        "min_r_ap_px": pc.DEFAULT_MIN_R_AP_PX,
        "ann_min_gap_px": pc.DEFAULT_ANN_MIN_GAP_PX,
        "gain": GAIN,
    }
    r_ap, r_in, r_out = pc._radii(FWHM_PX, scales)
    positions = np.column_stack([xs, ys])
    aperture = CircularAperture(positions, r=r_ap)
    annulus = CircularAnnulus(positions, r_in=r_in, r_out=r_out)
    bkg = ApertureStats(img, annulus, sigma_clip=SigmaClip(sigma=3.0, maxiters=5))
    phot = aperture_photometry(img, aperture)
    net = np.asarray(phot["aperture_sum"], dtype=float) - np.asarray(bkg.median, dtype=float) * aperture.area
    return net, r_ap


@pytest.fixture(scope="module")
def synthetic_result_dir(tmp_path_factory):
    """A minimal APEX result dir: raw FITS + step7 TSV + photometry_index.csv."""
    from astropy.io import fits

    root = tmp_path_factory.mktemp("apex_xcheck")
    data_dir = root  # raw frame lives beside the result dir (APEX default layout)
    result_dir = root / "result"
    step7 = result_dir / "step7_forced_phot"
    step7.mkdir(parents=True)

    img, xs, ys = _make_synthetic_frame()
    frame_name = "synthetic-0001-r.fits"
    fits.PrimaryHDU(data=img.astype(np.float32)).writeto(data_dir / frame_name, overwrite=True)

    net_flux, _r_ap = _measure_truth_photutils(img, xs, ys)
    # APEX-style instrumental magnitude + a plausible high SNR.
    flux_e = net_flux * GAIN
    mag_inst = ZP - 2.5 * np.log10(np.clip(flux_e, 1e-9, None))
    snr = np.sqrt(np.clip(flux_e, 0, None))  # shot-noise-dominated proxy, all >> 20

    tsv = pd.DataFrame(
        {
            "x": xs,
            "y": ys,
            "flux_net_adu": net_flux,
            "mag_inst": mag_inst,
            "mag_err": 1.0857 / np.clip(snr, 1e-9, None),
            "snr": snr,
            "sky": SKY,
            "apcorr": 1.0,
            "apcorr_applied": True,
            "is_saturated": False,
            "is_nonlinear": False,
            "bad_phot_flag": False,
            "off_frame_flag": False,
        }
    )
    tsv_path = step7 / f"photometry_{frame_name}.tsv"
    tsv.to_csv(tsv_path, sep="\t", index=False)

    pd.DataFrame(
        [{"file": frame_name, "filter": "r", "status": "ok", "fwhm_px": FWHM_PX, "path": str(tsv_path)}]
    ).to_csv(step7 / "photometry_index.csv", index=False)

    return {"result_dir": result_dir, "frame": frame_name, "tsv": tsv_path, "data_dir": data_dir}


# ── tests ────────────────────────────────────────────────────────────────────

def test_resolve_frame_image_finds_raw_frame(synthetic_result_dir):
    path = pc.resolve_frame_image(
        synthetic_result_dir["result_dir"], synthetic_result_dir["frame"]
    )
    assert path.exists()
    assert path.name == synthetic_result_dir["frame"]


def test_crosscheck_sep_agrees_with_photutils(synthetic_result_dir):
    res = pc.crosscheck_sep(
        pc.resolve_frame_image(synthetic_result_dir["result_dir"], synthetic_result_dir["frame"]),
        synthetic_result_dir["tsv"],
        fwhm_px=FWHM_PX,
        gain=GAIN,
        snr_min=20.0,
    )
    assert res["status"] == "ok"
    assert res["n_matched"] >= 20

    # sep vs photutils on identical apertures must agree to << 0.05 mag.
    assert abs(res["dmag_median"]) < 0.02, res
    assert res["dmag_mad"] < 0.01, res
    assert res["dmag_rms"] < 0.02, res
    assert res["robust_sigma"] < 0.01, res
    assert res["pearson_r"] > 0.999, res
    assert res["frac_within_0p05"] > 0.95, res

    # Positions land on real sources.
    assert res["sanity"]["positions_plausible"] is True
    assert res["sanity"]["star_over_sky_sigma"] > 5.0


def test_crosscheck_sep_writes_plots(synthetic_result_dir, tmp_path):
    out = tmp_path / "plots"
    res = pc.crosscheck_sep(
        pc.resolve_frame_image(synthetic_result_dir["result_dir"], synthetic_result_dir["frame"]),
        synthetic_result_dir["tsv"],
        fwhm_px=FWHM_PX,
        gain=GAIN,
        output_dir=out,
    )
    assert res["status"] == "ok"
    assert "plots" in res and len(res["plots"]) == 2
    for name in res["plots"]:
        assert (out / name).exists()


def test_run_crosscheck_orchestration(synthetic_result_dir, tmp_path):
    res = pc.run_crosscheck(
        synthetic_result_dir["result_dir"],
        tmp_path / "report",
        references=("sep",),
    )
    assert res["status"] == "ok"
    assert res["frame"] == synthetic_result_dir["frame"]
    assert res["positions_plausible"] is True
    assert res["references"]["sep"]["status"] == "ok"
    assert (tmp_path / "report" / "crosscheck.json").exists()


def test_crosscheck_iraf_skips_cleanly_when_unavailable(synthetic_result_dir, monkeypatch):
    # Force "no IRAF runtime" regardless of the host environment.
    monkeypatch.setattr(pc, "_iraf_runtime_available", lambda runtime_cmd: False)
    res = pc.crosscheck_iraf(
        pc.resolve_frame_image(synthetic_result_dir["result_dir"], synthetic_result_dir["frame"]),
        synthetic_result_dir["tsv"],
        fwhm_px=FWHM_PX,
        gain=GAIN,
    )
    assert res["status"] == "skipped"
    assert "reason" in res and res["reason"]


def test_load_apex_tsv_filters_bad_rows(synthetic_result_dir):
    df = pc.load_apex_tsv(synthetic_result_dir["tsv"], snr_min=20.0)
    assert len(df) >= 20
    assert (df["flux_net_adu"] > 0).all()
    assert (df["snr"] >= 20.0).all()
    assert "row_id" in df.columns


def test_import_is_qt_free():
    import sys

    # Importing the module must not drag in PyQt5 or apex.gui.
    before = {m for m in sys.modules if m.startswith("PyQt5") or m.startswith("apex.gui")}
    import importlib

    importlib.reload(pc)
    after = {m for m in sys.modules if m.startswith("PyQt5") or m.startswith("apex.gui")}
    assert not (after - before), f"crosscheck import pulled in GUI modules: {after - before}"
