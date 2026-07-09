"""Detector-calibration core tests (Step 0).

Synthesises frames with KNOWN bias/dark/flat, builds masters, calibrates, and
checks that the injected detector signature is recovered to the true science
frame.  This is the honesty gate for the paper's raw->science claim: it must be
green before any manuscript wording is changed back to "raw".
"""

from __future__ import annotations

import numpy as np
import pytest
from astropy.io import fits

from apex.analysis import calibration as cal
from apex.analysis.calibration import CalibrationOptions


# --------------------------------------------------------------------------
# Synthetic detector model:  raw = bias + dark_rate*exptime + flat * science
# --------------------------------------------------------------------------

SIZE = 128
DARK_EXP = 10.0
LIGHT_EXP = 10.0
FLAT_EXP = 10.0
BIAS_BASE = 300.0
DARK_RATE = 0.5           # DN / s
SKY = 100.0


def _true_bias() -> np.ndarray:
    # constant pedestal + a gentle column gradient (fixed pattern)
    col = np.linspace(0.0, 20.0, SIZE)[np.newaxis, :]
    return (BIAS_BASE + col) * np.ones((SIZE, SIZE), dtype=np.float64)


def _true_dark_rate() -> np.ndarray:
    rate = DARK_RATE * np.ones((SIZE, SIZE), dtype=np.float64)
    # a couple of hot pixels
    rate[10, 10] = 50.0
    rate[100, 64] = 80.0
    return rate


def _true_flat() -> np.ndarray:
    # radial vignette, normalised to unit median
    yy, xx = np.mgrid[0:SIZE, 0:SIZE]
    cy = cx = (SIZE - 1) / 2.0
    r2 = (xx - cx) ** 2 + (yy - cy) ** 2
    flat = 1.0 - 0.2 * (r2 / r2.max())     # ~0.8 at corners, 1.0 centre
    return flat / np.median(flat)


def _true_science() -> np.ndarray:
    img = SKY * np.ones((SIZE, SIZE), dtype=np.float64)
    # a few Gaussian "stars"
    for (y, x, flux) in [(40, 50, 5000.0), (80, 30, 2000.0), (60, 90, 8000.0)]:
        yy, xx = np.mgrid[0:SIZE, 0:SIZE]
        img += flux * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * 2.0 ** 2))
    return img


def _write(path, data, exptime, rng, read_noise=3.0):
    noisy = data + rng.normal(0.0, read_noise, size=data.shape)
    hdu = fits.PrimaryHDU(noisy.astype(np.float32))
    hdu.header["EXPTIME"] = float(exptime)
    hdu.header["FILTER"] = "V"
    hdu.writeto(path, overwrite=True)
    return str(path)


def _make_dataset(tmp_path, n=15, seed=1234):
    rng = np.random.default_rng(seed)
    B = _true_bias()
    Dr = _true_dark_rate()
    F = _true_flat()
    S = _true_science()

    bias_paths, dark_paths, flat_paths = [], [], []
    for i in range(n):
        bias_paths.append(_write(tmp_path / f"bias_{i}.fits", B, 0.0, rng))
        dark_paths.append(_write(tmp_path / f"dark_{i}.fits", B + Dr * DARK_EXP, DARK_EXP, rng))
        # flat: illuminated field through the flat response
        flat_signal = B + Dr * FLAT_EXP + F * 20000.0
        flat_paths.append(_write(tmp_path / f"flat_{i}.fits", flat_signal, FLAT_EXP, rng))

    raw_light = B + Dr * LIGHT_EXP + F * S
    light_path = _write(tmp_path / "light_0.fits", raw_light, LIGHT_EXP, rng)
    return bias_paths, dark_paths, flat_paths, light_path, (B, Dr, F, S)


# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------

def test_combine_median_rejects_outlier():
    base = np.full((16, 16), 100.0)
    frames = [base.copy() for _ in range(9)]
    spiked = base.copy()
    spiked[8, 8] = 60000.0            # cosmic-ray-like spike in one frame
    frames.append(spiked)
    med = cal.combine_frames(frames, method="median")
    assert med[8, 8] == pytest.approx(100.0, abs=1e-6)

    sc = cal.combine_frames(frames, method="sigmaclip_mean", sigma_low=3, sigma_high=3)
    assert sc[8, 8] == pytest.approx(100.0, abs=1.0)


def test_bias_dark_flat_roundtrip(tmp_path):
    opts = CalibrationOptions(combine_method="median", pedestal_mode="none")
    bias_p, dark_p, flat_p, light_p, (B, Dr, F, S) = _make_dataset(tmp_path)

    mbias, _ = cal.build_master_bias(bias_p, opts)
    mdark, dexp, _ = cal.build_master_dark(dark_p, opts, master_bias=mbias)
    mflat, fprov = cal.build_master_flat(flat_p, opts, master_bias=mbias,
                                         master_dark=mdark, dark_exp=dexp)

    # master sanity
    assert dexp == pytest.approx(DARK_EXP, abs=0.01)
    assert np.median(mflat) == pytest.approx(1.0, abs=1e-3)

    calibrated, header, qc = cal.calibrate_light_file(
        light_p, opts, master_bias=mbias, master_dark=mdark,
        dark_exp=dexp, master_flat=mflat,
    )

    # Recovered frame should match the true science frame (pedestal=none).
    # The residual scatter floor is the single light frame's own read noise
    # (~3 DN); what calibration must guarantee is zero *systematic* offset and
    # removal of the vignette structure.
    resid = calibrated.astype(np.float64) - S
    med_resid = float(np.nanmedian(resid))
    robust = float(np.nanmedian(np.abs(resid - med_resid)) * 1.4826)
    assert abs(med_resid) < 1.0        # DN — no systematic offset
    assert robust < 5.0                # DN — at the single-frame noise floor

    # Flat correction really happened: WITHOUT the flat, a corner sky pixel is
    # off by ~(F-1)*sky (the vignette), while WITH the flat it is within noise.
    no_flat, _, _ = cal.calibrate_light_file(
        light_p, opts, master_bias=mbias, master_dark=mdark, dark_exp=dexp,
    )
    cy = cx = 5  # corner sky region (no star)
    assert abs(float(no_flat[cy, cx]) - S[cy, cx]) > 10.0   # vignette present
    assert abs(float(calibrated[cy, cx]) - S[cy, cx]) < 8.0  # vignette removed

    # star flux preserved: peak near true peak within noise
    assert calibrated[60, 90] == pytest.approx(S[60, 90], abs=30.0)


def test_pedestal_is_additive_constant_and_recorded(tmp_path):
    bias_p, dark_p, flat_p, light_p, _ = _make_dataset(tmp_path)
    base = CalibrationOptions(pedestal_mode="none")
    fixed = CalibrationOptions(pedestal_mode="fixed", pedestal_value=250.0)

    mbias, _ = cal.build_master_bias(bias_p, base)
    mdark, dexp, _ = cal.build_master_dark(dark_p, base, master_bias=mbias)
    mflat, _ = cal.build_master_flat(flat_p, base, master_bias=mbias,
                                     master_dark=mdark, dark_exp=dexp)

    cal_none, _, _ = cal.calibrate_light_file(light_p, base, master_bias=mbias,
                                              master_dark=mdark, dark_exp=dexp,
                                              master_flat=mflat)
    cal_ped, hdr, qc = cal.calibrate_light_file(light_p, fixed, master_bias=mbias,
                                                master_dark=mdark, dark_exp=dexp,
                                                master_flat=mflat)
    diff = cal_ped.astype(np.float64) - cal_none.astype(np.float64)
    assert np.nanmedian(diff) == pytest.approx(250.0, abs=1e-3)
    assert np.nanstd(diff) == pytest.approx(0.0, abs=1e-3)   # purely additive
    assert int(hdr["PEDESTAL"]) == 250
    assert qc["pedestal"] == pytest.approx(250.0)


def test_flat_dead_pixel_becomes_nan(tmp_path):
    opts = CalibrationOptions(flat_min=0.01, pedestal_mode="none")
    bias_p, dark_p, flat_p, light_p, _ = _make_dataset(tmp_path)
    mbias, _ = cal.build_master_bias(bias_p, opts)
    data, header = cal.load_frame(light_p, opts)
    flat = np.ones_like(data)
    flat[5, 5] = 0.0                 # dead pixel (QE=0)
    out, _, qc = cal.calibrate_light(data, header, opts, master_bias=mbias,
                                     master_flat=flat)
    assert np.isnan(out[5, 5])
    assert qc["flat_bad_pct"] > 0.0


def test_dark_scaled_by_exposure_ratio(tmp_path):
    # light exposure twice the dark reference -> dark subtracted at k=2
    rng = np.random.default_rng(7)
    opts = CalibrationOptions(pedestal_mode="none", dark_scale=True)
    Dr = _true_dark_rate()
    mdark = (Dr * DARK_EXP)                       # "master dark" at 10s
    data = (Dr * (2 * DARK_EXP)).astype(np.float64)   # light at 20s, no sky
    hdr = fits.Header()
    hdr["EXPTIME"] = 2 * DARK_EXP
    out, _, qc = cal.calibrate_light(data, hdr, opts, master_dark=mdark,
                                     dark_exp=DARK_EXP)
    assert qc["dark_scale"] == pytest.approx(2.0, abs=1e-6)
    assert float(np.median(out)) == pytest.approx(0.0, abs=1e-6)


def _write_master(path, data, imagetyp, exptime=0.0, filt=""):
    hdu = fits.PrimaryHDU(np.asarray(data, dtype=np.float32))
    hdu.header["IMAGETYP"] = imagetyp
    hdu.header["EXPTIME"] = float(exptime)
    if filt:
        hdu.header["FILTER"] = filt
    hdu.writeto(path, overwrite=True)
    return str(path)


def test_prebuilt_master_used_directly(tmp_path):
    """A frame flagged IMAGETYP 'MASTER <kind>' is used verbatim, not re-stacked
    or re-reduced — mirrors AIPPI's manual-master auto-detection."""
    opts = CalibrationOptions()
    rng = np.random.default_rng(0)
    mbias = (500.0 + rng.normal(0, 1, (32, 32))).astype(np.float32)

    # master BIAS: returned verbatim
    pb = _write_master(tmp_path / "master_bias.fits", mbias, "Master Bias")
    out, prov = cal.build_master_bias([pb], opts)
    assert prov["master_input"] is True and prov["n_frames"] == 1
    np.testing.assert_allclose(out, mbias, atol=1e-2)

    # master DARK: verbatim + NOT bias-subtracted even when a master_bias is given
    mdark = (40.0 + rng.normal(0, 1, (32, 32))).astype(np.float32)
    pd = _write_master(tmp_path / "master_dark.fits", mdark, "Master Dark", exptime=30.0)
    outd, dexp, provd = cal.build_master_dark([pd], opts, master_bias=mbias)
    assert provd["master_input"] and provd["bias_subtracted"]
    assert dexp == pytest.approx(30.0)
    np.testing.assert_allclose(outd, mdark, atol=1e-2)      # bias was NOT subtracted

    # master FLAT: verbatim, near-unit median left as-is (not re-reduced)
    mflat = np.full((32, 32), 0.98, np.float32)
    pf = _write_master(tmp_path / "master_flat.fits", mflat, "Master Flat", filt="V")
    outf, provf = cal.build_master_flat([pf], opts, master_bias=mbias)
    assert provf["master_input"]
    assert float(np.median(outf)) == pytest.approx(0.98, abs=1e-3)

    # a master flat with a non-unit median is defensively renormalised to 1
    pf2 = _write_master(tmp_path / "master_flat2.fits",
                        np.full((32, 32), 1000.0, np.float32), "Master Flat", filt="V")
    outf2, _ = cal.build_master_flat([pf2], opts)
    assert float(np.median(outf2)) == pytest.approx(1.0, abs=1e-3)
