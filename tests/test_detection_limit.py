"""Tests for apex.analysis.detection_limit — predicted per-frame m50 + depth QC.

Two layers:

1. Pure unit tests (always run): formula algebra, PSF peak fraction,
   detection-fraction rolloff, and the frame_depth_qc flag logic.

2. Calibration-set validation (skipped when the data is not present):
   re-predicts the 50%-completeness magnitude of the 7 real-frame
   artificial-star injection runs (validation/paper/data_realframe_*/) from
   each frame's background noise (sep.Background rms x gain) and injection
   PSF peak fraction alone, and checks the residual RMS against the measured
   injection m50 stays at the ~0.1 mag level documented in
   validation/paper/논문작업/COMPLETENESS_REALFRAME_INVESTIGATION.md.
   Requires the untracked run products (empirical_psf.fits) and the original
   frames on E:\\APEX_validation — both site-local.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from apex.analysis.detection_limit import (
    detection_fraction_rolloff,
    estimate_peak_fraction_from_stars,
    frame_depth_qc,
    peak_fraction_from_psf,
    predict_frame_m50,
)
from apex.utils.constants import INSTRUMENTAL_ZMAG, PEAK_SN50_DETECTION

# PTC-measured gain of the calibration camera (Moravian C3-61000); the
# calibration frames carry no trustworthy gain header (EGAIN is 14x off),
# so the validation below pins the measured value, matching the runtime
# config used for those runs.
_CALIB_GAIN_E_PER_ADU = 0.689

_DATA_ROOT = Path(
    os.environ.get(
        "APEX_REALFRAME_DATA",
        Path(__file__).resolve().parents[1] / "validation" / "paper",
    )
)
_FRAME_ROOT = Path(
    os.environ.get("APEX_REALFRAME_FITS", r"E:\APEX_validation\reprocess")
)

# The 7 calibration runs: (run subdir, frame path relative to _FRAME_ROOT).
_CALIBRATION_RUNS = [
    ("data_realframe_M67i",           "M67/sci/pp_Messier67-0008-i.fit"),
    ("data_realframe_NGC6811R",       "NGC6811/sci/pp_NGC6811-0005-R.fit"),
    # sci/ 를 본다 — reprocess 의 reorg 단계가 calibrated/ 에서 sci/ 로 프레임을
    # **이동**(복사 아님)하므로 재처리 후 calibrated/ 는 비어 있다. 나머지 6개는
    # 원래 sci/ 를 보고 있었고 이 항목만 옛 경로에 남아 있었다 (2026-07-30).
    ("data_realframe_M13V",           "M13/sci/pp_messier13-0001-V.fit"),
    ("data_realframe_M67r_mid",       "M67/sci/pp_Messier67-0003-r.fit"),
    ("data_realframe_M67g_broad",     "M67/sci/pp_Messier67-0004-g.fit"),
    ("data_realframe_NGC6811R_broad", "NGC6811/sci/pp_NGC6811-0008-R.fit"),
    ("data_realframe_M13R_sharp",     "M13/sci/pp_messier13-0004-R.fit"),
]


# ── predict_frame_m50 ────────────────────────────────────────────────────────

def test_predict_m50_formula_value():
    # m50 = ZP - 2.5 log10(SN50 * sigma_e / p_peak), total-electron scale
    sigma_e, p_peak = 30.0, 0.0107
    expected = INSTRUMENTAL_ZMAG - 2.5 * np.log10(
        PEAK_SN50_DETECTION * sigma_e / p_peak
    )
    assert predict_frame_m50(sigma_e, p_peak) == pytest.approx(expected, abs=1e-12)


def test_predict_m50_sigma_scaling():
    # Doubling the background noise costs 2.5 log10(2) mag of depth.
    m_a = predict_frame_m50(10.0, 0.02)
    m_b = predict_frame_m50(20.0, 0.02)
    assert m_a - m_b == pytest.approx(2.5 * np.log10(2.0), abs=1e-12)


def test_predict_m50_peak_fraction_scaling():
    # A sharper PSF (larger p_peak) reaches fainter by 2.5 log10(ratio).
    m_broad = predict_frame_m50(10.0, 0.01)
    m_sharp = predict_frame_m50(10.0, 0.04)
    assert m_sharp - m_broad == pytest.approx(2.5 * np.log10(4.0), abs=1e-12)


def test_predict_m50_exptime_converts_to_count_rate_scale():
    # mag_inst is a count rate: same frame, longer exptime -> numerically
    # fainter m50 on the count-rate scale by 2.5 log10(t).
    m_total = predict_frame_m50(10.0, 0.02, exptime_s=1.0)
    m_rate = predict_frame_m50(10.0, 0.02, exptime_s=60.0)
    assert m_rate - m_total == pytest.approx(2.5 * np.log10(60.0), abs=1e-12)


@pytest.mark.parametrize(
    "sigma_e,p_peak,exptime",
    [
        (np.nan, 0.02, 1.0),
        (0.0, 0.02, 1.0),
        (-5.0, 0.02, 1.0),
        (10.0, np.nan, 1.0),
        (10.0, 0.0, 1.0),
        (10.0, 1.5, 1.0),
        (10.0, 0.02, 0.0),
    ],
)
def test_predict_m50_invalid_inputs_return_nan(sigma_e, p_peak, exptime):
    assert np.isnan(predict_frame_m50(sigma_e, p_peak, exptime_s=exptime))


# ── peak fraction estimators ─────────────────────────────────────────────────

def test_peak_fraction_from_psf_delta_and_uniform():
    delta = np.zeros((5, 5))
    delta[2, 2] = 3.7  # normalization must not matter
    assert peak_fraction_from_psf(delta) == pytest.approx(1.0)
    assert peak_fraction_from_psf(np.ones((3, 3))) == pytest.approx(1.0 / 9.0)


def test_peak_fraction_from_psf_gaussian_matches_analytic():
    # Wide Gaussian: peak fraction -> 1 / (2 pi sigma^2) as pixelization
    # becomes fine relative to the profile.
    sig = 3.0
    y, x = np.mgrid[-15:16, -15:16].astype(float)
    psf = np.exp(-(x**2 + y**2) / (2 * sig**2))
    assert peak_fraction_from_psf(psf) == pytest.approx(
        1.0 / (2 * np.pi * sig**2), rel=0.02
    )


def test_peak_fraction_from_psf_invalid():
    assert np.isnan(peak_fraction_from_psf(np.array([])))
    assert np.isnan(peak_fraction_from_psf(np.full((3, 3), np.nan)))
    assert np.isnan(peak_fraction_from_psf(np.zeros((3, 3))))


def test_estimate_peak_fraction_from_stars_median_and_selection():
    peak_e = np.array([10.0, 12.0, 14.0, 500.0, np.nan, -3.0])
    flux_e = np.array([1000.0, 1000.0, 1000.0, 100.0, 1000.0, 1000.0])
    # star 3 has ratio 5.0 (non-physical, dropped), 4 is NaN, 5 negative.
    p, n = estimate_peak_fraction_from_stars(peak_e, flux_e)
    assert n == 3
    assert p == pytest.approx(0.012)
    # mask excludes the first star -> median of {0.012, 0.014}
    mask = np.array([False, True, True, True, True, True])
    p2, n2 = estimate_peak_fraction_from_stars(peak_e, flux_e, mask)
    assert n2 == 2
    assert p2 == pytest.approx(0.013)


def test_estimate_peak_fraction_from_stars_empty():
    p, n = estimate_peak_fraction_from_stars(np.array([]), np.array([]))
    assert np.isnan(p) and n == 0


# ── detection-fraction rolloff ───────────────────────────────────────────────

def test_detection_rolloff_recovers_step_transition():
    rng = np.random.default_rng(42)
    mag = rng.uniform(10.0, 20.0, size=8000)
    # Linear ramp from 1 to 0 across mag 14.5-15.5 -> 50% at 15.0.
    p_det = np.clip((15.5 - mag) / 1.0, 0.0, 1.0)
    detected = rng.uniform(size=mag.size) < p_det
    m50 = detection_fraction_rolloff(mag, detected)
    assert m50 == pytest.approx(15.0, abs=0.15)


def test_detection_rolloff_ignores_bright_end_dip():
    # Saturated/shape-rejected bright stars can push the detected fraction
    # below 50% at the bright end; the rolloff must not lock onto that dip.
    rng = np.random.default_rng(3)
    mag = rng.uniform(9.0, 20.0, size=12000)
    p_det = np.clip((15.5 - mag) / 1.0, 0.0, 1.0)
    p_det[mag < 10.0] = 0.2  # bright-end dropout dip
    detected = rng.uniform(size=mag.size) < p_det
    m50 = detection_fraction_rolloff(mag, detected)
    assert m50 == pytest.approx(15.0, abs=0.15)


def test_detection_rolloff_no_crossing_returns_nan():
    mag = np.linspace(10, 12, 500)
    detected = np.ones(500, dtype=bool)  # never falls below 50%
    assert np.isnan(detection_fraction_rolloff(mag, detected))


def test_detection_rolloff_ignores_nonfinite_mags():
    mag = np.concatenate([np.linspace(10, 12, 500), [np.nan] * 50])
    detected = np.ones(550, dtype=bool)
    assert np.isnan(detection_fraction_rolloff(mag, detected))


# ── frame_depth_qc ───────────────────────────────────────────────────────────

def _synthetic_frame_qc(observed_shift: float = 0.0, tolerance: float = 0.5):
    """Build a self-consistent synthetic frame: stars detected exactly when
    brighter than the m50 the formula predicts (+ optional shift)."""
    rng = np.random.default_rng(7)
    sigma_adu, gain, p_peak = 40.0, 0.7, 0.015
    predicted = predict_frame_m50(sigma_adu * gain, p_peak)
    n = 6000
    mag = rng.uniform(predicted - 5.0, predicted + 2.0, size=n)
    m50_true = predicted + observed_shift
    p_det = np.clip((m50_true + 0.5 - mag) / 1.0, 0.0, 1.0)
    detected = rng.uniform(size=n) < p_det
    flux_e = 10 ** (0.4 * (INSTRUMENTAL_ZMAG - mag))
    peak_e = flux_e * p_peak
    bright = detected & (mag < predicted - 2.0)
    return frame_depth_qc(
        sky_sigma_adu=sigma_adu,
        gain_e_per_adu=gain,
        exptime_s=1.0,
        peak_e=peak_e,
        flux_e=flux_e,
        peak_star_mask=bright,
        mag_inst=mag,
        detected=detected,
        tolerance_mag=tolerance,
    )


def test_frame_depth_qc_consistent_frame_flags_ok():
    qc = _synthetic_frame_qc(observed_shift=0.0)
    assert qc["depth_qc_flag"] == "ok"
    assert qc["p_peak_frame"] == pytest.approx(0.015, rel=1e-6)
    assert abs(qc["depth_delta_mag"]) < 0.3
    assert qc["n_peak_stars"] > 0


def test_frame_depth_qc_shallow_frame_is_flagged():
    # Detection dies 1.2 mag brighter than sky+seeing predict -> anomaly.
    qc = _synthetic_frame_qc(observed_shift=-1.2)
    assert qc["depth_qc_flag"] == "depth_shallow"
    assert qc["depth_delta_mag"] < -0.5


def test_frame_depth_qc_uncomputable_stays_unflagged():
    qc = frame_depth_qc(
        sky_sigma_adu=np.nan,
        gain_e_per_adu=0.7,
        exptime_s=1.0,
        peak_e=np.array([]),
        flux_e=np.array([]),
        peak_star_mask=np.array([], dtype=bool),
        mag_inst=np.array([]),
        detected=np.array([], dtype=bool),
    )
    assert qc["depth_qc_flag"] == ""
    assert np.isnan(qc["predicted_m50"])
    assert np.isnan(qc["depth_delta_mag"])


# ── calibration-set validation (site-local data; skipped when absent) ────────

def _available_runs():
    runs = []
    for sub, rel in _CALIBRATION_RUNS:
        rd = _DATA_ROOT / sub / "artificial_star" / "benchmark_run"
        fp = _FRAME_ROOT / rel
        if (rd / "stars.csv").exists() and (rd / "empirical_psf.fits").exists() \
                and fp.exists():
            runs.append((sub, rd, fp))
    return runs


@pytest.mark.slow
def test_predicted_m50_reproduces_injection_calibration_set():
    sep = pytest.importorskip("sep")
    fits = pytest.importorskip("astropy.io.fits")
    pd = pytest.importorskip("pandas")

    runs = _available_runs()
    if len(runs) < len(_CALIBRATION_RUNS):
        pytest.skip(
            f"calibration data incomplete ({len(runs)}/{len(_CALIBRATION_RUNS)} "
            "runs with stars.csv + empirical_psf.fits + frame FITS)"
        )

    residuals = {}
    for sub, rd, fp in runs:
        stars = pd.read_csv(rd / "stars.csv")
        stars = stars[~stars["baseline_confounded"].astype(bool)]
        m50_measured = detection_fraction_rolloff(
            stars["magnitude_true"].to_numpy(float),
            stars["recovered"].to_numpy(bool),
        )
        assert np.isfinite(m50_measured), f"{sub}: no completeness rolloff"

        psf = fits.getdata(rd / "empirical_psf.fits").astype(float)
        p_peak = peak_fraction_from_psf(psf)
        frame = fits.getdata(str(fp)).astype(np.float64)
        sigma_e = float(np.median(sep.Background(frame).rms())) * _CALIB_GAIN_E_PER_ADU

        m50_predicted = predict_frame_m50(sigma_e, p_peak)
        residuals[sub] = m50_predicted - m50_measured

    res = np.array(list(residuals.values()))
    rms = float(np.sqrt(np.mean(res**2)))
    worst = max(residuals, key=lambda k: abs(residuals[k]))
    # Investigation record: residual RMS ~0.05 mag over the 7 runs; the gate
    # here allows drift up to the documented ~0.1 mag level.
    assert rms <= 0.12, f"residual RMS {rms:.3f} mag (worst {worst}: {residuals[worst]:+.3f})"
    assert max(abs(r) for r in res) <= 0.25, (
        f"worst-frame residual {residuals[worst]:+.3f} mag ({worst})"
    )
