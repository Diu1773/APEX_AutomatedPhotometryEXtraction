from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from apex.analysis.psf_diagnostics import (
    draw_psf_final_diagnostics,
    epsf_shape_metrics,
    load_psf_final_diagnostic_data,
)


def _circular_gaussian(size: int = 41, sigma: float = 4.0) -> np.ndarray:
    yy, xx = np.mgrid[:size, :size]
    center = 0.5 * (size - 1)
    model = np.exp(-0.5 * ((xx - center) ** 2 + (yy - center) ** 2) / sigma**2)
    return model / model.sum()


def test_load_psf_final_diagnostic_data_position_matches_catalogues(tmp_path):
    result_dir = tmp_path / "result"
    psf_dir = result_dir / "cmd_psf"
    aperture_dir = result_dir / "step7_forced_phot"
    psf_dir.mkdir(parents=True)
    aperture_dir.mkdir(parents=True)
    filename = "frame-B.fit"

    pd.DataFrame(
        {
            "x_fit": [10.1, 30.2, 50.0],
            "y_fit": [20.0, 40.1, 60.0],
            "flux_psf_e": [1000.0, 500.0, 900.0],
            "snr_psf": [80.0, 30.0, 50.0],
            "qfit": [0.1, 0.2, 0.3],
            "reduced_chi2": [1.0, 1.2, 1.1],
            "iter_found": [1, 2, 1],
            "flags_psf": [0, 0, 1],
        }
    ).to_csv(psf_dir / f"photometry_{filename}.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "x_fit": [10.0, 30.0, 70.0],
            "y_fit": [20.0, 40.0, 80.0],
            "flux_e": [800.0, 400.0, 300.0],
            "snr": [60.0, 20.0, 10.0],
            "sky": [100.0, 101.0, 99.0],
            "sky_std": [4.0, 4.2, 4.1],
            "detected_flag": [True, False, True],
            "bad_phot_flag": [False, False, False],
            "off_frame_flag": [False, False, False],
        }
    ).to_csv(aperture_dir / f"photometry_{filename}.tsv", sep="\t", index=False)

    matched = load_psf_final_diagnostic_data(result_dir, filename)

    assert len(matched) == 2
    assert matched["match_distance_px"].max() < 0.25
    assert matched["detected"].tolist() == [True, False]
    np.testing.assert_allclose(
        matched["delta_mag"],
        -2.5 * np.log10(np.array([1000.0 / 800.0, 500.0 / 400.0])),
    )
    assert np.isfinite(matched["nearest_neighbor_px"]).all()


def test_epsf_shape_metrics_are_small_for_circular_symmetric_model():
    metrics = epsf_shape_metrics(_circular_gaussian())

    assert metrics["ellipticity"] < 1e-6
    assert metrics["rotation_asymmetry"] < 1e-6
    assert metrics["centroid_offset_native_px"] < 1e-6


def test_draw_psf_final_diagnostics_builds_six_panel_summary():
    count = 120
    snr = np.geomspace(5.0, 200.0, count)
    data = pd.DataFrame(
        {
            "x": np.linspace(100.0, 900.0, count),
            "y": np.linspace(900.0, 100.0, count),
            "snr_aperture": snr,
            "delta_mag": -0.08 - 0.25 / np.sqrt(snr),
            "nearest_neighbor_px": np.linspace(4.0, 50.0, count),
        }
    )
    figure = Figure(figsize=(12, 7))
    epsf_reference = pd.DataFrame(
        {
            "x": [200.0, 400.0, 600.0],
            "y": [300.0, 500.0, 700.0],
            "selected": [True, False, True],
        }
    )

    summary = draw_psf_final_diagnostics(
        figure,
        data,
        _circular_gaussian(),
        filename="frame-B.fit",
        fwhm_px=5.0,
        pixel_scale_arcsec=0.4,
        core_center=(500.0, 500.0),
        core_radius_px=100.0,
        epsf_reference=epsf_reference,
    )

    assert summary["n_matched"] == count
    assert summary["high_snr_reference_n"] >= 20
    assert np.isfinite(summary["high_snr_reference_offset_mag"])
    assert np.isfinite(summary["epsf_rotation_asymmetry"])
    assert summary["epsf_reference_n"] == 2
    assert len(figure.axes) == 8  # Six panels plus two colorbars.


def test_final_diagnostics_warns_for_non_gaussian_high_snr_tail():
    count = 100
    delta = np.zeros(count)
    delta[::10] = 1.0
    data = pd.DataFrame(
        {
            "x": np.linspace(100.0, 900.0, count),
            "y": np.linspace(900.0, 100.0, count),
            "snr_aperture": np.full(count, 100.0),
            "delta_mag": delta,
            "nearest_neighbor_px": np.full(count, 30.0),
        }
    )

    summary = draw_psf_final_diagnostics(
        Figure(figsize=(12, 7)),
        data,
        _circular_gaussian(),
        filename="outlier-frame.fit",
        fwhm_px=5.0,
    )

    assert summary["status"] == "CHECK"
    assert summary["high_snr_robust_scatter_mag"] == 0.0
    assert summary["high_snr_rmse_mag"] > 0.2
    assert summary["high_snr_outlier_fraction_0p2mag"] == 0.1
    assert any("non-Gaussian tail" in warning for warning in summary["warnings"])
