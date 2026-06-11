from __future__ import annotations

import numpy as np
import pandas as pd

from apex.benchmark.artificial_stars import (
    extract_empirical_psf,
    inject_catalog,
    sample_injection_catalog,
)
from apex.benchmark.metrics import (
    mark_baseline_confounding,
    match_truth_to_detections,
)


def _add_gaussian(image, x, y, flux, sigma=1.5):
    yy, xx = np.mgrid[0:image.shape[0], 0:image.shape[1]]
    model = np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / sigma**2)
    image += flux * model / model.sum()


def test_empirical_psf_is_normalized_from_real_star_cutouts():
    image = np.full((128, 128), 100.0)
    positions = [(25.2, 25.4), (65.1, 25.2), (104.7, 25.3), (25.3, 85.1), (85.2, 90.4)]
    for x, y in positions:
        _add_gaussian(image, x, y, 20000.0)
    detections = pd.DataFrame(positions, columns=["x", "y"])

    result = extract_empirical_psf(
        image,
        detections,
        3.5,
        saturation_adu=65000.0,
        radius_fwhm=2.5,
        isolation_fwhm=3.0,
        min_stars=3,
    )

    assert result.n_used >= 3
    assert result.data.shape[0] % 2 == 1
    assert np.isclose(result.data.sum(), 1.0)
    center = result.data.shape[0] // 2
    assert np.unravel_index(np.argmax(result.data), result.data.shape) == (center, center)


def test_injection_is_deterministic_and_records_source_poisson_flux():
    image = np.full((96, 96), 250.0)
    psf = np.zeros((9, 9), dtype=float)
    psf[4, 4] = 1.0
    catalog = pd.DataFrame(
        {
            "injection_id": [1, 2],
            "x_true": [30.25, 65.75],
            "y_true": [40.5, 55.2],
            "magnitude_true": [18.0, 19.0],
        }
    )

    args = dict(
        gain_e_per_adu=2.0,
        zeropoint_mag=25.0,
    )
    first = inject_catalog(image, psf, catalog, rng=np.random.default_rng(123), **args)
    second = inject_catalog(image, psf, catalog, rng=np.random.default_rng(123), **args)

    assert np.array_equal(first[0], second[0])
    assert np.array_equal(first[2], second[2])
    assert np.allclose(first[1].sum() * 2.0, first[3]["flux_expected_e"].sum())
    assert np.allclose(first[2].sum() * 2.0, first[3]["flux_realized_e"].sum())
    expected = np.sum(10.0 ** (-0.4 * (catalog["magnitude_true"] - 25.0)))
    assert np.isclose(first[3]["flux_expected_e"].sum(), expected)


def test_sampling_respects_injection_separation():
    rng = np.random.default_rng(7)
    catalog = sample_injection_catalog(
        (256, 256),
        np.array([[128.0, 128.0]]),
        count=12,
        magnitude_min=16.0,
        magnitude_max=20.0,
        fwhm_px=3.0,
        psf_size=21,
        rng=rng,
        isolated_fraction=0.5,
        min_injected_sep_fwhm=4.0,
    )
    xy = catalog[["x_true", "y_true"]].to_numpy(float)
    distances = np.hypot(
        xy[:, None, 0] - xy[None, :, 0],
        xy[:, None, 1] - xy[None, :, 1],
    )
    distances += np.eye(len(xy)) * 1e9
    assert distances.min() >= 12.0
    assert set(catalog["placement_stratum"]) == {"isolated", "field"}


def test_matching_is_one_to_one_and_marks_existing_sources():
    truth = pd.DataFrame(
        {
            "injection_id": [1, 2],
            "x_true": [10.0, 10.8],
            "y_true": [10.0, 10.0],
        }
    )
    detections = pd.DataFrame({"x": [10.1], "y": [10.0]})
    matched = match_truth_to_detections(truth, detections, radius_px=1.0)
    assert matched["recovered"].sum() == 1

    baseline = pd.DataFrame({"x": [10.0], "y": [10.0]})
    marked = mark_baseline_confounding(matched, baseline, radius_px=0.3)
    assert bool(marked.loc[0, "baseline_confounded"])
    assert not bool(marked.loc[1, "baseline_confounded"])
