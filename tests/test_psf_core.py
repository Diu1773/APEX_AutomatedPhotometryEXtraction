import numpy as np

from apex.utils.psf_core import estimate_psf_core_cut, psf_core_keep_mask


def test_estimate_psf_core_cut_auto_density_excludes_center():
    rng = np.random.default_rng(123)
    core = rng.normal(loc=(500.0, 500.0), scale=20.0, size=(250, 2))
    field = rng.uniform(0.0, 1000.0, size=(600, 2))
    xy = np.vstack([core, field])

    cut = estimate_psf_core_cut(
        xy,
        (1000, 1000),
        5.0,
        enabled=True,
        center_mode="auto",
        radius_px=80.0,
        auto_min_density_ratio=1.2,
    )

    assert cut.enabled
    assert abs(cut.center_x - 500.0) < 40.0
    assert abs(cut.center_y - 500.0) < 40.0
    assert cut.n_excluded > 100

    keep = psf_core_keep_mask(xy, cut)
    assert keep.dtype == bool
    assert np.sum(~keep) == cut.n_excluded


def test_estimate_psf_core_cut_stays_off_for_sparse_inputs():
    xy = np.array([[10.0, 10.0], [20.0, 20.0], [30.0, 30.0]])

    cut = estimate_psf_core_cut(xy, (100, 100), 4.0, enabled=True, auto_min_sources=10)

    assert not cut.enabled
    assert cut.reason == "too_few_sources"
