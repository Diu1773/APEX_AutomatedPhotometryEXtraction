import numpy as np

import apex.analysis.psf_policy as psf_policy

from apex.analysis.psf_policy import (
    estimate_psf_flux_seeds,
    local_group_policy,
    merge_forced_catalog_seeds,
    measure_epsf_annulus_contamination,
    nearest_other_source_distance,
    plan_epsf_stars,
    plan_psf_fit_window,
    select_epsf_reference_stars,
    select_spatially_balanced,
)


def test_psf_fit_window_auto_reaches_energy_with_two_fwhm_floor():
    yy, xx = np.mgrid[-17:18, -17:18]
    psf = (1.0 + (xx**2 + yy**2) / 5.0**2) ** -2.0

    plan = plan_psf_fit_window(
        psf,
        8.0,
        mode="auto",
        target_energy_fraction=0.90,
    )

    assert plan.mode == "auto"
    assert plan.shape_px % 2 == 1
    assert plan.shape_px >= 17
    assert plan.energy_fraction >= 0.90
    assert plan.noise_equivalent_area_px > 0


def test_psf_fit_window_manual_uses_multiplier():
    psf = np.ones((35, 35), dtype=float)

    plan = plan_psf_fit_window(
        psf,
        8.67,
        mode="manual",
        manual_fwhm_mult=2.4,
    )

    assert plan.shape_px == 21
    assert plan.reason == "manual"


def test_epsf_star_plan_scales_sublinearly_and_honors_cap():
    assert plan_epsf_stars(100, 100).target == 20
    assert plan_epsf_stars(500, 500).target == 45
    assert plan_epsf_stars(2000, 2000).target == 89
    assert plan_epsf_stars(10000, 10000).target == 120
    assert plan_epsf_stars(2000, 2000, user_cap=30).target == 30
    assert plan_epsf_stars(100, 100, user_cap=1).target == 3
    assert plan_epsf_stars(100, 7).target == 7


def test_spatial_selection_round_robins_detector_cells():
    dense = np.column_stack([np.arange(20, dtype=float), np.arange(20, dtype=float)])
    sparse = np.array([[75.0, 20.0], [20.0, 75.0], [75.0, 75.0]])
    xy = np.vstack([dense, sparse])
    scores = np.arange(len(xy), 0, -1, dtype=float)

    selected = select_spatially_balanced(
        xy,
        scores,
        target=4,
        image_shape=(100, 100),
        grid_size=2,
    )

    assert selected.tolist() == [0, 20, 21, 22]


def test_epsf_nearest_distance_uses_all_detections_not_only_candidates():
    candidates = np.array([[10.0, 10.0], [80.0, 80.0]])
    all_sources = np.array([[10.0, 10.0], [12.0, 10.0], [50.0, 50.0], [80.0, 80.0]])

    distances = nearest_other_source_distance(candidates, all_sources)

    assert np.isclose(distances[0], 2.0)
    assert np.isclose(distances[1], np.hypot(30.0, 30.0))


def test_epsf_annulus_contamination_detects_neighbour_residual():
    rng = np.random.default_rng(12)
    image = rng.normal(0.0, 1.0, (101, 101))
    yy, xx = np.mgrid[:101, :101]
    for x_center in (30.0, 75.0):
        image += 500.0 * np.exp(-0.5 * ((xx - x_center) ** 2 + (yy - 50.0) ** 2) / 1.0**2)
    image += 120.0 * np.exp(-0.5 * ((xx - 36.0) ** 2 + (yy - 50.0) ** 2) / 1.0**2)

    score = measure_epsf_annulus_contamination(
        image,
        np.array([[30.0, 50.0], [75.0, 50.0]]),
        fwhm_px=2.0,
        background_rms=1.0,
    )

    assert score[0] > 5.0 * score[1]


def test_epsf_reference_selection_avoids_hidden_neighbour_and_cluster_core():
    rng = np.random.default_rng(4)
    image = rng.normal(0.0, 1.0, (100, 100))
    candidates = np.array([[50.0, 50.0], [20.0, 20.0], [80.0, 80.0]])
    all_sources = np.vstack([candidates, [[22.0, 20.0]]])

    selection = select_epsf_reference_stars(
        candidates,
        np.array([1000.0, 900.0, 800.0]),
        all_sources,
        image,
        target=1,
        image_shape=image.shape,
        grid_size=1,
        fwhm_px=2.0,
        isolation_fwhm_mult=2.0,
        background_rms=1.0,
        core_center=(50.0, 50.0),
        core_radius_px=10.0,
    )

    assert selection.selected_indices.tolist() == [2]
    assert not selection.core_safe[0]
    assert not selection.isolated[1]
    assert selection.core_safe[2]


def test_epsf_reference_target_is_a_cap_not_a_dirty_star_quota(monkeypatch):
    image = np.zeros((100, 100), dtype=float)
    candidates = np.column_stack([
        np.linspace(10.0, 90.0, 10),
        np.full(10, 50.0),
    ])
    monkeypatch.setattr(
        psf_policy,
        "measure_epsf_annulus_contamination",
        lambda *_args, **_kwargs: np.arange(10, dtype=float),
    )

    selection = select_epsf_reference_stars(
        candidates,
        np.linspace(1000.0, 2000.0, 10),
        candidates,
        image,
        target=10,
        image_shape=image.shape,
        grid_size=1,
        fwhm_px=1.0,
        isolation_fwhm_mult=2.0,
        background_rms=1.0,
        minimum_required=3,
    )

    assert len(selection.selected_indices) == 7
    assert selection.n_fallback_selected == 0


def test_epsf_reference_relaxes_morphology_when_post_pool_hits_min_but_misses_target(monkeypatch):
    points = np.array(
        [[x, y] for y in (10.0, 30.0, 50.0, 70.0) for x in (10.0, 30.0, 50.0, 70.0, 90.0)]
    )
    monkeypatch.setattr(
        psf_policy,
        "measure_epsf_annulus_contamination",
        lambda *_args, **_kwargs: np.zeros(len(points), dtype=float),
    )

    plan = plan_epsf_stars(200, len(points), minimum=3)
    assert plan.target == len(points)
    selection = select_epsf_reference_stars(
        points,
        np.linspace(1000.0, 2000.0, len(points)),
        points,
        np.zeros((100, 100), dtype=float),
        target=plan.target,
        image_shape=(100, 100),
        grid_size=2,
        fwhm_px=1.0,
        isolation_fwhm_mult=2.0,
        background_rms=1.0,
        minimum_required=3,
        morphology_ok=np.array([True, True, True] + [False] * (len(points) - 3)),
    )

    assert len(selection.selected_indices) == plan.target
    assert selection.n_morphology_relaxed_selected == len(points) - 3
    assert selection.n_fallback_selected == 0


def test_epsf_reference_clean_relaxed_outranks_dirty_strict(monkeypatch):
    points = np.array([[10.0, 10.0], [30.0, 10.0], [50.0, 10.0], [70.0, 10.0]])
    monkeypatch.setattr(
        psf_policy,
        "measure_epsf_annulus_contamination",
        lambda *_args, **_kwargs: np.array([10.0, 0.0, 0.0, 0.0]),
    )

    selection = select_epsf_reference_stars(
        points,
        np.ones(4),
        points,
        np.zeros((100, 100), dtype=float),
        target=2,
        image_shape=(100, 100),
        grid_size=1,
        fwhm_px=1.0,
        isolation_fwhm_mult=2.0,
        background_rms=1.0,
        minimum_required=2,
        morphology_ok=np.array([True, False, False, True]),
    )

    assert 0 not in selection.selected_indices  # dirty strict candidate
    assert 1 in selection.selected_indices      # clean relaxed candidate
    assert selection.n_morphology_relaxed_selected == 1
    assert selection.n_fallback_selected == 0


def test_local_group_policy_limits_optimizer_blast_radius():
    assert local_group_policy(1000, enabled=True, requested_max_size=20) == (3, 100)
    assert local_group_policy(5000, enabled=True, requested_max_size=4) == (3, 200)
    assert local_group_policy(20, enabled=True, requested_max_size=3) == (3, 3)
    assert local_group_policy(100, enabled=False, requested_max_size=3) == (1, 0)
    assert local_group_policy(100, enabled=True, requested_max_size=1) == (1, 0)


def test_forced_catalog_matching_anchors_detections_and_keeps_close_catalog_stars():
    merged = merge_forced_catalog_seeds(
        np.array([[10.0, 10.0], [30.0, 30.0]]),
        np.array([1, 2]),
        np.array([False, False]),
        np.array([[11.8, 10.0], [31.0, 30.0], [31.2, 30.0]]),
        np.array([100.0, 200.0, 50.0]),
        np.array([101.0, 102.0, 103.0]),
        match_radius_px=2.0,
    )

    assert merged.n_matched == 2
    assert merged.n_added == 1
    assert np.allclose(merged.xy[:2], [[11.8, 10.0], [31.0, 30.0]])
    assert merged.forced_mask.tolist() == [True, True, True]
    assert np.allclose(merged.xy[2], [31.2, 30.0])
    assert merged.flux_by_uid[1] == 100.0
    assert merged.flux_by_uid[2] == 200.0
    assert sorted(merged.flux_by_uid.values()) == [50.0, 100.0, 200.0]


def test_forced_star_cannot_steal_a_detected_neighbours_detection():
    """The tight-blend seed loss, reproduced at function level.

    One detection sits at a real, detected star R = (10, 10). A forced catalog
    star S = (17, 10) — a blend partner 7 px away, inside the 8.7 px match
    radius — has no detection of its own. Without context, S claims R's
    detection and the snap moves the seed to S: R's light is never modelled
    and nobody is told (measured on M13, 2026-08-14: 26 of 28 tight blends).
    With R passed as context it wins its own detection at distance zero, S
    falls through to the append path, and both stars get seeds.
    """
    detections = np.array([[10.0, 10.0]])
    forced = np.array([[17.0, 10.0]])

    broken = merge_forced_catalog_seeds(
        detections, np.array([1]), np.array([False]),
        forced, np.array([500.0]), np.array([np.nan]),
        match_radius_px=8.7,
    )
    # The failure mode this fix removes: one seed, relocated onto S.
    assert len(broken.xy) == 1 and np.allclose(broken.xy[0], [17.0, 10.0])

    fixed = merge_forced_catalog_seeds(
        detections, np.array([1]), np.array([False]),
        forced, np.array([500.0]), np.array([np.nan]),
        match_radius_px=8.7,
        context_xy=np.array([[10.0, 10.0]]),
    )
    assert len(fixed.xy) == 2
    assert np.allclose(fixed.xy[0], [10.0, 10.0]), "검출별이 제 검출을 지켜야 한다"
    assert np.allclose(fixed.xy[1], [17.0, 10.0]), "forced 별은 자기 씨앗을 새로 받아야 한다"
    assert fixed.n_matched == 0 and fixed.n_added == 1


def test_context_does_not_change_the_isolated_forced_case():
    """A forced star with its own detection must snap exactly as before."""
    merged = merge_forced_catalog_seeds(
        np.array([[10.0, 10.0], [40.0, 40.0]]),
        np.array([1, 2]),
        np.array([False, False]),
        np.array([[11.0, 10.0]]),
        np.array([300.0]),
        np.array([7.0]),
        match_radius_px=2.0,
        context_xy=np.array([[40.0, 40.0]]),
    )
    assert merged.n_matched == 1 and merged.n_added == 0
    assert np.allclose(merged.xy[0], [11.0, 10.0])
    assert merged.flux_by_uid[1] == 300.0


def test_matched_psf_flux_seed_recovers_total_flux():
    yy, xx = np.mgrid[-4:5, -4:5]
    kernel = np.exp(-(xx**2 + yy**2) / (2.0 * 1.2**2))
    kernel /= kernel.sum()

    def eval_psf(dx, dy):
        return np.exp(-(dx**2 + dy**2) / (2.0 * 1.2**2)) / kernel.sum()

    image = np.full((31, 31), 7.0, dtype=float)
    y_grid, x_grid = np.mgrid[:31, :31]
    image += 2500.0 * eval_psf(x_grid - 15.2, y_grid - 14.7)

    flux = estimate_psf_flux_seeds(
        image,
        np.array([[15.2, 14.7]]),
        eval_psf,
        fit_shape=9,
        fallback=np.array([100.0]),
    )

    assert np.isclose(flux[0], 2500.0, rtol=0.01)
