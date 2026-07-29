"""Adaptive resource and reference-star policies for PSF photometry."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree


@dataclass(frozen=True)
class EPSFStarPlan:
    """Per-frame ePSF reference-star selection budget."""

    n_detected: int
    n_candidates: int
    target: int
    user_cap: int
    grid_size: int

    @property
    def automatic(self) -> bool:
        return self.user_cap <= 0


@dataclass(frozen=True)
class EPSFReferenceSelection:
    """Diagnostics and selected indices for contamination-aware ePSF stars."""

    selected_indices: np.ndarray
    nearest_neighbor_px: np.ndarray
    contamination_score: np.ndarray
    core_distance_px: np.ndarray
    isolated: np.ndarray
    low_contamination: np.ndarray
    core_safe: np.ndarray
    quality_score: np.ndarray
    selection_tier: np.ndarray
    morphology_ok: np.ndarray | None = None

    @property
    def n_isolated(self) -> int:
        return int(np.count_nonzero(self.isolated))

    @property
    def n_low_contamination(self) -> int:
        return int(np.count_nonzero(self.low_contamination))

    @property
    def n_core_rejected(self) -> int:
        return int(np.count_nonzero(~self.core_safe))

    @property
    def n_fallback_selected(self) -> int:
        selected_tiers = self.selection_tier[self.selection_tier >= 0]
        return int(np.count_nonzero(selected_tiers > 0))

    @property
    def n_morphology_relaxed_selected(self) -> int:
        """Number of selected stars admitted by relaxing morphology cuts."""
        selected = self.selection_tier >= 0
        morphology = (
            self.morphology_ok
            if self.morphology_ok is not None
            else np.ones(len(self.selection_tier), dtype=bool)
        )
        return int(np.count_nonzero(selected & ~morphology))


@dataclass(frozen=True)
class ForcedSeedMerge:
    """Initial Step 8 sources after one-to-one forced-catalog matching."""

    xy: np.ndarray
    det_uids: np.ndarray
    forced_mask: np.ndarray
    flux_by_uid: dict[int, float]
    n_matched: int
    n_added: int


@dataclass(frozen=True)
class PSFFitWindowPlan:
    """Per-frame fit footprint selected from the measured PSF energy."""

    mode: str
    shape_px: int
    energy_fraction: float
    target_energy_fraction: float
    noise_equivalent_area_px: float
    reason: str


def _odd_size(value: float, *, minimum: int, maximum: int) -> int:
    size = max(int(minimum), min(int(maximum), int(round(float(value)))))
    if size % 2 == 0:
        size = size + 1 if size < maximum else size - 1
    return max(3, size)


def plan_psf_fit_window(
    native_psf: np.ndarray,
    fwhm_px: float,
    *,
    mode: str = "auto",
    manual_fwhm_mult: float = 2.4,
    target_energy_fraction: float = 0.90,
    minimum_fwhm_mult: float = 2.0,
    maximum_size_px: int = 31,
) -> PSFFitWindowPlan:
    """Choose a stable local fit window and report its captured PSF energy.

    Automatic mode finds the smallest centered square that reaches the target
    positive PSF energy. A two-FWHM floor avoids the strong source/sky
    degeneracy measured for undersized footprints.
    """

    psf = np.asarray(native_psf, dtype=float)
    requested_mode = str(mode).strip().lower()
    mode_safe = requested_mode if requested_mode in {"auto", "manual"} else "auto"
    fwhm_safe = max(float(fwhm_px), 1.0)
    max_size = max(9, int(maximum_size_px))
    if max_size % 2 == 0:
        max_size -= 1
    if psf.ndim == 2 and min(psf.shape) >= 3:
        support = min(max_size, psf.shape[0], psf.shape[1])
        if support % 2 == 0:
            support -= 1
        max_size = max(3, support)

    manual_size = _odd_size(
        max(1.0, float(manual_fwhm_mult)) * fwhm_safe,
        minimum=9,
        maximum=max_size,
    )
    positive = np.where(np.isfinite(psf) & (psf > 0), psf, 0.0)
    total = float(np.sum(positive))
    normalized = positive / total if total > 0 else positive
    normalized_sum_sq = float(np.sum(normalized * normalized))
    nea = 1.0 / normalized_sum_sq if normalized_sum_sq > 0 else np.nan

    def energy_for(size: int) -> float:
        if total <= 0 or psf.ndim != 2:
            return np.nan
        cy = psf.shape[0] // 2
        cx = psf.shape[1] // 2
        half = int(size) // 2
        cutout = normalized[
            max(0, cy - half):min(psf.shape[0], cy + half + 1),
            max(0, cx - half):min(psf.shape[1], cx + half + 1),
        ]
        return float(np.sum(cutout))

    if mode_safe == "manual" or total <= 0 or psf.ndim != 2:
        reason = "manual" if mode_safe == "manual" else "invalid_psf_fallback"
        return PSFFitWindowPlan(
            mode=mode_safe,
            shape_px=manual_size,
            energy_fraction=energy_for(manual_size),
            target_energy_fraction=float(target_energy_fraction),
            noise_equivalent_area_px=nea,
            reason=reason,
        )

    target = min(0.995, max(0.50, float(target_energy_fraction)))
    minimum_size = _odd_size(
        max(1.0, float(minimum_fwhm_mult)) * fwhm_safe,
        minimum=9,
        maximum=max_size,
    )
    selected = max_size
    for size in range(minimum_size, max_size + 1, 2):
        selected = size
        if energy_for(size) >= target:
            break
    fraction = energy_for(selected)
    reason = "target_reached" if np.isfinite(fraction) and fraction >= target else "size_cap"
    return PSFFitWindowPlan(
        mode="auto",
        shape_px=int(selected),
        energy_fraction=fraction,
        target_energy_fraction=target,
        noise_equivalent_area_px=nea,
        reason=reason,
    )


def merge_forced_catalog_seeds(
    detection_xy: np.ndarray,
    detection_uids: np.ndarray,
    detection_forced_mask: np.ndarray,
    forced_xy: np.ndarray,
    forced_flux: np.ndarray,
    forced_master_ids: np.ndarray,
    *,
    match_radius_px: float,
) -> ForcedSeedMerge:
    """Anchor matched detections and retain every unmatched catalog source.

    Matches are globally greedy in increasing separation and one-to-one. This
    prevents a displaced blend centroid from coexisting with a duplicate
    catalog seed while preserving multiple close catalog stars as separate
    fixed-position sources.
    """

    xy = np.asarray(detection_xy, dtype=float)
    uids = np.asarray(detection_uids, dtype=int)
    forced_mask = np.asarray(detection_forced_mask, dtype=bool)
    catalog_xy = np.asarray(forced_xy, dtype=float)
    catalog_flux = np.asarray(forced_flux, dtype=float)
    master_ids = np.asarray(forced_master_ids, dtype=float)
    if xy.ndim != 2 or xy.shape[1] != 2:
        raise ValueError("detection_xy must have shape (N, 2)")
    if catalog_xy.ndim != 2 or catalog_xy.shape[1] != 2:
        raise ValueError("forced_xy must have shape (M, 2)")
    if uids.shape != (len(xy),) or forced_mask.shape != (len(xy),):
        raise ValueError("detection arrays must have matching lengths")
    if catalog_flux.shape != (len(catalog_xy),) or master_ids.shape != (len(catalog_xy),):
        raise ValueError("forced-catalog arrays must have matching lengths")

    xy = xy.copy()
    uids = uids.copy()
    forced_mask = forced_mask.copy()
    radius = max(0.0, float(match_radius_px))
    pairs: list[tuple[float, int, int]] = []
    if len(xy) and len(catalog_xy) and radius > 0:
        tree = cKDTree(xy)
        for catalog_index, point in enumerate(catalog_xy):
            for detection_index in tree.query_ball_point(point, r=radius):
                distance = float(np.hypot(*(xy[detection_index] - point)))
                pairs.append((distance, catalog_index, int(detection_index)))

    matched_catalog: set[int] = set()
    matched_detection: set[int] = set()
    flux_by_uid: dict[int, float] = {}
    for _, catalog_index, detection_index in sorted(pairs):
        if catalog_index in matched_catalog or detection_index in matched_detection:
            continue
        xy[detection_index] = catalog_xy[catalog_index]
        forced_mask[detection_index] = True
        flux_by_uid[int(uids[detection_index])] = float(catalog_flux[catalog_index])
        matched_catalog.add(catalog_index)
        matched_detection.add(detection_index)

    appended_xy: list[list[float]] = []
    appended_uids: list[int] = []
    used_uids = set(map(int, uids))
    for catalog_index, point in enumerate(catalog_xy):
        if catalog_index in matched_catalog:
            continue
        master_id = master_ids[catalog_index]
        if np.isfinite(master_id):
            uid = -1000000 - int(master_id)
        else:
            uid = -2000000 - catalog_index
        while uid in used_uids:
            uid -= 1
        used_uids.add(uid)
        appended_xy.append([float(point[0]), float(point[1])])
        appended_uids.append(uid)
        flux_by_uid[uid] = float(catalog_flux[catalog_index])

    if appended_xy:
        xy = np.vstack([xy, np.asarray(appended_xy, dtype=float)])
        uids = np.concatenate([uids, np.asarray(appended_uids, dtype=int)])
        forced_mask = np.concatenate([
            forced_mask,
            np.ones(len(appended_xy), dtype=bool),
        ])

    return ForcedSeedMerge(
        xy=xy,
        det_uids=uids,
        forced_mask=forced_mask,
        flux_by_uid=flux_by_uid,
        n_matched=len(matched_catalog),
        n_added=len(appended_xy),
    )


def plan_epsf_stars(
    n_detected: int,
    n_candidates: int,
    *,
    user_cap: int = 0,
    minimum: int = 15,
    maximum: int = 120,
) -> EPSFStarPlan:
    """Scale the ePSF-star budget sub-linearly with the frame source count.

    ``user_cap=0`` selects the automatic ceiling. A positive value is a hard
    resource cap, not a request to use low-quality candidates to fill a quota.
    """

    n_detected = max(0, int(n_detected))
    n_candidates = max(0, int(n_candidates))
    minimum = max(3, int(minimum))
    maximum = max(minimum, int(maximum))
    user_cap = max(0, int(user_cap))

    if n_candidates == 0:
        target = 0
    else:
        automatic_target = int(round(2.0 * np.sqrt(max(1, n_detected))))
        target = min(n_candidates, max(minimum, automatic_target), maximum)
        if user_cap > 0:
            target = min(target, max(3, user_cap))

    if target >= 72:
        grid_size = 3
    elif target >= 8:
        grid_size = 2
    else:
        grid_size = 1

    return EPSFStarPlan(
        n_detected=n_detected,
        n_candidates=n_candidates,
        target=target,
        user_cap=user_cap,
        grid_size=grid_size,
    )


def select_spatially_balanced(
    xy: np.ndarray,
    scores: np.ndarray,
    *,
    target: int,
    image_shape: tuple[int, int],
    grid_size: int,
) -> np.ndarray:
    """Return deterministic candidate indices balanced across the detector.

    Candidates are ranked by quality within each grid cell and selected in
    rounds so a dense detector region cannot monopolize the ePSF model.
    """

    points = np.asarray(xy, dtype=float)
    quality = np.asarray(scores, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("xy must have shape (N, 2)")
    if quality.ndim != 1 or len(quality) != len(points):
        raise ValueError("scores must have shape (N,)")

    target = min(max(0, int(target)), len(points))
    if target == 0:
        return np.zeros(0, dtype=int)

    height, width = int(image_shape[0]), int(image_shape[1])
    if height <= 0 or width <= 0:
        raise ValueError("image_shape must be positive")
    grid_size = max(1, int(grid_size))

    finite_xy = np.isfinite(points[:, 0]) & np.isfinite(points[:, 1])
    safe_scores = np.where(np.isfinite(quality), quality, -np.inf)
    x_cell = np.clip((points[:, 0] / max(1.0, float(width))) * grid_size, 0, grid_size - 1).astype(int)
    y_cell = np.clip((points[:, 1] / max(1.0, float(height))) * grid_size, 0, grid_size - 1).astype(int)
    cell_id = y_cell * grid_size + x_cell

    queues: list[list[int]] = []
    for cid in range(grid_size * grid_size):
        members = np.flatnonzero(finite_xy & (cell_id == cid))
        ordered = sorted(members.tolist(), key=lambda idx: (-safe_scores[idx], idx))
        queues.append(ordered)

    selected: list[int] = []
    cursor = [0] * len(queues)
    while len(selected) < target:
        added = False
        for qidx, queue in enumerate(queues):
            pos = cursor[qidx]
            if pos >= len(queue):
                continue
            selected.append(queue[pos])
            cursor[qidx] += 1
            added = True
            if len(selected) >= target:
                break
        if not added:
            break

    return np.asarray(selected, dtype=int)


def nearest_other_source_distance(candidate_xy: np.ndarray, all_xy: np.ndarray) -> np.ndarray:
    """Distance from each candidate to the nearest other detected source."""
    candidates = np.asarray(candidate_xy, dtype=float)
    sources = np.asarray(all_xy, dtype=float)
    if candidates.ndim != 2 or candidates.shape[1] != 2:
        raise ValueError("candidate_xy must have shape (N, 2)")
    if sources.ndim != 2 or sources.shape[1] != 2:
        raise ValueError("all_xy must have shape (M, 2)")
    if len(candidates) == 0:
        return np.zeros(0, dtype=float)
    finite_sources = np.isfinite(sources[:, 0]) & np.isfinite(sources[:, 1])
    sources = sources[finite_sources]
    if len(sources) < 2:
        return np.full(len(candidates), np.inf, dtype=float)

    distances = cKDTree(sources).query(candidates, k=2, workers=1)[0]
    nearest = distances[:, 0].copy()
    self_match = nearest <= 1e-6
    nearest[self_match] = distances[self_match, 1]
    nearest[~np.isfinite(nearest)] = np.inf
    return nearest


def measure_epsf_annulus_contamination(
    image: np.ndarray,
    xy: np.ndarray,
    *,
    fwhm_px: float,
    background_rms: float,
    inner_fwhm: float = 2.0,
    outer_fwhm: float = 3.5,
) -> np.ndarray:
    """Measure local residual/background contamination around ePSF candidates.

    The score combines annulus noise excess, quadrant imbalance, and residual
    background relative to the candidate peak. Only relative ranking within a
    frame is used by the selection policy.
    """
    data = np.asarray(image, dtype=float)
    points = np.asarray(xy, dtype=float)
    if data.ndim != 2:
        raise ValueError("image must be two-dimensional")
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("xy must have shape (N, 2)")
    fwhm = max(float(fwhm_px), 1.0)
    rms = float(background_rms)
    if not np.isfinite(rms) or rms <= 0:
        sample = data.ravel()[::max(1, data.size // 65536)]
        sample = sample[np.isfinite(sample)]
        median = float(np.median(sample)) if sample.size else 0.0
        rms = float(1.4826 * np.median(np.abs(sample - median))) if sample.size else 1.0
    rms = max(rms, 1e-6)
    inner = max(1.0, float(inner_fwhm) * fwhm)
    outer = max(inner + 2.0, float(outer_fwhm) * fwhm)
    half = int(np.ceil(outer))
    height, width = data.shape
    scores = np.full(len(points), np.inf, dtype=float)

    for index, (x_center, y_center) in enumerate(points):
        if not np.isfinite(x_center) or not np.isfinite(y_center):
            continue
        x_mid = int(round(float(x_center)))
        y_mid = int(round(float(y_center)))
        x0, x1 = max(0, x_mid - half), min(width, x_mid + half + 1)
        y0, y1 = max(0, y_mid - half), min(height, y_mid + half + 1)
        if x1 - x0 < 5 or y1 - y0 < 5:
            continue
        yy, xx = np.mgrid[y0:y1, x0:x1]
        radius = np.hypot(xx - float(x_center), yy - float(y_center))
        ring = (radius >= inner) & (radius <= outer) & np.isfinite(data[y0:y1, x0:x1])
        values = data[y0:y1, x0:x1][ring]
        if values.size < 20:
            continue

        local_median = float(np.median(values))
        local_sigma = float(1.4826 * np.median(np.abs(values - local_median)))
        center_value = data[
            int(np.clip(y_mid, 0, height - 1)),
            int(np.clip(x_mid, 0, width - 1)),
        ]
        signal_scale = max(abs(float(center_value)) if np.isfinite(center_value) else 0.0, rms)

        quadrant_medians = []
        for x_positive, y_positive in ((False, False), (True, False), (False, True), (True, True)):
            quadrant = ring.copy()
            quadrant &= (xx >= x_center) if x_positive else (xx < x_center)
            quadrant &= (yy >= y_center) if y_positive else (yy < y_center)
            quadrant_values = data[y0:y1, x0:x1][quadrant]
            if quadrant_values.size:
                quadrant_medians.append(float(np.median(quadrant_values)))
        quadrant_span = (
            max(quadrant_medians) - min(quadrant_medians)
            if len(quadrant_medians) >= 2
            else 0.0
        )
        noise_excess = max(0.0, local_sigma / rms - 1.0)
        asymmetry = quadrant_span / signal_scale
        background_offset = abs(local_median) / signal_scale
        scores[index] = noise_excess + 5.0 * asymmetry + 2.0 * background_offset

    return scores


def _unit_rank(values: np.ndarray, *, lower_is_better: bool) -> np.ndarray:
    data = np.asarray(values, dtype=float)
    output = np.zeros(len(data), dtype=float)
    finite = np.flatnonzero(np.isfinite(data))
    if len(finite) == 0:
        return output
    order = finite[np.argsort(data[finite], kind="stable")]
    ranks = np.linspace(1.0, 0.0, len(order)) if lower_is_better else np.linspace(0.0, 1.0, len(order))
    output[order] = ranks
    return output


def select_epsf_reference_stars(
    candidate_xy: np.ndarray,
    candidate_flux: np.ndarray,
    all_xy: np.ndarray,
    image: np.ndarray,
    *,
    target: int,
    image_shape: tuple[int, int],
    grid_size: int,
    fwhm_px: float,
    isolation_fwhm_mult: float,
    background_rms: float,
    core_center: tuple[float, float] = (np.nan, np.nan),
    core_radius_px: float = np.nan,
    minimum_required: int = 3,
    morphology_ok: np.ndarray | None = None,
) -> EPSFReferenceSelection:
    """Select spatially balanced ePSF stars while avoiding crowded cutouts.

    The input is the pre-morphology candidate pool. Strict-morphology stars
    that are clean are selected first, followed by clean stars admitted by the
    morphology relaxation. Contamination tiers are used only to reach the
    hard minimum, preserving clean relaxed stars over dirty strict stars.
    """
    points = np.asarray(candidate_xy, dtype=float)
    flux = np.asarray(candidate_flux, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("candidate_xy must have shape (N, 2)")
    if flux.shape != (len(points),):
        raise ValueError("candidate_flux must have shape (N,)")
    if morphology_ok is None:
        morphology = np.ones(len(points), dtype=bool)
    else:
        morphology = np.asarray(morphology_ok, dtype=bool)
        if morphology.shape != (len(points),):
            raise ValueError("morphology_ok must have shape (N,)")

    target = min(max(0, int(target)), len(points))
    minimum_required = min(target, max(3, int(minimum_required)))
    nearest = nearest_other_source_distance(points, all_xy)
    contamination = measure_epsf_annulus_contamination(
        image,
        points,
        fwhm_px=fwhm_px,
        background_rms=background_rms,
    )
    isolation_radius = max(1.0, float(isolation_fwhm_mult) * max(float(fwhm_px), 1.0))
    isolated = nearest > isolation_radius

    finite_contamination = contamination[np.isfinite(contamination)]
    contamination_limit = (
        float(np.percentile(finite_contamination, 70.0))
        if finite_contamination.size
        else np.inf
    )
    low_contamination = np.isfinite(contamination) & (contamination <= contamination_limit)

    center_x, center_y = float(core_center[0]), float(core_center[1])
    has_core = (
        np.isfinite(center_x)
        and np.isfinite(center_y)
        and np.isfinite(core_radius_px)
        and float(core_radius_px) > 0
    )
    if has_core:
        core_distance = np.hypot(points[:, 0] - center_x, points[:, 1] - center_y)
        core_safe = core_distance >= float(core_radius_px)
        radial_safety = np.clip(
            (core_distance / float(core_radius_px) - 1.0) / 2.0,
            0.0,
            1.0,
        )
    else:
        core_distance = np.full(len(points), np.nan, dtype=float)
        core_safe = np.ones(len(points), dtype=bool)
        radial_safety = np.ones(len(points), dtype=float)

    separation_score = np.clip(nearest / isolation_radius, 0.0, 2.0)
    flux_rank = _unit_rank(flux, lower_is_better=False)
    clean_rank = _unit_rank(contamination, lower_is_better=True)
    quality = 4.0 * separation_score + flux_rank + 3.0 * clean_rank + 1.5 * radial_safety
    quality[~np.isfinite(points[:, 0]) | ~np.isfinite(points[:, 1])] = -np.inf

    selected: list[int] = []
    selected_mask = np.zeros(len(points), dtype=bool)
    selection_tier = np.full(len(points), -1, dtype=int)

    def take_from(mask: np.ndarray, count: int, contamination_tier: int) -> None:
        if count <= 0:
            return
        members = np.flatnonzero(mask & ~selected_mask)
        if members.size == 0:
            return
        chosen_local = select_spatially_balanced(
            points[members],
            quality[members],
            target=min(int(count), len(members)),
            image_shape=image_shape,
            grid_size=grid_size,
        )
        chosen = members[chosen_local]
        selected.extend(chosen.tolist())
        selected_mask[chosen] = True
        selection_tier[chosen] = contamination_tier

    clean = isolated & low_contamination & core_safe
    # Morphology is a quality preference, not a reason to use contaminated
    # references. This ordering is the important relaxation invariant.
    take_from(clean & morphology, target - len(selected), 0)
    take_from(clean & ~morphology, target - len(selected), 0)

    # Only the hard minimum may use contamination fallbacks. Within each
    # fallback tier, retain strict morphology preference for determinism.
    fallback_masks = (
        isolated & core_safe,
        core_safe,
        np.ones(len(points), dtype=bool),
    )
    for tier, tier_mask in enumerate(fallback_masks, start=1):
        if len(selected) >= minimum_required:
            break
        take_from(tier_mask & morphology, minimum_required - len(selected), tier)
        take_from(tier_mask & ~morphology, minimum_required - len(selected), tier)

    return EPSFReferenceSelection(
        selected_indices=np.asarray(selected, dtype=int),
        nearest_neighbor_px=nearest,
        contamination_score=contamination,
        core_distance_px=core_distance,
        isolated=isolated,
        low_contamination=low_contamination,
        core_safe=core_safe,
        quality_score=quality,
        selection_tier=selection_tier,
        morphology_ok=morphology,
    )


def local_group_policy(
    n_sources: int,
    *,
    enabled: bool,
    requested_max_size: int,
    hard_max_size: int = 3,
    max_fraction: float = 0.10,
    absolute_cap: int = 200,
) -> tuple[int, int]:
    """Return ``(max_group_size, grouped-source budget)`` for one frame."""

    n_sources = max(0, int(n_sources))
    requested_max_size = int(requested_max_size)
    if not enabled or n_sources < 2 or requested_max_size <= 1:
        return 1, 0

    max_size = min(max(2, requested_max_size), max(2, int(hard_max_size)))
    fraction = min(max(float(max_fraction), 0.0), 1.0)
    absolute_cap = max(0, int(absolute_cap))
    budget = int(np.ceil(fraction * n_sources))
    if absolute_cap > 0:
        budget = min(budget, absolute_cap)
    budget = max(max_size, budget)
    budget = min(n_sources, budget)
    return max_size, budget


def estimate_psf_flux_seeds(
    image: np.ndarray,
    xy: np.ndarray,
    eval_psf,
    *,
    fit_shape: int,
    fallback: np.ndarray | None = None,
) -> np.ndarray:
    """Estimate total-flux seeds by a local matched-PSF projection.

    Residual finders usually report a peak-pixel value. Passing that value as
    total PSF flux makes the first nonlinear fit unnecessarily difficult,
    especially for broad seeing. This projection converts each local residual
    peak to the PSF model's total-flux scale after removing a local background.
    """

    data = np.asarray(image, dtype=float)
    points = np.asarray(xy, dtype=float)
    if data.ndim != 2:
        raise ValueError("image must be two-dimensional")
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("xy must have shape (N, 2)")

    if fallback is None:
        fallback_values = np.ones(len(points), dtype=float)
    else:
        fallback_values = np.asarray(fallback, dtype=float)
        if fallback_values.shape != (len(points),):
            raise ValueError("fallback must have shape (N,)")
        fallback_values = np.where(
            np.isfinite(fallback_values) & (fallback_values > 0),
            fallback_values,
            1.0,
        )

    size = max(3, int(fit_shape))
    if size % 2 == 0:
        size += 1
    half = size // 2
    height, width = data.shape
    output = fallback_values.copy()

    for index, (x_cen, y_cen) in enumerate(points):
        if not np.isfinite(x_cen) or not np.isfinite(y_cen):
            continue
        x_int = int(round(float(x_cen)))
        y_int = int(round(float(y_cen)))
        x0, x1 = max(0, x_int - half), min(width, x_int + half + 1)
        y0, y1 = max(0, y_int - half), min(height, y_int + half + 1)
        if x1 - x0 < 3 or y1 - y0 < 3:
            continue

        yy, xx = np.mgrid[y0:y1, x0:x1]
        model = np.asarray(eval_psf(xx - x_cen, yy - y_cen), dtype=float)
        patch = data[y0:y1, x0:x1].astype(float, copy=True)
        finite = np.isfinite(model) & np.isfinite(patch)
        if np.count_nonzero(finite) < 5:
            continue

        edge = np.zeros(patch.shape, dtype=bool)
        edge[0, :] = edge[-1, :] = True
        edge[:, 0] = edge[:, -1] = True
        edge_values = patch[edge & np.isfinite(patch)]
        background = float(np.median(edge_values)) if edge_values.size else 0.0
        signal = patch - background
        denominator = float(np.sum(model[finite] ** 2))
        if denominator <= 0 or not np.isfinite(denominator):
            continue
        amplitude = float(np.sum(model[finite] * signal[finite]) / denominator)
        if np.isfinite(amplitude) and amplitude > 0:
            output[index] = amplitude

    return output


def psf_symmetric_mask(
    image: np.ndarray,
    xy: np.ndarray,
    *,
    background: float | np.ndarray = 0.0,
    neighbor_frac: float = 0.3,
) -> np.ndarray:
    """Keep only sources whose profile is symmetric like a real PSF.

    A star is isotropic, so at its peak both horizontal neighbours (or both
    vertical ones) sit at a sizable fraction of the peak. Point-like noise does
    not: an isolated hot pixel has no bright neighbour at all, a two-pixel pair
    has only one side, an L-shaped or diagonal cluster fails on both axes. The
    test compares neighbours *against the peak itself*, so it is scale
    invariant — faint and bright stars pass alike, and a spike of any
    brightness fails.

    This matters for ePSF reference selection: candidates are ranked by
    "bright and isolated", which is exactly what a cosmic ray maximises
    (all its flux in one pixel, no neighbours). On CMOS detectors 16-40% of
    detections can be 1-pixel spikes, and once they enter the ePSF the model
    comes out far too narrow (measured 2.75x on M67/QHY600, dropping PSF flux
    to 32% of aperture).

    Ported from the AstralImage cosmetic engine (``_extended_bright_cores``),
    which uses the same neighbour-pair test to protect stars from hot-pixel
    correction. Here it is used to *reject* non-stars rather than protect
    stars, but the geometric argument is identical.

    Parameters
    ----------
    image :
        Frame the sources were detected on.
    xy :
        ``(N, 2)`` array of source positions (x, y), pixel coordinates.
    background :
        Scalar or per-pixel background to subtract before the comparison.
        Neighbour ratios are meaningless on an un-subtracted sky.
    neighbor_frac :
        Minimum neighbour-to-peak ratio. 0.3 accepts a PSF as narrow as
        FWHM ~1.6 px while still rejecting single-pixel spikes.

    Returns
    -------
    np.ndarray
        Boolean mask, ``True`` where the source looks like a PSF. Sources too
        close to the border to test are kept (``True``) — the caller's other
        cuts handle edges.
    """
    xy = np.asarray(xy, dtype=float)
    if xy.ndim != 2 or xy.shape[0] == 0:
        return np.zeros(0, dtype=bool)

    resid = np.asarray(image, dtype=float) - background
    ny, nx = resid.shape
    keep = np.ones(len(xy), dtype=bool)

    xi = np.rint(xy[:, 0]).astype(int)
    yi = np.rint(xy[:, 1]).astype(int)
    inside = (xi >= 1) & (xi < nx - 1) & (yi >= 1) & (yi < ny - 1)
    if not np.any(inside):
        return keep

    xs, ys = xi[inside], yi[inside]
    peak = resid[ys, xs]
    thr = float(neighbor_frac) * peak

    left = resid[ys, xs - 1] >= thr
    right = resid[ys, xs + 1] >= thr
    up = resid[ys - 1, xs] >= thr
    down = resid[ys + 1, xs] >= thr

    # A non-positive peak carries no shape information; leave it to other cuts.
    symmetric = ((left & right) | (up & down)) | ~np.isfinite(peak) | (peak <= 0)
    keep[inside] = symmetric
    return keep
