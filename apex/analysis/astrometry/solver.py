"""Internal plate solver entry point.

Pipeline
--------
1. Source selection (caller already passed flux-weighted / isolation-aware
   subset from Step 4; we further trim to the brightest ``n_brightest_src``
   for quad construction).
2. Build a seed TAN WCS from (target RA, Dec, pixel scale) and project the
   Gaia catalogue into pixel coordinates.
3. Build quad features on both sides (``quad_matcher.build_quads``).
4. Match quad codes via 4-D kd-tree and verify each candidate with a RANSAC
   inlier count under the implied similarity transform.
5. Use the RANSAC inliers as seed pairs for ``astropy.wcs.utils.fit_wcs_from_points``.
6. Sigma-clip refine: re-project Gaia under the fitted WCS, rematch with a
   tight radius, refit. Iterate.
7. Independent verification: project the *full* Gaia catalogue and rematch
   against *all* detected sources. Reject if too few survive — this guards
   against the solution snapping onto a sparse bright subset.

The coarse (rotation × scale × translation) brute-force grid search that lived
in the previous solver is removed: quad matching is shift- and rotation-
invariant by construction, so we no longer need to scan that space.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.spatial import cKDTree
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from astropy.wcs.utils import fit_wcs_from_points

from apex.utils.constants import MAD_TO_SIGMA
from .quad_matcher import build_quads, match_quads, ransac_verify_topk


# ---------------------------------------------------------------------------
# Public types
# ---------------------------------------------------------------------------

@dataclass
class WCSSolveResult:
    converged: bool
    wcs: Optional[WCS]
    n_matches: int
    rms_arcsec: float
    rotation_deg: float
    scale_arcsec_per_px: float
    flip_x: bool
    flip_y: bool
    model: str = "TAN"
    sip_order: int = 0
    log: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# WCS helpers (ported from apex/utils/wcs_solver.py)
# ---------------------------------------------------------------------------

def _seed_wcs(ra: float, dec: float, scale_arcsec: float,
              img_shape: tuple[int, int]) -> WCS:
    """TAN WCS centred at (ra, dec) with given pixel scale (arcsec/px).

    CDELT[0] is negative so increasing x corresponds to decreasing RA, which
    is the standard TAN convention used elsewhere in APEX.
    """
    ny, nx = img_shape
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crpix = [nx / 2.0, ny / 2.0]
    w.wcs.crval = [ra, dec]
    w.wcs.cdelt = [-scale_arcsec / 3600.0, scale_arcsec / 3600.0]
    w.wcs.fix()
    return w


def _project_to_pixels(ra: np.ndarray, dec: np.ndarray, wcs: WCS) -> np.ndarray:
    sky = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    x, y = wcs.world_to_pixel(sky)
    return np.column_stack([np.asarray(x, float), np.asarray(y, float)])


def _wcs_cd_matrix(wcs: WCS) -> Optional[np.ndarray]:
    try:
        if getattr(wcs.wcs, "has_cd", lambda: False)():
            return np.asarray(wcs.wcs.cd, dtype=float)
    except Exception:
        pass
    try:
        pc = np.asarray(wcs.wcs.pc, dtype=float)
        cdelt = np.asarray(wcs.wcs.cdelt, dtype=float)
        if pc.shape == (2, 2) and cdelt.shape[0] >= 2:
            return pc * cdelt[:, np.newaxis]
    except Exception:
        pass
    return None


def _match_unique(
    src_xy: np.ndarray,
    cat_xy_t: np.ndarray,
    match_radius: float,
):
    """Greedy 1:1 nearest-neighbour matcher (refine path).

    Identical semantics to the helper used in the previous solver — see
    ``apex/utils/wcs_solver.py`` history. Centralised here so the refine
    path stays self-contained.
    """
    src_xy = np.asarray(src_xy, dtype=float)
    cat_xy_t = np.asarray(cat_xy_t, dtype=float)
    if len(src_xy) == 0 or len(cat_xy_t) == 0:
        empty = np.asarray([], dtype=int)
        return empty, empty, np.asarray([], dtype=float)

    valid_src = np.isfinite(src_xy).all(axis=1)
    valid_cat = np.isfinite(cat_xy_t).all(axis=1)
    if not valid_src.any() or not valid_cat.any():
        empty = np.asarray([], dtype=int)
        return empty, empty, np.asarray([], dtype=float)

    src_idx_all = np.flatnonzero(valid_src)
    cat_idx_all = np.flatnonzero(valid_cat)
    tree = cKDTree(cat_xy_t[valid_cat])
    dists, idxs = tree.query(src_xy[valid_src], k=1)
    ok = np.isfinite(dists) & (dists <= float(match_radius))
    if not ok.any():
        empty = np.asarray([], dtype=int)
        return empty, empty, np.asarray([], dtype=float)

    cand_src = src_idx_all[ok]
    cand_cat = cat_idx_all[np.asarray(idxs[ok], dtype=int)]
    cand_dist = np.asarray(dists[ok], dtype=float)
    order = np.argsort(cand_dist)

    used_cat: set[int] = set()
    keep_src: list[int] = []
    keep_cat: list[int] = []
    keep_dist: list[float] = []
    for pos in order:
        cat_i = int(cand_cat[pos])
        if cat_i in used_cat:
            continue
        used_cat.add(cat_i)
        keep_src.append(int(cand_src[pos]))
        keep_cat.append(cat_i)
        keep_dist.append(float(cand_dist[pos]))

    return (
        np.asarray(keep_src, dtype=int),
        np.asarray(keep_cat, dtype=int),
        np.asarray(keep_dist, dtype=float),
    )


def _sigma_clip_refine(
    src_xy: np.ndarray,
    cat_radec: np.ndarray,
    wcs: WCS,
    *,
    sigma: float = 3.0,
    n_iter: int = 4,
    min_pairs: int = 6,
    sip_degree: int | None = None,
    sip_min_pairs: int = 30,
):
    """Iterative re-match + refit. Returns (wcs, src_matched, cat_matched_radec).

    Each iteration:
      1. Project the catalogue under the current WCS.
      2. Compute the median nearest-neighbour distance to set an adaptive
         match radius (clipped to [3, 80] px).
      3. Greedy 1:1 match.
      4. Robust outlier cut at ``median + sigma · MAD_to_σ · MAD``.
      5. Refit WCS via ``fit_wcs_from_points``.
    """
    src_m = np.zeros((0, 2))
    cat_m = np.zeros((0, 2))
    for _ in range(n_iter):
        proj = _project_to_pixels(cat_radec[:, 0], cat_radec[:, 1], wcs)
        valid = np.isfinite(proj).all(axis=1)
        if not valid.any():
            break
        tree = cKDTree(proj[valid])
        d0, _ = tree.query(src_xy, k=1)
        finite_d = d0[np.isfinite(d0)]
        if len(finite_d) == 0:
            break
        best_r = float(np.clip(np.median(finite_d) * 3.0, 3.0, 80.0))
        src_idx, cat_idx, dists = _match_unique(src_xy, proj, best_r)
        if len(src_idx) < min_pairs:
            break

        med = float(np.median(dists))
        mad = float(np.median(np.abs(dists - med)))
        scatter = MAD_TO_SIGMA * mad if mad > 0 else float(np.std(dists))
        thresh = max(3.0, med + sigma * scatter)
        keep = dists <= thresh
        if int(keep.sum()) < min_pairs:
            break

        src_in = src_xy[src_idx[keep]]
        cat_in = cat_radec[cat_idx[keep]]
        src_m, cat_m = src_in, cat_in

        sky_in = SkyCoord(cat_in[:, 0] * u.deg, cat_in[:, 1] * u.deg)
        # SIP captures TAN residual + optical distortion once enough pairs are
        # locked in. Below sip_min_pairs the SIP fit is under-determined and
        # would oscillate, so stick to plain TAN for early iterations.
        use_sip = (
            sip_degree is not None and int(sip_degree) > 0
            and len(src_in) >= int(sip_min_pairs)
        )
        fit_kwargs = {"sip_degree": int(sip_degree)} if use_sip else {}
        try:
            wcs = fit_wcs_from_points(
                (src_in[:, 0], src_in[:, 1]),
                sky_in,
                proj_point=SkyCoord(
                    wcs.wcs.crval[0] * u.deg, wcs.wcs.crval[1] * u.deg
                ),
                projection="TAN",
                **fit_kwargs,
            )
        except Exception:
            # If SIP fit blew up, retry once with plain TAN before giving up.
            if use_sip:
                try:
                    wcs = fit_wcs_from_points(
                        (src_in[:, 0], src_in[:, 1]),
                        sky_in,
                        proj_point=SkyCoord(
                            wcs.wcs.crval[0] * u.deg, wcs.wcs.crval[1] * u.deg
                        ),
                        projection="TAN",
                    )
                except Exception:
                    break
            else:
                break

    return wcs, src_m, cat_m


def _fit_tan_wcs(src_xy: np.ndarray, sky: SkyCoord, proj_point: SkyCoord) -> WCS:
    return fit_wcs_from_points(
        (src_xy[:, 0], src_xy[:, 1]),
        sky,
        proj_point=proj_point,
        projection="TAN",
    )


def _fit_sip_wcs(
    src_xy: np.ndarray,
    sky: SkyCoord,
    proj_point: SkyCoord,
    sip_degree: int,
) -> WCS:
    return fit_wcs_from_points(
        (src_xy[:, 0], src_xy[:, 1]),
        sky,
        proj_point=proj_point,
        projection="TAN",
        sip_degree=int(sip_degree),
    )


def _pair_rms_arcsec(wcs: WCS, src_xy: np.ndarray, sky: SkyCoord) -> float:
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sky_pred = SkyCoord(wcs.pixel_to_world(src_xy[:, 0], src_xy[:, 1]))
            sep_arcsec = sky_pred.separation(sky).arcsec
        return float(np.sqrt(np.nanmean(sep_arcsec ** 2)))
    except Exception:
        return float("nan")


def _spatial_holdout_mask(
    src_xy: np.ndarray,
    *,
    fraction: float,
    min_holdout: int,
    min_train: int,
) -> np.ndarray | None:
    """Deterministic spatially spread train/holdout split for SIP validation."""
    n = int(len(src_xy))
    if n <= 0:
        return None
    frac = float(np.clip(fraction, 0.05, 0.50))
    n_hold = max(int(min_holdout), int(round(n * frac)))
    if n - n_hold < int(min_train):
        n_hold = n - int(min_train)
    if n_hold < int(min_holdout):
        return None

    xy = np.asarray(src_xy, dtype=float)
    cx = float(np.nanmedian(xy[:, 0]))
    cy = float(np.nanmedian(xy[:, 1]))
    radius = np.hypot(xy[:, 0] - cx, xy[:, 1] - cy)
    angle = np.arctan2(xy[:, 1] - cy, xy[:, 0] - cx)
    order = np.lexsort((radius, angle))
    pos = np.unique(np.rint(np.linspace(0, len(order) - 1, n_hold)).astype(int))
    hold_idx = [int(order[p]) for p in pos[:n_hold]]
    if len(hold_idx) < n_hold:
        used = set(hold_idx)
        for idx in order:
            i = int(idx)
            if i in used:
                continue
            hold_idx.append(i)
            used.add(i)
            if len(hold_idx) >= n_hold:
                break
    if len(hold_idx) < int(min_holdout):
        return None

    train_mask = np.ones(n, dtype=bool)
    train_mask[np.asarray(hold_idx, dtype=int)] = False
    if int(train_mask.sum()) < int(min_train):
        return None
    return train_mask


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------

def solve(
    sources_xy: np.ndarray,
    source_flux: Optional[np.ndarray],
    gaia_ra: np.ndarray,
    gaia_dec: np.ndarray,
    gaia_flux: Optional[np.ndarray],
    approx_ra: float,
    approx_dec: float,
    approx_scale_arcsec: float,
    img_shape: tuple[int, int],
    *,
    n_brightest_src: int = 200,
    n_brightest_cat: int = 500,
    min_matches: int = 8,
    refine_sigma: float = 3.0,
    quad_k_neighbor: int = 7,
    quad_neighbor_pool_factor: int = 3,
    quad_max_per_side: int = 1500,
    quad_code_tol: float = 0.020,
    quad_scale_ratio_tol: float = 0.25,
    ransac_inlier_radius_px: float = 4.0,
    ransac_max_trials: int = 1500,
    ransac_keep_candidates: int = 8,
    allow_reflection: bool = True,
    local_blind: bool = False,
    local_blind_radius_factor: float = 2.5,
    sip_degree: int = 3,
    sip_min_pairs: int = 30,
    sip_holdout_fraction: float = 0.25,
    sip_min_holdout_pairs: int = 8,
    sip_min_improvement: float = 0.10,
    gaia_pmra: Optional[np.ndarray] = None,
    gaia_pmdec: Optional[np.ndarray] = None,
    obs_epoch_jyear: Optional[float] = None,
    gaia_ref_epoch_jyear: float = 2016.0,
    log_fn=None,
    # Back-compat shims (ignored — the old solver took these):
    rotation_step_deg: float | None = None,
    scale_search_range: tuple[float, float] | None = None,
    n_scale_steps: int | None = None,
    match_radius_factor: float | None = None,
    translation_search: bool | None = None,
    max_translation_frac: float | None = None,
    translation_top_n: int | None = None,
) -> WCSSolveResult:
    """Solve WCS using astnet-style quad hashing + RANSAC + Gaia refine.

    Caller responsibilities
    -----------------------
    * Pass already-quality-filtered detections (Step 4 anchor/isolation logic
      handles this in step5_wcs_plate_solving.py).
    * Pass a Gaia subset that covers the target region with some margin. Normal
      mode re-filters to stars projected into the seed-WCS field; local-blind
      mode keeps a larger tangent-plane window around the target and lets
      quad/RANSAC solve translation.

    The legacy keyword arguments (``rotation_step_deg`` etc.) are accepted
    for call-site compatibility but ignored; they belonged to the brute-force
    grid search that has been removed.
    """
    del scale_search_range, n_scale_steps, match_radius_factor
    del translation_search, max_translation_frac, translation_top_n
    del rotation_step_deg

    log: list[str] = []

    def _log(msg: str) -> None:
        log.append(msg)
        if log_fn is not None:
            log_fn(msg)

    FAILED = WCSSolveResult(
        converged=False, wcs=None, n_matches=0,
        rms_arcsec=np.nan, rotation_deg=np.nan,
        scale_arcsec_per_px=np.nan, flip_x=False, flip_y=False,
        model="NONE", sip_order=0, log=log,
    )

    # ── Validate inputs ─────────────────────────────────────────────────────
    sources_xy = np.asarray(sources_xy, dtype=float)
    gaia_ra = np.asarray(gaia_ra, dtype=float)
    gaia_dec = np.asarray(gaia_dec, dtype=float)
    src_flux_arr = (
        np.asarray(source_flux, dtype=float)
        if source_flux is not None and len(source_flux) == len(sources_xy)
        else None
    )
    cat_flux_arr = (
        np.asarray(gaia_flux, dtype=float)
        if gaia_flux is not None and len(gaia_flux) == len(gaia_ra)
        else None
    )

    valid_src = np.isfinite(sources_xy).all(axis=1)
    sources_xy = sources_xy[valid_src]
    if src_flux_arr is not None:
        src_flux_arr = src_flux_arr[valid_src]

    valid_cat = np.isfinite(gaia_ra) & np.isfinite(gaia_dec)
    gaia_ra = gaia_ra[valid_cat]
    gaia_dec = gaia_dec[valid_cat]
    if cat_flux_arr is not None:
        cat_flux_arr = cat_flux_arr[valid_cat]

    if len(sources_xy) < min_matches or len(gaia_ra) < min_matches:
        _log(
            f"[WCS] Too few inputs: sources={len(sources_xy)} "
            f"gaia={len(gaia_ra)}"
        )
        return FAILED

    # ── Gaia proper-motion correction (linear, decade-scale) ───────────────
    # Gaia DR3 stores (pmra, pmdec) in mas/yr at reference epoch J2016.0, with
    # pmra already including the cos(dec) factor (i.e. it's the local east
    # angular velocity). Updating from ref_epoch to obs_epoch is a single
    # linear step:
    #   dec(t) = dec_0 + pmdec · dt
    #   ra(t)  = ra_0  + (pmra · dt) / cos(dec_0)
    # NaN PM entries are zero-filled by the caller, so unmeasured-PM stars
    # contribute as if stationary — safer than dropping them.
    if (
        obs_epoch_jyear is not None
        and gaia_pmra is not None and gaia_pmdec is not None
        and len(gaia_pmra) == len(gaia_ra)
        and len(gaia_pmdec) == len(gaia_ra)
    ):
        dt_yr = float(obs_epoch_jyear) - float(gaia_ref_epoch_jyear)
        if abs(dt_yr) > 1e-3:
            pmra_arr = np.asarray(gaia_pmra, dtype=float)[valid_cat]
            pmdec_arr = np.asarray(gaia_pmdec, dtype=float)[valid_cat]
            pmra_arr = np.where(np.isfinite(pmra_arr), pmra_arr, 0.0)
            pmdec_arr = np.where(np.isfinite(pmdec_arr), pmdec_arr, 0.0)
            mas_per_yr_to_deg = 1.0 / (3600.0 * 1000.0)
            cos_dec = np.cos(np.deg2rad(gaia_dec))
            cos_dec = np.where(np.abs(cos_dec) > 1e-6, cos_dec, 1e-6)
            d_ra = (pmra_arr * dt_yr * mas_per_yr_to_deg) / cos_dec
            d_dec = pmdec_arr * dt_yr * mas_per_yr_to_deg
            gaia_ra = gaia_ra + d_ra
            gaia_dec = gaia_dec + d_dec
            _log(
                f"[WCS] PM applied: dt={dt_yr:+.3f} yr  "
                f"|d_ra|_med={np.nanmedian(np.abs(d_ra * 3600.0 * 1000.0)):.1f} mas  "
                f"|d_dec|_med={np.nanmedian(np.abs(d_dec * 3600.0 * 1000.0)):.1f} mas"
            )

    ny, nx = img_shape

    # ── Seed WCS and project Gaia ──────────────────────────────────────────
    seed = _seed_wcs(approx_ra, approx_dec, approx_scale_arcsec, img_shape)
    cat_xy_full = _project_to_pixels(gaia_ra, gaia_dec, seed)

    diag_px = float(np.hypot(nx, ny))
    if bool(local_blind):
        cx, cy = nx / 2.0, ny / 2.0
        radius_px = 0.5 * diag_px * max(float(local_blind_radius_factor), 1.0)
        dist_c = np.hypot(cat_xy_full[:, 0] - cx, cat_xy_full[:, 1] - cy)
        in_field = np.isfinite(cat_xy_full).all(axis=1) & (dist_c <= radius_px)
        _log(
            f"[WCS] local-blind catalog window: radius={radius_px:.0f}px "
            f"factor={float(local_blind_radius_factor):.2f}"
        )
    else:
        # Small margin only — the brightness ranking of "in-field" catalogue
        # stars is what drives the quad-builder, so we must NOT pull in stars
        # far outside the actual FOV. Otherwise the brightest-N selection is
        # dominated by neighbours of the field that have no corresponding
        # detected source.
        field_pad = 0.10
        in_field = (
            np.isfinite(cat_xy_full).all(axis=1)
            & (cat_xy_full[:, 0] > -nx * field_pad)
            & (cat_xy_full[:, 0] < nx * (1.0 + field_pad))
            & (cat_xy_full[:, 1] > -ny * field_pad)
            & (cat_xy_full[:, 1] < ny * (1.0 + field_pad))
        )
    if not in_field.any():
        if bool(local_blind):
            _log("[WCS] No Gaia stars inside local-blind catalog window")
        else:
            _log("[WCS] No Gaia stars project into image field")
        return FAILED

    in_field_idx = np.flatnonzero(in_field)
    if cat_flux_arr is not None:
        flux_in = np.where(
            np.isfinite(cat_flux_arr[in_field_idx]),
            cat_flux_arr[in_field_idx],
            np.inf,
        )
        bright_order = in_field_idx[np.argsort(flux_in)]
    else:
        bright_order = in_field_idx

    cat_keep = bright_order[: int(n_brightest_cat)]
    cat_xy = cat_xy_full[cat_keep]
    cat_ra_kept = gaia_ra[cat_keep]
    cat_dec_kept = gaia_dec[cat_keep]

    # ── Pick brightest sources for quad construction ───────────────────────
    if src_flux_arr is not None:
        flux_score = np.where(np.isfinite(src_flux_arr), src_flux_arr, -np.inf)
        src_order = np.argsort(flux_score)[::-1]
    else:
        cx, cy = nx / 2.0, ny / 2.0
        d_c = np.hypot(sources_xy[:, 0] - cx, sources_xy[:, 1] - cy)
        src_order = np.argsort(d_c)
    src_xy = sources_xy[src_order[: int(n_brightest_src)]]

    _log(
        f"[WCS] quad-input: src={len(src_xy)} cat={len(cat_xy)} "
        f"(from {len(gaia_ra)} Gaia, {in_field.sum()} in field)"
    )

    if len(src_xy) < min_matches or len(cat_xy) < min_matches:
        _log("[WCS] Not enough stars for quad construction")
        return FAILED

    # ── Build quads on both sides ──────────────────────────────────────────
    scale_min_px = max(20.0, 0.03 * diag_px)
    scale_max_px = max(scale_min_px * 1.5, 0.85 * diag_px)

    src_quads = build_quads(
        src_xy,
        k_neighbor=int(quad_k_neighbor),
        neighbor_pool_factor=int(quad_neighbor_pool_factor),
        max_quads=int(quad_max_per_side),
        scale_min_px=scale_min_px,
        scale_max_px=scale_max_px,
    )
    cat_quads = build_quads(
        cat_xy,
        k_neighbor=int(quad_k_neighbor),
        neighbor_pool_factor=int(quad_neighbor_pool_factor),
        max_quads=int(quad_max_per_side),
        scale_min_px=scale_min_px,
        scale_max_px=scale_max_px,
    )
    _log(
        f"[WCS] quads built: src={len(src_quads)} cat={len(cat_quads)}  "
        f"(side range {scale_min_px:.0f}–{scale_max_px:.0f} px)"
    )
    if not src_quads or not cat_quads:
        _log("[WCS] Quad construction produced no usable quads")
        return FAILED

    # ── Match quad codes ───────────────────────────────────────────────────
    candidates = match_quads(
        src_quads,
        cat_quads,
        code_tol=float(quad_code_tol),
        scale_ratio_tol=float(quad_scale_ratio_tol),
        allow_reflection=bool(allow_reflection),
    )
    _log(f"[WCS] quad code matches: {len(candidates)}")
    if not candidates:
        _log("[WCS] No quad code matches — try widening quad_code_tol")
        return FAILED

    # ── RANSAC verify ──────────────────────────────────────────────────────
    ransac_candidates = ransac_verify_topk(
        src_xy,
        cat_xy,
        candidates,
        inlier_radius_px=float(ransac_inlier_radius_px),
        min_inliers=int(min_matches),
        max_trials=int(ransac_max_trials),
        try_reflection=bool(allow_reflection),
        n_keep=int(ransac_keep_candidates),
        log_fn=_log,
    )
    if not ransac_candidates:
        _log(
            f"[WCS] RANSAC failed (no candidate reached "
            f"min_matches={min_matches})"
        )
        return FAILED

    # ── Final fit all plausible candidates and choose by residual ──────────
    cat_radec = np.column_stack([cat_ra_kept, cat_dec_kept])
    proj_pt = SkyCoord(approx_ra * u.deg, approx_dec * u.deg)
    best_solution: dict | None = None

    for cand_idx, ransac in enumerate(ransac_candidates, start=1):
        _log(
            f"[WCS] Candidate {cand_idx}/{len(ransac_candidates)}: "
            f"ransac_inliers={ransac.n_inliers} rms_px={ransac.rms_px:.2f}"
        )
        src_in = src_xy[ransac.src_idx]
        ra_in = cat_ra_kept[ransac.cat_idx]
        dec_in = cat_dec_kept[ransac.cat_idx]
        if len(src_in) < min_matches:
            continue
        fit_proj_pt = proj_pt
        if bool(local_blind):
            fit_proj_pt = SkyCoord(
                float(np.nanmedian(ra_in)) * u.deg,
                float(np.nanmedian(dec_in)) * u.deg,
            )

        try:
            sky_in = SkyCoord(ra_in * u.deg, dec_in * u.deg)
            wcs = _fit_tan_wcs(src_in, sky_in, fit_proj_pt)
        except Exception as exc:
            _log(f"[WCS] Candidate {cand_idx}: initial fit failed: {exc}")
            continue

        # Sigma-clip refine over the full kept catalogue. Keep this TAN-only;
        # SIP is accepted below only after holdout validation.
        wcs, src_matched, cat_matched_radec = _sigma_clip_refine(
            sources_xy, cat_radec, wcs,
            sigma=refine_sigma, n_iter=4, min_pairs=min_matches,
        )
        if len(src_matched) < min_matches:
            src_matched = src_in
            cat_matched_radec = np.column_stack([ra_in, dec_in])
        if len(src_matched) < min_matches:
            _log(f"[WCS] Candidate {cand_idx}: too few refined pairs")
            continue

        sky_m = SkyCoord(
            cat_matched_radec[:, 0] * u.deg,
            cat_matched_radec[:, 1] * u.deg,
        )
        final_proj_pt = fit_proj_pt
        if bool(local_blind) and len(cat_matched_radec) > 0:
            final_proj_pt = SkyCoord(
                float(np.nanmedian(cat_matched_radec[:, 0])) * u.deg,
                float(np.nanmedian(cat_matched_radec[:, 1])) * u.deg,
            )
        try:
            tan_wcs = _fit_tan_wcs(src_matched, sky_m, final_proj_pt)
        except Exception as exc:
            _log(f"[WCS] Candidate {cand_idx}: TAN final fit failed: {exc}")
            continue

        final_wcs = tan_wcs
        final_model = "TAN"
        use_sip_final = int(sip_degree) > 0 and len(src_matched) >= int(sip_min_pairs)
        if use_sip_final:
            train_mask = _spatial_holdout_mask(
                src_matched,
                fraction=float(sip_holdout_fraction),
                min_holdout=int(sip_min_holdout_pairs),
                min_train=int(sip_min_pairs),
            )
            if train_mask is None:
                _log(
                    f"[WCS] Candidate {cand_idx}: SIP skipped "
                    f"(need train>={int(sip_min_pairs)} and "
                    f"holdout>={int(sip_min_holdout_pairs)})"
                )
            else:
                hold_mask = ~train_mask
                try:
                    sky_train = SkyCoord(
                        cat_matched_radec[train_mask, 0] * u.deg,
                        cat_matched_radec[train_mask, 1] * u.deg,
                    )
                    sky_hold = SkyCoord(
                        cat_matched_radec[hold_mask, 0] * u.deg,
                        cat_matched_radec[hold_mask, 1] * u.deg,
                    )
                    tan_train = _fit_tan_wcs(
                        src_matched[train_mask], sky_train, final_proj_pt,
                    )
                    sip_train = _fit_sip_wcs(
                        src_matched[train_mask], sky_train,
                        final_proj_pt, int(sip_degree),
                    )
                    tan_hold_rms = _pair_rms_arcsec(
                        tan_train, src_matched[hold_mask], sky_hold,
                    )
                    sip_hold_rms = _pair_rms_arcsec(
                        sip_train, src_matched[hold_mask], sky_hold,
                    )
                    improve = (
                        (tan_hold_rms - sip_hold_rms) / tan_hold_rms
                        if np.isfinite(tan_hold_rms) and tan_hold_rms > 0
                        and np.isfinite(sip_hold_rms)
                        else float("nan")
                    )
                    if np.isfinite(improve) and improve >= float(sip_min_improvement):
                        final_wcs = _fit_sip_wcs(
                            src_matched, sky_m, final_proj_pt, int(sip_degree),
                        )
                        final_model = f"SIP{int(sip_degree)}"
                        _log(
                            f"[WCS] Candidate {cand_idx}: SIP adopted "
                            f"degree={int(sip_degree)} holdout "
                            f"{tan_hold_rms:.3f}\"→{sip_hold_rms:.3f}\" "
                            f"improve={improve * 100.0:.1f}%"
                        )
                    else:
                        _log(
                            f"[WCS] Candidate {cand_idx}: SIP rejected "
                            f"holdout {tan_hold_rms:.3f}\"→{sip_hold_rms:.3f}\" "
                            f"improve={improve * 100.0 if np.isfinite(improve) else np.nan:.1f}%"
                        )
                except Exception as exc:
                    _log(
                        f"[WCS] Candidate {cand_idx}: SIP validation failed "
                        f"({exc}); keeping TAN"
                    )

        rms = _pair_rms_arcsec(final_wcs, src_matched, sky_m)

        cd = _wcs_cd_matrix(final_wcs)
        if cd is not None and np.isfinite(cd).all():
            rotation_deg = float(np.degrees(np.arctan2(-cd[0, 1], cd[1, 1])))
            det = float(cd[0, 0] * cd[1, 1] - cd[0, 1] * cd[1, 0])
            scale_arcsec_per_px = float(np.sqrt(abs(det)) * 3600.0)
            flip_x = bool(det > 0)  # CDELT[0]<0 convention -> det<0 normally
            flip_y = False
        else:
            rotation_deg = float("nan")
            scale_arcsec_per_px = float("nan")
            flip_x = False
            flip_y = False

        n_verified = 0
        verification_ok = True
        try:
            proj_all = _project_to_pixels(gaia_ra, gaia_dec, final_wcs)
            in_frame = (
                np.isfinite(proj_all).all(axis=1)
                & (proj_all[:, 0] > 0) & (proj_all[:, 0] < nx)
                & (proj_all[:, 1] > 0) & (proj_all[:, 1] < ny)
            )
            if in_frame.sum() >= min_matches:
                tight_r = max(float(ransac_inlier_radius_px), 2.5)
                v_src, _, _ = _match_unique(sources_xy, proj_all, tight_r)
                n_verified = int(len(v_src))
                if n_verified < min_matches:
                    verification_ok = False
        except Exception as exc:
            _log(f"[WCS] Candidate {cand_idx}: verification skipped ({exc})")

        _log(
            f"[WCS] Candidate {cand_idx}: RMS={rms:.3f}\" "
            f"matches={len(src_matched)} verified={n_verified}"
        )
        if not verification_ok:
            _log(f"[WCS] Candidate {cand_idx}: verification failed")
            continue

        score_rms = rms if np.isfinite(rms) else float("inf")
        score = (score_rms, -int(n_verified), -int(len(src_matched)), ransac.rms_px)
        candidate = {
            "score": score,
            "candidate_idx": cand_idx,
            "wcs": final_wcs,
            "model": final_model,
            "n_matches": int(len(src_matched)),
            "rms_arcsec": rms,
            "rotation_deg": rotation_deg,
            "scale_arcsec_per_px": scale_arcsec_per_px,
            "flip_x": flip_x,
            "flip_y": flip_y,
            "n_verified": int(n_verified),
        }
        if best_solution is None or score < best_solution["score"]:
            best_solution = candidate

    if best_solution is None:
        _log("[WCS] All RANSAC candidates failed final fit/verification")
        return FAILED

    final_wcs = best_solution["wcs"]
    rms = float(best_solution["rms_arcsec"])
    rotation_deg = float(best_solution["rotation_deg"])
    scale_arcsec_per_px = float(best_solution["scale_arcsec_per_px"])
    flip_x = bool(best_solution["flip_x"])
    flip_y = bool(best_solution["flip_y"])
    n_matches_final = int(best_solution["n_matches"])
    _log(
        f"[WCS] Selected candidate {best_solution['candidate_idx']}: "
        f"model={best_solution.get('model', 'TAN')}  "
        f"RMS={rms:.3f}\"  rot={rotation_deg:.2f}°  "
        f"scale={scale_arcsec_per_px:.4f}\"/px  "
        f"matches={n_matches_final} verified={best_solution['n_verified']}"
    )

    return WCSSolveResult(
        converged=True,
        wcs=final_wcs,
        n_matches=n_matches_final,
        rms_arcsec=rms,
        rotation_deg=rotation_deg,
        scale_arcsec_per_px=scale_arcsec_per_px,
        flip_x=flip_x,
        flip_y=flip_y,
        model=str(best_solution.get("model", "TAN")),
        sip_order=int(sip_degree) if str(best_solution.get("model", "")).startswith("SIP") else 0,
        log=log,
    )
