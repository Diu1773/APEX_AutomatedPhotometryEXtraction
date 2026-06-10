"""Astrometry.net-style geometric quad matching.

Replaces the brute-force (rotation × scale × translation) grid search used by
the previous internal solver. That score (nearest-neighbour count) plateaus on
crowded fields like M13 because both detected and Gaia distributions are dense,
so any plausible orientation accumulates a similar count and the global maximum
is not the true alignment.

A quad code is a 4-D shape descriptor that is invariant to translation,
rotation, and uniform scale. Two points A, B of the quad define the local
coordinate frame and the remaining two C, D contribute four numbers. In a
crowded field the chance that two unrelated 4-point configurations share the
same 4-D code is vanishingly small, so kd-tree matching in code space gives a
small candidate set that RANSAC can verify by counting inliers under the
implied similarity transform.

Reference: Lang et al. 2010, "Astrometry.net: Blind astrometric calibration of
arbitrary astronomical images", AJ 139, 1782.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from typing import Optional

import numpy as np
from scipy.spatial import cKDTree


# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

@dataclass
class Quad:
    """One geometric quad in a point set.

    indices : (A, B, C, D) — indices into the original xy array. A and B are
        the pair with maximum separation; C and D are the other two, ordered
        so that Cx ≤ Dx in the A→B canonical frame.
    code    : 4-vector (Cx, Cy, Dx, Dy) — the shape descriptor used for hash
        matching. Translation/rotation/scale invariant by construction.
    side_len: |AB| in input units (pixels). Larger quads constrain transforms
        more tightly, so we sort by this when picking which to try first.
    """

    indices: tuple[int, int, int, int]
    code: np.ndarray
    side_len: float


@dataclass
class RansacResult:
    transform: Optional[np.ndarray]   # 2x3 affine: cat_xy → src_xy
    src_idx: np.ndarray               # matched source indices
    cat_idx: np.ndarray               # matched catalog indices
    n_inliers: int
    rms_px: float
    log: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Quad construction
# ---------------------------------------------------------------------------

def _canonical_quad_code(p: np.ndarray):
    """Build the canonical quad code from 4 points.

    Frame convention: A → (0, 0), B → (1, 0).
      - A, B = pair with maximum separation in p.
      - In the AB frame, the remaining two points C, D have coords
        (Cx, Cy, Dx, Dy). Canonicalise C, D order so Cx ≤ Dx.
      - Resolve the A↔B ambiguity by requiring Cx + Dx ≤ 1 (geometric mid of
        C, D lies on the A-side of the AB midpoint). Swap A and B if not, and
        re-project (under the swap, x → 1 − x, y → −y).

    Returns ``(a_idx_local, b_idx_local, c_idx_local, d_idx_local, code,
    side_len)`` with the four indices referring to rows of ``p`` (0..3), or
    ``None`` if the quad is degenerate.
    """
    # Find the max-distance pair.
    best_d2 = -1.0
    a_loc = b_loc = -1
    for i, j in combinations(range(4), 2):
        dx = p[i, 0] - p[j, 0]
        dy = p[i, 1] - p[j, 1]
        d2 = dx * dx + dy * dy
        if d2 > best_d2:
            best_d2 = d2
            a_loc, b_loc = i, j
    if best_d2 <= 0:
        return None

    cd_loc = [k for k in range(4) if k != a_loc and k != b_loc]

    # Frame basis: u along AB (length |AB|² so that B → (1, 0));
    #              v perpendicular (rotated 90° CCW from AB).
    ab = p[b_loc] - p[a_loc]
    inv_l2 = 1.0 / best_d2
    ux, uy = ab[0] * inv_l2, ab[1] * inv_l2
    vx, vy = -ab[1] * inv_l2, ab[0] * inv_l2

    c_vec = p[cd_loc[0]] - p[a_loc]
    d_vec = p[cd_loc[1]] - p[a_loc]
    cx = ux * c_vec[0] + uy * c_vec[1]
    cy = vx * c_vec[0] + vy * c_vec[1]
    dx_ = ux * d_vec[0] + uy * d_vec[1]
    dy_ = vx * d_vec[0] + vy * d_vec[1]

    # Resolve A↔B ambiguity (canonical form: Cx + Dx ≤ 1).
    if cx + dx_ > 1.0:
        cx, dx_ = 1.0 - cx, 1.0 - dx_
        cy, dy_ = -cy, -dy_
        a_loc, b_loc = b_loc, a_loc

    # Order C, D so that Cx ≤ Dx.
    if cx > dx_:
        cx, dx_ = dx_, cx
        cy, dy_ = dy_, cy
        c_orig, d_orig = cd_loc[1], cd_loc[0]
    else:
        c_orig, d_orig = cd_loc[0], cd_loc[1]

    code = np.array([cx, cy, dx_, dy_], dtype=float)
    side_len = float(np.sqrt(best_d2))
    return a_loc, b_loc, c_orig, d_orig, code, side_len


# Combination lookup for the batched canonical-code path. Index k selects the
# k-th of the 6 unordered pairs of 4 points; (_AB_I, _AB_J) is that pair and
# (_CD0, _CD1) are the two remaining points in ascending index order — matching
# the scalar ``[k for k in range(4) if k not in (a, b)]`` ordering exactly.
_AB_I = np.array([0, 0, 0, 1, 1, 2])
_AB_J = np.array([1, 2, 3, 2, 3, 3])
_CD0  = np.array([2, 1, 1, 0, 0, 0])
_CD1  = np.array([3, 3, 2, 3, 2, 1])


def _canonical_quad_codes_batch(pts: np.ndarray):
    """Vectorised equivalent of ``_canonical_quad_code`` over many quads.

    ``pts`` has shape ``(M, 4, 2)``. Returns
    ``(a_loc, b_loc, c_loc, d_loc, code, side_len, valid)`` where the first
    four are ``(M,)`` int local indices (0..3), ``code`` is ``(M, 4)``,
    ``side_len`` is ``(M,)`` and ``valid`` is an ``(M,)`` bool mask (degenerate
    quads with a zero max separation are ``False``). The algorithm mirrors the
    scalar version exactly so results are bit-for-bit comparable; it exists
    purely to replace the per-quad Python loop that dominated the profile.
    """
    pts = np.asarray(pts, dtype=float)
    m = pts.shape[0]
    if m == 0:
        e = np.zeros(0, dtype=int)
        return (e, e, e, e, np.zeros((0, 4)), np.zeros(0), np.zeros(0, bool))

    row = np.arange(m)
    # Squared distance for each of the 6 pairs → pick the max-separation pair.
    diff = pts[:, _AB_I, :] - pts[:, _AB_J, :]            # (M, 6, 2)
    d2 = np.einsum("mkc,mkc->mk", diff, diff)             # (M, 6)
    arg = np.argmax(d2, axis=1)                           # (M,)
    best_d2 = d2[row, arg]

    a_loc = _AB_I[arg]
    b_loc = _AB_J[arg]
    cd0 = _CD0[arg]
    cd1 = _CD1[arg]

    pA = pts[row, a_loc]
    pB = pts[row, b_loc]
    pC = pts[row, cd0]
    pD = pts[row, cd1]

    valid = best_d2 > 0.0
    inv_l2 = np.where(valid, 1.0 / np.where(valid, best_d2, 1.0), 0.0)
    ab = pB - pA
    ux = ab[:, 0] * inv_l2
    uy = ab[:, 1] * inv_l2
    vx = -ab[:, 1] * inv_l2
    vy = ab[:, 0] * inv_l2

    cvec = pC - pA
    dvec = pD - pA
    cx = ux * cvec[:, 0] + uy * cvec[:, 1]
    cy = vx * cvec[:, 0] + vy * cvec[:, 1]
    dx_ = ux * dvec[:, 0] + uy * dvec[:, 1]
    dy_ = vx * dvec[:, 0] + vy * dvec[:, 1]

    # Resolve A↔B ambiguity (canonical form Cx + Dx ≤ 1): x→1−x, y→−y, swap A,B.
    swap_ab = (cx + dx_) > 1.0
    cx = np.where(swap_ab, 1.0 - cx, cx)
    dx_ = np.where(swap_ab, 1.0 - dx_, dx_)
    cy = np.where(swap_ab, -cy, cy)
    dy_ = np.where(swap_ab, -dy_, dy_)
    a_out = np.where(swap_ab, b_loc, a_loc)
    b_out = np.where(swap_ab, a_loc, b_loc)

    # Order C, D so Cx ≤ Dx.
    swap_cd = cx > dx_
    cx2 = np.where(swap_cd, dx_, cx)
    dx2 = np.where(swap_cd, cx, dx_)
    cy2 = np.where(swap_cd, dy_, cy)
    dy2 = np.where(swap_cd, cy, dy_)
    c_out = np.where(swap_cd, cd1, cd0)
    d_out = np.where(swap_cd, cd0, cd1)

    code = np.column_stack([cx2, cy2, dx2, dy2])
    side_len = np.sqrt(np.where(valid, best_d2, 0.0))
    return (
        a_out.astype(int), b_out.astype(int), c_out.astype(int), d_out.astype(int),
        code, side_len, valid,
    )


def build_quads(
    xy: np.ndarray,
    *,
    k_neighbor: int = 7,
    neighbor_pool_factor: int = 3,
    max_quads: int = 500,
    scale_min_px: float = 20.0,
    scale_max_px: float = 5000.0,
    require_compact: bool = True,
) -> list[Quad]:
    """Build quad features from a point set.

    For each star, sample ``k_neighbor`` neighbours from a wider nearest-neighbour
    pool and form all C(k, 3) quads consisting of the star itself plus three of
    those neighbours. Sampling from a pool rather than taking only the closest
    neighbours avoids crowded-field core patterns where every quad is tiny and
    repetitive.

    Filters:
      * |AB| in [scale_min_px, scale_max_px] — discards quads that span the
        whole image (rarely useful, scale-degenerate) and tiny quads (noisy).
      * If ``require_compact`` is true, both C and D must satisfy
        0 ≤ Cx, Dx ≤ 1 — i.e. they lie inside the AB strip. This is the
        astnet "all four points inside the smallest circle through the
        extremal pair" condition in code-space form. Rejects elongated quads
        whose code is dominated by C/D extrapolating beyond AB.

    Sort by side length descending so the most discriminating quads come
    first. Cap at ``max_quads``.
    """
    xy = np.asarray(xy, dtype=float)
    n = len(xy)
    if n < 4:
        return []

    k = min(int(k_neighbor), n - 1)
    if k < 3:
        return []
    pool_factor = max(1, int(neighbor_pool_factor))
    k_query = min(max(k, k * pool_factor), n - 1)

    tree = cKDTree(xy)
    # k+1 because the query point itself is returned as the first neighbour.
    _, nbrs = tree.query(xy, k=k_query + 1)

    # Phase 1: gather unique 4-point index sets (cheap Python/bookkeeping).
    # Phase 2: compute every canonical code in one batched numpy call. The
    # per-quad scalar code computation used to dominate the profile, so it is
    # hoisted out of the loop.
    seen: set[tuple[int, int, int, int]] = set()
    idx4_list: list[tuple[int, int, int, int]] = []

    for i in range(n):
        pool = [int(j) for j in nbrs[i] if int(j) != i]
        if len(pool) > k:
            # Evenly sample rank-ordered neighbours: keep local constraints while
            # also injecting longer baselines that survive dense globular cores.
            pick = np.linspace(0, len(pool) - 1, k)
            pick = np.unique(np.rint(pick).astype(int))
            neigh = [pool[int(j)] for j in pick[:k]]
            if len(neigh) < k:
                for j in pool:
                    if j not in neigh:
                        neigh.append(j)
                    if len(neigh) >= k:
                        break
        else:
            neigh = pool[:k]
        if len(neigh) < 3:
            continue
        for trio in combinations(neigh, 3):
            idx4 = tuple(sorted((i, *trio)))
            if idx4 in seen:
                continue
            seen.add(idx4)
            idx4_list.append(idx4)

    if not idx4_list:
        return []

    quads_idx = np.asarray(idx4_list, dtype=int)            # (M, 4)
    pts_batch = xy[quads_idx]                               # (M, 4, 2)
    a_loc, b_loc, c_loc, d_loc, codes, sides, valid = _canonical_quad_codes_batch(pts_batch)

    keep = valid & (sides >= scale_min_px) & (sides <= scale_max_px)
    if require_compact:
        cx = codes[:, 0]; cy_ = codes[:, 1]; dx_ = codes[:, 2]; dy_ = codes[:, 3]
        keep &= (cx >= 0.0) & (cx <= 1.0) & (dx_ >= 0.0) & (dx_ <= 1.0)
        keep &= (np.abs(cy_) <= 1.0) & (np.abs(dy_) <= 1.0)

    kept = np.flatnonzero(keep)
    if kept.size == 0:
        return []

    rows = np.arange(quads_idx.shape[0])
    g_a = quads_idx[rows, a_loc]
    g_b = quads_idx[rows, b_loc]
    g_c = quads_idx[rows, c_loc]
    g_d = quads_idx[rows, d_loc]

    quads: list[Quad] = [
        Quad(
            indices=(int(g_a[r]), int(g_b[r]), int(g_c[r]), int(g_d[r])),
            code=codes[r],
            side_len=float(sides[r]),
        )
        for r in kept
    ]

    quads.sort(key=lambda q: q.side_len, reverse=True)
    return quads[:max_quads]


# ---------------------------------------------------------------------------
# Quad matching
# ---------------------------------------------------------------------------

def match_quads(
    src_quads: list[Quad],
    cat_quads: list[Quad],
    *,
    code_tol: float = 0.01,
    scale_ratio_tol: float = 0.15,
    k_nearest: int = 5,
    allow_reflection: bool = True,
) -> list[tuple[Quad, Quad]]:
    """Match source quads to catalogue quads via 4-D code kd-tree lookup.

    For each src quad, the ``k_nearest`` closest catalogue quads in code space
    are considered. A pair survives if:
      * 4-D code L2 distance < ``code_tol``
      * |log(src_side / cat_side)| < ``scale_ratio_tol``

    If ``allow_reflection`` is True, also search for the mirror code
    (Cx, −Cy, Dx, −Dy) — this covers image / sky handedness flips without
    needing to build a second set of source quads.
    """
    if not src_quads or not cat_quads:
        return []

    cat_codes = np.array([q.code for q in cat_quads], dtype=float)
    cat_scales = np.array([q.side_len for q in cat_quads], dtype=float)
    tree = cKDTree(cat_codes)
    k_eff = max(1, min(int(k_nearest), len(cat_quads)))

    src_codes = np.array([q.code for q in src_quads], dtype=float)     # (Ns, 4)
    src_sides = np.array([q.side_len for q in src_quads], dtype=float) # (Ns,)
    log_src = np.where(src_sides > 0, np.log(np.where(src_sides > 0, src_sides, 1.0)), np.nan)
    log_cat = np.where(cat_scales > 0, np.log(np.where(cat_scales > 0, cat_scales, 1.0)), np.nan)

    signs = [np.array([1.0, 1.0, 1.0, 1.0])]
    if allow_reflection:
        signs.append(np.array([1.0, -1.0, 1.0, -1.0]))

    pairs: list[tuple[Quad, Quad]] = []
    for sign in signs:
        # One batched kd-tree query for ALL source quads (cKDTree.query is
        # vectorised) instead of a Python per-quad loop — this was the main
        # remaining hotspot after build_quads was vectorised.
        dists, idxs = tree.query(src_codes * sign, k=k_eff)
        # cKDTree.query returns 1-D arrays when k_eff == 1; normalise to (Ns, k).
        if dists.ndim == 1:
            dists = dists.reshape(-1, 1)
            idxs = idxs.reshape(-1, 1)
        # Scale-ratio test vectorised across (Ns, k_eff).
        scale_ok = np.abs(log_src[:, None] - log_cat[idxs]) <= scale_ratio_tol
        keep = np.isfinite(dists) & (dists <= code_tol) & np.isfinite(log_src)[:, None] & scale_ok
        si_arr, kk_arr = np.nonzero(keep)
        for si, kk in zip(si_arr.tolist(), kk_arr.tolist()):
            ci = int(idxs[si, kk])
            pairs.append((src_quads[si], cat_quads[ci]))
    return pairs


# ---------------------------------------------------------------------------
# RANSAC verification
# ---------------------------------------------------------------------------

def _similarity_from_pair(
    src_xy: np.ndarray, cat_xy: np.ndarray, sq: Quad, cq: Quad,
) -> Optional[np.ndarray]:
    """2-D similarity transform (rotation + uniform scale + translation,
    possibly reflected) mapping cat → src using the A, B endpoints of the
    quad pair. Returns a 2×3 matrix or None if degenerate.

    For two corresponding pairs (cA, cB) → (sA, sB) the similarity is unique
    up to reflection; we try both orientations (det > 0 and det < 0) and let
    the inlier count decide.
    """
    sA = src_xy[sq.indices[0]]
    sB = src_xy[sq.indices[1]]
    cA = cat_xy[cq.indices[0]]
    cB = cat_xy[cq.indices[1]]

    v_src = sB - sA
    v_cat = cB - cA
    len_src2 = v_src[0] ** 2 + v_src[1] ** 2
    len_cat2 = v_cat[0] ** 2 + v_cat[1] ** 2
    if len_src2 <= 0 or len_cat2 <= 0:
        return None

    inv_cat2 = 1.0 / len_cat2
    cos_t = (v_src[0] * v_cat[0] + v_src[1] * v_cat[1]) * inv_cat2
    sin_t = (v_src[1] * v_cat[0] - v_src[0] * v_cat[1]) * inv_cat2
    a = cos_t
    b = -sin_t
    c = sin_t
    d = cos_t
    tx = sA[0] - (a * cA[0] + b * cA[1])
    ty = sA[1] - (c * cA[0] + d * cA[1])
    return np.array([[a, b, tx], [c, d, ty]], dtype=float)


def _similarity_from_pair_reflected(
    src_xy: np.ndarray, cat_xy: np.ndarray, sq: Quad, cq: Quad,
) -> Optional[np.ndarray]:
    """Same as ``_similarity_from_pair`` but with a reflection (det < 0).

    We achieve this by negating the y-component of the catalogue vector before
    solving — equivalent to mirroring across the x axis in catalogue space.
    """
    sA = src_xy[sq.indices[0]]
    sB = src_xy[sq.indices[1]]
    cA = cat_xy[cq.indices[0]]
    cB = cat_xy[cq.indices[1]]

    v_src = sB - sA
    v_cat = cB - cA
    v_cat_m = np.array([v_cat[0], -v_cat[1]], dtype=float)
    len_src2 = v_src[0] ** 2 + v_src[1] ** 2
    len_cat2 = v_cat_m[0] ** 2 + v_cat_m[1] ** 2
    if len_src2 <= 0 or len_cat2 <= 0:
        return None

    inv_cat2 = 1.0 / len_cat2
    cos_t = (v_src[0] * v_cat_m[0] + v_src[1] * v_cat_m[1]) * inv_cat2
    sin_t = (v_src[1] * v_cat_m[0] - v_src[0] * v_cat_m[1]) * inv_cat2
    # M_reflected = M_rot · diag(1, -1)
    a = cos_t
    b = sin_t
    c = sin_t
    d = -cos_t
    tx = sA[0] - (a * cA[0] + b * cA[1])
    ty = sA[1] - (c * cA[0] + d * cA[1])
    return np.array([[a, b, tx], [c, d, ty]], dtype=float)


def _apply_affine(xy: np.ndarray, M: np.ndarray) -> np.ndarray:
    return (M[:, :2] @ xy.T).T + M[:, 2]


def _count_unique_inliers(
    src_xy: np.ndarray,
    cat_t: np.ndarray,
    src_tree: cKDTree,
    inlier_radius_px: float,
):
    """Greedy 1:1 nearest-neighbour matching of cat→src within radius.

    Returns ``(src_idx, cat_idx, dists)`` arrays. Same convention as the
    helper from the previous solver, but without the duplicate-resolution
    overhead of the larger ``_match_unique`` because we are inside the RANSAC
    inner loop.
    """
    dists, idxs = src_tree.query(cat_t, k=1)
    finite = np.isfinite(dists)
    keep = finite & (dists <= inlier_radius_px)
    if not keep.any():
        empty = np.array([], dtype=int)
        return empty, empty, np.array([], dtype=float)

    cand_src = idxs[keep].astype(int)
    cand_cat = np.flatnonzero(keep)
    cand_dist = dists[keep].astype(float)
    order = np.argsort(cand_dist)

    used_src: set[int] = set()
    keep_s: list[int] = []
    keep_c: list[int] = []
    keep_d: list[float] = []
    for pos in order:
        s = int(cand_src[pos])
        if s in used_src:
            continue
        used_src.add(s)
        keep_s.append(s)
        keep_c.append(int(cand_cat[pos]))
        keep_d.append(float(cand_dist[pos]))
    return (
        np.asarray(keep_s, dtype=int),
        np.asarray(keep_c, dtype=int),
        np.asarray(keep_d, dtype=float),
    )


def ransac_verify(
    src_xy: np.ndarray,
    cat_xy: np.ndarray,
    candidates: list[tuple[Quad, Quad]],
    *,
    inlier_radius_px: float = 3.0,
    min_inliers: int = 12,
    max_trials: int = 500,
    try_reflection: bool = True,
    log_fn=None,
) -> RansacResult:
    """For each candidate quad pair compute the implied similarity, count
    inliers, and keep the best.

    Order the candidates by source-quad side length descending so the most
    constraining transforms are tried first. Stop after ``max_trials`` to keep
    runtime bounded even when ``match_quads`` returns thousands of pairs.
    """
    src_xy = np.asarray(src_xy, dtype=float)
    cat_xy = np.asarray(cat_xy, dtype=float)
    log: list[str] = []

    top = ransac_verify_topk(
        src_xy,
        cat_xy,
        candidates,
        inlier_radius_px=inlier_radius_px,
        min_inliers=min_inliers,
        max_trials=max_trials,
        try_reflection=try_reflection,
        n_keep=1,
        log_fn=log_fn,
    )
    if top:
        return top[0]

    src_xy = np.asarray(src_xy, dtype=float)
    cat_xy = np.asarray(cat_xy, dtype=float)
    log: list[str] = []
    return RansacResult(
        transform=None,
        src_idx=np.array([], dtype=int),
        cat_idx=np.array([], dtype=int),
        n_inliers=0,
        rms_px=np.inf,
        log=log,
    )


def _inlier_overlap(a: RansacResult, b: RansacResult) -> float:
    if len(a.src_idx) == 0 or len(b.src_idx) == 0:
        return 0.0
    pa = set(zip(a.src_idx.tolist(), a.cat_idx.tolist()))
    pb = set(zip(b.src_idx.tolist(), b.cat_idx.tolist()))
    if not pa or not pb:
        return 0.0
    return len(pa & pb) / max(1, min(len(pa), len(pb)))


def ransac_verify_topk(
    src_xy: np.ndarray,
    cat_xy: np.ndarray,
    candidates: list[tuple[Quad, Quad]],
    *,
    inlier_radius_px: float = 3.0,
    min_inliers: int = 12,
    max_trials: int = 500,
    try_reflection: bool = True,
    n_keep: int = 6,
    distinct_overlap: float = 0.85,
    log_fn=None,
) -> list[RansacResult]:
    """Return several plausible RANSAC transforms, not just the top count.

    Dense fields can produce a wrong transform with a high nearest-neighbour
    inlier count. The solver can now final-fit and score several distinct
    candidates by astrometric residual instead of committing to the first high
    count transform.
    """
    src_xy = np.asarray(src_xy, dtype=float)
    cat_xy = np.asarray(cat_xy, dtype=float)
    log: list[str] = []

    if not candidates:
        return []

    candidates = sorted(candidates, key=lambda p: p[0].side_len, reverse=True)
    src_tree = cKDTree(src_xy)
    viable: list[RansacResult] = []
    best_any = RansacResult(
        transform=None,
        src_idx=np.array([], dtype=int),
        cat_idx=np.array([], dtype=int),
        n_inliers=0,
        rms_px=np.inf,
        log=log,
    )

    n_tried = 0
    for sq, cq in candidates:
        if n_tried >= max_trials:
            break
        n_tried += 1

        for builder in (
            (_similarity_from_pair, False),
            (_similarity_from_pair_reflected, True),
        ):
            fn, is_refl = builder
            if is_refl and not try_reflection:
                continue
            M = fn(src_xy, cat_xy, sq, cq)
            if M is None:
                continue
            cat_t = _apply_affine(cat_xy, M)
            src_i, cat_i, dists = _count_unique_inliers(
                src_xy, cat_t, src_tree, inlier_radius_px,
            )
            n_in = int(len(src_i))
            rms = float(np.sqrt(np.mean(dists ** 2))) if n_in else np.inf
            result = RansacResult(
                transform=M,
                src_idx=src_i,
                cat_idx=cat_i,
                n_inliers=n_in,
                rms_px=rms,
                log=log,
            )
            if (
                n_in > best_any.n_inliers
                or (n_in == best_any.n_inliers and rms < best_any.rms_px)
            ):
                best_any = result
            if n_in >= min_inliers:
                viable.append(result)

    viable.sort(key=lambda r: (-r.n_inliers, r.rms_px))
    kept: list[RansacResult] = []
    for cand in viable:
        if len(kept) >= max(1, int(n_keep)):
            break
        if any(_inlier_overlap(cand, old) >= float(distinct_overlap) for old in kept):
            continue
        kept.append(cand)

    msg = (
        f"[quad] RANSAC tried={n_tried} best_inliers={best_any.n_inliers} "
        f"rms_px={best_any.rms_px:.2f} kept={len(kept)}"
    )
    if log_fn is not None:
        log_fn(msg)
    log.append(msg)

    if kept:
        return kept

    if best_any.n_inliers > 0:
        log.append(
            f"RANSAC failed best_inliers={best_any.n_inliers} < {min_inliers}"
        )
    return []
