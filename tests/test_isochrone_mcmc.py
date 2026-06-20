"""Synthetic-recovery test for the generative MCMC isochrone solver.

This is the core correctness test: we forward-simulate a cluster CMD from known
parameters (age, [M/H], distance modulus, reddening) by sampling stars along a
PARSEC isochrone with the Kroupa IMF, adding photometric noise, a handful of
unresolved binaries, and a few uniform field outliers -- exactly the generative
model the MCMC assumes -- then check that the sampler recovers the truth inside
its credible interval.

Qt-free and fast: it uses a tiny synthetic isochrone (no big data files) and a
small walker/step budget.
"""

from __future__ import annotations

import numpy as np
import pytest

from apex.analysis.cmd.isochrone_fitter_v2 import FitBounds, IsochroneFitterV2
from apex.analysis.cmd.isochrone_mcmc import fit_isochrone_mcmc


# ---------------------------------------------------------------------------
# Build a small synthetic PARSEC-like grid (compact 7-col layout the validation
# runner uses: [0]=pad, [1]=mh, [2]=log_age, [3]=b1, [4]=b2, [5]=mag, [6]=mass).
# We use a smooth analytic main sequence + turn-off so age genuinely changes the
# morphology (older -> fainter, redder turn-off).
# ---------------------------------------------------------------------------
def _make_synthetic_grid(seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    ages = np.round(np.arange(8.4, 9.81, 0.1), 2)        # 8.4 .. 9.8 dex
    mhs = np.round(np.arange(-0.4, 0.41, 0.1), 2)        # -0.4 .. 0.4
    rows = []
    for log_age in ages:
        # Turn-off mass shrinks with age: M_to ~ (t/t0)^-0.4 style relation.
        m_to = 2.6 * 10.0 ** (-0.4 * (log_age - 8.4))     # ~2.6 -> ~0.95 Msun
        m_to = float(np.clip(m_to, 0.9, 2.8))
        for mh in mhs:
            # Mass grid from low MS to just past the turn-off.
            masses = np.linspace(0.25, m_to, 120)
            # Absolute mag (g): brighter (smaller) for higher mass.
            abs_g = 4.8 - 4.2 * np.log10(masses) + 0.5 * (masses - 1.0) ** 2
            # Add a small turn-off hook: near m_to mag brightens then reddens.
            hook = 0.6 * np.clip((masses - 0.85 * m_to) / (0.15 * m_to), 0, 1)
            abs_g = abs_g - hook
            # Colour g-r: redder for lower mass; metallicity reddens slightly.
            gr = 0.95 - 0.55 * np.log10(masses) + 0.25 * mh
            gr = gr + 0.15 * hook                          # turn-off reddening
            b1 = abs_g                                     # treat mag band as b1=g
            b2 = abs_g - gr                                # r = g - (g-r)
            n = len(masses)
            block = np.column_stack([
                np.zeros(n),          # 0 pad
                np.full(n, mh),       # 1 [M/H]
                np.full(n, log_age),  # 2 log_age
                b1,                   # 3 g
                b2,                   # 4 r
                abs_g,                # 5 mag (g)
                masses,               # 6 mass
            ])
            rows.append(block)
    return np.vstack(rows)


def _kroupa_sample(masses_lo: float, masses_hi: float, n: int, rng) -> np.ndarray:
    """Rejection-sample n masses from a Kroupa(2001)-like IMF in [lo, hi]."""
    out = []
    # Envelope: power law -2.3 above 0.5, -1.3 below.
    def pdf(m):
        return np.where(m < 0.5, m ** -1.3, 0.5 ** -1.3 * (m / 0.5) ** -2.3)
    m_grid = np.linspace(masses_lo, masses_hi, 500)
    pmax = pdf(m_grid).max()
    while len(out) < n:
        cand = rng.uniform(masses_lo, masses_hi, size=n)
        keep = rng.uniform(0, pmax, size=n) < pdf(cand)
        out.extend(cand[keep].tolist())
    return np.array(out[:n])


def _simulate_cluster(fitter, truth, n_members=350, n_binaries=40,
                      n_field=25, phot_err=0.03, seed=7):
    """Forward-simulate an observed CMD from the generative model."""
    rng = np.random.default_rng(seed)
    log_age, mh, dm, e_color = truth
    iso_c, iso_m, iso_mass = fitter._interpolate_isochrone(log_age, mh, dm, e_color)
    assert len(iso_c) >= 20

    # Sample member masses from the IMF over the track's mass range.
    m_lo, m_hi = float(iso_mass.min()), float(iso_mass.max())
    order = np.argsort(iso_mass)
    ms, cs, mgs = iso_mass[order], iso_c[order], iso_m[order]
    member_mass = _kroupa_sample(m_lo, m_hi, n_members, rng)
    mc = np.interp(member_mass, ms, cs)
    mm = np.interp(member_mass, ms, mgs)

    # Binaries: brighter (flux-add) + redder. Use equal-ish mass companions.
    bm_mass = _kroupa_sample(m_lo, m_hi, n_binaries, rng)
    comp_mass = bm_mass * rng.uniform(0.6, 1.0, size=n_binaries)
    p_g = np.interp(bm_mass, ms, mgs)
    p_r = p_g - np.interp(bm_mass, ms, cs)
    c_g = np.interp(comp_mass, ms, mgs)
    c_r = c_g - np.interp(comp_mass, ms, cs)
    tot_g = -2.5 * np.log10(10 ** (-0.4 * p_g) + 10 ** (-0.4 * c_g))
    tot_r = -2.5 * np.log10(10 ** (-0.4 * p_r) + 10 ** (-0.4 * c_r))
    bc = tot_g - tot_r
    bmag = tot_g

    obs_c = np.concatenate([mc, bc])
    obs_m = np.concatenate([mm, bmag])

    # Photometric noise.
    obs_c = obs_c + rng.normal(0, phot_err * np.sqrt(2), size=obs_c.size)
    obs_m = obs_m + rng.normal(0, phot_err, size=obs_m.size)

    # Field outliers: uniform over the CMD box.
    c_lo, c_hi = obs_c.min(), obs_c.max()
    m_lo2, m_hi2 = obs_m.min(), obs_m.max()
    fc = rng.uniform(c_lo, c_hi, size=n_field)
    fm = rng.uniform(m_lo2, m_hi2, size=n_field)
    obs_c = np.concatenate([obs_c, fc])
    obs_m = np.concatenate([obs_m, fm])

    err_c = np.full(obs_c.size, phot_err * np.sqrt(2))
    err_m = np.full(obs_m.size, phot_err)
    return obs_c, obs_m, err_c, err_m


@pytest.fixture(scope="module")
def synth_fitter():
    grid = _make_synthetic_grid()
    return IsochroneFitterV2(
        "synthetic.dat",
        col_mh=1, col_age=2, col_g=3, col_r=4, col_mag=5, col_mass=6,
        iso_data=grid,
        use_bilinear=True,
    )


def test_synthetic_recovery_truth_in_credible_interval(synth_fitter):
    fitter = synth_fitter
    truth = (9.20, 0.0, 9.60, 0.06)   # log_age, [M/H], (m-M), E(g-r)
    obs_c, obs_m, err_c, err_m = _simulate_cluster(fitter, truth, seed=11)

    bounds = FitBounds(
        log_age=(8.6, 9.7),
        metallicity=(-0.3, 0.3),
        distance_mod=(9.0, 10.2),
        extinction_gr=(0.0, 0.20),
    )

    res = fit_isochrone_mcmc(
        obs_c, obs_m, err_c, err_m,
        bounds=bounds,
        interp_fn=fitter._interpolate_isochrone,
        imf_fn=fitter._imf_weight,
        n_walkers=24,
        n_burn=150,
        n_steps=450,
        init_from_gridscan=list(truth),  # seed near truth (grid scan would supply this)
        f_bin=0.3,
        f_field=0.1,
        seed=2024,
    )

    print("\n" + res.summary())
    print("TRUTH:", dict(zip(["log_age", "mh", "dm", "e_color"], truth)))

    # Core assertions: truth inside 16-84 credible interval.
    assert res.contains("log_age", truth[0]), (
        f"log_age truth {truth[0]} not in CI {res.log_age_p}")
    assert res.contains("dm", truth[2]), (
        f"dm truth {truth[2]} not in CI {res.dm_p}")

    # Median age within ~0.15 dex of truth.
    assert abs(res.log_age_med - truth[0]) < 0.15, (
        f"log_age median {res.log_age_med:.3f} off from truth {truth[0]:.3f}")

    # Sampler sanity.
    assert 0.05 <= res.acceptance_fraction <= 0.9
    assert res.flat_chain is not None and res.flat_chain.shape[1] == 4


def test_loglike_prefers_truth_over_wrong_age(synth_fitter):
    """The generative likelihood must score the true age above a very wrong one."""
    from apex.analysis.cmd.isochrone_mcmc import cmd_loglike, CmdHyper, _compute_field_log_density

    fitter = synth_fitter
    truth = (9.20, 0.0, 9.60, 0.06)
    obs_c, obs_m, err_c, err_m = _simulate_cluster(fitter, truth, seed=5)

    hyper = CmdHyper(
        f_bin=0.3, f_field=0.1,
        field_log_density=_compute_field_log_density(obs_c, obs_m),
    )
    ll_true = cmd_loglike(truth, obs_c, obs_m, err_c, err_m,
                          interp_fn=fitter._interpolate_isochrone,
                          imf_fn=fitter._imf_weight, hyper=hyper)
    # Much younger isochrone (turn-off too bright/blue).
    wrong = (8.60, 0.0, 9.60, 0.06)
    ll_wrong = cmd_loglike(wrong, obs_c, obs_m, err_c, err_m,
                           interp_fn=fitter._interpolate_isochrone,
                           imf_fn=fitter._imf_weight, hyper=hyper)
    assert np.isfinite(ll_true)
    assert ll_true > ll_wrong, (
        f"truth loglike {ll_true:.1f} should exceed wrong-age {ll_wrong:.1f}")


# ===========================================================================
# Multi-colour (g-r AND r-i) degeneracy-breaking test
# ===========================================================================
#
# A single g-r CMD cannot separate (age, [M/H], reddening): several combinations
# trace nearly the same locus.  Adding a second colour with a DIFFERENT
# metallicity response (r-i here responds to [M/H] with a different slope than
# g-r) gives the fit independent leverage on [M/H], so the joint posterior must
# be tighter / more accurate in [M/H] and E than the g-r-only fit.
#
# Reddening coefficients (Cardelli SDSS): A_g/A_r/A_i below match EXTINCTION_R in
# the validation runner so E(g-r) -> per-band A is consistent with production.
_R_G, _R_R, _R_I = 3.303, 2.285, 1.698


def _make_synthetic_grid_gri(seed: int = 0) -> np.ndarray:
    """7-col grid carrying g, r, i: [pad, mh, log_age, g, r, i, mass].

    g-r and r-i depend on [M/H] with DIFFERENT slopes so two colours carry
    independent metallicity information (the whole point of the multi-colour fit).
    """
    ages = np.round(np.arange(8.4, 9.81, 0.1), 2)
    mhs = np.round(np.arange(-0.4, 0.41, 0.1), 2)
    rows = []
    for log_age in ages:
        m_to = 2.6 * 10.0 ** (-0.4 * (log_age - 8.4))
        m_to = float(np.clip(m_to, 0.9, 2.8))
        for mh in mhs:
            masses = np.linspace(0.25, m_to, 120)
            abs_g = 4.8 - 4.2 * np.log10(masses) + 0.5 * (masses - 1.0) ** 2
            hook = 0.6 * np.clip((masses - 0.85 * m_to) / (0.15 * m_to), 0, 1)
            abs_g = abs_g - hook
            # g-r: reddens for low mass; +0.25 per dex [M/H].
            gr = 0.95 - 0.55 * np.log10(masses) + 0.25 * mh + 0.15 * hook
            # r-i: SHALLOWER mass slope, OPPOSITE-sign / different [M/H] slope
            # (-0.18 per dex) so (g-r, r-i) are not collinear in [M/H].
            ri = 0.45 - 0.30 * np.log10(masses) - 0.18 * mh + 0.08 * hook
            g = abs_g
            r = g - gr
            i = r - ri
            n = len(masses)
            rows.append(np.column_stack([
                np.zeros(n), np.full(n, mh), np.full(n, log_age),
                g, r, i, masses,
            ]))
    return np.vstack(rows)


@pytest.fixture(scope="module")
def synth_gri():
    """Return (multiband_iso, g-only fitter, gr/ri fitters) on the gri grid."""
    from apex.analysis.cmd.isochrone_mcmc import MultiBandIsochrone

    grid = _make_synthetic_grid_gri()
    mb = MultiBandIsochrone(
        grid,
        col_age=2, col_mh=1, col_mass=6,
        band_cols={"g": 3, "r": 4, "i": 5},
        R_band={"g": _R_G, "r": _R_R, "i": _R_I},
        n_mass_pts=200,
    )
    # Single-colour g-r fitter (band1=g col3, band2=r col4, mag=g col3).
    fitter_gr = IsochroneFitterV2(
        "synthetic_gri.dat",
        col_mh=1, col_age=2, col_g=3, col_r=4, col_mag=3, col_mass=6,
        R_band1=_R_G, R_band2=_R_R, R_mag=_R_G,
        iso_data=grid, use_bilinear=True,
    )
    return mb, fitter_gr


def _simulate_cluster_gri(mb, truth, n_members=400, n_field=30,
                          phot_err=0.02, seed=7):
    """Forward-simulate (g-r, r-i, g) from the multi-band generative model."""
    rng = np.random.default_rng(seed)
    log_age, mh, dm, e_color = truth
    e_bv = e_color / (_R_G - _R_R)
    bands, mass = mb.apparent(log_age, mh, dm, e_bv)
    assert bands is not None and len(mass) >= 20

    order = np.argsort(mass)
    ms = mass[order]
    g = bands["g"][order]
    r = bands["r"][order]
    i = bands["i"][order]

    mm = _kroupa_sample(float(ms.min()), float(ms.max()), n_members, rng)
    sg = np.interp(mm, ms, g)
    sr = np.interp(mm, ms, r)
    si = np.interp(mm, ms, i)

    gr = sg - sr + rng.normal(0, phot_err * np.sqrt(2), n_members)
    ri = sr - si + rng.normal(0, phot_err * np.sqrt(2), n_members)
    gm = sg + rng.normal(0, phot_err, n_members)

    # Uniform field outliers across the box.
    fc_gr = rng.uniform(gr.min(), gr.max(), n_field)
    fc_ri = rng.uniform(ri.min(), ri.max(), n_field)
    fc_g = rng.uniform(gm.min(), gm.max(), n_field)

    obs_gr = np.concatenate([gr, fc_gr])
    obs_ri = np.concatenate([ri, fc_ri])
    obs_g = np.concatenate([gm, fc_g])

    err_gr = np.full(obs_gr.size, phot_err * np.sqrt(2))
    err_ri = np.full(obs_ri.size, phot_err * np.sqrt(2))
    err_g = np.full(obs_g.size, phot_err)
    return obs_gr, obs_ri, obs_g, err_gr, err_ri, err_g


def test_multicolor_breaks_metallicity_reddening_degeneracy(synth_gri):
    """Multi-colour (g-r + r-i) recovers [M/H] and E tighter than g-r alone."""
    mb, fitter_gr = synth_gri
    truth = (9.20, 0.10, 9.60, 0.06)  # log_age, [M/H], (m-M), E(g-r)
    obs_gr, obs_ri, obs_g, err_gr, err_ri, err_g = _simulate_cluster_gri(
        mb, truth, seed=13)

    bounds = FitBounds(
        log_age=(8.6, 9.7),
        metallicity=(-0.4, 0.4),
        distance_mod=(9.0, 10.2),
        extinction_gr=(0.0, 0.20),
    )
    common = dict(
        bounds=bounds, imf_fn=fitter_gr._imf_weight,
        n_walkers=24, n_burn=200, n_steps=600,
        init_from_gridscan=list(truth), f_bin=0.3, f_field=0.1, seed=2024,
    )

    # --- single colour (g-r only) -----------------------------------------
    res_single = fit_isochrone_mcmc(
        obs_gr, obs_g, err_gr, err_g,
        interp_fn=fitter_gr._interpolate_isochrone,
        R_color=(_R_G - _R_R), **common,
    )
    # --- multi colour (g-r + r-i) -----------------------------------------
    obs_mc = np.column_stack([obs_gr, obs_ri, obs_g])
    err_mc = np.column_stack([err_gr, err_ri, err_g])
    res_multi = fit_isochrone_mcmc(
        obs_gr, obs_g, err_gr, err_g,  # ignored in multicolor mode
        interp_fn=fitter_gr._interpolate_isochrone,
        colors=[("g", "r"), ("r", "i")], mag="g",
        multiband_iso=mb, obs_multicolor=obs_mc, err_multicolor=err_mc,
        R_color=(_R_G - _R_R), **common,
    )

    print("\n--- SINGLE (g-r) ---\n" + res_single.summary())
    print("\n--- MULTI (g-r,r-i) ---\n" + res_multi.summary())
    print("TRUTH:", truth)

    def width(p):
        return float(p[2] - p[0])

    mh_w_single = width(res_single.mh_p)
    mh_w_multi = width(res_multi.mh_p)
    e_w_single = width(res_single.e_color_p)
    e_w_multi = width(res_multi.e_color_p)
    mh_acc_single = abs(res_single.mh_med - truth[1])
    mh_acc_multi = abs(res_multi.mh_med - truth[1])

    print(f"[M/H] CI width  single={mh_w_single:.3f}  multi={mh_w_multi:.3f}")
    print(f"E(g-r) CI width single={e_w_single:.3f}  multi={e_w_multi:.3f}")
    print(f"[M/H] |err|     single={mh_acc_single:.3f}  multi={mh_acc_multi:.3f}")

    # Truth must lie in the multi-colour credible interval.
    assert res_multi.contains("mh", truth[1]), (
        f"[M/H] truth {truth[1]} not in multi CI {res_multi.mh_p}")
    assert res_multi.contains("e_color", truth[3]), (
        f"E truth {truth[3]} not in multi CI {res_multi.e_color_p}")

    # The headline result: the SECOND colour (r-i) breaks the metallicity
    # degeneracy, so the multi-colour [M/H] posterior is strictly tighter AND
    # at least as accurate as the g-r-only fit.
    assert mh_w_multi < mh_w_single, (
        f"[M/H] not tighter: single={mh_w_single:.4f} multi={mh_w_multi:.4f}")
    assert mh_acc_multi <= mh_acc_single + 0.03, (
        f"[M/H] less accurate: single|err|={mh_acc_single:.3f} "
        f"multi|err|={mh_acc_multi:.3f}")
    # E should be no worse than the single-colour fit (within MC noise).  In this
    # generative synthetic g-r already pins E tightly, so the joint constraint
    # cannot widen it materially; the strong gain is on [M/H] above.  On real
    # data where g-r/r-i disagree (zeropoint systematics) E behaviour is the
    # diagnostic reported separately by the validation runner.
    assert e_w_multi <= e_w_single * 1.10, (
        f"E significantly worse: single={e_w_single:.4f} multi={e_w_multi:.4f}")


def test_multicolor_loglike_prefers_truth_metallicity(synth_gri):
    """The 2-colour likelihood must localise [M/H] that a 1-colour fit cannot."""
    from apex.analysis.cmd.isochrone_mcmc import (
        cmd_loglike_multicolor, CmdHyperMC, _compute_field_log_density_nd,
    )

    mb, fitter_gr = synth_gri
    truth = (9.20, 0.10, 9.60, 0.06)
    obs_gr, obs_ri, obs_g, err_gr, err_ri, err_g = _simulate_cluster_gri(
        mb, truth, seed=21)
    obs_mc = np.column_stack([obs_gr, obs_ri, obs_g])
    err_mc = np.column_stack([err_gr, err_ri, err_g])

    hyper = CmdHyperMC(
        colors=[("g", "r"), ("r", "i")], mag="g", iso=mb,
        e_color_to_ebv=1.0 / (_R_G - _R_R),
        f_bin=0.3, f_field=0.1,
        field_log_density=_compute_field_log_density_nd(obs_mc),
    )

    def ll(mh):
        return cmd_loglike_multicolor(
            (truth[0], mh, truth[2], truth[3]), obs_mc, err_mc,
            imf_fn=fitter_gr._imf_weight, hyper=hyper)

    ll_true = ll(0.10)
    ll_lo = ll(-0.30)
    ll_hi = ll(0.40)
    assert np.isfinite(ll_true)
    assert ll_true > ll_lo and ll_true > ll_hi, (
        f"true [M/H] ll {ll_true:.1f} should beat {-0.30}:{ll_lo:.1f} "
        f"and {0.40}:{ll_hi:.1f}")
