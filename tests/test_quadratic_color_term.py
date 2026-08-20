"""The quadratic colour term is measured every run and applied only if told to.

It used to be fitted per field and adopted whenever the scatter did not get
worse. A quadratic has one more free parameter, so "not worse" costs almost
nothing — but the real problem is what the curvature turned out to be. Measured
on one camera (Moravian C3-61000) across five clusters, 2026-08-21:

    R band ct2:  0.000, 0.000, +0.123, +0.042, -0.244
    B band ct2: -0.088, -0.197,  0.000

A colour term describes how a filter and detector differ from the standard
passband. It belongs to the instrument, so it cannot depend on which cluster is
in the frame. Lending one field's curvature to another made 10 of 24 pairs
*worse* than using none at all, by up to 17 mmag — and the single most
significant curvature in the set (NGC 6811 R, |ct2|/sigma = 46) was the one that
transferred worst, so no sharper significance test would have caught it.

Applying it per field put up to 0.086 mag of colour-dependent error into the
magnitudes (M3 V). Supplying the *instrument's* value instead beat even the
field's own fit: on M13 B, robust scatter went 0.0547 (linear) -> 0.0408 (that
field's own quadratic) -> 0.0373 (a constant shared across fields).
"""

from __future__ import annotations

import numpy as np
import pytest

from apex.analysis.cmd.zeropoint_runner import (
    _quad_coefficient_sigma,
    parse_quadratic_color_terms,
    robust_weighted_polyfit,
)


# ── reading the setting ────────────────────────────────────────────────────

def test_empty_means_no_quadratic_anywhere():
    """The default. A curvature nobody measured is not a curvature."""
    assert parse_quadratic_color_terms("") == {}
    assert parse_quadratic_color_terms(None) == {}
    assert parse_quadratic_color_terms("   ") == {}


def test_per_filter_values_are_read():
    assert parse_quadratic_color_terms("B=-0.10,g=-0.06") == {"B": -0.10, "g": -0.06}
    assert parse_quadratic_color_terms("B = -0.10 ; V = 0.0") == {"B": -0.10, "V": 0.0}


def test_filter_keys_keep_their_case():
    """Johnson is upper, SDSS is lower, and `R` and `r` are different filters."""
    parsed = parse_quadratic_color_terms("R=-0.24,r=+0.12")
    assert parsed == {"R": -0.24, "r": 0.12}


def test_a_typo_is_dropped_rather_than_stored_under_no_filter():
    """`=0.1` would register a curvature under the empty name and never apply."""
    assert parse_quadratic_color_terms("=0.1") == {}
    assert parse_quadratic_color_terms("B=,V=abc,R=-0.05") == {"R": -0.05}
    assert parse_quadratic_color_terms("nonsense") == {}


# ── reporting what the run would have fitted ───────────────────────────────

def _curved(n=400, ct2=-0.10, ct=0.15, zp=-3.9, noise=0.02, seed=3):
    rng = np.random.default_rng(seed)
    x = rng.uniform(-0.2, 1.4, n)
    y = zp + ct * x + ct2 * x * x + rng.normal(0.0, noise, n)
    w = np.full(n, 1.0 / noise**2)
    return x, y, w


def test_sigma_is_small_when_the_curvature_is_well_constrained():
    x, y, w = _curved(ct2=-0.10, noise=0.02)
    coeffs, _n, _s = robust_weighted_polyfit(x, y, w=w, degree=2)
    sigma = _quad_coefficient_sigma(x, y, w, coeffs)
    assert np.isfinite(sigma) and sigma > 0
    assert abs(coeffs[0] - (-0.10)) < 4 * sigma


def test_sigma_grows_when_the_colour_baseline_shrinks():
    """A quadratic over a narrow span is barely constrained, and says so."""
    rng = np.random.default_rng(11)
    sigmas = []
    for span in (1.6, 0.4):
        x = rng.uniform(0.0, span, 400)
        y = -3.9 + 0.15 * x - 0.10 * x * x + rng.normal(0.0, 0.02, 400)
        w = np.full(400, 1.0 / 0.02**2)
        coeffs, _n, _s = robust_weighted_polyfit(x, y, w=w, degree=2)
        sigmas.append(_quad_coefficient_sigma(x, y, w, coeffs))
    assert sigmas[1] > 3 * sigmas[0]


def test_sigma_survives_degenerate_input():
    assert np.isnan(_quad_coefficient_sigma([1.0], [1.0], [1.0], [0.0, 0.0, 0.0]))
    x = np.full(50, 0.5)                       # no colour baseline at all
    sigma = _quad_coefficient_sigma(x, np.zeros(50), np.ones(50), [0.0, 0.0, 0.0])
    assert np.isnan(sigma) or sigma > 1e3


# ── the property that made the old rule wrong ──────────────────────────────

def test_an_extra_parameter_almost_never_makes_the_residual_worse():
    """Why "scatter not worse" could not be a test.

    Fit straight data — no curvature at all — and the quadratic still matches
    or beats the line on its own residuals, because it has a parameter free to
    absorb noise. A rule that adopts on "not worse" adopts on nothing.
    """
    rng = np.random.default_rng(5)
    wins = 0
    trials = 40
    for seed in range(trials):
        rng = np.random.default_rng(seed)
        x = rng.uniform(-0.2, 1.4, 300)
        y = -3.9 + 0.15 * x + rng.normal(0.0, 0.02, 300)   # strictly linear
        w = np.full(300, 1.0 / 0.02**2)
        r1 = y - np.polyval(np.polyfit(x, y, 1, w=np.sqrt(w)), x)
        r2 = y - np.polyval(np.polyfit(x, y, 2, w=np.sqrt(w)), x)
        if float(r2 @ r2) <= float(r1 @ r1) + 1e-12:
            wins += 1
    assert wins == trials, (
        f"{trials - wins} of {trials} straight-line datasets had the quadratic "
        "fit worse — least squares cannot do that, so the comparison is broken"
    )


def test_the_config_row_reaches_the_runner(tmp_path):
    """A key present only in the map is silently dropped by the two-stage loader."""
    import json

    from apex.config.parameters_cmd import read_params

    cfg = tmp_path / "apex_config.json"
    cfg.write_text(json.dumps({
        "io": {"data_dir": str(tmp_path), "result_dir": str(tmp_path)},
        "cmd": {"zp": {"quadratic_color_term": "B=-0.11"}},
    }), encoding="utf-8")
    P = read_params(cfg).P
    assert getattr(P, "zp_quadratic_color_term", None) == "B=-0.11"
    assert parse_quadratic_color_terms(P.zp_quadratic_color_term) == {"B": -0.11}


def test_the_default_config_applies_no_curvature(tmp_path):
    import json

    from apex.config.parameters_cmd import read_params

    cfg = tmp_path / "apex_config.json"
    cfg.write_text(json.dumps({
        "io": {"data_dir": str(tmp_path), "result_dir": str(tmp_path)},
    }), encoding="utf-8")
    P = read_params(cfg).P
    assert parse_quadratic_color_terms(getattr(P, "zp_quadratic_color_term", "")) == {}
