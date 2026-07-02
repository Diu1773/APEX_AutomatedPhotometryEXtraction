"""Tests for the exact standard-color solve used by Step-10 calibration.

The transformation model per filter is ``std_f = inst_f + zp_f + ct_f * C``
with ``C`` a *standard* color. solve_standard_colors must recover the true
standard colors from instrumental magnitudes alone (no catalog input), both
for a shared color pair (B,V -> B-V) and for a chained pair (R -> V-R that
depends on V whose color term uses B-V).
"""

import numpy as np
import pytest

pytest.importorskip("PyQt5")

from apex.gui.workflow.cmd.step10_zeropoint_calibration import solve_standard_colors


RNG = np.random.default_rng(20260703)


def _make_truth(n=200):
    """True standard mags along a fake main sequence."""
    v = RNG.uniform(12.0, 18.0, n)
    bv = RNG.uniform(0.2, 1.6, n)
    vr = 0.55 * bv + 0.05  # arbitrary smooth relation
    b = v + bv
    r = v - vr
    return b, v, r


def _inst_from_truth(std, zp, ct, color_true):
    """Invert std = inst + zp + ct*C  ->  inst = std - zp - ct*C."""
    return std - zp - ct * color_true


FIT = {
    "B": {"zp": -3.90, "ct": +0.12, "color_col": "B_V"},
    "V": {"zp": -3.40, "ct": +0.06, "color_col": "B_V"},
    "R": {"zp": -3.80, "ct": -0.11, "color_col": "V_R"},
}


def test_solve_recovers_shared_pair_exactly():
    b, v, _ = _make_truth()
    bv_true = b - v
    inst = {
        "B": _inst_from_truth(b, FIT["B"]["zp"], FIT["B"]["ct"], bv_true),
        "V": _inst_from_truth(v, FIT["V"]["zp"], FIT["V"]["ct"], bv_true),
    }
    colors = solve_standard_colors(inst, {k: FIT[k] for k in ("B", "V")})
    assert "B_V" in colors
    assert np.allclose(colors["B_V"], bv_true, atol=1e-6)


def test_solve_recovers_chained_pairs():
    b, v, r = _make_truth()
    bv_true, vr_true = b - v, v - r
    inst = {
        "B": _inst_from_truth(b, FIT["B"]["zp"], FIT["B"]["ct"], bv_true),
        "V": _inst_from_truth(v, FIT["V"]["zp"], FIT["V"]["ct"], bv_true),
        "R": _inst_from_truth(r, FIT["R"]["zp"], FIT["R"]["ct"], vr_true),
    }
    colors = solve_standard_colors(inst, FIT)
    assert np.allclose(colors["B_V"], bv_true, atol=1e-5)
    assert np.allclose(colors["V_R"], vr_true, atol=1e-5)


def test_solve_beats_raw_instrumental_color():
    """The solved color must sit on the standard scale; the raw instrumental
    color differs by the zp difference and the color-term stretch."""
    b, v, _ = _make_truth()
    bv_true = b - v
    inst = {
        "B": _inst_from_truth(b, FIT["B"]["zp"], FIT["B"]["ct"], bv_true),
        "V": _inst_from_truth(v, FIT["V"]["zp"], FIT["V"]["ct"], bv_true),
    }
    raw = inst["B"] - inst["V"]
    solved = solve_standard_colors(inst, {k: FIT[k] for k in ("B", "V")})["B_V"]
    # raw color is offset (zp difference ~0.5) and stretched (1 - dct) vs truth
    assert np.median(np.abs(raw - bv_true)) > 0.3
    assert np.median(np.abs(solved - bv_true)) < 1e-6


def test_missing_band_yields_nan_only_for_affected_pair():
    b, v, r = _make_truth(50)
    bv_true, vr_true = b - v, v - r
    inst = {
        "B": _inst_from_truth(b, FIT["B"]["zp"], FIT["B"]["ct"], bv_true),
        "V": _inst_from_truth(v, FIT["V"]["zp"], FIT["V"]["ct"], bv_true),
        "R": _inst_from_truth(r, FIT["R"]["zp"], FIT["R"]["ct"], vr_true),
    }
    inst["R"][:10] = np.nan  # R missing for the first 10 stars
    colors = solve_standard_colors(inst, FIT)
    assert np.isnan(colors["V_R"][:10]).all()
    assert np.isfinite(colors["V_R"][10:]).all()
    assert np.isfinite(colors["B_V"]).all()  # B-V unaffected by missing R


def test_unfitted_band_pair_excluded():
    b, v, _ = _make_truth(50)
    bv_true = b - v
    inst = {
        "B": _inst_from_truth(b, FIT["B"]["zp"], FIT["B"]["ct"], bv_true),
        "V": _inst_from_truth(v, FIT["V"]["zp"], FIT["V"]["ct"], bv_true),
    }
    # R has a fit entry whose pair references an unfitted band: V_R must be absent
    fit = {k: FIT[k] for k in ("B", "V")}
    fit["R"] = {"zp": -3.8, "ct": -0.1, "color_col": "V_R"}
    colors = solve_standard_colors(inst, fit)
    assert "B_V" in colors
    assert "V_R" not in colors
