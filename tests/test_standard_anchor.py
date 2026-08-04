"""Unit tests for apex.analysis.cmd.standard_anchor (no network)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from apex.analysis.cmd.standard_anchor import (
    MAX_SANE_OFFSET_MAG,
    anchor_qc_frame,
    apply_anchor,
    compute_anchor,
    resolve_band_magnitude,
)

RNG = np.random.default_rng(42)


def _make_field(n=80, offsets=None, scatter=0.02):
    """Synthetic wide table + standard catalog on the same sky positions.

    The standard catalog stores Vmag + colours (the MMJ93 layout), so the
    U/B chains are exercised, with known per-band offsets injected into the
    APEX side: mag_std = standard - offset  =>  compute_anchor must recover
    +offset.
    """
    offsets = offsets or {"U": -0.131, "B": +0.051, "V": +0.008}
    ra = 132.8 + RNG.uniform(-0.1, 0.1, n)
    dec = 11.8 + RNG.uniform(-0.1, 0.1, n)
    v = RNG.uniform(12.0, 16.0, n)
    bv = RNG.uniform(0.4, 1.4, n)
    ub = RNG.uniform(-0.1, 1.0, n)
    std = pd.DataFrame({
        "RAJ2000": ra, "DEJ2000": dec,
        "Vmag": v, "B-V": bv, "U-B": ub,
    })
    u_true = v + bv + ub
    b_true = v + bv
    wide = pd.DataFrame({
        "ID": np.arange(n),
        "ra_deg": ra, "dec_deg": dec,
        "mag_std_U": u_true - offsets["U"] + RNG.normal(0, scatter, n),
        "mag_std_B": b_true - offsets["B"] + RNG.normal(0, scatter, n),
        "mag_std_V": v - offsets["V"] + RNG.normal(0, scatter, n),
    })
    return wide, std, offsets


def test_recovers_injected_offsets():
    wide, std, offsets = _make_field()
    res = compute_anchor(wide, std, "TEST/CAT", match_radius_arcsec=1.0,
                         min_stars=20)
    assert res.n_matched == len(wide)
    got = {b.band: b for b in res.bands}
    for band, true_off in offsets.items():
        assert got[band].applied
        assert got[band].offset == pytest.approx(true_off, abs=0.01)
        assert got[band].sigma < 0.05


def test_apply_anchor_shifts_only_mag_std():
    wide, std, offsets = _make_field()
    wide["mag_inst_U"] = wide["mag_std_U"] - 5.0
    res = compute_anchor(wide, std, "TEST/CAT")
    out = apply_anchor(wide, res)
    got = {b.band: b for b in res.bands}
    for band in offsets:
        shift = out[f"mag_std_{band}"] - wide[f"mag_std_{band}"]
        assert np.allclose(shift, got[band].offset)
    assert np.array_equal(out["mag_inst_U"], wide["mag_inst_U"])


def test_min_stars_guard_blocks_application():
    wide, std, _ = _make_field(n=10)
    res = compute_anchor(wide, std, "TEST/CAT", min_stars=20)
    assert all(not b.applied for b in res.bands)
    assert any("not applied" in w for w in res.warnings)


def test_insane_offset_not_applied():
    wide, std, _ = _make_field(offsets={"U": 2.5, "B": 0.05, "V": 0.0})
    res = compute_anchor(wide, std, "TEST/CAT")
    got = {b.band: b for b in res.bands}
    assert not got["U"].applied
    assert abs(got["U"].offset) > MAX_SANE_OFFSET_MAG
    assert got["B"].applied
    out = apply_anchor(wide, res)
    assert np.array_equal(out["mag_std_U"], wide["mag_std_U"])  # untouched


def test_resolve_band_direct_column_wins():
    std = pd.DataFrame({"Umag": [10.0], "Vmag": [9.0], "B-V": [0.5], "U-B": [0.1]})
    assert resolve_band_magnitude(std, "U").iloc[0] == 10.0
    # chain path when the direct column is absent
    std2 = std.drop(columns=["Umag"])
    assert resolve_band_magnitude(std2, "U").iloc[0] == pytest.approx(9.6)
    assert resolve_band_magnitude(std2, "R") is None  # no V-R column


def test_qc_frame_shape():
    wide, std, _ = _make_field()
    res = compute_anchor(wide, std, "TEST/CAT")
    qc = anchor_qc_frame(res)
    assert set(qc["band"]) == {"U", "B", "V"}
    assert {"offset_mag", "robust_sigma", "n_matched", "applied"} <= set(qc.columns)


def test_no_position_rows_warns():
    wide = pd.DataFrame({"mag_std_V": [12.0], "ra_deg": [np.nan], "dec_deg": [np.nan]})
    std = pd.DataFrame({"RAJ2000": [10.0], "DEJ2000": [10.0], "Vmag": [12.0]})
    res = compute_anchor(wide, std, "TEST/CAT")
    assert res.n_matched == 0
    assert res.warnings
