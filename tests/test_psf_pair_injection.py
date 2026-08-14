"""Deliberate blends, so the crowded regime is not a property of the field.

`pair_fraction` places a companion at a chosen separation from a share of the
injections. Worth recording why it is *not* needed: the prediction that a
sparse field could not populate the tight-crowding bins was wrong. An injection
only has to sit beside one real star, and the stratified sampler seeks such
spots — the LCO 1 m frame (median nearest neighbour 13 FWHM, nothing inside
1.5 FWHM) filled the tightest bin with 102 of 400 injections, against 86 on
M13 (2026-08-14).

What pairs buy is control and reach: a separation chosen rather than whatever
the field offers, blends in frames with too few stars to sit beside, and
separations closer than `min_real_sep_fwhm` permits against real stars.

These tests pin both halves: the new mode produces blends, and the old mode is
untouched.
"""

from __future__ import annotations

import numpy as np
import pytest

from apex.benchmark.psf_artificial_stars import sample_stratified_injections

FWHM = 5.0
SHAPE = (2048, 2048)


def _sparse_field(seed: int = 3, n: int = 400) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.uniform(60.0, SHAPE[0] - 60.0, size=(n, 2))


def _sample(**kwargs):
    return sample_stratified_injections(
        SHAPE, _sparse_field(), count=120, fwhm_px=FWHM,
        rng=np.random.default_rng(11), psf_size=25, **kwargs)


def test_default_leaves_the_previous_behaviour_intact():
    table = _sample()
    assert len(table) == 120
    assert not table["is_pair_companion"].any()
    assert (table["pair_id"] == 0).all()
    # Crowding still comes from the real field, and the mirror column agrees.
    assert np.allclose(table["nearest_any_sep_fwhm"],
                       table["nearest_real_sep_fwhm"])


def test_pairs_put_stars_at_the_requested_separation():
    table = _sample(pair_fraction=0.5, pair_separations_fwhm=(1.0,))
    companions = table[table["is_pair_companion"]]
    assert len(companions) == 30, "count must stay the requested count"

    for _, companion in companions.iterrows():
        primary = table[(table["pair_id"] == companion["pair_id"])
                        & (~table["is_pair_companion"])].iloc[0]
        separation = np.hypot(companion["x_true"] - primary["x_true"],
                              companion["y_true"] - primary["y_true"])
        assert separation == pytest.approx(1.0 * FWHM, rel=1e-6)


def test_pairs_populate_the_tight_crowding_bin_a_sparse_field_cannot():
    without = _sample()
    with_pairs = _sample(pair_fraction=0.5, pair_separations_fwhm=(0.8, 1.2))

    tight = "0.75-1.5 FWHM"
    before = int((without["crowding_bin"] == tight).sum())
    after = int((with_pairs["crowding_bin"] == tight).sum())
    assert after > before, f"쌍을 심었는데 최혼잡 구간이 안 늘었다 ({before} -> {after})"
    assert with_pairs["nearest_any_sep_fwhm"].median() < \
        without["nearest_any_sep_fwhm"].median()


def test_companions_still_respect_the_real_star_exclusion():
    """A companion may sit on its primary, never on a real star."""
    real = _sparse_field()
    table = sample_stratified_injections(
        SHAPE, real, count=60, fwhm_px=FWHM, rng=np.random.default_rng(5),
        psf_size=25, min_real_sep_fwhm=0.75, pair_fraction=1.0,
        pair_separations_fwhm=(0.8,))
    assert (table["nearest_real_sep_fwhm"] >= 0.75).all()


def test_a_separation_list_without_a_positive_value_is_refused():
    with pytest.raises(ValueError):
        _sample(pair_fraction=0.5, pair_separations_fwhm=(0.0, -1.0))
