"""Split-R-hat as the second convergence opinion for the CMD fit.

The autocorrelation criterion is the wrong tool for this posterior. Measured on
M67, `tau` grew from 44 at 400 steps to 202 at 4,000 — it never settles, so the
verdict is always "run longer" no matter how good the answer is. Meanwhile the
estimates were stable to three decimals across 32x2000, 64x4000 and 64x8000
(age 2.805 / 2.803 / 2.802, [M/H] -0.504 / -0.502 / -0.502).

Split-R-hat is reported as a *mixing-efficiency* number, not as a verdict:
it assumes independent chains, and emcee walkers propose from each other. On
this posterior the walkers spread along the degeneracy ridge and each barely
moves, which inflates R-hat to 3-4 even though four different seeds agree to
0.02 dex in [M/H]. These tests pin the statistic itself against cases whose
answer is known by construction; whether a given value should block a fit is a
separate question that the seed test answers better.
"""

from __future__ import annotations

import numpy as np
import pytest

from apex.analysis.cmd.isochrone_mcmc import split_rhat


def _converged(steps=800, walkers=32, dim=3, seed=0):
    """Independent draws from one distribution — halves must agree."""
    rng = np.random.default_rng(seed)
    return rng.normal(size=(steps, walkers, dim))


def test_a_stationary_chain_gives_rhat_near_one():
    rhat = split_rhat(_converged())
    assert rhat.shape == (3,)
    assert np.all(rhat < 1.05), rhat
    assert np.all(rhat > 0.95), rhat


def test_walkers_stuck_in_different_places_are_caught():
    """Half the ensemble sits at a different mode — the classic failure."""
    chain = _converged(seed=1)
    chain[:, : chain.shape[1] // 2, 0] += 40.0
    rhat = split_rhat(chain)
    assert rhat[0] > 1.05, f"a split ensemble should be flagged, got {rhat[0]}"


def test_a_drifting_chain_is_caught():
    """Still moving means the first half and the second half disagree."""
    steps, walkers = 800, 32
    drift = np.linspace(0.0, 25.0, steps)[:, None, None]
    chain = _converged(steps, walkers, 1, seed=2) + drift
    assert split_rhat(chain)[0] > 1.05


def test_one_slow_parameter_does_not_hide_behind_good_ones():
    chain = _converged(dim=3, seed=3)
    chain[:, :, 2] += np.linspace(0.0, 25.0, chain.shape[0])[:, None]
    rhat = split_rhat(chain)
    assert np.all(rhat[:2] < 1.05)
    assert rhat[2] > 1.05
    assert float(np.max(rhat)) == pytest.approx(rhat[2])


@pytest.mark.parametrize("chain", [None,
                                   np.zeros((3, 4, 2)),      # too few steps
                                   np.zeros((10, 4))])       # wrong rank
def test_unusable_input_returns_none_rather_than_guessing(chain):
    assert split_rhat(chain) is None


def test_a_frozen_parameter_does_not_produce_a_fake_verdict():
    """Zero within-chain variance would divide by zero; must not be a number."""
    chain = np.zeros((400, 16, 1))
    rhat = split_rhat(chain)
    assert rhat is None or not np.isfinite(rhat[0])
