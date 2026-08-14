"""The injection kernel's fit must be judged by its shape, not its parameters.

`fit_frame_moffat` measures the frame's own profile so artificial stars can be
injected with something that is neither engine's model. Its acceptance test used
to read `1.21 < alpha < 9.9`, which quietly assumes the PSF has measurable
wings. It does not always: gamma and alpha are degenerate for a near-Gaussian
core — as alpha grows the Moffat becomes a Gaussian and gamma grows with it, so
both run to whatever bound is set while the profile stays where it is. On the
LCO 0.4 m frame that rejected 39 of 40 stars and stopped the benchmark before
it started, even though every one of those fits reproduced the measured FWHM to
within 2 %.

These tests pin the replacement: a fit is accepted when the profile it implies
matches the frame, and rejected when it does not, whatever alpha happens to be.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO / "validation"))

from run_psf_artificial_stars import fit_frame_moffat  # noqa: E402


def _frame(fwhm: float, beta: float, shape=(600, 600), n_stars=30, seed=1):
    """A synthetic frame of well-separated Moffat stars of known width."""
    rng = np.random.default_rng(seed)
    gamma = fwhm / (2.0 * np.sqrt(2.0 ** (1.0 / beta) - 1.0))
    image = rng.normal(100.0, 2.0, size=shape)
    rows = []
    step = shape[0] // 6
    for index in range(n_stars):
        x = 60 + (index % 5) * step
        y = 60 + (index // 5) * step
        if x >= shape[1] - 60 or y >= shape[0] - 60:
            continue
        half = int(round(4 * fwhm))
        yy, xx = np.mgrid[-half:half + 1, -half:half + 1].astype(float)
        star = 6000.0 * (1.0 + (xx ** 2 + yy ** 2) / gamma ** 2) ** (-beta)
        image[y - half:y + half + 1, x - half:x + half + 1] += star
        rows.append({"x": float(x), "y": float(y), "snr": 400.0,
                     "is_saturated": False})
    return image, pd.DataFrame(rows)


def _run(tmp_path, fwhm, beta):
    image, table = _frame(fwhm, beta)
    tsv = tmp_path / "step7.tsv"
    table.to_csv(tsv, sep="\t", index=False)
    return fit_frame_moffat(image, tsv, fwhm)


def test_a_wingless_near_gaussian_psf_is_accepted(tmp_path):
    """beta 30 is the case that used to fail: no wings, alpha unconstrained."""
    gamma, alpha, n_used = _run(tmp_path, fwhm=5.25, beta=30.0)
    assert n_used >= 8
    implied = 2.0 * gamma * np.sqrt(2.0 ** (1.0 / alpha) - 1.0)
    assert implied == pytest.approx(5.25, rel=0.10)


def test_an_ordinary_moffat_still_recovers_its_parameters(tmp_path):
    gamma, alpha, n_used = _run(tmp_path, fwhm=7.0, beta=3.0)
    assert n_used >= 8
    assert alpha == pytest.approx(3.0, rel=0.25)
    implied = 2.0 * gamma * np.sqrt(2.0 ** (1.0 / alpha) - 1.0)
    assert implied == pytest.approx(7.0, rel=0.10)


def test_a_frame_whose_stars_are_the_wrong_width_is_refused(tmp_path):
    """Told the wrong FWHM, the guard must reject rather than inject garbage."""
    image, table = _frame(fwhm=5.0, beta=3.0)
    tsv = tmp_path / "step7.tsv"
    table.to_csv(tsv, sep="\t", index=False)
    with pytest.raises(RuntimeError, match="부족"):
        fit_frame_moffat(image, tsv, fwhm_px=15.0)
