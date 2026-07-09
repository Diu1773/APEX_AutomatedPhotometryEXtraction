"""Cosmetic correction tests — cosmic-ray + hot pixel, star-protected."""

from __future__ import annotations

import numpy as np
import pytest

astroscrappy = pytest.importorskip("astroscrappy")

from apex.analysis import cosmetic as cos


def _star(img, x, y, flux, sigma=1.4):
    yy, xx = np.mgrid[0:img.shape[0], 0:img.shape[1]]
    img += flux * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * sigma ** 2))


def test_clean_removes_cr_and_hot_preserves_star():
    rng = np.random.default_rng(0)
    img = np.full((128, 128), 300.0)
    _star(img, 64, 64, 8000.0)                    # bright star, well away from artefacts
    img += rng.normal(0, 6.0, img.shape)
    truth_peak = img[64, 64]

    dirty = img.copy()
    dirty[20, 20] += 30000.0                       # cosmic ray (isolated)
    dirty[20, 21] += 25000.0                       # 2-px CR streak
    dirty[100, 100] += 20000.0                     # hot pixel

    cleaned, mask, n = cos.clean_frame(
        dirty, gain=1.5, readnoise=6.0, satlevel=65535.0)

    # artefacts flagged + removed
    assert mask[20, 20] and mask[20, 21] and mask[100, 100]
    assert cleaned[20, 20] < 1000.0 and cleaned[100, 100] < 1000.0
    # star core protected (not flagged) and its flux preserved
    assert not mask[64, 64]
    assert cleaned[64, 64] == pytest.approx(truth_peak, rel=0.02)


def test_star_protect_mask_shields_cores_not_spikes():
    img = np.full((64, 64), 100.0)
    _star(img, 32, 32, 5000.0, sigma=1.6)
    img[10, 10] += 20000.0                          # isolated hot pixel
    m = cos.star_protect_mask(img)
    assert m is not None
    assert m[32, 32]            # star core protected
    assert not m[10, 10]        # isolated spike NOT protected


def test_hot_pixel_mask_from_dark():
    dark = np.full((32, 32), 5.0)
    dark[8, 8] = 500.0
    m = cos.hot_pixel_mask(dark, sigma=6.0)
    assert m is not None and m[8, 8] and m.sum() == 1


def test_dead_pixels_preserved_as_nan():
    img = np.full((48, 48), 300.0)
    img[5, 5] = np.nan
    img[20, 20] += 30000.0
    cleaned, mask, _ = cos.clean_frame(img, gain=1.5, readnoise=6.0)
    assert np.isnan(cleaned[5, 5])          # dead pixel stays NaN
    assert mask[20, 20]                     # CR still removed
