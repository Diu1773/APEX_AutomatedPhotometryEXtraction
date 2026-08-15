"""The aperture correction must be read at r_ap, not near it.

Step 7 builds a 14-point growth curve spanning `0.4*r_ap` to `1.15*r_ref` and
takes the enclosed fraction at the science radius. It used to take it at
whichever grid point was *closest* — but the grid spacing scales with FWHM
(1.5 px at FWHM 8), and r_ap sits on the steep rising part of the curve where a
fraction of a pixel is real flux. Measured across 45 frames of M67 and M13 on
2026-08-16, the nearest-point read made apcorr 0.08-0.10 mag too small, and
almost always in the same direction: at the FWHM this instrument delivers the
nearest grid point lands *outside* r_ap, so the curve is read where more light
is already enclosed.

Nothing failed. `apcorr` is a per-frame scalar and the published magnitudes
come through Step 10, which fits a per-frame zero point and absorbs any uniform
offset exactly — so the CMD never saw it. What was wrong was the reported
quantity itself, and anything read off instrumental magnitudes (depth limits).
"""

from __future__ import annotations

import numpy as np
import pytest


def _enclosed(radii: np.ndarray, curve: np.ndarray, r_ap: float) -> float:
    """The lookup as Step 7 now does it (mirrors forced_photometry.py)."""
    finite = np.isfinite(curve)
    if finite.sum() >= 2:
        return float(np.interp(r_ap, radii[finite], curve[finite]))
    return 0.0


def _nearest(radii: np.ndarray, curve: np.ndarray, r_ap: float) -> float:
    index = int(np.argmin(np.abs(radii - r_ap)))
    return float(curve[index])


def _moffat_growth(radii: np.ndarray, fwhm: float, beta: float = 2.5) -> np.ndarray:
    """Enclosed fraction of a Moffat, normalised at the outermost radius."""
    gamma = fwhm / (2.0 * np.sqrt(2.0 ** (1.0 / beta) - 1.0))
    enclosed = 1.0 - (1.0 + (radii / gamma) ** 2) ** (1.0 - beta)
    return enclosed / enclosed[-1]


def _apex_grid(fwhm: float, n: int = 14):
    r_ap = max(4.0, 0.8 * fwhm)
    r_ref = max(r_ap + 2.0, 2.4 * fwhm)
    return np.linspace(max(2.0, r_ap * 0.4), r_ref * 1.15, n), r_ap


def test_the_value_is_taken_at_r_ap_itself():
    radii = np.array([2.0, 4.0, 6.0, 8.0])
    curve = np.array([0.20, 0.50, 0.80, 1.00])
    assert _enclosed(radii, curve, 5.0) == pytest.approx(0.65)


def test_a_grid_point_that_lands_on_r_ap_is_unchanged():
    radii = np.array([2.0, 4.0, 6.0, 8.0])
    curve = np.array([0.20, 0.50, 0.80, 1.00])
    assert _enclosed(radii, curve, 4.0) == pytest.approx(0.50)


@pytest.mark.parametrize("fwhm", [6.0, 7.0, 8.0])
def test_the_old_lookup_understated_apcorr_at_realistic_seeing(fwhm):
    """The regression this closes, on a curve with a known answer."""
    radii, r_ap = _apex_grid(fwhm)
    curve = _moffat_growth(radii, fwhm)
    interpolated, nearest = _enclosed(radii, curve, r_ap), _nearest(radii, curve, r_ap)
    # Nearest reads further out, so more light is enclosed and apcorr is smaller.
    assert nearest > interpolated
    shift = 2.5 * np.log10((1.0 / interpolated) / (1.0 / nearest))
    assert 0.02 < shift < 0.30, f"FWHM {fwhm}: {shift * 1000:.0f} mmag"


def test_the_grid_spacing_is_what_makes_it_matter():
    """Spacing scales with FWHM; the offset from r_ap does too."""
    offsets = []
    for fwhm in (4.0, 6.0, 8.0):
        radii, r_ap = _apex_grid(fwhm)
        offsets.append(abs(radii[int(np.argmin(np.abs(radii - r_ap)))] - r_ap))
    assert offsets[2] > offsets[0]


def test_a_denser_grid_and_interpolation_agree():
    """Both fix the same thing, so they must land in the same place."""
    fwhm = 7.0
    coarse, r_ap = _apex_grid(fwhm)
    dense, _ = _apex_grid(fwhm, n=1401)
    from_coarse = _enclosed(coarse, _moffat_growth(coarse, fwhm), r_ap)
    from_dense = _nearest(dense, _moffat_growth(dense, fwhm), r_ap)
    assert from_coarse == pytest.approx(from_dense, abs=0.01)


def test_a_curve_with_holes_interpolates_across_the_hole():
    """The gap is skipped, not treated as zero flux."""
    radii = np.array([2.0, 4.0, 6.0, 8.0])
    curve = np.array([0.20, np.nan, 0.80, 1.00])
    assert _enclosed(radii, curve, 5.0) == pytest.approx(0.65)


@pytest.mark.parametrize("curve", [
    np.array([np.nan, 0.5, np.nan]),     # one usable point
    np.array([np.nan, np.nan, np.nan]),  # none
])
def test_too_few_points_means_no_correction_rather_than_a_guess(curve):
    """One point extrapolated across every radius is worse than saying nothing."""
    radii = np.array([2.0, 4.0, 6.0])
    assert _enclosed(radii, curve, 5.5) == 0.0
