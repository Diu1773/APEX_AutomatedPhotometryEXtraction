"""Tests for apex.utils.astro_utils — airmass, BJD, filter helpers."""
from __future__ import annotations

import numpy as np
import pytest

from apex.utils.astro_utils import (
    AIRMASS_FORMULAS,
    DEFAULT_AIRMASS_FORMULA,
    MAX_REASONABLE_AIRMASS,
    MIN_REASONABLE_AIRMASS,
    airmass_from_alt,
    compute_bjd_tdb,
    compute_bjd_tdb_array,
    hardie_airmass,
    is_reasonable_airmass,
    kasten_young_airmass,
    normalize_filter_name,
    secz_airmass,
    young_irvine_airmass,
)


# ---------------------------------------------------------------------------
# normalize_filter_name
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    ("V",       "V"),
    ("B",       "B"),
    ("sdss-r",  "r"),
    ("sdss_g",  "g"),
    ("Rc",      "R"),
    ("IC",      "I"),
    ("halpha",  "Ha"),
    ("OIII",    "OIII"),
    ("[OIII]",  "OIII"),
    ("clear",   "Clear"),
    ("CLR",     "Clear"),
    ("luminance","L"),
    (None,      ""),
    ("",        ""),
    ("UNKNOWN_FILTER", "UNKNOWN_FILTER"),   # pass-through
])
def test_normalize_filter_name(raw, expected):
    assert normalize_filter_name(raw) == expected


# ---------------------------------------------------------------------------
# Airmass formulas — known values
# ---------------------------------------------------------------------------

class TestKastenYoung:
    def test_zenith(self):
        # At zenith (alt=90°) airmass should be ≈1
        am = kasten_young_airmass(90.0)
        assert abs(am - 1.0) < 0.01

    def test_45_deg(self):
        # alt=45° → airmass ≈ 1.41 (sec(45°) ≈ √2)
        am = kasten_young_airmass(45.0)
        assert 1.38 < am < 1.45

    def test_30_deg(self):
        am = kasten_young_airmass(30.0)
        assert 1.90 < am < 2.05

    def test_negative_alt_is_nan(self):
        am = kasten_young_airmass(-5.0)
        assert not np.isfinite(am)

    def test_zero_alt_is_nan(self):
        am = kasten_young_airmass(0.0)
        assert not np.isfinite(am)

    def test_array_input(self):
        alts = np.array([90.0, 60.0, 45.0, -10.0])
        ams = kasten_young_airmass(alts)
        assert ams.shape == (4,)
        assert abs(ams[0] - 1.0) < 0.01
        assert not np.isfinite(ams[3])

    def test_monotone_decreasing_with_altitude(self):
        alts = np.linspace(10, 89, 50)
        ams = kasten_young_airmass(alts)
        assert np.all(np.diff(ams) < 0)


class TestSeczAirmass:
    def test_zenith(self):
        assert abs(secz_airmass(90.0) - 1.0) < 0.001

    def test_60_deg(self):
        # sec(30°) = 2/√3 ≈ 1.155
        am = secz_airmass(60.0)
        assert abs(am - 1.0 / np.cos(np.deg2rad(30.0))) < 1e-9

    def test_negative_returns_nan(self):
        assert not np.isfinite(secz_airmass(-1.0))


class TestYoungIrvine:
    def test_zenith(self):
        am = young_irvine_airmass(90.0)
        assert abs(am - 1.0) < 0.01

    def test_lower_than_secz_at_horizon(self):
        # Young & Irvine corrects downward at high airmass
        yi = young_irvine_airmass(20.0)
        sz = secz_airmass(20.0)
        assert yi < sz


class TestHardie:
    def test_zenith(self):
        am = hardie_airmass(90.0)
        assert abs(am - 1.0) < 0.01

    def test_finite_at_45(self):
        am = hardie_airmass(45.0)
        assert np.isfinite(am) and am > 1.0


class TestAirmassFromAlt:
    def test_default_formula_is_kasten_young(self):
        assert DEFAULT_AIRMASS_FORMULA == "Kasten & Young (1989)"
        am_default = airmass_from_alt(45.0)
        am_ky = kasten_young_airmass(45.0)
        assert am_default == am_ky

    def test_all_named_formulas_give_finite_at_45(self):
        for name in AIRMASS_FORMULAS:
            am = airmass_from_alt(45.0, formula=name)
            assert np.isfinite(am), f"{name} returned non-finite at alt=45°"

    def test_unknown_formula_falls_back_to_default(self):
        am = airmass_from_alt(60.0, formula="nonexistent")
        assert np.isfinite(am)


# ---------------------------------------------------------------------------
# is_reasonable_airmass
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("val,expected", [
    (1.0,  True),
    (1.5,  True),
    (5.9,  True),
    (6.0,  True),
    (0.9,  False),   # below min
    (6.1,  False),   # above max
    (np.nan, False),
    ("x",  False),
    (0.0,  False),
])
def test_is_reasonable_airmass(val, expected):
    assert is_reasonable_airmass(val) == expected


def test_is_reasonable_airmass_custom_bounds():
    assert is_reasonable_airmass(10.0, max_airmass=12.0)
    assert not is_reasonable_airmass(10.0, max_airmass=8.0)


# ---------------------------------------------------------------------------
# BJD_TDB
# ---------------------------------------------------------------------------

# Crab Nebula coords, Palomar Observatory
_RA  = 83.8221
_DEC = 22.0146
_LAT = 33.3563
_LON = -116.8648
_ALT = 1706.0
_JD  = 2459000.5   # arbitrary valid JD


def test_compute_bjd_tdb_returns_finite():
    bjd = compute_bjd_tdb(_JD, _RA, _DEC, _LAT, _LON, _ALT)
    assert np.isfinite(bjd)


def test_compute_bjd_tdb_close_to_jd():
    # Barycentric correction should be within ±0.01 day (≈ ±14 min)
    bjd = compute_bjd_tdb(_JD, _RA, _DEC, _LAT, _LON, _ALT)
    assert abs(bjd - _JD) < 0.01


def test_compute_bjd_tdb_nan_input_returns_nan():
    assert not np.isfinite(compute_bjd_tdb(np.nan, _RA, _DEC, _LAT, _LON))
    assert not np.isfinite(compute_bjd_tdb(_JD, np.nan, _DEC, _LAT, _LON))


def test_compute_bjd_tdb_array_shape():
    jd_arr = np.array([_JD, _JD + 1.0, np.nan])
    out = compute_bjd_tdb_array(jd_arr, _RA, _DEC, _LAT, _LON, _ALT)
    assert out.shape == (3,)
    assert np.isfinite(out[0])
    assert np.isfinite(out[1])
    assert not np.isfinite(out[2])   # NaN in → NaN out


def test_compute_bjd_tdb_array_consistent_with_scalar():
    jds = np.array([_JD, _JD + 0.5])
    vec = compute_bjd_tdb_array(jds, _RA, _DEC, _LAT, _LON, _ALT)
    for i, jd in enumerate(jds):
        scalar = compute_bjd_tdb(jd, _RA, _DEC, _LAT, _LON, _ALT)
        assert abs(vec[i] - scalar) < 1e-9, f"index {i}: vec={vec[i]}, scalar={scalar}"
