"""`pixel_scale_arcsec` unset must produce a message, not a TypeError.

The parameter models parse it with `_as_float_or_none`, so when the user has
not configured a telescope/camera the attribute exists with value None. A
`getattr(P, "pixel_scale_arcsec", np.nan)` therefore returns None — the default
never fires — and `float(None)` raised inside the WCS worker, killing step 5
instead of telling the user what to set.
"""

from __future__ import annotations

import types

import numpy as np
import pytest

from apex.analysis.wcs_solve import configured_pixel_scale


def test_unset_is_nan_not_a_crash():
    """The exact shape the config produces when the key is absent from the TOML."""
    P = types.SimpleNamespace(pixel_scale_arcsec=None)
    assert np.isnan(configured_pixel_scale(P))


def test_attribute_missing_entirely_is_nan():
    assert np.isnan(configured_pixel_scale(types.SimpleNamespace()))


def test_configured_value_passes_through():
    P = types.SimpleNamespace(pixel_scale_arcsec=0.83)
    assert configured_pixel_scale(P) == pytest.approx(0.83)
    # strings survive a hand-edited TOML
    assert configured_pixel_scale(
        types.SimpleNamespace(pixel_scale_arcsec="1.25")) == pytest.approx(1.25)


@pytest.mark.parametrize("bad", [0.0, -1.0, "", "abc", float("nan")])
def test_unusable_values_are_nan(bad):
    """Callers guard with `> 0` / `isfinite`, so every unusable value has to
    land on NaN for those guards to fire."""
    P = types.SimpleNamespace(pixel_scale_arcsec=bad)
    assert np.isnan(configured_pixel_scale(P))


def test_guard_expression_used_by_callers_rejects_the_unset_case():
    """The two guard shapes in wcs_solve must both reject an unset scale."""
    scale = configured_pixel_scale(types.SimpleNamespace(pixel_scale_arcsec=None))
    assert (not np.isfinite(scale)) or scale <= 0      # 1683 / 4092 form
    assert not (scale > 0)                             # 2656 form


def test_real_params_with_no_pixel_scale(tmp_path):
    """End to end through the actual config loader: a TOML without the key."""
    from apex.config.parameters_lc import read_params

    path = tmp_path / "parameters.toml"
    path.write_text('[io]\ndata_dir = "."\n', encoding="utf-8")
    params = read_params(path)
    assert params.P.pixel_scale_arcsec is None         # present, but None
    with pytest.raises(TypeError):                     # what the old code did
        float(getattr(params.P, "pixel_scale_arcsec", np.nan))
    assert np.isnan(configured_pixel_scale(params.P))  # what it does now
