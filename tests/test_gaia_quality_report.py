"""A skipped quality cut must announce itself.

`gaia_quality_mask` is permissive when a column is absent, which is correct —
a catalog fetched before Gaia published RUWE should still be usable. Doing it
silently is what cost a night. On M67 the same code, same target and same
night kept 564 calibrators in one run and 503 in the next, with the zero point
moving only 7 mmag so nothing looked wrong. The cause: ESA TAP timed out, APEX
fell back to VizieR, and the VizieR table carries `ruwe` while the ESA query
does not fetch it. One run applied RUWE <= 1.4 and the other skipped it.

Measured on the two real catalogs from those runs: ESA 1900/1900 pass (cut not
applied), VizieR 1749/1900 (151 rejected). These tests pin the reporting so a
future refactor cannot make the skip silent again, and pin the mask itself so
adding the report did not change any decision.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from apex.utils.gaia_quality import gaia_quality_mask, gaia_quality_report


def _catalog(**extra) -> pd.DataFrame:
    base = {
        "gaia_G": [12.0, 14.0, 16.0, 18.0],
        "gaia_BP_RP": [0.5, 0.8, 1.1, 1.4],
    }
    base.update(extra)
    return pd.DataFrame(base)


def test_missing_ruwe_column_is_reported_not_silent():
    """The ESA case: no `ruwe` column, so the cut cannot run."""
    mask, report = gaia_quality_report(_catalog())
    assert mask.all(), "a cut that cannot run must reject nothing"
    ruwe = report["cuts"]["ruwe"]
    assert ruwe["applied"] is False
    assert "ruwe" in ruwe["reason"], "the reason must name the missing column"
    assert report["n_passed"] == 4


def test_present_ruwe_column_applies_the_cut_and_counts_it():
    """The VizieR case: `ruwe` present, so the cut runs and is counted."""
    mask, report = gaia_quality_report(_catalog(ruwe=[1.0, 1.2, 1.5, 3.0]))
    assert list(mask) == [True, True, False, False]
    ruwe = report["cuts"]["ruwe"]
    assert ruwe["applied"] is True
    assert ruwe["threshold"] == 1.4
    assert ruwe["n_rejected"] == 2
    assert report["n_passed"] == 2 and report["n_rejected"] == 2


def test_non_finite_ruwe_is_kept_not_rejected():
    """Permissive within a column too: an unmeasured RUWE is not evidence."""
    mask, report = gaia_quality_report(_catalog(ruwe=[1.0, np.nan, 2.0, np.nan]))
    assert list(mask) == [True, True, False, True]
    assert report["cuts"]["ruwe"]["n_finite"] == 2
    assert report["cuts"]["ruwe"]["n_rejected"] == 1


def test_missing_excess_columns_are_reported():
    _mask, report = gaia_quality_report(_catalog(ruwe=[1.0, 1.0, 1.0, 1.0]))
    excess = report["cuts"]["bp_rp_excess"]
    assert excess["applied"] is False
    assert "phot_bp_rp_excess_factor" in excess["reason"]


def test_excess_cut_runs_when_all_three_columns_present():
    df = _catalog(ruwe=[1.0] * 4,
                  phot_bp_rp_excess_factor=[1.2, 1.25, 1.3, 9.0])
    mask, report = gaia_quality_report(df)
    excess = report["cuts"]["bp_rp_excess"]
    assert excess["applied"] is True
    assert excess["n_rejected"] == int((~mask).sum())
    assert excess["n_rejected"] >= 1, "the 9.0 outlier must be caught"


def test_mask_helper_matches_the_report_mask():
    """Adding the report must not have changed any decision."""
    for df in (
        _catalog(),
        _catalog(ruwe=[1.0, 1.5, np.nan, 2.0]),
        _catalog(ruwe=[1.0] * 4, phot_bp_rp_excess_factor=[1.2, 1.3, 5.0, 1.1]),
    ):
        mask, _ = gaia_quality_report(df)
        assert np.array_equal(mask, gaia_quality_mask(df))


def test_report_counts_are_self_consistent():
    df = _catalog(ruwe=[1.0, 1.5, 2.0, 1.1],
                  phot_bp_rp_excess_factor=[1.2, 1.3, 1.25, 9.0])
    mask, report = gaia_quality_report(df)
    assert report["n_input"] == len(df)
    assert report["n_passed"] == int(mask.sum())
    assert report["n_passed"] + report["n_rejected"] == report["n_input"]
