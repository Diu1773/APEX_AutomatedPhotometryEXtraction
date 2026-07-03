import numpy as np
import pandas as pd

from apex.utils.gaia_quality import (
    gaia_corrected_excess_factor,
    gaia_cstar_sigma,
    gaia_quality_mask,
)


def _locus(x):
    """Reference locus f(BP-RP) re-derived independently for the test."""
    if x < 0.5:
        return 1.154360 + 0.033772 * x + 0.032277 * x**2
    if x < 4.0:
        return 1.162004 + 0.011464 * x + 0.049255 * x**2 - 0.005879 * x**3
    return 1.057572 + 0.140537 * x


def test_cstar_is_zero_on_the_locus():
    bp_rp = np.array([0.0, 0.4, 0.5, 1.0, 2.5, 3.9, 4.0, 5.0])
    excess = np.array([_locus(x) for x in bp_rp])
    cstar = gaia_corrected_excess_factor(bp_rp, excess)
    assert np.allclose(cstar, 0.0, atol=1e-12)


def test_cstar_sigma_grows_with_magnitude():
    g = np.array([6.0, 14.0, 17.0, 19.0, 21.0])
    s = gaia_cstar_sigma(g)
    assert np.all(np.diff(s) > 0)
    assert abs(s[0] - 0.0059898) < 1e-4  # bright limit ~ the constant term


def test_quality_mask_rejects_contaminated_and_high_ruwe():
    bp_rp = np.array([1.0, 1.0, 1.0, 1.0])
    g = np.array([18.0, 18.0, 18.0, 18.0])
    clean = _locus(1.0)
    sigma = float(gaia_cstar_sigma(18.0))
    df = pd.DataFrame(
        {
            "gaia_BP_RP": bp_rp,
            "gaia_G": g,
            # rows: clean / contaminated (10 sigma) / clean / clean
            "phot_bp_rp_excess_factor": [clean, clean + 10 * sigma, clean, clean],
            # rows: ok / ok / bad astrometry / missing
            "ruwe": [1.0, 1.0, 2.5, np.nan],
        }
    )
    mask = gaia_quality_mask(df)
    assert mask.tolist() == [True, False, False, True]  # NaN ruwe is permissive


def test_quality_mask_permissive_when_columns_missing():
    df = pd.DataFrame({"gaia_BP_RP": [1.0, 2.0], "gaia_G": [15.0, 16.0]})
    assert gaia_quality_mask(df).all()
