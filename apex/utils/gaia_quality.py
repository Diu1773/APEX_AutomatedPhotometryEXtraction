"""Gaia DR3 photometric-quality helpers (Qt-free).

Implements the corrected BP/RP flux-excess metric C* and its G-dependent
scatter from Riello et al. 2021 (A&A 649, A3, Table 2 and Eq. 18), plus a
combined calibrator-quality mask (C* + RUWE). In crowded fields the BP/RP
prism windows of faint stars are contaminated by neighbours, biasing BP-RP
and therefore every magnitude transformed from it — C* is the published
detector for exactly that failure mode.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    "gaia_corrected_excess_factor",
    "gaia_cstar_sigma",
    "gaia_quality_mask",
]


def gaia_corrected_excess_factor(bp_rp, excess_factor) -> np.ndarray:
    """Return C* = phot_bp_rp_excess_factor minus the normal-source locus.

    Piecewise polynomial locus f(BP-RP) from Riello+2021 Table 2; a
    well-behaved isolated source has C* ~ 0 regardless of colour.
    """
    x = np.asarray(bp_rp, dtype=float)
    e = np.asarray(excess_factor, dtype=float)
    locus = np.full_like(x, np.nan)
    m1 = np.isfinite(x) & (x < 0.5)
    m2 = np.isfinite(x) & (x >= 0.5) & (x < 4.0)
    m3 = np.isfinite(x) & (x >= 4.0)
    locus[m1] = 1.154360 + 0.033772 * x[m1] + 0.032277 * x[m1] ** 2
    locus[m2] = (
        1.162004
        + 0.011464 * x[m2]
        + 0.049255 * x[m2] ** 2
        - 0.005879 * x[m2] ** 3
    )
    locus[m3] = 1.057572 + 0.140537 * x[m3]
    return e - locus


def gaia_cstar_sigma(g_mag) -> np.ndarray:
    """1-sigma scatter of C* for well-behaved sources at magnitude G.

    Riello+2021 Eq. 18: sigma(G) = 0.0059898 + 8.817481e-12 * G^7.618399.
    """
    g = np.asarray(g_mag, dtype=float)
    return 0.0059898 + 8.817481e-12 * np.power(g, 7.618399)


def gaia_quality_mask(
    df: pd.DataFrame,
    *,
    ruwe_max: float = 1.4,
    cstar_nsigma: float = 3.0,
    bp_rp_col: str = "gaia_BP_RP",
    g_col: str = "gaia_G",
    excess_col: str = "phot_bp_rp_excess_factor",
    ruwe_col: str = "ruwe",
) -> np.ndarray:
    """Boolean mask of rows passing the Gaia quality cuts.

    Permissive on missing data: a cut only rejects rows where the needed
    columns are present and finite, so catalogs built before these columns
    were fetched behave exactly as today.
    """
    n = len(df)
    mask = np.ones(n, dtype=bool)

    if ruwe_col in df.columns:
        ruwe = pd.to_numeric(df[ruwe_col], errors="coerce").to_numpy(float)
        mask &= (~np.isfinite(ruwe)) | (ruwe <= float(ruwe_max))

    if excess_col in df.columns and bp_rp_col in df.columns and g_col in df.columns:
        cstar = gaia_corrected_excess_factor(
            pd.to_numeric(df[bp_rp_col], errors="coerce"),
            pd.to_numeric(df[excess_col], errors="coerce"),
        )
        sigma = gaia_cstar_sigma(pd.to_numeric(df[g_col], errors="coerce"))
        bad = (
            np.isfinite(cstar)
            & np.isfinite(sigma)
            & (np.abs(cstar) > float(cstar_nsigma) * sigma)
        )
        mask &= ~bad
    return mask
