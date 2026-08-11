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
    "gaia_quality_report",
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


def gaia_quality_report(
    df: pd.DataFrame,
    *,
    ruwe_max: float = 1.4,
    cstar_nsigma: float = 3.0,
    bp_rp_col: str = "gaia_BP_RP",
    g_col: str = "gaia_G",
    excess_col: str = "phot_bp_rp_excess_factor",
    ruwe_col: str = "ruwe",
) -> tuple[np.ndarray, dict]:
    """Quality mask plus a record of which cuts were actually applied.

    Being permissive on missing columns is right — a catalog fetched before
    these columns existed should still work — but doing it *silently* cost a
    night of investigation. On M67 the same code, same target and same night
    kept 564 calibrators in one run and 503 in the next; the zero point moved
    only 7 mmag, so nothing looked broken, and the cause was invisible from the
    outputs. It was this: ESA TAP timed out and APEX fell back to VizieR. Both
    serve the same Gaia DR3 (`gaiadr3.gaia_source` vs VizieR `I/355/gaiadr3`);
    the difference is that APEX's VizieR query SELECTs `RUWE` and its ESA query
    does not — the ESA table has the column, the query simply omits it. So one
    run applied the RUWE cut and the other skipped it, removing 84 of 910
    calibrators (9.2 %) — matching the observed drop of 59-61 per band.

    Which cuts run must therefore be visible in the outputs, not inferred from
    a catalog-provenance JSON three directories away. The returned dict says
    for each cut whether it was applied, why not if it was skipped, and how
    many rows it rejected.
    """
    n = len(df)
    mask = np.ones(n, dtype=bool)
    report: dict = {"n_input": int(n), "cuts": {}}

    if ruwe_col in df.columns:
        ruwe = pd.to_numeric(df[ruwe_col], errors="coerce").to_numpy(float)
        rejected = np.isfinite(ruwe) & (ruwe > float(ruwe_max))
        mask &= ~rejected
        report["cuts"]["ruwe"] = {
            "applied": True, "threshold": float(ruwe_max),
            "n_rejected": int(rejected.sum()),
            "n_finite": int(np.isfinite(ruwe).sum()),
        }
    else:
        report["cuts"]["ruwe"] = {
            "applied": False, "reason": f"column '{ruwe_col}' not in catalog"}

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
        report["cuts"]["bp_rp_excess"] = {
            "applied": True, "nsigma": float(cstar_nsigma),
            "n_rejected": int(bad.sum()),
        }
    else:
        missing = [c for c in (excess_col, bp_rp_col, g_col)
                   if c not in df.columns]
        report["cuts"]["bp_rp_excess"] = {
            "applied": False, "reason": f"missing column(s): {missing}"}

    report["n_passed"] = int(mask.sum())
    report["n_rejected"] = int(n - mask.sum())
    return mask, report


def gaia_quality_mask(df: pd.DataFrame, **kwargs) -> np.ndarray:
    """Boolean mask of rows passing the Gaia quality cuts.

    Permissive on missing data: a cut only rejects rows where the needed
    columns are present and finite, so catalogs built before these columns
    were fetched behave exactly as today. Use :func:`gaia_quality_report` when
    the caller should record *which* cuts ran — see its docstring for why that
    matters.
    """
    mask, _ = gaia_quality_report(df, **kwargs)
    return mask
