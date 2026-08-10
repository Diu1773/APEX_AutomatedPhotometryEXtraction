"""APEX SYSREM vs PySysRem's core update loop on the identical matrix.

The scorecard's A axis: APEX carries its own Tamuz, Mazeh & Zucker (2005)
implementation with a missing-data/weighting contract. PySysRem
(github.com/stephtdouglas/PySysRem) is the independently written implementation
the exoplanet literature cites, but it is not importable — its output section
uses Python-2 syntax (`print >> f`, sysrem.py line 118) — and its I/O is
coupled to per-star light-curve files. Its *algorithmic core* (sysrem.py lines
73-101: five components, ten alternating weighted-least-squares updates of the
star and epoch coefficients, subtract the outer product) is pure numpy, so it
is reproduced below verbatim in structure, clearly marked. What is being
compared is their update equations against APEX's, on one matrix.

Identical-input protocol — and a trap found the hard way: pre-centring with
APEX's *weighted mean* makes PySysRem's first update degenerate. Its c-init
computes exactly that weighted mean per star, which is then exactly zero, and
the following a-update divides 0/0 into NaN; nanmean then hides the wreckage
(first attempt: component-1 correlation 0.09, delta RMS 48 mmag — garbage,
not disagreement). So the matrix is centred per star by the MEDIAN — PySysRem's
own native prep — and fed identically to both. APEX subtracts its weighted mean
on top; on median-centred data that is a per-star shift of a few mmag and part
of APEX's documented contract. This asymmetry of robustness (APEX's frame-side
c-init survives either centring, PySysRem's star-side init does not) is itself
an audit observation.

Known convergence difference, reported rather than hidden: PySysRem runs a
fixed 10 inner iterations; APEX iterates to tol=1e-6 (max 50). A second
variant runs the PySysRem loop with 50 iterations to separate "different
convergence point" from "different mathematics".

Input: YZ Boo two-night workspace, filter g, via the same Qt-free service the
GUI uses.
"""

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))

from apex.analysis.light_curve.photometry_source_service import (  # noqa: E402
    load_filter_photometry_timeseries,
)
from apex.analysis.light_curve.sysrem import sysrem  # noqa: E402

RESULT_DIR = Path(r"E:\APEX_validation\reprocess\YZBoo_2n\result")
OUT = Path(__file__).parent / "sysrem_compare.json"

table, _src = load_filter_photometry_timeseries(RESULT_DIR, "g")
print(f"long table: {len(table)} rows, "
      f"{table['star_id'].nunique()} stars x {table['frame'].nunique()} frames")

pivot_mag = table.pivot_table(index="star_id", columns="frame",
                              values="mag", aggfunc="first")
pivot_err = table.pivot_table(index="star_id", columns="frame",
                              values="mag_err", aggfunc="first")
# Stars observed in at least 80 % of frames keep the matrix well-conditioned
# on both sides (PySysRem has no explicit missing-data model beyond weights).
coverage = pivot_mag.notna().mean(axis=1)
keep = coverage >= 0.8
mag = pivot_mag[keep].to_numpy(float)
err = pivot_err[keep].to_numpy(float)
print(f"matrix after >=80% coverage cut: {mag.shape[0]} stars x {mag.shape[1]} frames")

# --- identical pre-centring: per-star MEDIAN (PySysRem's native prep) ------
valid_mag = np.isfinite(mag)
r0 = mag - np.nanmedian(np.where(valid_mag, mag, np.nan), axis=1, keepdims=True)
observed = valid_mag

# --- APEX ------------------------------------------------------------------
apex_res = sysrem(r0, err, n_iter=5)
resid_apex = apex_res.residuals


# --- PySysRem core (sysrem.py lines 73-101, structure verbatim) ------------
def pysysrem_core(residuals, errors, n_components=5, n_inner=10):
    residuals = np.where(np.isfinite(residuals), residuals, 0.0)
    errors = np.where(np.isfinite(errors) & (errors > 0), errors, 1e9)
    stars_dim, epoch_dim = residuals.shape
    systematics = []
    for _ in range(n_components):
        c = np.zeros(stars_dim)
        a = np.ones(epoch_dim)
        for _ in range(n_inner):
            for s in range(stars_dim):
                err_squared = errors[s] ** 2
                c[s] = (np.sum(a * residuals[s] / err_squared)
                        / np.sum(a**2 / err_squared))
            for ep in range(epoch_dim):
                err_squared = errors[:, ep] ** 2
                a[ep] = (np.sum(c * residuals[:, ep] / err_squared)
                         / np.sum(c**2 / err_squared))
        syserr = np.outer(c, a)
        residuals = residuals - syserr
        systematics.append(syserr)
    return residuals, systematics


apex_systems = [np.outer(comp.a, comp.c) for comp in apex_res.components]
in_rms = float(np.sqrt(np.nanmean(r0[observed] ** 2))) * 1000
# Agreement must be judged where the data has weight: the unweighted RMS is
# dominated by faint stars whose cells carry almost no weight, and there the
# weighted low-rank problem is indifferent between solutions.
bright = observed & (np.where(np.isfinite(err), err, np.inf) < 0.02)
print(f"input residual RMS {in_rms:.2f} mmag | bright cells (err<20 mmag): "
      f"{int(bright.sum())}/{int(observed.sum())}")

metrics = {"input_rms_mmag": in_rms,
           "n_bright_cells": int(bright.sum()), "n_observed": int(observed.sum())}
for label, n_inner in (("verbatim_10iter", 10), ("extended_50iter", 50)):
    resid_py, sys_py = pysysrem_core(r0, err, n_inner=n_inner)
    finite_frac = float(np.isfinite(resid_py[observed]).mean())
    assert finite_frac > 0.999, f"PySysRem core produced NaNs ({finite_frac})"
    delta_all = (resid_apex - resid_py)[observed]
    delta_bright = (resid_apex - resid_py)[bright]
    # Deflation order is not unique when component strengths are comparable —
    # match components by best absolute correlation, not by index.
    corr = np.zeros((len(apex_systems), len(sys_py)))
    for i, sa in enumerate(apex_systems):
        for j, sp in enumerate(sys_py):
            corr[i, j] = abs(np.corrcoef(sa[observed], sp[observed])[0, 1])
    best_match = corr.max(axis=1)
    metrics[label] = {
        "apex_resid_rms_mmag": float(np.sqrt(np.nanmean(resid_apex[observed]**2))) * 1000,
        "pysysrem_resid_rms_mmag": float(np.sqrt(np.nanmean(resid_py[observed]**2))) * 1000,
        "delta_rms_all_mmag": float(np.sqrt(np.nanmean(delta_all**2))) * 1000,
        "delta_rms_bright_mmag": float(np.sqrt(np.nanmean(delta_bright**2))) * 1000,
        "delta_max_bright_mmag": float(np.nanmax(np.abs(delta_bright))) * 1000,
        "component_best_match_corr": [round(float(v), 4) for v in best_match],
    }
    m = metrics[label]
    print(f"[{label}] resid RMS APEX {m['apex_resid_rms_mmag']:.1f} / "
          f"PySysRem {m['pysysrem_resid_rms_mmag']:.1f} mmag | "
          f"delta RMS all {m['delta_rms_all_mmag']:.2f} / "
          f"bright {m['delta_rms_bright_mmag']:.3f} mmag "
          f"(max {m['delta_max_bright_mmag']:.3f})")
    print(f"          component best-match |corr|: {m['component_best_match_corr']}")

metrics["matrix"] = {"stars": int(mag.shape[0]), "frames": int(mag.shape[1])}
OUT.write_text(json.dumps(metrics, indent=1), encoding="utf-8")
print(f"saved -> {OUT}")
