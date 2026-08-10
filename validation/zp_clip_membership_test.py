"""Why did the zero-point fit keep 564 calibrators one run and 503 the next?

STATUS: hypothesis REFUTED, cause UNRESOLVED. Read this before reusing it.

The Phase 3 clean run reproduced M67's calibrated magnitudes to 1.5 mmag (MAD)
and its colours to 0.5 mmag, but `zp_fit_coefficients.csv` reported a different
number of surviving calibrators: g 564 -> 503, r 641 -> 582, i 692 -> 633, and
the g colour term moved 0.212 -> 0.228.

What is established:

* The measurement-level selection is BIT-IDENTICAL between the two runs.
  `frame_zeropoint_cut_summary.csv` matches on every field (g: 8740 total,
  8470 ref_ok, 8456 delta_ok, 6055 snr_ok, 6000 kept) and `median_snr_ref`
  agrees to six decimals (132.050814). Nothing upstream of the final fit
  selected differently.
* The calibrator table is the same 913 IDs, and step10's fitting code has not
  changed since 2026-08-06.
* So the swing arises inside the final robust polynomial fit alone.

The hypothesis was mechanical: `N` counts inliers after *iterative* sigma
clipping (`robust_weighted_polyfit`: five rounds of "refit, recompute MAD
scatter, keep |r - med| <= 3 sigma"), and clipping membership is a discrete
cascade, so a small perturbation might swing it while leaving the fit stable.

The experiment below injects noise the size of the measured Step 7 difference
(2 mmag) and refits. It does NOT reproduce the swing: N moves 0.9 % where the
runs differ by 10.8 %. The hypothesis is refuted as stated.

It is also not a clean test of the production fit: this reproduction returns
N = 757 and 32 mmag scatter where step10 reports 564 and 24 mmag, so it is
missing the weights and/or the per-star aggregation the real path uses.
Attributing the swing needs step10's own fit instrumented, not this stand-in.

Why it was not pursued further: the science output is unaffected at the level
that matters (1.5 mmag on magnitudes, 0.5 mmag on colours, zero point within
7 mmag against a fit whose own scatter is 24 mmag). `n_fit_calibrators` is a
diagnostic field, not a product. It is recorded as unresolved rather than
explained away.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).absolute().parents[1]   # validation/ is an E: junction
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from apex.gui.workflow.cmd.step10_zeropoint_calibration import (  # noqa: E402
    robust_weighted_polyfit,
)

CAL = Path(r"E:\APEX_validation\reprocess\M67\result\cmd_zeropoint"
           r"\gaia_sdss_calibrator_by_ID.csv")
BAND, COLOR = "g", ("g", "r")
PERTURB_MMAG = 2.0          # the maximum Step 7 difference measured on M67
N_TRIALS = 40
SEED = 2024

cal = pd.read_csv(CAL)
ref_col = f"ref_{BAND}"
inst_col = f"mag_inst_{BAND}"
c1, c2 = (f"ref_{COLOR[0]}", f"ref_{COLOR[1]}")
need = [ref_col, inst_col, c1, c2]
missing = [c for c in need if c not in cal.columns]
if missing:
    raise SystemExit(f"calibrator table lacks {missing}; columns are "
                     f"{sorted(cal.columns)[:20]}")

ok = np.isfinite(cal[need]).all(axis=1)
cal = cal[ok]
y0 = (cal[ref_col] - cal[inst_col]).to_numpy(float)   # the zero point + colour
x = (cal[c1] - cal[c2]).to_numpy(float)
print(f"calibrators with finite {BAND}: {len(cal)}")

coeff0, n0, s0 = robust_weighted_polyfit(x, y0, degree=2)
print(f"unperturbed: N={n0}  scatter={s0*1000:.2f} mmag  "
      f"zp(at colour 0)={coeff0[-1]:.6f}")

rng = np.random.default_rng(SEED)
ns, zps, scatters = [], [], []
for _ in range(N_TRIALS):
    y = y0 + rng.normal(0.0, PERTURB_MMAG / 1000.0, size=y0.size)
    coeff, n, s = robust_weighted_polyfit(x, y, degree=2)
    if coeff is None:
        continue
    ns.append(n)
    zps.append(float(coeff[-1]))
    scatters.append(s)

ns = np.array(ns)
zps = np.array(zps)
print(f"\nwith {PERTURB_MMAG:.0f} mmag noise injected, {len(ns)} trials:")
print(f"  N      : min {ns.min()}  median {int(np.median(ns))}  max {ns.max()}"
      f"   (spread {ns.max()-ns.min()}, {100*(ns.max()-ns.min())/n0:.1f} % of N)")
print(f"  zp     : spread {(zps.max()-zps.min())*1000:.2f} mmag "
      f"(sd {zps.std()*1000:.2f})")
print(f"  scatter: {np.mean(scatters)*1000:.2f} +- {np.std(scatters)*1000:.2f} mmag")

observed_swing = 564 - 503
print(f"\nobserved run-to-run swing on M67 g: {observed_swing} stars "
      f"({100*observed_swing/564:.1f} %)")
print("verdict:", "consistent — membership is chaotic, the fit is not"
      if ns.max() - ns.min() >= observed_swing * 0.5 else
      "NOT reproduced by this mechanism — look further")
