"""Does the RUWE cut remove the binaries step12 is trying to model?

step12 fits M13's CMD with `f_bin = 0.3` and a four-ratio binary ridge: it
assumes ~30 % of the stars are unresolved binaries sitting up to 0.75 mag above
the main sequence, and it uses them. RUWE > 1.4 is, by construction, a detector
for unresolved binaries — the astrometric fit is poor because two stars are
wobbling around a common centre. Filtering on RUWE and then fitting with
f_bin = 0.3 would be assuming a population that has just been removed.

A first, crude check compared median B magnitude inside a colour slice and
found the cut stars only 0.06 mag fainter — no signal. That test cannot answer
the question: binaries are displaced *perpendicular to the main sequence*, and
a median magnitude inside a wide colour slice averages that displacement away.

This measures the displacement properly. Fit a main-sequence ridge line
(robust, iterative) to the blue main sequence, then compare the distribution of
vertical offsets from that ridge for stars the cut keeps versus stars it
removes. If RUWE preferentially removes binaries, the removed stars sit
systematically ABOVE the ridge (brighter at fixed colour).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).absolute().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from apex.utils.gaia_quality import (  # noqa: E402
    gaia_corrected_excess_factor,
    gaia_cstar_sigma,
)

CMD = Path(r"E:\APEX_validation\cstar_test\M13_cmd\result\cmd_zeropoint"
           r"\median_by_ID_filter_wide_cmd.csv")
OUT = REPO / "validation" / "ruwe_binary_conflict.json"
MS_COLOR = (-0.15, 0.35)     # blue main sequence in B-V for M13
RUWE_MAX, CSTAR_NSIG = 1.4, 3.0


def ridge(colour: np.ndarray, mag: np.ndarray, deg: int = 2,
          iters: int = 5) -> np.ndarray:
    """Robust main-sequence ridge: fit, clip the bright half, refit.

    Clipping is one-sided on purpose. Binaries lie above the sequence, so a
    symmetric clip would let them pull the ridge up and hide the very effect
    being measured.
    """
    m = np.isfinite(colour) & np.isfinite(mag)
    coeffs = np.polyfit(colour[m], mag[m], deg)
    for _ in range(iters):
        resid = mag - np.polyval(coeffs, colour)
        s = 1.4826 * np.nanmedian(np.abs(resid[m] - np.nanmedian(resid[m])))
        keep = m & (resid > -1.5 * s)          # drop the bright tail
        if keep.sum() < deg + 5 or keep.sum() == m.sum():
            break
        m = keep
        coeffs = np.polyfit(colour[m], mag[m], deg)
    return coeffs


d = pd.read_csv(CMD)
bv = (pd.to_numeric(d["mag_std_B"], errors="coerce")
      - pd.to_numeric(d["mag_std_V"], errors="coerce")).to_numpy(float)
B = pd.to_numeric(d["mag_std_B"], errors="coerce").to_numpy(float)
ruwe = pd.to_numeric(d["ruwe"], errors="coerce").to_numpy(float)
cstar = gaia_corrected_excess_factor(d["gaia_BP_RP"], d["phot_bp_rp_excess_factor"])
sigma = gaia_cstar_sigma(d["gaia_G"])

ms = np.isfinite(bv) & np.isfinite(B) & (bv > MS_COLOR[0]) & (bv < MS_COLOR[1])
coeffs = ridge(bv[ms], B[ms])
# Negative offset = brighter than the ridge = where binaries live.
offset = B - np.polyval(coeffs, bv)

cuts = {
    "RUWE > 1.4": np.isfinite(ruwe) & (ruwe > RUWE_MAX),
    "C* > 3 sigma": (np.isfinite(cstar) & np.isfinite(sigma)
                     & (np.abs(cstar) > CSTAR_NSIG * sigma)),
}

print(f"main-sequence sample: {int(ms.sum())} stars, "
      f"B-V in {MS_COLOR}, ridge deg 2")
print(f"{'cut':16s}{'removed':>9s}{'kept med':>10s}{'removed med':>13s}"
      f"{'shift':>9s}{'p(KS)':>9s}")

from scipy import stats  # noqa: E402

results = {}
for name, bad in cuts.items():
    a = offset[ms & ~bad]
    b = offset[ms & bad]
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if b.size < 10:
        print(f"{name:16s}  too few removed ({b.size})")
        continue
    med_a, med_b = float(np.median(a)), float(np.median(b))
    ks = stats.ks_2samp(a, b)
    print(f"{name:16s}{b.size:9d}{med_a:10.3f}{med_b:13.3f}"
          f"{med_b - med_a:+9.3f}{ks.pvalue:9.2g}")
    results[name] = {
        "n_kept": int(a.size), "n_removed": int(b.size),
        "median_offset_kept": med_a, "median_offset_removed": med_b,
        "shift_mag": med_b - med_a, "ks_p": float(ks.pvalue),
    }

print("\nnegative shift = removed stars are BRIGHTER than the ridge "
      "= the binary sequence")
OUT.write_text(json.dumps(
    {"ms_colour_range": MS_COLOR, "ridge_coeffs": coeffs.tolist(),
     "n_ms": int(ms.sum()), "cuts": results}, indent=1), encoding="utf-8")
print(f"saved -> {OUT}")
