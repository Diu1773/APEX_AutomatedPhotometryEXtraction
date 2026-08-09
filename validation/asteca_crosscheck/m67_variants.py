"""What else, besides binaries, is bending the M67 fit?

Binaries were confirmed by ablation: switching ASteCA's off moves it to APEX's
corner. But both tools still miss the literature values, so something is shared.
Each row below is one hypothesis, run through the same forward model on the
same stars, so the shifts are comparable to each other and to the baseline.

  DR        differential reddening across the field, ASteCA's own parameter,
            pinned to 0 in the first comparison
  bright    the brightest 5 % — M67's blue stragglers sit above and blueward of
            the turn-off (the two bluest members here are g-r = -0.26 and +0.26
            against a median of +0.57), and any bright-end non-linearity lands
            in the same place
  faint     the faintest 10 % — star counts are dominated by the lower main
            sequence while the age information sits at the turn-off
"""

import json
import time

import numpy as np
from astropy.table import Table
from scipy.optimize import differential_evolution

import asteca

WIDE = (r"E:\APEX_validation\reprocess\M67\result\cmd_zeropoint"
        r"\median_by_ID_filter_wide_cmd.csv")
ISO = r"C:\ast_v\iso_m67\parsec_sdss_subset.dat"
OUT = r"C:\ast_v\m67_variants.json"
SEED, ZSUN = 2024, 0.0152
DM_WINDOW = (9.631, 9.751)

t = Table.read(WIDE, format="ascii.csv")


def c(name):
    return np.asarray(t[name], dtype=float)


base = (
    (np.abs(c("pmra") + 10.979) <= 0.6)
    & (np.abs(c("pmdec") + 2.916) <= 0.6)
    & (np.abs(c("parallax") - 1.153) <= 0.15)
)
for b in "gri":
    base &= np.isfinite(c(f"mag_std_{b}")) & np.isfinite(c(f"mag_std_err_{b}"))

g_all = c("mag_std_g")
BRIGHT_CUT = float(np.percentile(g_all[base], 5))
FAINT_CUT = float(np.percentile(g_all[base], 90))
print(f"members {base.sum()}   bright cut g<{BRIGHT_CUT:.2f}   "
      f"faint cut g>{FAINT_CUT:.2f}", flush=True)

isochs = asteca.Isochrones(
    model="PARSEC", isochs_path=ISO, mag="gmag",
    color=("gmag", "rmag"), color2=("rmag", "imag"),
    magnitude_effl=4750.0, color_effl=(4750.0, 6220.0),
    color2_effl=(6220.0, 7630.0), verbose=0)
mets = np.array(isochs.met_age_dict["met"], dtype=float)
logas = np.array(isochs.met_age_dict["loga"], dtype=float)


def fit(label, mask, *, binaries=True, dr_free=False):
    sel = np.asarray(mask)
    cluster = asteca.Cluster(
        ra=c("ra_deg")[sel], dec=c("dec_deg")[sel],
        mag=c("mag_std_g")[sel], e_mag=c("mag_std_err_g")[sel],
        color=(c("mag_std_g") - c("mag_std_r"))[sel],
        e_color=np.hypot(c("mag_std_err_g"), c("mag_std_err_r"))[sel],
        color2=(c("mag_std_r") - c("mag_std_i"))[sel],
        e_color2=np.hypot(c("mag_std_err_r"), c("mag_std_err_i"))[sel],
        verbose=0)
    synth = asteca.Synthetic(isochs, seed=SEED, verbose=0)
    synth.calibrate(cluster)
    lkl = asteca.Likelihood(cluster)
    n_synth = int(sel.sum())
    binfrac = {} if binaries else {"alpha": 0.0, "beta": 0.0}

    bounds = [(mets.min(), mets.max()),
              (max(logas.min(), 9.0), min(logas.max(), np.log10(8e9))),
              (0.0, 1.0), DM_WINDOW]
    if dr_free:
        bounds.append((0.0, 0.5))

    def cost(theta):
        met, loga, av, dm = theta[:4]
        dr = float(theta[4]) if dr_free else 0.0
        try:
            return float(lkl.get(synth.generate(
                {"met": float(met), "loga": float(loga), "Av": float(av),
                 "DR": dr, "Rv": 3.1, "dm": float(dm), **binfrac},
                N_stars=n_synth)))
        except Exception:
            return 1e9

    t0 = time.time()
    r = differential_evolution(cost, bounds, seed=SEED, maxiter=60, popsize=12,
                               tol=0.01, polish=False, disp=False)
    out = {"label": label, "n_stars": n_synth, "binaries": binaries,
           "age_gyr": float(10 ** r.x[1] / 1e9),
           "MH": float(np.log10(r.x[0] / ZSUN)),
           "e_bv": float(r.x[2] / 3.1), "dm": float(r.x[3]),
           "DR": float(r.x[4]) if dr_free else 0.0,
           "likelihood": float(r.fun), "elapsed_s": round(time.time() - t0, 1)}
    print(f"{label:26s} n={n_synth:4d}  age {out['age_gyr']:5.2f}  "
          f"[M/H] {out['MH']:+.3f}  E(B-V) {out['e_bv']:.3f}  "
          f"dm {out['dm']:.3f}  DR {out['DR']:.3f}", flush=True)
    return out


bright = base & (g_all >= BRIGHT_CUT)
faint = base & (g_all <= FAINT_CUT)
both = base & (g_all >= BRIGHT_CUT) & (g_all <= FAINT_CUT)

runs = [
    fit("baseline", base),
    fit("+DR free", base, dr_free=True),
    fit("-brightest 5%", bright),
    fit("-faintest 10%", faint),
    fit("-both ends", both),
    fit("-both ends +DR", both, dr_free=True),
    fit("no binaries (ref)", base, binaries=False),
]
json.dump({"bright_cut_g": BRIGHT_CUT, "faint_cut_g": FAINT_CUT, "runs": runs},
          open(OUT, "w"), indent=1)
print(f"\nAPEX  age 2.80  [M/H] -0.502  E(B-V) 0.174  dm 9.633")
print(f"lit.  age ~4.0  [M/H]  0.00   E(B-V) 0.040  dm ~9.70")
print(f"\nsaved -> {OUT}")
