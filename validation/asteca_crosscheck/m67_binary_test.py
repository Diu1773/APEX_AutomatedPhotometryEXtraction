"""ASteCA on the very same M67 catalogue APEX fitted.

APEX's own fit on this catalogue did not recover the literature values
(age 2.80 Gyr, [M/H] -0.50, E(B-V) 0.173, distance modulus railed at the lower
edge of its parallax window, MCMC acceptance 0.126). Two explanations fit that
equally well and they call for opposite responses:

  * the data cannot constrain it — g/r/i alone leaves the metallicity and
    extinction degenerate, and no tool will do better;
  * APEX's likelihood or sampler is at fault.

Giving a second, independently written forward model the identical stars and
the identical isochrone grid separates them.

What this does and does not compare: ASteCA 0.7 ships the pieces (synthetic
cluster generation, Poisson likelihood ratio of Tremmel et al. 2013) but not
the genetic optimiser of the 2015 paper, so the search below is mine. This is
therefore a comparison of *forward models and likelihoods* on fixed data, not
of two complete published tools.
"""

import json
import time
from pathlib import Path

import numpy as np
from astropy.table import Table          # the venv has astropy, not pandas
from scipy.optimize import differential_evolution

import asteca

WIDE = Path(r"E:\APEX_validation\reprocess\M67\result\cmd_zeropoint"
            r"\median_by_ID_filter_wide_cmd.csv")
APEX_FIT = Path(r"E:\APEX_validation\reprocess\M67\result\cmd_isochrone"
                r"\isochrone_fit.json")
ISO = r"C:\ast_v\iso_m67\parsec_sdss_subset.dat"
OUT = Path(r"C:\ast_v\m67_asteca_result.json")
SEED = 2024

apex = json.loads(APEX_FIT.read_text(encoding="utf-8"))
meta = apex["member_meta"]

df = Table.read(WIDE, format="ascii.csv")

# Reproduce APEX's member cut exactly, so both tools see the same stars.
pm_c, pm_s = meta["pm_center"], meta["pm_sigma"]
plx_c, plx_s = meta["parallax_center"], meta["parallax_sigma"]
nsig, plx_nsig = meta["nsig"], meta["plx_nsig"]
def col(name):
    return np.asarray(df[name], dtype=float)

keep = (
    (np.abs(col("pmra") - pm_c[0]) <= nsig * pm_s[0])
    & (np.abs(col("pmdec") - pm_c[1]) <= nsig * pm_s[1])
    & (np.abs(col("parallax") - plx_c) <= plx_nsig * plx_s)
)
for band in ("g", "r", "i"):
    keep &= np.isfinite(col(f"mag_std_{band}")) & np.isfinite(col(f"mag_std_err_{band}"))
members = df[keep]

def m(name):
    return np.asarray(members[name], dtype=float)
print(f"members: {len(members)}  (APEX fitted {apex['n_stars']})", flush=True)

cluster = asteca.Cluster(
    ra=m("ra_deg"), dec=m("dec_deg"),
    mag=m("mag_std_g"), e_mag=m("mag_std_err_g"),
    color=m("mag_std_g") - m("mag_std_r"),
    e_color=np.hypot(m("mag_std_err_g"), m("mag_std_err_r")),
    color2=m("mag_std_r") - m("mag_std_i"),
    e_color2=np.hypot(m("mag_std_err_r"), m("mag_std_err_i")),
    verbose=1,
)

isochs = asteca.Isochrones(
    model="PARSEC", isochs_path=ISO,
    mag="gmag", color=("gmag", "rmag"), color2=("rmag", "imag"),
    magnitude_effl=4750.0, color_effl=(4750.0, 6220.0),
    color2_effl=(6220.0, 7630.0), verbose=1,
)

synth = asteca.Synthetic(isochs, seed=SEED, verbose=1)
synth.calibrate(cluster)
likelihood = asteca.Likelihood(cluster)

mets = np.array(isochs.met_age_dict["met"], dtype=float)
logas = np.array(isochs.met_age_dict["loga"], dtype=float)
print(f"grid: met {mets.min():.5f}-{mets.max():.5f}, "
      f"loga {logas.min():.2f}-{logas.max():.2f}", flush=True)

# Same windows APEX searched: age 1-8 Gyr, and its parallax-derived dm window.
LOGA = (max(logas.min(), 9.0), min(logas.max(), np.log10(8e9)))
DM = tuple(meta["dm_window"])
AV = (0.0, 1.0)
BOUNDS = [(mets.min(), mets.max()), LOGA, AV, DM]
calls = {"n": 0}


# Match the synthetic sample to the observed one: the default of 100 synthetic
# stars against 409 observed makes the likelihood noisy enough to matter.
N_SYNTH = int(len(members))


def make_cost(binaries: bool):
    """Cost with ASteCA's default binary fraction, or with binaries switched off.

    APEX fits a single isochrone track and does not model unresolved binaries.
    Those sit up to 0.75 mag above the main sequence, so a model without them
    can only absorb that width by reddening or lowering the metallicity — which
    is the direction APEX errs in. Turning ASteCA's binaries off tests whether
    that alone accounts for the gap.
    """
    binfrac = {} if binaries else {"alpha": 0.0, "beta": 0.0}

    def cost(theta):
        met, loga, av, dm = theta
        calls["n"] += 1
        try:
            synth_clust = synth.generate(
                {"met": float(met), "loga": float(loga), "Av": float(av),
                 "DR": 0.0, "Rv": 3.1, "dm": float(dm), **binfrac},
                N_stars=N_SYNTH)
            return float(likelihood.get(synth_clust))
        except Exception:
            return 1e9
    return cost


ZSUN = 0.0152
fits = {}
for label, binaries in (("with_binaries", True), ("no_binaries", False)):
    calls["n"] = 0
    t0 = time.time()
    r = differential_evolution(make_cost(binaries), BOUNDS, seed=SEED,
                               maxiter=60, popsize=12, tol=0.01,
                               polish=False, disp=False)
    fits[label] = {"met_Z": float(r.x[0]), "MH": float(np.log10(r.x[0] / ZSUN)),
                   "age_gyr": float(10 ** r.x[1] / 1e9), "e_bv": float(r.x[2] / 3.1),
                   "distance_mod": float(r.x[3]), "likelihood": float(r.fun),
                   "elapsed_s": round(time.time() - t0, 1), "n_calls": calls["n"]}
    f = fits[label]
    print(f"ASteCA {label:14s} age {f['age_gyr']:.2f} Gyr  [M/H] {f['MH']:+.3f}  "
          f"E(B-V) {f['e_bv']:.3f}  dm {f['distance_mod']:.3f}  "
          f"[{f['elapsed_s']:.0f}s]", flush=True)

res = type("R", (), {"x": [fits["with_binaries"]["met_Z"], 0, 0, 0],
                     "fun": fits["with_binaries"]["likelihood"]})()
met = fits["with_binaries"]["met_Z"]
loga = np.log10(fits["with_binaries"]["age_gyr"] * 1e9)
av = fits["with_binaries"]["e_bv"] * 3.1
dm = fits["with_binaries"]["distance_mod"]
elapsed = fits["with_binaries"]["elapsed_s"]

result = {
    "n_members": int(len(members)),
    "asteca": {
        "met_Z": float(met),
        "log_age": float(loga),
        "age_gyr": float(10 ** loga / 1e9),
        "Av": float(av),
        "e_bv": float(av / 3.1),
        "distance_mod": float(dm),
        "distance_pc": float(10 ** (dm / 5.0 + 1.0)),
        "likelihood": float(res.fun),
        "n_calls": calls["n"],
        "elapsed_s": round(elapsed, 1),
    },
    "apex": apex["summary"],
    "settings": {"bounds": {"met": BOUNDS[0], "loga": LOGA, "Av": AV, "dm": DM},
                 "seed": SEED, "isochrones": ISO,
                 "optimiser": "scipy differential_evolution (not ASteCA's own)"},
}
OUT.write_text(json.dumps(result, indent=1), encoding="utf-8")

print(f"\nASteCA  age {result['asteca']['age_gyr']:.2f} Gyr  Z {met:.5f}  "
      f"E(B-V) {av/3.1:.3f}  dm {dm:.3f}   [{elapsed:.0f}s, {calls['n']} evals, "
      f"N_synth={N_SYNTH}]")
a = apex["summary"]
print(f"APEX    age {a['age_gyr'][1]:.2f} Gyr  [M/H] {a['metallicity'][1]:+.3f}  "
      f"E(B-V) {a['e_bv'][1]:.3f}  dm {a['distance_mod'][1]:.3f}")
print(f"lit.    age ~4 Gyr        [M/H] ~0.0      E(B-V) ~0.04       dm ~9.7")
print(f"\nsaved -> {OUT}")
