"""Single- vs multi-colour MCMC comparison for real cluster CMDs (headless).

Runs the SAME generative MCMC twice on one cluster: once on g-r only and once on
(g-r, r-i, g), with identical walker/step budgets, then prints a side-by-side
posterior comparison.  Standalone (Qt-free); used for the multi-colour
degeneracy-breaking validation on M67 / M37.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from apex.analysis.cmd.isochrone_fitter_v2 import FitBounds, IsochroneFitterV2
from apex.analysis.cmd.isochrone_mcmc import MultiBandIsochrone, fit_isochrone_mcmc
from validation.cmd_step12_realdata import (
    EXTINCTION_R,
    extract_cmd_arrays,
    extract_multicolor_arrays,
    load_compact_isochrone,
    load_multiband_isochrone,
    resolve_wide_csv,
)


def _fmt(p):
    if p is None:
        return "n/a"
    lo, mid, hi = p
    return f"{mid:+.3f} (+{hi - mid:.3f}/-{mid - lo:.3f})"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("result_dir")
    ap.add_argument("--cluster", default="cluster")
    ap.add_argument("--age-bounds", nargs=2, type=float, required=True)
    ap.add_argument("--mh-bounds", nargs=2, type=float, required=True)
    ap.add_argument("--dm-bounds", nargs=2, type=float, required=True)
    ap.add_argument("--ecolor-bounds", nargs=2, type=float, required=True)
    ap.add_argument("--n-walkers", type=int, default=40)
    ap.add_argument("--n-burn", type=int, default=1000)
    ap.add_argument("--n-steps", type=int, default=3500)
    ap.add_argument("--data-snr-min", type=float, default=10.0)
    ap.add_argument("--fit-snr-min", type=float, default=5.0)
    ap.add_argument("--max-stars", type=int, default=4000)
    ap.add_argument("--mag", default="g")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    result_dir = Path(args.result_dir)
    wide_csv = resolve_wide_csv(result_dir)
    df = pd.read_csv(wide_csv)

    iso_file = ROOT / "isochrone" / "PARSEC" / "sdss" / "iso_data_cmd39.dat"
    cache_dir = ROOT / "validation" / "reports" / "iso_cache"

    age_b = tuple(args.age_bounds)
    mh_b = tuple(args.mh_bounds)
    bounds = FitBounds(
        log_age=age_b, metallicity=mh_b,
        distance_mod=tuple(args.dm_bounds),
        extinction_gr=tuple(args.ecolor_bounds),
    )
    mag = args.mag

    common = dict(
        bounds=bounds, n_walkers=args.n_walkers,
        n_burn=args.n_burn, n_steps=args.n_steps, seed=args.seed, progress=False,
        f_bin=0.3, f_field=0.1,
    )

    # ---- single colour (g-r) ------------------------------------------------
    color = ("g", "r")
    R_color = EXTINCTION_R["g"] - EXTINCTION_R["r"]
    obs_c, obs_m, ce, me, _, dmeta = extract_cmd_arrays(
        df, color, mag, data_snr_min=args.data_snr_min,
        max_stars=args.max_stars, seed=args.seed,
    )
    iso_data, _ = load_compact_isochrone(iso_file, color, mag, age_b, mh_b, cache_dir=cache_dir)
    fitter = IsochroneFitterV2(
        iso_file, col_mh=1, col_age=2, col_g=3, col_r=4, col_mag=5, col_mass=6,
        R_band1=EXTINCTION_R["g"], R_band2=EXTINCTION_R["r"], R_mag=EXTINCTION_R[mag],
        iso_data=iso_data, use_bilinear=False,
    )
    grid = fitter.fit_grid_scan(obs_c, obs_m, ce, me, bounds=bounds, snr_min=args.fit_snr_min)
    gb = grid.best_fit
    init = [gb.log_age, gb.metallicity, gb.distance_mod, gb.extinction_gr]

    ce_c = np.clip(ce, 0.01, None)
    me_c = np.clip(me, 0.01, None)
    snr = np.minimum(1.0 / ce_c, 1.0 / me_c)
    mask = (snr > args.fit_snr_min) & np.isfinite(obs_c) & np.isfinite(obs_m)

    t0 = time.perf_counter()
    res_single = fit_isochrone_mcmc(
        obs_c[mask], obs_m[mask], ce_c[mask], me_c[mask],
        interp_fn=fitter._interpolate_isochrone, imf_fn=fitter._imf_weight,
        init_from_gridscan=init, R_color=R_color, **common,
    )
    t_single = time.perf_counter() - t0

    # ---- multi colour (g-r, r-i, g) -----------------------------------------
    colors = [("g", "r"), ("r", "i")]
    bands = ["g", "r", "i"]
    obs_mc, err_mc, axes, mc_meta = extract_multicolor_arrays(
        df, colors, mag, data_snr_min=args.data_snr_min,
        max_stars=args.max_stars, seed=args.seed,
    )
    mb_data, band_col, _ = load_multiband_isochrone(iso_file, bands, age_b, mh_b, cache_dir=cache_dir)
    mb = MultiBandIsochrone(
        mb_data, col_age=2, col_mh=1, col_mass=mb_data.shape[1] - 1,
        band_cols=band_col, R_band={b: EXTINCTION_R[b] for b in bands},
    )
    t0 = time.perf_counter()
    res_multi = fit_isochrone_mcmc(
        obs_c[mask], obs_m[mask], ce_c[mask], me_c[mask],  # ignored
        interp_fn=fitter._interpolate_isochrone, imf_fn=fitter._imf_weight,
        init_from_gridscan=init, R_color=R_color,
        colors=colors, mag=mag, multiband_iso=mb,
        obs_multicolor=obs_mc, err_multicolor=err_mc, **common,
    )
    t_multi = time.perf_counter() - t0

    # ---- report -------------------------------------------------------------
    def block(tag, r, nst, secs):
        return {
            "tag": tag,
            "n_stars": int(nst),
            "log_age": list(r.log_age_p),
            "age_gyr": list(r.age_gyr_p),
            "metallicity": list(r.mh_p),
            "distance_mod": list(r.dm_p),
            "e_color": list(r.e_color_p),
            "e_bv": list(r.e_bv_p) if r.e_bv_p else None,
            "acceptance": float(r.acceptance_fraction),
            "autocorr": [float(x) for x in r.autocorr_time] if r.autocorr_time is not None else None,
            "convergence_ok": bool(r.convergence_ok),
            "elapsed_sec": float(secs),
        }

    out = {
        "cluster": args.cluster,
        "bands_used_single": dmeta.get("columns"),
        "bands_used_multi": mc_meta.get("columns"),
        "n_single": int(mask.sum()),
        "n_multi": int(obs_mc.shape[0]),
        "single": block("single g-r", res_single, mask.sum(), t_single),
        "multi": block("multi g-r,r-i", res_multi, obs_mc.shape[0], t_multi),
    }

    print(f"\n================ {args.cluster}: SINGLE vs MULTI ================")
    print(f"  N stars: single={out['n_single']}  multi={out['n_multi']}")
    print(f"  band cols single={out['bands_used_single']}")
    print(f"  band cols multi ={out['bands_used_multi']}")
    for key, label in [("log_age", "log(Age)"), ("age_gyr", "Age[Gyr]"),
                       ("metallicity", "[M/H]"), ("distance_mod", "(m-M)"),
                       ("e_color", "E(g-r)"), ("e_bv", "E(B-V)")]:
        s = out["single"][key]
        m = out["multi"][key]
        print(f"  {label:9s}  single={_fmt(s):28s}  multi={_fmt(m)}")

    def width(p):
        return None if p is None else (p[2] - p[0])
    print("  --- CI width (p84-p16): tighter = better-constrained ---")
    for key, label in [("metallicity", "[M/H]"), ("e_color", "E(g-r)"),
                       ("log_age", "log(Age)")]:
        ws = width(out["single"][key]); wm = width(out["multi"][key])
        verdict = "MULTI TIGHTER" if wm < ws else "single tighter"
        print(f"  {label:9s}  single={ws:.4f}  multi={wm:.4f}  -> {verdict}")
    print(f"  acceptance single={out['single']['acceptance']:.3f} "
          f"multi={out['multi']['acceptance']:.3f}")
    print(f"  conv_ok single={out['single']['convergence_ok']} "
          f"multi={out['multi']['convergence_ok']}")
    print(f"  elapsed single={t_single:.0f}s multi={t_multi:.0f}s")

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"  wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
