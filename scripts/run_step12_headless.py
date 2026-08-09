"""Run Step-12 (isochrone MCMC fit) headless, without the GUI.

Reads the Step-10 wide CMD table and calls the SAME service the GUI Step-12
tab calls (``apex.analysis.cmd.isochrone_fit_service.fit_cluster_isochrone``),
so a headless fit and a GUI fit execute identical code. Qt is never imported.

    .venv-deploy/Scripts/python.exe scripts/run_step12_headless.py \
        --params E:/APEX_validation/psf_crossinstrument/m67_lco/parameters.toml \
        --colors B-V --mag V --age-gyr 1 12 --parallax-prior

Outputs (under <result_dir>/cmd_isochrone/):
    isochrone_fit.json   fitted parameters + uncertainties + run settings
    isochrone_cmd.png    observed CMD with the best-fit track
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


def _console(value) -> str:
    enc = getattr(sys.stdout, "encoding", None) or "utf-8"
    return str(value).encode(enc, errors="backslashreplace").decode(enc)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--params", default=str(REPO / "parameters.toml"))
    ap.add_argument("--colors", default="B-V",
                    help="colour(s), e.g. 'B-V' or 'g-r,r-i'")
    ap.add_argument("--mag", default="V", help="magnitude band for the CMD y-axis")
    ap.add_argument("--iso", default=None, help="isochrone grid (default: config/bundled)")
    ap.add_argument("--age-gyr", nargs=2, type=float, default=[0.5, 13.0],
                    metavar=("MIN", "MAX"))
    ap.add_argument("--mh-prior", nargs=2, type=float, default=None,
                    metavar=("MEAN", "SIGMA"))
    ap.add_argument("--ebv-prior", nargs=2, type=float, default=None,
                    metavar=("MEAN", "SIGMA"))
    ap.add_argument("--parallax-prior", action="store_true",
                    help="derive a (m-M)0 prior from Gaia parallaxes of members")
    ap.add_argument("--dm-prior", nargs=2, type=float, default=None,
                    metavar=("MEAN", "SIGMA"),
                    help="explicit (m-M)0 Gaussian prior — use for distant "
                         "clusters (globulars) where Gaia parallax is unusable")
    ap.add_argument("--mh-bounds", nargs=2, type=float, default=[-1.0, 0.5],
                    metavar=("MIN", "MAX"),
                    help="[M/H] bounds (widen for globulars, e.g. -2.2 0.5)")
    ap.add_argument("--no-membership", action="store_true")
    ap.add_argument("--max-stars", type=int, default=500)
    ap.add_argument("--walkers", type=int, default=32)
    ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--burn", type=int, default=600)
    ap.add_argument("--err-floor", type=float, default=0.02)
    # The binary and field fractions differ between clusters — an old open
    # cluster is not a sparse young one — so they must not stay compiled-in
    # constants. Defaults match IsochroneFitConfig; the cross-check against
    # ASteCA (validation/asteca_crosscheck/) is what made their fixedness
    # visible.
    ap.add_argument("--f-bin", type=float, default=None,
                    help="binary fraction in the CMD mixture (default 0.3)")
    ap.add_argument("--f-field", type=float, default=None,
                    help="field-contamination fraction (default 0.1)")
    ap.add_argument("--seed", type=int, default=2024)
    # The convergence question this posterior actually poses is "does the
    # answer depend on where the walkers started?", not "is the
    # autocorrelation time small". tau never settles here (44 at 400 steps,
    # 202 at 4,000) and split-R-hat is inflated because emcee's walkers are
    # not independent chains. Repeating the fit from different seeds answers
    # it directly: four seeds agreed to 0.02 dex in [M/H].
    ap.add_argument("--seed-check", type=int, default=0, metavar="N",
                    help="repeat the fit with N extra seeds and report the "
                         "spread of the medians (0 = off)")
    ap.add_argument("--mag-max", type=float, default=None,
                    help="drop stars fainter than this in the magnitude band "
                         "(diagnostic: isolates the bright end where the "
                         "Gaia-transformation colour drift is flat)")
    args = ap.parse_args()

    import numpy as np
    import pandas as pd

    from apex.config.parameters_cmd import read_params
    from apex.analysis.cmd.isochrone_data import EXTINCTION_R, default_iso_file
    from apex.analysis.cmd.isochrone_fit_service import (
        IsochroneFitConfig,
        fit_cluster_isochrone,
    )
    from apex.utils.step_paths_cmd import step10_zp_dir, step12_iso_dir

    params = read_params(args.params)
    result_dir = Path(params.P.result_dir)

    wide = step10_zp_dir(result_dir) / "median_by_ID_filter_wide_cmd.csv"
    if not wide.exists():
        print(f"[error] Step-10 wide CMD table not found: {wide}")
        return 1
    df = pd.read_csv(wide)
    print(f"wide CMD: {len(df)} rows | columns: "
          f"{[c for c in df.columns if c.startswith('mag_std')]}")

    if args.mag_max is not None:
        col = f"mag_std_{args.mag}"
        before = len(df)
        v = pd.to_numeric(df[col], errors="coerce")
        df = df[v.notna() & (v <= float(args.mag_max))].copy()
        print(f"mag cut: {col} <= {args.mag_max} -> {len(df)}/{before} rows")

    colors = [tuple(part.split("-", 1)) for part in args.colors.split(",")]
    for b1, b2 in colors:
        for b in (b1, b2):
            if f"mag_std_{b}" not in df.columns:
                print(f"[error] column mag_std_{b} missing for colour {b1}-{b2}")
                return 1

    b1, b2 = colors[0]
    r_color = EXTINCTION_R.get(b1, 3.303) - EXTINCTION_R.get(b2, 2.285)
    a_min = max(args.age_gyr[0], 1e-3)
    a_max = max(args.age_gyr[1], a_min * 1.01)
    age_bounds = (9.0 + float(np.log10(a_min)), 9.0 + float(np.log10(a_max)))

    ecolor_prior = None
    if args.ebv_prior:
        ecolor_prior = (args.ebv_prior[0] * r_color,
                        max(args.ebv_prior[1] * abs(r_color), 1e-3))

    # 그리드는 색에서 자동 선택한다(Johnson vs SDSS). config 의
    # isochrone.file_path 를 쓰지 않는 이유: 설정이 다른 대상에서 복사되면
    # 측광 시스템이 안 맞는 그리드를 조용히 쓰게 된다.
    iso_file = Path(args.iso) if args.iso else default_iso_file(colors[0])
    if not iso_file.exists():
        print(f"[error] isochrone grid not found: {iso_file}")
        return 1

    cfg = IsochroneFitConfig(
        colors=colors, mag_band=args.mag, iso_file=iso_file,
        age_bounds=age_bounds,
        mh_bounds=(float(args.mh_bounds[0]), float(args.mh_bounds[1])),
        dm_bounds=(5.0, 18.0), ecolor_bounds=(0.0, 1.0),
        mh_prior=tuple(args.mh_prior) if args.mh_prior else None,
        ecolor_prior=ecolor_prior,
        dm_prior=tuple(args.dm_prior) if args.dm_prior else None,
        parallax_distance_prior=bool(args.parallax_prior),
        use_membership=not args.no_membership,
        max_stars=args.max_stars, n_walkers=args.walkers,
        n_steps=args.steps, n_burn=args.burn,
        err_floor=args.err_floor, seed=args.seed,
        **({} if args.f_bin is None else {"f_bin": float(args.f_bin)}),
        **({} if args.f_field is None else {"f_field": float(args.f_field)}),
    )
    print(f"fit: colors={colors} mag={args.mag} age={a_min}-{a_max} Gyr "
          f"iso={iso_file.name} membership={cfg.use_membership} "
          f"parallax_prior={cfg.parallax_distance_prior}")

    t0 = time.perf_counter()
    last = [-1.0]

    def _progress(frac, msg):
        if frac - last[0] >= 0.1 or frac >= 1.0:
            last[0] = frac
            print(f"  [{frac*100:3.0f}%] {_console(msg)}", flush=True)

    out = fit_cluster_isochrone(df, cfg, make_figures=False, progress_cb=_progress)
    elapsed = time.perf_counter() - t0

    seed_check = None
    if args.seed_check > 0:
        import dataclasses

        keys = ("age_gyr", "metallicity", "e_bv", "distance_mod")
        runs = [{k: out.summary[k][1] for k in keys if k in out.summary}]
        for i in range(int(args.seed_check)):
            alt = dataclasses.replace(cfg, seed=int(args.seed) + 1000 * (i + 1))
            print(f"  [seed-check {i+1}/{args.seed_check}] seed={alt.seed}",
                  flush=True)
            r = fit_cluster_isochrone(df, alt, make_figures=False)
            runs.append({k: r.summary[k][1] for k in keys if k in r.summary})
        seed_check = {"n_seeds": len(runs), "runs": runs, "spread": {}}
        print("  seed-check spread (max - min over seeds):")
        for k in keys:
            vals = [r[k] for r in runs if k in r]
            if not vals:
                continue
            spread = float(max(vals) - min(vals))
            seed_check["spread"][k] = spread
            print(f"    {k:14s} {spread:.4f}   "
                  f"[{min(vals):.4f}, {max(vals):.4f}]", flush=True)

    out_dir = step12_iso_dir(result_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    payload = {
        "summary": out.summary,
        "n_stars": int(out.n_stars),
        "member_meta": out.member_meta,
        "warnings": list(out.warnings),
        "settings": {
            "colors": [f"{a}-{b}" for a, b in colors],
            "mag_band": args.mag,
            "age_gyr": [a_min, a_max],
            "iso_file": str(iso_file),
            "mh_prior": args.mh_prior,
            "ebv_prior": args.ebv_prior,
            "parallax_prior": bool(args.parallax_prior),
            "membership": cfg.use_membership,
            "mag_max": args.mag_max,
            "walkers": args.walkers, "steps": args.steps, "burn": args.burn,
            "err_floor": args.err_floor, "seed": args.seed,
            "f_bin": cfg.f_bin, "f_field": cfg.f_field,
        },
        "elapsed_s": elapsed,
        "seed_check": seed_check,
        "wide_table": str(wide),
    }
    (out_dir / "isochrone_fit.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
    )

    # CMD figure — matplotlib only (no Qt), built from the returned plot arrays.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6.2, 7.0))
        ax.plot(out.obs_color, out.obs_mag, ".", ms=3, alpha=0.55,
                color="#1f77b4", label=f"stars (N={out.n_stars})")
        ax.plot(out.iso_color, out.iso_mag, "-", lw=1.6, color="crimson",
                label="best-fit isochrone")
        ax.invert_yaxis()
        ax.set_xlabel(out.color_label or f"{b1} - {b2}")
        ax.set_ylabel(out.mag_label or args.mag)
        s = out.summary
        bits = []
        for k, lab in (("log_age", "log age"), ("age_gyr", "age [Gyr]"),
                       ("mh", "[M/H]"), ("dm", "(m-M)0"), ("ebv", "E(B-V)")):
            if k in s and isinstance(s[k], (int, float)):
                bits.append(f"{lab}={s[k]:.3f}")
        ax.set_title(" · ".join(bits) if bits else "isochrone fit", fontsize=9)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out_dir / "isochrone_cmd.png", dpi=140)
        plt.close(fig)
    except Exception as exc:
        print(f"[warn] figure failed: {_console(exc)}")

    print(f"\n[done] {elapsed:.1f}s  stars={out.n_stars}")
    for k, v in sorted(out.summary.items()):
        if isinstance(v, (int, float, str)):
            print(f"  {k:22s} = {v}")
    for w in out.warnings:
        print(f"  [warn] {_console(w)}")
    print(f"  -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
