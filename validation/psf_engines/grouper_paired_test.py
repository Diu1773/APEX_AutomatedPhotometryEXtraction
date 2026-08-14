"""Did enlarging the simultaneous group actually help, or is 17 stars just noisy?

The three-engine table showed the most crowded bin improve from 0.034 to 0.024
mag when the grouper was widened, while the isolated bins did not move at all.
That localisation is suggestive, but each bin holds only 17-31 stars in one
trial, and a scatter estimated from 17 numbers carries roughly 18 % uncertainty
on its own — enough to manufacture a 29 % "improvement" out of nothing.

The two runs are not independent samples, though: the same implanted stars were
measured twice, changing only the grouper. That makes the comparison paired, so
the question can be asked star by star — did *this* star's error shrink? — which
removes the star-to-star brightness and position scatter that dominates the
unpaired estimate.

Reported per crowding bin:

* the two scatters, as in the engine table;
* how many stars improved versus worsened, with an exact two-sided sign test;
* the median change in |error|, bootstrapped for a confidence interval.

The sign test assumes nothing about the error distribution, which matters here
because photometric errors in blended stars are not Gaussian.

Result (2026-08-14, two trials, M13 B frame): the widening reached 6× more stars
— 95 solved jointly before, 564 after — and changed nothing. Overall scatter
0.039 both ways; of the stars whose value actually moved, 85 improved and 77
worsened, p = 0.58, median change in |error| −0.0009 mag with a 95 % interval
straddling zero. The trial-1 improvement in the most crowded bin did not
survive trial 2: only 2 of that bin's 23 stars had moved at all, which is enough
to shift a robust scatter computed from seventeen numbers. That is the reason
this script exists — the engine table alone would have reported a win.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).absolute().parent
sys.path.insert(0, str(HERE))

from compare_recovery import remove_zeropoint, robust_scatter  # noqa: E402
from compare_three_engines import apex_moffat_scored  # noqa: E402

ISOLATED = ("3-6 FWHM", "6-inf FWHM")


def load(work: Path, run: str, trial: int, frame: str, truth: pd.DataFrame,
         gates: dict, radius: float) -> pd.DataFrame | None:
    tsv = work / f"{run}_trial{trial}" / "result" / "cmd_psf" / f"photometry_{frame}.tsv"
    if not tsv.exists():
        return None
    scored = apex_moffat_scored(truth, tsv, gates, radius)
    scored, _ = remove_zeropoint(scored, min_snr=50.0, isolated_bins=ISOLATED)
    return scored


def bootstrap_median(x: np.ndarray, n: int = 10000, seed: int = 0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    draws = rng.choice(x, size=(n, x.size), replace=True)
    lo, hi = np.percentile(np.median(draws, axis=1), [2.5, 97.5])
    return float(lo), float(hi)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", default=r"E:\APEX_validation\psf_engines\M13_ast2x200_moffat")
    ap.add_argument("--frame", default="pp_messier13-0005-B.fit")
    ap.add_argument("--trials", type=int, default=2)
    ap.add_argument("--baseline", default="apexhybrid")
    ap.add_argument("--variant", default="apexgroup")
    ap.add_argument("--match-radius-px", type=float, default=1.5)
    ap.add_argument("--postfit-snr-min", type=float, default=3.0)
    ap.add_argument("--postfit-qfit-max", type=float, default=3.0)
    ap.add_argument("--postfit-chi2-max", type=float, default=25.0)
    ap.add_argument("--output", default=str(HERE / "grouper_paired_test.csv"))
    args = ap.parse_args()

    work = Path(args.work)
    gates = {"snr": args.postfit_snr_min, "qfit": args.postfit_qfit_max,
             "chi2": args.postfit_chi2_max}
    truth_all = pd.read_csv(work / "truth.csv")

    parts = []
    used = []
    for trial in range(1, args.trials + 1):
        truth = truth_all[truth_all["trial"] == trial].reset_index(drop=True)
        if truth.empty:
            continue
        base = load(work, args.baseline, trial, args.frame, truth, gates,
                    args.match_radius_px)
        var = load(work, args.variant, trial, args.frame, truth, gates,
                   args.match_radius_px)
        if base is None or var is None:
            print(f"[건너뜀] trial {trial}: 두 실행 중 하나가 없다")
            continue
        used.append(trial)
        both = np.isfinite(base["delta_mag"]) & np.isfinite(var["delta_mag"])
        parts.append(pd.DataFrame({
            "trial": trial,
            "crowding_bin": truth["crowding_bin"],
            "base": base["delta_mag"],
            "var": var["delta_mag"],
        })[both.to_numpy()])

    if not parts:
        raise SystemExit("짝지을 회차가 없다")
    d = pd.concat(parts, ignore_index=True)
    d["abs_base"] = d["base"].abs()
    d["abs_var"] = d["var"].abs()
    d["change"] = d["abs_var"] - d["abs_base"]

    print(f"\n회차 {used} · 짝지은 별 {len(d)}개 "
          f"({args.baseline} → {args.variant}, 그룹 묶기만 다름)")
    head = (f"\n{'혼잡 구간':>16}{'건드림/전체':>12}{'기준 산포':>10}{'변경 산포':>10}"
            f"{'좋아짐':>8}{'나빠짐':>8}{'부호검정 p':>12}"
            f"{'|오차| 중앙변화 (95% 구간)':>30}")
    print(head)
    print("-" * 104)

    # Most stars come out bit-identical: widening the radius from 1.5 to 2.5
    # FWHM only pulls in neighbours that sit in that annulus, and a star with no
    # neighbour there is fitted exactly as before. Counting those untouched
    # stars as "no worse" would dilute any real effect into significance noise,
    # so the sign test runs on the stars the change actually reached.
    d["touched"] = d["var"].sub(d["base"]).abs() > 1e-4

    rows = []
    for label in list(d["crowding_bin"].unique()) + ["전체"]:
        part = d if label == "전체" else d[d["crowding_bin"] == label]
        n_all = len(part)
        part = part[part["touched"]]
        if len(part) < 3:
            print(f"{label:>16}{n_all:>5}  건드린 별 {len(part)}개 — 검정 불가")
            continue
        better = int((part["change"] < 0).sum())
        worse = int((part["change"] > 0).sum())
        # Exact binomial on the stars that moved at all; ties carry no sign.
        p = (stats.binomtest(better, better + worse, 0.5).pvalue
             if better + worse else float("nan"))
        med = float(np.median(part["change"]))
        lo, hi = bootstrap_median(part["change"].to_numpy())
        s_b = robust_scatter(part["base"].to_numpy(float))
        s_v = robust_scatter(part["var"].to_numpy(float))
        print(f"{label:>16}{len(part):>5}{s_b:>10.3f}{s_v:>10.3f}"
              f"{better:>8}{worse:>8}{p:>12.3f}"
              f"{med:>+16.4f}  [{lo:+.4f}, {hi:+.4f}]")
        rows.append({"crowding_bin": label, "n_pairs": len(part),
                     "scatter_baseline": s_b, "scatter_variant": s_v,
                     "n_better": better, "n_worse": worse, "sign_test_p": p,
                     "median_abs_change": med, "ci_lo": lo, "ci_hi": hi})

    print("\n음수 = 그룹 묶기가 오차를 줄였다. p 는 '좋아짐/나빠짐이 반반' 가설의 확률.")
    pd.DataFrame(rows).to_csv(args.output, index=False)
    print(f"saved -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
