"""Is the unsubtracted light actually the cause of the 41 mmag deficit?

The residual comparison found APEX leaving a median 23 % of a star's own flux
at the tightest blends where ALLSTAR leaves 1.5 %. Tempting as it is to declare
that the cause of the magnitude gap in the same bin, the two findings could be
unrelated symptoms: the measurement aperture (1 FWHM) overlaps the neighbour's
position in this bin by construction, so "light left at the star" may be a
neighbour ALLSTAR cleaned up and APEX did not, without the star's own fitted
flux being any worse for it. This script carries that assumption and tests the
link directly, three ways.

**Link test** — split the tightest-bin stars by how much light APEX left at
their position and compare the paired magnitude error (|APEX| − |ALLSTAR|) in
each class. If the deficit lives only in the big-leftover stars, the mechanism
is confirmed; if the clean-leftover stars show the same deficit, the residual
is a red herring for the magnitudes.

**Whose light** — measure the residual again in a half-FWHM aperture at the
star and at its nearest real neighbour separately, and check whether each
engine's output actually contains a fitted row for that neighbour and at what
flux relative to step 7's estimate. Distinguishes "the star was under-fitted"
from "the neighbour was under-fitted" from "the neighbour was never fitted".

**Variant sweep** — the same leftover statistic on the runs that already
exist: the wider-grouping run (every star solved jointly) and the raised
iteration-ceiling run. If those leave the same light behind, budget and
iteration count are excluded as the mechanism *for the pixels* too, not just
for the scatter.

Sign convention worth keeping in mind while reading the output: a star whose
own flux was under-fitted comes out too faint (delta_mag > 0); a star whose
fit absorbed an unsubtracted neighbour comes out too bright (delta_mag < 0).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

HERE = Path(__file__).absolute().parent
sys.path.insert(0, str(HERE))

from compare_recovery import (  # noqa: E402
    daophot_flux_adu, load_daophot, match_to_truth, remove_zeropoint,
)
from gap_decompose import ISOLATED, allstar_rows, apex_rows  # noqa: E402
from residual_compare import aperture_stats, sky_zeroed  # noqa: E402

TIGHT = "0.75-1.5 FWHM"


def boot_median(x: np.ndarray, n: int = 10000, seed: int = 0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    draws = np.median(rng.choice(x, size=(n, x.size), replace=True), axis=1)
    return float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", default=r"E:\APEX_validation\psf_engines\M13_ast2x200_moffat")
    ap.add_argument("--frame", default="pp_messier13-0005-B.fit")
    ap.add_argument("--apex-run", default="apexhybrid")
    ap.add_argument("--trials", type=int, default=2)
    ap.add_argument("--fwhm-px", type=float, default=7.052391)
    ap.add_argument("--output", default=str(HERE / "leftover_link_test.csv"))
    args = ap.parse_args()

    work = Path(args.work)
    fwhm = args.fwhm_px
    gates = {"snr": 3.0, "qfit": 3.0, "chi2": 25.0}
    truth_all = pd.read_csv(work / "truth.csv")

    scored_a, scored_d, per_star = [], [], []
    for trial in range(1, args.trials + 1):
        truth = truth_all[truth_all["trial"] == trial].reset_index(drop=True)
        tsv = (work / f"{args.apex_run}_trial{trial}" / "result" / "cmd_psf"
               / f"photometry_{args.frame}.tsv")
        dao_csv = work / f"daophot_trial{trial}.csv"
        scored_a.append(apex_rows(truth, tsv, gates, 1.5))
        scored_d.append(allstar_rows(truth, dao_csv, 25.0, 60.0, 1.5))

        # Residual images, blank-field level at zero, as in residual_compare.
        apex_res = sky_zeroed(
            work / f"{args.apex_run}_trial{trial}" / "result" / "cmd_psf"
            / f"residual_{args.frame}")
        dao_res = sky_zeroed(work / f"daophot_work_trial{trial}" / "sub.fits")

        # The nearest real star: second neighbour in the step-7 list, because
        # the implanted star itself is the first at zero distance.
        s7 = pd.read_csv(work / f"trial_{trial:04d}" / "result"
                         / "step7_forced_phot" / f"photometry_{args.frame}.tsv",
                         sep="\t")
        sx = pd.to_numeric(s7["x_fit"], errors="coerce").to_numpy()
        sy = pd.to_numeric(s7["y_fit"], errors="coerce").to_numpy()
        sf = pd.to_numeric(s7.get("flux_net_adu", s7.get("flux")),
                           errors="coerce").to_numpy()
        ok = np.isfinite(sx) & np.isfinite(sy) & np.isfinite(sf)
        sx, sy, sf = sx[ok], sy[ok], sf[ok]
        s7_tree = cKDTree(np.column_stack([sx, sy]))

        # Full (ungated) engine outputs, to ask whether the neighbour was
        # fitted at all and at what flux.
        af = pd.read_csv(tsv, sep="\t")
        for c in ("x_fit", "y_fit", "flux_psf_e", "gain_e_per_adu", "mag_psf"):
            af[c] = pd.to_numeric(af[c], errors="coerce")
        af = af[np.isfinite(af["mag_psf"])]
        apex_tree = cKDTree(af[["x_fit", "y_fit"]].to_numpy(float))
        apex_flux_adu = (af["flux_psf_e"] / af["gain_e_per_adu"]).to_numpy(float)

        dao = load_daophot(dao_csv)
        dao["mag"] = pd.to_numeric(dao["mag"], errors="coerce")
        dao = dao[np.isfinite(dao["mag"])]
        dao_tree = cKDTree(dao[["x", "y"]].to_numpy(float))
        dao_flux_adu = daophot_flux_adu(dao["mag"].to_numpy(float), 25.0, 60.0)

        for _, star in truth.iterrows():
            x, y = float(star["x_true"]), float(star["y_true"])
            flux = float(star["flux_realized_adu"])
            nb_d, nb_i = s7_tree.query([x, y], k=2)
            nx, ny_, nf = sx[nb_i[1]], sy[nb_i[1]], sf[nb_i[1]]

            a1, _ = aperture_stats(apex_res, x, y, 1.0 * fwhm)
            d1, _ = aperture_stats(dao_res, x, y, 1.0 * fwhm)
            a_half, _ = aperture_stats(apex_res, x, y, 0.5 * fwhm)
            a_nb, _ = aperture_stats(apex_res, nx, ny_, 0.5 * fwhm)
            d_half, _ = aperture_stats(dao_res, x, y, 0.5 * fwhm)
            d_nb, _ = aperture_stats(dao_res, nx, ny_, 0.5 * fwhm)

            ad, ai = apex_tree.query([nx, ny_], k=1)
            dd, di = dao_tree.query([nx, ny_], k=1)
            per_star.append({
                "trial": trial, "crowding_bin": star["crowding_bin"],
                "target_snr": star["target_snr"],
                "sep_px": float(nb_d[1]), "nb_flux_adu": float(nf),
                "star_flux_adu": flux,
                "apex_left_1fwhm": a1 / flux, "dao_left_1fwhm": d1 / flux,
                "apex_left_star_half": a_half / flux,
                "apex_left_nb_half": a_nb / max(nf, 1.0),
                "dao_left_star_half": d_half / flux,
                "dao_left_nb_half": d_nb / max(nf, 1.0),
                "apex_has_nb_row": bool(ad <= 1.0),
                "dao_has_nb_row": bool(dd <= 1.0),
                "apex_nb_flux_ratio": (float(apex_flux_adu[ai] / nf)
                                       if ad <= 1.0 and nf > 0 else np.nan),
                "dao_nb_flux_ratio": (float(dao_flux_adu[di] / nf)
                                      if dd <= 1.0 and nf > 0 else np.nan),
            })

    A = remove_zeropoint(pd.concat(scored_a, ignore_index=True),
                         min_snr=50.0, isolated_bins=ISOLATED)[0]
    D = remove_zeropoint(pd.concat(scored_d, ignore_index=True),
                         min_snr=50.0, isolated_bins=ISOLATED)[0]
    S = pd.DataFrame(per_star)
    S["apex_dmag"] = A["delta_mag"].to_numpy(float)
    S["dao_dmag"] = D["delta_mag"].to_numpy(float)
    S.to_csv(args.output, index=False)

    tight = S[S["crowding_bin"] == TIGHT].copy()
    both = np.isfinite(tight["apex_dmag"]) & np.isfinite(tight["dao_dmag"])
    print(f"최혼잡 구간 {len(tight)}개 · 그중 두 엔진 게이트 통과(짝지음 가능) {int(both.sum())}개\n")

    # ── 1. Link test ─────────────────────────────────────────────────────
    print("■ 연결 검사 — 남은 빛이 큰 별에서만 등급이 나빠지는가")
    print(f"{'남은 빛 (1·FWHM)':>18}{'N':>5}{'|APEX|-|ALLSTAR| 중앙':>22}"
          f"{'95% 구간 (mmag)':>20}{'APEX 부호중앙':>13}")
    print("-" * 82)
    classes = [("<= 0.1 (깨끗)", tight["apex_left_1fwhm"] <= 0.1),
               ("0.1-0.5", (tight["apex_left_1fwhm"] > 0.1)
                & (tight["apex_left_1fwhm"] <= 0.5)),
               ("> 0.5 (크게 남음)", tight["apex_left_1fwhm"] > 0.5)]
    for label, mask in classes:
        part = tight[mask & both]
        if len(part) < 4:
            print(f"{label:>18}{len(part):>5}   — 표본 부족")
            continue
        diff = (part["apex_dmag"].abs() - part["dao_dmag"].abs()).to_numpy(float)
        lo, hi = boot_median(diff)
        print(f"{label:>18}{len(part):>5}{np.median(diff) * 1000:>+18.1f} mmag"
              f"   [{lo * 1000:>+6.1f}, {hi * 1000:>+6.1f}]"
              f"{np.median(part['apex_dmag']) * 1000:>+11.1f}")
    ok = both.to_numpy()
    r1, p1 = spearmanr(tight["apex_left_1fwhm"][ok],
                       tight["apex_dmag"][ok].abs())
    r2, p2 = spearmanr(tight["apex_left_1fwhm"][ok],
                       (tight["apex_dmag"].abs() - tight["dao_dmag"].abs())[ok])
    print(f"\n  상관 (남은빛 vs |APEX 오차|)        r={r1:+.3f} p={p1:.4f}")
    print(f"  상관 (남은빛 vs |APEX|-|ALLSTAR|)  r={r2:+.3f} p={p2:.4f}")

    # ── 2. Whose light ───────────────────────────────────────────────────
    big = tight[tight["apex_left_1fwhm"] > 0.5]
    print(f"\n■ 누구의 빛인가 — 크게 남은 {len(big)}개, 반경 0.5·FWHM 로 분리")
    print(f"  별 자리 남은 빛   중앙 {np.median(big['apex_left_star_half']):+.3f} (별 밝기 대비)")
    print(f"  이웃 자리 남은 빛 중앙 {np.median(big['apex_left_nb_half']):+.3f} (이웃 밝기 대비)")
    print(f"  [ALLSTAR 같은 별] 별 {np.median(big['dao_left_star_half']):+.3f}"
          f" · 이웃 {np.median(big['dao_left_nb_half']):+.3f}")
    print(f"  APEX 출력에 이웃 행 있음   {int(big['apex_has_nb_row'].sum())}/{len(big)}"
          f" · 이웃 플럭스/step7 중앙 {np.nanmedian(big['apex_nb_flux_ratio']):.2f}")
    print(f"  ALLSTAR 출력에 이웃 행 있음 {int(big['dao_has_nb_row'].sum())}/{len(big)}"
          f" · 이웃 플럭스/step7 중앙 {np.nanmedian(big['dao_nb_flux_ratio']):.2f}")

    # ── 3. Variant sweep ─────────────────────────────────────────────────
    print("\n■ 변형 실행의 같은 잔차 — 예산·반복수가 픽셀에서도 무죄인지")
    variants = [("hybrid (기준)", args.apex_run), ("group (전부 동시적합)", "apexgroup"),
                ("iter (반복 30)", "apexiter"), ("fix (씨앗 문맥 수정)", "apexfix"),
                ("free (위치잠금 해제)", "apexfree")]
    tight_truth = []
    for trial in range(1, args.trials + 1):
        t = truth_all[truth_all["trial"] == trial].reset_index(drop=True)
        tight_truth.append(t[t["crowding_bin"] == TIGHT])
    for label, run in variants:
        vals = []
        for trial, tt in zip(range(1, args.trials + 1), tight_truth):
            res_path = (work / f"{run}_trial{trial}" / "result" / "cmd_psf"
                        / f"residual_{args.frame}")
            if not res_path.exists():
                continue
            res = sky_zeroed(res_path)
            for _, star in tt.iterrows():
                s, _ = aperture_stats(res, float(star["x_true"]),
                                      float(star["y_true"]), 1.0 * fwhm)
                vals.append(s / float(star["flux_realized_adu"]))
        if vals:
            print(f"  {label:>22}: 남은 빛 중앙 {np.median(vals):+.4f} (N={len(vals)})")
        else:
            print(f"  {label:>22}: 잔차 파일 없음")

    print(f"\nsaved -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
