"""Take the remaining APEX-vs-ALLSTAR gap apart, one factor at a time.

Three explanations for APEX's larger scatter have been tested and killed: the
PSF model (ePSF, analytic Moffat and the DAOPHOT-style hybrid all land within
0.002 mag of each other), the interpolation order used to sample the model, and
simultaneous group fitting (six times more stars solved jointly changed
nothing). APEX sits at 0.039 mag where ALLSTAR sits at 0.029 on the same stars,
and nothing about how the model is built or fitted accounts for it.

Before running more variants it is worth asking a cheaper question of the data
already on disk: **does each engine know how noisy it is?** Both report a
per-star magnitude error, so the normalised residual

    chi = (measured - true) / sigma_reported

has a robust scatter near 1 whenever the error model is complete. The answer
splits the search:

* both near 1 — the noise models are honest and APEX is genuinely extracting
  less information per star; the cause is upstream of the fit (how many pixels
  carry signal, what background was removed);
* APEX above 1 — APEX carries an error its own model does not know about, which
  is a bug-shaped problem rather than an information one;
* APEX near 1 but its *reported* sigma is larger — APEX is correctly reporting
  that it took on more noise, and the question becomes where that noise entered.

The same normalisation is then cut by brightness, by crowding, by local
background and by position in the frame, because a gap that is flat in every cut
means something global (a constant added variance) while a gap concentrated in
one cut names its own cause.

Nothing here re-runs an engine; it reads the products the earlier comparisons
already wrote.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).absolute().parent
sys.path.insert(0, str(HERE))

from compare_recovery import (  # noqa: E402
    load_daophot,
    daophot_flux_adu, match_to_truth, remove_zeropoint, robust_scatter,
)
from compare_three_engines import QUALITY  # noqa: E402

ISOLATED = ("3-6 FWHM", "6-inf FWHM")


def apex_rows(truth: pd.DataFrame, tsv: Path, gates: dict,
              radius_px: float) -> pd.DataFrame:
    """APEX measurements matched to truth, carrying its reported error."""
    d = pd.read_csv(tsv, sep="\t")
    for c in ("x_fit", "y_fit", "mag_psf", "mag_psf_err", "snr_psf",
              "n_pixels_fit", "psf_nea_px", "reduced_chi2", *QUALITY):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    keep = ((d["flags_psf"] == 0) & np.isfinite(d["mag_psf"])
            & (d["snr_psf"] >= gates["snr"]) & (d["qfit"] <= gates["qfit"])
            & (d["reduced_chi2"] <= gates["chi2"]))
    d = d[keep]
    table = pd.DataFrame({"x": d["x_fit"].to_numpy(), "y": d["y_fit"].to_numpy(),
                          "mag": d["mag_psf"].to_numpy()})
    idx = match_to_truth(truth, table, radius_px)
    out = truth.copy()
    take = lambda col: np.where(idx >= 0, d[col].to_numpy()[idx], np.nan)  # noqa: E731
    out["delta_mag"] = (np.where(idx >= 0, table["mag"].to_numpy()[idx], np.nan)
                        + 2.5 * np.log10(truth["flux_realized_adu"].to_numpy(float)))
    out["sigma_reported"] = take("mag_psf_err")
    out["n_pixels_fit"] = take("n_pixels_fit")
    out["fit_chi2"] = take("reduced_chi2")
    out["sky_fitted"] = np.nan  # APEX subtracts the background before fitting
    return out


def allstar_rows(truth: pd.DataFrame, csv: Path, zmag: float, itime: float,
                 radius_px: float) -> pd.DataFrame:
    dao = load_daophot(csv)
    for c in ("mag", "merr", "msky", "chi", "sharpness"):
        if c in dao.columns:
            dao[c] = pd.to_numeric(dao[c], errors="coerce")
    dao = dao[np.isfinite(dao["mag"])]
    idx = match_to_truth(truth, dao, radius_px)
    measured = np.where(idx >= 0,
                        daophot_flux_adu(dao["mag"].to_numpy(float)[idx],
                                         zmag, itime), np.nan)
    out = truth.copy()
    out["delta_mag"] = -2.5 * np.log10(
        measured / truth["flux_realized_adu"].to_numpy(float))
    take = lambda col: np.where(idx >= 0, dao[col].to_numpy()[idx], np.nan)  # noqa: E731
    out["sigma_reported"] = take("merr")
    out["n_pixels_fit"] = np.nan  # not reported by ALLSTAR
    out["fit_chi2"] = take("chi") ** 2
    out["sky_fitted"] = take("msky")
    return out


def cut_by(frame: pd.DataFrame, column: str, label: str,
           edges: np.ndarray | None = None, n_bins: int = 4) -> pd.DataFrame:
    """Group a frame into bins of one variable and report scatter and chi."""
    v = pd.to_numeric(frame[column], errors="coerce")
    ok = np.isfinite(v) & np.isfinite(frame["delta_mag"])
    if ok.sum() < 12:
        return pd.DataFrame()
    if edges is None:
        edges = np.nanpercentile(v[ok], np.linspace(0, 100, n_bins + 1))
        edges[-1] += abs(edges[-1]) * 1e-6 + 1e-9
    which = np.digitize(v, edges) - 1
    rows = []
    for k in range(len(edges) - 1):
        part = frame[ok & (which == k)]
        if len(part) < 6:
            continue
        rows.append({
            "cut": label,
            "bin": f"{edges[k]:.4g}–{edges[k + 1]:.4g}",
            "n": len(part),
            "scatter": robust_scatter(part["delta_mag"].to_numpy(float)),
            "sigma_med": float(np.nanmedian(part["sigma_reported"])),
            "chi": robust_scatter(
                (part["delta_mag"] / part["sigma_reported"]).to_numpy(float)),
        })
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", default=r"E:\APEX_validation\psf_engines\M13_ast2x200_moffat")
    ap.add_argument("--frame", default="pp_messier13-0005-B.fit")
    ap.add_argument("--apex-run", default="apexhybrid")
    ap.add_argument("--trials", type=int, default=2)
    ap.add_argument("--zmag", type=float, default=25.0)
    ap.add_argument("--itime", type=float, default=60.0)
    ap.add_argument("--match-radius-px", type=float, default=1.5)
    ap.add_argument("--postfit-snr-min", type=float, default=3.0)
    ap.add_argument("--postfit-qfit-max", type=float, default=3.0)
    ap.add_argument("--postfit-chi2-max", type=float, default=25.0)
    ap.add_argument("--output", default=str(HERE / "gap_decompose.csv"))
    args = ap.parse_args()

    work = Path(args.work)
    gates = {"snr": args.postfit_snr_min, "qfit": args.postfit_qfit_max,
             "chi2": args.postfit_chi2_max}
    truth_all = pd.read_csv(work / "truth.csv")

    per_engine: dict[str, list[pd.DataFrame]] = {"APEX": [], "ALLSTAR": []}
    for trial in range(1, args.trials + 1):
        truth = truth_all[truth_all["trial"] == trial].reset_index(drop=True)
        tsv = (work / f"{args.apex_run}_trial{trial}" / "result" / "cmd_psf"
               / f"photometry_{args.frame}.tsv")
        csv = work / f"daophot_trial{trial}.csv"
        if truth.empty or not tsv.exists() or not csv.exists():
            print(f"[건너뜀] trial {trial}")
            continue
        per_engine["APEX"].append(
            apex_rows(truth, tsv, gates, args.match_radius_px))
        per_engine["ALLSTAR"].append(
            allstar_rows(truth, csv, args.zmag, args.itime, args.match_radius_px))

    frames: dict[str, pd.DataFrame] = {}
    for name, parts in per_engine.items():
        if not parts:
            raise SystemExit("채점할 회차가 없다")
        frame = pd.concat(parts, ignore_index=True)
        frame, zp = remove_zeropoint(frame, min_snr=50.0, isolated_bins=ISOLATED)
        frames[name] = frame
        print(f"  {name:>8} 영점 {zp['zeropoint_offset_mag']:+.3f} mag "
              f"({zp['n_offset_stars']}개 밝은 고립성)")

    # Only stars both engines kept, so the comparison is not a selection effect.
    both = np.ones(len(frames["APEX"]), dtype=bool)
    for frame in frames.values():
        both &= np.isfinite(frame["delta_mag"]).to_numpy()
    for name in frames:
        frames[name] = frames[name][both].reset_index(drop=True)
    print(f"\n두 엔진이 모두 살린 별 {int(both.sum())}개로 분해한다.\n")

    print(f"{'엔진':>10}{'실제 산포':>12}{'보고 오차 중앙':>16}{'chi 산포':>12}"
          f"{'판정':>28}")
    print("-" * 80)
    summary = []
    for name, frame in frames.items():
        s = robust_scatter(frame["delta_mag"].to_numpy(float))
        sig = float(np.nanmedian(frame["sigma_reported"]))
        chi = robust_scatter(
            (frame["delta_mag"] / frame["sigma_reported"]).to_numpy(float))
        verdict = ("오차모형 정직" if 0.75 <= chi <= 1.35
                   else ("과소보고 — 모르는 오차원 있음" if chi > 1.35
                         else "과대보고 — 오차를 부풀림"))
        print(f"{name:>10}{s:>12.4f}{sig:>16.4f}{chi:>12.2f}{verdict:>28}")
        summary.append({"engine": name, "scatter": s, "sigma_median": sig,
                        "chi_scatter": chi})

    cuts = [("target_snr", "밝기 (주입 SNR)", None),
            ("nearest_real_sep_fwhm", "혼잡 (최근접 이웃)", None),
            ("preinjection_local_rms_adu", "국소 배경 잡음", None),
            ("radius_arcmin", "프레임 내 위치", None)]
    parts = []
    for col, label, edges in cuts:
        if col not in frames["APEX"].columns:
            continue
        print(f"\n■ {label}")
        print(f"{'구간':>18}{'N':>5}"
              f"{'APEX 산포':>11}{'ALLSTAR 산포':>13}{'비':>7}"
              f"{'APEX chi':>10}{'ALLSTAR chi':>13}")
        a = cut_by(frames["APEX"], col, label, edges)
        b = cut_by(frames["ALLSTAR"], col, label, edges)
        if a.empty or b.empty:
            print("   (표본 부족)")
            continue
        merged = a.merge(b, on=["cut", "bin"], suffixes=("_apex", "_dao"))
        merged["ratio"] = merged["scatter_apex"] / merged["scatter_dao"]
        for _, r in merged.iterrows():
            print(f"{r['bin']:>18}{r['n_apex']:>5}"
                  f"{r['scatter_apex']:>11.4f}{r['scatter_dao']:>13.4f}"
                  f"{r['ratio']:>7.2f}{r['chi_apex']:>10.2f}{r['chi_dao']:>13.2f}")
        parts.append(merged)

    if parts:
        out = pd.concat(parts, ignore_index=True)
        out.to_csv(args.output, index=False)
        print(f"\nsaved -> {args.output}")

    paired = paired_effect(frames["APEX"], frames["ALLSTAR"])
    paired.to_csv(Path(args.output).with_name("gap_paired_effect.csv"), index=False)
    pd.DataFrame(summary).to_csv(
        Path(args.output).with_name("gap_decompose_summary.csv"), index=False)
    return 0


def paired_effect(apex: pd.DataFrame, dao: pd.DataFrame) -> pd.DataFrame:
    """How much better is ALLSTAR on the very same star?

    A ratio of two scatters is a weak statistic here: it throws away the pairing
    and, with fewer than a couple hundred stars, its confidence interval reaches
    from "no difference" to "half again as noisy". Both engines measured the
    same star on the same pixels, so the sharper question is whether *this
    star's* error is smaller in one engine, which cancels everything the two
    share — most importantly the leftover light already sitting at each
    injection site, which drives both engines equally.

    Reported per crowding bin: the median of |APEX error| - |ALLSTAR error| in
    mmag with a bootstrapped interval, a Wilcoxon signed-rank p, and that p
    multiplied by the number of bins, because four bins invite one accident.
    """
    from scipy.stats import wilcoxon

    d = (apex["delta_mag"].abs() - dao["delta_mag"].abs()).to_numpy(float)
    bins = [b for b in apex["crowding_bin"].unique() if (apex["crowding_bin"] == b).sum() >= 10]

    def boot_median(x: np.ndarray, n: int = 20000, seed: int = 3) -> tuple[float, float]:
        rng = np.random.default_rng(seed)
        draws = np.median(rng.choice(x, size=(n, x.size), replace=True), axis=1)
        return tuple(np.percentile(draws, [2.5, 97.5]))

    print("\n■ 같은 별에서 누가 더 정확했나 "
          "(양수 = ALLSTAR 가 그 별을 더 정확히 쟀다)")
    print(f"{'혼잡 구간':>16}{'N':>5}{'중앙 차 (mmag)':>16}{'95% 구간':>22}"
          f"{'p':>9}{'구간수 보정 p':>14}")
    print("-" * 84)
    rows = []
    for label in list(bins) + ["전체"]:
        part = d if label == "전체" else d[(apex["crowding_bin"] == label).to_numpy()]
        if part.size < 10:
            continue
        lo, hi = boot_median(part)
        p = float(wilcoxon(part).pvalue)
        adj = min(1.0, p * len(bins)) if label != "전체" else float("nan")
        print(f"{label:>16}{part.size:>5}{np.median(part) * 1000:>+16.1f}"
              f"   [{lo * 1000:>+6.1f}, {hi * 1000:>+6.1f}]{p:>9.3f}"
              f"{adj if np.isfinite(adj) else float('nan'):>14.3f}")
        rows.append({"crowding_bin": label, "n": int(part.size),
                     "median_diff_mmag": float(np.median(part)) * 1000,
                     "ci_lo_mmag": lo * 1000, "ci_hi_mmag": hi * 1000,
                     "wilcoxon_p": p, "p_adjusted": adj})
    print("\n가장 한산한 구간에서 0 에 붙고 가장 겹친 구간에서만 벌어지면,\n"
          "격차는 엔진의 전반적 정밀도가 아니라 겹친 별 처리에 있다.")
    return pd.DataFrame(rows)


if __name__ == "__main__":
    raise SystemExit(main())
