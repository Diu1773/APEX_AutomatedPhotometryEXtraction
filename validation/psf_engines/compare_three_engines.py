"""Score three PSF engines on one set of implanted stars.

The two-engine comparison left an open question: APEX's larger scatter came
from the *model*, not from how it picks reference stars — the empirical ePSF is
flexible enough to absorb the reference stars' pixel noise, while DAOPHOT's
analytic Moffat is smooth by construction. APEX has a Moffat mode of its own
(`psf_build_mode`), unlocked 2026-08-14. If the model is really the cause, then
APEX-Moffat should sit closer to DAOPHOT than APEX-ePSF does.

Three measurements of the same implanted stars, in the same frames:

* **APEX ePSF** — the production default, read from the injection benchmark's
  own recovery table.
* **APEX Moffat** — the same step 8 with `build_mode=moffat`, re-run on the
  same injected frames.
* **IRAF ALLSTAR** — DAOPHOT, from the earlier run.

Scoring reuses `compare_recovery`'s functions unchanged, so all three go
through one definition of matching, zero-point removal and quality gating.
APEX-Moffat and ALLSTAR are both scored the "flat table" way (match to truth by
position); APEX-ePSF keeps its benchmark's own delta_mag, exactly as the
two-engine comparison did.

The stars are implanted with an independently fitted Moffat, which is neither
engine's model — established 2026-08-14 as the neutral injection.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).absolute().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))

from compare_recovery import (  # noqa: E402
    daophot_flux_adu, match_to_truth, remove_zeropoint, summarise,
)

QUALITY = ("flags_psf", "snr_psf", "qfit", "reduced_chi2")


def apex_epsf_scored(truth: pd.DataFrame, recovery: pd.DataFrame,
                     gates: dict) -> pd.DataFrame:
    """APEX's own recovery table, gated by its own post-fit policy."""
    key = ["x_true", "y_true"]
    carry = key + ["delta_mag"] + [c for c in QUALITY if c in recovery.columns]
    out = truth[key + ["target_snr", "crowding_bin"]].merge(
        recovery[carry].drop_duplicates(subset=key), on=key, how="left")
    kept = (
        (pd.to_numeric(out["flags_psf"], errors="coerce") == 0)
        & (pd.to_numeric(out["snr_psf"], errors="coerce") >= gates["snr"])
        & (pd.to_numeric(out["qfit"], errors="coerce") <= gates["qfit"])
        & (pd.to_numeric(out["reduced_chi2"], errors="coerce") <= gates["chi2"])
    )
    out.loc[~kept, "delta_mag"] = np.nan
    return out


def apex_moffat_scored(truth: pd.DataFrame, tsv: Path, gates: dict,
                       radius_px: float) -> pd.DataFrame:
    """APEX step 8 run with build_mode=moffat, matched to truth by position.

    Scored the same way ALLSTAR is: the engine's own quality flags decide what
    it keeps, then position matching decides which implanted star it is.
    """
    d = pd.read_csv(tsv, sep="\t")
    for c in ("x_fit", "y_fit", "mag_psf", "mag_psf_err", *QUALITY):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    keep = (
        (d["flags_psf"] == 0) & np.isfinite(d["mag_psf"])
        & (d["snr_psf"] >= gates["snr"]) & (d["qfit"] <= gates["qfit"])
        & (d["reduced_chi2"] <= gates["chi2"])
    )
    d = d[keep]
    table = pd.DataFrame({"x": d["x_fit"], "y": d["y_fit"],
                          "mag": d["mag_psf"], "merr": d["mag_psf_err"]})
    idx = match_to_truth(truth, table, radius_px)
    # Instrumental magnitudes already; flux conversion only matters for the
    # IRAF scale, so delta_mag is a magnitude difference directly.
    measured = np.where(idx >= 0, table["mag"].to_numpy()[idx], np.nan)
    true_mag = -2.5 * np.log10(truth["flux_realized_adu"].to_numpy(float))
    out = truth[["x_true", "y_true", "target_snr", "crowding_bin"]].copy()
    out["delta_mag"] = measured - true_mag
    return out


def allstar_scored(truth: pd.DataFrame, csv: Path, zmag: float, itime: float,
                   radius_px: float) -> pd.DataFrame:
    dao = pd.read_csv(csv)
    dao["mag"] = pd.to_numeric(dao["mag"], errors="coerce")
    dao = dao[np.isfinite(dao["mag"])]
    idx = match_to_truth(truth, dao, radius_px)
    measured = np.where(idx >= 0,
                        daophot_flux_adu(dao["mag"].to_numpy(float)[idx],
                                         zmag, itime), np.nan)
    out = truth[["x_true", "y_true", "target_snr", "crowding_bin"]].copy()
    out["delta_mag"] = -2.5 * np.log10(
        measured / truth["flux_realized_adu"].to_numpy(float))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", default=r"E:\APEX_validation\psf_engines\M13_ast2x200_moffat")
    ap.add_argument("--frame", default="pp_messier13-0005-B.fit")
    ap.add_argument("--trials", type=int, default=2)
    ap.add_argument("--zmag", type=float, default=25.0)
    ap.add_argument("--itime", type=float, default=60.0)
    ap.add_argument("--match-radius-px", type=float, default=1.5)
    ap.add_argument("--postfit-snr-min", type=float, default=3.0)
    ap.add_argument("--postfit-qfit-max", type=float, default=3.0)
    ap.add_argument("--postfit-chi2-max", type=float, default=25.0)
    ap.add_argument("--offset-min-snr", type=float, default=50.0)
    ap.add_argument("--isolated-bins", default="3-6 FWHM|6-inf FWHM")
    ap.add_argument("--output", default=str(HERE / "recovery_three_engines.csv"))
    args = ap.parse_args()

    work = Path(args.work)
    gates = {"snr": args.postfit_snr_min, "qfit": args.postfit_qfit_max,
             "chi2": args.postfit_chi2_max}
    isolated = tuple(args.isolated_bins.split("|"))
    truth_all = pd.read_csv(work / "truth.csv")
    recovery_all = pd.read_csv(work / "recovery.csv")

    parts: list[pd.DataFrame] = []
    for trial in range(1, args.trials + 1):
        truth = truth_all[truth_all["trial"] == trial].reset_index(drop=True)
        recovery = recovery_all[recovery_all["trial"] == trial]
        moffat_tsv = (work / f"apexmoffat_trial{trial}" / "result" / "cmd_psf"
                      / f"photometry_{args.frame}.tsv")
        dao_csv = work / f"daophot_trial{trial}.csv"
        if truth.empty or not moffat_tsv.exists() or not dao_csv.exists():
            print(f"[건너뜀] trial {trial}: 입력 없음")
            continue

        engines = {
            "APEX ePSF": apex_epsf_scored(truth, recovery, gates),
            "APEX Moffat": apex_moffat_scored(truth, moffat_tsv, gates,
                                              args.match_radius_px),
            "IRAF ALLSTAR": allstar_scored(truth, dao_csv, args.zmag,
                                           args.itime, args.match_radius_px),
        }
        scored = {}
        for name, frame in engines.items():
            frame, zp = remove_zeropoint(frame, min_snr=args.offset_min_snr,
                                         isolated_bins=isolated)
            scored[name] = frame
            print(f"  trial {trial} {name:>13}: 영점 {zp['zeropoint_offset_mag']:+.3f} mag "
                  f"({zp['n_offset_stars']}개)")

        # Same-star comparison needs every engine to have kept the star.
        both = np.ones(len(truth), dtype=bool)
        for frame in scored.values():
            both &= np.isfinite(frame["delta_mag"]).to_numpy()

        rows = []
        for name, frame in scored.items():
            rows += summarise(frame, name)
            for row in summarise(frame[both], name):
                if row["scope"] == "overall":
                    row["scope"], row["label"] = "matched", "all three kept"
                    rows.append(row)
                elif row["scope"] == "crowding":
                    row["scope"] = "matched-crowding"
                    rows.append(row)
        table = pd.DataFrame(rows)
        table["trial"] = trial
        parts.append(table)

    if not parts:
        raise SystemExit("채점된 회차가 없다")

    combined = pd.concat(parts, ignore_index=True)
    pooled = (combined.groupby(["engine", "scope", "label"], as_index=False)
              .agg(n_truth=("n_truth", "sum"), n_recovered=("n_recovered", "sum"),
                   bias_mag=("bias_mag", "median"),
                   scatter_mag=("scatter_mag", "median")))
    pooled["completeness"] = pooled["n_recovered"] / pooled["n_truth"]
    pooled.to_csv(args.output, index=False)

    order = ["APEX ePSF", "APEX Moffat", "IRAF ALLSTAR"]
    head = (f"\n{'engine':>14}{'scope':>18}{'label':>16}{'N':>6}{'rec':>6}"
            f"{'compl':>8}{'bias':>9}{'scatter':>9}")
    print(head)
    print("-" * len(head))
    for scope in ("overall", "matched", "matched-crowding"):
        part = pooled[pooled["scope"] == scope]
        for label in sorted(part["label"].unique()):
            for name in order:
                r = part[(part["label"] == label) & (part["engine"] == name)]
                if r.empty:
                    continue
                r = r.iloc[0]
                print(f"{name:>14}{scope:>18}{label:>16}{r['n_truth']:>6}"
                      f"{r['n_recovered']:>6}{r['completeness']:>8.2f}"
                      f"{r['bias_mag']:>9.3f}{r['scatter_mag']:>9.3f}")
        print()
    print(f"saved -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
