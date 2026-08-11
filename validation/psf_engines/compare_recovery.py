"""Score APEX step 8 and IRAF ALLSTAR against the same implanted stars.

Neither engine is truth, so comparing them to each other can only say they
disagree, never which one is right. Artificial stars fix that: a star of known
flux is added to the real frame, both engines measure the same pixels, and the
error each one makes is a fact rather than a difference.

Everything the two engines are scored on comes from the same injection table —
the same positions, the same fluxes, the same crowding bins. APEX's numbers are
read from the recovery table its own benchmark wrote; DAOPHOT's are matched to
the truth positions here, with the same match radius.

Two things this has to get right, or it produces a false winner
--------------------------------------------------------------
**The magnitude scales are not the same.** ALLSTAR's magnitudes inherit the
zero point of the aperture `phot` used, and no aperture correction is applied
by the chain, so every ALLSTAR magnitude is faint by the light that aperture
missed — about 0.3 mag for a 1·FWHM aperture. APEX's ePSF flux is a total
flux. Subtracting one from the other measures the aperture correction, not the
engine. A constant offset is therefore removed from each engine before any bias
is quoted, measured on the bright isolated stars where neither engine is under
stress, which is the same convention `validation/psf_archive/` already uses.
The raw offset is reported too, so it can be checked.

**"Recovered" does not mean the same thing for the two engines.** APEX's
benchmark adds the implanted positions to the step-7 forced catalogue, so step 8
fits at every one of them and returns a number for all of them by construction.
ALLSTAR rejects stars it cannot fit. Reading APEX's 100 % against ALLSTAR's
fraction as a completeness win would be comparing "always answers" with
"answers when confident". `--apex-mode` selects which question is being asked:
`forced` keeps every fitted position, `gated` keeps only those that pass APEX's
own post-fit quality policy, which is the setting comparable to a rejecting
engine.

Reading the output
------------------
`completeness` is the fraction of implanted stars the engine returned a usable
magnitude for, under the mode above. `bias` is the median of measured minus
true magnitude after the constant offset is removed, so a positive bias means
the star was reported fainter than it was — the signature of lost flux.
`scatter` is a robust (MAD-based) spread, because a handful of catastrophic
blends would otherwise set the number for every bin they land in.

The bins that matter are the crowding bins. A PSF engine's whole claim is that
it separates blends, so its advantage should appear where neighbours are close
and vanish where they are not. An engine that wins everywhere equally is
probably being credited for a zero-point difference that survived the offset
removal.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

REPO = Path(__file__).absolute().parents[2]


def daophot_flux_adu(mag: np.ndarray, zmag: float, itime: float) -> np.ndarray:
    """Invert IRAF's magnitude definition, mag = zmag - 2.5*log10(counts/itime)."""
    return itime * np.power(10.0, (zmag - mag) / 2.5)


def match_to_truth(truth: pd.DataFrame, table: pd.DataFrame,
                   radius_px: float) -> np.ndarray:
    """Row index in `table` nearest each implanted star, or -1 beyond `radius_px`."""
    if table.empty:
        return np.full(len(truth), -1)
    tree = cKDTree(table[["x", "y"]].to_numpy(float))
    distance, index = tree.query(truth[["x_true", "y_true"]].to_numpy(float), k=1)
    return np.where(distance <= radius_px, index, -1)


def robust_scatter(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if values.size < 3:
        return float("nan")
    return float(1.4826 * np.median(np.abs(values - np.median(values))))


def remove_zeropoint(frame: pd.DataFrame, *, min_snr: float,
                     isolated_bins: tuple[str, ...]) -> tuple[pd.DataFrame, dict]:
    """Subtract the engine's constant offset, measured where it is not stressed.

    Bright, isolated implanted stars are where both engines should be at their
    best, so whatever offset remains there is a magnitude-scale difference
    (ALLSTAR's missing aperture correction) rather than a failure to measure.
    Falling back to the whole sample when that subset is empty keeps a small
    run from silently reporting an uncorrected scale.
    """
    delta = frame["delta_mag"]
    reference = frame[
        np.isfinite(delta)
        & (pd.to_numeric(frame["target_snr"], errors="coerce") >= min_snr)
        & frame["crowding_bin"].astype(str).isin(isolated_bins)
    ]
    basis = "bright isolated"
    if len(reference) < 3:
        reference = frame[np.isfinite(delta)]
        basis = "whole sample (too few bright isolated stars)"
    offset = float(np.median(reference["delta_mag"])) if len(reference) else 0.0

    out = frame.copy()
    out["delta_mag_raw"] = delta
    out["delta_mag"] = delta - offset
    return out, {"zeropoint_offset_mag": offset, "n_offset_stars": int(len(reference)),
                 "offset_basis": basis}


def summarise(frame: pd.DataFrame, engine: str) -> list[dict]:
    """Completeness, bias and scatter overall and per bin."""
    rows: list[dict] = []
    groups: list[tuple[str, str, pd.DataFrame]] = [("overall", "all", frame)]
    for column, scope in (("crowding_bin", "crowding"), ("target_snr", "snr")):
        if column in frame.columns:
            for label, group in frame.groupby(column, dropna=False):
                groups.append((scope, str(label), group))

    for scope, label, group in groups:
        recovered = group[np.isfinite(group["delta_mag"])]
        rows.append({
            "engine": engine, "scope": scope, "label": label,
            "n_truth": int(len(group)),
            "n_recovered": int(len(recovered)),
            "completeness": float(len(recovered) / len(group)) if len(group) else float("nan"),
            "bias_mag": float(np.median(recovered["delta_mag"])) if len(recovered) else float("nan"),
            "scatter_mag": robust_scatter(recovered["delta_mag"].to_numpy(float)),
        })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--apex-recovery", required=True)
    ap.add_argument("--daophot", required=True)
    ap.add_argument("--zmag", type=float, default=25.0)
    ap.add_argument("--itime", type=float, default=60.0)
    ap.add_argument("--match-radius-px", type=float, default=1.5)
    ap.add_argument("--trial", type=int, default=None,
                    help="여러 회차로 주입한 경우 이 회차만 채점한다. "
                         "회차마다 프레임이 다르므로 DAOPHoT 결과도 같은 회차의 것을 줄 것")
    ap.add_argument("--apex-mode", default="gated", choices=("forced", "gated"),
                    help="forced = 적합한 모든 위치, gated = APEX 자체 품질정책을 "
                         "통과한 것만. ALLSTAR 는 버리는 엔진이므로 gated 가 대응된다")
    # Defaults are APEX's own postfit_* settings.
    ap.add_argument("--postfit-snr-min", type=float, default=3.0)
    ap.add_argument("--postfit-qfit-max", type=float, default=3.0)
    ap.add_argument("--postfit-chi2-max", type=float, default=25.0)
    ap.add_argument("--offset-min-snr", type=float, default=50.0)
    ap.add_argument("--isolated-bins", default="3-6 FWHM|6-inf FWHM")
    ap.add_argument("--output", default=None)
    args = ap.parse_args()

    truth = pd.read_csv(args.truth)
    apex = pd.read_csv(args.apex_recovery)
    dao = pd.read_csv(args.daophot)

    # Each trial implants a different set of stars into its own copy of the
    # frame, so a DAOPHOT run belongs to exactly one trial. Matching across
    # trials would pair a measurement with a star that was never in that image.
    if args.trial is not None:
        for name, frame in (("truth", truth), ("apex", apex)):
            if "trial" not in frame.columns:
                raise SystemExit(f"--trial 을 줬는데 {name} 표에 trial 열이 없다")
        truth = truth[truth["trial"] == args.trial].reset_index(drop=True)
        apex = apex[apex["trial"] == args.trial].reset_index(drop=True)
        if truth.empty:
            raise SystemExit(f"trial {args.trial} 의 주입별이 없다")
        print(f"trial {args.trial}: 주입 {len(truth)}개\n")

    # APEX's benchmark already scored itself against this truth table; reusing
    # its delta_mag keeps the comparison on APEX's own definition rather than a
    # second one invented here.
    key = ["x_true", "y_true"]
    # `psf_recovered` only records that step 8 returned a fit — it is True for
    # every forced position, so it cannot stand in for a rejection decision.
    # APEX's actual post-fit policy is these four thresholds, the same ones
    # `load_lightcurve_frame_photometry` applies before a magnitude is used.
    quality = ["flags_psf", "snr_psf", "qfit", "reduced_chi2"]
    carry = key + ["delta_mag"] + [c for c in quality if c in apex.columns]
    apex_scored = truth[key + ["target_snr", "crowding_bin"]].merge(
        apex[carry].drop_duplicates(subset=key), on=key, how="left")

    if args.apex_mode == "gated":
        missing = [c for c in quality if c not in apex_scored.columns]
        if missing:
            raise SystemExit(f"--apex-mode gated 인데 회수표에 없는 열: {missing}")
        kept = (
            (pd.to_numeric(apex_scored["flags_psf"], errors="coerce") == 0)
            & (pd.to_numeric(apex_scored["snr_psf"], errors="coerce") >= args.postfit_snr_min)
            & (pd.to_numeric(apex_scored["qfit"], errors="coerce") <= args.postfit_qfit_max)
            & (pd.to_numeric(apex_scored["reduced_chi2"], errors="coerce") <= args.postfit_chi2_max)
        )
        # A star APEX fitted but its own policy would discard is not recovered
        # under this question; blank it so it counts as a miss exactly the way
        # an ALLSTAR rejection does.
        apex_scored.loc[~kept, "delta_mag"] = np.nan

    dao = dao.copy()
    dao["mag"] = pd.to_numeric(dao["mag"], errors="coerce")
    dao = dao[np.isfinite(dao["mag"])]
    index = match_to_truth(truth, dao, args.match_radius_px)
    measured = np.where(index >= 0,
                        daophot_flux_adu(dao["mag"].to_numpy(float)[index],
                                         args.zmag, args.itime),
                        np.nan)
    true_flux = truth["flux_realized_adu"].to_numpy(float)
    dao_scored = truth[["x_true", "y_true", "target_snr", "crowding_bin"]].copy()
    # Same sign convention as APEX's benchmark: measured minus true, so a
    # positive value means the star was reported too faint.
    dao_scored["delta_mag"] = -2.5 * np.log10(measured / true_flux)

    isolated = tuple(args.isolated_bins.split("|"))
    apex_scored, apex_zp = remove_zeropoint(
        apex_scored, min_snr=args.offset_min_snr, isolated_bins=isolated)
    dao_scored, dao_zp = remove_zeropoint(
        dao_scored, min_snr=args.offset_min_snr, isolated_bins=isolated)

    # A rejecting engine measures its scatter on whatever it chose to keep. If
    # ALLSTAR discards the blends it cannot fit and APEX keeps them, ALLSTAR's
    # spread is computed on an easier sample and the comparison rewards it for
    # being selective rather than for being accurate. Restricting both engines
    # to the stars both kept removes that advantage; the difference in what
    # each one keeps is already reported as completeness.
    both = (np.isfinite(apex_scored["delta_mag"])
            & np.isfinite(dao_scored["delta_mag"])).to_numpy()
    apex_matched = apex_scored[both].copy()
    dao_matched = dao_scored[both].copy()

    rows = (summarise(apex_scored, f"APEX step8 ({args.apex_mode})")
            + summarise(dao_scored, "IRAF ALLSTAR"))
    # Broken out by crowding as well, because "equivalent overall" can hide the
    # regime that matters: a PSF engine earns its place in the blended bin, and
    # an engine that is equal on average while losing there is not equal.
    for frame, engine in ((apex_matched, f"APEX step8 ({args.apex_mode})"),
                          (dao_matched, "IRAF ALLSTAR")):
        for row in summarise(frame, engine):
            if row["scope"] == "overall":
                row["scope"], row["label"] = "matched", "both kept"
                rows.append(row)
            elif row["scope"] == "crowding":
                row["scope"] = "matched-crowding"
                rows.append(row)
    table = pd.DataFrame(rows)

    print(f"등급 영점 제거 — 기준: SNR >= {args.offset_min_snr:.0f} 이고 "
          f"{' 또는 '.join(isolated)}")
    for name, zp in (("APEX step8", apex_zp), ("IRAF ALLSTAR", dao_zp)):
        print(f"  {name:>13}: {zp['zeropoint_offset_mag']:+.3f} mag "
              f"(별 {zp['n_offset_stars']}개, {zp['offset_basis']})")
    print()

    header = (f"{'engine':>13}{'scope':>10}{'label':>16}{'N':>5}{'rec':>5}"
              f"{'compl':>8}{'bias':>9}{'scatter':>9}")
    print(header)
    print("-" * len(header))
    for scope in ("overall", "matched", "matched-crowding", "crowding", "snr"):
        part = table[table["scope"] == scope].sort_values(["label", "engine"])
        for _, r in part.iterrows():
            print(f"{r['engine']:>13}{r['scope']:>10}{r['label']:>16}"
                  f"{r['n_truth']:>5}{r['n_recovered']:>5}"
                  f"{r['completeness']:>8.2f}{r['bias_mag']:>9.3f}"
                  f"{r['scatter_mag']:>9.3f}")
        print()

    output = Path(args.output) if args.output else (
        REPO / "validation" / "psf_engines" / "recovery_comparison.csv")
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, index=False)
    (output.parent / f"{output.stem}_inputs.json").write_text(json.dumps({
        "truth": str(args.truth), "apex_recovery": str(args.apex_recovery),
        "daophot": str(args.daophot), "zmag": args.zmag, "itime": args.itime,
        "match_radius_px": args.match_radius_px, "apex_mode": args.apex_mode,
        "zeropoint": {"APEX step8": apex_zp, "IRAF ALLSTAR": dao_zp,
                      "min_snr": args.offset_min_snr,
                      "isolated_bins": list(isolated)},
    }, indent=1), encoding="utf-8")
    print(f"saved -> {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
