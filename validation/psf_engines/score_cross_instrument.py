"""Score APEX, ALLSTAR and photutils on one target, with the gates off.

Gates are a policy layer, not part of an engine: ALLSTAR has none, APEX has
four, and comparing APEX's filtered output against ALLSTAR's raw output
measures conservatism rather than photometry. Everything here is scored on
whatever each engine measured, which is the only apples-to-apples version.

Writes one row per engine to `cross_instrument_scores.csv` in the work
directory so a night's runs can be collected without re-reading logs.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

HERE = Path(__file__).absolute().parent
sys.path.insert(0, str(HERE))

from compare_recovery import (  # noqa: E402
    daophot_flux_adu, load_daophot, match_to_truth, remove_zeropoint,
    robust_scatter,
)

ISOLATED = ("3-6 FWHM", "6-inf FWHM")
TIGHT = "0.75-1.5 FWHM"


def _scored(truth: pd.DataFrame, x, y, mag) -> pd.DataFrame:
    idx = match_to_truth(truth, pd.DataFrame({"x": x, "y": y}), 1.5)
    out = truth.copy()
    out["delta_mag"] = (np.where(idx >= 0, np.asarray(mag)[idx], np.nan)
                        + 2.5 * np.log10(truth["flux_realized_adu"].to_numpy(float)))
    return out


def apex(work: Path, frame: str, trials: int) -> pd.DataFrame:
    parts = []
    for trial in range(1, trials + 1):
        truth = pd.read_csv(work / "truth.csv")
        truth = truth[truth["trial"] == trial].reset_index(drop=True)
        table = pd.read_csv(work / f"trial_{trial:04d}" / "result" / "cmd_psf"
                            / f"photometry_{frame}.tsv", sep="\t")
        for column in ("x_fit", "y_fit", "mag_psf"):
            table[column] = pd.to_numeric(table[column], errors="coerce")
        table = table[np.isfinite(table["mag_psf"])].reset_index(drop=True)
        parts.append(_scored(truth, table["x_fit"], table["y_fit"],
                             table["mag_psf"].to_numpy()))
    return remove_zeropoint(pd.concat(parts, ignore_index=True),
                            min_snr=50.0, isolated_bins=ISOLATED)[0]


def photutils(work: Path, trials: int, suffix: str = "") -> pd.DataFrame | None:
    parts = []
    for trial in range(1, trials + 1):
        path = work / f"photutils{suffix}_trial{trial}.csv"
        if not path.exists():
            return None
        truth = pd.read_csv(work / "truth.csv")
        truth = truth[truth["trial"] == trial].reset_index(drop=True)
        table = pd.read_csv(path)
        table = table[np.isfinite(table["mag_psf"])].reset_index(drop=True)
        parts.append(_scored(truth, table["x_fit"], table["y_fit"],
                             table["mag_psf"].to_numpy()))
    return remove_zeropoint(pd.concat(parts, ignore_index=True),
                            min_snr=50.0, isolated_bins=ISOLATED)[0]


def allstar(work: Path, trials: int, exptime: float) -> pd.DataFrame | None:
    parts = []
    for trial in range(1, trials + 1):
        path = work / f"daophot_trial{trial}.csv"
        if not path.exists():
            return None
        truth = pd.read_csv(work / "truth.csv")
        truth = truth[truth["trial"] == trial].reset_index(drop=True)
        table = load_daophot(path)
        table["mag"] = pd.to_numeric(table["mag"], errors="coerce")
        table = table[np.isfinite(table["mag"])]
        idx = match_to_truth(truth, table, 1.5)
        flux = np.where(idx >= 0,
                        daophot_flux_adu(table["mag"].to_numpy(float)[idx],
                                         25.0, exptime), np.nan)
        row = truth.copy()
        row["delta_mag"] = -2.5 * np.log10(
            flux / truth["flux_realized_adu"].to_numpy(float))
        parts.append(row)
    return remove_zeropoint(pd.concat(parts, ignore_index=True),
                            min_snr=50.0, isolated_bins=ISOLATED)[0]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", required=True)
    ap.add_argument("--frame", required=True)
    ap.add_argument("--exptime", type=float, required=True)
    ap.add_argument("--trials", type=int, default=2)
    ap.add_argument("--label", default="")
    args = ap.parse_args()

    work = Path(args.work)
    total = len(pd.read_csv(work / "truth.csv"))
    engines: dict[str, pd.DataFrame] = {"APEX": apex(work, args.frame, args.trials)}
    for name, frame in (("photutils", photutils(work, args.trials)),
                        ("ALLSTAR", allstar(work, args.trials, args.exptime))):
        if frame is not None:
            engines[name] = frame

    label = args.label or work.name
    print(f"\n=== {label} · 게이트 없음 · 인공별 {total} ===\n")
    print(f'{"엔진":>11}{"건진수":>8}{"완전도":>8}{"편향":>11}{"산포":>9}')
    print("-" * 49)
    rows = []
    for name, frame in engines.items():
        values = frame["delta_mag"].to_numpy(float)
        kept = int(np.isfinite(values).sum())
        rows.append({"target": label, "engine": name, "n_injected": total,
                     "n_recovered": kept, "completeness": kept / total,
                     "bias_mag": float(np.nanmedian(values)),
                     "scatter_mag": robust_scatter(values[np.isfinite(values)])})
        print(f'{name:>11}{kept:>8}{kept / total:>8.2f}'
              f'{np.nanmedian(values) * 1000:>+9.1f}mm'
              f'{robust_scatter(values[np.isfinite(values)]):>9.4f}')

    common = np.ones(total, dtype=bool)
    for frame in engines.values():
        common &= np.isfinite(frame["delta_mag"]).to_numpy()
    bins = engines["APEX"]["crowding_bin"].to_numpy()
    print(f"\n공통 별 {int(common.sum())}개")
    print(f'{"엔진":>11}{"산포":>9}{"최혼잡 산포":>13}')
    print("-" * 34)
    for name, frame in engines.items():
        values = frame["delta_mag"].to_numpy(float)
        tight = common & (bins == TIGHT)
        rows[[r["engine"] for r in rows].index(name)].update({
            "n_common": int(common.sum()),
            "scatter_common": robust_scatter(values[common]),
            "scatter_tight": robust_scatter(values[tight]),
        })
        print(f'{name:>11}{robust_scatter(values[common]):>9.4f}'
              f'{robust_scatter(values[tight]):>13.4f}')

    print("\n짝지은 |APEX|-|X| (양수 = 상대가 우세)")
    rng = np.random.default_rng(11)
    for name in ("ALLSTAR", "photutils"):
        if name not in engines:
            continue
        diff = (engines["APEX"]["delta_mag"].abs()
                - engines[name]["delta_mag"].abs()).to_numpy(float)[common]
        interval = np.percentile(
            np.median(rng.choice(diff, (20000, diff.size), replace=True), axis=1),
            [2.5, 97.5])
        p = float(wilcoxon(diff).pvalue)
        print(f"   vs {name:>10}: 중앙 {np.median(diff) * 1000:+6.1f} mmag "
              f"[{interval[0] * 1000:+6.1f}, {interval[1] * 1000:+6.1f}] p={p:.4f}")
        rows[0][f"paired_vs_{name}_mmag"] = float(np.median(diff)) * 1000
        rows[0][f"paired_vs_{name}_p"] = p

    # photutils fit-window sweep, if the runner produced one.
    sweep = sorted(work.glob("photutils_fs*_trial1.csv"))
    if sweep:
        print("\nphotutils 적합창 민감도 (공통 별 산포)")
        for path in sweep:
            suffix = path.name[len("photutils"):-len("_trial1.csv")]
            frame = photutils(work, args.trials, suffix)
            if frame is None:
                continue
            values = frame["delta_mag"].to_numpy(float)
            ok = np.isfinite(values) & common
            print(f"   창 {suffix[3:]:>3}: 산포 {robust_scatter(values[ok]):.4f} "
                  f"(N={int(ok.sum())})")

    out = work / "cross_instrument_scores.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nsaved -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
