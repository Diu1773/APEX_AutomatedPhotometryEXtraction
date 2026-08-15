"""Does DAOPHOT's profile-error term earn its place as a default?

APEX weights each pixel by `background² + counts/gain` — photon and read noise
and nothing else. DAOPHOT adds a term proportional to the counts themselves
(`proferr`, 5 % by default), on the argument that the PSF model is never exactly
right and the error it makes is largest where the star is brightest. The effect
is to stop a bright core from dominating the fit, which is exactly the pixel a
blended neighbour has to compete with.

Whether that helps is an empirical question, and it is not obviously the same
answer on every telescope: the term's justification is *model* error, so an
instrument whose PSF the model fits well has less to gain. This pairs the same
implanted stars measured with the term off and on, per instrument.

Run `run_apex_variant.py --tag prof --set psf.profile_error_frac=0.05` first;
everything upstream of Step 8 is shared with the baseline by construction.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

HERE = Path(__file__).absolute().parent
sys.path.insert(0, str(HERE))

from compare_recovery import (  # noqa: E402
    match_to_truth, remove_zeropoint, robust_scatter,
)

ISOLATED = ("3-6 FWHM", "6-inf FWHM")
BIN_ORDER = ("0.75-1.5 FWHM", "1.5-3 FWHM", "3-6 FWHM", "6-inf FWHM")


def _score(work: Path, frame: str, trials: int, tag: str | None,
           zeropoint: str = "global") -> pd.DataFrame:
    """Δmag per implanted star, zeropoint removed the same way for both arms.

    ``per-trial`` exists for the QHY600 field, whose frames carry a ~90 mmag
    zeropoint swing between trials that is still unexplained. Left in, it is
    the largest term in |Δmag| and buries whatever the fit weighting did; it is
    removed identically from both arms, so it cannot manufacture a difference.
    """
    parts = []
    for trial in range(1, trials + 1):
        truth = pd.read_csv(work / "truth.csv")
        truth = truth[truth["trial"] == trial].reset_index(drop=True)
        folder = (work / f"trial_{trial:04d}" if tag is None
                  else work / f"apex{tag}_trial{trial}")
        table = pd.read_csv(folder / "result" / "cmd_psf"
                            / f"photometry_{frame}.tsv", sep="\t")
        for column in ("x_fit", "y_fit", "mag_psf"):
            table[column] = pd.to_numeric(table[column], errors="coerce")
        table = table[np.isfinite(table["mag_psf"])].reset_index(drop=True)
        index = match_to_truth(truth, pd.DataFrame({"x": table["x_fit"],
                                                    "y": table["y_fit"]}), 1.5)
        row = truth.copy()
        row["delta_mag"] = (
            np.where(index >= 0, table["mag_psf"].to_numpy()[index], np.nan)
            + 2.5 * np.log10(truth["flux_realized_adu"].to_numpy(float)))
        if zeropoint == "per-trial":
            row = remove_zeropoint(row, min_snr=50.0,
                                   isolated_bins=ISOLATED)[0]
        parts.append(row)
    joined = pd.concat(parts, ignore_index=True)
    if zeropoint == "per-trial":
        return joined
    return remove_zeropoint(joined, min_snr=50.0, isolated_bins=ISOLATED)[0]


def _stats(values: np.ndarray) -> tuple[int, float, float]:
    good = values[np.isfinite(values)]
    return (int(good.size),
            float(np.median(good)) if good.size else float("nan"),
            robust_scatter(good))


def _paired(base: np.ndarray, test: np.ndarray) -> tuple[int, float, float]:
    """Wilcoxon on |Δmag| where both arms measured the same star.

    Absolute error rather than signed: an arm that swaps a bright bias for an
    equal faint one has not improved, and the signed test would not notice.
    """
    both = np.isfinite(base) & np.isfinite(test)
    if both.sum() < 10:
        return int(both.sum()), float("nan"), float("nan")
    a, b = np.abs(base[both]), np.abs(test[both])
    shift = float(np.median(b - a)) * 1000.0
    if np.allclose(a, b):
        return int(both.sum()), shift, 1.0
    return int(both.sum()), shift, float(wilcoxon(a, b).pvalue)


def _report(label: str, base: pd.DataFrame, test: pd.DataFrame,
            total: int) -> dict:
    print(f"\n=== {label} ===")
    print(f'{"구간":>15}{"주입":>6}{"기본 n":>8}{"산포":>9}'
          f'{"prof n":>8}{"산포":>9}{"|오차| 변화":>12}{"p":>9}')
    print("-" * 76)
    rows = []
    scopes = [("전체", np.ones(len(base), bool))]
    scopes += [(name, (base["crowding_bin"] == name).to_numpy())
               for name in BIN_ORDER if (base["crowding_bin"] == name).any()]
    for name, mask in scopes:
        a = base.loc[mask, "delta_mag"].to_numpy(float)
        b = test.loc[mask, "delta_mag"].to_numpy(float)
        n_a, _, s_a = _stats(a)
        n_b, _, s_b = _stats(b)
        n_pair, shift, p = _paired(a, b)
        print(f"{name:>15}{int(mask.sum()):>6}{n_a:>8}{s_a:>9.4f}"
              f"{n_b:>8}{s_b:>9.4f}{shift:>+11.2f}m{p:>9.3g}")
        rows.append({"label": label, "scope": name, "injected": int(mask.sum()),
                     "n_base": n_a, "scatter_base": s_a,
                     "n_prof": n_b, "scatter_prof": s_b,
                     "n_paired": n_pair, "abs_error_shift_mmag": shift,
                     "wilcoxon_p": p})
    n_a = int(np.isfinite(base["delta_mag"]).sum())
    n_b = int(np.isfinite(test["delta_mag"]).sum())
    print(f"완전도 {n_a/total:.3f} -> {n_b/total:.3f}   "
          f"(주입 {total}, |오차| 변화가 음수면 prof 가 낫다)")
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True, action="append",
                    help="주입 실행 폴더. 기기마다 하나씩, 여러 번 쓸 수 있다")
    ap.add_argument("--frame", required=True, action="append")
    ap.add_argument("--label", default=[], action="append")
    ap.add_argument("--trials", type=int, default=2)
    ap.add_argument("--tag", default="prof")
    ap.add_argument("--base-tag", default=None,
                    help="비교 기준을 회차의 step8 이 아니라 다른 변형으로 둔다. "
                         "M13 은 prof 가 hybrid 위에 얹혀 있어 --base-tag hybrid")
    ap.add_argument("--zeropoint", choices=("global", "per-trial"),
                    default="global",
                    help="per-trial 은 회차마다 영점을 따로 뺀다. QHY600 처럼 "
                         "회차 간 영점이 흔들리는 데이터에서만 쓴다")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    if len(args.work) != len(args.frame):
        raise SystemExit("--work 와 --frame 의 개수가 같아야 한다")
    rows = []
    for position, (work_text, frame) in enumerate(zip(args.work, args.frame)):
        work = Path(work_text)
        label = (args.label[position] if position < len(args.label)
                 else work.name)
        total = int((pd.read_csv(work / "truth.csv")["trial"]
                     <= args.trials).sum())
        base = _score(work, frame, args.trials, args.base_tag,
                      args.zeropoint)
        test = _score(work, frame, args.trials, args.tag,
                      args.zeropoint)
        rows.extend(_report(label, base, test, total))

    if args.output:
        pd.DataFrame(rows).to_csv(args.output, index=False)
        meta = args.output.with_suffix(".inputs.json")
        meta.write_text(json.dumps(
            {"work": args.work, "frame": args.frame, "trials": args.trials,
             "tag": args.tag, "base_tag": args.base_tag,
             "zeropoint": args.zeropoint,
             "profile_error_frac_tested": 0.05,
             "match_radius_px": 1.5, "zeropoint": "min_snr=50, isolated bins"},
            indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"\n표 -> {args.output}\n명세 -> {meta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
