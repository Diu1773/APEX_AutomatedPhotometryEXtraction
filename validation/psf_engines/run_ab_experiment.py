"""Run the whole APEX-vs-ALLSTAR artificial-star experiment unattended.

The preliminary comparison rested on 25 implanted stars, which left most
magnitude × crowding cells with two or three measurements. This raises the count
to a level where the cells mean something, and does the bookkeeping that the
larger sample forces.

Why several trials rather than one big injection
------------------------------------------------
Implanted stars are themselves neighbours. Adding 400 stars to a frame that
holds 1,600 raises the source density by a quarter, which changes the very
crowding the experiment is trying to measure. Splitting the same total across
trials keeps each frame close to the original — 200 into 1,600 is a 12 %
perturbation — at the cost of one extra Step 4 + Step 8 pass per trial.

Trials are not interchangeable afterwards. Each one implants a different set of
stars into its own copy of the frame, so a DAOPHOT run belongs to exactly one
trial and scoring has to be done per trial before the tables are stacked.
Matching across trials would pair a measurement with a star that was never in
that image.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

REPO = Path(__file__).absolute().parents[2]
VENV = REPO / ".venv-deploy" / "Scripts" / "python.exe"
HERE = Path(__file__).absolute().parent


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def run(argv: list[str], what: str) -> None:
    log(f"시작 — {what}")
    started = time.perf_counter()
    completed = subprocess.run(argv, cwd=str(REPO), text=True)
    if completed.returncode != 0:
        raise SystemExit(f"실패 ({completed.returncode}): {what}")
    log(f"완료 — {what} · {time.perf_counter() - started:.0f}s")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", default="M13")
    ap.add_argument("--frame", default="pp_messier13-0005-B.fit")
    ap.add_argument("--fwhm", type=float, required=True)
    ap.add_argument("--trials", type=int, default=2)
    ap.add_argument("--injections", type=int, default=200)
    ap.add_argument("--fitrad-fwhm", type=float, default=1.7,
                    help="APEX 의 포위에너지 90 %% 창에 맞춘 값. 이 교란요인은 "
                         "2026-08-12 에 결론을 안 바꾼다고 확인됐다")
    ap.add_argument("--phase3", default=r"E:\APEX_validation\phase3")
    ap.add_argument("--work", default=r"E:\APEX_validation\psf_engines")
    args = ap.parse_args()

    target = Path(args.phase3) / args.target
    work = Path(args.work) / f"{args.target}_ast{args.trials}x{args.injections}"
    log(f"{args.target} · {args.frame} · {args.trials}회 x {args.injections}개 "
        f"· 산출 {work}")

    run([str(VENV), "-X", "utf8", "validation/run_psf_artificial_stars.py",
         "--source-fits", str(target / "sci" / args.frame),
         "--baseline-result-dir", str(target / "result"),
         "--parameter-file", str(target / "apex_config.json"),
         "--output-dir", str(work),
         "--trials", str(args.trials), "--injections", str(args.injections)],
        "APEX 주입 + step4 + step8")

    scored: list[pd.DataFrame] = []
    for trial in range(1, args.trials + 1):
        tdir = work / f"trial_{trial:04d}"
        injected = tdir / "data" / args.frame
        positions = (tdir / "result" / "step7_forced_phot"
                     / f"photometry_{args.frame}.tsv")
        if not injected.exists() or not positions.exists():
            log(f"[건너뜀] trial {trial}: 주입 프레임이나 step7 표가 없다")
            continue

        dao_csv = work / f"daophot_trial{trial}.csv"
        run([str(VENV), "-X", "utf8", str(HERE / "daophot_allstar.py"),
             "--frame", str(injected), "--positions", str(positions),
             "--output", str(dao_csv),
             "--workdir", str(work / f"daophot_work_trial{trial}"),
             "--fwhm", str(args.fwhm),
             "--fitrad-fwhm", str(args.fitrad_fwhm)],
            f"IRAF ALLSTAR trial {trial}")

        out = work / f"comparison_trial{trial}.csv"
        run([str(VENV), "-X", "utf8", str(HERE / "compare_recovery.py"),
             "--truth", str(work / "truth.csv"),
             "--apex-recovery", str(work / "recovery.csv"),
             "--daophot", str(dao_csv), "--trial", str(trial),
             "--output", str(out)],
            f"채점 trial {trial}")
        table = pd.read_csv(out)
        table["trial"] = trial
        scored.append(table)

    if not scored:
        raise SystemExit("채점된 회차가 없다")

    combined = pd.concat(scored, ignore_index=True)
    # Completeness and the robust spread are both ratios/medians, so the trials
    # are pooled by re-deriving them from the counts rather than by averaging
    # per-trial values, which would weight a thin cell like a full one.
    pooled = (combined.groupby(["engine", "scope", "label"], as_index=False)
              .agg(n_truth=("n_truth", "sum"),
                   n_recovered=("n_recovered", "sum"),
                   bias_mag=("bias_mag", "median"),
                   scatter_mag=("scatter_mag", "median"),
                   trials=("trial", "nunique")))
    pooled["completeness"] = pooled["n_recovered"] / pooled["n_truth"]

    combined.to_csv(work / "comparison_all_trials.csv", index=False)
    pooled.to_csv(HERE / "recovery_comparison_pooled.csv", index=False)

    header = (f"{'engine':>20}{'scope':>10}{'label':>16}{'N':>6}{'rec':>6}"
              f"{'compl':>8}{'bias':>9}{'scatter':>9}")
    print("\n=== 회차 통합 ===")
    print(header)
    print("-" * len(header))
    for scope in ("overall", "matched", "crowding", "snr"):
        part = pooled[pooled["scope"] == scope].sort_values(["label", "engine"])
        for _, r in part.iterrows():
            print(f"{r['engine']:>20}{r['scope']:>10}{r['label']:>16}"
                  f"{r['n_truth']:>6}{r['n_recovered']:>6}"
                  f"{r['completeness']:>8.2f}{r['bias_mag']:>9.3f}"
                  f"{r['scatter_mag']:>9.3f}")
        print()
    log(f"saved -> {HERE / 'recovery_comparison_pooled.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
