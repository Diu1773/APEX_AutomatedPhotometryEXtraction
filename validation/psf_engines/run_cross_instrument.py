"""One unattended pass of the three-engine benchmark on a new instrument.

The chain is injection -> APEX step 8 -> DAOPHOT/ALLSTAR -> photutils ->
scoring, and it has been assembled by hand twice now (M13, NGC 5985). Doing it
a third time from the shell is how a step gets forgotten at 2 a.m., so it lives
here with the per-target numbers as arguments and every stage logged.

Each stage is skipped when its product already exists, so a re-run after a
failure resumes rather than repeating the expensive parts.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).absolute().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))

PYTHON = str(REPO / ".venv-deploy" / "Scripts" / "python.exe")


def run(label: str, command: list[str], log: Path) -> bool:
    print(f"[{time.strftime('%H:%M:%S')}] {label} …", flush=True)
    started = time.perf_counter()
    with log.open("w", encoding="utf-8", errors="replace") as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT,
                                text=True)
    elapsed = time.perf_counter() - started
    ok = result.returncode == 0
    print(f"[{time.strftime('%H:%M:%S')}] {label} "
          f"{'완료' if ok else '실패'} ({elapsed / 60:.1f}분) -> {log.name}", flush=True)
    if not ok:
        print(log.read_text(encoding="utf-8", errors="replace")[-1500:], flush=True)
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--name", required=True, help="label used for the work dir")
    ap.add_argument("--source-fits", required=True)
    ap.add_argument("--baseline-result-dir", required=True)
    ap.add_argument("--parameter-file", required=True)
    ap.add_argument("--frame", required=True, help="basename of the frame")
    ap.add_argument("--fwhm-px", type=float, required=True)
    ap.add_argument("--gain", type=float, required=True)
    ap.add_argument("--background-rms", type=float, required=True)
    ap.add_argument("--pixel-scale-arcsec", type=float, required=True)
    ap.add_argument("--exptime", type=float, required=True)
    ap.add_argument("--trials", type=int, default=2)
    ap.add_argument("--injections", type=int, default=200)
    ap.add_argument("--fit-shape", type=int, default=0,
                    help="photutils fit window; 0 = 2.4x FWHM rounded odd")
    ap.add_argument("--fit-shape-sweep", type=int, nargs="*", default=[],
                    help="extra photutils windows, to show the choice is not "
                         "what decided the result")
    ap.add_argument("--work-root", default=r"E:\APEX_validation\psf_engines")
    args = ap.parse_args()

    work = Path(args.work_root) / args.name
    work.mkdir(parents=True, exist_ok=True)
    logs = work / "logs"
    logs.mkdir(exist_ok=True)
    fit_shape = args.fit_shape or (int(round(2.4 * args.fwhm_px)) // 2 * 2 + 1)
    print(f"작업 {work}\n  FWHM {args.fwhm_px:.2f}px · gain {args.gain} · "
          f"배경RMS {args.background_rms:.2f} · 적합창 {fit_shape}", flush=True)

    if not (work / "truth.csv").exists():
        ok = run("인공별 주입 + APEX step8", [
            PYTHON, "-X", "utf8", str(REPO / "validation" / "run_psf_artificial_stars.py"),
            "--source-fits", args.source_fits,
            "--baseline-result-dir", args.baseline_result_dir,
            "--output-dir", str(work),
            "--parameter-file", args.parameter_file,
            "--trials", str(args.trials), "--injections", str(args.injections),
            "--fwhm-px", str(args.fwhm_px), "--gain-e-per-adu", str(args.gain),
            "--background-rms-adu", str(args.background_rms),
            "--pixel-scale-arcsec", str(args.pixel_scale_arcsec),
            "--inject-kernel", "moffat", "--flux-scale", "off",
            "--use-grouper", "off", "--fit-window-mode", "manual",
            "--fit-shape-fwhm-mult", "2.4",
            "--forced-match-radius-fwhm", "1.25",
            "--postfit-qfit-noise-max", "3.0",
        ], logs / "01_inject.log")
        if not ok:
            return 1
    else:
        print("인공별 주입: 이미 있음, 건너뜀", flush=True)

    for trial in range(1, args.trials + 1):
        out = work / f"daophot_trial{trial}.csv"
        if out.exists():
            print(f"DAOPHOT 회차 {trial}: 이미 있음, 건너뜀", flush=True)
            continue
        run(f"DAOPHOT ALLSTAR 회차 {trial}", [
            PYTHON, "-X", "utf8", str(HERE / "daophot_allstar.py"),
            "--frame", str(work / f"trial_{trial:04d}" / "data" / args.frame),
            "--positions", str(work / f"trial_{trial:04d}" / "result"
                               / "step7_forced_phot" / f"photometry_{args.frame}.tsv"),
            "--output", str(out),
            "--workdir", str(work / f"daophot_work_trial{trial}"),
            "--fwhm", str(args.fwhm_px), "--fitrad-fwhm", "1.7",
            "--calgorithm", "none", "--varorder", "0", "--nclean", "0",
        ], logs / f"02_daophot_t{trial}.log")

    for window in [fit_shape, *args.fit_shape_sweep]:
        tag = "" if window == fit_shape else f"_fs{window}"
        run(f"photutils (적합창 {window})", [
            PYTHON, "-X", "utf8", str(HERE / "photutils_engine.py"),
            "--work", str(work), "--frame", args.frame, "--apex-run", "trial",
            "--trials", str(args.trials), "--fit-shape", str(window),
            "--gain", str(args.gain), "--background-rms", str(args.background_rms),
            "--grouper-fwhm", "1.5", "--fwhm-px", str(args.fwhm_px),
        ], logs / f"03_photutils{tag or '_default'}.log")
        if tag:
            for trial in range(1, args.trials + 1):
                src = work / f"photutils_trial{trial}.csv"
                if src.exists():
                    src.replace(work / f"photutils{tag}_trial{trial}.csv")

    # Re-run the default last so the scored files are the chosen window.
    if args.fit_shape_sweep:
        run(f"photutils 재실행 (선택 창 {fit_shape})", [
            PYTHON, "-X", "utf8", str(HERE / "photutils_engine.py"),
            "--work", str(work), "--frame", args.frame, "--apex-run", "trial",
            "--trials", str(args.trials), "--fit-shape", str(fit_shape),
            "--gain", str(args.gain), "--background-rms", str(args.background_rms),
            "--grouper-fwhm", "1.5", "--fwhm-px", str(args.fwhm_px),
        ], logs / "03_photutils_final.log")

    manifest = {
        "name": args.name, "frame": args.frame,
        "fwhm_px": args.fwhm_px, "gain_e_per_adu": args.gain,
        "background_rms_adu": args.background_rms,
        "pixel_scale_arcsec": args.pixel_scale_arcsec,
        "exptime_s": args.exptime, "photutils_fit_shape": fit_shape,
        "photutils_fit_shape_sweep": args.fit_shape_sweep,
        "trials": args.trials, "injections": args.injections,
    }
    (work / "cross_instrument_inputs.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n입력 명세 -> {work / 'cross_instrument_inputs.json'}", flush=True)
    print("다음: score_cross_instrument.py 로 채점", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
