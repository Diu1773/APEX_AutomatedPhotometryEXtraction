"""Run PSF photometry across the Phase 3 targets and rebuild their CMDs.

The Phase 3 reprocessing ran steps 1-7 and then step 10, skipping step 8,
because step 10 reads star IDs straight from the forced-aperture tables and
does not need PSF output. On M13 that shortcut cost 12-16 % in zero-point
scatter, 18-22 % of the calibrators, and a 76 -> 58 mmag widening of the giant
branch, so the remaining targets are measured the same way.

Two things are protected. The aperture zero-point directory is copied to
`cmd_zeropoint_APERTURE` before step 10 re-runs, because step 10 switches to
PSF silently once a valid signature exists and the aperture products are the
source of the published numbers. And each target is skipped when its PSF
signature is already present, so an interrupted batch resumes instead of
recomputing hours of fitting.

Nothing here touches raw frames, calibration, detection, WCS, or the master
catalogue. Step 8 writes a new directory and step 10 reads it.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
VENV_PY = REPO / ".venv-deploy" / "Scripts" / "python.exe"
PHASE3 = Path(r"E:\APEX_validation\phase3")

# Ordered cheapest first so a partial night still finishes whole targets.
TARGETS = ("NGC6811", "M67", "M3", "M5")


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def run(script: str, config: Path) -> int:
    return subprocess.run(
        [str(VENV_PY), "-X", "utf8", str(REPO / "scripts" / script),
         "--params", str(config)],
        cwd=str(REPO),
    ).returncode


def run_target(target: str, *, force: bool) -> bool:
    root = PHASE3 / target
    config = root / "apex_config.json"
    result = root / "result"
    if not config.exists():
        log(f"[SKIP] {target}: apex_config.json 없음")
        return False

    frames = sorted((result / "step7_forced_phot").glob("photometry_*.tsv"))
    if not frames:
        log(f"[SKIP] {target}: step7 산출물 없음")
        return False

    signature = result / "cmd_psf" / "psf_output_signature.json"
    if signature.exists() and not force:
        log(f"[SKIP] {target}: PSF 서명이 이미 있다 ({signature.name})")
    else:
        log(f"[STEP8] {target}: {len(frames)} 프레임 PSF 측광 시작")
        started = time.perf_counter()
        if run("run_step8_headless.py", config) != 0:
            log(f"[FAIL] {target}: step8 실패")
            return False
        log(f"[STEP8] {target}: 완료 {time.perf_counter() - started:.0f}s")

    # The aperture products are the published Phase 3 numbers. Preserve them
    # before step 10 overwrites cmd_zeropoint with the PSF solution.
    aperture = result / "cmd_zeropoint"
    preserved = result / "cmd_zeropoint_APERTURE"
    if aperture.exists() and not preserved.exists():
        shutil.copytree(aperture, preserved)
        log(f"[KEEP] {target}: 조리개 영점 보존 -> {preserved.name}")

    log(f"[STEP10] {target}: PSF 로 영점·CMD 재구성")
    if run("run_step10_headless.py", config) != 0:
        log(f"[FAIL] {target}: step10 실패")
        return False

    coefficients = result / "cmd_zeropoint" / "zp_fit_coefficients.csv"
    source = "unknown"
    if coefficients.exists():
        text = coefficients.read_text(encoding="utf-8", errors="replace")
        source = "psf" if ",psf," in text else "aperture"
    log(f"[DONE] {target}: step10 photometry_source={source}")
    if source != "psf":
        log(f"[WARN] {target}: step10 이 PSF 를 채택하지 않았다 — 서명을 확인할 것")
    return source == "psf"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", default=",".join(TARGETS))
    parser.add_argument("--force", action="store_true",
                        help="PSF 서명이 있어도 step8 을 다시 돌린다")
    args = parser.parse_args()

    targets = [t.strip() for t in args.targets.split(",") if t.strip()]
    log(f"대상 {len(targets)}개: {', '.join(targets)}")

    results: dict[str, bool] = {}
    started = time.perf_counter()
    for target in targets:
        results[target] = run_target(target, force=args.force)

    log(f"전체 {time.perf_counter() - started:.0f}s")
    for target, ok in results.items():
        log(f"  {target}: {'PSF 채택' if ok else '실패/미채택'}")
    return 0 if all(results.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
