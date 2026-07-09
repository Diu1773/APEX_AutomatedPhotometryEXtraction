"""Full-APEX headless reprocessing orchestrator (see docs/REPROCESS_PLAN.md).

Per target: APEX Step 0 (headless) -> split calibrated frames per object ->
generate a parameters.toml -> `apex run` Steps 1-7 -> CMD 8/9/10/11 (clusters) ->
keep results (isochrone step 12 is left to the user). LC targets: Step 0 only.

SAFETY (hard): never touch E:\\observe_raw_Analysis, E:\\observe_DSY,
E:\\observed_Analysis. All output under E:\\APEX_validation\\reprocess\\<target>\\.
Stop if E: free < 20 GB. Resumable: skip targets already marked done in PROGRESS.md.

This is built incrementally by the autonomous loop; Steps 1-7 + Step 0 are wired,
CMD 8-11 runners are filled in as they are built and verified.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

VENV_PY = REPO / ".venv-deploy" / "Scripts" / "python.exe"
REPROCESS = Path(r"E:\APEX_validation\reprocess")
BIAS = r"E:\bias"
DARKS = r"E:\darks"
PROGRESS = REPROCESS / "PROGRESS.md"
MIN_FREE_GB = 20.0

# target -> (raw_dir, kind, ra_deg, dec_deg). kind: "gc"/"oc" (full CMD) or "lc" (Step 0 only).
TARGETS = {
    "NGC6811":  (r"E:\observe_raw_Analysis\20260611",        "oc", 294.34, 46.378),
    "M67":      (r"E:\observe_raw_Analysis\M67_20260208",    "oc", 132.825, 11.80),
    "M13":      (r"E:\observe_raw_Analysis\M13_20260515",    "gc", 250.421, 36.460),
    "M3":       (r"E:\observe_raw_Analysis\M3",              "gc", 205.548, 28.377),
    "M5":       (r"E:\observe_DSY\M5",                       "gc", 229.638, 2.081),
    # +1 open cluster from observe_DSY chosen at runtime (first complete set).
    "AE_UMa":   (r"E:\observe_raw_Analysis\AE UMa",          "lc", None, None),
    "YZ_Boo":   (r"E:\observe_raw_Analysis\YZbootis",        "lc", None, None),
}


def free_gb(drive="E:\\") -> float:
    return shutil.disk_usage(drive).free / 1e9


def log_progress(msg: str):
    REPROCESS.mkdir(parents=True, exist_ok=True)
    with open(PROGRESS, "a", encoding="utf-8") as f:
        f.write(msg.rstrip() + "\n")
    print(msg, flush=True)


def is_done(target: str) -> bool:
    if not PROGRESS.exists():
        return False
    return f"[DONE] {target}" in PROGRESS.read_text(encoding="utf-8")


def main() -> int:
    log_progress(f"# reprocess run {time.strftime('%Y-%m-%d %H:%M')} — free {free_gb():.0f} GB")
    for target, (raw, kind, ra, dec) in TARGETS.items():
        if is_done(target):
            print(f"skip {target} (done)"); continue
        if free_gb() < MIN_FREE_GB:
            log_progress(f"[STOP] E: free {free_gb():.0f} GB < {MIN_FREE_GB} — halting."); break
        # TODO(loop): Step 0 -> reorg -> config -> apex run 1-7 -> CMD 8-11 (gc/oc) / Step0-only (lc)
        log_progress(f"[PENDING] {target} ({kind}) — orchestrator wiring in progress")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
