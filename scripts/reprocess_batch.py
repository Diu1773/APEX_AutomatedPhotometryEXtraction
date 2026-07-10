"""Full-APEX headless reprocessing orchestrator (see docs/REPROCESS_PLAN.md).

Per target: APEX Step 0 (headless) -> split calibrated frames per object ->
generate a parameters.toml -> `apex run` Steps 1-7 -> CMD 8/9/10/11 (clusters).
Isochrone (step 12) left to the user. LC targets: Step 0 only.

SAFETY (hard): never touch E:\\observe_raw_Analysis, E:\\observe_DSY,
E:\\observed_Analysis. All output under E:\\APEX_validation\\reprocess\\<target>\\.
Stop if E: free < 20 GB. Resumable: skip targets marked [DONE] in PROGRESS.md.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
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
TEMPLATE = REPO / "parameters.toml"
MIN_FREE_GB = 20.0

# target -> (raw_dir, kind, ra_deg, dec_deg). kind: gc/oc (full CMD) or lc (Step 0 only).
TARGETS = {
    "NGC6811": (r"E:\observe_raw_Analysis\20260611",     "oc", 294.34, 46.378),
    "M67":     (r"E:\observe_raw_Analysis\M67_20260208", "oc", 132.825, 11.80),
    "M13":     (r"E:\observe_raw_Analysis\M13_20260515", "gc", 250.421, 36.460),
    "M3":      (r"E:\observe_raw_Analysis\M3",           "gc", 205.548, 28.377),
    "M5":      (r"E:\observe_DSY\M5",                    "gc", 229.638, 2.081),
    "AE_UMa":  (r"E:\observe_raw_Analysis\AE UMa",       "lc", None, None),
    "YZ_Boo":  (r"E:\observe_raw_Analysis\YZbootis",     "lc", None, None),
}


def free_gb(drive="E:\\") -> float:
    return shutil.disk_usage(drive).free / 1e9


def log(msg: str):
    REPROCESS.mkdir(parents=True, exist_ok=True)
    with open(PROGRESS, "a", encoding="utf-8") as f:
        f.write(msg.rstrip() + "\n")
    print(msg, flush=True)


def is_done(target: str) -> bool:
    return PROGRESS.exists() and f"[DONE] {target}" in PROGRESS.read_text(encoding="utf-8")


def run_step0(target: str, raw_dir: str) -> Path:
    """APEX Step 0 headless -> reprocess/<target>/calibrated/. Skips if already present."""
    out = REPROCESS / target
    if list((out / "calibrated").rglob("pp_*.fit")):
        print(f"[{target}] Step 0 output exists — skip"); return out
    r = subprocess.run([str(VENV_PY), "-X", "utf8",
                        str(REPO / "scripts" / "_reprocess_step0.py"), raw_dir, target],
                       cwd=str(REPO))
    if r.returncode != 0:
        raise RuntimeError(f"Step 0 failed for {target}")
    return out


def reorg_per_object(target: str) -> Path:
    """Move this target's calibrated frames into reprocess/<target>/sci/ (one object)."""
    cal = REPROCESS / target / "calibrated"
    sci = REPROCESS / target / "sci"
    sci.mkdir(parents=True, exist_ok=True)
    stem = re.split(r"[_\d]", target)[0]  # NGC6811 -> NGC, M13 -> M ; refine per object below
    frames = list(cal.rglob(f"pp_{target}-*.fit")) + list(cal.rglob(f"pp_{target.replace('_','')}-*.fit"))
    if not frames:  # fall back: object token from filenames
        frames = [p for p in cal.rglob("pp_*.fit")
                  if target.replace("_", "").lower() in p.name.replace(" ", "").lower()]
    for p in frames:
        dst = sci / p.name
        if not dst.exists():
            shutil.move(str(p), str(dst))
    print(f"[{target}] reorg: {len(list(sci.glob('pp_*.fit')))} frames in sci/")
    return sci


def gen_config(target: str, sci: Path, ra: float, dec: float) -> Path:
    """Write reprocess/<target>/parameters.toml from the template with new io/target."""
    txt = TEMPLATE.read_text(encoding="utf-8")
    result = REPROCESS / target / "result"
    esc = lambda p: str(p).replace("\\", "\\\\")   # TOML needs doubled backslashes
    # lambda replacements avoid re.sub interpreting backslashes in the replacement
    txt = re.sub(r'(?m)^data_dir\s*=.*$', lambda m: f'data_dir = "{esc(sci)}"', txt)
    txt = re.sub(r'(?m)^result_dir\s*=.*$', lambda m: f'result_dir = "{esc(result)}"', txt)
    if ra is not None:
        txt = re.sub(r'(?m)^ra_deg\s*=.*$', lambda m: f'ra_deg = {ra}', txt)
        txt = re.sub(r'(?m)^dec_deg\s*=.*$', lambda m: f'dec_deg = {dec}', txt)
    cfg = REPROCESS / target / "parameters.toml"
    cfg.write_text(txt, encoding="utf-8")
    return cfg


def apex_run(cfg: Path, steps="1-7", mode="cmd"):
    r = subprocess.run([str(VENV_PY), "-m", "apex.cli", "run", "--mode", mode,
                        "--steps", steps, "--config", str(cfg), "--force"], cwd=str(REPO))
    if r.returncode != 0:
        raise RuntimeError(f"apex run {steps} failed ({cfg})")


def run_one(target: str) -> bool:
    raw, kind, ra, dec = TARGETS[target]
    log(f"\n### {target} ({kind}) — start, E: free {free_gb():.0f} GB")
    run_step0(target, raw)
    if kind == "lc":
        log(f"[DONE-Step0] {target} (lc, preprocess only)"); return True
    sci = reorg_per_object(target)
    if not list(sci.glob("pp_*.fit")):
        log(f"[SKIP] {target}: no calibrated frames after reorg"); return False
    cfg = gen_config(target, sci, ra, dec)
    apex_run(cfg, "1-7", "cmd")
    log(f"[STEP1-7] {target}: forced photometry done -> {REPROCESS/target/'result'}")
    # CMD 8-11 runners wired once built; then mark [DONE]
    log(f"[PENDING-CMD] {target}: steps 8-11 await runner build")
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="run a single target")
    a = ap.parse_args()
    log(f"# reprocess {time.strftime('%Y-%m-%d %H:%M')} — E: free {free_gb():.0f} GB")
    names = [a.only] if a.only else list(TARGETS)
    for t in names:
        if is_done(t):
            print(f"skip {t} (done)"); continue
        if free_gb() < MIN_FREE_GB:
            log(f"[STOP] E: free {free_gb():.0f} GB < {MIN_FREE_GB} — halt."); break
        try:
            run_one(t)
        except Exception as e:
            log(f"[ERROR] {t}: {type(e).__name__}: {e}")
            break
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
