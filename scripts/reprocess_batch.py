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


def _object_aliases(target: str) -> set[str]:
    """Filename object-name aliases for a target. A single raw night can image
    several objects (e.g. 20260611 has NGC6811 AND NGC3231), and calibrated
    frames are named by the object header (pp_<object>-NNNN-FILT.fit), which is
    often the long form: M3 -> 'messier3', M67 -> 'Messier67'. Match on that."""
    t = target.replace("_", "").lower()
    aliases = {t}
    m = re.match(r"m(\d+)$", t)
    if m:
        aliases |= {f"messier{m.group(1)}"}
    return aliases


_FILTER_SUFFIX = re.compile(r"[-_ ](g|r|i|z|u|b|v|ha|o3|oiii|sii|l)$", re.IGNORECASE)


def _norm_object(o: str) -> str:
    """Normalize a FITS OBJECT value to a comparable token: strip a trailing
    filter suffix ('M5-g' -> 'm5'), drop spaces, lowercase ('NGC 6811' handled
    by callers via aliases)."""
    o = _FILTER_SUFFIX.sub("", str(o).strip())
    return o.replace(" ", "").lower()


def _obj_token(name: str) -> str:
    """Filename object token: pp_messier3-0001-B.fit -> 'messier3'."""
    stem = name[3:] if name.lower().startswith("pp_") else name
    return stem.split("-")[0].strip().lower()


# Per-target filter allow-list. None = keep every filter. M3's SDSS g/i/r frames
# carry i-band fringe residuals that SEP detects as thousands of spurious sources
# (master catalogue blows up to 16714 vs ~2000), so M3 is restricted to the clean
# Johnson B/V/R set — the same system as M13/NGC 6811.
FILTER_ALLOW = {
    "M3": {"b", "v", "r"},
}
_FILT_RE = re.compile(r"-(\d+)-([A-Za-z0-9]+?)(?:_\d+)?\.fits?$", re.IGNORECASE)


def _frame_filter(name: str) -> str:
    """Filter token from a calibrated filename: pp_messier3-0001-B.fit -> 'b'."""
    m = _FILT_RE.search(name)
    return m.group(2).lower() if m else ""


def _frame_object(path: Path) -> str:
    """Normalized OBJECT header of a calibrated frame ('' if unreadable)."""
    try:
        from astropy.io import fits
        return _norm_object(fits.getheader(str(path)).get("OBJECT", ""))
    except Exception:
        return ""


def reorg_per_object(target: str) -> Path:
    """Move this target's calibrated frames into reprocess/<target>/sci/. A frame
    belongs to the target if its FITS OBJECT header (filter suffix stripped) OR
    its filename token matches the target aliases {m<n>, messier<n>} / NGC key.
    Header-driven so it works when filenames carry no object (M5: pp_0002.fit,
    OBJECT='M5-g'). Double-calibrated frames (raw already pp_ -> pp_pp_*) are
    excluded: their OBJECT still matches, but re-calibrating a calibrated frame
    is wrong, so only the clean single-pass set is kept."""
    cal = REPROCESS / target / "calibrated"
    sci = REPROCESS / target / "sci"
    sci.mkdir(parents=True, exist_ok=True)
    aliases = _object_aliases(target)
    allow = FILTER_ALLOW.get(target)          # None = all filters
    all_frames = [p for p in cal.rglob("pp_*.fit")
                  if not p.name.lower().startswith("pp_pp_")]
    frames, seen_obj, n_filt_drop = [], {}, 0
    for p in all_frames:
        obj = _frame_object(p)
        tok = _obj_token(p.name)
        key = obj or tok
        seen_obj[key] = seen_obj.get(key, 0) + 1
        if not (obj in aliases or tok in aliases):
            continue
        if allow is not None and _frame_filter(p.name) not in allow:
            n_filt_drop += 1
            continue
        frames.append(p)
    n_dpp = sum(1 for p in cal.rglob("pp_pp_*.fit"))
    if allow is not None:
        print(f"[{target}] filter allow-list {sorted(allow)}: dropped {n_filt_drop} off-filter frames")
    for p in frames:
        dst = sci / p.name
        if not dst.exists():
            shutil.move(str(p), str(dst))
    n = len(list(sci.glob("pp_*.fit")))
    excluded = sorted(k for k in seen_obj if k not in aliases)
    print(f"[{target}] reorg: {n} frames -> sci/ (aliases={sorted(aliases)}; "
          f"excluded_objects={excluded}; double_pp_skipped={n_dpp})")
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


def run_step10(cfg: Path) -> Path:
    """CMD Step 10 (ZP calibration -> CMD table) headless. Step 8/9 not needed
    on the forced-aperture path (Step 10 reads ID from Step 7 TSVs + Step 6
    master catalog). Returns the produced CMD table path."""
    r = subprocess.run([str(VENV_PY), "-X", "utf8",
                        str(REPO / "scripts" / "run_step10_headless.py"),
                        "--params", str(cfg)], cwd=str(REPO))
    if r.returncode != 0:
        raise RuntimeError(f"step10 failed ({cfg})")
    result = cfg.parent / "result"
    return result / "cmd_zeropoint" / "median_by_ID_filter_wide_cmd.csv"


def run_one(target: str) -> bool:
    raw, kind, ra, dec = TARGETS[target]
    log(f"\n### {target} ({kind}) — start, E: free {free_gb():.0f} GB")
    # Resumable: if the CMD is already built, skip BEFORE Step 0 so a re-run does
    # not re-calibrate (reorg has already emptied calibrated/ into sci/).
    new_cmd = REPROCESS / target / "result" / "cmd_zeropoint" / "median_by_ID_filter_wide_cmd.csv"
    if kind != "lc" and new_cmd.exists():
        log(f"[DONE] {target}: CMD table already present -> {new_cmd}")
        return True
    run_step0(target, raw)
    if kind == "lc":
        log(f"[DONE-Step0] {target} (lc, preprocess only)"); return True
    sci = reorg_per_object(target)
    if not list(sci.glob("pp_*.fit")):
        log(f"[SKIP] {target}: no calibrated frames after reorg"); return False
    cfg = gen_config(target, sci, ra, dec)
    apex_run(cfg, "1-7", "cmd")
    log(f"[STEP1-7] {target}: forced photometry done -> {REPROCESS/target/'result'}")
    new_cmd = run_step10(cfg)
    if not new_cmd.exists():
        log(f"[SKIP] {target}: step10 produced no CMD table"); return False
    # NOTE: no comparison against E:\observed_Analysis (AIPPI-preprocessed) — that
    # is the author's own tool, not an independent reference, so it is not a
    # validation. Independent accuracy checks (IRAF / ccdproc / PS1 / Gaia) are
    # done per-figure in validation/paper/. Filter systems also differ per target
    # (M13 Johnson vs its AIPPI SDSS run), so a band cross-match is not meaningful.
    log(f"[DONE] {target}: full-APEX CMD complete ({new_cmd}); isochrone step 12 left to user")
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="run a single target")
    ap.add_argument("--lc", action="store_true",
                    help="also run LC targets (Step 0 preprocess). Deferred by "
                         "default: LC calibrated output is ~180 GB and needs "
                         "E:\\observed_Analysis deleted first.")
    a = ap.parse_args()
    log(f"# reprocess {time.strftime('%Y-%m-%d %H:%M')} — E: free {free_gb():.0f} GB")
    names = [a.only] if a.only else list(TARGETS)
    for t in names:
        if is_done(t):
            print(f"skip {t} (done)"); continue
        if TARGETS[t][1] == "lc" and not a.lc:
            print(f"skip {t} (lc deferred — pass --lc to run)"); continue
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
