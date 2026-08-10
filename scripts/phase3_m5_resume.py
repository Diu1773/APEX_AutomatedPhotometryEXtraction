"""Finish M5 in the Phase 3 tree, starting from its preserved calibrated frames.

M5 is the one target whose Step 0 cannot be re-run: its raw light frames are no
longer on this machine. Verified three ways before concluding —
`E:\\observe_DSY\\M5\\M5_20250308\\light5\\` is empty (zero files of any
extension), the whole night directory holds 24 files and all of them are darks
or flats, and `E:\\observed_Analysis\\M5\\light\\` has no FITS either. The
batch therefore aborted with "no light frames" and stopped before the remaining
steps.

Nothing scientific is lost by starting from the preserved calibrated frames.
Step 0 is proven bit-identical run to run — O1's gate compared 30/30 M67 frames
at max |delta| = 0, and this clean run reproduced 21/21, 30/30 and 15/15 pixel
data exactly on the three targets that did re-run. Re-running M5's Step 0 would
write the same bytes it already has. What cannot be claimed for M5 is that its
calibrated frames were *regenerated* in this run, and the results document says
so rather than letting a copied file pass as a recomputed one.

So: copy the preserved frames into the Phase 3 tree, then run exactly what the
batch would have run — Steps 1-7 and CMD Step 10 — through the batch's own
functions, so the code path is identical to the other four targets.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import scripts.reprocess_batch as batch  # noqa: E402

TARGET = "M5"
PHASE3 = Path(r"E:\APEX_validation\phase3")
BASELINE_SCI = Path(r"E:\APEX_validation\reprocess\M5\sci")

# Rebind the batch's output root exactly as `--out` does, so gen_config and
# run_step10 write into the Phase 3 tree and not over the baseline.
batch.REPROCESS = PHASE3
batch.PROGRESS = PHASE3 / "PROGRESS.md"
batch.OUT_DRIVE = str(PHASE3.anchor) or batch.OUT_DRIVE

sci = PHASE3 / TARGET / "sci"
sci.mkdir(parents=True, exist_ok=True)

frames = sorted(BASELINE_SCI.glob("pp_*.fit"))
if not frames:
    raise SystemExit(f"no preserved calibrated frames in {BASELINE_SCI}")

copied = 0
for src in frames:
    dst = sci / src.name
    if dst.exists() and dst.stat().st_size == src.stat().st_size:
        continue
    shutil.copy2(src, dst)
    copied += 1
print(f"[M5] calibrated frames in place: {len(frames)} ({copied} copied, "
      f"rest already present) — Step 0 INHERITED, not recomputed", flush=True)

_raw, _kind, ra, dec = batch.TARGETS[TARGET]
cfg = batch.gen_config(TARGET, sci, ra, dec)
print(f"[M5] config: {cfg}", flush=True)

batch.apex_run(cfg, "1-7", "cmd")
print(f"[M5] Steps 1-7 done -> {PHASE3/TARGET/'result'}", flush=True)

cmd_table = batch.run_step10(cfg)
if not cmd_table.exists():
    raise SystemExit("[M5] step10 produced no CMD table")
print(f"[M5] DONE — CMD table: {cmd_table}", flush=True)
