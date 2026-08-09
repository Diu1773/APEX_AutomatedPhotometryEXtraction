"""Is the shared chain limited by the code or by the volume it reads from?

Every stage of Steps 1-7 reads the full frame set — 21 x 61 MB for the NGC 6811
fixture — and none of them scale: detect gains 10 % at 2 workers and then
flattens, wcs is fastest at 1, forced photometry is 5.5x *slower* at 12. On 16
logical cores that pattern says the bottleneck is not CPU.

The fixture lives on a USB-attached portable SSD. USB storage has a shallow
command queue and a single pipe, so concurrent readers do not buy bandwidth,
they add per-request overhead. If that is the mechanism, the same run on the
internal NVMe should scale where the USB volume does not — and the per-stage
ceilings are then a property of *this storage*, not of APEX.

Runs Step 7 alone (the worst-affected stage) at 1 and 12 workers, once with
everything on the USB volume and once with everything copied to the internal
disk. Steps 1-6 outputs are copied, not recomputed, so both sides start from
byte-identical inputs.

Usage:  python -X utf8 benchmark/storage_scaling.py [--workers 1,12]
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from apex.benchmark import resources  # noqa: E402

CFG = Path(r"E:\APEX_validation\reprocess\NGC6811\apex_config_20260807.json")
USB_SCI = Path(r"E:\APEX_validation\reprocess\NGC6811\sci")
USB_TREE = Path(r"E:\APEX_validation\bench\arm_serial_1")
LOCAL_ROOT = Path(r"C:\apex_bench")
OUT = REPO / "benchmark" / "perf" / "20260809" / "storage_scaling.json"

PHASES = ["load_s", "wcs_fwhm_s", "phot_apcorr_s", "save_s", "elapsed_s"]


def stage_local() -> tuple[Path, Path]:
    """Copy frames and the Steps 1-6 outputs to the internal disk once."""
    sci, tree = LOCAL_ROOT / "sci", LOCAL_ROOT / "result"
    if not sci.exists():
        print(f"copying frames -> {sci} …", flush=True)
        shutil.copytree(USB_SCI, sci)
    if tree.exists():
        shutil.rmtree(tree)
    print(f"copying Steps 1-6 outputs -> {tree} …", flush=True)
    shutil.copytree(USB_TREE, tree,
                    ignore=shutil.ignore_patterns("step7_forced_phot"))
    return sci, tree


def run_step7(label: str, sci: Path, tree: Path, workers: int) -> dict:
    import pandas as pd

    env = {**os.environ, "APEX_MAX_WORKERS": str(workers)}
    cmd = [sys.executable, "-X", "utf8", "-m", "apex.cli", "run", "--mode", "cmd",
           "--config", str(CFG), "--steps", "7", "--data-dir", str(sci),
           "--result-dir", str(tree), "--force"]
    started = time.perf_counter()
    rc = subprocess.call(cmd, cwd=str(REPO), env=env)
    wall = time.perf_counter() - started

    index = pd.read_csv(tree / "step7_forced_phot" / "photometry_index.csv")
    ok = index[index["status"] == "ok"]
    record = {"storage": label, "workers": workers, "returncode": rc,
              "wall_s": round(wall, 1), "frames": int(len(ok)),
              "machine": resources.machine_state()}
    for phase in PHASES:
        if phase in ok:
            record[phase] = round(float(ok[phase].sum()), 1)
    return record


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", default="1,12")
    ap.add_argument("--keep", action="store_true",
                    help="keep the copied frames on the internal disk")
    args = ap.parse_args()
    counts = [int(w) for w in args.workers.split(",")]

    records = []
    try:
        local_sci, local_tree = stage_local()
        for workers in counts:
            for label, sci, tree in (("usb", USB_SCI, USB_TREE),
                                     ("nvme", local_sci, local_tree)):
                r = run_step7(label, sci, tree, workers)
                records.append(r)
                print(f"{label:>5} w={workers:<3} wall={r['wall_s']:7.1f}s  "
                      + "  ".join(f"{p}={r.get(p, 0):.0f}s" for p in PHASES),
                      flush=True)
                resources.save_metrics(
                    {"label": "storage_scaling", "runs": records}, OUT)
    finally:
        if not args.keep and LOCAL_ROOT.exists():
            print(f"removing {LOCAL_ROOT} …", flush=True)
            shutil.rmtree(LOCAL_ROOT, ignore_errors=True)

    print(f"\n{'storage':>8}{'w':>4}{'wall':>10}" +
          "".join(f"{p:>15}" for p in PHASES))
    for r in records:
        print(f"{r['storage']:>8}{r['workers']:>4}{r['wall_s']:>9.1f}s" +
              "".join(f"{r.get(p, 0):>14.1f}s" for p in PHASES))
    print(f"\nsaved -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
