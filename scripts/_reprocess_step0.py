"""APEX Step 0 headless for one target: scan <raw_dir> + E:\\bias + E:\\darks,
build masters, calibrate lights -> E:\\APEX_validation\\reprocess\\<name>\\.
Drives the real GUI workers (parity). Called by scripts/reprocess_batch.py.

    .venv-deploy\\Scripts\\python.exe scripts/_reprocess_step0.py <raw_dir> <name>
"""
from __future__ import annotations
import os, sys, time
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from pathlib import Path
REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from collections import Counter
from PyQt5.QtWidgets import QApplication
from apex.gui.workflow.step0_detector_calibration import _ScanWorker, _CalibrationWorker
from apex.analysis.calibration import CalibrationOptions
from apex.analysis import calibration_scan as scan

BIAS = r"E:\bias"
DARKS = r"E:\darks"


def main() -> int:
    raw_dir = sys.argv[1]
    name = sys.argv[2]
    out = Path(r"E:\APEX_validation\reprocess") / name
    app = QApplication(sys.argv)
    t0 = time.time()
    sw = _ScanWorker([raw_dir, BIAS, DARKS])
    box = {}
    sw.done.connect(lambda fs: box.__setitem__("f", fs))
    sw.run()
    fr = box["f"]
    print(f"[{name}] scanned {len(fr)}  by-type={dict(Counter(f.ftype for f in fr))}  "
          f"nights={scan.nights(fr)}", flush=True)
    if not any(f.ftype == "light" for f in fr):
        print(f"[{name}] no light frames — abort", flush=True)
        return 2
    opts = CalibrationOptions(combine_method="median", pedestal_mode="none",
                              cosmetic_enable=False)
    w = _CalibrationWorker(fr, "__all__", out, opts)
    w.logline.connect(lambda m: print("  ", m, flush=True)
                      if ("master" in m or "warn" in m or "Calibrating" in m) else None)
    box2 = {}
    w.finished_ok.connect(lambda s: box2.__setitem__("s", s))
    w.failed.connect(lambda m: print("  [FAILED]", m, flush=True))
    w._run()
    s = box2.get("s", {})
    print(f"[{name}] calibrated {s.get('n_calibrated')} lights / {s.get('n_nights')} nights "
          f"in {time.time()-t0:.0f}s -> {out}", flush=True)
    return 0 if s.get("n_calibrated") else 3


if __name__ == "__main__":
    raise SystemExit(main())
