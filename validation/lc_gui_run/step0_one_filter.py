"""Step 0 보정을 한 필터만 돌린다 (LC 완주 검증용).

    python step0_one_filter.py <raw_dir> <out_name> <filter>
"""
from __future__ import annotations
import os, sys, time
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from pathlib import Path

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from collections import Counter
from PyQt5.QtWidgets import QApplication
from apex.gui.workflow.step0_detector_calibration import _ScanWorker, _CalibrationWorker
from apex.analysis.calibration import CalibrationOptions

BIAS = r"E:\bias"
DARKS = r"E:\darks"


def main() -> int:
    raw_dir, name, want_filt = sys.argv[1], sys.argv[2], sys.argv[3]
    out = Path(r"E:\APEX_validation\reprocess") / name
    app = QApplication(sys.argv)
    t0 = time.time()

    sw = _ScanWorker([raw_dir, BIAS, DARKS])
    box = {}
    sw.done.connect(lambda fs: box.__setitem__("f", fs))
    sw.run()
    frames = box["f"]

    # 라이트만 필터로 거른다. bias/dark 는 필터가 없고, flat 은 매칭 단계에서
    # 해당 필터만 쓰이므로 그대로 둔다.
    kept = [f for f in frames
            if f.ftype != "light" or str(getattr(f, "filt", "")) == want_filt]
    n_light = sum(1 for f in kept if f.ftype == "light")
    print(f"[{name}] scan {len(frames)} → 사용 {len(kept)} "
          f"(light {want_filt} {n_light}장) by-type={dict(Counter(f.ftype for f in kept))}",
          flush=True)
    if not n_light:
        print("no light frames — abort", flush=True)
        return 2

    opts = CalibrationOptions(combine_method="median", pedestal_mode="none")
    w = _CalibrationWorker(kept, "__all__", out, opts)
    w.logline.connect(lambda m: print("  ", m, flush=True)
                      if ("master" in m or "warn" in m or "Calibrating" in m) else None)
    seen = {"n": 0}

    def _tick(done, total, msg):
        if done and done % 10 == 0 and done != seen["n"]:
            seen["n"] = done
            el = time.time() - t0
            print(f"   {done}/{total}  {el/60:.1f}분 경과  장당 {el/max(done,1):.1f}s",
                  flush=True)

    w.progress.connect(_tick)
    box2 = {}
    w.finished_ok.connect(lambda s: box2.__setitem__("s", s))
    w.failed.connect(lambda m: print("  [FAILED]", m, flush=True))
    w._run()
    s = box2.get("s", {})
    el = time.time() - t0
    print(f"[{name}] calibrated {s.get('n_calibrated')} lights in {el:.0f}s "
          f"({el/max(s.get('n_calibrated') or 1,1):.1f}s/장) -> {out}", flush=True)
    return 0 if s.get("n_calibrated") else 3


if __name__ == "__main__":
    raise SystemExit(main())
