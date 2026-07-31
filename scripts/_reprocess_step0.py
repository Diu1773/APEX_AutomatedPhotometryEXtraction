"""APEX Step 0 headless for one target: scan <raw_dir> + E:\\bias + E:\\darks,
build masters, calibrate lights -> E:\\APEX_validation\\reprocess\\<name>\\.
Calls the same Qt-free core the GUI Step-0 window drives (parity, no Qt).
Called by scripts/reprocess_batch.py.

    .venv-deploy\\Scripts\\python.exe scripts/_reprocess_step0.py <raw_dir> <name>
"""
from __future__ import annotations
import sys, time
from pathlib import Path
REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from collections import Counter
from apex.analysis import calibration_scan as scan
from apex.analysis.calibration import CalibrationOptions
from apex.analysis.calibration_run import ALL_NIGHTS, run_calibration

BIAS = r"E:\bias"
DARKS = r"E:\darks"
# Korea (KASI/KNUE sites): the observing-night split needs a local reference
# when a frame header carries no SITELONG.
TZ_OFFSET_HOURS = 9.0


def main() -> int:
    raw_dir = sys.argv[1]
    name = sys.argv[2]
    out = Path(r"E:\APEX_validation\reprocess") / name
    t0 = time.time()
    fr = []
    for root in (raw_dir, BIAS, DARKS):
        fr.extend(scan.scan_folder(root, tz_offset_hours=TZ_OFFSET_HOURS,
                                   warn=lambda m: print("  ", m, flush=True)))
    print(f"[{name}] scanned {len(fr)}  by-type={dict(Counter(f.ftype for f in fr))}  "
          f"nights={scan.nights(fr)}", flush=True)
    if not any(f.ftype == "light" for f in fr):
        print(f"[{name}] no light frames — abort", flush=True)
        return 2
    # cosmetic(우주선·핫픽셀 제거)은 CalibrationOptions 기본값(True)을 따른다.
    # 예전에는 여기서 명시적으로 껐는데, 그 탓에 CMOS 프레임의 1픽셀 스파이크가
    # ePSF 참조별로 뽑혀 PSF 가 무너졌다(2026-07-29 M67/QHY600 실측).
    opts = CalibrationOptions(combine_method="median", pedestal_mode="none")

    def _log(m):
        if "master" in m or "warn" in m or "Calibrating" in m:
            print("  ", m, flush=True)

    try:
        s = run_calibration(fr, ALL_NIGHTS, out, opts, log=_log)
    except Exception as exc:
        print("  [FAILED]", exc, flush=True)
        return 3
    print(f"[{name}] calibrated {s.get('n_calibrated')} lights / {s.get('n_nights')} nights "
          f"in {time.time()-t0:.0f}s -> {out}", flush=True)
    if s.get("n_temp_mismatch"):
        print(f"[{name}] WARNING: {s['n_temp_mismatch']} light(s) had no dark within "
              f"{opts.temp_match_tol_c:g}°C", flush=True)
    return 0 if s.get("n_calibrated") else 3


if __name__ == "__main__":
    raise SystemExit(main())
