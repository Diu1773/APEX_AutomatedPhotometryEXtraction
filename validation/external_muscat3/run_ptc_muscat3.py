# -*- coding: utf-8 -*-
"""세 카메라를 같은 조건으로 재서 헤더값과 견준다."""
import glob, os, sys
sys.path.insert(0, r"C:\Users\KNUE\galaxy-jobs\apex-lco-m67-muscat3-20260907-001\src")
from apex.analysis.detector_ptc import characterize_detector, scan_calibration_frames, measure_levels
from apex.utils.io_utils import read_fits_header

base = r"C:\Users\KNUE\galaxy-jobs\apex-lco-m67-muscat3-20260907-001\inputs\calib"
CAMS = (("gp", "ep04"), ("rp", "ep02"), ("ip", "ep03"))
TOL = 0.06          # 하늘 플랫은 밝기가 계속 흘러서 2 % 로는 짝이 안 잡힌다

print(f"{'카메라':<8}{'밴드':<5}{'헤더 GAIN':>10}{'PTC gain':>20}{'비':>8}"
      f"{'헤더 RN':>9}{'PTC RN':>9}{'쌍':>4}{'R^2':>8}")
print('-' * 82)
for band, cam in CAMS:
    bias = sorted(glob.glob(os.path.join(base, "bias", band, "*.fits.fz")))
    flat = sorted(glob.glob(os.path.join(base, "skyflat", band, "*.fits.fz")))
    h = read_fits_header(flat[0])
    hg, hr = float(h.get("GAIN", 0)), float(h.get("RDNOISE", 0))
    try:
        r = characterize_detector(bias, flat, signal_floor=0.0,
                                  box_radius=300, tolerance=TOL)
    except Exception as e:
        print(f"{cam:<8}{band:<5}{hg:>10.3f}{'실패: '+str(e)[:40]:>20}")
        continue
    print(f"{cam:<8}{band:<5}{hg:>10.3f}"
          f"{r.gain_eff:>13.4f} +-{r.gain_eff_err:<6.4f}"
          f"{r.gain_eff/hg:>8.3f}{hr:>9.1f}{r.read_noise_eff:>9.2f}"
          f"{r.n_pairs:>4}{r.r_squared:>8.4f}")
    for line in r.log:
        if 'WARNING' in line:
            print(f"        경고: {line.split('WARNING:')[-1].strip()[:70]}")
