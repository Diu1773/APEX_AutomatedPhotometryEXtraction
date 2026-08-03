# -*- coding: utf-8 -*-
"""Per-step APEX-vs-ccdproc calibration cross-check — data generator.

Recreates (and this time keeps in the repo) the generator behind
``data/calib_crosscheck_ngc6811.json``: the original script was lost after its
one 2026-07-09 run, leaving the figure unable to access per-step difference
*arrays*.  This rerun therefore also writes ``data/calib_crosscheck_maps.npz``
holding a max-pooled |APEX - ccdproc| map per step, which the reworked
fig12 panel (a) draws (user 2026-08-03: per-step plots instead of the
log-lollipop whose five bit-identical points sat on a fake 1e-9 floor).

Inputs (all real Moravian C3-61000 frames, 2x2 binning, night 2026-06-11 for
darks/flats/light; the bias library is the observatory's -10 C set):
  light  E:\\observe_raw_Analysis\\20260611\\NGc6811\\NGC6811-0001-B.fit (60 s)
  darks  same night, header EXPTIME == 60 s, first 8 by name
  flats  same night, flat_B_0001..0005
  bias   E:\\bias\\bias-0091..0098
Every input file name is recorded in the JSON so the subset is reproducible.

Step semantics (mirrors the manuscript 3.2 wording):
  * master bias / dark / flat  — each side builds its own master from the same
    raw frames with its own primitives (APEX ``apex.analysis.calibration`` in a
    float32 stack; ccdproc ``combine`` forced to the same stack dtype).
  * bias/dark/flat application — isolated op: both sides apply the SAME master
    to the SAME input, so any difference is the arithmetic itself.
  * full pipeline — true end-to-end: APEX ``calibrate_light_file`` vs the
    chained ccdproc reduction, each from its own masters.

Run:  .venv-deploy\\Scripts\\python -X utf8 validation\\paper\\calib_crosscheck_ngc6811.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData
import ccdproc

from apex.analysis import calibration as cal

NIGHT = Path(r"E:\observe_raw_Analysis\20260611")
LIGHT = NIGHT / "NGc6811" / "NGC6811-0001-B.fit"
BIAS_DIR = Path(r"E:\bias")
DATA = REPO / "validation" / "paper" / "data"
POOL = 8      # |delta| 최대값 풀링 크기 (표시용 지도)


def pick_darks_60s(n: int = 8) -> list[Path]:
    out = []
    for p in sorted((NIGHT / "dark").glob("*.fit")):
        try:
            if float(fits.getheader(p).get("EXPTIME", -1)) == 60.0:
                out.append(p)
        except OSError:
            continue
        if len(out) == n:
            break
    if len(out) != n:
        raise SystemExit(f"60 s dark {len(out)}/{n}")
    return out


def maxpool_abs(d: np.ndarray, k: int = POOL) -> np.ndarray:
    a = np.abs(np.nan_to_num(d, nan=0.0))
    h, w = a.shape
    a = a[: h - h % k, : w - w % k]
    return a.reshape(h // k, k, w // k, k).max(axis=(1, 3)).astype(np.float32)


def stats(d: np.ndarray) -> dict:
    f = d[np.isfinite(d)]
    med = float(np.median(f))
    return {"max_abs": float(np.max(np.abs(f))),
            "robust_sigma": float(1.4826 * np.median(np.abs(f - med)))}


def main() -> int:
    biases = sorted(BIAS_DIR.glob("bias-*.fit"))[:8]
    darks = pick_darks_60s(8)
    flats = [NIGHT / "flats" / f"flat_B_{i:04d}.fit" for i in range(1, 6)]
    for p in [LIGHT, *biases, *darks, *flats]:
        if not p.exists():
            raise SystemExit(f"없음: {p}")

    # cosmetic(L.A.Cosmic+핫픽셀 수리)은 끈다. 2026-07 부터 기본 켬이라 켠 채로
    # 비교하면 수리된 15만 픽셀이 산술 차이로 잡힌다(최대 1.2e4 DN — 실측).
    # 이 비교의 주장은 bias/dark/flat *산술*의 동등성이고, cosmetic 단계는
    # 3.2 의 주입 시험이 따로 검증한다.
    opts = cal.CalibrationOptions(combine_method="median", pedestal_mode="none",
                                  cosmetic_enable=False)

    # ── APEX 쪽 마스터 ──
    mbias_a, _ = cal.build_master_bias(biases, opts)
    mdark_a, dexp, _ = cal.build_master_dark(darks, opts, master_bias=mbias_a)
    mflat_a, _ = cal.build_master_flat(flats, opts, master_bias=mbias_a,
                                       master_dark=mdark_a, dark_exp=dexp)

    # ── ccdproc 쪽 마스터 (같은 raw, ccdproc 프리미티브, 같은 스택 dtype) ──
    def ccd(p: Path) -> CCDData:
        d, _ = cal.load_frame(p, opts)          # 같은 로더(오버스캔 규약 공유)
        return CCDData(np.asarray(d, dtype=np.float32), unit="adu")

    mbias_c = ccdproc.combine([ccd(p) for p in biases], method="median",
                              dtype=np.float32)
    dark_ccs = [ccdproc.subtract_bias(ccd(p), mbias_c) for p in darks]
    mdark_c = ccdproc.combine(dark_ccs, method="median", dtype=np.float32)
    flat_ccs = []
    for p in flats:
        f = ccdproc.subtract_bias(ccd(p), mbias_c)
        f = ccdproc.subtract_dark(f, mdark_c, scale=True,
                                  dark_exposure=dexp * u.s,
                                  data_exposure=float(fits.getheader(p)["EXPTIME"]) * u.s)
        flat_ccs.append(f.divide(float(np.nanmedian(f.data)) * u.dimensionless_unscaled))
    mflat_c = ccdproc.combine(flat_ccs, method="median", dtype=np.float32)
    mflat_c = mflat_c.divide(float(np.nanmedian(mflat_c.data)) * u.dimensionless_unscaled)

    # ── 적용 단계 (고립 연산: 양쪽에 같은 입력·같은 마스터) ──
    light64, lhdr = cal.load_frame(LIGHT, opts)
    lexp = float(lhdr["EXPTIME"])
    L = CCDData(np.asarray(light64, dtype=np.float64), unit="adu")
    M_BIAS = CCDData(np.asarray(mbias_a, dtype=np.float64), unit="adu")
    M_DARK = CCDData(np.asarray(mdark_a, dtype=np.float64), unit="adu")

    apex_b = np.asarray(light64, dtype=np.float64) - mbias_a
    cc_b = ccdproc.subtract_bias(L, M_BIAS).data
    apex_d = apex_b - mdark_a * (lexp / dexp)
    cc_d = ccdproc.subtract_dark(CCDData(apex_b, unit="adu"), M_DARK, scale=True,
                                 dark_exposure=dexp * u.s,
                                 data_exposure=lexp * u.s).data
    safe = np.where(mflat_a < opts.flat_min, 1.0, mflat_a)
    apex_f = apex_d / safe
    cc_f = ccdproc.flat_correct(
        CCDData(apex_d, unit="adu"),
        CCDData(np.asarray(safe, dtype=np.float64), unit="adu"),
        norm_value=1.0).data

    # ── 완전 파이프라인 (각자 자기 마스터로 끝까지) ──
    apex_full, _h, _qc = cal.calibrate_light_file(
        LIGHT, opts, master_bias=mbias_a, master_dark=mdark_a,
        dark_exp=dexp, master_flat=mflat_a)
    lc = ccd(LIGHT)
    lc = ccdproc.subtract_bias(lc, mbias_c)
    lc = ccdproc.subtract_dark(lc, mdark_c, scale=True,
                               dark_exposure=dexp * u.s, data_exposure=lexp * u.s)
    safe_c = CCDData(np.where(mflat_c.data < opts.flat_min, 1.0,
                              mflat_c.data).astype(np.float32), unit="adu")
    cc_full = ccdproc.flat_correct(lc, safe_c, norm_value=1.0).data

    diffs = {
        "master_bias": mbias_a - mbias_c.data.astype(np.float64),
        "master_dark": mdark_a - mdark_c.data.astype(np.float64),
        "master_flat": mflat_a - mflat_c.data.astype(np.float64),
        "bias_subtract": apex_b - cc_b,
        "dark_subtract": apex_d - cc_d,
        "flat_correct": apex_f - cc_f,
        "full_pipeline": apex_full.astype(np.float64) - cc_full.astype(np.float64),
    }

    steps = {k: stats(v) for k, v in diffs.items()}
    maps = {k: maxpool_abs(v) for k, v in diffs.items()}

    out = {
        "light": LIGHT.name, "light_exp_s": lexp,
        "light_median": float(np.nanmedian(light64)),
        "reference": f"astropy ccdproc {ccdproc.__version__}",
        "night": "2026-06-11", "camera": "Moravian C3-61000 (2x2)",
        "n_bias": len(biases), "n_dark": len(darks), "n_flat": len(flats),
        "dark_exp": dexp,
        "inputs": {"bias": [p.name for p in biases],
                   "dark": [p.name for p in darks],
                   "flat": [p.name for p in flats]},
        "steps": steps,
    }
    (DATA / "calib_crosscheck_ngc6811.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")
    np.savez_compressed(DATA / "calib_crosscheck_maps.npz",
                        pool=POOL, shape=np.array(light64.shape), **maps)

    for k in diffs:
        print(f"  {k:14s} max|Δ|={steps[k]['max_abs']:.2e}  σ={steps[k]['robust_sigma']:.2e}")
    print("wrote JSON + maps.npz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
