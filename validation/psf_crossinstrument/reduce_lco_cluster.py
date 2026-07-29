# -*- coding: utf-8 -*-
"""LCO 성단 raw 다중프레임을 APEX식으로 보정해 sci/ 로 낸다 (QHY600 단일앰프).

BANZAI QHY600 레시피 (fig13 에서 e91 대조로 검증된 것):
    e- = (raw*GAIN - BIASLVL - master_bias - master_dark*t_ratio) / master_flat

BIASLVL 은 raw 에 없고 e91 헤더에만 있다. 같은 마스터를 쓰는 프레임은 값이
같으므로(M67 2024-10-29 B·V 둘 다 376.363) 받은 e91 에서 읽어 전 프레임에 쓴다.
받은 e91 1장으로 보정 결과를 대조해 어긋나면 즉시 멈춘다.

출력 단위가 전자이므로 APEX 설정은 gain=1.0 을 쓴다.

실행:
    .venv-deploy\\Scripts\\python -X utf8 validation/psf_crossinstrument/reduce_lco_cluster.py m67_lco
    .venv-deploy\\Scripts\\python -X utf8 validation/psf_crossinstrument/reduce_lco_cluster.py m45_wide
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

BASE = Path(r"E:\APEX_validation\psf_crossinstrument")

# 대상별 출력 파일명 접두 (APEX 는 pp_<object>-NNNN-FILT.fit 형태를 쓴다)
PREFIX = {"m67_lco": "m67lco", "m45_wide": "m45lco"}

KEEP = [
    "DATE-OBS", "EXPTIME", "FILTER", "OBJECT", "RA", "DEC", "CAT-RA", "CAT-DEC",
    "AIRMASS", "INSTRUME", "TELESCOP", "SITEID", "LATITUDE", "LONGITUD", "HEIGHT",
    "PIXSCALE", "RDNOISE", "MJD-OBS", "DAY-OBS", "OBSTYPE",
]


def _sec(s: str) -> tuple[int, int, int, int] | None:
    """FITS 섹션 '[x1:x2,y1:y2]'(1-based 양끝 포함) → 0-based 슬라이스 경계."""
    if not s or str(s).strip().upper() in ("N/A", "UNKNOWN", ""):
        return None
    m = re.findall(r"-?\d+", str(s))
    if len(m) != 4:
        return None
    x1, x2, y1, y2 = (int(v) for v in m)
    return x1 - 1, x2, y1 - 1, y2


def _sci(path: Path):
    """2-D 데이터와 헤더. float32 로 읽는다 — M45 풀프레임(61.7 MP)에서 마스터
    4장을 float64 로 들면 2 GB 라 다른 실행과 겹치면 메모리가 모자란다. 값은
    전자 단위 10^4 규모여서 float32(유효 7자리)로 0.001 e- 까지 남는다."""
    with fits.open(path) as h:
        for hd in h:
            if hd.data is not None and hd.data.ndim == 2:
                return hd.data.astype(np.float32), hd.header.copy()
    raise RuntimeError(f"no 2-D HDU: {path}")


def main(key: str) -> None:
    d = BASE / key
    raw_dir, mst_dir, ref_dir = d / "raw", d / "masters", d / "banzai_ref"
    out_dir = d / "sci"
    out_dir.mkdir(parents=True, exist_ok=True)

    # 마스터
    mbias, _ = _sci(next(mst_dir.glob("*bias*.fits.fz")))
    mdark, mdh = _sci(next(mst_dir.glob("*dark*.fits.fz")))
    dexp = float(mdh.get("EXPTIME", 300.0))
    flats = {}
    for p in mst_dir.glob("*skyflat*"):
        f = re.search(r"-([A-Za-z0-9]+)\.fits", p.name)
        if f:
            flats[f.group(1)] = _sci(p)[0]
    print(f"masters: bias·dark(t={dexp:.0f}s)·flat{sorted(flats)}")

    # BIASLVL — 받은 e91 중 아무거나 (같은 마스터면 동일)
    ref_frames = sorted(ref_dir.glob("*e91*"))
    _rd, rh = _sci(ref_frames[0])
    biaslvl = float(rh["BIASLVL"])
    lvls = {float(_sci(p)[1]["BIASLVL"]) for p in ref_frames}
    print(f"BIASLVL = {biaslvl:.3f} e-  (받은 e91 {len(ref_frames)}장의 값: {sorted(lvls)})")
    if max(lvls) - min(lvls) > 1.0:
        raise SystemExit("BIASLVL 이 프레임마다 1 e- 넘게 다르다 — 상수 적용 불가")

    # 프레임별 보정
    raws = sorted(raw_dir.glob("*e00*"))
    counter: dict[str, int] = {}
    checked = False
    for p in raws:
        raw, H = _sci(p)
        filt = str(H["FILTER"]).strip()
        if filt not in flats:
            print(f"  SKIP {p.name}: flat 없음 ({filt})")
            continue
        gain = float(H["GAIN"])
        lexp = float(H["EXPTIME"])
        mf = flats[filt]

        # 읽기 모드에 따라 오버스캔이 있다 — QHY600 full_frame 은 상단 30행에
        # BIASSEC 이 있고 TRIMSEC 으로 잘라야 마스터와 크기가 맞는다.
        # central30x30 크롭 모드는 BIASSEC='N/A' 라 이 블록을 건너뛴다.
        bs, ts = _sec(H.get("BIASSEC", "")), _sec(H.get("TRIMSEC", ""))
        over = 0.0
        if bs is not None:
            x1, x2, y1, y2 = bs
            over = float(np.median(raw[y1:y2, x1:x2]))
            raw = raw - np.float32(over)
        if ts is not None:
            x1, x2, y1, y2 = ts
            raw = raw[y1:y2, x1:x2]
        if raw.shape != mbias.shape:
            raise SystemExit(
                f"{p.name}: 트림 후 {raw.shape} 가 마스터 {mbias.shape} 와 다르다"
            )

        apex = (raw * gain - biaslvl - mbias - mdark * (lexp / dexp)) / np.where(
            np.abs(mf) < 1e-6, np.nan, mf
        )

        # 첫 프레임 중 e91 이 있는 것으로 보정 산술 대조
        ref = ref_dir / p.name.replace("-e00", "-e91")
        if ref.exists():
            e91, _ = _sci(ref)
            diff = apex - e91
            fin = np.isfinite(diff)
            med = float(np.median(diff[fin]))
            mad = float(1.4826 * np.median(np.abs(diff[fin] - med)))
            sky = float(np.median(e91[fin]))
            print(f"  [대조] {p.name[-18:]} vs BANZAI: Δmed={med:+.4f} e- "
                  f"robustσ={mad:.4f} (sky {sky:.1f} e-)")
            if abs(med) > 1.0 or mad > 1.0:
                raise SystemExit("BANZAI 대조 실패 — 보정 산술을 확인할 것")
            checked = True

        n_nan = int(np.sum(~np.isfinite(apex)))
        apex[~np.isfinite(apex)] = np.nanmedian(apex)
        ov_txt = f" over={over:.1f}ADU" if bs is not None else ""

        counter[filt] = counter.get(filt, 0) + 1
        name = f"pp_{PREFIX[key]}-{counter[filt]:04d}-{filt}.fit"
        hdr = fits.Header()
        for k in KEEP:
            if k in H:
                hdr[k] = H[k]
        hdr["BUNIT"] = ("electron", "calibrated frame in electrons")
        hdr["EGAIN"] = (1.0, "data already in electrons")
        hdr["SATURATE"] = (float(rh.get("SATURATE", 46200.0)), "e-, from BANZAI e91")
        hdr["CALSRC"] = ("APEX-arith vs BANZAI", "calibration provenance")
        fits.PrimaryHDU(apex.astype(np.float32), hdr).writeto(out_dir / name, overwrite=True)
        print(f"  {name}  median={np.nanmedian(apex):7.2f} e-  nan_filled={n_nan}{ov_txt}")

    if not checked:
        raise SystemExit("e91 대조를 한 번도 못 했다 — banzai_ref 파일명 확인")
    print(f"\n{key}: {sum(counter.values())}장 → {out_dir}  ({counter})")


if __name__ == "__main__":
    for k in (sys.argv[1:] or ["m67_lco"]):
        main(k)
