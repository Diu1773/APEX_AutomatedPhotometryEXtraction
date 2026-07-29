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

# 레포 루트를 경로에 넣어 apex 패키지(우주선 제거)를 쓴다.
# resolve() 를 쓰면 안 된다 — validation/ 은 E:\APEX_validation_output 으로의
# 정션이라 링크를 따라가면 레포 밖(E:\)으로 나가 apex 를 못 찾는다.
# absolute() 는 링크를 따라가지 않는다.
_REPO = Path(__file__).absolute().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

BASE = Path(r"E:\APEX_validation\psf_crossinstrument")

# 대상별 출력 파일명 접두 (APEX 는 pp_<object>-NNNN-FILT.fit 형태를 쓴다)
PREFIX = {"m67_lco": "m67lco", "m45_wide": "m45lco", "m67_ubv": "m67ubv"}

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


def main(key: str, clean_cr: bool = True) -> None:
    from apex.analysis.cosmetic import (
        HAS_ASTROSCRAPPY,
        clean_frame,
        hot_pixel_mask,
    )

    if clean_cr and not HAS_ASTROSCRAPPY:
        raise SystemExit("astroscrappy 가 없다 — --no-cr 로 끄거나 설치할 것")

    d = BASE / key
    raw_dir, mst_dir, ref_dir = d / "raw", d / "masters", d / "banzai_ref"
    out_dir = d / "sci"
    out_dir.mkdir(parents=True, exist_ok=True)

    # 마스터
    mbias, _ = _sci(next(mst_dir.glob("*bias*.fits.fz")))
    mdark, mdh = _sci(next(mst_dir.glob("*dark*.fits.fz")))
    dexp = float(mdh.get("EXPTIME", 300.0))
    # 마스터 플랫의 필터 접미사 구분자가 세대마다 다르다:
    #   2024 QHY600 : ...bin1x1-B.fits.fz   (하이픈)
    #   2015 SBIG   : ..._bin1x1_U.fits.fz  (언더스코어)
    flats = {}
    for p in mst_dir.glob("*skyflat*"):
        f = re.search(r"[-_]([A-Za-z0-9]+)\.fits", p.name)
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

        # 오버스캔·트림 시점은 기기 세대마다 다르다. 실측한 세 경우:
        #   QHY600 central30x30 : BIASSEC='N/A', raw·마스터 모두 트림 전 크기
        #   QHY600 full_frame   : 상단 30행 BIASSEC, 마스터가 전부 **트림 후** 크기
        #   SBIG STL-6303       : BIASSEC='UNKNOWN', **bias 만 트림 전**,
        #                         dark·flat 은 트림 후 → bias 를 빼고 잘라야 한다
        # 그래서 순서를 고정하지 않고, 마스터를 적용하기 직전에 크기를 맞춘다.
        bs, ts = _sec(H.get("BIASSEC", "")), _sec(H.get("TRIMSEC", ""))
        over = 0.0
        if bs is not None:
            x1, x2, y1, y2 = bs
            over = float(np.median(raw[y1:y2, x1:x2]))
            raw = raw - np.float32(over)

        def _align(arr: np.ndarray, target: tuple[int, int], what: str) -> np.ndarray:
            if arr.shape == target:
                return arr
            if ts is not None:
                x1, x2, y1, y2 = ts
                cut = arr[y1:y2, x1:x2]
                if cut.shape == target:
                    return cut
            raise SystemExit(
                f"{p.name}: {what} 적용 전 크기 {arr.shape} → 목표 {target} 로 못 맞춤"
            )

        # 마스터 bias 의 단위도 세대마다 다르다:
        #   QHY600(2024) : 전자 단위, median≈0. BIASLVL 이 별도 스칼라 레벨
        #                  → raw*gain - BIASLVL - bias
        #   SBIG(2015)   : **ADU 단위**, median≈1048. BIASLVL = median*gain
        #                  → (raw - bias)*gain   (BIASLVL 을 또 빼면 이중 차감)
        # BIASLVL 이 median(bias)*gain 과 일치하는지로 판별한다.
        mb_med = float(np.median(mbias))
        bias_in_adu = (
            abs(biaslvl) > 1.0
            and abs(mb_med * gain - biaslvl) < 0.05 * abs(biaslvl)
        )
        if bias_in_adu:
            x = (_align(raw, mbias.shape, "bias") - mbias) * np.float32(gain)
        else:
            x = raw * np.float32(gain) - np.float32(biaslvl)
            x = _align(x, mbias.shape, "bias") - mbias
        x = _align(x, mdark.shape, "dark") - mdark * np.float32(lexp / dexp)
        x = _align(x, mf.shape, "flat")
        apex = x / np.where(np.abs(mf) < 1e-6, np.nan, mf)

        # 첫 프레임 중 e91 이 있는 것으로 보정 산술 대조
        ref = ref_dir / p.name.replace("-e00", "-e91")
        if ref.exists():
            e91, _ = _sci(ref)
            diff = apex - e91
            fin = np.isfinite(diff)
            med = float(np.median(diff[fin]))
            mad = float(1.4826 * np.median(np.abs(diff[fin] - med)))
            sky = float(np.median(e91[fin]))
            # 게이트는 하늘 대비 상대값으로 본다. 절대 e- 로 잡으면 하늘이
            # 밝은 세트(M45 sky 399 e-)와 어두운 세트(U 300 s, sky 73 e-)를
            # 같은 잣대로 재게 되고, 읽기잡음이 큰 옛 CCD 가 부당하게 걸린다.
            rel_med = abs(med) / max(abs(sky), 1e-6)
            rel_mad = mad / max(abs(sky), 1e-6)
            print(f"  [대조] {p.name[-18:]} vs BANZAI: Δmed={med:+.4f} e- ({rel_med*100:.2f}%) "
                  f"robustσ={mad:.4f} ({rel_mad*100:.2f}%) (sky {sky:.1f} e-)")
            if rel_med > 0.05 or rel_mad > 0.15:
                raise SystemExit("BANZAI 대조 실패 — 보정 산술을 확인할 것")
            checked = True

        n_nan = int(np.sum(~np.isfinite(apex)))
        apex[~np.isfinite(apex)] = np.nanmedian(apex)

        # 우주선·핫픽셀 제거 — APEX Step 0 과 같은 함수(astroscrappy = L.A.Cosmic,
        # objlim 으로 실제 천체를 보호)를 쓴다. BANZAI e91 에는 이 처리가 들어
        # 있지만 우리는 raw 에서 직접 보정하므로 여기서 해야 한다.
        # 안 하면 CMOS 프레임의 점광원 잡음이 EPSF 참조별로 뽑혀 PSF 가 무너진다
        # (2026-07-29 실측: M67/QHY600 은 FWHM<=1.5px 검출이 40% 였고 EPSF 가
        # 실제 별보다 2.75배 좁게 만들어져 PSF 플럭스가 구경의 32% 로 떨어졌다.
        # Sinistro CCD 는 0% 라 무사했다).
        n_cr = 0
        if clean_cr:
            hmask = hot_pixel_mask(mdark, 6.0) if mdark is not None else None
            apex, _cmask, n_cr = clean_frame(
                apex,
                gain=1.0,                       # 이 배열은 이미 전자 단위
                readnoise=float(H.get("RDNOISE", 6.5)),
                satlevel=float(rh.get("SATURATE", 65535.0)),
                sigclip=4.5, objlim=5.0, hot_mask=hmask,
            )

        ov_txt = f" over={over:.1f}ADU" if bs is not None else ""
        cr_txt = f" cr={n_cr}" if clean_cr else " cr=off"

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
        print(f"  {name}  median={np.nanmedian(apex):7.2f} e-  nan_filled={n_nan}{ov_txt}{cr_txt}")

    if not checked:
        raise SystemExit("e91 대조를 한 번도 못 했다 — banzai_ref 파일명 확인")
    print(f"\n{key}: {sum(counter.values())}장 → {out_dir}  ({counter})")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    clean = "--no-cr" not in sys.argv
    for k in (args or ["m67_lco"]):
        main(k, clean_cr=clean)
