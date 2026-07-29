# -*- coding: utf-8 -*-
"""LCO 교차기기 PSF 테스트용 sci 프레임 준비.

두 공개 LCO 프레임(archive.lco.global)을 APEX식 보정 산술로 처리해
파이프라인 입력용 FITS 로 저장한다. 보정 산술은 Fig13 에서 BANZAI 와
비트 수준 대조가 끝난 것과 동일하다 (Sinistro robust σ 0.32%, QHY +0.06 e-).

  Sinistro (1m elp, fa16, Fairchild CCD 4-amp, 0.389"/px) — NGC 5985, rp 120 s
      4앰프 조립(오버스캔→트림→앰프별 gain→DETSEC 배치) 후 (raw-b-d*r)/f
  QHY600 (0.4m coj, sq36 CMOS 단일앰프, 0.74"/px) — Proxima Cen 필드, V 20 s
      (raw*GAIN - BIASLVL - mbias - mdark*ratio)/mflat   (BANZAI 레시피)

출력 단위는 전자(e-)이므로 APEX 설정은 gain=1.0 을 쓴다.

실행: .venv-deploy\\Scripts\\python -X utf8 validation/psf_crossinstrument/prepare_lco_sci.py
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
from astropy.io import fits

EXT = Path(r"E:\APEX_validation\external")
OUT = Path(r"E:\APEX_validation\psf_crossinstrument")

# 원본 헤더에서 sci 헤더로 넘길 키 (측광·야간분류·좌표 시드에 필요한 것만)
KEEP = [
    "DATE-OBS", "EXPTIME", "FILTER", "OBJECT", "RA", "DEC", "CAT-RA", "CAT-DEC",
    "AIRMASS", "INSTRUME", "TELESCOP", "SITEID", "LATITUDE", "LONGITUD", "HEIGHT",
    "PIXSCALE", "RDNOISE", "MJD-OBS", "UTSTART", "DAY-OBS",
]


def _sec(s: str) -> tuple[int, int, int, int]:
    m = re.findall(r"-?\d+", s)
    return tuple(int(x) for x in m)  # type: ignore[return-value]


def _sci_of(fz: Path, shape=None):
    with fits.open(fz) as h:
        for hd in h:
            if hd.data is not None and hd.data.ndim == 2 and (
                shape is None or hd.data.shape == shape
            ):
                return hd.data.astype(np.float64), hd.header
    raise RuntimeError(f"no 2-D HDU in {fz}")


def _write_sci(out_path: Path, data: np.ndarray, src_header, extra: dict):
    hdr = fits.Header()
    for k in KEEP:
        if k in src_header:
            hdr[k] = src_header[k]
    hdr["BUNIT"] = ("electron", "calibrated frame in electrons")
    hdr["EGAIN"] = (1.0, "data already in electrons")
    for k, (v, c) in extra.items():
        hdr[k] = (v, c)
    hdr["CALSRC"] = ("APEX-arith vs BANZAI (fig13)", "calibration provenance")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=data.astype(np.float32), header=hdr).writeto(
        out_path, overwrite=True
    )
    print(f"  wrote {out_path}  {data.shape}  median={np.nanmedian(data):.2f} e-")


def prepare_sinistro() -> None:
    L = EXT / "LCO_sinistro"
    print("[SINISTRO] assemble 4-amp raw -> electrons")
    with fits.open(L / "raw_e00.fits.fz") as h:
        mosaic = np.zeros((4096, 4096), np.float64)
        for i in range(1, 5):
            hd = h[i].header
            d = h[i].data.astype(np.float64)
            x1, x2, y1, y2 = _sec(hd["DATASEC"])
            bx1, bx2, by1, by2 = _sec(hd["BIASSEC"])
            gain = float(hd["GAIN"])
            over = np.median(d[by1 - 1:by2, bx1 - 1:bx2])
            sci = (d[y1 - 1:y2, x1 - 1:x2] - over) * gain
            dx1, dx2, dy1, dy2 = _sec(hd["DETSEC"])
            if dx1 > dx2:
                sci = np.fliplr(sci)
                dx1, dx2 = dx2, dx1
            if dy1 > dy2:
                sci = np.flipud(sci)
                dy1, dy2 = dy2, dy1
            mosaic[dy1 - 1:dy2, dx1 - 1:dx2] = sci
        raw_hdr = h[1].header.copy()
        for k in KEEP:  # 앰프 헤더에 없는 키는 프라이머리에서
            if k not in raw_hdr and k in h[0].header:
                raw_hdr[k] = h[0].header[k]
        lexp = float(h[0].header.get("EXPTIME") or raw_hdr.get("EXPTIME") or 120.0)

    mb, _ = _sci_of(L / "master_bias.fits.fz", (4096, 4096))
    md, mdh = _sci_of(L / "master_dark.fits.fz", (4096, 4096))
    mf, _ = _sci_of(L / "master_flat.fits.fz", (4096, 4096))
    dexp = float(mdh.get("EXPTIME", 1.0))
    ratio = lexp / dexp if dexp else 1.0
    apex = (mosaic - mb - md * ratio) / np.where(np.abs(mf) < 1e-6, np.nan, mf)

    # 자기검증: BANZAI e91 대비 (fig13 수치 재현 확인)
    e91, _ = _sci_of(L / "banzai_e91.fits.fz", (4096, 4096))
    d = apex - e91
    fin = np.isfinite(d)
    med = np.median(d[fin])
    mad = 1.4826 * np.median(np.abs(d[fin] - med))
    print(f"  vs BANZAI: Δmedian={med:+.3f} e-  robustσ={mad:.3f} e-  (fig13: +1.95 / 0.74)")
    assert abs(med) < 4.0 and mad < 2.0, "BANZAI 대조 실패 — 보정 산술 확인"

    # NaN(플랫 0 나눗셈)은 파이프라인에 들어가기 전에 중앙값으로 대치
    n_nan = int(np.sum(~np.isfinite(apex)))
    apex[~np.isfinite(apex)] = np.nanmedian(apex)
    print(f"  NaN filled: {n_nan}")

    _write_sci(
        OUT / "sinistro" / "sci" / "pp_ngc5985-0001-rp.fit",
        apex,
        raw_hdr,
        {
            "SATURATE": (126000.0, "e-, from BANZAI e91 header"),
        },
    )


def prepare_qhy600() -> None:
    Q = EXT / "LCO_qhy600_0m4"
    print("[QHY600] single-CCD reduce -> electrons")
    raw, raw_hdr = _sci_of(Q / "raw_e00.fits.fz", (2400, 2400))
    mb, _ = _sci_of(Q / "master_bias.fits.fz", (2400, 2400))
    md, mdh = _sci_of(Q / "master_dark.fits.fz", (2400, 2400))
    mf, _ = _sci_of(Q / "master_flat.fits.fz", (2400, 2400))
    e91, e91h = _sci_of(Q / "banzai_e91.fits.fz", (2400, 2400))

    gain = float(raw_hdr.get("GAIN", 0.79))
    biaslvl = float(e91h["BIASLVL"])
    lexp = float(raw_hdr.get("EXPTIME", 20.0))
    dexp = float(mdh.get("EXPTIME", 300.0))
    ratio = lexp / dexp if dexp else 1.0
    apex = (raw * gain - biaslvl - mb - md * ratio) / np.where(
        np.abs(mf) < 1e-6, np.nan, mf
    )

    d = apex - e91
    fin = np.isfinite(d)
    med = np.median(d[fin])
    mad = 1.4826 * np.median(np.abs(d[fin] - med))
    print(f"  vs BANZAI: Δmedian={med:+.4f} e-  robustσ={mad:.4f} e-  (fig13: +0.063 / 0.080)")
    assert abs(med) < 0.5 and mad < 0.5, "BANZAI 대조 실패 — 보정 산술 확인"

    n_nan = int(np.sum(~np.isfinite(apex)))
    apex[~np.isfinite(apex)] = np.nanmedian(apex)
    print(f"  NaN filled: {n_nan}")

    _write_sci(
        OUT / "qhy600" / "sci" / "pp_proximafield-0001-V.fit",
        apex,
        raw_hdr,
        {
            "SATURATE": (47400.0, "e-, from BANZAI e91 header"),
        },
    )


if __name__ == "__main__":
    prepare_sinistro()
    prepare_qhy600()
    print("done.")
