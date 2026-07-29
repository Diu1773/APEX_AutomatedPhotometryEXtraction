# -*- coding: utf-8 -*-
"""LCO 교차기기 PSF 테스트용 parameters.toml 생성.

검증된 reprocess M13 설정을 탬플릿으로 삼아 기기별 값만 치환한다
(reprocess_batch.gen_config 와 같은 방식). 프레임 단위가 전자(e-)이므로
gain 1.0, RDNOISE 는 BANZAI 헤더값.

실행: .venv-deploy\\Scripts\\python -X utf8 validation/psf_crossinstrument/gen_configs.py
"""
from __future__ import annotations

import re
from pathlib import Path

TEMPLATE = Path(r"E:\APEX_validation\reprocess\M13\parameters.toml")
BASE = Path(r"E:\APEX_validation\psf_crossinstrument")

# (하위폴더, 치환 dict)
CAMERAS = {
    # LCO 1m elp / Sinistro fa16 (Fairchild CCD, 4-amp) — NGC 5985 field, rp 120 s
    # 0.389"/px, 15 um px → focal = 206.265*15/0.389 = 7953 mm
    "sinistro": {
        "target_name": "NGC5985_field",
        "ra_deg": 234.90503,      # 15:39:37.207
        "dec_deg": 59.33209,      # +59:19:55.54
        "telescope_focal_mm": 7953.0,
        "camera_pixel_um": 15.0,
        "binning": 1,
        "gain_e_per_adu": 1.0,    # 프레임이 이미 전자 단위
        "rdnoise_e": 8.4,
        "saturation_adu": 126000.0,
        "datamax_adu": 120000.0,
        "datamin_adu": -100.0,
        "pixel_scale_arcsec": 0.389,
        "guess_arcsec": 2.0,
        "filter_key": "rp",
        "site_lat": 30.679833, "site_lon": -104.015173, "site_alt": 2030.0, "site_tz": -6.0,
    },
    # LCO 0.4m coj / sq36 = QHY600 CMOS (단일앰프) — Proxima Cen field, V 20 s
    # 0.74"/px, 3.76 um px → focal = 206.265*3.76/0.74 = 1048 mm
    "qhy600": {
        "target_name": "ProximaCen_field",
        "ra_deg": 217.36826,      # 14:29:28.382
        "dec_deg": -62.67383,     # -62:40:25.80
        "telescope_focal_mm": 1048.0,
        "camera_pixel_um": 3.76,
        "binning": 1,
        "gain_e_per_adu": 1.0,
        "rdnoise_e": 3.08,
        "saturation_adu": 47400.0,
        "datamax_adu": 45000.0,
        "datamin_adu": -100.0,
        "pixel_scale_arcsec": 0.74,
        "guess_arcsec": 2.5,
        "filter_key": "V",
        "site_lat": -31.2728196, "site_lon": 149.0708466, "site_alt": 1130.0, "site_tz": 10.0,
        # 남천 -62°: 로컬 astrometry.net 인덱스(2mass-04 시리즈 48청크 중 31개)가
        # 남천을 안 덮어 blind 풀이 실패 → 전천 D50 을 가진 ASTAP 으로 지정.
        "wcs_engine": "astap",
    },
}


def _sub(txt: str, key: str, value, section_hint: str | None = None) -> str:
    """^key = ... 행을 치환한다. 문자열 값은 따옴표로 감싼다."""
    if isinstance(value, str):
        rep = f'{key} = "{value}"'
    elif isinstance(value, bool):
        rep = f"{key} = {str(value).lower()}"
    else:
        rep = f"{key} = {value}"
    pat = re.compile(rf"(?m)^{re.escape(key)}\s*=.*$")
    n = len(pat.findall(txt))
    if n == 0:
        raise KeyError(f"template lacks key: {key}")
    return pat.sub(lambda m: rep, txt, count=0 if section_hint == "ALL" else 1), n


def main() -> None:
    base_txt = TEMPLATE.read_text(encoding="utf-8-sig")
    esc = lambda p: str(p).replace("\\", "\\\\")
    for cam, C in CAMERAS.items():
        txt = base_txt
        report = []

        def one(key, value, allmatch=False):
            nonlocal txt
            txt, n = _sub(txt, key, value, "ALL" if allmatch else None)
            report.append(f"{key}={value!r}({n})")

        one("data_dir", esc(BASE / cam / "sci"))
        one("result_dir", esc(BASE / cam / "result"))
        one("name", C["target_name"])
        one("ra_deg", C["ra_deg"])
        one("dec_deg", C["dec_deg"])
        one("telescope_focal_mm", C["telescope_focal_mm"])
        one("camera_pixel_um", C["camera_pixel_um"])
        one("binning", C["binning"])
        one("gain_e_per_adu", C["gain_e_per_adu"])
        one("rdnoise_e", C["rdnoise_e"])
        one("saturation_adu", C["saturation_adu"])
        one("datamax_adu", C["datamax_adu"])
        one("datamin_adu", C["datamin_adu"])
        one("pixel_scale_arcsec", C["pixel_scale_arcsec"])
        one("guess_arcsec", C["guess_arcsec"])
        one("global_ref_filter", C["filter_key"])
        one("filter_keep", C["filter_key"])
        one("lat_deg", C["site_lat"])
        one("lon_deg", C["site_lon"])
        one("alt_m", C["site_alt"])
        one("tz_offset_hours", C["site_tz"])
        # WCS 솔버 스케일 범위 — 탬플릿(0.393"/px 기기)값이 새면 0.74"/px 는 범위 밖
        one("astnet_local_scale_low", round(C["pixel_scale_arcsec"] * 0.85, 3))
        one("astnet_local_scale_high", round(C["pixel_scale_arcsec"] * 1.15, 3))
        # 솔버 엔진 명시 (탬플릿에 키가 없으므로 [wcs] 절 머리에 삽입)
        if C.get("wcs_engine"):
            txt = txt.replace("[wcs]\n", f'[wcs]\nengine = "{C["wcs_engine"]}"\n', 1)
            report.append(f"engine={C['wcs_engine']!r}(insert)")
        # IRAF 도구 절의 epadu/readnoise 도 전자 단위로 (검증에 안 쓰이지만 일관성)
        txt, _ = _sub(txt, "epadu", 1.0, "ALL")
        txt, _ = _sub(txt, "readnoise", C["rdnoise_e"], "ALL")

        out = BASE / cam / "parameters.toml"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(txt, encoding="utf-8")
        print(f"[{cam}] wrote {out}")
        print("   " + " | ".join(report[:8]))
        print("   " + " | ".join(report[8:]))


if __name__ == "__main__":
    main()
