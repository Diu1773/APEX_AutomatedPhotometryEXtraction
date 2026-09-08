"""LCO 아카이브 워크스페이스의 apex_config.json 을 헤더에서 읽어 만든다.

**값을 손으로 적지 않는다.** 기기 상수는 프레임 헤더에 있고, 그것이 정본이다.
이 스크립트는 과학 프레임 한 장을 열어 필요한 값을 꺼내고, 나머지는
`parameters.example.json` 의 기본값을 그대로 쓴다.

**밴드 폴더 이름을 박지 않는다.** 처음에는 `gp·rp·ip` 를 박아 두었는데, 그러면
기기가 바뀔 때마다 이 파일을 고쳐야 한다 — 하드코딩을 고치면서 하드코딩을 늘리는
일이다(2026-09-08 교정 C-201). `inputs/science/` 밑에 있는 폴더를 그대로 훑는다.
kb26 은 `B·V·rp·ip·zs`, MuSCAT3 는 `gp·rp·ip` 인데 코드는 같다.

**헤더에서 가져오는 것**

    GAIN·RDNOISE   프레임마다 다르므로 값을 박지 않고 noise_use_fits_header 를
                   켠다. MuSCAT3 는 카메라마다 1.90·1.88·1.80 으로 다르다.
                   PTC 로 재서 헤더가 1 % 안에서 맞는 것을 확인했다
                   (validation/external_muscat3/PTC_MuSCAT3.md).
    SATURATE       포화 한계 (ADU)
    PIXSCALE       화소 크기 (초각). CD 행렬과도 대조한다.
    CRVAL1/2       그 밤의 겨눔 좌표 — 대상 좌표의 출발점
    BIASSEC        오버스캔 자리. 있으면 오버스캔을 켜고 폭을 여기서 정한다.

실행:
    python -X utf8 scripts/make_lco_config.py --job <작업폴더> [--dry-run]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).absolute().parents[1]))

from apex.utils.io_utils import read_fits_header  # noqa: E402

REPO = Path(__file__).absolute().parents[1]


def _first_science_frame(job: Path) -> Path:
    """`inputs/science/` 밑을 훑어 첫 과학 프레임을 찾는다.

    밴드 폴더 이름을 목록으로 갖고 있지 않다 — 기기마다 다르기 때문이다.
    폴더가 없이 파일이 바로 놓인 경우도 받는다.
    """
    root = job / "inputs" / "science"
    if not root.is_dir():
        raise SystemExit(f"과학 프레임 폴더가 없다: {root}")
    hits = sorted(root.rglob("*.fits*"))
    if not hits:
        raise SystemExit(f"과학 프레임을 못 찾았다: {root}")
    bands = sorted({p.parent.name for p in hits if p.parent != root})
    if bands:
        print(f"밴드 폴더: {' · '.join(bands)}")
    return hits[0]


def _parse_biassec(value) -> tuple[str, int] | None:
    """``[2049:2080,1:2048]`` -> ("right", 32). 세로 오버스캔도 알아본다."""
    try:
        body = str(value).strip().strip("[]")
        xs, ys = body.split(",")
        x0, x1 = (int(v) for v in xs.split(":"))
        y0, y1 = (int(v) for v in ys.split(":"))
    except Exception:
        return None
    if (x1 - x0) < (y1 - y0):
        return ("right" if x0 > 1 else "left", x1 - x0 + 1)
    return ("top" if y0 > 1 else "bottom", y1 - y0 + 1)


def build(job: Path) -> dict:
    frame = _first_science_frame(job)
    h = read_fits_header(frame)
    print(f"헤더를 읽은 프레임: {frame.name}")

    pixscale = float(h.get("PIXSCALE") or 0.0)
    cd11 = float(h.get("CD1_1") or 0.0)
    cd_scale = abs(cd11) * 3600.0 if cd11 else 0.0
    if pixscale and cd_scale and abs(pixscale - cd_scale) / pixscale > 0.02:
        print(f"  주의: PIXSCALE {pixscale:.4f} 와 CD 행렬 {cd_scale:.4f} 가 2 % 넘게 다르다")
    scale = pixscale or cd_scale
    if not scale:
        raise SystemExit("화소 크기를 헤더에서 못 읽었다 (PIXSCALE·CD1_1 둘 다 없음)")

    data = json.loads((REPO / "parameters.example.json").read_text(encoding="utf-8"))

    data.setdefault("io", {})
    data["io"]["data_dir"] = str(job / "inputs")
    data["io"]["result_dir"] = str(job / "results")
    # 파일 이름이 ogg2m001-ep04-20210317-0059-e00 꼴이다. 밤은 하이픈 사이의 8 자리.
    data["io"]["night_parse_mode"] = "regex"
    data["io"]["night_parse_regex"] = r".*-(20\d{6})-"

    data.setdefault("target", {})
    data["target"]["name"] = str(h.get("OBJECT") or "M67").strip()
    data["target"]["ra_deg"] = float(h.get("CRVAL1"))
    data["target"]["dec_deg"] = float(h.get("CRVAL2"))

    inst = data.setdefault("instrument", {})
    # gain·rdnoise 는 프레임마다 헤더에서 읽는다 (카메라마다 다르다).
    inst["noise_use_fits_header"] = True
    inst["gain_e_per_adu"] = float(h.get("GAIN") or 0) or inst.get("gain_e_per_adu")
    inst["rdnoise_e"] = float(h.get("RDNOISE") or 0) or inst.get("rdnoise_e")
    sat = float(h.get("SATURATE") or 0)
    if sat:
        inst["saturation_adu"] = sat
        # MAXLIN 이 따로 없으면 포화를 선형 한계로도 쓴다. 헤더가 둘을 같은
        # 값으로 적고 있으므로 여기서 임의로 낮추지 않는다.
        inst["datamax_adu"] = float(h.get("MAXLIN") or sat)
    inst["binning"] = int(str(h.get("CCDSUM") or "1 1").split()[0])

    # 화소 크기는 초점거리·화소크기 조합 대신 잰 값을 직접 준다. 헤더의
    # CCDXPIXE 는 단위가 [m] 이라고 적혀 있으나 실제로는 화소 개수라 못 쓴다.
    data.setdefault("match", {})["pixel_scale_arcsec"] = round(scale, 4)

    cal = data.setdefault("calibration", {})
    cal["enabled"] = True
    over = _parse_biassec(h.get("BIASSEC"))
    if over:
        edge, width = over
        cal["overscan"] = {"enable": True, "edge": edge,
                           "width": width, "trim": True}
        print(f"  오버스캔: {edge} 쪽 {width} 열 (BIASSEC {h.get('BIASSEC')})")
    else:
        cal.setdefault("overscan", {})["enable"] = False

    print(f"  대상    : {data['target']['name']} "
          f"({data['target']['ra_deg']:.4f}, {data['target']['dec_deg']:.4f})")
    print(f"  화소    : {scale:.4f} 초각/px")
    print(f"  GAIN    : {inst['gain_e_per_adu']} e-/ADU  (프레임마다 헤더에서 다시 읽음)")
    print(f"  RDNOISE : {inst['rdnoise_e']} e-")
    print(f"  포화    : {inst.get('saturation_adu')} ADU")
    return data


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="LCO 워크스페이스 설정을 헤더에서 만든다")
    ap.add_argument("--job", required=True, help="작업 폴더 (inputs/ 를 담고 있는)")
    ap.add_argument("--dry-run", action="store_true", help="쓰지 않고 보여만 준다")
    a = ap.parse_args(argv)

    job = Path(a.job)
    data = build(job)
    out = job / "apex_config.json"
    if a.dry_run:
        print(f"\n(--dry-run) {out} 에 쓰지 않았다")
        return 0
    out.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n",
                   encoding="utf-8")
    print(f"\n썼다: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
