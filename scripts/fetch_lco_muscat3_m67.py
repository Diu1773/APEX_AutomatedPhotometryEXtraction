"""LCO MuSCAT3 의 M67 을 원본·보정·BANZAI 처리본까지 받는다.

**왜 이것인가.** APEX 는 M67 을 Moravian C3-61000 으로 이미 처리해 두었다
(`validation/paper/data_realframe_M67g_broad` · `M67r_mid` · `M67i`). MuSCAT3 은
2 m 망원경(ogg 2m0a)에 붙어 **gp·rp·ip·zs 네 밴드를 동시에** 찍고, 그중 셋이
APEX 의 g·r·i 와 겹친다. **같은 성단 · 같은 필터 · 다른 기기**라서 표준성야보다
직접적인 대조가 된다(STATUS.md 의 목표 둘).

**기기가 실제로 얼마나 다른가** (헤더에서 직접 읽은 값).

    항목        MuSCAT3 (ep04)      Moravian C3-61000
    GAIN        1.9 e-/ADU          0.689
    RDNOISE     12.0 e-             2.1
    화소        0.266 초각/px       0.393 (M13)
    포화        64,000 ADU          -
    오버스캔    있음 (BIASSEC)      없음

**어느 밤을 쓰나.** M67 관측은 네 밤에 흩어져 있고 크기가 아주 다르다.

    밤          밴드당 과학    BIAS   DARK   SKYFLAT
    2021-03-14        1         -      -       -
    2021-03-17       60        64     20     없음
    2021-03-19        6        64     20      12
    2021-03-20       10        32     10      12

**03-17 을 쓴다.** 프레임이 열 배 많고 BANZAI 처리본도 같은 밤에 60 장씩 있어
짝이 맞기 때문이다. 다만 **그 밤에는 skyflat 이 없어서 다음 날 아침(03-18) 것을
쓴다.** 관측소에서 흔히 하는 일이고 LCO 자신의 BANZAI 도 가장 가까운 좋은 보정을
쓰지만, **논문에는 「플랫은 이웃 밤 것」이라고 적어야 한다.** (2026-09-07 결정)

**BANZAI 처리본도 받는다.** level 91 이 LCO 자신의 파이프라인이 낸 결과다. 같은
밤 같은 밴드로 60 장씩 있으므로 **APEX 의 Step 0 을 관측소의 파이프라인과 화소
단위로 견줄 수 있다** — 원고의 Fig 13 이 쓰는 바로 그 방식이다.

**크기.** 한 장 3.4 MB 로 재서 확인했다.

    과학 gp·rp·ip (네 밤 전부)   228 장   0.78 GB
    BIAS  03-17 · 세 카메라      192 장   0.65 GB
    DARK  03-17 · 세 카메라       60 장   0.20 GB
    SKYFLAT 03-18 · 세 카메라     36 장   0.12 GB
    BANZAI 03-17 · 세 밴드       180 장   (level 91 은 더 크다)

인증이 필요 없고 표준 라이브러리만 쓴다. 이미 받은 파일은 건너뛰므로 끊겨도 다시
돌리면 이어진다.

실행:
    python -X utf8 scripts/fetch_lco_muscat3_m67.py --out D:/LCO_M67_MuSCAT3
    python -X utf8 scripts/fetch_lco_muscat3_m67.py --out ... --part calib
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).absolute().parent))

from fetch_lco import download, query  # noqa: E402

# 과학 프레임은 네 밤을 다 받아 둔다(작다). 실제로 쓸 밤은 03-17 이다.
SCIENCE_WINDOW = ("2021-03-01", "2021-04-01")
SCIENCE_NIGHT = "2021-03-17"

# 보정: bias·dark 는 같은 밤, skyflat 은 그 밤에 없어 다음 날 아침 것을 쓴다.
CALIB_PLAN = (
    ("BIAS", "2021-03-17", "2021-03-18"),
    ("DARK", "2021-03-17", "2021-03-18"),
    ("SKYFLAT", "2021-03-18", "2021-03-19"),
)
CALIB_NIGHT = {"BIAS": "2021-03-17", "DARK": "2021-03-17", "SKYFLAT": "2021-03-18"}

# MuSCAT3 은 카메라마다 밴드가 고정이다. APEX 의 g·r·i 와 겹치는 셋이 앞의 셋.
CAMERAS = (("ep04", "gp"), ("ep02", "rp"), ("ep03", "ip"), ("ep05", "zs"))
DEFAULT_BANDS = ("gp", "rp", "ip")


def _night_of(row: dict) -> str:
    return str(row.get("observation_day") or "")


def _run(label: str, out: Path, rows: list) -> tuple[int, int]:
    if not rows:
        print(f"  [{label}] 맞는 프레임이 없다.", flush=True)
        return 0, 0
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    n_bytes = download(rows, out)
    dt = time.perf_counter() - t0
    speed = f" · {n_bytes/1e6/dt:.1f} MB/s" if dt > 0 and n_bytes else ""
    print(f"  [{label}] {len(rows)} 장 · {n_bytes/1e9:.2f} GB · {dt:.0f} 초{speed}",
          flush=True)
    return len(rows), n_bytes


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="MuSCAT3 의 M67 을 받는다 (03-17 밤)")
    ap.add_argument("--out", required=True, help="받을 폴더")
    ap.add_argument("--part", nargs="+", default=["science", "calib", "banzai"],
                    choices=["science", "calib", "banzai"],
                    help="받을 부분 (기본: 셋 다)")
    ap.add_argument("--bands", nargs="+", default=list(DEFAULT_BANDS),
                    help=f"밴드 (기본: {' '.join(DEFAULT_BANDS)})")
    a = ap.parse_args(argv)

    out = Path(a.out)
    want = set(a.bands)
    cams = [(c, b) for c, b in CAMERAS if b in want]
    parts = set(a.part)

    total_n = total_b = 0
    t0 = time.perf_counter()

    if "science" in parts:
        print(f"== M67 원본 (level 0, {SCIENCE_WINDOW[0]} ~ {SCIENCE_WINDOW[1]}) ==",
              flush=True)
        for cam, band in cams:
            rows = query("M67", 0, 300, cam, None,
                         start=SCIENCE_WINDOW[0], end=SCIENCE_WINDOW[1])
            n, b = _run(f"{cam} {band}", out / "science" / band, rows)
            total_n += n
            total_b += b

    if "calib" in parts:
        print(f"\n== 보정 프레임 (bias·dark {CALIB_NIGHT['BIAS']} · "
              f"skyflat {CALIB_NIGHT['SKYFLAT']}) ==", flush=True)
        for ctype, start, end in CALIB_PLAN:
            night = CALIB_NIGHT[ctype]
            for cam, band in cams:
                rows = query(None, 0, 500, cam, None, config_type=ctype,
                             start=start, end=end)
                # 질의 창이 하루를 넘겨 잡히므로 밤으로 한 번 더 거른다.
                rows = [r for r in rows if _night_of(r) == night]
                n, b = _run(f"{ctype} {cam} {band}",
                            out / "calib" / ctype.lower() / band, rows)
                total_n += n
                total_b += b

    if "banzai" in parts:
        print(f"\n== BANZAI 처리본 (level 91, {SCIENCE_NIGHT} 밤) ==", flush=True)
        for cam, band in cams:
            rows = query("M67", 91, 300, cam, None,
                         start=SCIENCE_WINDOW[0], end=SCIENCE_WINDOW[1])
            rows = [r for r in rows if _night_of(r) == SCIENCE_NIGHT]
            n, b = _run(f"{cam} {band}", out / "banzai" / band, rows)
            total_n += n
            total_b += b

    dt = time.perf_counter() - t0
    print(f"\n합계 {total_n} 장 · {total_b/1e9:.2f} GB · {dt/60:.1f} 분", flush=True)
    print(f"받은 곳: {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
