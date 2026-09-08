"""LCO 0.4 m (kb26) 의 M67 을 원본·보정·BANZAI 처리본까지 받는다 — 세 번째 기기.

**왜 또 다른 기기인가.** APEX 는 지금 M67 을 두 기기로 처리해 두었다 — 사장님의
Moravian C3-61000(CMOS)과 LCO 2 m 의 MuSCAT3 다. 둘 다 잘 맞았지만
(g −0.015 · r −0.019 · i −0.015), **둘 다 중형 이상이고 SDSS 필터만 겹쳤다.**
그리고 MuSCAT3 의 rp 에서 1 차 색항이 적합한 별에서만 맞는 고장이 나왔는데
(`validation/COLOR_TERM_STABILITY.md`), 그것이 그 기기만의 일인지 알 수 없다.

**kb26 이 채우는 빈칸은 셋이다.**

    검출기      SBIG STL-6303E (Kodak KAF-6303E) — **CCD**
                앞의 둘은 Moravian 이 CMOS(IMX455), MuSCAT3 이 별개의 CCD 다.
    구경        0.4 m — **APEX 가 실제로 겨냥하는 급**이다. 앞의 둘은 2 m 와
                사장님 망원경이었다.
    필터        Johnson **B·V** 가 있다. 사장님의 M13·M3·NGC 6811 이 B/V/R 이므로,
                지금까지 못 해 본 **Johnson 계열의 기기 사이 대조**가 열린다.

**어느 밤을 쓰나 — 2019-05-03 (dayobs).** LCO 0.4 m 망원경들이 M67 을 찍은 밤이
열일곱인데 대부분 밴드당 한두 장이다. 이 밤만 다르다.

    밴드      장수   노출
    B          6     40 s
    V          6     40 s
    rp         6     40 s
    ip         6     40 s
    zs         6     60 s
    up         3    450 s

**밤은 반드시 `dayobs` 로 고른다.** 날짜 구간으로 고르면 관측소의 시차 때문에
그 밤이 통째로 빠지거나 이웃 밤이 섞인다(`fetch_lco.query` 의 설명 참조).
실제로 이 자료도 UTC 로는 05-03 과 05-04 에 걸쳐 있다.

**다섯 밴드 모두 색지수가 붙는다** (`FILTER_COLOR_PREF`).

    B → B−V     V → B−V     r → r−i     i → r−i     z → r−z

gp 가 없어서 r 은 첫 후보(g−r)를 못 쓰고 둘째 후보 r−i 로 간다. **이것도 볼거리다** —
같은 밴드를 다른 색축으로 맞추면 색항이 어떻게 달라지는지 여기서 보인다.

**보정은 같은 밤 것이 다 있다** — MuSCAT3 때는 skyflat 이 없어 이웃 밤 것을 썼는데
여기는 그럴 필요가 없다.

    BIAS      34 장   0 s
    DARK      10 장   300 s
    SKYFLAT   B·V·rp·ip·zs 각 3~5 장

**다크 노출이 과학 프레임보다 훨씬 길다**(300 s 대 40 s). Step 0 이 노출 시간으로
비례 축소해 빼므로 문제는 없지만, 읽기잡음이 그만큼 확대되어 들어온다.
**결과를 볼 때 이 점을 기억해야 한다.**

**BANZAI 처리본(level 91)도 받는다.** 관측소 자신의 파이프라인이 낸 결과이므로
APEX 의 Step 0 을 화소 단위로 견줄 수 있다 — MuSCAT3 에서 상관 0.9998 을 낸 그
방식이다.

**크기.** 한 장 6~8 MB 로 잡아 전부 1.1 GB 안팎이다.

인증이 필요 없고 표준 라이브러리만 쓴다. 이미 받은 파일은 건너뛰므로 끊겨도 다시
돌리면 이어진다.

실행:
    python -X utf8 scripts/fetch_lco_kb26_m67.py --out D:/LCO_M67_kb26
    python -X utf8 scripts/fetch_lco_kb26_m67.py --out ... --part calib
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).absolute().parent))

from fetch_lco import download, query  # noqa: E402

INSTRUMENT = "kb26"
NIGHT = "2019-05-03"          # 과학·보정이 모두 이 밤에 있다
CALIB_TYPES = ("BIAS", "DARK", "SKYFLAT")

#: APEX 가 색지수를 붙일 수 있는 다섯. up 은 3 장뿐이고 450 s 라 기본에서 뺐다.
DEFAULT_BANDS = ("B", "V", "rp", "ip", "zs")

#: 플랫은 필터마다 받아야 한다. bias·dark 는 필터와 무관하므로 한 번만 받는다.
FILTER_FREE = ("BIAS", "DARK")


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
    ap = argparse.ArgumentParser(
        description=f"LCO {INSTRUMENT} 의 M67 을 받는다 ({NIGHT} 밤)")
    ap.add_argument("--out", required=True, help="받을 폴더")
    ap.add_argument("--part", nargs="+", default=["science", "calib", "banzai"],
                    choices=["science", "calib", "banzai"],
                    help="받을 부분 (기본: 셋 다)")
    ap.add_argument("--bands", nargs="+", default=list(DEFAULT_BANDS),
                    help=f"밴드 (기본: {' '.join(DEFAULT_BANDS)})")
    a = ap.parse_args(argv)

    out = Path(a.out)
    bands = list(a.bands)
    parts = set(a.part)
    total_n = total_b = 0
    t0 = time.perf_counter()

    if "science" in parts:
        print(f"== M67 원본 (level 0 · {INSTRUMENT} · {NIGHT} 밤) ==", flush=True)
        for band in bands:
            rows = query("M67", 0, 300, INSTRUMENT, band, dayobs=NIGHT)
            n, b = _run(band, out / "science" / band, rows)
            total_n += n
            total_b += b

    if "calib" in parts:
        print(f"\n== 보정 프레임 ({NIGHT} 밤, 과학과 같은 밤) ==", flush=True)
        for ctype in CALIB_TYPES:
            if ctype in FILTER_FREE:
                rows = query(None, 0, 800, INSTRUMENT, None,
                             config_type=ctype, dayobs=NIGHT)
                n, b = _run(ctype, out / "calib" / ctype.lower(), rows)
                total_n += n
                total_b += b
                continue
            for band in bands:
                rows = query(None, 0, 800, INSTRUMENT, band,
                             config_type=ctype, dayobs=NIGHT)
                n, b = _run(f"{ctype} {band}",
                            out / "calib" / ctype.lower() / band, rows)
                total_n += n
                total_b += b

    if "banzai" in parts:
        print(f"\n== BANZAI 처리본 (level 91 · {NIGHT} 밤) ==", flush=True)
        for band in bands:
            rows = query("M67", 91, 300, INSTRUMENT, band, dayobs=NIGHT)
            n, b = _run(band, out / "banzai" / band, rows)
            total_n += n
            total_b += b

    dt = time.perf_counter() - t0
    print(f"\n합계 {total_n} 장 · {total_b/1e9:.2f} GB · {dt/60:.1f} 분", flush=True)
    print(f"받은 곳: {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
