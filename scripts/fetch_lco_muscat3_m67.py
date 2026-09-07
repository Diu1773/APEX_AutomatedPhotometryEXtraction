"""LCO MuSCAT3 의 M67 한 밤을 원본부터 보정 프레임까지 통째로 받는다.

**왜 이것인가.** APEX 는 M67 을 Moravian C3-61000 으로 이미 처리해 두었다
(`validation/paper/data_realframe_M67g_broad` · `M67r_mid` · `M67i`). MuSCAT3 은
2 m 망원경(ogg 2m0a)에 붙어 **gp·rp·ip·zs 네 밴드를 동시에** 찍고, 그중 셋이
APEX 의 g·r·i 와 겹친다. **같은 성단 · 같은 필터 · 다른 기기**라서 표준성야보다
직접적인 대조가 된다(STATUS.md 의 목표 둘).

**무엇을 받나.** 2021-03-20 밤 하나. 재서 확인한 크기는 아래와 같다(한 장 3.4 MB).

    M67 raw   gp 76 · rp 76 · ip 76 · zs 77  =  305 장  1.04 GB
    BIAS      네 카메라 각 64                =  256 장  0.87 GB
    DARK      네 카메라 각 20                =   80 장  0.27 GB
    SKYFLAT   네 카메라 각 12                =   48 장  0.16 GB
                                               ---------------
                                               689 장  2.34 GB

보정 프레임까지 받는 이유는 **Step 0(원본 → 과학 프레임)을 여기서도 돌리기
위해서다.** 2026-09-02 갤럭시북 실행은 보정 프레임을 못 찾아 Step 0 을 건너뛰었고,
그래서 「원본부터 끝까지」라고 말할 수 없었다.

인증이 필요 없고 표준 라이브러리만 쓴다. 이미 받은 파일은 건너뛰므로 끊겨도 다시
돌리면 이어진다.

실행:
    python -X utf8 scripts/fetch_lco_muscat3_m67.py --out D:/LCO_M67_MuSCAT3
    python -X utf8 scripts/fetch_lco_muscat3_m67.py --out ... --only-science
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).absolute().parent))

from fetch_lco import download, query  # noqa: E402

NIGHT_START, NIGHT_END = "2021-03-20", "2021-03-21"
SCIENCE_START, SCIENCE_END = "2021-03-01", "2021-04-01"

# MuSCAT3 은 카메라마다 밴드가 고정이다.
CAMERAS = (("ep04", "gp"), ("ep02", "rp"), ("ep03", "ip"), ("ep05", "zs"))
CALIB_TYPES = ("BIAS", "DARK", "SKYFLAT")


def _run(label: str, out: Path, rows: list) -> tuple[int, int, float]:
    if not rows:
        print(f"  [{label}] 맞는 프레임이 없다.", flush=True)
        return 0, 0, 0.0
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    n_bytes = download(rows, out)
    dt = time.perf_counter() - t0
    print(f"  [{label}] {len(rows)} 장 · {n_bytes/1e9:.2f} GB · {dt:.0f} 초"
          + (f" · {n_bytes/1e6/dt:.1f} MB/s" if dt > 0 else ""), flush=True)
    return len(rows), n_bytes, dt


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="MuSCAT3 의 M67 한 밤을 통째로 받는다")
    ap.add_argument("--out", required=True, help="받을 폴더")
    ap.add_argument("--only-science", action="store_true",
                    help="보정 프레임은 건너뛰고 M67 만 받는다")
    ap.add_argument("--only-calib", action="store_true",
                    help="M67 은 건너뛰고 보정 프레임만 받는다")
    ap.add_argument("--bands", nargs="+", default=None,
                    help="받을 밴드 (기본: 넷 다). 예: gp rp ip")
    a = ap.parse_args(argv)

    out = Path(a.out)
    want = set(a.bands) if a.bands else {b for _, b in CAMERAS}
    cams = [(c, b) for c, b in CAMERAS if b in want]

    total_n = total_b = 0
    t0 = time.perf_counter()

    if not a.only_calib:
        print(f"== M67 원본 ({SCIENCE_START} ~ {SCIENCE_END}) ==", flush=True)
        for cam, band in cams:
            rows = query("M67", 0, 300, cam, None,
                         start=SCIENCE_START, end=SCIENCE_END)
            n, b, _ = _run(f"{cam} {band}", out / "science" / band, rows)
            total_n += n
            total_b += b

    if not a.only_science:
        print(f"\n== 보정 프레임 ({NIGHT_START} 밤) ==", flush=True)
        for ctype in CALIB_TYPES:
            for cam, band in cams:
                rows = query(None, 0, 400, cam, None, config_type=ctype,
                             start=NIGHT_START, end=NIGHT_END)
                n, b, _ = _run(f"{ctype} {cam}",
                               out / "calib" / ctype.lower() / band, rows)
                total_n += n
                total_b += b

    dt = time.perf_counter() - t0
    print(f"\n합계 {total_n} 장 · {total_b/1e9:.2f} GB · {dt/60:.1f} 분", flush=True)
    print(f"받은 곳: {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
