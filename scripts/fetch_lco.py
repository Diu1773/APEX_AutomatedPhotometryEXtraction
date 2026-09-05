"""LCO 공개 아카이브에서 프레임을 골라 내려받는다.

다른 기기 자료로 APEX 를 끝까지 돌리기 위한 도구다(STATUS.md 의 목표 둘).
인증이 필요 없고 표준 라이브러리만 쓰므로 다른 기계에 그대로 옮겨 써도 된다.

먼저 무엇이 있는지 본다:

    python -X utf8 scripts/fetch_lco.py --target M67 --list

기기와 필터를 골라 받는다:

    python -X utf8 scripts/fetch_lco.py --target M67 --instrument sq30 \
        --filter rp --level 0 --limit 20 --out E:/APEX_validation/external/LCO_M67

`--level 0` 이 raw 이고 `91` 이 LCO 의 BANZAI 가 처리한 것이다. 둘 다 받아 두면
APEX 의 Step 0 을 BANZAI 와 화소 단위로 견줄 수 있다(Fig 13 이 그 방식이다).

이미 받은 파일은 건너뛰므로 중간에 끊겨도 다시 돌리면 이어진다.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

API = "https://archive-api.lco.global/frames/"
UA = {"User-Agent": "APEX-fetch/1.0 (photometry pipeline validation)"}


def _get(url: str) -> dict:
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))


def query(target: str | None, level: int, limit: int, instrument: str | None,
          filt: str | None, config_type: str | None = None,
          start: str | None = None, end: str | None = None) -> list[dict]:
    """Walk the paginated frame list, newest first, up to `limit` rows."""
    params = {"public": "true", "reduction_level": str(level),
              "limit": "100"}  # 익명 사용자는 한 쪽에 100 이 상한이다
    if target:
        params["target_name"] = target
    if instrument:
        params["instrument_id"] = instrument
    if filt:
        params["primary_optical_element"] = filt
    if config_type:
        params["configuration_type"] = config_type
    if start:
        params["start"] = start
    if end:
        params["end"] = end
    url = API + "?" + urllib.parse.urlencode(params)
    rows: list[dict] = []
    while url and len(rows) < limit:
        page = _get(url)
        rows.extend(page.get("results", []))
        url = page.get("next")
    return rows[:limit]


def show(rows: list[dict]) -> None:
    if not rows:
        print("맞는 프레임이 없다.")
        return
    tally = Counter()
    for r in rows:
        tally[(r.get("instrument_id"), r.get("telescope_id"), r.get("site_id"),
               r.get("primary_optical_element"), round(float(r.get("exposure_time") or 0)))] += 1
    print(f"{'기기':>7} {'망원경':>7} {'관측소':>5} {'필터':>10} {'노출(s)':>8} {'개수':>5}")
    print("-" * 50)
    for (inst, tel, site, fl, exp), n in sorted(tally.items(), key=lambda kv: -kv[1]):
        print(f"{str(inst):>7} {str(tel):>7} {str(site):>5} {str(fl):>10} {exp:>8} {n:>5}")
    print(f"\n합계 {len(rows)} 프레임 (질의 상한까지)")


def download(rows: list[dict], out: Path) -> int:
    out.mkdir(parents=True, exist_ok=True)
    total = 0
    for i, r in enumerate(rows, 1):
        name = r.get("filename") or f"{r.get('id')}.fits.fz"
        dest = out / name
        if dest.exists() and dest.stat().st_size > 0:
            print(f"[{i}/{len(rows)}] 이미 있음 {name}")
            continue
        url = r.get("url")
        if not url:
            print(f"[{i}/{len(rows)}] 내려받기 주소 없음 {name}")
            continue
        tmp = dest.with_suffix(dest.suffix + ".part")
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=300) as resp, tmp.open("wb") as fh:
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
        tmp.replace(dest)
        mb = dest.stat().st_size / 1e6
        total += dest.stat().st_size
        print(f"[{i}/{len(rows)}] {name}  {mb:.1f} MB", flush=True)
    print(f"\n받은 양 {total/1e9:.2f} GB → {out}")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", help="천체 이름 (예: M67). 보정 프레임을 받을 때는 생략한다")
    ap.add_argument("--config-type", dest="config_type",
                    help="EXPOSE(기본) · BIAS · DARK · SKYFLAT")
    ap.add_argument("--start", help="관측 시작일 (예: 2021-03-01)")
    ap.add_argument("--end", help="관측 종료일 (예: 2021-04-01)")
    ap.add_argument("--level", type=int, default=0, help="0=raw, 91=BANZAI 처리본")
    ap.add_argument("--instrument", help="기기 코드 (예: sq30)")
    ap.add_argument("--filter", dest="filt", help="필터 (예: rp, V)")
    ap.add_argument("--limit", type=int, default=200, help="최대 프레임 수")
    ap.add_argument("--list", action="store_true", help="받지 않고 목록만 본다")
    ap.add_argument("--out", type=Path, help="저장 폴더 (받을 때 필요)")
    a = ap.parse_args(argv)

    if not a.target and not a.config_type:
        ap.error("--target 이나 --config-type 중 하나는 있어야 한다")
    rows = query(a.target, a.level, a.limit, a.instrument, a.filt,
                 a.config_type, a.start, a.end)
    if a.list or not a.out:
        show(rows)
        if not a.list:
            print("\n받으려면 --out 을 준다.")
        return 0
    show(rows)
    print()
    return download(rows, a.out)


if __name__ == "__main__":
    raise SystemExit(main())
