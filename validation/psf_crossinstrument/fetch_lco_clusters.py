# -*- coding: utf-8 -*-
"""LCO 공개 아카이브에서 성단 raw + BANZAI 마스터를 받는다 (재개형).

교차기기 CMD 검증용 두 세트:

  M45  sq32 (0.4m tfn, QHY600 **풀프레임** 2.0°x1.3°) 2025-01-14 · B10+V10 · 120 s
       → 광시야 EPSF 공간변화 + 성단 CMD. 우리 Moravian(0.43°)의 4.6배 시야.
  M67  sq32 (같은 기기 central 30'x30') 2024-10-29 · B5+V5 · 60 s
       → 우리 Moravian M67(2026-02-08)과 **같은 성단 다른 기기** CMD 대조.

이미 받은 파일은 크기가 맞으면 건너뛴다(중단해도 재실행하면 이어받음).

실행: .venv-deploy\\Scripts\\python -X utf8 validation/psf_crossinstrument/fetch_lco_clusters.py
"""
from __future__ import annotations

import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path(r"E:\APEX_validation\psf_crossinstrument")
API = "https://archive-api.lco.global/frames/"

SETS = {
    "m45_wide": dict(target="M45", inst="sq32", day="2025-01-14"),
    "m67_lco": dict(target="M67", inst="sq32", day="2024-10-29"),
    # U+B+V 가 같은 밤에 있는 유일한 M67 세트 (LCO 0.4m / kb74 = SBIG STL-6303).
    # U-B 색이 있어야 나이-금속함량 축퇴가 풀린다 — B-V 만으로는 [M/H] 가
    # 사전분포를 1.7σ 밀어내고 금속결핍 쪽으로 rail 한다(REPORT_M67_CROSS.md §2).
    "m67_ubv": dict(target="M67", inst="kb74", day="2015-04-01",
                    filters=("U", "B", "V")),
}

DEFAULT_FILTERS = ("B", "V")


def q(**kw) -> dict:
    kw.setdefault("public", "true")
    with urllib.request.urlopen(API + "?" + urllib.parse.urlencode(kw), timeout=90) as r:
        return json.loads(r.read())


def frame_by_id(fid) -> dict:
    with urllib.request.urlopen(f"{API}{fid}/", timeout=60) as r:
        return json.loads(r.read())


def fetch(url: str, dest: Path) -> tuple[int, bool]:
    """(bytes, downloaded?) — 이미 완전한 파일이면 건너뛴다."""
    req = urllib.request.Request(url, headers={"Range": "bytes=0-0"})
    with urllib.request.urlopen(req, timeout=90) as r:
        cr = r.headers.get("Content-Range", "")
        total = int(cr.split("/")[-1]) if "/" in cr else 0
    if dest.exists() and total and dest.stat().st_size == total:
        return total, False
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    t0 = time.time()
    with urllib.request.urlopen(url, timeout=600) as r, tmp.open("wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)
    tmp.replace(dest)
    dt = max(time.time() - t0, 1e-6)
    size = dest.stat().st_size
    print(f"    {dest.name}  {size/1e6:.1f}MB  {size/dt/1e6:.1f}MB/s")
    return size, True


def run_set(key: str, spec: dict) -> None:
    base = OUT / key
    filters = tuple(spec.get("filters", DEFAULT_FILTERS))
    print(f"\n=== {key}: {spec['target']} / {spec['inst']} / {spec['day']} "
          f"/ {'+'.join(filters)} ===")

    sci = [
        r for r in q(target_name=spec["target"], INSTRUME=spec["inst"],
                     DAY_OBS=spec["day"], RLEVEL="00", OBSTYPE="EXPOSE",
                     limit="100")["results"]
        if r["FILTER"] in filters
    ]
    e91 = [
        r for r in q(target_name=spec["target"], INSTRUME=spec["inst"],
                     DAY_OBS=spec["day"], RLEVEL="91", OBSTYPE="EXPOSE",
                     limit="100")["results"]
        if r["FILTER"] in filters
    ]
    counts = {f: sum(1 for r in sci if r["FILTER"] == f) for f in filters}
    print(f"  raw e00: " + " ".join(f"{f}{n}" for f, n in counts.items())
          + f"  |  e91 대조: {len(e91)} (필터당 1장만 받음)")

    masters: dict[str, dict] = {}
    for f in filters:
        ref = next((r for r in e91 if r["FILTER"] == f), None)
        if ref is None:
            continue
        for rid in ref.get("related_frames", []):
            m = frame_by_id(rid)
            if m.get("OBSTYPE") in ("BIAS", "DARK", "SKYFLAT"):
                masters[m["basename"]] = m

    jobs: list[tuple[str, Path]] = []
    for r in sci:
        jobs.append((r["url"], base / "raw" / r["filename"]))
    seen_filt = set()
    for r in e91:  # 보정 산술 대조용 필터당 1장
        if r["FILTER"] not in seen_filt:
            seen_filt.add(r["FILTER"])
            jobs.append((r["url"], base / "banzai_ref" / r["filename"]))
    for m in masters.values():
        jobs.append((m["url"], base / "masters" / m["filename"]))

    print(f"  파일 {len(jobs)}개 (raw {len(sci)} · e91 {len(seen_filt)} · master {len(masters)})")
    got = skipped = 0
    for i, (url, dest) in enumerate(jobs, 1):
        try:
            _sz, dl = fetch(url, dest)
        except Exception as exc:  # 개별 실패는 건너뛰고 계속 — 재실행하면 이어받는다
            print(f"    [{i}/{len(jobs)}] FAIL {dest.name}: {exc}")
            continue
        if dl:
            got += 1
        else:
            skipped += 1
    total_mb = sum(p.stat().st_size for p in base.rglob("*") if p.is_file()) / 1e6
    print(f"  완료: 새로 받음 {got} · 기존 {skipped} · 총 {total_mb:.0f}MB → {base}")


if __name__ == "__main__":
    keys = sys.argv[1:] or list(SETS)
    for k in keys:
        run_set(k, SETS[k])
    print("\ndone.")
