# -*- coding: utf-8 -*-
"""STYLE_HARNESS.md §3 「인용 범위」 중 기계로 확인 가능한 부분만 검사한다.

- 본문 인용 키가 references.bib 에 있는가
- 저자 수와 et al. 표기가 맞는가 (2인은 "A & B", 1인은 단독)
- 인용 매크로 밖에 저자명을 손으로 적어 두지 않았는가

인용이 그 문장을 실제로 뒷받침하는지(범위 초과)는 사람이 원문을 봐야 한다.
"""
from __future__ import annotations

import re
from pathlib import Path

HERE = Path(__file__).parent
BIB = HERE.parent / "references.bib"
MD = HERE / "MANUSCRIPT_ko.md"

md = MD.read_text(encoding="utf-8")
bib = BIB.read_text(encoding="utf-8")

keys = set()
for m in re.finditer(r"\\cite[a-z]*(?:\[[^\]]*\])*\{([^}]+)\}", md):
    keys.update(k.strip() for k in m.group(1).split(","))

ents = {}
for m in re.finditer(r"@\w+\{([^,]+),(.*?)\n\}", bib, re.S):
    ents[m.group(1).strip()] = m.group(2)

print(f"본문 인용 키 {len(keys)}개 / bib 항목 {len(ents)}개")
missing = sorted(k for k in keys if k not in ents)
print("bib 에 없는 키:", missing or "없음")

print("\n■ 저자 수 대비 표기")
for k in sorted(keys & set(ents)):
    au = re.search(r"author\s*=\s*[{\"](.*)", ents[k])
    if not au:
        print(f"  [!] {k}: author 필드 없음")
        continue
    raw = au.group(1)
    depth, buf = 0, []
    for ch in raw:                       # 중괄호 균형까지 읽는다
        if ch == "{":
            depth += 1
        elif ch == "}":
            if depth == 0:
                break
            depth -= 1
        buf.append(ch)
    n = len(re.split(r"\s+and\s+", "".join(buf).strip()))
    if n == 1:
        print(f"  [단독] {k}: 저자 1인 — 'et al.' 금지")
    elif n == 2:
        print(f"  [2인 ] {k}: 'A & B' 로 표기해야 함")

print("\n■ 인용 매크로 밖에 손으로 적은 저자명")
found = False
for m in re.finditer(r"(?<!\{)\b([A-Z][A-Za-z]+)\s+(?:et\s+al\.|&\s+[A-Z][A-Za-z]+)", md):
    ctx = md[max(0, m.start() - 30):m.end() + 20].replace("\n", " ")
    if "\\cite" in md[max(0, m.start() - 60):m.start()]:
        continue
    print(f"  {m.group(0)}   …{ctx}…")
    found = True
if not found:
    print("  없음")
