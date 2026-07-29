# -*- coding: utf-8 -*-
"""국문 논문 말뭉치 대조 — 하네스 규칙이 실제 국문 학술문체와 맞는지 잰다.

배경. 하네스 규칙 일부가 근거 없이 들어갔다(부록 Z). 특히 E4 등급인 것들:
  - 긴 대시 0 목표 (EASWA 에는 34건 있다)
  - 대구 'A가 아니라 B' 절당 1회 (숫자에 근거가 없다)
사용자 지적: *"뭐 국문논문이라도 찾아서 보던가"*

말뭉치 셋으로 비교한다.
  A. 국문 기술논문  — refs_ko/ 의 관측·원격탐사 논문. **저자가 Claude 가 아니다**
  B. EASWA         — 사용자 본인 교육논문(미완)
  C. APEX 교정 전   — Claude 가 쓴 것

A 는 「국문 학술문체가 실제로 어떻게 쓰이는가」의 기준이고, B 는 「이 사용자가
어떻게 쓰는가」다. 둘 다 0 인데 C 만 있으면 규칙화해도 된다.

    python -X utf8 corpus_compare.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).parent
REFS_KO = HERE / "refs_ko"
EASWA = Path(r"C:\Users\bmffr\Desktop\2026_Cosmos_ERP\EASWA_논문_전체수정본_v11.md")

# 하네스가 검사하거나 검사할 뻔한 것들
PATTERNS = [
    ("긴 대시 —",        r"—"),
    ("대구 A가 아니라 B", r"[가-힣\w)]+(?:이|가) 아니라 "),
    ("몫으로 남",        r"몫으로 남"),
    ("~로 이어진다",      r"[으]?로 이어진다"),
    ("핵심은",           r"핵심은"),
    ("효과적",           r"효과적"),
    ("주목할 점은",       r"주목할 점은"),
    ("부족·필요한 것은",   r"(?:부족|필요)한 것은"),
    ("음을 보인 사례",     r"음을 보[인여](?:주는)? 사례"),
    ("~로 볼 수 있다",    r"[으]?로 볼 수 있다"),
    ("~을 시사한다",      r"[을를] 시사한다"),
    ("3항 이상 가운뎃점",  r"[가-힣A-Za-z0-9]{2,12}(?:·[가-힣A-Za-z0-9]{2,12}){2,}"),
]


def load_ko_corpus() -> tuple[str, list[str]]:
    """refs_ko 의 PDF 에서 한글 본문을 뽑아 붙인다. 스캔본(한글 3,000자 미만)은 뺀다."""
    try:
        from pypdf import PdfReader
    except ImportError:
        sys.exit("pypdf 가 필요하다:  .venv-deploy\\Scripts\\python -m pip install pypdf")
    texts, used = [], []
    for p in sorted(REFS_KO.glob("*.pdf")):
        try:
            r = PdfReader(str(p))
            t = "\n".join((pg.extract_text() or "") for pg in r.pages)
        except Exception:
            continue
        if len(re.findall(r"[가-힣]", t)) < 3000:
            continue                      # 스캔본
        if p.stem == "astro_edu":
            continue                      # 교육논문이라 B 와 같은 장르
        texts.append(t)
        used.append(p.name)
    return "\n".join(texts), used


def main() -> int:
    ko, used = load_ko_corpus()
    if not ko:
        sys.exit("refs_ko 에 쓸 수 있는 국문 논문이 없다")
    easwa = EASWA.read_text(encoding="utf-8") if EASWA.exists() else ""
    before = HERE / "_apex_before.md"
    apex = before.read_text(encoding="utf-8") if before.exists() else \
        (HERE / "MANUSCRIPT_ko.md").read_text(encoding="utf-8")

    corp = [("국문 기술논문", ko), ("EASWA(본인)", easwa), ("APEX(Claude)", apex)]
    print(f"■ 말뭉치")
    for name, t in corp:
        print(f"   {name:14s} {len(t):8,}자")
    print(f"   국문 기술논문 구성: {', '.join(used)}\n")

    head = f"{'패턴':20s}" + "".join(f"{n:>16s}" for n, _ in corp) + "   판정"
    print(head)
    print("─" * len(head))
    for label, rx in PATTERNS:
        d = []
        for _, t in corp:
            n = len(re.findall(rx, t))
            d.append(n / len(t) * 10000 if t else 0.0)
        ko_d, ea_d, ap_d = d
        base = max(ko_d, ea_d)
        if base == 0 and ap_d > 0:
            v = "규칙화 가능 — Claude 만 쓴다"
        elif base == 0 and ap_d == 0:
            v = "전부 0 — 근거 없음"
        elif ap_d > base * 1.8:
            v = f"경고만 (Claude {ap_d/base:.1f}배)"
        else:
            v = "규칙에서 뺄 것 — 국문 학술문체다"
        print(f"{label:20s}" + "".join(f"{x:14.2f}/만" for x in d) + f"   {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
