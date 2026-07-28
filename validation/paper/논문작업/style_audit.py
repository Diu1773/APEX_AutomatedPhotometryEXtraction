# -*- coding: utf-8 -*-
"""STYLE_HARNESS.md §6 의 내부 검사를 기계로 돌린다.

사람이 판정할 것(논리 전개·문단 역할·인용 범위)은 여기서 다루지 않는다.
문자열로 잡히는 것만 잡는다 — 대신 예외 없이 전부 잡는다.

    python -X utf8 style_audit.py [원고.md]
"""
from __future__ import annotations

import re
import sys
from collections import OrderedDict
from pathlib import Path

HERE = Path(__file__).parent
DEFAULT = HERE / "MANUSCRIPT_ko.md"

# ─────────────────────────────────────────────────────────────────────────
# 검사 항목. (심각도, 이름, 정규식, 설명)
# 심각도는 하네스 §7 의 S1~S4 를 따른다.
# ─────────────────────────────────────────────────────────────────────────
CHECKS: "OrderedDict[str, tuple]" = OrderedDict()


def check(sev, name, pattern, why):
    CHECKS[name] = (sev, re.compile(pattern), why)


# S2 — 근거 없는 장면·인물·경험 (하네스 §3)
check("S2", "가상 인물·장면",
      r"어느 학생|한 학생|어떤 연구자|한 연구자|사용자는[^.]{0,20}느[낀끼]|"
      r"학생조차|초심자가[^.]{0,15}[쥐잡]",
      "관찰·조사·인용이 없는 사용 장면은 쓸 수 없다")
check("S2", "측정 없는 수행 서술",
      r"한나절|버튼 하나|클릭 몇 번|쓸 만한|손에 [쥐잡]|금세|눈 깜짝",
      "수행시간·결과 품질은 조사 자료 없이 주장할 수 없다")

# S3 — 영어 에세이 번역투 (하네스 §4.1)
check("S3", "번역투 은유",
      r"막힘없이|인상에 남는다|손에 [쥐쥔]|짐은 증거|물리에 정박|"
      r"갉아먹|거짓말하는|가로질러|믿음으로 받아들|소유하지 않는다|"
      r"벽을 없애|문턱을 낮[추춘]|잡음 속에서|눈이 되어|숨을 쉬",
      "은유·의인화는 기술적 서술로 바꾼다")

# S3 — AI 해설 상투구 (하네스 §4.2)
check("S3", "AI 해설구",
      r"핵심은|요점은|중요한 점은|주목할 점은|이 절이 묻는|바로 이 점이|"
      r"다시 말해|흥미롭게도|무엇보다|단순한 .{1,12}를 넘어|"
      r"의미는 분명하다|말해 준다\.",
      "결과 뒤에 의미 부여 문장을 매번 붙이지 않는다")

# S3 — 대구·슬로건 (하네스 §4.3)
check("S3", "대구 'A가 아니라 B'",
      # '아니라' 뒤에 쉼표가 오는 경우를 놓치던 허점을 막았다(2026-07-28).
      r"[가-힣\w)]+[이가] 아니라[,\s]|에 그치지 않고|"
      r"뿐만 아니라 .{1,20}까지|에서 .{1,12}(으)?로 전환",
      "절당 1회를 넘기지 않는다. 남길 자리는 두 문장 대조로 바꾼다")

# S2/S3 — 저자의 자기평가 (하네스 §4.4)
check("S2", "자기평가·홍보어",
      r"정직[한하]|가장 강[한하]|튼튼[한하]|인상적|직관적|검증된 (파이프라인|도구|"
      r"소프트웨어)|확립된|실제 하늘에서|우수[한함]|성공적|강력[한하]|"
      r"완벽[한하]|손색이 없",
      "평가는 수치가 허용하는 범위까지만 쓴다")

# S4 — 형식 번역투 (하네스 §4.5)
check("S4", "긴 대시", r"—", "문장을 나누거나 괄호를 쓴다")

# S4 — 억지 번역어 (2026-07-28 추가)
#
# 이 검사가 없어서 서론을 "문체 1건"으로 통과시켰다. 실제로는 용어가 전부
# 지어낸 한국어였다. 사용자 지적: *"광곡선이 아니라 광도 곡선 또는 그냥 light
# curve라고 해; 원시 프레임도 그냥 raw 프레임이라고 하던가, 억지로 한국어로
# 번역을 하지말라고"*
#
# 판별 기준은 「한국어 표준 표기가 실제로 있는가」다. 있으면 그것을 쓰고
# (aperture=구경), 없는데 지어냈으면 영문 또는 통용 표기로 되돌린다.
check("S4", "억지 번역어",
      r"광곡선|원시 프레임|원시 CCD|원시 영상|천체측정 해|짝짓기|짝지은|"
      r"되찾음|되찾기|되찾은|하늘값|밀집장|성장곡선|광자전달곡선|"
      r"마스터 천체 목록|오버스캔 트림|하늘 미리보기|축소 파이프라인|"
      r"축소 엔진|축소 도구|검출 문턱(?!값)|고정 문턱(?!값)",
      "표준 표기가 없으면 영문을 쓴다. 한국어를 새로 만들지 않는다")

# 「3항 병렬」검사는 뺐다 (2026-07-28).
# 하네스 §4.5 가 막으려는 것은 "리듬 맞추려고 억지로 세 개를 채우는" 습관인데,
# `A·B·C` 정규식은 그걸 재지 못한다. 4~5항 나열의 앞 세 개를 잘라 세고,
# 반복되는 고정 용어(`분석·시각화·기준값` 6회)를 6건으로 센다.
# 실측: 가운뎃점 나열의 항 수 분포가 APEX 2항 81.0%/3항 17.2%,
#       EASWA 2항 83.6%/3항 12.5% 로 사실상 같다. 둘 다 3항 편중이 없다.
# 지표가 습관이 아니라 나열 길이에 반응하므로 부적합하다(전역 철칙 12).
# 이 습관을 정말 재려면 항 수 분포를 봐야 한다 — 필요하면 그때 따로 만든다.

# S4 — 같은 용어의 괄호 병기 반복 (첫 등장만 허용)
REPEAT_GLOSS = re.compile(r"([가-힣]{2,12})\(([A-Za-z][A-Za-z \-]{2,30})\)")


def scan(text: str):
    lines = text.splitlines()
    hits = OrderedDict((k, []) for k in CHECKS)
    for i, ln in enumerate(lines, 1):
        if ln.startswith(("#", "|", "```", "<!--")):
            continue
        for name, (sev, rx, why) in CHECKS.items():
            for m in rx.finditer(ln):
                s = max(0, m.start() - 26)
                hits[name].append((i, m.group(0), ln[s:m.end() + 26].strip()))
    # 괄호 병기 반복
    seen, repeat = {}, []
    for i, ln in enumerate(lines, 1):
        for m in REPEAT_GLOSS.finditer(ln):
            key = m.group(1) + "|" + m.group(2).strip().lower()
            if key in seen:
                repeat.append((i, m.group(0), f"{seen[key]}행에서 이미 병기"))
            else:
                seen[key] = i
    return hits, repeat


def main():
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT
    text = src.read_text(encoding="utf-8")
    hits, repeat = scan(text)

    n_para = len([p for p in text.split("\n\n") if len(p) > 80])
    print(f"■ {src.name} — {len(text):,}자 / 문단 {n_para}개\n")

    total = 0
    for name, (sev, _rx, why) in CHECKS.items():
        h = hits[name]
        total += len(h)
        mark = "OK " if not h else f"{sev} "
        print(f"[{mark}] {name:16s} {len(h):4d}건   {why}")
        for line, tok, ctx in h[:4]:
            print(f"        {line:>5}: …{ctx}…")
        if len(h) > 4:
            print(f"        … 외 {len(h) - 4}건")
    print(f"[{'OK ' if not repeat else 'S4 '}] {'괄호 병기 반복':16s} {len(repeat):4d}건   "
          f"전문용어는 처음 한 번만 병기한다")
    for line, tok, ctx in repeat[:4]:
        print(f"        {line:>5}: {tok}  ({ctx})")
    total += len(repeat)

    print(f"\n합계 {total}건.  S1(사실·인용 오류)은 기계로 못 잡는다 — 사람이 본다.")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
