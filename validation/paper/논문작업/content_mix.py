# -*- coding: utf-8 -*-
"""원고를 「천문 / 도구검증 / 소프트웨어공학 / 배경·목적」으로 나눠 센다.

투고 전략(D-P01)을 정하려면 "천문이 얼마고 공학이 얼마인지"를 알아야 한다.
절 단위 분류는 Claude 가 본문을 읽고 손으로 매겼다. 자동 판정이 아니므로
근거를 한 줄씩 붙였다. 이견이 있으면 CLASS 표만 고치면 다시 계산된다.

  ASTRO  천문 내용. 하늘·천체·측광량이 주어이고, 결론이 그 대상에 대한 것.
  TOOL   도구 검증. 천문 방법으로 재지만 결론은 소프트웨어 산출값에 대한 것.
  ENG    소프트웨어 공학. 구조·설정·재현성·성능.
  FRAME  배경·목적·대상·한계. 어느 저널이든 필요한 틀.
"""
import re
from pathlib import Path

MD = Path(__file__).parent / "MANUSCRIPT_ko.md"

CLASS = {
    "1. 서론":               ("FRAME", "도구 계보와 연구 목적. 천문 결과 없음"),
    "2.1":  ("FRAME", "설계 요구사항 세 가지. 왜 이렇게 만들었나"),
    "2.2":  ("ENG",   "계층 구조와 Qt 분리"),
    "2.3":  ("ENG",   "외부 라이브러리 대 자체 구현 구분. 표 1"),
    "2.4":  ("ASTRO", "검출·WCS·카탈로그·강제 측광 절차 자체. 측광 방법 서술"),
    "2.5":  ("ASTRO", "CMD/LC 분석 단계. PSF·SYSREM·PDM"),
    "2.6":  ("FRAME", "사용자가 각 단계에서 내리는 결정. 대상 독자 논거"),
    "2.7":  ("ENG",   "TOML 설정·JSON 상태·캐시 서명"),
    "3.1":  ("FRAME", "검증 접근과 자료 명세. 방법론 틀"),
    "3.2":  ("TOOL",  "합성 주입으로 보정 되돌리기"),
    "3.3":  ("ASTRO", "PTC 로 gain·읽기잡음·암전류 실측. 검출기 물리량"),
    "3.4":  ("TOOL",  "ccdproc 대비 비트동일. 구현 일치 확인"),
    "3.5":  ("TOOL",  "BANZAI 대비 교차기기. 구현 일치 확인"),
    "3.6":  ("ASTRO", "실측 프레임 인공별 주입 완전도. 깊이는 하늘·시상이 정한다"),
    "3.7":  ("TOOL",  "astrometric solution 정확도"),
    "3.8":  ("ASTRO", "측광 오차 모형. sigma_m = 1.0857/SNR 검증"),
    "3.9":  ("ASTRO", "구경·문턱값·하늘밝기·시상 민감도. 관측 조건 의존성"),
    "3.10": ("TOOL",  "SEP 대비 합성 참값 일치"),
    "3.11": ("TOOL",  "IRAF/DAOPHOT 대비 실측 일치. 파라미터 정합표"),
    "3.12": ("ASTRO", "M5·M13 코어 밀집도 대 측광법. 성단 관측 결과"),
    "3.13": ("ASTRO", "PS1 교차대조. Gaia BP faint 결함 규명"),
    "3.14": ("ASTRO", "두 측광계 CMD 재현. 주계열 능선 일치"),
    "3.15": ("TOOL",  "SYSREM·PDM 모듈 정확성"),
    "3.16": ("TOOL",  "프레임 QC 게이트 동작"),
    "4. 과학 적용": ("ASTRO", "YZ Boo 광도 곡선과 주기. 유일한 과학 산출물 절"),
    "5.1":  ("FRAME", "검증 결과의 범위"),
    "5.2":  ("FRAME", "미검증 구성요소"),
    "5.3":  ("FRAME", "일반성의 한계. 단일 기기 문제"),
    "5.4":  ("FRAME", "적용 지침"),
    "6. 결론": ("FRAME", "결론"),
}

LABEL = {"ASTRO": "천문", "TOOL": "도구검증", "ENG": "소프트웨어공학", "FRAME": "배경·목적·한계"}


def main():
    s = MD.read_text(encoding="utf-8")
    parts = re.split(r"\n(#{2,3}) ", s)
    rows = []
    for i in range(1, len(parts) - 1, 2):
        body = parts[i + 1]
        title = body.split("\n")[0].strip()
        if not re.match(r"^\d", title):
            continue
        n = len(body)
        if n < 60:           # 하위절만 있는 상위 절 머리
            continue
        # 번호 토큰으로 정확히 맞춘다. startswith 로 하면 "3.10" 이 "3.1" 에 걸린다.
        num = re.match(r"^(\d+(?:\.\d+)?)", title).group(1)
        key = next((k for k in CLASS
                    if re.match(r"^(\d+(?:\.\d+)?)", k).group(1) == num), None)
        if key is None:
            print(f"  [분류 없음] {title}")
            continue
        cat, why = CLASS[key]
        rows.append((title, n, cat, why))

    tot = sum(n for _, n, _, _ in rows)
    print(f"■ 본문 {tot:,}자 (초록·후기 제외)\n")
    for t, n, cat, why in rows:
        print(f"  {LABEL[cat]:8s} {n:5,}  {100*n/tot:4.1f}%  {t[:34]:36s} {why}")

    print()
    agg = {}
    for _, n, cat, _ in rows:
        agg[cat] = agg.get(cat, 0) + n
    for cat in ("ASTRO", "TOOL", "ENG", "FRAME"):
        v = agg.get(cat, 0)
        print(f"  {LABEL[cat]:16s} {v:6,}자  {100*v/tot:5.1f}%  {'█'*int(100*v/tot/2)}")
    print(f"\n  천문 + 도구검증 = {100*(agg.get('ASTRO',0)+agg.get('TOOL',0))/tot:.1f}%")
    print(f"  순수 소프트웨어공학 = {100*agg.get('ENG',0)/tot:.1f}%")


if __name__ == "__main__":
    main()
