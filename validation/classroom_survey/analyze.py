"""수업 설문 응답을 논문 표로 만든다.

원자료는 `responses.tsv`(익명), 문항 전문은 `items.tsv`.
본문 수치를 하드코딩하지 않고 여기서 뽑는다 — 자료가 바뀌면 표도 바뀐다.

    python -X utf8 validation/classroom_survey/analyze.py
"""
from __future__ import annotations

import csv
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
LIKERT_MAX = 5

# 짧은 표 라벨. 원문은 items.tsv 에 그대로 있다.
LABEL = {
    "AX1": "단계형 구성이 분석 절차 이해에 도움",
    "AX2": "이미지에서 측광·CMD·광도곡선까지 따라갈 수 있었다",
    "AX3": "자동 측광과 결과 시각화가 유용",
    "AX4": "구경 측광으로 별 밝기를 재는 과정을 이해",
    "AX5": "여러 단계를 한 프로그램에서 해서 편리",
    "AX6": "실제 관측자료 탐구 활동에 활용 가능",
    "AI1": "천문 이미지 자료를 확인·이해하는 데 도움",
    "AI2": "화면 구성과 조작이 쉬웠다",
    "AI3": "이미지 밝기·품질·전처리 필요성을 확인",
    "AI4": "원본이 전처리를 거쳐 분석 자료가 되는 과정을 이해",
    "AI5": "전문 프로그램을 처음 접하는 학생도 접근하기 쉬웠다",
    "AI6": "실제 관측자료 수업·탐구에 활용 가능",
    "OV1": "직접 찍은 자료를 분석하는 활동이 흥미를 높였다",
    "OV2": "관측·전처리·측광·분석·해석의 연결을 이해",
    "OV3": "사용법 설명이 활동을 수행하기에 충분했다",
    "OV4": "학교 수업·동아리·연구 활동에 활용 가능",
    "OV5": "전체 수업 경험에 만족",
}
GROUPS = [("APEX", "AX"), ("AstralImage", "AI"), ("전체 활동", "OV")]


def load() -> list[dict[str, str]]:
    with (HERE / "responses.tsv").open(encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def numeric(rows: list[dict[str, str]], key: str) -> list[float]:
    out = []
    for r in rows:
        try:
            out.append(float(r[key]))
        except (KeyError, TypeError, ValueError):
            pass
    return out


def main() -> None:
    rows = load()
    n = len(rows)
    print(f"응답 {n} 명 · 척도 1~{LIKERT_MAX}\n")

    prior = numeric(rows, "prior")
    novice = sum(1 for v in prior if v <= 2)
    print("배경")
    print(f"  사전 분석 경험 (1 = 없음): 중앙 {statistics.median(prior):.0f} · "
          f"{novice}/{len(prior)} 명이 2 이하")
    used = [r["programs"] for r in rows]
    both = sum(1 for u in used if u.startswith("둘"))
    print(f"  두 프로그램 모두 사용: {both}/{n} 명")
    if both < n:
        print(f"  ! 한 명은 AstralImage 만 썼다고 답했으나 APEX 문항에도 응답했다 — "
              f"APEX 표의 유효 N 을 {n} 로 볼지 {both} 로 볼지 본문에 밝힌다")

    for title, prefix in GROUPS:
        keys = [k for k in LABEL if k.startswith(prefix)]
        print(f"\n{title}")
        print(f"  {'항목':<48} {'N':>2} {'평균':>5} {'중앙':>4} {'최소':>4} {'최대':>4}")
        for k in sorted(keys):
            v = numeric(rows, k)
            print(f"  {LABEL[k]:<48} {len(v):>2} {statistics.mean(v):>5.2f} "
                  f"{statistics.median(v):>4.1f} {min(v):>4.0f} {max(v):>4.0f}")
        allv = [x for k in keys for x in numeric(rows, k)]
        print(f"  {'— 묶음 평균':<48} {len(allv):>2} {statistics.mean(allv):>5.2f}")

    print("\n자유 서술 (응답한 사람만)")
    for key, title in (("AI_free", "AstralImage"), ("AX_free", "APEX"),
                       ("OV_free", "수업 전체")):
        texts = [r[key] for r in rows if r.get(key, "").strip()]
        print(f"\n  {title} — {len(texts)}/{n} 명")
        for t in texts:
            print(f"    · {t}")


if __name__ == "__main__":
    main()
