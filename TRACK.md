# TRACK: apex-dev — APEX 소프트웨어 v1.0 릴리즈

> 이 파일이 이 트랙 상태의 **유일한 진실**이다. 세션을 열면 여기부터 읽는다.
> 세션을 끝내기 전에 `## 지금` / `## 다음 3개` / `## 함정`을 반드시 다시 쓴다.
> 논문 트랙은 별도: [TRACK_PAPER.md](TRACK_PAPER.md)

## 완료 정의

측광 파이프라인이 실관측 데이터에서 GUI·헤드리스 양쪽으로 완주하고,
회귀 테스트 전량 통과 상태로 v1.0 릴리즈 태그가 붙는다.

## 오라클 — 이 트랙의 검증 명령

```bash
.venv-deploy/Scripts/python.exe -m pytest tests/ -q
```

- **통과 기준:** 614 passed, 0 failed
- **소요:** 약 7분
- **마지막 실행:** 2026-07-27 — 614 passed, 6 warnings
- CI: `.github/workflows/tests.yml` 외 3개

## 지금

- 마지막 커밋: `95af5b2` 2026-07-27 — docs(harness): add harness protocol section to AGENTS.md
- 하네스 설치됨 (2026-07-27): AGENTS.md에 세션 종료 프로토콜, Stop hook이 이 파일 갱신을 검사
- 미푸쉬: 0 · 미커밋: 0 (2026-07-27 하네스 구축 세션에서 정리)
- 진행 중: 프레임별 예측 검출한계 QC 게이트 완성됨. 원고 §3 반영 대기
- 막힌 것: 없음

## 다음 3개

1. `ref_union_min_frames` 기본값 검토 — ALLFRAME(stetson1994) 선례가 확보돼 근거는 생겼으나,
   **마스터에 유령이 실제로 있는지 먼저 재는 것이 순서**
2. 4000프레임 재검출 — 원본 계수(`n_raw_detections`)가 이제 기록되므로 재실행하면 실제 숫자가 남음
3. `source_quality.py:115` All-NaN slice 경고 처리 — roundness 계산에서 발생 (테스트 4건)

## 사용자 판단 필요

- `ref_union_min_frames` 기본값을 올릴지 — 유령 측정 결과를 보고 결정
- `validation/`을 최종적으로 어떻게 보관할지 (현재 git 제외, 디스크 20 GB 상주)

## 함정

- **`validation/` 은 커밋하지 않는다.** 20,088 MB / 10,990 파일. `.gitignore`에 있음.
  논문 자산(`validation/paper/figures`, `captions`, `논문작업`, `fig*.py`)만 추적한다.
- **실행은 `run.bat` 또는 `.venv-deploy\Scripts\python`.** 시스템 python 아님.
- **`step10_zeropoint_calibration.py:2199, 2804`** — pandas `DataFrameGroupBy.apply` 폐기 예정.
  pandas를 올리면 **영점보정 결과가 조용히 달라질 수 있다.** 올리기 전에 회귀 스냅샷을 뜰 것.
- 초과검출 게이트 임계는 실측 산포(실프레임 194장/24그룹, 최대 1.53배, 99분위 1.49)에서
  나왔다. REVIEW 1.8 / FAIL 2.5. **임의 숫자가 아니므로 함부로 바꾸지 말 것.**
- `n_sources`가 아니라 `n_raw_detections`를 봐야 한다. 전자는 필터·상한 뒤 값이라
  과검출이 정상처럼 보인다.

## 최근 세션 원문

`C:\Users\bmffr\Desktop\Main\harvest\apex\` (15건, 2026-05~07)
