# 성능 벤치마크·최적화 개발계획 (2026-08-07)

`APEX_ACTION_PLAN.md` 의 P2 "performance baseline" 항목과
`APEX_LC_PERFORMANCE_REVIEW.md` 의 "측정 계측부터" 순서를 실행 계획으로 편다.
논문은 이 작업의 배경이고, 산출물은 **수치가 나오는 벤치마크 체계**와
**parity 가 증명된 최적화**다. 작업 주체는 이 세션(Claude)이 Codex 로부터
인수했다 (2026-08-07 사용자 지시).

**전제 (불변)**
- 처리 순서: **공용체인 → CMD → LC** (사용자 결정)
- 검증된 compiled 패키지(SEP·photutils·Astropy WCS·SciPy)를 속도 이유로
  재구현하지 않는다 (사용자 지시 원문)
- baseline(`validation/BASELINE_2026-08-07.md`, `E:\APEX_validation\reprocess`)은
  보존하고, 모든 최적화는 그 결과에 대한 parity 로 판정한다
- 이전 cache 와 새 코드를 섞은 산출물은 논문 결과로 쓰지 않는다 —
  최종 clean run 만 쓴다
- `APEX_ACTION_PLAN.md` P0 의 자격증명 회전·저장소 경계는 **사용자 결정
  대기**로 이 계획 밖에 둔다

---

## Phase 0 — 계측 기반 (도구 없이는 수치도 없다) · 예상 반나절

| ID | 무엇을 만든다 | 그래서 무엇이 나온다 |
|---|---|---|
| T0.1 | `psutil` 을 선택 의존성 그룹으로 추가 (`pyproject [project.optional-dependencies] bench`) | 자식 프로세스 포함 peak RSS 를 잴 수단. 런타임 의존성은 늘지 않음 |
| T0.2 | `apex/benchmark/resources.py` — `measure()` 컨텍스트: wall time·peak RSS(자식 합산, 폴링 스레드)·worker 수·git commit·패키지 버전을 JSON 으로 | 모든 벤치가 같은 형식의 기계가독 기록을 남김 |
| T0.3 | `apex bench` CLI 서브커맨드 — 고정 fixture 정의(M67 30 프레임·NGC 6811 21 프레임·YZBoo 1밤)와 실행·기록. 결과는 `benchmark/runs/perf_<날짜>/` 누적 | 명령 한 줄로 재실행 가능한 벤치. Codex P2 "machine-readable results" 충족 |
| T0.4 | LC 로더 카운터 — `photometry_loader` 에 `frames_loaded`·`rows_loaded`·`cache_hit/miss` 계수기 (LC 리뷰 1단계 요구 그대로) | LC 최적화의 효과를 "몇 번 읽었나"로 직접 증명 |
| T0.5 | parity 게이트 명령 — 자기일관성 검사 + `compare_reprocess_runs` + 단계별 요약(검출수·WCS 수용·finite mask·apcorr)을 한 번에 | 최적화마다 즉흥 스크립트 없이 한 명령으로 판정 |

**오라클**: 계측 오버헤드 < 2 % (같은 fixture 계측 on/off 비교) ·
전체 pytest 통과 · 계측 JSON 이 커밋됨

## Phase 1 — before 수치 (현재 코드 그대로) · 예상 반나절 + 야간 배치

| ID | 측정 | 그래서 무엇이 나온다 |
|---|---|---|
| B1 | Step 0 peak RSS vs 프레임 수 (M67 부분집합 5/10/20/30) | "메모리가 프레임 수에 비례한다"가 그래프로 — streaming 의 근거이자 before 값 |
| B2 | 공용체인 worker 훑기 — NGC 6811 steps 4–7, workers 1/2/4/8/12/16, **결과 비트동일 확인 포함** | 논문 §5 스케일링 표의 원자료 + resource-aware worker 의 설계 근거 |
| B3 | LC cold/warm 로드 프로파일 — YZBoo step 9 빌드, T0.4 카운터로 | "프레임×별 반복 읽기"의 실제 횟수·시간 — LC 최적화의 before |
| B4 | **재현성 잡음 바닥** — 같은 입력 3회 반복: `n_sources` 편차·등급 MAD | parity 허용오차의 근거. 알려진 비결정성(M13 R 940→943)의 크기를 정량화 |

**오라클**: 네 측정 전부 `benchmark/runs/perf_baseline_*/` 에 JSON + 표로 커밋

## Phase 2 — 최적화 (한 번에 하나, parity 게이트) · 공용체인 먼저

| ID | 변경 | parity 게이트 (전부 통과해야 채택) | 예상 |
|---|---|---|---|
| O1 | **Streaming/two-pass Step 0** — 프레임 전량 float64 보유를 스트리밍 누적으로 | calibrated 프레임 **비트동일** · 마스터 비트동일 · peak RSS 가 프레임 수와 무관함을 B1 대비 수치로 | 1일 |
| O2 | **단계별 worker** — B2 실측: detect/wcs 는 2–4 에서 포화, **forcedphot 은 병렬이 4× 느림**(222.6→950 s @w=12) → 단계별 기본값 {detect/wcs: 2–4, forcedphot: 1} + RAM 승인 제어 | 통계 게이트(아래) 통과 · B2 대비 시간 개선을 수치로 · RAM 압박 시나리오 완주 | 반나절 + step7 병렬 경로 규명 별건 |
| O3 | **LC 공용 cache** — sid_map 공유 + frame×star compact matrix + bounded preload (LC 리뷰 P0 세 건, **마지막에**) | 광도값 자릿수 동일 · T0.4 카운터로 읽기 횟수 감소 수치 · cold build 시간 before/after | 1–2일 |

각 게이트에 **전체 pytest 통과**가 포함된다. 게이트 실패 = 그 최적화 폐기
(baseline 이 남아 있으므로 되돌림은 git revert 하나).

> **게이트 정정 (2026-08-07, B4 실측)**: step 4 이후 산출물은 비트동일이
> 성립하지 않는다 — 같은 입력·같은 worker 로도 소수 별 좌표가 1.5 px 까지
> 튀고 마스터 목록이 ±1 흔들린다. **근원 규명 완료: SEP 의 블렌드
> 디블렌딩**(`sep_flag=1` 인 밝고 붐비는 천체 약 2.7 %의 flux/shape 만
> 흔들리고 x·y·peak·FWHM 은 비트동일; RANSAC 가설은 반증 — 난수 미사용). 게이트는 실측 잡음 바닥 기반 통계형으로 바꾼다: 매칭 mag
> 중앙 MAD < 0.1 mmag · 바뀐 측정 < 5 % · 최대 |Δ| < 30 mmag · n_detect ±5 ·
> 목록 ±2. **step 0 마스터·calibrated 프레임은 비트동일 게이트 유지.**
> 수치 근거: benchmark/perf/20260807/RESULTS.md

### O2 설계 확정 (2026-08-07)

일반 도구의 "auto"는 거의 전부 CPU 수만 본다 — `multiprocessing.Pool`
기본값·pytest `-n auto`·`make -j` 전부. RAM 까지 보는 자동화는
dask(워커당 `memory_limit` + 80 % 일시정지·스필)나 배치 시스템(작업당
메모리 신청)처럼 큰 시스템에만 있고, **실행 중 동적 조절은 진동·복잡성
때문에 이 규모에서 과하다.** 일반 도구가 RAM auto 를 못 하는 이유는
작업당 메모리를 미리 모르기 때문인데, APEX 는 워크로드가 균질(같은 크기
프레임의 반복)해서 작업당 peak 를 실측할 수 있다 — 그래서 시작 시
승인 제어가 신뢰 가능하다.

```
workers = clamp(1, cpu상한(기존 75 %), 가용RAM × 0.6 ÷ 작업당 peak)
```

- 작업당 peak 는 추측이 아니라 **Phase 1 (B1/B2) 실측 계수**다
- I/O 무거운 단계는 별도 상한 — E: 외장 디스크에서 동시 읽기는 역효과이며,
  B2 스윕이 꺾이는 지점이 그 상한이다
- **auto 가 고른 값과 근거(가용 RAM·계수)를 `pipeline_run.json` 에 기록**한다
  (Codex P2 "worker counts" 항목). "auto 였다"로 끝나면 재현이 안 된다
- 결과는 worker 수와 무관해야 하고(비트동일 게이트), 따라서 auto 는 순수
  성능 손잡이다. 사용자 override(`max_workers`)는 유지한다
- 1차 방어선은 worker 수가 아니라 **작업 내부 메모리 바운드**
  (row-band combine — ccdproc `mem_limit` 와 같은 패턴)이고, worker 상한은
  2차다

CMD 는 별도 최적화 항목이 없다 — MCMC 는 계산 자체가 지배적이고
(LC 리뷰 P2 판정과 동일 논리) 재현성 seed 는 이미 있다. CMD 차례에는
Phase 3 clean run 의 CMD 산출물 검증으로 갈음한다.

## Phase 3 — 최종 clean run + 수치 고정 · 예상 반나절 + 야간

1. 채택된 코드로 raw 부터 전 대상 재실행 — **계측 켠 채로**
   (`reprocess_batch --out` 새 폴더, 씨앗 상수 가드 동작)
2. `BASELINE_2026-08-07.md` 와 같은 형식으로 기록: commit·환경·시간·
   **peak RSS(이번엔 있음)**·자기일관성
3. before/after 표 완성: RSS(B1↔O1)·worker 스케일링(B2↔O2)·
   LC 읽기 횟수(B3↔O3)
4. parity 확인 후 이 산출물이 논문 그림·수치의 **유일한 출처**가 된다

### O2 진행 상황 (2026-08-08)

| 항목 | 상태 |
|---|---|
| 단계별 상한 `{detect: 4, wcs: 4, crop: 4, forcedphot: 1}` | 구현·테스트 완료 (`46431eb`) |
| detect·forcedphot 상한의 실측 근거 | 유효 (B2) |
| **wcs 상한의 실측 근거** | **없음** — B2 의 wcs 행이 캐시 공유로 무효 |
| 구·신 정책 A/B 수치 | 재측정 중 (`benchmark/o2_ab.py`) |
| RAM 승인 제어 + 결정 근거 기록 | 미구현 |

**계측을 두 번 고쳤다.** 둘 다 수치의 해석을 바꾼 버그다.

1. **자식 프로세스 계측** (`4771192`) — Windows venv 의 `python.exe` 는 약 5 MB
   런처이고 진짜 인터프리터는 자식이다. `measure_command` 의 `peak_rss_mb` 가
   어느 실행이든 5.0 MB 를 보고했다. B2/B4 표는 처음부터 자식 합산 필드를
   썼으므로 영향 없다.
2. **`--result-dir` 가 캐시를 안 옮김** (`d8397ba`) — `read_params` 가
   `cache_dir` 을 설정의 `result_dir` 기준으로 굳힌 뒤에 `RunContext.build` 가
   `result_dir` 만 덮었다. 결과 폴더가 매번 새것이어도 검출·WCS 캐시는
   설정 트리에 공유된다 — 이 계획이 금지한 "이전 cache 와 새 코드를 섞은
   산출물" 을 헤드리스 경로가 구조적으로 만들고 있었다.

**추가 측정 (B2′)**: 캐시가 분리된 조건에서 worker 훑기를 다시 돌려
detect·wcs 의 진짜 병렬 거동을 잰다. forcedphot 은 재측정 불필요(캐시 없음,
게다가 B2 의 w=1 은 콜드였는데도 최속이라 결론이 과소평가 쪽).

## 병행 추적 (게이트 아님, 별건)

- **NGC 6811 7 프레임 화소 차이** 규명 — cosmetic 양쪽 켬인데 다름.
  B4 재현성 측정이 단서를 줄 것
- ~~검출 비결정성 원인 추적~~ → **규명 완료 (2026-08-07)**: SEP 디블렌딩.
  남은 것은 완화 여부 — `deblend_cont` 조정이나 블렌드 천체 제외가 과학적으로
  타당한지는 별도 판단이며, 현재는 그대로 두고 논문에 정직하게 서술한다

## 일정 요약

| 구간 | 예상 |
|---|---|
| Phase 0 계측 | 반나절 |
| Phase 1 before | 반나절 + 야간 배치 |
| Phase 2 (O1→O2→O3) | 3–4일 (게이트 포함) |
| Phase 3 clean run | 반나절 + 야간 |
| **합계** | **약 1주** |
