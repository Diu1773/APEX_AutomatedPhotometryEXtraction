# Phase 1 before 수치 — 2026-08-07 야간 배치

계측: `apex bench` (`8bf6d51`), psutil 폴링 0.1 s, 자식 프로세스 합산.
환경: Python 3.12.3 · 논리코어 16 · numpy 2.4.4 (각 JSON 의 `env` 블록이 정본).

## B1 — Step 0 peak RSS vs 프레임 수 (M67 raw 부분집합, in-process)

| lights | wall (s) | peak RSS (MB) |
|---:|---:|---:|
| 5 | 259.7 | 5,041 |
| 10 | 337.6 | 5,161 |
| 20 | 516.4 | 5,263 |
| 30 | 643.3 | 5,204 |

**가설이 뒤집혔다.** "메모리가 light 프레임 수에 비례한다"는 계획의 전제였는데,
실측 peak RSS 는 **5.0–5.26 GB 로 평평**하다. light 는 이미 한 장씩 처리되고
있고, 지배 항은 **보정 풀 로딩** — 매 실행이 bias 102 장·dark 196 장·flat 18 장을
스캔해 마스터를 만들며, 그 결합 입력(한 그룹의 전체 프레임 리스트)이 통째로
메모리에 올라간다. 첫 재처리에서 봤던 1.71 GiB 할당 실패(sigmaclip 작업 배열)도
같은 계열이다.

wall time 은 lights 에 선형: **약 187 s 고정(스캔+마스터) + 15.4 s/frame**.

**O1 streaming 의 표적이 바뀐다** — light 경로가 아니라 **마스터 결합의 입력
경로**다. 그룹 프레임을 전량 로드하는 대신 행 밴드 단위로 파일에서 직접 읽어
누적하면(PTC 의 `HDU.section` 기법과 동일 계열) peak 를 수백 MB 로 묶을 수 있다.
결합 연산 자체는 이미 밴드 단위(`_COMBINE_CHUNK_BYTES`)라 산술은 그대로다.

**O1 판정선(before 확정)**: peak RSS 5.0–5.3 GB → 목표 < 1 GB, calibrated
프레임·마스터 비트동일.

한계: n ≤ 30, 보정 풀 고정(M67 fixture). light 수백 장 규모에서 다른 항이
나타날 가능성은 clean run 계측이 확인한다.

## B2 — worker 훑기 · B4 — 재현성 (1차 시도 실패 → 재실행)

첫 배치는 여섯 worker 값 전부 `Config not found` 로 즉사했다. 원인은 벤치가
가리킨 `parameters_20260807.toml` 이 존재하지 않았던 것 — 트리 승격 때
`cp … 2>/dev/null` 이 조용히 실패했다 (clean run 의 설정은 TOML 이 아니라
JSON 정본 `apex_config.json` 이었다). 교훈: **오류를 버리는 복사는 실패도
버린다.** 다섯 대상의 JSON 을 `apex_config_20260807.json` 으로 제대로 복사하고
벤치 경로를 고쳐 재실행했다. 결과는 `sweep_NGC6811.json`·`repro_NGC6811.json`
(재실행본)이 정본이다.
