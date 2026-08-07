# Native 구현의 탄생 이력과 근거 등급 (git 고고학)

작성 2026-08-07. `APEX_COMPONENT_AUDIT.md` 가 **현재 상태의 지도**라면, 이 문서는
**시간축** — 각 자체 구현이 언제, 어떤 이유로 태어났는지를 git 이력에서 복원해
근거 등급을 붙인다. 논문 §2.3 의 "왜 자체 구현했나" 열과 Section 15 심사 방어의
원천 자료다.

**처리 순서 (사용자 결정 2026-08-07): 공용체인 → CMD → LC.**
LC 후반(주기분석)은 사용자 이해가 선행되어야 하므로 마지막에 다룬다.

## 증거 등급

- **(a) 당시 기록** — 커밋 메시지·문서에 동기가 남아 있어 그대로 인용 가능
- **(b) 현재 증명** — git 은 침묵하지만 지금의 검증·측정·문헌 준거로 방어됨
- **(c) 증거 없음** — 기억뿐. 논문에서 빼거나 측정으로 (b) 승격해야 함

주의: `34bb792` ("chore: commit untracked source so CI …", 2026-06-11)와
`5ce4a8b` ("backup", 2026-05-08)는 **벌크 커밋**이다 — 여기서 태어난 파일의
실제 작성 시점은 그보다 앞이고, git 에 동기가 없다.

## 공용체인 (Step 0–7)

| 모듈 | 탄생 | 당시 커밋이 말하는 동기 | 유형 | 등급 |
|---|---|---|---|---|
| `calibration.py`·`cosmetic.py`·`overscan.py` | 2026-07-09 `3258ca7` | detector calibration 통합 (raw→science) | 시설 파이프라인 부재 환경 대응 | **(a)+(b)** — ccdproc 비트동일·IRAF·BANZAI 검증 |
| quad solver — `quad_matcher.py` | 2026-06-10 `79005b9`·`3a2190f` | **"perf: vectorize quad canonical-code" · "batch the kd-tree query"** | 속도 (git 에 실재) | **(a)** |
| quad solver — `solver.py` | `34bb792` 벌크 | 없음 | 의존성 제거 (실행파일·인덱스·네트워크 불요) | **(b)** — §3.6 ASTAP·astrometry.net 대조 |
| `wcs_solve.py` | 2026-04-28, 24커밋 | "parallelize internal solver"(06-01) · "header→target→blind hint chain"(06-10) · "source-set ladder"(06-11) | 견고성 사다리 + 병렬화 | **(a)** |
| `detection.py`·`forced_photometry.py` | 06-21 헤드리스화 | **"perf(mem): float32 로 프레임당 RAM 절반"**(`c5cfe62` 07-10) · "per-frame m50 + depth QC gate"(`d92be0a`) | 메모리 + 깊이 게이트 | **(a)** — float32 절반은 산술적으로 자명 |
| `frame_qc.py` | 07-02 `f12280a` | "physically calibrated **depth-cost gate**"(`5e2a6df`) · **"fail frames that detect far more sources than the group"**(`0f01faa` 07-26 — 초과검출 게이트, 논문 신규성) | 신규 QC 방법 | **(a)** |
| `photometric_qc.py` | 07-02 `c72b331` | "post-photometry transparency" | 투명도 QC | **(a)** |
| `extinction.py` (SVD 투영 소광) | `34bb792` 벌크 | 없음 | z_j–k1 축퇴 제거 (수학적 필요) | **(b)** — CLAUDE.md 수학노트 + 검증 |
| `fast_stats.py` (Bottleneck 래퍼) | `5ce4a8b` 벌크 | 없음 (docstring: "optional Bottleneck acceleration", numpy 폴백) | 속도 (nan-통계) | **(b)** — 2026-08-07 실측: 결합 축(axis 0, 30×140×4788 float32)에서 nanmedian **2.3×**·nanmean **3.8×**, 전체 축약은 이득 없음(0.9×). **측정 중 `bn.nanstd` 가 float32 전체 프레임에서 +21 % 어긋나는 잠재 결함 발견·수정** (axis=None 시 float64 승격, 회귀 테스트 7건). 프로덕션 경로는 원래 무사 — `finite_*` 는 float64 승격, 결합은 축 길이 ≤41 |
| `detector_ptc.py` | 2026-08-06 `aa85439` | 제조사 미공표·헤더 EGAIN 14배 오류 | 필요성 | **(a)** — `DETECTOR_CONSTANTS.md` 완비 |

## CMD (Step 8–12)

| 모듈 | 탄생 | 당시 커밋이 말하는 동기 | 유형 | 등급 |
|---|---|---|---|---|
| `isochrone_mcmc.py` | 06-21 `ea32c63`, 8커밋 | "EEP interpolation — **fixes age recovery**"(`a5ebdd1`) · "full inter-color covariance"(`2fa1f3a`) · "Gaia parallax prior"(`ff480ad`) · "reproducible MCMC (seed)"(`556fbed`) | 과학적 필요 (축퇴·재현성) | **(a)** — 메모리 `project_isochrone_fit_flaw` 에 전체 서사 |
| `standard_anchor.py` | 08-04 `45d6f0a` | U/B/V 영점을 한 표준계에 앵커해야 U−B 가 축퇴를 품 | 과학적 필요 | **(a)** — 소거실험 문서화 |
| `isochrone_fitter_v2.py` | 04-28 벌크 | 없음 | — | (b) — MCMC 로 대체·검증됨 |

## LC (Step 8–11) — 마지막 처리

| 모듈 | 탄생 | 당시 커밋이 말하는 동기 | 유형 | 등급 |
|---|---|---|---|---|
| `sysrem.py` | `34bb792` 벌크 | 없음 (docstring: Tamuz+2005 MNRAS 구현 명시) | 부재 (표준 라이브러리에 없음) | **(b)** |
| `period_analysis_service.py` | 04-28, 8커밋 | "numerical accuracy — WLS·SYSREM matrix·error propagation·BLS"(`4861df4`) · **"speed up bootstrap FAP"**(`63344b3`) · "iterative pre-whitening (**Period04/Breger**)"(`018369d`) · "resolve the multi-night **1-day alias** via reference-night selection"(`05c883f`) | 부재(PDM) + 속도 + 문헌 준거 | **(a)** |
| `period_alias_service.py` | 07-11 `b6c77f7` | 1일 별칭 서비스 (WIP 표기) | 별칭 진단 | (a) — WIP 이력 주의 |
| `global_ensemble.py`·`loader.py`·`merge/id_match.py` | 04-28 벌크 | 없음 | 공용 자료구조 | (b) — YZ Boo 재현으로 간접 |

**주기분석 이해용 메모**: 사용자가 아직 이 부분을 이해하지 못한다고 밝혔으나,
git 서사 자체가 학습 자료다 — Period04/Breger 준거(pre-whitening), WLS/오차전파
수정, 기준밤 선택에 의한 1일 별칭 해소가 커밋 단위로 남아 있다. 이해 정리는
이 커밋들을 따라가면 된다.

## 요약 — 논문이 지금 쓸 수 있는 것과 필요한 것

**(c) 로 남는 항목은 없다.** 마지막 후보였던 `fast_stats` 는 2026-08-07 실측으로
(b) 승격했고, 그 측정이 잠재 결함(`bn.nanstd` float32 전체 축약 +21 %)까지
적발해 수정으로 이어졌다 — 근거 정렬 작업이 코드 품질 문제를 직접 잡아낸
사례다. 나머지 자체 구현은 전부 당시 기록 또는 현재 검증으로 방어된다.

**속도를 근거로 쓸 수 있는 모듈** (git 에 perf 커밋 실재):
quad vectorize/batch(`79005b9`·`3a2190f`) · 내장 솔버 병렬화(`d572aae`) ·
forced-phot float32(`c5cfe62`) · bootstrap FAP(`63344b3`).
단, **수치 주장(몇 배)은 여전히 미측정** — 구조 서술("vectorized", "parallel")로
쓰거나, 필요 시 worker 훑기·마이크로벤치로 채운다.
`APEX_PERFORMANCE_AUDIT.md` 의 "committed baseline 없음" 판정과 일치한다.

**벌크 커밋 3인방** (`solver.py`·`extinction.py`·`sysrem.py`): git 동기는 없지만
셋 다 (b) — 외부 엔진 대조·수학적 필요·문헌 구현 명시로 방어된다. 사후 합리화가
필요한 모듈은 없다.
