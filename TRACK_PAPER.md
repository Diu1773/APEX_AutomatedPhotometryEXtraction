# TRACK: apex-paper — APEX 논문 RASTI 투고

> 개발 트랙은 별도: [TRACK.md](TRACK.md). 오라클이 다르므로 분리했다.
> **이 트랙은 야간 배치 대상이 아니다.** 주장·구조·해석은 사람 판단이다.

## 완료 정의

**RASTI** 에 투고 완료. (투고처 확정 2026-07-28 01:42 — D-001, `Main/DECISIONS.md`)

## 오라클 — 추적가능성

기계로 돌리는 테스트는 없다. 대신 이 조건을 만족해야 완료다.

- 본문의 모든 수치가 `validation/paper/fig*.py` 출력으로 **역추적 가능**
- 표·그림·본문 숫자가 서로 어긋나지 않음
- 인용은 실제 서지에서 확인된 것만 (**DOI 조작 금지** — 기존 규약)

권장 기법: 수치를 원고에 하드코딩하지 말고 스크립트가 매크로 파일을 생성하게 한다.
데이터가 바뀌면 본문이 자동 갱신되고, 본문과 표가 어긋나는 사고가 원천 차단된다.

## 지금

- 마지막 커밋: `4b6b079` 2026-07-28 00:07 — paper(tex): discussion reflects real-frame depth + QC gate
- **§3 은 이미 쓰였다.** 이 파일이 9커밋 동안 "§3만 비어 있음"으로 낡아 있었다 (F-019).
  채운 커밋: `f240804` 19:20 · `6a8bcaa` 19:33 · `4b6b079` 00:07

  2026-07-28 01:xx 실측 — `validation/paper/MANUSCRIPT.tex` 43,254자 / section 6개:

  | 절 | 분량 | 하위절 |
  |---|---:|---:|
  | 1. Introduction | 5,052자 | 0 |
  | 2. Design and implementation | 6,640자 | 4 |
  | **3. Validation** | **18,060자** | **11** |
  | 4. Science application (YZ Boo) | 936자 | 0 |
  | 5. Discussion | 6,302자 | 4 |
  | 6. Conclusion | 1,968자 | 0 |

  운영 QC 절(`sec:depthqc`)은 수치가 박힌 완성 산문이다 — 3성단 **57 프레임**,
  예측 대 실현 깊이 **RMS 78 mmag**, 전 프레임이 **±0.5 mag** 게이트 안,
  `S/N50 = 4.05` 고정, 마스터 한계에 근접한 **9 프레임**은 circularity guard 로 제외.

- **2026-07-28 01:50 세션 발견 (`6e7860d`) — 투고본이 한국어 초안보다 6개 하위절 뒤진다.**
  `MANUSCRIPT.tex` 검증 하위절 11개 vs `MANUSCRIPT_ko.md` 16개.
  투고본에 없는 것: Step 0 검출기 보정 · 검출기 특성화 · ccdproc 교차검증 ·
  교차기기 보정(LCO Sinistro/QHY600 vs BANZAI) · astrometry · 시계열 코어.
  그 세션의 표현: *"raw-to-science in one tool 주장의 증거 기반 전체"* 가 빠져 있다 —
  **개별 그림 수정보다 이게 우선이다.**
- fig1 LCO 확장은 **기각** (`6e7860d`): 논문 스스로 교차기기 일반성을 범위 밖으로
  선언했으므로 반쪽 일반성 주장만 부른다. 같은 카메라 안에서 강화한다 (조건 확대 또는 M5/M3).
  그림 재작업 계획: `validation/paper/논문작업/FIGURE_REBUILD_PLAN.md` (R1~R6 반복 실패 회고 포함)
- 막힌 것: 없음. 다만 **§4 가 936자로 유독 얇다** — 다른 절의 1/6 이다.
  의도된 축약인지 미완인지는 사람 판단.

## 사용자 의견

- **투고처 RASTI vs PASP 최종 확정** — A안 — RASTI / B. PASP  
  `07-28 01:42` 대시보드 결정. 카드 원문: `Main/DECISIONS.md`
- **APEX 논문 §3 을 지금 쓸지, 재검출 결과를 기다릴지** — 논문써 근데 내가 fig들 엄청 지적했던거 알아?
  `2026-07-27 23:37` 대시보드 결정. 카드 원문: `Main/DECISIONS.md`
  → ⚠️ **결정 시점에 §3 은 이미 쓰여 있었다.** 이 파일이 낡아서 "비어 있음"으로 보였다.
  → 살아 있는 요구는 **fig 쪽**이다. 지적 원문 179건을 캐냈다:
  `Main/harvest/corrections_fig.md` · 요약은 `Main/TASTE.md` T-001~003.
  핵심 한 줄: **"그림은 자기가 무엇으로 만들어졌는지 그림 안에서 말해야 한다"**
  (데이터 출처 · prior · 파라미터 범위 · Fig.N 캡션). 2026-06-21 · 07-01 · 07 세 번 나왔다.

## 다음 3개

1. **ko→tex 이식: 빠진 검증 하위절 6개** — Step 0 보정 · 검출기 특성화 · ccdproc ·
   교차기기 보정 · astrometry · 시계열 코어. `MANUSCRIPT_ko.md` 가 원본이다.
   (`6e7860d` 판단: 개별 그림 수정보다 우선. LCO 데이터의 자리도 여기다)
2. **fig 재작업 — `FIGURE_REBUILD_PLAN.md` 의 실행 순서대로.** T-001 기준
   (그림이 자기 데이터 출처를 그림 안에서 말하는가) 점검을 겸한다.
   원문 179건: `Main/harvest/corrections_fig.md`
3. `references.bib` 43건 최종 검증 (`/ars-citation-check`)

## 사용자 판단 필요

## 함정

- **신규성 근거는 확정됐다.** Dragonfly(하한만) · ZTF(개수 기반 없음) · iPTF/PTFIDE(상한 없음)
  — 세 파이프라인 전부 초과검출을 보지 않는다. 이 공백이 게이트의 신규성이다.
- 문헌 확인은 IOPscience/ar5iv에서 직접. PDF가 안 열리면 ar5iv를 쓴다(iPTF 사례).
- 그림은 `validation/paper/figures/` 49개가 git 추적 중. `data*/`는 제외(627 MB).

## 핵심 문헌

| 키 | 역할 |
|---|---|
| `smtn002` | Rubin 단일방문 깊이. 한계flux ∝ σ_sky·FWHM — APEX `depth_ref`와 동일 스케일링 |
| `stetson1994` | ALLFRAME. "상당 비율 프레임에 나타나는 것만 인정" = min-frames 선례 |
| `danieli2020` | Dragonfly. `NOBJ < 1000` 게이트로 프레임 40~60% 폐기 |

## 최근 세션 원문

`C:\Users\bmffr\Desktop\Main\harvest\apex\` (논문 관련: 2026-07-09/12 PAPER.md 세션들)
