# TRACK: apex-paper — APEX 논문 RASTI 투고

> 개발 트랙은 별도: [TRACK.md](TRACK.md). 오라클이 다르므로 분리했다.
> **이 트랙은 야간 배치 대상이 아니다.** 주장·구조·해석은 사람 판단이다.

## 완료 정의

RASTI(또는 PASP)에 투고 완료.

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

- 막힌 것: 없음. 다만 **§4 가 936자로 유독 얇다** — 다른 절의 1/6 이다.
  의도된 축약인지 미완인지는 사람 판단.

## 사용자 의견

- **APEX 논문 §3 을 지금 쓸지, 재검출 결과를 기다릴지** — 논문써 근데 내가 fig들 엄청 지적했던거 알아?
  `2026-07-27 23:37` 대시보드 결정. 카드 원문: `Main/DECISIONS.md`
  → ⚠️ **결정 시점에 §3 은 이미 쓰여 있었다.** 이 파일이 낡아서 "비어 있음"으로 보였다.
  → 살아 있는 요구는 **fig 쪽**이다. 지적 원문 179건을 캐냈다:
  `Main/harvest/corrections_fig.md` · 요약은 `Main/TASTE.md` T-001~003.
  핵심 한 줄: **"그림은 자기가 무엇으로 만들어졌는지 그림 안에서 말해야 한다"**
  (데이터 출처 · prior · 파라미터 범위 · Fig.N 캡션). 2026-06-21 · 07-01 · 07 세 번 나왔다.

## 다음 3개

1. **fig 49장을 T-001 기준으로 점검** — 각 그림이 「무슨 데이터로 돌렸는지」를
   그림 안에 담고 있는가. 원문은 `Main/harvest/corrections_fig.md`.
   담기지 않은 그림 목록부터 뽑는다. 이건 기계 판정 가능하다 (캡션 문자열 검사).
2. 포지셔닝 문단 확정 — ①깊이예측·③union마스터는 기존 계보 편입, **②검출수 상한이 신규**
3. `references.bib` 43건 최종 검증 (`/ars-citation-check`)

## 사용자 판단 필요

### D-001 · 투고처 RASTI vs PASP 최종 확정
- 무엇: APEX 논문의 투고처를 RASTI(RAS Techniques and Instruments)와 PASP(Publications of the ASP) 중 확정
- 지금: 원고는 이미 `mnras.cls`(MNRAS/RASTI 공용 OUP 템플릿)로 43,254자·section 6개·그림 49장까지 완성. 투고처만 미정
- 선택지: A. RASTI / B. PASP
- 근거:
  · **범위 적합성** — RASTI의 공식 scope가 "software for processing data (including pipelines), data analysis and modelling"을 명시적으로 포함한다 (academic.oup.com/rasti/pages/why-publish, 2026-07 확인). APEX는 검증 섹션이 전체 원고의 42%(18,060자/11개 하위절)를 차지하는 파이프라인 검증 논문 — RASTI가 명시한 범위와 정확히 겹친다. PASP는 범용 천문학지에 계측·소프트웨어 코너가 있는 형태로, 과학 성과(§4, YZ Boo)의 비중이 더 요구될 가능성이 있다.
  · **포맷 비용** — 원고가 이미 `mnras.cls`로 작성됨. RASTI는 MNRAS와 동일한 OUP LaTeX 클래스를 그대로 받아 재포맷이 0. PASP는 AASTeX 필수 — 43k자·49그림 전체를 재조판해야 한다.
  · **비용** — RASTI APC는 £1,339(비회원 기준, RAS 회원 20% 할인, LMIC 국가 전액 면제) 단일 트랙(완전 OA). PASP는 OA $3,490, 또는 구독 트랙 페이지당 $129 + 부록 $175 — 49개 그림이 딸린 원고는 페이지 수가 많아 구독 트랙도 결코 저렴하지 않다 (iopscience.iop.org/journal/1538-3873/page/publication-charges, 2026-07 확인). RASTI가 확실히 저렴하다.
  · **지표/인지도** — PASP는 IF 7.22(2024)·SJR 2.517·Q1로 인지도가 높은 종합 천문학지다. RASTI는 SJR 0.856, 2021년 창간이라 전통 Impact Factor가 아직 없다 (다만 Web of Science ESCI·ADS·DOAJ·Scopus에는 색인됨). 이 항목은 PASP 우위다.
  · **심사 소요 기간** — 양쪽 모두 공식 게재까지 소요일수를 찾지 못함 (찾지 못함).
- 권고: **RASTI.** scope 명문화 일치 + 재포맷 비용 0 + APC 절반 이하, 세 축이 전부 RASTI를 가리킨다. PASP의 우위는 인지도(IF/SJR)뿐인데, 이건 저자가 감내할 수 있는 트레이드오프다. 단, 공저자(지도교수)가 인지도를 우선순위로 두면 PASP로 뒤집을 근거가 된다 — 그 판단은 사람 몫.
- 되돌리기: 만약 RASTI 심사에서 "과학 성과 비중 부족"으로 반려되면, §4(YZ Boo, 현재 936자로 다른 절의 1/6)를 확장한 뒤 PASP로 재투고를 검토한다.

- 신규성 주장 강도 — "세 파이프라인 전부 초과검출을 안 본다"를 어디까지 밀지

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
