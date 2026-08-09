# 선행연구 위치잡기 — 도구 수준 대안과 노벨티 판정

작성 2026-08-09. `APEX_COMPONENT_AUDIT.md` 는 **프리미티브 수준**(astropy·
photutils·ccdproc·SEP)에서 "왜 패키지를 통째로 쓰지 않았나"를 이미 답했다.
비어 있던 것은 **도구 수준** — 완제품 파이프라인과의 비교다. 심사자가
"왜 ASteCA 를 안 쓰나", "AutoPhOT 과 뭐가 다른가"를 물었을 때 답할 근거다.

판정 규칙은 `APEX_COMPONENT_AUDIT.md` 의 것을 그대로 쓴다 —
**조립되었다는 이유만으로 novel 이 되지 않는다.**

---

## 1. 도구 수준 지형도

| 도구 | 출처 | 범위 | APEX 와 겹치는 구간 | 확인 |
|---|---|---|---|---|
| **AstroImageJ** | Collins+ 2017, AJ 153, 77 | GUI 데스크톱(Java/ImageJ). **bias·dark·flat·비선형 보정 포함**, 시계열 다중구경 차등측광, 앙상블, 동시 추세제거·트랜싯 적합. 측성은 **astrometry.net 웹 인터페이스** | **Step 0–7 + LC 전체.** 형태가 가장 가깝다 | 확인(본문) |
| **PhotometryPipeline** | Mommert 2017 (ascl:1703.004) | "small to medium-sized observatories" 대상 자동 파이프라인. SExtractor + SCAMP 로 정합·측광, 온라인 카탈로그로 측광 보정 | Step 4–7 + **대상 사용자층이 동일** | 확인 |
| **AutoPhOT** | Brennan & Fraser 2022, **A&A 667, A62** | 초신성·트랜지언트용. 구경+PSF 측광, 템플릿 차감, 인공성 주입 한계등급 | Step 4–7 | 확인 |
| **PhoPS** | Erece & Kilic, arXiv:2607.27414 (2026-07-29 제출, **아직 arXiv**) | 측성+측광 자동화. **관측 시각으로 전파한 Gaia DR3 인덱스를 동적 생성** → 사전 설치 인덱스 불요. RANSAC 기반 위치의존 영점 | Step 5 + 7 | 확인 |
| **ASteCA** | Perren+ 2015, **A&A 576, A6** (현행 v0.7) | 성단 구조·멤버십·합성 CMD·우도. **입력이 이미 만들어진 측광 카탈로그** — 이미지 처리를 하지 않는다 | **Step 12 만** | 확인 |
| `isochrones` | Morton 2015 (ascl:1503.010) | 항성모델 격자 인터페이스. **개별 별**의 질량·나이·거리·소광 추정 | 성단 CMD 동시 피팅의 직접 대체물 아님 (격자 계층으로는 대체 가능) | 확인 |
| `PySysRem` | Tamuz+ 2005 알고리즘의 파이썬 구현 | SYSREM | LC 추세제거 | 확인 |
| `wotan` | Hippke+ 2019, AJ 158, 143 | 시계열 추세제거 종합. **트랜싯 탐색용 단일 광도곡선** 겨냥 | 다중밤 앙상블과 겨냥이 다름 | 확인 |
| `astropy.timeseries` | — | Lomb–Scargle, BLS | **PDM 은 없다** | 확인 |
| VaST · Astrokit | Sokolovsky+ · — | 변광성 탐색·고정밀 차등측광 | LC | 미확인 |

---

## 2. APEX 구성요소별 최근접 대안과 판정

판정 3단: **주장 가능** / **워크플로 기여** / **주장 불가**

| APEX 구성요소 | 최근접 도구 수준 대안 | 판정 | 근거 |
|---|---|---|---|
| Step 0 검출기 보정 (raw→science) | ccdproc (프리미티브) · AstroImageJ (도구) | **워크플로 기여** | 알고리즘은 표준. ccdproc 과 **비트동일** 증명이 자산이지 우월성이 아니다. 다만 **AutoPhOT·PhoPS·ASteCA 는 raw 를 다루지 않는다** — 파이프라인 경계로서는 실재하는 차이 |
| 내장 WCS 솔버 (인덱스 불요) | **PhoPS** 가 같은 문제를 같은 방향으로 푼다 (동적 Gaia 인덱스) | **주장 불가로 하향** | 2026-07 이전이라면 "인덱스 없는 측성"이 카드였다. 지금은 선행연구가 있으므로 **인용하고 차이를 서술**해야 한다 |
| 검출·구경측광 | SEP·photutils (프리미티브) · AutoPhOT·PP (도구) | **주장 불가** | 커널은 위임. APEX 기여는 오케스트레이션·QC |
| 프레임 QC 게이트 (depth-cost, 초과검출) | BANZAI 등 시설 QC. 물리 보정된 depth-cost 게이트는 흔치 않음 | **주장 가능(조건부)** | `APEX_NATIVE_PROVENANCE.md` 가 "논문 신규성"으로 표시. **다만 문헌 대조가 아직 없다** — 확인 필요 |
| 구경보정 워크플로 | photutils 프리미티브 + 각 파이프라인의 자체 정책 | **워크플로 기여** | 감사 문서 문구 그대로: "does not replace photutils" |
| PTC 검출기 특성화 | LSST `cp_pipe` PTC 태스크. 방법은 Janesick 표준 | **주장 불가 (방법)** / **주장 가능 (결과)** | 구현은 수십 줄. **가치는 발견에 있다 — 헤더 EGAIN 이 14배 틀렸고 측광에 미치는 영향을 정량화** |
| CMD 이소크론 피팅 | **ASteCA** — 같은 문제(나이·금속량·거리·소광 동시) | **주장 불가로 하향** | ASteCA 는 합성 CMD + 유전 알고리즘, APEX 는 EEP 보간 + 전색상 공분산 + Gaia 시차 사전분포 + MCMC. **접근이 다르나 문제는 같다.** 인용·비교 필수 |
| U/B/V 영점 재앵커 | — | **주장 가능 (결과)** | 이건 엔진이 아니라 **발견**이다. rail 해소 소거실험이 근거 |
| SYSREM | `PySysRem` | **주장 불가** | 알고리즘도 구현도 선행 존재. APEX 고유는 결측·가중·타깃 제외 계약 |
| PDM | **없음** (astropy 는 LS·BLS 만) | **워크플로 기여** | 선택한 스택에 부재하는 것은 사실. 다만 Stellingwerf 1978 구현이며 novel 아님 |
| 다중밤 앙상블·1일 별칭 해소 | AstroImageJ 앙상블 · wotan(겨냥 다름) | **미확인** | AIJ 가 다중밤 성단 규모를 다루는지 확인 필요 |

---

## 3. 남는 것 — 논문이 실제로 주장할 수 있는 것

개별 구성요소에 알고리즘 노벨티는 **거의 없다.** 그건 문제가 아니다 —
AutoPhOT 도 photutils·astropy 위에 서 있고 A&A Section 15 에 실렸다.
남는 것은 셋이다.

**(1) 파이프라인 경계** — 확인된 사실로 뒷받침된다.
AutoPhOT·PhoPS·ASteCA 는 **이미 보정된 데이터에서 시작**한다. PhotometryPipeline
과 AstroImageJ 는 보정을 다루지만 각각 소행성 측광·트랜싯 시계열에 특화돼
있다. **raw → science → (성단 CMD ∧ 다중밤 광도곡선)** 을 한 도구에서
끝내는 조합은 조사 범위에서 확인되지 않았다.

**(2) 발견** — 구현자가 누구든 무관하게 성립한다.
- 헤더 EGAIN 이 실측 대비 14배 틀림 (0.0495 vs 0.689) 과 그 측광 영향
- B 필터 faint 편차의 범인이 Gaia BP 감광이고, 밴드마다 범인이 바뀜
- 재현성의 근원이 SEP 디블렌딩 (2.7 %, 중앙 MAD 0.0 mmag, 꼬리 15.9 mmag)
- U/B/V 영점을 한 표준계에 앵커해야 이소크론 rail 이 풀림

**(3) 검증의 폭** — AutoPhOT 이 쓴 전략과 같고, APEX 가 더 넓다.
ccdproc 비트동일 · IRAF 교차검증 · PS1/Gaia 대조 · 문헌 주기 재현 ·
LCO 2기기 vs BANZAI 교차기기.

---

## 4. 반드시 인용해야 할 것 (누락 시 심사에서 지적됨)

| 문헌 | 왜 |
|---|---|
| Perren, Vázquez & Piatti 2015, A&A 576, A6 (ASteCA) | Step 12 와 같은 문제를 푸는 선행 도구 |
| Collins+ 2017, AJ 153, 77 (AstroImageJ) | 형태가 가장 가까운 GUI 데스크톱 도구 |
| Mommert 2017 (PhotometryPipeline) | "small to medium-sized observatories" 를 먼저 표방 |
| Brennan & Fraser 2022, A&A 667, A62 (AutoPhOT) | 같은 절의 직전 사례이자 문체 준거 |
| Erece & Kilic 2026, arXiv:2607.27414 (PhoPS) | 인덱스 없는 측성을 먼저 발표 |
| Tamuz+ 2005 · Stellingwerf 1978 · Hippke+ 2019 | 알고리즘 출처 (일부는 이미 인용) |

---

## 4-b. AstroImageJ 본문 확인 결과 (2026-08-09)

기능 목록 22 개를 본문에서 직접 확인했다. **AIJ 는 예상보다 넓다** — 5 번이
"Data Processor (DP) facility for image calibration including bias, dark, flat,
and nonlinearity correction" 이므로 **raw 보정을 한다.** 6 번이 시계열
다중구경 차등측광, 10 번이 동시 추세제거 적합이다. 즉 **광도곡선 쪽에서는
raw → science → LC 를 이미 한 도구가 덮는다.**

두 가지가 갈린다.

**(1) 측성이 네트워크 의존이다.** 기능 11 번은
*"Plate solving and addition of World Coordinate System (WCS) headers to images
seamlessly using the **astrometry.net web interface**"* 다. 웹 인터페이스이므로
인터넷 없는 돔에서는 성립하지 않는다. **APEX 내장 솔버의 오프라인 알리바이는
AIJ 에 대해 살아 있고**, git 상 `wcs_solve.py` 2026-04-28 · quad solver
06-10 으로 PhoPS arXiv 제출(07-29)보다 앞선다 — "독립 개발, 동시기 발표"로
서술 가능하다.

**(2) 성단 경로가 없다.** 22 개 기능 어디에도 마스터 카탈로그·다중필터
교차매칭·CMD·이소크론이 없다. 초록은 *"streamlined for time-series
differential photometry … especially exoplanet transits"*, 사용처는 KELT
트랜싯 후속관측 팀이다.

### 그래서 "파이프라인 경계" 주장을 이렇게 좁힌다

| 갈래 | 판정 |
|---|---|
| **광도곡선** | **novelty 주장 불가.** AstroImageJ 가 raw → LC 를 이미 덮는다. APEX 는 비교·위치잡기만 한다 (차이: Python vs Java, 헤드리스 CLI·재현성, 다중 밤 병합, 성단 규모 마스터 카탈로그 — 마지막 둘은 AIJ 에 부재) |
| **성단 CMD** | **주장 가능.** raw → science → 마스터 카탈로그 → CMD → 이소크론을 한 도구에서 끝내는 사례가 조사 범위에서 확인되지 않았다. ASteCA 는 카탈로그에서 시작하고, AIJ 는 CMD 경로가 없다 |
| **오프라인 측성** | **조건부 주장 가능.** AIJ 는 웹, ASTAP·astrometry.net 은 인덱스 DB 필요, PhoPS 는 동적 생성이나 2026-07 발표. 인용 후 "독립 개발" 서술 |

**함의**: 논문의 무게 중심을 **성단 CMD 갈래**에 두고, LC 갈래는 "같은 도구
안에서 같은 측광 산출물로 이어진다"는 통합 논거로만 쓰는 것이 안전하다.

또 하나 — AIJ 도 *"verified the accuracy of AIJ against IRAF, IDL, and MaxIm
DL"* 로 같은 검증 전략을 쓴다. APEX 의 ccdproc·IRAF·PS1/Gaia 대조는 이 절의
표준 관행이며, 그 자체로 novelty 는 아니지만 **빠지면 감점**이다.

## 5. 아직 확인 못 한 것

논문 문장으로 쓰기 전에 확인해야 한다. **기억으로 쓰지 말 것.**

1. ~~AstroImageJ 가 성단 규모 CMD 를 다루는가~~ → **4-b 에서 해소.**
   기능 목록에 CMD·마스터 카탈로그·이소크론이 없다. 다만 raw 보정과 LC 는
   덮으므로 주장 범위를 성단 갈래로 좁혔다.
2. **PhotometryPipeline 이 raw 보정(bias/dark/flat)을 하는가.** 지금 확인된
   것은 정합·측광·보정(zero-point)까지다.
3. **depth-cost 게이트·초과검출 게이트에 문헌 선례가 있는가.** 지금 "논문
   신규성" 표기의 유일한 근거는 내부 커밋 메시지다.
4. **ASteCA 와 APEX step12 의 실제 결과 비교.** 같은 성단(M67·NGC 6811)에
   둘 다 돌려 비교하면 §5 가 크게 강해진다 — ASteCA 입력이 카탈로그이므로
   **APEX Step 7 산출물을 그대로 먹일 수 있다.** 비용이 낮다.

---

## 6. ASteCA 직접 비교 — 실현 가능성 확인 (2026-08-09)

`APEX_PRIOR_ART.md` §5-4 의 "APEX Step 7 산출물을 ASteCA 에 그대로 먹인다"가
실제로 되는지 확인했다. **된다.**

| 항목 | 확인 결과 |
|---|---|
| 설치 | `pip install asteca` 로 0.7.0 설치됨 (Python ≥ 3.12). **APEX venv 를 건드리지 않도록 별도 venv** `C:\ast_v` 에 두었다 (경로가 길면 Windows 파일명 제한에 걸린다) |
| 입력 API | `Cluster(ra, dec, mag, e_mag, color, e_color, plx, e_plx, pmra, e_pmra, pmde, e_pmde, …)` — numpy 배열을 직접 받는다. 파일 포맷 변환이 필요 없다 |
| APEX 산출물 적합성 | `cmd_zeropoint/median_by_ID_filter_wide_cmd.csv` 가 **그 컬럼을 전부 갖고 있다.** M67 기준 1,005 별, `mag_cal_g/r/i` + 오차, `ra_deg`·`dec_deg`, 그리고 `parallax`·`pmra`·`pmdec` + 오차(ASteCA 멤버십 판정용) |
| 이소크론 격자 | ASteCA 는 PARSEC·MIST·BASTI 를 읽는다. **APEX 가 PARSEC CMD 3.9 격자를 이미 갖고 있다** — SDSS 396 MB(`umag gmag rmag imag zmag`), Johnson 434 MB(`UXmag BXmag Bmag Vmag Rmag Imag …`). 컬럼명(`Zini`·`logAge`·`Mini`)도 ASteCA 기대값과 같다. **다운로드 불필요이고, 두 도구가 동일한 이론 모델을 쓴다** |

**방법론 차이가 실재하며 논문에 쓸 만하다:**

- **ASteCA** — 합성 성단 생성(IMF·질량·이항성·소광법칙) 후 관측 CMD 와
  Poisson likelihood ratio(Tremmel+ 2013)로 비교
- **APEX** — 이소크론에 직접 MCMC. EEP 보간 + 전색상 공분산 + Gaia 시차
  사전분포

같은 성단·같은 측광·같은 이론 격자에 **접근이 다른 두 도구**를 태우는 비교는
§5 에 그대로 들어간다.

### 걸림돌 — APEX 쪽 step12 결과가 재처리 트리에 없다

재처리는 Step 0–7 + CMD10(영점)까지만 돌았다. `E:\APEX_validation\reprocess\
M67\result\` 에 `cmd_zeropoint` 는 있으나 이소크론 산출물이 없다. 메모리에
기록된 M67 복귀값([M/H] +0.06 · 나이 3.83 · E 0.011)은 `observed_Analysis`
트리의 것이고 그 입력은 이미 없다.

**따라서 비교 순서는**: (1) 재처리 트리에 APEX step12 실행 →
(2) 같은 `median_by_ID_filter_wide_cmd.csv` 로 ASteCA 실행 → (3) 대조.
(1) 이 MCMC 라 비용이 있다.
