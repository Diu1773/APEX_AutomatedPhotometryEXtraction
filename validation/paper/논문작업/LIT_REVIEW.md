# APEX 논문 — Literature Search Report (ARS lit-review, 2026-07-12)

Phase 0 설정 = 사용자 브리프로 확정(포지셔닝·타겟·검증원칙 변경 금지). Phase 1 산출물.
**모든 항목은 WebSearch/ADS/DOI로 실존 확인함. 미확인 인용 없음.**

## Paper Configuration Record (요약)

| 항목 | 값 |
|---|---|
| 유형 | Tools/software paper (IMRaD 변형: Design–Validation–Demonstration) |
| 분야 | 천문 기기·방법 (astro-ph.IM) |
| 투고 | RASTI 1지망 (스코프: software·data analysis 명시, **AI 사용 공개 요구 — 우리 thesis와 정합**) / PASP 2지망 |
| 언어 | 영어 (국문본은 검토용) |
| 타겟 독자 | 도구 장벽에 막힌 학부생·입문 연구자·소형망원경 관측자 |
| 기여 | 접근성(가이드 raw→science GUI+CLI) + 다층 검증 |

## Search Strategy

- 경로: ADS·arXiv·저널원문(A&A/MNRAS/PASP/JOSS)·WebSearch 교차확인. 세션 내 원문 열람: AutoPhOT(ar5iv 전문), AstroImageJ(전문), DES Balrog(HTML 전문+Fig3), AutoPhOT Fig11/12/15/16 이미지.
- 축: A 비교도구 / B 기반원전 / C 검증방법론·그림관례 / D 접근성·정책.
- 스크리닝: 관련성(우리 §에 실사용) + 실존검증 통과만 포함. 초기후보 ~45 → 최종 40.

## Screening Results
- 기존 bib 28키 전부 유지(전부 실사용 중) · 신규 검증 통과 12건 추가 → 최종 40.

---

## Annotated Bibliography

### 축 A — 비교 도구 (§1.2 갭 논증, §5.4 위치짓기)

**Brennan & Fraser (2022)** A&A 667, A62. doi:10.1051/0004-6361/202243067 — AutoPhOT. **주 준거.**
- 검증방식: §7에서 SExtractor/DAOPHOT/PSFEx와 같은 별 Δmag·Δerr vs 기기등급 패널(Fig14), 구경↔PSF 자체비교(Fig13), 한계등급=인공별 주입 후 **읽어내기**(Fig12, 곡선피팅 없음), erf 검출확률은 부록D 식으로만.
- 채택: 검증 논문 구조(방법 파이프라인순+집중 검증절), Fig12형 완전도, Fig13형 구경↔PSF, Fig14형 독립엔진 Δ-vs-mag. 버림: transient 특화 부분(template subtraction).

**Collins et al. (2017)** AJ 153, 77. doi:10.3847/1538-3881/153/2/77 — AstroImageJ.
- 검증방식: 정식 검증절 없음(~90% 기능기술). "IRAF·IDL·MaxIm DL과 대조했다" 산문 1문단 + 정밀도 사례 인용. 비교표 없음.
- 채택: GUI 도구 논문의 선례(유저스터디 불요, 산문검증 허용선). 우리는 이보다 검증을 강하게 → 차별점 근거.

**Tsiaras et al. (2019)** — HOPS (기존 bib tsiaras2019). §1 훅·갭. 단일시야 transit 특화 GUI 선례.

**Mommert (2017)** Astronomy & Computing 18, 47–53. doi:10.1016/j.ascom.2017.02.002, arXiv:1702.00834 — PHOTOMETRYPIPELINE.
- 검증방식: 보정등급 정확도 ≤0.03 mag, 천체측정 ~0.3″ vs 카탈로그 — **파이프라인 정확도를 카탈로그 대조 수치로 요약하는 관례**.
- 채택: §1.2(CLI 자동 파이프라인, GUI 아님·소행성 특화) + §3.7 천체측정 정확도 보고 형식(카탈로그 대비 ″ 단위).

**Campagnolo (2019)** PASP 131, 024501. arXiv:1811.01408 — ASTROPOP. §1.2: 파이썬 측광+편광 파이프라인, 특정 기기 지향, GUI 아님.

**Garcia et al. (2022)** MNRAS 509, 4817. doi:10.1093/mnras/stab3113 — prose. §1.2: 모듈러 파이프라인 *프레임워크* — 스크립팅 전제 → "파이썬 문턱" 논증에 정확히 부합.

**Karpov, STDPipe** ASCL:2112.006 (+ STDWeb arXiv:2411.16470). §1.2: transient 특화 루틴 모음, 역시 라이브러리형. (저널논문 없음 — ASCL 인용)

**Hroch (2014)** ASCL — Munipack (기존 hroch2014) / **Benn (2012)** JAAVSO — VStar (기존 benn2012). §1.2 변광성 GUI들.

→ **갭 확정**: GUI형(AIJ·HOPS·MuniWin·VStar)=단일시야/시계열 특화, 파이프라인형(PP·ASTROPOP·prose·STDPipe·AutoPhOT)=스크립팅/CLI 전제. **raw→성단CMD·다중밤LC를 가이드 GUI로 잇는 도구 부재** — §1.2 산문 논증 완성 (표 불요, AIJ·AutoPhOT 선례).

### 축 B — 기반 원전 (§2 인용 의무)

**Stetson (1987)** PASP 99, 191 — DAOPHOT (기존). 인공별 시험 원류 + IRAF 교차검증 대상.
**Bertin & Arnouts (1996)** A&AS 117, 393 — SExtractor (기존). SEP의 알고리즘 원전.
**Barbary (2016)** JOSS 1(6), 58. doi:10.21105/joss.00058 — **SEP. 신규 필수** (우리 주 검출엔진인데 bib 누락이었음!).
**photutils / astropy 3부작** (기존). **Lang et al. (2010)** AJ 139, 1782 — astrometry.net (기존). **ASTAP** (기존 astap).
**Craig et al., ccdproc** ASCL:1510.007, doi:10.5281/zenodo.1069648 — **신규 필수** (§3.4 비트동일 대조인데 bib 누락!).
**McCully et al. (2018)** Proc. SPIE 10707, arXiv:1811.04163 — BANZAI. **신규 필수** (§3.5 교차기기 대조인데 bib 누락!).
**van Dokkum (2001)** — L.A.Cosmic (기존). **Tamuz et al. (2005)** — SYSREM (기존). **Stellingwerf (1978)** — PDM (기존).
**VanderPlas (2018)** ApJS 236, 16. doi:10.3847/1538-4365/aab766 — Lomb–Scargle 이해. **신규** (§2.4 LS 인용 보강).
**Tody (1986)** Proc. SPIE 627, 733 / **Tody (1993)** ASP Conf 52, 173 — IRAF 원전. **신규** (§1에서 IRAF를 논하며 원전 미인용은 결례).

### 축 C — 검증 방법론·그림 관례 (§3 사양의 근거)

**Anbajagane et al. (2025)** OJA, arXiv:2501.05683 — DES Y6 Balrog (기존 anbajagane2025).
- 1.46억 주입 → 구간비율 자체가 매끈, **함수피팅 없음**(Fig3: 검출률 vs 등급, 밴드별 + 90% 세로선). 크라우딩은 2D맵/별도 런.
- 채택: "표본 크면 곡선=데이터" 원칙, 깊이 세로선 표기.

**Masci (2011)** IPAC 노트 / **Kashyap et al. (2010)** ApJ 719, 900 (기존) — erf 검출확률=S/N 함수. §3.6 본문 인용 전용(그림 아님).

**Stetson (1990)** PASP 102, 932. doi:10.1086/132719 — DAOGROW 성장곡선. **신규 필수** — 우리 자작 apcorr의 방법론 원전. §2.3 apcorr 서술 + §3 "성장곡선 기반 구경보정은 Stetson 1990 확립 기법, 우리는 그 자동화 워크플로" 논증.

**Janesick (2007)** *Photon Transfer: DN→λ*, SPIE Press PM170. **신규 필수** — §3.3 PTC gain·RN 실측의 표준 방법론 원전 (현재 alarcon2023 검출기 대조만 있고 방법 원전 없음).

**Riello et al. (2021)** (기존) — Gaia BP faint 결함, §3.13 정직사례 근거. **Pancino et al. (2022)** (기존) — 표준화 맥락. **Chambers et al. (2016)** (기존) — PS1.

→ 그 밖의 관례 확인: 독립엔진 교차검증 = AutoPhOT Fig14형(Δmag vs mag, 잔차 평평성) — **우리 fig 이미 부합**. 구경↔PSF = AutoPhOT Fig13형 — **부합**(우리는 +이웃거리 축 확장). 보정 대조 = 차영상+픽셀히스토그램(BANZAI/ccdproc 실무 관례) — **부합**.

### 축 D — 접근성·정책 (§1 훅·§6·disclosure)

**Fitzgerald et al. (2014)** PASA 31, e037. arXiv:1407.6586 — 20년간 고교급 천문 연구 프로젝트 리뷰. **신규** — "자료 획득은 쉬워졌으나 **자료 처리·분석이 병목**"이라는 우리 §1 논증의 문헌 근거 (ARiC 프로젝트들의 성패 요인 분석).
**RASTI 저널 정책** (academic.oup.com/rasti/pages/general-instructions) — AI 저자 불인정 + **사용 공개 요구**. §back matter AI disclosure가 저널 요구와 정합함을 명시할 근거 (본문 인용은 editorial: RASTI 1,1 doi:10.1093/rasti/rzac002 — 저널 창간사·스코프).

---

## Literature Matrix (요약)

| 소스 | 포지셔닝(A) | 원전(B) | 검증관례(C) | 접근성(D) | 품질 |
|---|---|---|---|---|---|
| Brennan&Fraser 2022 | **main** | | **main** | x | High |
| Collins 2017 | main | | x | x | High |
| Mommert 2017 / Campagnolo 2019 / Garcia 2022 / STDPipe | main | | | | Med-High |
| Stetson 1987/1990 | | **main** | main | | High |
| Barbary 2016·photutils·astropy·Lang·ASTAP·ccdproc·BANZAI·vanDokkum·Tamuz·Stellingwerf·VanderPlas | | main | x | | High |
| Anbajagane 2025 · Masci · Kashyap | | | main | | High |
| Janesick 2007 | | main | main | | High |
| Riello·Pancino·Chambers | | x | main | | High |
| Tody 1986/1993 · Fitzgerald 2014 · RASTI editorial | x | x | | **main** | High |

## Identified Gaps (→ 우리 기여문)
1. **raw→성단 CMD/다중밤 LC 전 구간 가이드 GUI 부재** (GUI형=특화, 파이프라인형=스크립팅).
2. **GUI 측광 도구의 체계적 검증 부재** (AIJ 산문 1문단이 현 수준) — 우리 §3가 갭 자체를 메움.
3. **AI 보조 개발 과학 SW의 검증 프레임 부재** — RASTI 정책은 '공개'까지만; 우리는 공개+검증 결합.
4. 교육 문헌(Fitzgerald 2014)이 지적한 **분석 병목**을 직접 겨냥한 도구 논문 부재.

---

## ★ 확정표: 컴포넌트 × 채택 검증방식 × 근거 레퍼런스 × 그림 형식 (fig 재작성 사양서)

| § | 컴포넌트 | 채택 검증방식 | 근거 레퍼런스 | 그림 형식 (관례) | 현재 fig 판정 |
|---|---|---|---|---|---|
| 3.2 | 검출기 보정 | 합성 inject-recover + 실프레임 전후 | 표준 축소(ccdproc 문서 관례) | 전/후 영상+프로파일, 잔차 통계 | **유지** (fig10→새1) |
| 3.3 | 검출기 특성화 | PTC (flat-pair 분산·평균) | **Janesick 2007** | 분산 vs 신호 log-log + gain 적합 | **유지**, Janesick 인용 추가 (fig11→새2) |
| 3.4 | ccdproc 교차 | 픽셀 대 픽셀 차이 | **Craig ccdproc** | 단계별 Δ 히스토그램/막대 | **유지**, ccdproc 정식 인용 추가 (fig12→새3) |
| 3.5 | 교차기기 | 독립 파이프라인 산출 대조 | **McCully 2018 (BANZAI)** | 차영상 + 픽셀분포 | **유지**, BANZAI 정식 인용 추가 (fig13→새4) |
| 3.6 | 검출 완전도 | 대량 MC 주입→구간비율, 깊이 read-off, 곡선피팅 없음 + 컷아웃 | Stetson 1987 · **AutoPhOT Fig12** · **Balrog Fig3** · Masci/Kashyap(erf, 본문) | 완전도 vs 등급(점구름+pooled선+m50선) + 주입 컷아웃 줄 | **이미 재작성 완료** (21k, fig1→새5) |
| 3.7 | 천체측정 해 | 별-Gaia 잔차 통계 + 내장솔버 vs ASTAP/astnet 해 일치 | Lang 2010 · **Mommert 2017**(″단위 보고) | 잔차 산포/히스토그램(px·″) + 솔버간 Δ | **신규 제작 필요** (F_wcs) |
| 3.8 | 오차 모델 | 주입 참값 pull 분포 + σ=1.0857/SNR 추종 | 통계 표준(+AutoPhOT Fig14의 Δerr 관례) | pull 히스토그램+단위정규, RMS vs SNR | **유지** (fig2→새6) |
| 3.9 | 민감도 | 파라미터·조건 스윕 | (관례 자유 — DES 계열 robustness) | 지표 vs 파라미터 다패널 | **유지** (fig3→새7) |
| 3.10–3.11 | 독립엔진 교차 | 같은 별·같은 파라미터 Δmag vs mag, MAD/RMS | **AutoPhOT Fig14** · Stetson 1987 | Δmag vs mag 산포+잔차 평평성 (T2 파라미터 매칭표) | **유지** (fig4·5→새8·9), 형식 이미 부합 |
| 3.12 | 밀집장 구경↔PSF | 내부 일치 vs 이웃거리 | **AutoPhOT Fig13**(구경↔PSF 관례) | Δ(구경−PSF) vs 분리, 구간 중앙값 | **유지** (fig9→새10) |
| 3.13 | 참조카탈로그 | PS1 독립 대조로 참조결함 분리 | Riello 2021 · Chambers 2016 | Δ vs mag 밴드별 (정직사례) | **유지** (fig7→새11) |
| 3.14 | CMD 재현 | 능선 대조 (피팅 없음) | Pancino 2022(표준화 맥락) | CMD 오버레이+능선 잔차 | **유지** (fig8→새12) |
| 3.15 | 시계열 코어 | 주입 주기신호 복원 + PDM↔LS 일치 + 전/후 rms | **Tamuz 2005**(rms 개선 관례) · Stellingwerf 1978 · VanderPlas 2018 | 전/후 detrend LC + 주기도 복원 | **신규 제작 필요** (F_ts, YZ Boo 재실행 후) |
| 3.16 | 프레임 QC | 주입결함 판정 + blind spot 정직보고 | (관례 자유) | 판정 매트릭스 | **유지** (fig6→새13) |

**판정 요약**: "다 엉터리"는 아니다 — **문제였던 건 fig1(완전도)뿐이었고 이미 레퍼런스 준거로 재작성 완료.** 나머지 11개는 관례에 부합(4개는 인용만 보강). **없는 그림 2개**(천체측정 F_wcs, 시계열 F_ts)가 진짜 구멍.

## references.bib 갱신 목록

**추가 (12)**: barbary2016(SEP), ccdproc, mccully2018(BANZAI), stetson1990(DAOGROW), janesick2007(PTC), vanderplas2018(LS), tody1986, tody1993(IRAF), mommert2017(PP), campagnolo2019(ASTROPOP), garcia2022(prose), fitzgerald2014(교육). — karpov STDPipe(ASCL)·RASTI editorial(rasti2022)은 본문 필요시.
**유지 (28)**: 기존 전부.

## Recommended Sources by Section

| § | 핵심 소스 |
|---|---|
| §1 훅·갭 | Tsiaras, Tody×2, Stetson 1987, Bertin, photutils/astropy, Collins, Hroch, Benn, Mommert, Campagnolo, Garcia, Fitzgerald 2014 |
| §2 설계 | photutils, astropy×3, Barbary, Lang, ASTAP, ccdproc, vanDokkum, Stetson 1990, Tamuz, Stellingwerf, VanderPlas |
| §3 검증 | 확정표 열 참조 |
| §5 논의 | Brennan&Fraser, Collins, Anbajagane, Riello |
| Back matter | RASTI 정책, CRediT |
