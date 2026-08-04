# 레퍼런스 지식베이스

논문에 인용하는 문헌이 **실제로 무엇인지** 기록한다. 서지정보가 아니라, 그 문헌이
무엇을 했고 **무엇을 뒷받침할 수 있는지**를 적는다.

원문 PDF는 `refs/` 에 있다 (git 제외 — 저작권). 받은 곳은 각 항목에 URL 로 적었다.

## 왜 만들었나

2026-07-29, HOPS 를 서론에서 경쟁 도구로 계속 인용해 놓고 **그것이 1쪽짜리 학회 초록이라는
사실을 그날 처음 확인**했다. 사용자 지적: *"역시 래퍼런스 연구가 중요하다니깐"*.
서지만 맞으면 된다고 여긴 탓이다. **인용의 범위는 서지가 아니라 원문이 정한다.**

## 기록 규칙

- `원문 확인` = 원문을 읽었다 · `서지 확인` = 종류·저널만 확인 · `미확인` = 둘 다 안 함
- **미확인 문헌의 내용을 근거로 문장을 쓰지 않는다.** 서지 나열까지만 허용
- 「APEX 가 주장할 수 있는 것」에 없는 것을 본문에 쓰면 인용 범위 초과다

---

# 빠른 답

## 경쟁 도구 한눈에

| 도구 | 취지 | 문헌 종류 | 인터페이스 | 정확도 수치 | 범위 |
|---|---|---|---|---|---|
| **HOPS** | 교육 + 시민과학 + 연구 | **EPSC-DPS 초록 1쪽** | GUI | **없음** | 단일 시야 식(transit) |
| **AstroImageJ** | **교육 + 연구 (APEX 와 동일)** | AJ 정식 논문 27쪽 | GUI | **"확인했다" 한 문장, 수치 없음** | 단일 시야 시계열 |
| **Munipack** | 변광성 | `@MISC` 소프트웨어 등록 | GUI | 없음 | 시계열 |
| **VStar** | 변광성 | JAAVSO | GUI | 없음 | 시계열 분석 |
| PhotometryPipeline | 자동 측광 | Astronomy & Computing | 스크립트 | 있음 | 소행성·일반 |
| ASTROPOP | 편광·측광 | **PASP** | 스크립트 | 있음 | 일반 |
| prose | 모듈형 영상처리 | MNRAS | 스크립트 | 있음 | 일반 |
| **AutoPhOT** | 과도천체 자동화 | A&A | 스크립트 | **있음 (본문 12%)** | 과도천체 |
| **APEX** | 교육 + 연구 | (준비 중) | **GUI** | **있음 (본문 40%)** | **성단 CMD + 다중야간 LC** |

**한 줄 결론.** GUI 도구 넷 중 정식 저널 논문은 AstroImageJ 하나뿐이고, 제대로 논문화된
파이프라인 넷은 **전부 스크립트 기반**이다. APEX 의 자리는 *GUI 쪽에 수치가 있는 검증이
없다*는 공백이며, 차별점은 **대상 독자가 아니라 범위와 검증의 형태**다.

## 투고처 한눈에 (2026-07-29 확인)

| 저널 | IF | 색인 | 상시 여부 | APEX 적합도 |
|---|---:|---|---|---|
| **PASP** | **6.8** | **SCIE** | 상시 | 높음. **ASTROPOP 선례** |
| Astronomy and Computing | 1.8 | SCIE | 상시 | 높음 |
| RASTI | 2.5 | **ESCI (SCIE 아님)** | 상시 | 높음이나 색인이 걸림 |
| SPIE JATIS | — | SCIE | 상시 | 낮음 (기기 시스템 중심) |
| SPIE ATI 논문집 | — | **비 SCIE (학술대회)** | **격년, 다음 2028** | 높음. IRAF 원문이 여기 |
| A&A | ~6 | SCIE | 상시 | 높음이나 문턱 높음. AutoPhOT 게재지 |

- SCI 와 SCIE 는 사실상 통합되어 국내 서류에서 "SCI(E)" 로 함께 쓴다. **ESCI 는 별도 색인**이며
  2023년부터 IF 를 받지만 SCIE 는 아니다. 인정 여부는 기관마다 다르므로 확인이 필요하다.
- SPIE ATI 2026 은 7월 5–10일 코펜하겐에서 열렸고 **이미 끝났다**. 소프트웨어 트랙은
  "Software and Cyberinfrastructure for Astronomy IX".

---

# 문헌별 상세

## 경쟁 도구

### collins2017 · AstroImageJ

| | |
|---|---|
| 무엇 | ImageJ 위에 얹은 천문 전용 GUI. 영상 보정부터 광도곡선 적합까지 |
| 문헌 | The Astronomical Journal (SCIE). 27쪽 확장판 + 동료심사 축약판(2017b) 별도 |
| 원문 | `refs/aij.pdf` · https://arxiv.org/pdf/1701.04817 |
| 확인 상태 | **원문 확인** (2026-07-29) |
| 목차 | 1 서론 / 2 개요와 기본 기능(툴바·영상 표시·측정 표) / 3 Data Processor 영상 보정 / 4 초정밀 측광과 광도곡선(단일 구경·다중 구경 차등·Multi-Plot·적합과 추세제거·비교별 앙상블) / 5 업데이터 — **검증 절 없음** |
| 서론 구조 | 문제 제기 → AIJ 소개 → **기능 목록 20개** → **사용 실적** → 정확도 한 문장 → 정밀도 → 대화형 기능 → 실시간 모드 → 문서 안내 |
| 자기 포지션(초록) | *"research grade image calibration and analysis tools with a GUI driven approach"* · *"enables new users, even at the level of **undergraduate student, high school student, or amateur astronomer**"* |
| 정확도 서술 | 본문 **한 문장**. *"We and the KELT follow-up team have verified the accuracy of AIJ against … IRAF, IDL, and MaxIm DL."* 수치·그림·절 없음. 27쪽에 compare/validate/benchmark 0회, agreement 1회 |
| 정밀도 근거 | 별도 논문(Collins et al. 2017a). WASP-12b·Qatar-1b 식 잔차 RMS 183·255 ppm, 타이밍 잔차 30초 이내. 이는 **정밀도(산포)**이지 기준 대비 **정확도(일치)**가 아니다 |
| 사용 실적 | KELT 후속관측 팀 약 30명 중 대부분, 행성 10편 발표 + 8편 진행. 학부 실습실, 고교 수업 |
| 주장 가능 | ① 단일 시야 시계열 차등 측광 중심 ② **APEX 와 같은 대상 독자를 명시한 선례** ③ 정확도 비교의 수치와 조건이 논문에 제시되지 않았다 |
| 주장 불가 | "AIJ 는 검증하지 않았다" — 했다고 밝히고 있다 |
| **함정** | **APEX 의 대상 독자 서술은 신규성이 아니다.** AIJ 가 2017년에 같은 말을 했다 |
| 통찰 | AIJ 의 신뢰 근거는 **"많이 쓰인다"**(KELT·행성 10편)다. APEX 에는 그 근거가 없으므로 **측정으로 대신**해야 한다. 검증 40% 의 정당화가 여기 있다 |

### tsiaras2019 · HOPS

| | |
|---|---|
| 무엇 | HOlomon Photometric Software. python GUI, 오픈소스, 3-OS |
| 문헌 | **EPSC-DPS 2019 학회 초록. 본문 1쪽이 전부** |
| 원문 | `refs/hops2019.pdf` · https://meetingorganizer.copernicus.org/EPSC-DPS2019/EPSC-DPS2019-1594-1.pdf |
| 확인 상태 | **원문 확인** (2026-07-29, 전문) |
| 기능(원문 명시) | (a) reduction — **마스터 bias/dark/flat 계산과 과학 프레임 보정** (b) 프레임 선별 (c) 정렬(자오선 반전 대응) (d) 측광 — 구경과 PSF, 대상·비교별 대화형 선택 (e) 식 맞추기 — MCMC. PyLightcurve 사용 |
| 자기 포지션 | 시민과학 겨냥. *"과학적 데이터 분석에 기여할 수 있고 동시에 교육 도구로도"* |
| 주장 가능 | 위 기능 목록. 교육과 연구를 함께 겨냥한 **선례** |
| 주장 불가 | **정확도·성능·수행시간·사용성.** 초록에 측정값이 하나도 없다 |
| **함정** | "논문"이라 부르면 안 된다. **HOPS 도 검출기 보정을 한다** — raw→science 를 배타적 차별점으로 쓰면 틀린다 |

### hroch2014 · Munipack / benn2012 · VStar

| | |
|---|---|
| 문헌 | Munipack `@MISC`(소프트웨어 등록, 저널 논문 아님) · VStar JAAVSO |
| 확인 상태 | **서지 확인** |
| 주장 가능 | 변광성 시계열 분석에 초점을 둔다는 것까지 |

### mommert2017 · PhotometryPipeline / campagnolo2019 · ASTROPOP / garcia2022 · prose

| | |
|---|---|
| 문헌 | Astronomy and Computing · **PASP** · MNRAS. 전부 SCIE 정식 논문 |
| 인터페이스 | **스크립트/명령행** |
| 확인 상태 | **서지 확인** |
| 주장 가능 | 자동화되어 있으나 스크립트로 구동한다는 것 |
| 부수 사실 | **ASTROPOP 이 PASP 에 실렸다** — PASP 가 이 종류를 받는다는 직접 증거 |

### brennan2022 · AutoPhOT

| | |
|---|---|
| 무엇 | 과도천체(transient) 자동 측광 파이프라인 |
| 문헌 | A&A 667, A62 (SCIE, open access CC BY). 16쪽 |
| 원문 | `refs/autophot2022.pdf` |
| 확인 상태 | **원문 확인** (2026-07-29, 서론·목차) |
| 목차 | 1 서론 / 2 전처리(영상 보정·스태킹·대상 식별·메타데이터·WCS·우주선·FWHM) / 3 측광(구경·PSF) / 4 보정(영점·색항·대기소광) / 5 영상 차감 / 6 한계등급 / **7 검증(2쪽 = 12%)** / 8 결론 + **부록 4개** |
| 서론 구조 | **IRAF 계보와 그 붕괴**(NOIRLab 개발 중단·64비트 불가·PyRAF 지원 종료) → 다른 패키지(SExtractor·A-PHOT·PP) → 자기 소개 → 대상 문제(이질적 자료, 관측자마다 파라미터가 달라 어긋남) → 가용성과 논문 구성 |
| 서론이 닫는 방식 | *"이 논문의 목적은 AutoPhOT 패키지를 **간략히 소개**하는 것"* |
| 주장 가능 | 「스스로 만들지 않은 기준에 대해 검증한다」는 방식이 이 논문과 같다는 것 |
| 대조점 | APEX 서론은 **"검증"** 으로 닫는다. 검증 비중 12% 대 40% |

### tody1986 · tody1993 · IRAF

| | |
|---|---|
| 문헌 | **Proc. SPIE 627**(20쪽) · ASPC 52(11쪽) |
| 원문 | `refs/iraf1986.pdf` · `refs/iraf1993.pdf` · https://iraf-community.github.io/doc/iraf.pdf · `/iraf92.pdf` |
| 확인 상태 | **원문 확인** (2026-07-29, 목차·서론) |
| 1986 목차 | 1 서론 / 2 시스템 구조 / 3 CL 명령어 / 4 응용 소프트웨어 / 5 프로그래밍 환경 / 6 가상 OS / 7 호스트 인터페이스 |
| 1986 서론 | 4문단. 프로젝트 연혁(1981 시작) → 채택 현황(STScI, 40개 기관) → 문서 범위 → 목적과 구성. **천문 내용 0** |
| **함정** | 경쟁 도구가 없던 시절의 글이라 **서론 모델로 쓰면 안 된다** |
| 부수 사실 | IRAF 원문이 **Proc. SPIE** 다. 이 바닥 소프트웨어 논문의 원조가 SPIE 논문집 |

---

## 방법·기준 문헌

### stetson1987 · DAOPHOT

| | |
|---|---|
| 문헌 | PASP 99, 191 (32쪽) · `refs/stetson1987.pdf` · https://articles.adsabs.harvard.edu/pdf/1987PASP...99..191S |
| 확인 상태 | **원문 확인** (2026-07-29, 초록·구성) |
| 무엇 | 밀집장 항성 측광 프로그램. 겹친 별상(blended stellar images)의 정확한 측광을 위한 수학적 알고리즘에 초점 |
| 특기 | 초록이 스스로 *"known shortcomings of the current program"* 을 논한다고 밝힌다. 도구 논문이 자기 한계를 쓰는 선례 |
| 주장 가능 | 밀집장 측광의 표준 참조. §3.11 교차확인의 기준 엔진 |

### bertin1996 · SExtractor

| | |
|---|---|
| 문헌 | A&AS 117, 393 (12쪽) · `refs/bertin1996.pdf` |
| 확인 상태 | **원문 확인** (2026-07-29, 초록) |
| 무엇 | 천체 검출·디블렌딩·측정·분류. 신경망 기반 별/은하 분리 |
| **함정** | 초록이 *"particularly suited to the analysis of large **extragalactic** surveys"* 라고 밝힌다. **밀집 성단 측광용이 아니다** — 그 자리는 DAOPHOT 이다. 검출 엔진으로만 인용해야 한다 |
| 주장 불가 | 밀집장 항성 측광 성능. SEP(barbary2016)와 혼동 금지 — **SEP 는 Barbary 2016 이다** |

### stetson1994 · ALLFRAME

| | |
|---|---|
| 문헌 | PASP 106, 250 (31쪽) · `refs/stetson1994.pdf` |
| 확인 상태 | **원문 확인** (2026-07-29, 초록) |
| 무엇 | 구상성단 M15 코어(중심 2분각)의 CMD·색-색도를 CFHT + HST 자료로. **ALLFRAME 을 처음 상세히 기술** |
| 방법 | 모든 프레임의 기하·측광 정보를 **동시에** 사용 |
| 주장 가능 | 밀집 성단 코어 측광의 선례. §3.12(M5·M13 코어)의 문헌 근거 |

### mccully2018 · BANZAI

| | |
|---|---|
| 문헌 | Proc. SPIE (9쪽) · `refs/mccully2018.pdf` · arXiv:1811.04163 |
| 확인 상태 | **원문 확인** (2026-07-29, 초록) |
| 무엇 | LCOGT 로봇 망원경망의 실시간 처리 파이프라인. 하룻밤 수천 장 |
| 수행 범위(원문) | 기기 특성 제거(불량화소 마스킹, bias·dark 제거, flat 보정), astrometric fitting, 소스 카탈로그 추출 |
| 주장 가능 | §3.5 교차기기 비교의 기준. **독립적으로 개발·운영되는 출판 파이프라인**이라는 것 |

### riello2021 · Gaia EDR3 측광

| | |
|---|---|
| 문헌 | A&A 649, A3 (35쪽) · `refs/riello2021.pdf` · arXiv:2012.01916 |
| 확인 상태 | **원문 확인** (2026-07-29, BP 어두운 쪽 서술 검색) |
| **핵심 문장** | *"**BP tends to be systematically brighter towards the faint end**: it would therefore make sense to include a restriction on G_BP in the archive query"* |
| 기전 | 어두운 쪽에서 BP/RP 스펙트럼의 S/N 이 낮고 **잔여 배경(residual background)** 이 남는다. 보정 계수가 Table 5 에 두 등급 구간으로 주어진다 |
| 주장 가능 | §3.13 의 "Gaia BP 참조의 알려진 결함" 이 **문헌으로 확인된 사실**이라는 것. 원문이 아카이브 질의에서 G_BP 제한을 권고할 정도다 |
| 주의 | 부호를 쓸 때 확인할 것 — 원문은 BP 가 어두운 쪽에서 **더 밝게**(등급이 작게) 치우친다고 말한다 |

### janesick2007 · Photon Transfer

| | |
|---|---|
| 문헌 | **단행본** (SPIE Press, *Photon Transfer: DN → λ*) |
| 확인 상태 | **미확인 — 책이라 내려받을 수 없다** |
| 주장 가능 | PTC 방법의 표준 참조로 서지 인용까지. **책의 특정 내용을 근거로 문장을 쓰지 않는다** |

---

## 확인이 남은 것

서지만 확인했거나 미확인이다. **내용을 근거로 문장을 쓰기 전에 원문을 봐야 한다.**

`stetson1987`(DAOPHOT) · `bertin1996`(SExtractor) · `barbary2016`(SEP) · `photutils` ·
`astropy2013/2018/2022` · `stetson1990`(성장 곡선) · `stetson1994`(ALLFRAME 밀집장) ·
`janesick2007`(PTC) · `vandokkum2001`(L.A.Cosmic) · `mccully2018`(BANZAI) ·
`tonry2018`(ATLAS) · `kessler2015`(DES) · `smtn002`(Rubin 한계등급) · `tamuz2005`(SYSREM) ·
`stellingwerf1978`(PDM) · `lang2010`(astrometry.net) · `chambers2016`(PS1) ·
`riello2021`(Gaia BP) · `pancino2022`(표준별) · `yang2018`(YZ Boo) · `astap` · `ccdproc` 외

**우선순위** — 본문에서 그 내용을 근거로 삼는 순서대로: `stetson1987` · `bertin1996` ·
`janesick2007` · `stetson1994` · `mccully2018` · `riello2021`.

### 2026-08-04 · Appendix B instrument specifications

The following manufacturer/observatory pages are recorded as specification
sources, separate from the measured detector constants in Fig. 3:

| BibTeX key | source | status | use in manuscript |
|---|---|---|---|
| `qhy600_specs` | QHYCCD, QHY600PH Series, https://www.qhyccd.com/astronomical-camera-qhy600/ | source page checked | Appendix B: IMX455, 3.76 μm, 9576×6388, 16-bit |
| `lco_sinistro_specs` | Las Cumbres Observatory, Sinistro instrument, https://lco.global/observatory/instruments/sinistro/ | source page checked | Appendix B: CCD486, 4096×4097, 15 μm, four amplifiers |
| `moravian_c361000_specs` | Moravian Instruments, C3-61000, https://www.gxccd.com/shop?action=product&cat=28&lang=405&page=1348&subcat=0 | source page checked | Appendix B: IMX455, 9576×6388, 3.76 μm |

These entries support hardware specifications only. Gain, read noise,
dark current, and the 2×2 pixel scale used in the main validation remain the
measured or configured values reported in the manuscript, not values inferred
from the product pages.
