# APEX Paper — Master Plan v2 (2026-07-10)

*v1(2026-07-05 ars-plan)을 대체. v1의 뼈대(RASTI, techniques-paper 구조, 정직성
스파인, ars 스킬시퀀스)는 계승하고, 랩미팅 피드백 + 포지셔닝 재논의(2026-07-10)
+ Step 0 통합 + all-APEX 재처리를 반영해 전면 개정. 논문 작업의 단일 앵커.*

---

## Locked decisions (v2)

| 항목 | 값 | 근거 |
|---|---|---|
| **Venue** | **RASTI** (1지망) / PASP (2지망) | 소프트웨어·기법 검증 전문지, AUTOPHOT 인접 |
| **포지셔닝** | photutils/astropy 위 **가이드형 GUI+CLI 통합 워크플로** — "HOPS가 transit에 해준 것을 성단 CMD·변광성에" | photutils와 정면대결 ❌, 그 위의 통합/가이드 레이어 ✅ |
| **"IRAF 대체"의 의미** | IRAF의 *알고리즘*이 아니라 **진입 경험**(설치 일주일·낡은 매뉴얼·막히면 끝)을 대체 | 알고리즘 경쟁은 지는 프레임 (photutils가 이미 공식 후계) |
| **타겟층** | 도구 장벽에 막힌 학생·입문연구자·자기주도 학습자·소규모 관측소 | §1 motivation으로 **서술만** — tools paper이므로 만족도/사용성 조사 안 함 (AstroImageJ·AutoPhOT 선례) |
| **주장 수준** | **Level 2**: "검증된 연구급 CMD·광도곡선·주기를 짜깁기·코딩 없이" | Level 1(교육)은 목표 미달, Level 3(ALLSTAR 정면)은 과욕 |
| **측광 스코프** | **B안**: 구경+apcorr = 주력(full 검증) / PSF(photutils ePSF) = 제공+일치검증 | crowded 코어는 PSF 제공하되 "가이드로 된다" 수준, mmag 경쟁 안 함 |
| **검증 철학** | 검증 깊이 = novelty만큼. **라이브러리 래퍼는 인용+end-to-end 산출물 검증, 자작 수치코어 4개만 직접 검증** | novelty는 통합/워크플로이지 알고리즘이 아님 |
| **AIPPI** | **전면 배제** — 저자 자작 툴은 독립 레퍼런스 아님. 그림·본문·배치 어디에도 검증 근거 사용 금지 | 순환논리 (기존 Fig10 패널c 제거 필요) |
| **이소크론(Step 12)** | **명시적 scope-out** — 도구는 제공, 정확도 주장 제외 | 퇴화 역문제 (age-[M/H]-거리-적색화) |
| **데이터** | **all-APEX 재처리본만** (Step0→1-7→10, E:\APEX_validation\reprocess\) — AIPPI 전처리 기반 결과 사용 금지 | raw→science 주장의 데이터 무결성 |
| **AI 개발 공개** | balanced frame: §1 한 문장 + disclosure 전문 | v1 계승 |

## Thesis (한 문장)

> 표준 알고리즘(photutils/astropy)을 가이드형 GUI+헤드리스 워크플로로 통합한
> 측광 파이프라인은, 비전문가를 raw 프레임에서 성단 CMD·변광성 주기까지
> 스크립팅 없이 데려가면서도, 다층·독립·재현가능한 검증을 통과해 연구급
> 정확도에 도달할 수 있다.

## 경쟁 지형 (§1·§5 도구 비교표의 원천)

| 툴 | 정체 | 빈틈 (APEX가 채우는 것) |
|---|---|---|
| IRAF/DAOPHOT | 레거시 툴킷 (EOL ~2019) | 설치·학습 절벽, 문서 부패 |
| photutils | 라이브러리 (IRAF 공식 후계) | 코드 작성 필수, 워크플로/GUI/보정 없음 |
| AutoPhOT | CLI 파이프라인 | 초신성/transient 특화 |
| AstroImageJ | GUI | transit 특화, 보정파이프라인·CMD 없음 |
| HOPS | GUI | transit 특화 |
| MuniWin/VStar | GUI | 좁은 범위·높은 학습장벽 (저자 실경험) |
| **APEX** | **가이드 GUI+CLI, raw→science** | **성단 CMD + 다중야간 변광성 주기, 통합·재현가능** |

---

## 파이프라인 전 스텝 맵 × 검증 매트릭스 (Table 1의 원천)

**자작 수치코어 4개 (직접 검증 대상):**
① 내부 quad 플레이트솔버 (Step5; Lang+2010식 quad 해시+RANSAC+`fit_wcs_from_points`)
② 성장곡선 apcorr (Step7; apcorr=1/f_enc, 프레임별 밝은 고립별)
③ 벡터화 Stellingwerf PDM (LC11)
④ SYSREM (LC10; Tamuz+2005)
나머지 과학코어는 전부 신뢰 라이브러리 얇은 래퍼.

| Step | 계산 | 감싼 것/자작 | 검증 | 상태 |
|---|---|---|---|---|
| 0 보정 | bias/dark/flat | 자작 numpy | 합성 truth·ccdproc 비트일치·IRAF·PTC·LCO 교차기기 | ✅ Fig10-13 (Fig10 AIPPI패널 제거) |
| 1 스캔 | 헤더/타겟 | astroquery SIMBAD | 부기 — 서술만 | — |
| 2 크롭 | 픽셀 슬라이스 | 자작 | 부기 (WCS 보존 note) | — |
| 3 스카이QC | imexam 뷰어 | Qt/numpy | 부기 (프레임QC는 Fig6) | — |
| 4 검출 | segm/DAO/sep 3엔진 | **photutils/sep** | 합성 완전도 주입→복원 | ✅ Fig1 (라벨 개정) |
| 5 WCS | 플레이트솔브 | **자작 quad**① + astrometry.net/ASTAP 옵션 | **Gaia 잔차(이미 산출됨) + 내부↔외부 솔버 일치** | 🔶 신규 그림 |
| 6 마스터 | Gaia 교차매치·union dedup | astropy match_to_catalog_sky | CMD/Gaia 일치로 간접 | ✅ |
| 7 구경측광 | forced + apcorr② | photutils aperture + 자작② | 합성 truth + IRAF/sep + apcorr 복원 | ✅ Fig2,4,5 + 🔶 IRAF 매칭 재실행 |
| 8 PSF | ePSF 반복피팅 | **photutils** EPSFBuilder/IterativePSF | 구경↔PSF 일치 (M5·M13 crowded) | ✅ Fig9 |
| 9 ID편집 | 큐레이션 | pandas/Qt | 부기 | — |
| 10 ZP | ZP+색항 robust fit | numpy lstsq (자작 sigma-clip) | Gaia·PS1 독립 대조 | ✅ Fig7,8 (all-APEX 재생성) |
| 11 CMD플롯 | 뷰어 | Qt | 부기 | — |
| 12 이소크론 | emcee+PARSEC | emcee | **scope-out** | ⛔ |
| LC8 선택 | 타겟/비교별 | pandas/Qt | 큐레이션 | — |
| LC9 LC빌드 | 시계열+BJD_TDB+airmass | **astropy** light_travel_time | BJD 1점 sanity (외부 계산기 대조) | 🔶 |
| LC10 detrend | 차등·ensemble·SYSREM④ | 자작④ | **합성 주입신호 보존 + 계통 제거 복원** | 🔶 신규 |
| LC11 주기 | LS/BLS + PDM③ | astropy LS/BLS + 자작③ | **PDM↔LS 일치 + 문헌 주기 복원** | 🔶 신규 |

## Step7 측광 방법론 (코드 검증 완료 — §2 서술의 근거)

작은 구경(1.0×FWHM, `method="exact"`) flux → 로컬 스카이 annulus(6–9×FWHM,
sigma-clip) → **프레임별 성장곡선 apcorr**(SNR≥40·고립·비포화 참조별 ≤250개,
apcorr=1/enclosed_fraction; NGC6811 V 실측 1.4295) 적용 →
mag_inst = zmag − 2.5·log₁₀(F·apcorr/t_exp) → Step10 프레임별
ZP = median(m_Gaia − m_inst) + 색항. **aperture → apcorr → ZP 순서 표준 확인.**
apcorr는 프레임당 단일 스칼라(전역 PSF 가정) — §5 한계에 한 줄 명시.
Step10은 mag_psf를 우선 탐색하므로 Step8 실행 시 PSF mag가 자동 사용됨.

---

## 목차 (v2 — v1 6장 구조 계승, §2에 파이프라인 방법 통합)

```
§1 Introduction (~900w)
   진입장벽 서사(비인격화, HOPS 훅→IRAF 절벽→GUI 대안 실패) → 도구 지형
   2축 갭 → APEX 한 문단 + 기여 목록 + AI-assisted 한 문장 → 로드맵
   기여: (i) raw→science 통합 가이드 워크플로 (ii) CMD+시계열 2모드
        (iii) 컴포넌트별 검증 스위트 (iv) 공개 SW + 헤드리스 재현 코어
§2 Design & Implementation (~1500w)
   아키텍처(GUI/Qt-free 코어 분리, 2모드 12/11스텝, 공유 0-7) ·
   설계철학("신뢰 라이브러리 래퍼 + 자작 코어 4개" — Table 1) ·
   스텝별 방법 서술+수식:
     보정 ((raw−B)−k·D)/F · 검출(photutils segm/DAO/sep) ·
     천체측정(quad 솔버, Gaia QC 잔차) · 마스터(union dedup, 하이브리드 ID) ·
     구경+apcorr(위 방법론 절) · PSF(ePSF) ·
     표준화(m_std=m_inst+ZP+c·color[+c₂·color²]; Riello+2021, Pancino+2022) ·
     시계열(BJD_TDB Eastman+2010 · SYSREM Tamuz+2005 · PDM Stellingwerf 1978 ·
     LS VanderPlas 2018 · BLS Kovács+2002 · prewhitening Lenz&Breger 2005)
   설치·배포(pip, MIT) · GUI 스크린샷 1컷 + 워크플로 다이어그램
§3 Validation (~2800w) — 증거 코어, 컴포넌트별
   3.1 보정: 합성 truth·ccdproc 비트일치·IRAF·PTC(gain 0.681, 헤더 EGAIN 16×
       오류 발견 사례)·LCO 2카메라 교차기기
   3.2 검출·완전도: 주입→복원 m50 (로지스틱=경험적 요약임을 명시)
   3.3 천체측정: Gaia 잔차 분포 + 내부 quad 솔버 vs astrometry.net 일치 [신규]
   3.4 구경측광: 오차모델 pull(CCD식=이론 예측) · sep · IRAF(파라미터 매칭표
       T2) · apcorr 복원 · 파라미터 민감도
   3.5 PSF·crowded: 구경↔PSF 일치 (M5·M13, 코어 밀도 34×)
   3.6 표준화·CMD: Gaia·PS1 독립 대조 (BP faint bias 사례 = 정직성 스파인)
   3.7 시계열 코어: SYSREM 주입복원 · PDM↔LS · BJD sanity [신규]
   3.8 프레임 QC: 주입 결함→판정 (+투명도 blind spot 정직 보고)
§4 End-to-end science demonstrations (~900w)
   4.1 성단 CMD 갤러리: all-APEX 5성단 (NGC6811·M67·M13·M3·M5; 단일카메라 명시)
   4.2 변광성 주기: YZ Boo(P=0.104092d; 주력) · AE UMa(0.086017d; prewhitening
       필요시 future work 게이트 유지 — v1 결정 계승)
   4.3 워크드 예제: raw→CMD A4 1p + raw→주기 A4 1p (부록; 접근성 시연 =
       설문의 대체물)
§5 Discussion & Limitations (~1000w)
   확립된 것 / 도구 비교표(T5) / 한계: 단일기기·멀티앰프 미지원(Sinistro
   사례)·이소크론 scope-out·sub-resolution blending·apcorr 전역 스칼라 /
   실용 권고(구경 1.0-1.2×FWHM, 2단 QC, faint 레퍼런스 선택)
§6 Conclusion (~350w)
Back matter: Data&Code Availability · AI disclosure · CRediT · 부록 A4 예제 2장
```

## 그림 플랜 (최종 세트)

**컨벤션 (전 그림 공통):**
1. 겹친 곡선은 **수식 명시 + 정체 라벨** 필수 — 3분류:
   「이론 예측(자유파라미터 0)」/「이론기반 적합(적합값±σ)」/「경험적 요약」
2. 데이터 출처 배너: 「대상+기기+날짜」 (예: NGC 6811 · C3-61000 · 2026-06-11)
3. AIPPI 흔적 금지. all-APEX 재처리본만.
4. 스타일: `apex_paper_style.py` (Okabe-Ito, PDF+PNG)

| # | 그림 | 원본 | 작업 |
|---|---|---|---|
| F1 | 워크플로/아키텍처 다이어그램 + GUI 1컷 | 신규 | 제작 |
| F2 | 완전도 m50 (주입→복원) | Fig1 | 라벨 개정("경험적 요약"); erf 이론유도 격상은 선택 |
| F3 | 오차모델 pull (σ=1.0857/SNR 이론 예측) | Fig2 | 유지 |
| F4 | 파라미터·관측조건 민감도 | Fig3 | 유지 |
| F5 | sep 대조 (합성 truth) | Fig4 | 유지 |
| F6 | **IRAF 대조 (all-APEX NGC6811, 파라미터 매칭)** + T2 | fig_apex_iraf | **매칭 재실행**(ap 1.0·ann 6·dann 3; E: 경합 시 C: 스크래치 경유) — NGC457은 보조 |
| F7 | **WCS 검증: Gaia 잔차 + 내부↔외부 솔버 일치** | 신규 | 제작 (지표 이미 산출됨) |
| F8 | 프레임 QC (주입결함→판정) | Fig6 | 유지 |
| F9 | 구경↔PSF crowded (M5·M13) | Fig9 | 유지 (all-APEX 재측광 갱신 가능) |
| F10 | 검출기보정 before/after + cosmetic | Fig10 | **AIPPI 패널(c) 제거** → ccdproc/합성 대체 |
| F11 | PTC gain·RN·dark | Fig11 | 수식+"이론기반 적합" 라벨 보강 |
| F12 | ccdproc per-step 비트일치 | Fig12 | 유지 |
| F13 | LCO 교차기기 (QHY600·Sinistro) | Fig13 | 유지 |
| F14 | Gaia·PS1 대조 (BP faint bias) | Fig7 | **all-APEX 재생성** (PS1 위치 재매칭) |
| F15 | CMD 재현 (all-APEX vs Gaia + PS1) | fig_apex_cmd+Fig8 | 병합 재생성 (ridge RMS 16mmag 확보) |
| F16 | **5성단 CMD 갤러리 (all-APEX)** | 신규 | M3·M5 배치 완료 후 |
| F17 | **SYSREM 주입복원 + PDM↔LS 일치** | 신규 | LC 재처리 후 |
| F18 | **YZ Boo(·AE UMa) 위상접힘 LC + 주기 vs 문헌** | 신규(fig_lc_yzboo 참고) | LC 재처리 후 |

**표:** T1 스텝×라이브러리×검증(§2) · T2 IRAF/DAOPHOT 파라미터 1:1 매칭(§3.4)
· T3 데이터 출처(사이트+기기+기간) · T4 검증결과 요약(컴포넌트당 지표 한 줄)
· T5 도구 비교 매트릭스(§5)

## 인용할 이론·방법 (ars-lit-review 입력)

측광·검출: Stetson 1987 (DAOPHOT) · Bertin&Arnouts 1996 (SExtractor) ·
Barbary 2016 (sep) · Bradley+ (photutils) · Howell (CCD equation/오차모델)
천체측정: Lang+2010 (astrometry.net quads) · Gaia DR3 (Vallenari+2023)
표준화: Riello+2021 (Gaia 밴드변환) · Pancino+2022 (Gaia→Johnson B) ·
Chambers+2016 (PS1)
시계열: Stellingwerf 1978 (PDM) · Tamuz+2005 (SYSREM) · VanderPlas 2018 (LS) ·
Kovács+2002 (BLS) · Lenz&Breger 2005 (Period04) · Eastman+2010 (BJD_TDB) ·
Yang+2018 (YZ Boo)
보정·검출기: van Dokkum 2001 (L.A.Cosmic/astroscrappy) · ccdproc (Craig+) ·
Alarcón+2023 (IMX455 특성) · BANZAI (McCully+)
지형: Tody 1986/1993 (IRAF) · Collins+2017 (AstroImageJ) · Brennan&Fraser 2022
(AutoPhOT) · Astropy Collaboration 2013/2018/2022 · emcee (Foreman-Mackey+2013)

## 남은 작업 (우선순위)

1. **F7 WCS 그림** — 지표 이미 있음, 그림화만 (반나절)
2. **F6 IRAF 매칭 재실행 + T2 표** (C: 스크래치 경유, 1시간)
3. **F16 5성단 갤러리** — M3·M5 배치 완료 대기 (반나절)
4. **F10 AIPPI 패널 제거 + F2 라벨 개정 + F11 수식 라벨** (반나절)
5. **LC 재처리** (`reprocess_batch.py --lc`, E: 236GB 확보됨) → F17·F18 (하루+)
6. **F14/F15 all-APEX 재생성** (PS1 위치 재매칭 포함)
7. F1 다이어그램 + T1-T5 표 + MANUSCRIPT.md 목차 개정 (이 플랜 반영)
8. A4 워크드 예제 2장 (부록)
9. ars 스킬시퀀스 (v1 계승): ars-lit-review → ars-outline → 초안 →
   ars-citation-check → ars-disclosure → ars-format-convert → ars-reviewer

## 현재 데이터 상태 (2026-07-10)

- all-APEX CMD: NGC6811(2084) · M67(1288) · M13(1543) ✅ / M3·M5 배치 진행 중
- LC 원본: AE UMa 597f · YZ Boo 879f (Step0 대기; E:\observed_Analysis FITS
  삭제로 236GB 확보, D:\APEX_backup이 유일 백업 — HDD 유지)
- 검증 자산: Fig1-13 + fig_apex_cmd(ridge 16mmag) + fig_apex_iraf(MAD 0.0097 —
  단 구경 0.8 vs 1.0 불일치, 매칭 재실행 예정)
- IRAF 인프라: apex/benchmark/iraf_crosscheck_cli.py (WSL pyraf 2.2.3 동작)
