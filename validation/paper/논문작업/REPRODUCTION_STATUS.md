# 레퍼런스 검증-프레임 재현 현황 — /goal 자율 배치 (2026-07-23)

**목표**(사용자 /goal): "각 주요 논문들의 검증 프레임을 그대로 할 수 있는 만큼 레프로듀싱."
= 참조 논문들이 파이프라인을 검증한 방식(그림·방법)을 APEX **실측 데이터**로 재현.

원칙: 실측 데이터·정규 코드경로만(철칙 #5), 설명 안 되는 수치는 추적(#12), 논란 영역(Gaia-CMD
보정계·이소크론) 회피, 각 그림 message-first 제목·내부 ID 노출 금지(FIGURE_SPEC).

## ⚠️ "재현"의 정의 — 기기간 수치 비교가 아님 (2026-07-23 확정)

레퍼런스들은 **전부 다른 기기**다: Honeycutt=초기 CCD, Collier Cameron(WASP)=e2v CCD 광각,
Kovács(HATNet)=CCD, Ofek/Masci(ZTF)=대형 CCD 모자이크. APEX=Sony IMX455 **CMOS**+소형망원경
(Moravian C3-61000). 따라서:

- **비교 불가(기기 특정·절대값)**: floor 절대값(mmag), 한계등급/깊이, 광자곡선 위치. 구경·사이트·
  scintillation(∝D^{-2/3})·노출·검출기(CCD↔CMOS)에 좌우됨. **"APEX가 WASP만큼 정밀"류 문구 금지.**
- **재현 대상(기기 무관)**: (1) 검증 **방법/프로토콜**(주입-회수·앙상블·에러모델 검증·독립엔진 대조),
  (2) **보편 구조/법칙**(RMS-vs-mag 모양=photon 추종+floor 평평; 완전도 m50 ∝ σ·FWHM², 상수 C만 기기 흡수).
- **"재현"의 뜻** = 커뮤니티 표준 검증 배터리를 **APEX 자기 데이터·자기 기기**에 그대로 걸어 통과함을
  보이는 것 (AutoPhOT도 자기 데이터로, DES Balrog도 DES로). **수치는 APEX(C3-61000 CMOS) 기준으로만 보고.**
- **기기 의존도 순위**(강→약): ① IRAF 교차(같은 프레임 두 코드, 완전 무관) > ② 완전도 스케일링(법칙 무관,
  절대 m50만 기기별) > ③ 측성(arcsec 반쯤 물리, 광학 의존) > ④ 정밀도 floor(절대값 기기특정, 모양만 재현).
- 진짜 **기기간 일반성**은 별도 축에서: LCO Sinistro CCD + QHY600 CMOS vs BANZAI 교차기기(sub-percent).
  레퍼런스 재현과 역할이 다름 — 섞지 말 것. 원고: 전 성단 동일카메라라 다중성단≠다중기기도 명시.

## ✅ 완료 (이번 배치, 4개 — 모두 실측·커밋됨)

| 그림 | 재현 대상 논문 | 데이터 | 핵심 수치 | 파일 |
|---|---|---|---|---|
| **완전도(F5)** | Haynes 2002 / DES Balrog / DAOPHOT ADDSTAR / AutoPhOT App.D | 실측 3프레임(M67i·NGC6811R·M13V) 주입-회수 + 합성 verification | m50 17.7/15.6/14.9 (sky+seeing 추적), M67실측=합성깊이, 컷아웃 사다리 | `fig_completeness_realvssynth.py` |
| **측성(F6)** | Ofek 2019 / Masci 2019 (ZTF) | 실측 66프레임 3성단 WCS QC | Gaia 잔차 중앙값 0.26″, 66/66 solved | `fig_astrometry.py` |
| **독립엔진 교차(F10)** | Schechter 1993 (DoPHOT vs DAOPHOT) / AutoPhOT Fig14 | NGC6811 V, APEX vs IRAF phot, 499별 고정좌표 | MAD 9.7 mmag, r=0.99989, faint까지 평평 | `fig_iraf_crosscheck.py` |
| **정밀도 floor** | Honeycutt 1992 / Collier Cameron 2006 (WASP) / Kovács 2005 (TFA) | 실측 M67 10프레임 r, 1073별 | floor g/r/i = 5.4/5.2/6.7 mmag, 에러모델 일치 | `fig_precision_floor.py` |

각 그림의 정직성 포인트:
- **완전도**: 합성 vs 실측 2.7 mag 격차를 끝까지 추적 → M13 프레임의 넓은 seeing(7.6px)+높은
  하늘(1315 ADU)이 주범, "합성이 낙관적"이란 일반화는 오해라 제거. verification(합성=기계정확성)→
  validation(실측=성능) 사다리로 재구성. 상세: `COMPLETENESS_REALFRAME_INVESTIGATION.md`.
- **정밀도**: 경험적 RMS와 파이프라인 보고 mag_err가 photon-noise 영역에서 일치(에러모델 검증) +
  bright end에서 계통 floor 노출. Honeycutt 앙상블 ZP 보정(표준 후처리, 재구현 아님).

## ⏸ 보류 (데이터/영역 제약)

| 후보 | 재현 대상 | 막힌 이유 |
|---|---|---|
| 주기 복원(PDM/LS 주기도) | Stellingwerf 1978 / Graham 2013 | LC 데이터(AE UMa·YZ Boo) 전처리 deferred(100GB, observed_Analysis 삭제 대기). 합성신호 주입 verification은 가능하나 실측 아님 |
| 광도 보정해(ZP·색항·잔차) | Padmanabhan 2008 / Schlafly 2012 | Gaia-synthetic 보정계는 이소크론/B필터 디버깅 이력 있는 **논란 영역** — 사용자 판단 필요 |
| 소광(Bouguer, zp vs airmass) | Hardie 1962 | frame_zeropoint에 airmass 있으나 ZP 보정 레퍼런스에 얽힘. airmass span 미확인 |
| Curve-of-growth(구경보정) | Stetson / Howell | apcorr_summary는 프레임당 단일값만 — 다중구경 재실행 필요 |

## 데이터 인벤토리 (E:\APEX_validation\reprocess\)
- 다중프레임 forced phot: M13(V 15), M67(g/r/i 10×3=37), NGC6811(21). M5는 step7 없음.
- WCS QC: M13/M67/NGC6811 `step5_wcs/frame_wcs_qc.csv` (66프레임).
- IRAF 교차: `benchmark/runs/ngc6811_iraf_allapex_v1`, `ngc457_iraf_crosscheck_g0016_v1`.
- 실측주입: `validation/paper/data_realframe_M13V/` (M13 V, 3000별, m50=14.9).
- 보정(논란): `cmd_zeropoint/` (zp_fit_coefficients, gaia_sdss_calibrator, gaia_cmd_comparison).

## 사용자 판단 필요
1. **보류 항목 중 진행할 것**: (a) 광도보정 잔차 그림(논란 영역 — 진행하려면 어느 레퍼런스로
   보정할지 결정 필요), (b) LC 주기복원(전처리 재개 필요), (c) 소광 Bouguer, (d) curve-of-growth 재실행.
2. **원고 §3.6 반영**: 완전도를 실측-프레임 주입으로 전환(리뷰 권고 해소) — 프로즈 갱신할지.
3. **정밀도 floor 그림**을 §3 어느 절에 넣을지(신규 절 vs 기존 시계열 절).
