# 완전도(completeness) 실측-프레임 주입 조사 — 2026-07-23

## 배경
리뷰(LIT_REVIEW_R2 axis-2, Haynes 2002 MNRAS 334,262)가 "합성 프레임 주입은 비표준·낙관 편향"
이라 지적 → §3.6 완전도를 **실측 프레임 주입**(DAOPHOT ADDSTAR / DES Balrog / AutoPhOT 표준)으로
전환하라는 권고. 이를 실행하며 수치를 추적한 기록.

## 실행
`apex/benchmark/validate.py:run_artificial_star_suite(reference_frame=<실측 FITS>)` 로 M13 V 프레임
(`E:\APEX_validation\reprocess\M13\calibrated\20260515\pp_messier13-0001-V.fit`, 3194×4788, 60s)에
인공별 주입. 주입 등급범위 11–17 (module const monkeypatch), trials=60, 50별/trial, bootstrap=300.
출력: `validation/paper/data_realframe_M13V/`.

## 결과 (수치)
- 실측 M13 완전도: **m90≈14.36, m50≈14.84, m10≈15.33** (기기등급, count-rate ZP=25).
- 완전도 곡선이 mag 14.5→15.5에서 **절벽**(0.85→0.02) — 매끈한 photon-noise erf 롤오프 아님.

## 추적 (철칙 #12 — 설명 안 되는 수치)
격차 원인을 끝까지 분해:

| 항목 | 합성 프레임 | 실측 M13 V |
|---|---|---|
| 배경 median | 150 ADU | 1315 ADU |
| 배경 노이즈(SEP bkg rms) | 11.3 ADU | 43.6 ADU |
| **별 FWHM** | (좁음, 이상적) | **7.64 px (중앙값, 실측)** |
| 주입 PSF FWHM | — | **7.40 px** ← 프레임에서 측정한 경험적 PSF |
| m50 | ≈17.5 | 14.84 |

**핵심 규명**:
1. 주입 flux 스케일은 **정확** — mag 15 주입별 flux_expected=7861 e 는 실측 mag-15 별과 동일
   (실검출 중앙값이 mag 14.87). ZP 불일치 아님.
2. 절벽의 물리적 원인 = **넓은 seeing(FWHM 7.6px)**. 총 flux 7861 e가 7.4px PSF로 퍼지면
   **peak = 84 e**, 검출 임계(3.2σ ≈ 96 e) 아래 → 검출 실패. mag 14.5(peak 134 e)에서
   ~50%, mag 14(peak 216 e)에서 회수. peak가 임계를 넘느냐가 절벽을 만든다.
3. 주입 PSF(7.40px) = 실측 별 FWHM(7.64px) → **주입 정확, 버그 아님**. 이 프레임이 실제로
   얕다(넓은 seeing + 높은 하늘).
4. **독립 교차검증**: 파이프라인 자체 clean 검출(875개, 품질필터)의 faint 롤오프 —
   90%가 mag<15.33, 50%가 <14.53. 주입 완전도(m50=14.84, m10=15.33)와 **정확히 일치.**
   → 주입-회수 완전도가 경험적 검출 한계와 독립적으로 일치 = 강력한 validation.
   (필터 없는 raw sep.extract가 flux 446 e까지 잡던 건 노이즈 잡음 — 파이프라인이 걸러냄.)

## 판단
- 실측 주입은 **올바르게 작동** → §3.6을 실측 프레임 주입으로 전환 가능(리뷰 권고 해소).
- 그러나 합성 m50≈17.5 vs 실측 14.84의 **2.7 mag "격차"는 이 M13 프레임의 나쁜 seeing이 주범**
  (합성은 이상적 PSF+저노이즈). "합성이 2.7 mag 낙관적"이라 라벨하면 **오해** — 대부분이
  프레임별 seeing 차이지 일반적 "합성 낙관"이 아님. Haynes식 정직한 재현은 seeing·하늘을 실측에
  **맞춘 합성 대조군**이 필요(향후 확장으로 남김, 오늘은 안 함).
- **채택 프레이밍**: verification(합성 = 기계·수치 정확성, known truth) → validation(실측 = 실제
  성능, 경험적 검출한계와 교차확인). 원고의 Oberkampf&Trucano/Portillo 3-rung 사다리와 정확히 매핑.
  Δm50 화살표·"낙관적" 문구 제거.

## 다중프레임 확장 (2026-07-23, 사용자 지적 — "M13 하나면 대표성 없다")
M13-0001-V(7.64px)는 M13 프레임 중에서도 나쁜 seeing(M13 중앙값 6.94, 최선 5.9). 대표성 우려 →
좋은-seeing 개방성단 프레임 2개에 동일 주입 추가:

| 프레임 | sky(ADU) | FWHM(px) | m50 | 깊이 |
|---|---|---|---|---|
| M67 i (`sci/pp_Messier67-0008-i.fit`) | 27 | 5.2 | **17.65** | 깊음 |
| NGC6811 R (`sci/pp_NGC6811-0005-R.fit`) | 1315 | 5.3 | **15.64** | 중간 |
| M13 V (`calibrated/…/pp_messier13-0001-V.fit`) | 1315 | 7.6 | **14.90** | 얕음 |
| 합성(verification) | 150 | 3.4 | 17.59 | (이상적) |

**결과가 물리와 정확히 일치**: M67→NGC6811 2.0mag차 = 하늘 50배(√ 노이즈 ~1.9mag), NGC6811→M13
0.74mag차 = seeing 5.3→7.6px. **주입이 프레임 품질(sky+seeing)을 정확히 추적** → M13은 이상치 아님,
그냥 제일 얕은 프레임. **결정적**: 어두운-하늘 샤프 프레임(M67 i, m50 17.65) = 합성 깊이(17.59)에 도달
→ "합성이 낙관적"이란 우려는 **확정적으로 반증**(실측 좋은 프레임이 합성을 따라잡음). 데이터:
`data_realframe_{M67i,NGC6811R,M13V}/`. 실측 주입별 컷아웃 사다리(M13, `injection_cutouts.npz`)를
그림 하단에 병기(`make_injection_cutouts.py`, inject_flux_catalog 정규경로).

## 산출물
- `fig_completeness_realvssynth.py` → 4곡선(실측 3 + 합성) 다중프레임 사다리 + 컷아웃 하단 스트립.
- 데이터: `validation/paper/data_realframe_{M67i,NGC6811R,M13V}/` (실측 3), `data/artificial_star/` (합성).
- **파이프라인 반영 (2026-07-24)**: S/N50 법칙을 APEX 본체의 프레임별 depth QC 게이트로 구현.
  `apex/analysis/detection_limit.py` (`predict_frame_m50` 등, Qt-free) + step7 `frame_stats.csv`에
  `predicted_m50`/`observed_m50`(마스터카탈로그 검출률 50% 롤오프)/`depth_delta_mag`/`depth_qc_flag`
  (기본 허용 0.5 mag, `depth_qc_tolerance_mag`). 상수 `PEAK_SN50_DETECTION=4.05`는
  `apex/utils/constants.py`, gain은 런타임 config. 검증: `tests/test_detection_limit.py` —
  7런 재예측 잔차 RMS 0.048 mag; 실 NGC6811 3프레임 step7 재실행에서 predicted vs 주입 실측
  −0.008/+0.109 mag, observed 롤오프와도 전부 |Δ|<0.17 (모두 "ok"). '예측력 실증' 절(실별
  검출수 6% 일치)이 predicted↔observed 비교의 근거.
- **실데이터·GUI 검증 (2026-07-24 확장)**: M67 9프레임(g/r/i, 60s) 추가 — 전 프레임
  |predicted−observed| ≤ 0.073 mag, 주입 캘리브레이션 3프레임과 +0.111/+0.049/−0.038.
  GUI(CMD Step7)에서 동일 3프레임 실행: 수치 헤드리스와 비트동일, TSV·frame_stats 정상 저장,
  허용오차 0.05로 조이면 depth_shallow/depth_deep 플래그+로그 경고 발동 확인.
  파라미터는 `[photometry.depth_qc] tolerance_mag/min_snr`로 노출(cmd/lc 공통), step7
  캐시 서명에 포함(변경 시 자동 재계산 확인). 관측 롤오프는 최대완전도 bin부터 스캔
  (포화 밝은별 dip의 허위 조기교차 방지).

## 후속 정정 (2026-07-24, 사용자 의심 "어떻게 M67이 딱 맞지?" 적중)
- 합성 프레임의 S/N50을 직접 계산: **5.19** — 실측 7프레임 법칙(4.05±0.18) 위에 **없음**.
- → **"M67 실측=합성 깊이라서 합성은 낙관 아님" 주장 철회**: m50 일치(17.59 vs 17.65)는
  합성의 높은 노이즈(11.2 vs 5.1 e-)를 샤프한 PSF(p_peak 0.064 vs 0.025)가 상쇄한
  **파라미터 우연**. 조건의 대응이 아님.
- 대신 **새 체계 발견**: S/N50이 커널 FWHM에 단조 의존 (5.3px→4.3, 7–9px→3.8-3.9,
  합성 3.4px→5.2; Spearman ρ=-0.982). 기전 = 검출의 최소 픽셀면적 요구(minarea=5 등) —
  샤프한 프로파일은 적은 픽셀에 flux 집중 → 같은 peak S/N에서 면적 기준 미달.
  minarea-5 나이브 모델은 방향·순서 재현, 크기는 과소(매치드필터·매칭반경 미포함) — 정량 모델 미완.
- 정직한 결론: **S/N 기준으로 합성은 낙관이 아니라 약간 보수적**(더 높은 peak S/N 요구).
  캡션에서 "M67이 합성 따라잡음" 문구 제거, FWHM 트렌드 서술로 교체.

## 완전도 곡선의 예측력 실증 (2026-07-24, 사용자 질문 "실제 검출수 예측 되는거네?")
NGC6811 R soft(480s) 프레임: 마스터 카탈로그 실별 1,972개(off-frame 제외)에 주입 완전도
C_inj(m)을 적용한 예측 검출수 **Σ C_inj(m_i) = 1,629** vs 실제 detected_flag **1,728** —
**6% 일치**. (등급 변환: mag_inst는 count-rate라 주입 스케일로 −2.5log10(480) 시프트,
TSV flux_e는 apcorr 반영 총전자 확인.)
- bin별: 밝은 평지 일치(실별 0.95 vs 주입 0.97 — 실별은 블렌드·플래그로 약간 낮음).
  전이구간(14.6)에서 실별 회수율이 높게 보이는 건 **측정등급 binning의 Eddington 편향**
  (검출된 별은 위로 요동한 것들) — 적분은 bin 이동이 상쇄돼 강건.
- 논문 활용: "주입 완전도가 실별 검출수를 6% 이내로 예측" = (b) 붕괴와 독립적인 3번째 검증축.
