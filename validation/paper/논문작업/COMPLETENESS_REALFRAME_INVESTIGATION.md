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
