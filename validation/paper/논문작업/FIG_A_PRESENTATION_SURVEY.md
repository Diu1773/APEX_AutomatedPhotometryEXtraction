# fig (a) — 등급공간 완전도의 논문 배치 조사 (2026-07-24)

사용자: "(a)를 어떻게 논문에 녹일지 좀 더 생각. AutoPhOT 말고 다른 논문들도 봐야 할 듯."
웹조사로 확인한 4계열 관행 + APEX 적용 옵션. (그림 번호·세부는 원고 인용 전 citation-check 필수.)

## 관행 4계열

**① 도구 논문 (AutoPhOT, SExtractor 계열)** — 예시 프레임 1개의 완전도 점 + 한계등급 인용.
(a)의 현재 형태가 이 계열. 역할 = "도구가 내놓는 산출물의 예시".

**② 서베이 전달함수 논문 (DES Balrog Y3=Everett+2022 ApJS / Y6=arXiv:2501.05683)** —
완전도를 서베이 전체에 걸쳐 집계, **깊이의 분포**(필드별 m50 지도/히스토그램)가 산출물.
Balrog Y6 Fig3: griz별 검출확률 vs 등급 + 90% 완전도 수직선. 개별 프레임 곡선이 아니라
**조건별 변동 자체를 상품화**.

**③ 운영 모니터링 (DES-SN DiffImg=Kessler+2015 AJ 150,172 / ATLAS=Tonry+2018 PASP 130,988)** —
완전도를 1회성 그림이 아니라 **상시 QC 지표**로: DES-SN은 fake SN을 매 이미지에 상시 주입해
단일이미지 검출효율을 모니터, ATLAS는 **알려진 소행성**(외부 정답지!)의 회수율로 망원경·밤별
효율 곡선(Fig3, 효율 vs V)을 운영 중 측정. → 우리 "1728 실별 체크"가 정확히 ATLAS 방식의 축소판.

**④ 성단/CMD AST 전통 (DOLPHOT/HST·JWST 계열)** — 완전도를 등급·색의 함수로 재서
**CMD 위에 50%/90% 완전도 선**으로 얹거나 광도함수 보정에 직접 소비. 로지스틱 피팅 m_c 관례.

## APEX 적용 옵션

- **A (현상 유지, ① 계열)**: §3.6 = 현재 2패널(예시 3프레임 + S/N붕괴). (a)는 pedagogy,
  (b)가 주장. 최소 노력. 리스크: "(a)가 왜 3개뿐?"류 질문은 캡션 방어.
- **B (② 계열 추가)**: 법칙으로 **66프레임 전체의 예측 m50 분포**(frame_stats의 σ·FWHM만으로)
  를 스트립/히스토그램으로 — 주입 7프레임을 마커로 겹침. "도구가 데이터셋 전체의 깊이를
  산출"을 보여줌. 비용 낮음(법칙 이미 검증됨).
- **C (③ 계열, QC 게이트와 결합)** ★추천: predicted-m50 QC 게이트(백로그 task_5bb4af3e)
  구현 후 **"예측 vs 실현 깊이, 66프레임"** 산점 = Kessler/ATLAS식 운영 검증 그림.
  (a)는 §3.6에 예시로 남고, 서베이-스케일 깊이 이야기는 QC 절(§3.9 인근)로 이동.
  실별 회수 체크(6% 일치)도 이 절에 한 문장으로.
- ④는 §3.6이 아니라 CMD 과학 절에서 (50% 선을 CMD에 얹기) — 별도 후보.

## 판단 보류
사용자 "좀 더 생각해봐야" — **A/B/C 선택은 미결**. C는 QC 게이트 구현이 선행.

Sources (원고 인용 전 재검증):
- https://arxiv.org/abs/2501.05683 (Balrog Y6) · https://iopscience.iop.org/article/10.3847/1538-4365/ac26c1 (Balrog Y3)
- https://iopscience.iop.org/article/10.1088/1538-3873/aabadf (Tonry 2018 ATLAS)
- Kessler+2015 DiffImg (DES-SN fakes 상시주입 모니터링)
- https://iopscience.iop.org/article/10.3847/1538-4365/ad2600 (JWST ERS DOLPHOT)

## 원문 그림 확인 (2026-07-24, reffigs/에 사본)
Kessler+2015 (arXiv:1507.05137, ar5iv 렌더)에서 **우리 설계의 1:1 선례 3장** 확인:
- **Fig 6** (`kessler2015_fig6_perepoch_maghist.png`): 에폭 1개의 fake 등급분포 + 검출된 것
  음영 + m_eff=1/2 점선 — **우리 (a)의 프레임당 판독과 동일 논리** (deep/shallow 필드 병렬).
- **Fig 7** (`..._fig7_m50_distribution.png`): "각 항목 = 에폭 1개"인 **m_eff(0.5)의 분포
  히스토그램**(밴드×deep/shallow) — **우리 realized_m50.csv 66프레임의 정확한 선례.**
  C 그림의 실현-깊이 축은 이 형식(분포)으로 제시 가능. 단 Kessler는 예측 없이 모니터만 —
  **우리의 법칙 기반 predicted-vs-realized는 Kessler를 넘어서는 기여점.**
- **Fig 8** (`..._fig8_eff_vs_SN_splits.png`): 검출효율 ε vs **계산 S/N** + ZP/PSF/σ_sky
  절반 분할 곡선 — **우리 (b) S/N 붕괴의 직접 선례.** PSF 분할에서 곡선이 살짝 갈라짐 =
  우리 FWHM 트렌드와 동일 현상 (부호는 축 정의 차이: Kessler는 총-flux S/N이라 broad가
  오른쪽, 우리는 peak S/N이라 sharp가 오른쪽 — 같은 물리).
- **Balrog Y6 Fig 3** (`balrog_y6_fig3_detrate_griz.png`): 서베이 전체 집계 검출률 vs mag,
  90% 수직선. 전이가 우리보다 완만한 건 전 footprint 조건 분산이 합쳐진 것 — 집계형의 특징.
→ 결론: C 그림 = (i) per-frame m50 분포(Kessler Fig7 형식) + (ii) predicted vs realized 산점
  (우리 추가 기여) 2패널이 관행과 기여를 모두 만족. (b)에 FWHM-분할 곡선 옵션(Kessler Fig8 방식).
