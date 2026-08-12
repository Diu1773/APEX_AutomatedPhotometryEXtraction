# APEX step 8 과 IRAF DAOPHOT ALLSTAR 의 측정 파라미터

이 표가 없으면 두 엔진의 비교는 논문에 쓸 수 없다. DAOPHOT 열은 실행이
IRAF `lpar` 로 덤프한 값(`iraf_parameters.json`)이고 — 하네스가 설정하려
**의도한** 값이 아니라 태스크가 발화할 때 IRAF 가 **실제로 들고 있던** 값이다.
APEX 열은 그 step 8 산출을 만든 워크스페이스 설정에서 읽는다.

대상 프레임: `pp_messier13-0005-B.fit` · 좌표 1701 개 · ALLSTAR 적합 1659 개 · 554 s

## 검출기 상수 — 두 엔진이 같아야 하는 것

자료의 성질이다. 다르면 비교가 엔진이 아니라 입력의 불일치를 잰다.

| 항목 | APEX step 8 | IRAF DAOPHOT | 비고 |
|---|---|---|---|
| Gain | 0.68 | `0.68` | 동일 — APEX 설정에서 두 엔진에 같이 준다 |
| Read noise | 2.35 | `2.35` | 동일 |
| Good-data minimum | 0.1 | `0.1` | 동일 |
| Good-data maximum | 55000 | `55000.0` | 동일 |
| FWHM | 7.05239 | `7.052391` | 동일 — 프레임에서 측정한 값 |
| Sky sigma | 29.1113 | `29.11127622` | 동일 — sigma_clipped_stats |
| Exposure | 60 | `60.0` | 동일 — 헤더 |

## 방법 파라미터 — 각 엔진 저자의 선택

DAOPHOT 에 APEX 값을 강제하면 APEX 의 답 쪽으로 조율하는 것이 된다.
표준 지침(Massey & Davis 1992)을 따르되 **차이는 전부 적는다.**

| 항목 | APEX step 8 | IRAF DAOPHOT | 비고 |
|---|---|---|---|
| PSF model radius | 4.0·FWHM 상자 | `29.209564` | DAOPHOT 4·FWHM+1 (Massey & Davis 1992) · APEX epsf_size_fwhm_mult=4 의 상자 |
| Fit radius | auto · 포위에너지 0.9 | `11.9890647` | **차이** — APEX 는 포위에너지 90 % 자동창(약 1.7·FWHM). DAOPHOT 표준은 1·FWHM 이나 여기서는 APEX 에 맞춰 1.7 로 돌렸다(교란 제거 확인) |
| Analytic function | epsf | `"auto"` | **차이** — DAOPHOT 은 auto 선택(실제 채택: moffat15), APEX 는 경험적 ePSF |
| Spatial variation | per_frame | `0` | **차이** — DAOPHOT varorder=0(일정), APEX 는 프레임당 ePSF 하나. varorder=1 로도 돌려 결론이 안 바뀜을 확인 |
| Recenter during fit | yes (적합 내) | `yes` | 동일 취지 — 두 엔진 모두 적합 중 위치를 미세조정 |
| Fit sky | yes | `yes` | 동일 취지 |
| PSF-star cleaning | 5 | `0` | **차이 — DAOPHOT 에 불리하다.** nclean=0 이라 PSF 별 정제 반복을 껐다. APEX 는 epsf_maxiters=5 로 ePSF 를 반복 정련한다 |
| Max fit iterations | 8 | `50` | 각자 기본값 |
| PSF reference stars | 자동 (오염인지 필터) | `60` | pstselect 로 자동 선정 · APEX 는 오염인지 필터로 선정 |
| Initial aperture | Step 7 구경보정 경로 | `"7.0524"` | DAOPHOT 은 phot 의 초기 등급용. ALLSTAR 등급의 영점이 여기서 오므로 밝고 한산한 별로 상수 오프셋을 제거한 뒤에만 bias 를 비교한다 |
| Sky annulus | 해당 없음 (적합 내 하늘) | `28.209564` | **구조가 다르다** — DAOPHOT 은 고리에서 하늘을 재고, APEX step 8 은 하늘을 적합의 자유변수로 함께 푼다. 구멍 반경이라는 개념이 없다 |
| Sky annulus width | 해당 없음 | `14.104782` | 위와 같음 |
| Centering in `phot` | 강제 (좌표 고정) | `"none"` | **none 으로 껐다.** 좌표를 주는 강제측광이라 재중심을 켜면 희미한 별의 중심이 이웃으로 끌려간다(주입별 25 중 7 만 생존한 실측으로 발각) |

## 표에서 읽어야 할 것

- **`nclean = 0` 은 DAOPHOT 에 불리한 설정이었고, 켜서 확인했다.** PSF 별 정제
  반복을 끄면 `pstselect` 가 고른 60 개가 그대로 모형에 들어간다(`psf.pst` 60 →
  `psf.opst` 60, 기각 0). 사람이 눈으로 걸러내던 단계의 자동 대체물을 안 쓴 것이다.
  `--nclean 5` 로 다시 돌려도 **DAOPHOT 이 여전히 정밀도에서 앞선다** — 같은 별
  기준 혼잡 사분위 산포가 오히려 더 조밀해진다(0.75–1.5 FWHM 0.060 → 0.038,
  6–inf FWHM 0.031 → 0.019, ALLSTAR). 즉 `nclean = 0` 결과가 이 방향으로
  보수적이라는 추정이 실측으로 확인됐다.
- **적합창·공간변화·PSF 정제 셋 다 따로 실험해서 닫았다.** `--fitrad-fwhm 1.7`,
  `--varorder 1`, `--nclean 5` 로 각각 재실행했고 결론이 바뀌지 않았다.
- **초기 구경이 ALLSTAR 등급의 영점을 정한다.** 그래서 bias 를 비교하기
  전에 밝고 한산한 별에서 상수 오프셋을 제거한다(측정값 +0.296 mag,
  1·FWHM 구경보정 크기와 일치).

전체 파라미터(설명 포함)는 `*_iraf_parameters.json` 에 있다. 이 표는
그 중 비교에 영향을 주는 것만 뽑은 것이다.
