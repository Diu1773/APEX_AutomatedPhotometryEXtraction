# PSF 측광 분석 및 현황 정리

## 1. 결론 요약

이 데이터(지상 관측, FWHM 5~8px)에서 **step7 aperture photometry로 CMD가 충분히 나온다.**
PSF는 선택적 향상 도구이며, 아래 조건이 충족될 때 의미가 있다.

| 조건 | PSF 효과 | 현재 데이터 |
|---|---|---|
| 이웃이 aperture 안에 들어오는 별 | 유의미 | 일부 |
| 중심부 블렌드 분리 | 유의미 | **지상 한계로 불가** |
| 고립된 별 | 거의 없음 | 대부분 |

---

## 2. photutils가 안 된 이유

### 속도 문제
`IterativePSFPhotometry + SourceGrouper(max_size=25)`:
- 그룹당 75-param LM 최적화 (25별 × 3파라미터)
- 프레임당 13분+ → 실용 불가

### 근본적 접근 문제
- photutils는 grouped LM으로 동시 피팅 시도
- 그룹 크기가 크면 O(group_size³) 비용
- 정확도는 있지만 이 데이터 규모에서 속도 한계

---

## 3. 현재 구현: ALLSTAR (Newton step)

### 알고리즘
DAOPHOT ALLSTAR (Stetson 1987) 방식으로 재구현:

```
각 별마다:
  1. 이웃 별 모델 빼기 (neighbor subtraction)
  2. PSF 편미분 계산: ∂PSF/∂x, ∂PSF/∂y
  3. 3×3 정규방정식 한 번 풀기 → (dx, dy, dflux)
  4. 위치 + flux 업데이트

outer iteration:
  → 잔차에서 DAOStarFinder로 재검출
  → 새 소스 추가 후 재피팅
  → max_iter(=3)번 반복
```

### scipy LM 대비 개선
| | scipy LM | Newton step |
|---|---|---|
| 별당 비용 | 20~40회 반복 | 3×3 행렬 1회 |
| 안정성 | max_dflux 수십억 발생 | 안정적 |
| 프레임당 시간 | 60~75s (no group) | ~260s (group on) |

### 파라미터
- `fit_engine = "allstar"` — Newton step 사용
- `build_mode = "epsf"` — per-frame EPSFBuilder
- `group max size = 1` — 그룹핑 없음 (권장)
- `max_iter = 3` — outer iteration 횟수

---

## 4. 그룹핑 (Grouping)

### 구현 상태
`_build_groups` + `_allstar_newton_group` 구현 완료.
1.5×FWHM 이내 별들을 3N×3N 선형 시스템으로 동시 피팅.

### 이 데이터에서 효과 없는 이유
1. **코어 소스 미검출**: step4에서 블렌드로 묶인 소스는 ALLSTAR 입력에 없음 → 그룹핑해도 피팅 대상 없음
2. **잔차 재검출도 실패**: 블렌드 잔차는 점광원 형태 아님 → DAOStarFinder sharpness/roundness 기준 탈락
3. **step4 분리 검출된 별**: 이미 1 FWHM 이상 떨어진 별들 → sequential neighbor subtraction으로 충분

**결론: 이 데이터에서 group max=1이 실용적 최적.**

---

## 5. 중심부 측광 한계

### 왜 안 되는가
```
지상 2.5" seeing + M5 코어 밀도:
  → 별 간격 < 1 FWHM → PSF 완전 겹침
  → step4 deblend 한계로 개별 소스 분리 불가
  → PSF 피팅 입력 없음
  → 잔차도 블롭 형태 → 재검출 불가
```

어떤 PSF 알고리즘을 써도 **검출 단계에서 걸러지는 것이 근본 원인**.

### IRAF도 동일한 한계
DAOPHOT ALLSTAR + GROUP도 이 상황에서 동일한 잔차 패턴을 보인다.

---

## 6. 해결 방안

### 6-1. 단기 (현재 가능)
- **PSF skip**: step7 aperture 결과로 CMD 분석 → 이미 충분
- **group max=1 Newton ALLSTAR**: aperture 대비 minor contamination 개선, ~260s/frame

### 6-2. 중기 (구현 가능)
**Core ROI Forced Linear PSF Solver**

```
1. step6 master catalog의 x,y → 각 프레임으로 변환
2. 코어 ROI 선택 (클러스터 중심 반경 기준)
3. x,y 고정, flux만 미지수
4. sparse A 행렬 빌드: A[pixel, star] = PSF_i(pixel)
5. scipy.sparse.linalg.lsqr로 동시 solve
6. 잔차에서 양의 peak만 재검출
```

- 장점: 위치 고정이라 검출 문제 우회, 전역 동시 solve로 블렌드 처리
- 단점: 정확한 좌표 필요 (HST 카탈로그 or deep stack)

### 6-3. 장기 (데이터 수준)
- **AO 관측**: FWHM < 1" → 코어 분리 가능
- **HST 데이터**: 이미 ACS Survey of Globular Clusters 존재
- **Lucky imaging**: 순간 좋은 시상 프레임 선택

---

## 7. 현재 권장 워크플로우

```
step1~7 (aperture) → CMD 확인 → 문제 없으면 완료

PSF 적용이 필요한 경우:
  - CMD scatter가 crowding으로 크게 나올 때
  - 특정 magnitude에서 systematic bias 확인될 때

PSF 적용 시:
  - fit_engine = allstar, group max = 1
  - build_mode = epsf (per-frame)
  - max_iter = 3
```

---

## 8. 코드 위치

| 기능 | 파일 | 함수/클래스 |
|---|---|---|
| Newton single step | step8_psf_photometry.py | `_allstar_newton_one` |
| Newton group step | step8_psf_photometry.py | `_allstar_newton_group` |
| Group 빌드 | step8_psf_photometry.py | `_build_groups` |
| ALLSTAR outer loop | step8_psf_photometry.py | `_allstar_fit` |
| PSF evaluator | step8_psf_photometry.py | `_make_psf_evaluator` |
| Moffat 빌드 | step8_psf_photometry.py | `_build_moffat_psf` |
