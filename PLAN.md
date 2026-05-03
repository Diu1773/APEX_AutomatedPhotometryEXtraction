# APEX Pipeline Restructure + Forced Aperture Photometry

## Context

현재 문제: M13 같은 crowded field에서 IDMatch match rate가 ~50%로 측광값 절반이 버려짐.
근본 원인: detection-first 구조 — detection 안 된 별은 아예 측광 안 됨.
해결: master catalog 위치를 프레임에 투영해 강제 측광(forced aperture photometry) 추가.

코드 확인으로 발견한 추가 문제:
- CMD Step 5(Aperture Phot)가 WCS(Step 7) 이전에 실행됨 → ID 모르는 채 측광
- IDMatch가 transform을 메모리에만 계산하고 저장 안 함 → Forced Phot에서 재사용 불가
- RefBuild에 crowding_flag 없음

---

## 새 구조

### CMD (12단계)
| UI# | 이름 | 변경 | 디렉토리 |
|-----|------|------|----------|
| 1 | File Selection | 동일 | step1_file_selection/ |
| 2 | Crop | 동일 | step2_crop/ |
| 3 | Sky QC | 동일 | step3_sky_preview/ |
| 4 | Detection | 동일 | step4_detection/ |
| 5 | WCS | 이동 (기존 7→5) | step6_wcs/ ← 디렉토리명 유지 |
| 6 | MasterBuild | 확장 (기존 RefBuild) | step7_refbuild/ ← 디렉토리명 유지 |
| 7 | **Forced Aperture Phot** | **신규** | step_forced_phot/ |
| 8 | PSF Phot | 이동 (기존 6→8), skip 유지 | cmd_psf/ |
| 9 | Master ID Editor | 번호 조정 (9→9) | cmd_selection/ |
| 10 | Zeropoint | 번호 조정 | cmd_zeropoint/ |
| 11 | CMD Plot | 번호 조정 | cmd_plot/ |
| 12 | Isochrone | 번호 조정 | cmd_isochrone/ |

### LC (11단계)
| UI# | 이름 | 변경 | 디렉토리 |
|-----|------|------|----------|
| 1–4 | File Selection~Detection | 동일 | |
| 5 | WCS | 이동 (기존 6→5) | step6_wcs/ |
| 6 | MasterBuild | 확장 | step7_refbuild/ |
| 7 | **Forced Aperture Phot** | **신규** | step_forced_phot/ |
| 8 | Target Selection | 번호 조정 | lc_selection/ |
| 9 | LC Builder | 번호 조정 | (기존 step10) |
| 10 | Detrend | 번호 조정 | |
| 11 | Period Analysis | 번호 조정 | |

> step5_aperture/ 와 step8_idmatch/ 디렉토리명은 유지 (기존 데이터 호환)
> 해당 step은 main flow에서 제거하지만 파일 코드는 보존

---

## 파일별 변경 사항

### 신규 파일
- `apex/gui/workflow/step_forced_aperture_phot.py`
  - ForcedPhotWorker(QThread) + ForcedPhotWindow(StepWindowBase)
  - step_names key: `"forced_aperture_phot"`

### 수정 파일

**`apex/gui/main_window.py`**
- CMD step_names: Aperture Phot 제거, WCS↔PSF 순서 변경, Forced Phot 삽입
- LC step_names: Aperture Phot 제거, WCS 이동, Forced Phot 삽입
- `_open_step_window()` dispatch 인덱스 전면 갱신

**`apex/utils/step_paths.py`**
- 추가: `STEP_FORCED_PHOT_DIRNAME = "step_forced_phot"`
- 추가: `def step_forced_phot_dir(result_dir) -> Path`

**`apex/gui/workflow/step7_ref_build.py`** (MasterBuild 확장)
- 추가: 레퍼런스 프레임 기반 per-frame transform 계산 및 저장
  ```python
  # step7_refbuild/transform_{fname}.json
  {"affine_matrix": [...], "residual_rms_px": float, "n_anchors": int}
  ```
- 추가: master catalog에 `crowding_flag` 컬럼
  ```python
  master_df["crowding_flag"] = master_df["neighbor_dist_px"] < 2.5 * fwhm_ref
  ```

**`apex/gui/workflow/cmd/step6_psf_photometry.py`** (입력 경로 변경)
- 기존: `step5_aperture_dir(result_dir) / "photometry_index.csv"`
- 변경: `step_forced_phot_dir(result_dir) / "photometry_index.csv"`
- PSF worker가 aperture photometry를 초기 flux 추정으로 읽는 모든 경로 동일하게 변경

**`apex/gui/workflow/cmd/step11_zeropoint_calibration.py`**
- photometry_index 탐색 우선순위 변경:
  1. `step_forced_phot_dir / "photometry_index.csv"` (신규 우선)
  2. `step6_psf_dir / "photometry_index.csv"` (PSF 있으면)
  3. `step5_aperture_dir / "photometry_index.csv"` (fallback)
- idmatch 읽기 → forced phot의 master_id 컬럼으로 대체
  - `idmatch_{fname}.csv`의 source_id 조인 로직 → `photometry_{fname}.tsv`의 master_id 직접 사용

**`apex/gui/workflow/cmd/step10_master_id_editor.py`**
- `step8_idmatch_dir / "idmatch_{fname}.csv"` 읽기 → `step_forced_phot_dir / "photometry_{fname}.tsv"` 읽기
- master_id 컬럼 기반으로 편집 UI 유지

**`apex/gui/workflow/lc/step10_lightcurve_builder.py`**
- `step5_photometry_dir / "photometry_index.csv"` → `step_forced_phot_dir / "photometry_index.csv"`
- `photometry_{fname}.tsv` 읽기 경로 동일하게 변경

### 완전 삭제 (레거시 없이)
- `apex/gui/workflow/step5_aperture_photometry.py`
- `apex/gui/workflow/step5_aperture_worker.py`
- `apex/gui/workflow/step8_star_id_matching.py`

> project_state.json backward compat 불필요 (기존 데이터 재사용 안 함)

---

## Forced Aperture Phot 구현 상세

### 입력
- `step7_refbuild/ref_catalog_{filter}.tsv` — master 위치 (ra_deg, dec_deg, crowding_flag)
- `step7_refbuild/transform_{fname}.json` — per-frame affine transform
- `step6_wcs/` — WCS 헤더
- FITS 이미지 (cropped 우선)

### 처리 흐름 (ForcedPhotWorker.run)
```python
for fname in file_list:
    master_df = load_master_catalog(filter)
    wcs = load_wcs(fname)

    # 1. WCS로 master RA/Dec → 픽셀 좌표 초기값
    master_xy = wcs.all_world2pix(master_df[["ra_deg","dec_deg"]], 0)

    # 2. affine transform으로 정밀 보정 (transform_{fname}.json)
    master_xy_corrected = apply_affine(master_xy, transform)

    # 3. 검출된 별은 재중심 허용 (±max_recenter_shift px)
    #    미검출 별은 위치 고정
    detected_flag = recenter(master_xy_corrected, data)

    # 4. 소형 고정 aperture로 측광
    r_ap = params.forced_r_ap_px  # 별도 파라미터, 기존 r_ap와 다름
    phot = aperture_photometry(data, CircularAperture(positions, r_ap))

    # 5. apcorr: 밝은 detected 별에서만 계산 → 전체 적용
    apcorr = compute_apcorr(phot[detected_flag & bright])

    # 6. 출력 컬럼
    # master_id, source_id, x_pred, y_pred, x_fit, y_fit,
    # detected_flag, forced_flag, centroid_shift_px,
    # flux, flux_err, mag_inst, mag_err, snr,
    # sky, apcorr, bad_phot_flag
```

### 출력
- `step_forced_phot/photometry_{fname}.tsv` — per-frame
- `step_forced_phot/photometry_index.csv` — 프레임 요약
- `step_forced_phot/apcorr_summary.csv`

### parameters.toml 추가 항목
```toml
[forced_phot]
r_ap_scale = 0.8        # FWHM 배수 (기존 aperture_scale보다 작게)
max_recenter_shift = 2.0
min_snr_forced = 3.0    # forced별 upper limit 기준
```

---

## validate_step 로직

**ForcedPhotWindow.validate_step():**
```python
return (step_forced_phot_dir(result_dir) / "photometry_index.csv").exists()
```

**PSFPhotometryWindow.validate_step():** 변경 없음 (skip or index exists)

---

## 검증 방법

1. **syntax check**: `python -m compileall apex main.py`
2. **import smoke test**: `python scripts/smoke_steps.py`
3. **CMD 전체 실행**: M67 데이터로 Step 1→12 순차 실행, photometry_index.csv 생성 확인
4. **Forced phot 결과 확인**:
   - `step_forced_phot/photometry_index.csv` 존재
   - `detected_flag=True` 비율 ≈ 기존 match_rate
   - `forced_flag=True` 별도 포함 → 전체 별 수 증가 확인
5. **PSF skip 경로**: PSF skip 후 Zeropoint까지 진행 가능한지 확인
6. **LC 경로**: LC 모드 Step 1→7(Forced Phot)→9(LC Builder) 실행

---

## 주의사항

### Step 넘버링 실수 방지
- `main_window.py`의 `step_names` 리스트 인덱스와 `_open_step_window` dispatch 인덱스가 반드시 일치해야 함
- dispatch에서 매직 넘버(`step_index == 5`) 대신 상수 정의 권장:
  ```python
  _WCS_STEP = 4
  _MASTER_BUILD_STEP = 5
  _FORCED_PHOT_STEP = 6
  _PSF_STEP = 7  # CMD only
  ```
- step_names 리스트 변경 시 dispatch 인덱스 동시 업데이트 필수 (테스트로 검증)

### Tools 호환성 점검 필요
다음 tools가 step 경로를 직접 참조할 수 있음 — 구현 전 확인:
- `apex/gui/tools/` 하위 모든 파일에서 `step5_aperture_dir`, `step8_idmatch_dir` 참조 여부 grep
- 참조 발견 시 `step_forced_phot_dir`로 교체 또는 fallback 로직 추가

---

## 구현 순서

1. `step_paths.py` — `step_forced_phot_dir` 추가 (5분)
2. `main_window.py` — step_names + dispatch 갱신 (30분)
3. `step7_ref_build.py` — transform 저장 + crowding_flag (2시간)
4. `step_forced_aperture_phot.py` — 신규 worker + window (4시간)
5. `step6_psf_photometry.py` — 입력 경로 변경 (30분)
6. `step11_zeropoint_calibration.py` — 입력 경로 변경 (1시간)
7. `step10_master_id_editor.py` — 입력 경로 변경 (1시간)
8. `lc/step10_lightcurve_builder.py` — 입력 경로 변경 (1시간)
