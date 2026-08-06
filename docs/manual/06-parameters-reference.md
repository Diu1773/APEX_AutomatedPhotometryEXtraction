# 6. 파라미터 레퍼런스 (`parameters.toml`)

[← 분석 도구](05-tools.md) · [매뉴얼 목차](index.md) · [다음: 문제 해결 →](07-troubleshooting.md)

APEX의 모든 수치 설정은 작업 폴더의 **`parameters.toml`** 한 파일에 모여 있습니다. 대부분 각
단계의 **`Parameters`** 버튼에서 GUI로 바꿀 수 있고, 저장하면 이 파일에도 반영됩니다.

> **단위 규칙:** 키 이름의 접미사가 단위입니다 — `_px`(픽셀), `_arcsec`(초각), `_adu`, `_s`(초), `_deg`(도), `_mm`/`_um`.
> **편집 대상은 `parameters.toml`** 입니다. `parameters.example.toml`은 기본값 원본이니 건드리지 마세요.

아래는 자주 만지는 섹션부터 정리한 레퍼런스입니다. (각 단계에서 어떤 키를 쓰는지는 해당
단계 문서의 "파라미터" 표도 함께 보세요.)

---

## `[io]` — 입출력·경로

| 키 | 의미 | 기본값 |
| --- | --- | --- |
| `data_dir` | FITS 입력 폴더 | `data/example` |
| `result_dir` | 결과 출력 루트 | `<data_dir>/result` |
| `cache_dir` | 캐시 폴더 | `cache` |
| `filename_prefix` | 스캔 시 파일명 접두어 필터 | `""` |
| `night_gap_hours` | (LC) 밤을 가르는 JD 간격(시간) | 8.0 |
| `night_parse_mode` | 밤 파싱 방식(regex 등) | `regex` |
| `airmass_formula` | 공기질량 공식 | `Kasten & Young (1989)` |
| `airmass_update_header` / `_mode` | 헤더 AIRMASS 재기록 정책 | `false` / `overwrite` |

## `[target]` — 대상

| 키 | 의미 | 기본값 |
| --- | --- | --- |
| `name` | 대상 이름 | `M13` |
| `ra_deg` / `dec_deg` | 대상 중심 좌표(도) — WCS 힌트·카탈로그 기준 | 250.42 / 36.46 |

## `[instrument]` — 장비

| 키 | 의미 | 기본값 |
| --- | --- | --- |
| `telescope_focal_mm` | 망원경 초점거리(mm) | 3947.0 |
| `camera_pixel_um` | 픽셀 크기(µm) | 3.76 |
| `binning` | 비닝 | 2 |
| `gain_e_per_adu` | 게인(e⁻/ADU) — 노이즈 모델 | 0.68 |
| `rdnoise_e` | 리드노이즈(e⁻) | 2.5 |
| `saturation_adu` | 포화값(ADU) | 65000 |
| `datamin_adu` / `datamax_adu` | 측광 유효 하한/상한 | 0.1 / 55000 |
| `zp_initial` | 미리보기 영점 | 25.0 |
| `noise_use_fits_header` | 헤더 게인/리드노이즈 우선 사용 | `false` |

## `[crop]` — 크롭 (헤드리스)

| 키 | 의미 | 기본값 |
| --- | --- | --- |
| `enable` | 크롭 사용 여부 | `false` |
| `x0/y0/x1/y1` | 자를 픽셀 사각형 | (주석) |

> GUI로 그린 `step2_crop/crop_rect.json`이 이 값보다 우선합니다.

## `[parallel]` — 병렬·재실행

| 키 | 의미 | 기본값 |
| --- | --- | --- |
| `mode` | 워커 선택 정책 | `auto` |
| `max_workers` | 0이면 자동 결정 | 0 |
| `resume_mode` | 호환되는 기존 결과 재사용 | `true` |
| `force_redetect` / `force_rephot` | 검출/측광 강제 재계산 | `false` |

---

## `[detection]` — 소스 검출 (Step 4)

| 키 | 의미 | 기본값 |
| --- | --- | --- |
| `engine` | 검출 엔진(sep/segm/dao/peak) | `sep` |
| `sigma` | 검출 임계 시그마(기준) | 3.2 |
| `sigma_by_filter` | 필터별 시그마(자동 채움) | `{}` |
| `minarea_pix` | 최소 면적(px) | 3 |
| `keep_max` | 프레임당 최대 검출 수 | 6000 |
| `mode` | 검출 모드 프리셋 | `normal` |

**`[detection.deblend]`**(디블렌딩): `enable`=true, `nthresh`=64, `contrast`=0.004 …
**`[detection.peak]`**(피크 보조, 기본 off), **`[detection.dao]`**(DAO 정제, 기본 off).

## `[background]` — 배경 추정

| 키 | 의미 | 기본값 |
| --- | --- | --- |
| `enable` / `in_detect` | 배경 추정 사용 / 검출 시 사용 | `true` / `true` |
| `box` | 배경 박스(px) | 64 |
| `method` | 추정 방식 | `median` |
| `downsample` | 다운샘플 | 4 |

## `[fwhm]` — FWHM·미리보기 QC

| 키 | 의미 | 기본값 |
| --- | --- | --- |
| `guess_arcsec` | FWHM 시드(초각) | 2.5 |
| `px_min/px_max` | FWHM 허용 픽셀 범위 | 3.0 / 10.0 |
| `arcsec_min/max` | FWHM 허용 초각 범위 | 1.0 / 4.0 |
| `elong_max` | 최대 신장도 | 1.3 |

## `[qc]` — 프레임 품질 게이트

| 키 | 의미 | 기본값 |
| --- | --- | --- |
| `gate_enable` | QC 게이트 사용 | `true` |
| `sky_sigma_max_e` | sky σ 상한(e⁻) | 25.0 |
| `keep_positions_if_fail` | 실패해도 위치 유지 | `true` |

---

## `[photometry]` — 강제 조리개 측광 (Step 7)

| 키 | 의미 | 기본값 |
| --- | --- | --- |
| `mode` | 측광 모드 | `apcorr` |
| `recenter` | 검출 소스로 재중심 | `true` |
| `min_snr_for_mag` | 등급 계산 최소 SNR | 3.0 |
| `use_qc_pass_only` | QC 통과 프레임만 | `false` |

**`[photometry.scales]`** — FWHM 대비 배율: `aperture_scale`=1.0, `annulus_scale`=4.0, `dannulus_scale`=2.0, `center_cbox_scale`=1.5.
**`[photometry.radii]`** — 최소 반경·클리핑: `min_r_ap_px`=4.0, `sigma_clip`=3.0, `max_iter`=5.
**`[photometry.apcorr]`** — 조리개 보정: `apply`=true, `min_n`=12, `min_snr`=40.0, `max_sources`=250, `optimize_scales`=true.

---

## `[wcs]` / `[wcs_refine]` / `[wcs_qc]` — WCS 솔빙 (Step 5)

**`[wcs]`** 주요 키:
| 키 | 의미 | 기본값 |
| --- | --- | --- |
| `do_plate_solve` | 플레이트 솔빙 수행 | `true` |
| `astap_exe` | ASTAP 실행파일 경로 | `C:/Program Files/astap/astap_cli.exe` |
| `astap_database` | ASTAP 별 DB | `D80` |
| `astap_search_radius_deg` | 탐색 반경(도) | 8.0 |
| `astnet_local_enable` | astrometry.net 로컬 사용 | `true` |
| `astnet_local_use_wsl` | WSL 사용 | `true` |
| `refine_enable` | CRPIX/SIP 정제 | `true` |
| `max_workers` | WCS 병렬 워커 | 4 |
| `require_qc_pass` | Step 4 QC 통과 프레임만 | `true` |

**`[wcs_refine]`**: `enable`=true, `max_match`=600, `match_r_fwhm`=2.0, `min_match`=25, `max_sep_arcsec`=2.5.
**`[wcs_qc]`**: `match_radius_arcsec`=2.5, `min_match_n`=50, `max_rms_px`=2.5, `min_inlier_rate`=0.5.

## `[gaia]` / `[simbad]` — 카탈로그 조회

| 키 | 의미 | 기본값 |
| --- | --- | --- |
| `mag_max` | Gaia 등급 상한 | 25.0 |
| `wcs_mag_max` | WCS용 등급 상한 | 18.0 |
| `match_tol_arcsec` | 매칭 허용(초각) | 2.0 |
| `timeout_s` / `retry` | 타임아웃 / 재시도 | 30.0 / 2 |
| `pmem_method` | 멤버십 방법 | `gmm3d` |
| `pmem_ruwe_max` | RUWE 상한 | 2.0 |
| `simbad.timeout_s` | SIMBAD 타임아웃 | 20.0 |

---

## `[master]` / `[idmatch]` / `[refbuild]` — 마스터 카탈로그 (Step 6)

**`[master]`**: `n_master`=1000, `preserve_ids`=true, `keep_max`=12000, `filter_keep`="r".
**`[idmatch]`**: `tol_px`=1.0, `mode`="crowded", `two_pass_enable`=true, `transform_mode`="similarity".
**`[refbuild]`**: `sat_drop_pct`=20, `elong_drop_pct`=20, `per_date`=true, `ref_cat_min_sources`=50, `wcs_match_radius_arcsec`=2.0, `wcs_min_match_rate`=0.2.

---

## CMD 전용 — `[cmd]` / `[psf]` / `[isochrone]`

**`[cmd]`**(영점·표준화, Step 10):
| 키 | 의미 | 기본값 |
| --- | --- | --- |
| `snr_calib_min` | CMD 표시·보정 SNR 하한 | 50.0 |
| `apply_extinction` | k·X 소광 적용 | `false` |
| `extinction_mode` | absorb/two_step | `absorb` |
| `frame_zp_min_n` | 프레임 영점 최소 기준성 | 5 |
| `membership_mode` | 멤버십 컷 | `normal` |

`[cmd.zp]`: `clip_sigma`=3.0, `fit_iters`=5, `slope_absmax`=1.0. `[cmd.color]`: `clip_sigma`=3.0, `slope_absmax`=2.0.

**`[psf]`**(PSF 측광, Step 8): `epsf_oversampling`=2, `epsf_size_px`=25, `n_stars_max`=30, `max_iter`=3, `redetect_sigma`=4.5, `parallel_workers`=6, `build_mode`="epsf".

**`[isochrone]`**(아이소크론, Step 12):
| 키 | 의미 | 기본값 |
| --- | --- | --- |
| `file_path` | 아이소크론 파일 경로 | `""` |
| `age_init` | log Age 초기값 | 9.7 |
| `mh_init` | [Fe/H] 초기값 | -0.1 |
| `eg_r_init` | E(color) 초기값 | 0.0033 |
| `dm_init` | 거리지수 초기값 | 9.46 |
| `fit_fraction` | 피팅에 쓸 비율 | 0.6 |

---

## LC 전용 — `[light_curve.*]` / `[site]`

**`[light_curve.color_index_by_filter]`** / **`[light_curve.color_term_by_filter]`** — 필터별 색지수·색항(디트렌드 ΔC에 사용).

**`[site]`**(관측지 — 공기질량·BJD 계산):
| 키 | 의미 | 기본값 |
| --- | --- | --- |
| `lat_deg` / `lon_deg` | 위도/경도(도) | 36.607 / 127.360 |
| `alt_m` | 고도(m) | 81.0 |
| `tz_offset_hours` | 표준시 오프셋 | 9.0 |

---

## 도구 전용 — `[extinction_fit]` / `[tools.*]`

**`[extinction_fit]`**(소광 피팅 도구): `order`=1, `min_points`=5, `clip_sigma`=3.0, `fit_iters`=5, `use_color_terms`=false.

**`[tools.iraf.params]`** / **`[tools.iraf.filters.*]`** — IRAF/DAOPHOT 도구 파라미터(DATAPARS/FINDPARS/CENTERPARS/FITSKYPARS/PHOTPARS). 픽셀스케일 `pix_scale`=0.392, 게인 `epadu`=0.68 등.

**`[tools.qa_report]`** — QA 리포트: `min_n_frames`=5, `min_snr`=20.0, `max_chi2_nu`=5.0, `exclude_saturated`=true.

**`[tools.iraf_compare]`** — IRAF 비교: `match_tol_px`=1.5, `cutout_scale`=2.5.

---

## 그 밖의 섹션 (보통 손대지 않음)

- `[alignment]` — 전역 정렬 기준 필터/인덱스
- `[overlay]` — ID 라벨·시프트 벡터 표시
- `[hud5x]` — 5배 HUD 측광 미리보기 배율
- `[transform]` — src→ref 변환 저장
- `[parameters]` — 재중심 한계(`max_recenter_shift`=2.0), 센터링 이상치(`centroid_outlier_px`=1.0) 등
- `[cross_frame]` — RANSAC 프레임 간 매칭
- `[ui]` — 로그 길이·진행바·캔버스 픽셀

---

## 전체 기본값을 한 번에 보고 싶다면

저장소의 **`parameters.example.toml`** 이 모든 키의 기본값과 단위를 담은 원본입니다.
생성 문서 **`docs/parameter-inventory.md`** 에는 스키마 기반 자동 인벤토리가 있습니다.

```powershell
# 인벤토리 재생성(개발자용)
python scripts\inventory_foundation.py --write-docs
```

---

[← 분석 도구](05-tools.md) · [매뉴얼 목차](index.md) · [다음: 문제 해결 →](07-troubleshooting.md)
