# 설정이 코드에 닿는가 — 두 방향 감사

2026-08-16. 계기는 하늘 클리핑 설정이 Step 7 에 안 닿던 건(같은 모양 세 번째)을
고친 뒤 사용자가 *"다른놈들도 문제없는지 체크해 다시"* 라고 지시한 것이다.
로더를 읽어서가 아니라 **실제 워크스페이스를 걸어서** 재측정했다.

## 설정이 죽는 세 자리

| 이름 | 무슨 일이 일어나나 | 증상 |
|---|---|---|
| 매핑 없음 | 키가 `parameter_map` 에 없어 `raw` 에도 못 들어간다 | 파일에 적어도 로더가 통째로 무시 |
| **중간 유실** | `raw` 에는 들어오는데 `SimpleNamespace(...)` 호출이 그 이름을 인자로 안 적어 버려진다 | `P` 에 속성이 아예 없어 모든 `getattr(P, 이름, 기본값)` 이 기본값을 집는다 |
| 아무도 안 읽음 | `P` 까지 도착했는데 어느 모듈도 안 읽는다 | 값은 살아 있고 쓰이지 않는다 |

세 경우 모두 **실행은 성공하고 로그는 아무 말도 안 한다.** 정확도를 실제로
망친 건 가운데 것이다.

측정: `python -X utf8 -m apex.config.config_audit <apex_config.json> [--mode lc]`

## 이번에 고친 것

### 1. 기기 광학 제원이 P 까지 못 갔다 — 유일하게 수치가 틀릴 뻔한 건

`instrument.telescope_focal_mm` 과 `camera_pixel_um` 은 매핑돼 있는데
생성자가 안 받아서 `P` 에 없었다. `InstrumentConfig` 는 이 둘을
`getattr(params.P, "telescope_focal_mm", 3947.0)` 로 읽으므로 **어느 기기든
3947 mm / 3.76 µm — CDK500 + C3-61000 — 로 고정**됐다.

| 워크스페이스 | 설정 초점거리 | 설정 픽셀 | 올바른 픽셀스케일 | 고치기 전 실제 | 배율 |
|---|---|---|---|---|---|
| phase3 다섯 성단 | 3947 mm | 3.76 µm | 0.3930 "/px | 0.3930 | 1.00 (우연히 일치) |
| psf_crossinstrument/qhy600 | 1048 mm | 3.76 µm | 0.7400 "/px | 0.1965 | **0.27** |
| psf_crossinstrument/sinistro | 7953 mm | 15.0 µm | 0.3890 "/px | 0.1965 | **0.51** |
| psf_crossinstrument/m67_ubv | 7949 mm | 3.76 µm | 0.2335 "/px | 0.1965 | 0.84 |

**저장된 산출물은 무사하다.** 헤드리스 경로는 `InstrumentConfig` 를 아예 만들지
않고 `match.pixel_scale_arcsec` 를 파일에서 그대로 읽는데, 교차기기 설정들은
거기에 올바른 값(0.74 · 0.389 · 0.2335)을 직접 적어뒀다. `InstrumentConfig` 는
`apex/gui/main_window.py` 두 곳에서만 생성되고, 거기서 **올바른 값을 기본값으로
덮어썼다**. 즉 GUI 로 교차기기 워크스페이스를 열었을 때만 틀렸다.

고친 뒤 네 워크스페이스 모두 `InstrumentConfig` 의 픽셀스케일이 설정값과 일치한다.

### 2. 구경 배수가 양쪽에서 끊겨 있었다

Step 7 은 `forced_r_ap_scale` / `forced_ref_ap_scale` 을 읽었고 그 이름은
아무도 설정하지 않았다. 동시에 설정의 `photometry.apcorr.small_scale` /
`.large_scale` 은 `raw` 까지 왔다가 버려졌다. 양쪽 다 동작했고 서로 만난 적이
없다. 디스크의 설정 파일 50 개가 **어떤 실행도 쓴 적 없는 구경 반지름**을
싣고 있었다.

배선하면서 **산출물이 안 바뀌게** 맞췄다.

- 코드 기본값과 50 개 설정 파일, 그리고 두 예시 템플릿(`parameters.example.json`
  이 정본이고 `.toml` 은 주석용 사본)을 모두 `small_scale = 0.8` ·
  `large_scale = 2.4` 로 통일했다 — 지금까지의 모든 실행이 실제로 쓴 반지름이다.
  예시 JSON 을 같이 안 고치면 **새로 만드는 워크스페이스만 다른 구경**을 갖는다.
- 파일이 주장하던 1.0 / 3.0 은 `parameters.example.toml` 에서 복사된 값이고,
  같은 블록의 `small_scale_min = 0.8` · `large_scale_min = 2.4` 가 보여주듯
  코드는 자기 범위의 **바닥**에 앉아 있었다.
- 검증: `validation/apcorr/reproduce_apcorr.py` 로 M13 6 프레임 — 배선 후에도
  저장된 Step 7 결과와 최대 차이 **3.4e-7**.

**1.0 / 3.0 으로 올리는 건 별개의 결정이다** (아래 「사용자 판단 필요」).

CMD 와 LC 가 어긋나 있던 것도 같이 맞췄다. Step 7 은 두 모드 공용인데 LC 는
1.0 / 3.0 을 선언했고 CMD 는 필드 자체가 없었다 — 배선하는 순간 같은 엔진이
모드에 따라 다른 구경을 잴 뻔했다.

### 3. Step 7 캐시가 구경 변경을 못 알아챘다

캐시 무효화 서명 `_FORCED_SIGNATURE_PARAMS` 36 개 중 9 개가 없는 이름이라
`getattr(P, k, None)` 이 매번 `None` 을 해시에 넣었다. 즉 **구경이나 하늘
클리핑을 바꿔도 낡은 캐시가 재사용**됐다. 네 개(`forced_*_scale`,
`phot_sigma_clip`, `phot_max_iter`)를 실제 이름으로 바꾸고 서명 버전을
2 → 3 으로 올렸다. 남은 다섯은 `ref_cat_*` 로, 어떤 설정 파일도 쓰지 않아
낡을 수가 없다(아래 잔여 항목).

같은 고아 이름을 읽던 `apex/benchmark/runner.py`,
`apex/benchmark/photometry_crosscheck.py`, `apex/gui/tools/iraf_photometry.py`,
그리고 Step 7 GUI 툴팁도 함께 정렬했다.

### 4. 구경보정 on/off 스위치가 아무 일도 안 했다

`photometry.apcorr.apply` 는 `P` 까지 도착하는데 **어느 모듈도 읽지 않았다**.
Step 7 은 `flux_corr = flux_arr * apcorr` 를 무조건 실행했으므로
`apply = false` 로 꺼도 보정이 그대로 적용됐다. 50 개 설정 파일이 전부
`true` 라서 배선해도 산출물은 안 바뀐다 — 확인하고 이었다. 캐시 서명에도
넣어 껐다 켤 때 재계산되게 했다.

## 세 번째 범주 — P 까지 갔는데 아무도 안 읽는 것

M13 설정이 실제로 적은 매핑 키 364 개 중 **66 개**가 여기 해당한다.

| 구역 | 개수 | 예 |
|---|---|---|
| `idmatch` | 19 | `adaptive_retry_threshold` · `force` · `loose_radius_arcsec` · `match_r_fwhm` |
| `overlay` | 10 | `label_fontsize` · `max_labels` · `shift_min_px` (그림 주석 표시) |
| `photometry` | 8 | `apcorr.apply`(위에서 배선함) · `apcorr.scale_min/max/step` · `mode` |
| `gaia` | 6 | `pmem_method` · `pmem_ruwe_max` · `derived_enable` (구성원 확률) |
| `cross_frame` | 4 | `ransac_tol_px` · `ransac_max_iter` · `match_tol_px` |
| 나머지 | 19 | `alignment.global_align` · `background.method` · `qc.gate_enable` · `psf.save_model_image` 등 |

`idmatch` 와 `overlay` 는 GUI 편의 기능의 흔적으로 보이고, `gaia.pmem_*` 와
`cross_frame.ransac_*` 는 기능이 있는지부터 확인이 필요하다. 이 범주는 유실이
아니라 **미구현 또는 폐기**이므로 하나씩 판정해야 하며 이번에 손대지 않았다.

## 남은 것 — 못박아 뒀고, 늘어나면 테스트가 깨진다

`tests/test_config_settings_reach_code.py` 가 맵의 **모든 키를 넣은 설정**을
만들어 사라지는 이름을 전부 대조한다. 새 이름이 끼면 실패한다.

### CMD — 27 개

| 묶음 | 이름 | 판정 |
|---|---|---|
| 구경보정 배수 **탐색기** | `apcorr_optimize_scales` · `apcorr_small_scale_min/max` · `apcorr_large_scale_min/max` · `apcorr_max_pairs` · `apcorr_min_gap_fwhm` | 기능이 구현된 적 없다. 범위 안을 탐색하는 코드가 없고 고정 한 쌍만 쓴다. 지우거나 구현하거나 — 결정 필요 |
| 마스터 카탈로그 | `N_master` · `master_keep_max` · `master_filter_keep` · `master_flux_quantile` · `master_iso_min_sep_pix` | 코드는 같은 개념을 `ref_cat_max_sources` 라는 **다른 이름**으로 읽는다. 이름 통일 필요 |
| Step 10 색항 | `color_clip_sigma` · `color_fit_iters` · `color_slope_absmax` | 색항 적합 파라미터. 배선하면 step10 결과가 움직일 수 있어 확인 필요 |
| 프레임 QC 게이트 | `gate_nsrc_min` · `gate_sky_sigma_max_e` · `keep_positions_if_qc_fail` | 게이트 임계값. 배선 시 통과/탈락 프레임이 바뀔 수 있다 |
| 낡은 중복 | `ra_deg` · `dec_deg` | 무해. `target.ra_deg → target_ra_deg` 가 실제 경로이고 그건 도착한다 |
| 기타 흔적 | `clip_min_adu` · `clip_max_adu` · `cmd_max_sources` · `night_gap_hours` · `save_src2ref_tforms` · `ui_canvas_px` · `wcs_match_radius_arcsec` | 개별 확인 필요 |

### LC — 47 개

CMD 와 겹치는 것 외에 LC 만의 큰 덩어리 둘:

- **`wcs_qc_*` 열 개 전부** (`min_match_n` · `max_rms_px` · `max_p99_px` ·
  `min_inlier_rate` · `require_wcs_ok` 등) 와 `idmatch_wcs_qc_*` 여섯 개.
  WCS 품질 게이트 임계값이 통째로 안 닿는다.
- **PSF 관련 여섯**: `psf_profile_error_frac` · `psf_interp_order` ·
  `psf_forced_position_lock` · `psf_final_pass_max_iter` ·
  `psf_grouper_budget_frac` · `psf_grouper_budget_cap`.
  `psf_profile_error_frac` 은 2026-08-16 에 CMD 쪽에서 기본값을 정한 그 설정이다.

전체 목록은 테스트 파일의 `LC_DROPPED` 에 있다.

### 매핑 자체가 없는 키

M13 설정 기준 193 개. 대부분은 GUI 편의 항목이거나 `hud5x.*` 처럼 `P._raw` 로만
읽는 것들이라 정상이다. 이 감사에서는 분류하지 않았다.

## 사용자 판단 필요

**구경 반지름을 0.8 / 2.4 로 둘 것인가, 1.0 / 3.0 으로 올릴 것인가.**

- 지금까지의 모든 측정값은 0.8 / 2.4 다. 파일이 주장하던 1.0 / 3.0 은 쓰인 적 없다.
- 물리적으로는 둘 다 방어 가능하다. `r_ap = 0.8·FWHM` 은 가우시안 기준 포함광량
  83 %, `1.0·FWHM` 은 더 흔한 선택이다. `r_ref` 는 2.4 보다 3.0 이 날개를
  더 담아 표준에 가깝다.
- **올리면 다섯 성단 전부 Step 7 → 8 → 10 재처리**가 필요하고, 오늘 밤 측정한
  구경보정 대조 수치가 전부 다시 움직인다.
- 지금은 산출물을 지키는 쪽(0.8 / 2.4)으로 고정해 뒀다.
