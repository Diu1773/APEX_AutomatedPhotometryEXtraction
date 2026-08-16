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

## 구조 수리 — 설정 정의를 한 군데로 (2026-08-16)

위 넷을 고친 뒤에도 **끊길 자리 자체는 그대로**였다. 같은 설정이 네 군데에
따로 적혀 있었기 때문이다.

| 어디 | 무엇을 적고 있었나 | 결과 |
|---|---|---|
| `parameter_map.py` | 점표기 경로 → 평평한 이름 (370행) | 정본 |
| `parameters_cmd.py` 생성자 | 이름 + 타입 + 기본값 (452줄) | 손으로 동기화 |
| `parameters_cmd.py:75` 로컬 맵 | 경로 → 이름 (295행) | **즉시 덮어써져 버려짐** |
| `schema.py` | pydantic 모델 36개 (1762줄) | **아무도 안 씀** |

버려지던 로컬 사본은 화석이었다 — 정본보다 CMD 는 81행, LC 는 123행이 빠져
있었다. 편집해도 아무 일이 안 일어나는 파일이 두 개 있었던 셈이다.

**고친 방식.** 맵의 행이 설정의 전부를 담게 했다.

```python
(('instrument', 'telescope_focal_mm'), 'telescope_focal_mm', 'float', 3947.0)
#  어디서 오는가              뭐라 부르는가        무슨 타입   없으면 뭘 쓰나
```

`build_settings(raw, TOML_KEY_MAP)` 가 이 행들을 그대로 값으로 바꾸므로,
**새 설정은 한 줄이고 반만 추가되는 일이 불가능하다.** 337행이 이 형태로
바뀌었고, 두 로더의 손으로 쓴 인자 559줄과 버려지던 사본 565줄이 사라졌다.

두 부분만 있는 행은 맵이 알지만 직접 세우지 않는 설정이다 — 정규화가
필요하거나(`str(...).strip().lower()`), CMD·LC 가 다른 기본값을 쓰거나(5개),
다른 값에서 계산되거나, `P._raw` 로만 읽거나, **아직 코드에 안 닿는 것**이다.
마지막 부류가 아래 목록이고, 이제는 맵을 보면 눈에 띈다.

**증명.** 실제 워크스페이스 51 개 × 두 모드의 **속성 35,751 개를 이행 전후로
전부 비교해 차이 0**. 맵만 바꾼 시점, 로더를 바꾼 시점, 주석을 정리한 시점,
화석을 지운 시점 네 번 모두 확인했다.

**`schema.py` 정정.** 앞서 "아무도 안 쓴다"고 적었는데 `apex/`·`scripts/` 만
보고 내린 판정이라 틀렸다 — `tests/test_parameters_foundation.py:318` 이
`Parameters.from_toml` 을 부른다. 제품 코드 경로에는 없지만 테스트 하나가
붙들고 있다. 또 그건 타입 변환이 아니라 **검증**(범위·열거형) 계층이라 키 맵이
대체하지 못한다. 지우지 않았다 — 채택할지 버릴지는 사람 판단.

## 배선한 12 개 — 동작은 그대로 (2026-08-16)

설정 파일이 적어 둔 값이 **코드가 이미 쓰던 상수와 같은** 것들만 골라
배선했다. 그래서 오늘 산출물은 안 바뀌고 설정만 실재하게 된다.

| 모드 | 설정 | 값 |
|---|---|---|
| CMD·LC | `io.night_gap_hours` | 8.0 |
| LC | `detection.keep_max` · `background.downsample` | 6000 · 4 |
| LC | `wcs_qc.clip_sigma` · `require_wcs_ok` · `max_rms_px` | 3.0 · true · 2.5 |
| LC | `wcs_qc.min_inlier_rate` · `max_edge_ratio` · `max_center_offset_arcsec` | 0.5 · 0.0 · 0.0 |
| LC | `refbuild.master_union` · `union_min_frames` | true · 1 |

확인: 51 개 워크스페이스 재측정에서 **기존 값 변화 0 · 새로 도착 12 개**,
도착값은 전부 코드 상수와 동일.

## LC 가 설정을 무시하는 자리 — 배선하면 결과가 바뀐다

CMD 는 읽고 LC 는 안 읽는 설정들이 있다. Step 7·8 과 WCS 해는 두 모드가
**같은 코드**를 쓰므로, 어느 모드로 열었느냐에 따라 다른 숫자로 돈다.

| 설정 | 설정 파일 요구 | LC 가 실제로 쓰는 값 |
|---|---|---|
| `wcs_qc.min_match_n` | 50 | **20** |
| `wcs_qc.min_match_rate` | 0.05 | **0.20** |
| `wcs_qc.max_p99_px` | 5.2 | **5.0** |
| `wcs_qc.match_radius_arcsec` | 2.5 | **2.0** |
| `idmatch.wcs_qc_*` 7 개 | 값 있음 | 전부 미도달 |
| `psf.profile_error_frac` 등 6 개 | 값 있음 | 전부 미도달 |

`wcs_qc` 는 **프레임 합격/불합격 기준**이라 배선하면 LC 광도곡선에 들어가는
프레임이 달라진다. 배관 수리가 아니라 과학 결정이므로 손대지 않고 못박아 뒀다
(`LC_DROPPED`). AE UMa·YZ Boo 재현 작업과 직결된다.

## 세 번째 표면 — 파라미터 창 (2026-08-16)

사용자: *"GUI 파라매터 창들이랑 연결되는 거 맞지? 근본적인 해결책도?"*

체인은 **위젯이 `getattr(P, 이름, 상수)` 를 보여주고 → 사용자가 고치면 `P` 에
쓰고 → `save_toml` 이 키 맵을 돌며 파일에 적는** 구조다. 어느 고리가 끊겨도
창은 "Parameters saved." 를 띄우고 아무것도 안 쓴다. 재현:

```
불러온 뒤 P 에 있나 : False
파일의 값           : 3.0
save_toml 결과      : True   ← 창은 "Parameters saved." 를 띄운다
저장 후 파일의 값   : 3.0    ← 9.9 로 바꿨는데 그대로
다시 불러온 값      : <속성 없음>
```

**모드별로 따로 봐야 한다.** 처음엔 `CMD 맵 ∪ LC 맵` 으로 검사해 30 개로 봤는데,
그건 너무 느슨해서 가장 큰 건을 놓쳤다 — **Step 6(마스터 카탈로그)은 공용
단계인데 설정 17 개가 전부 LC 전용 맵에만** 있었다. CMD 로 열면 창이 보여주고
하나도 저장하지 못한다. 모드별 재측정: **CMD 27 개 · LC 3 개.**

| 창 | 개수 | 조치 |
|---|---|---|
| Step 6 마스터 카탈로그 | 17 | 공용 맵으로 이동 |
| 소광 보정 도구 | 17 + 2 | 새 행 추가 · 공용 이동 |
| 대기질량 도구 (`night_parse_*` 포함) | 7 | 공용 이동 — 설정 파일 `io` 절에 값이 있는데 CMD 가 안 읽고 있었다 |
| Step 4·5 (`detect_mode` 등) | 3 | 공용 이동 + LC 로더에 같은 줄 |
| 그 외 창 | 6 | 새 행 추가 |

전부 **코드가 이미 쓰던 상수를 행 기본값으로** 넣었다. 51 개 워크스페이스
재측정에서 **기존 값 변화 0 · 사라진 것 0 · 새로 도착 84 개.**

### 근본 해결 — 창을 맵 행에서 생성한다

같은 설정이 창에서 한 번 더 선언되던 것을 없앴다. 행에 표시 정보(라벨·범위·
단위)를 얹고 `specs_from_map()` 이 창을 만든다. 창은 **보여줄 행의 이름과
순서만** 선언한다.

```python
_STEP6_ATTRS = ("ref_select_sat_pct", "ref_per_date", ...)
_STEP6_SPECS = specs_from_map(_STEP6_ATTRS)
```

맵에 행이 없는 이름을 창에 올리면 `specs_from_map` 이 **즉시 예외**를 던진다 —
저장 못 하는 위젯을 만들 방법이 사라진다. Step 6 로 시범 적용했고, 생성된
19 개 항목이 손으로 쓰던 것과 **완전히 동일**함을 확인했다(라벨·타입·범위·
스텝·소수점·단위·기본값 전부).

회귀: `tests/test_gui_settings_persist.py` 가 모드별로 검사한다.

### 남은 결정 하나

`extfit_min_points` — LC 맵은 파일 값 5 를 주고 CMD 코드는 10 을 쓴다. 배선하면
CMD 의 소광 적합이 받아들이는 점 수가 달라진다. 예외 목록에 근거와 함께 남겼다.

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
