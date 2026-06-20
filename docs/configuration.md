# APEX Configuration Guide

APEX는 TOML 설정을 사용합니다.

- `parameters.example.toml`: repository에 저장되는 기본값
- `parameters.toml`: 사용자 환경의 실제 runtime 설정
- `apex/config/parameters_cmd.py`: CMD mode mapping
- `apex/config/parameters_lc.py`: LC mode mapping

## First Launch

루트 launcher 또는 배포된 `APEX.exe`는 `parameters.toml`이 없으면
`parameters.example.toml`을 복사합니다. Mode-specific entry point를
직접 실행할 때는 먼저 root launcher를 한 번 실행하거나 runtime
파일을 준비하는 것이 좋습니다. 설정 변경은 runtime 파일에 적용하고,
기본값을 변경할 때만 example 파일과 mode mapping을 함께 수정하십시오.

경로는 TOML 파일 위치가 아니라 application의 working directory와
workflow code에서 해석될 수 있으므로, 재현 가능한 분석에는 절대
경로 또는 project-relative 경로를 일관되게 사용하는 것이 좋습니다.

## Core Sections

### `[io]`

- `data_dir`: FITS 입력 경로
- `result_dir`: workflow 출력 루트
- `cache_dir`: 공통 cache 경로
- `filename_prefix`: FITS scan prefix filter
- `night_*`: LC night grouping과 filename parsing
- `airmass_*`: airmass 계산 및 header update 정책

### `[target]`

- `name`: 대상 이름
- `ra_deg`, `dec_deg`: target center, WCS hint, catalog query 기준

Step 1의 SIMBAD resolution으로 갱신할 수 있습니다.

When the selected FITS files already contain usable sky coordinates, Step 1
also offers `Use Header RA/Dec`. This accepts the representative header
coordinates as `target.ra_deg` and `target.dec_deg`, so name resolution is not
required before continuing.

### `[instrument]`

- `telescope_focal_mm`
- `camera_pixel_um`
- `binning`
- `gain_e_per_adu`
- `rdnoise_e`
- `saturation_adu`
- `datamin_adu`, `datamax_adu`

gain과 read noise는 측광 오차 모델에 직접 영향을 줍니다.
`noise_use_fits_header`가 켜져 있으면 frame header 값을 우선 사용할 수
있습니다.

### `[parallel]`

- `mode`: worker 선택 정책
- `max_workers`: `0`이면 자동 결정
- `resume_mode`: 호환되는 기존 결과 재사용
- `force_redetect`, `force_rephot`: 해당 step 강제 재계산

### `[detection]`, `[background]`, `[fwhm]`

Step 3-4의 background estimation, detection threshold, deblending, FWHM/QC
범위를 제어합니다. Detection parameter를 바꾸면 Step 4 이후 결과를
재생성해야 합니다.

Step 4 materializes `detection.sigma_by_filter` from the active FITS `FILTER`
headers before detection. Missing per-filter values inherit
`detection.sigma`. Canonical keys are case-sensitive where the photometric
systems differ: Johnson-Cousins uses `R` and `I`, while Sloan uses `r` and
`i`. Stale keys from a different dataset are not applied to the active files.

### `[photometry]`

Step 7 aperture photometry와 reference-star selection을 제어합니다.

- `[photometry.scales]`: FWHM 대비 aperture/annulus scale
- `[photometry.radii]`: 최소 pixel radius와 clipping
- `[photometry.apcorr]`: aperture correction 및 scale optimization

### `[wcs]`, `[wcs_refine]`, `[wcs_qc]`

외부 solver 경로, Gaia refinement, WCS acceptance threshold를 제어합니다.
Internal solver의 interactive parameter는 Step 5 dialog에서 관리되는
값도 있으므로, UI와 TOML 중 실제 ownership을 코드에서 확인해야
합니다.

Internal solver는 다음 정보에 민감합니다.

- 정확한 pixel scale 추정
- header 또는 Step 1 target center
- 충분한 Step 4 source
- field를 덮는 Gaia catalog
- quad code/RANSAC tolerance

### `[gaia]`, `[simbad]`

catalog query limit, retry, magnitude threshold, membership와 derived field를
제어합니다. 큰 field에서 `gaia.wcs_mag_max`를 너무 어둡게 설정하면
query와 matching 비용이 크게 증가할 수 있습니다.

### `[master]`, `[idmatch]`, `[match]`

Step 6 master catalog 크기와 frame 간 source matching/QC를 제어합니다.
match tolerance 변경은 source ID와 downstream photometry를 바꿀 수
있습니다.

### `[cmd]`, `[psf]`, `[isochrone]`

CMD calibration, PSF fitting, isochrone input/default fit 값을 제어합니다.
PSF는 선택 단계이며 crowded field가 아닌 경우 Step 7 aperture 결과를
직접 사용할 수 있습니다.

Step 12 derives available bands from the isochrone file header. PARSEC CMD
tables and APEX-normalized BaSTI tables are supported. The configured
isochrone file must use the same photometric system as the calibrated CMD.
Johnson-Cousins `B`, `V`, `R`, `I` selections are rejected for an SDSS
`ugriz` file instead of reusing columns with similar positions.

Run `python tools/download_basti_isochrones.py` to download and normalize the
near-solar BaSTI Johnson-Cousins grid under `isochrone/BaSTI/johnson`.

Gaia-based Johnson `B` calibration uses the Pancino et al. (2022) dwarf
relation. This is intentionally different from the legacy approximate `B`
conversion, which compressed `B-V` and changed the lower-main-sequence slope.

The Step 12 CMD display and auto-fit start at SNR 20. Below this level,
flux-boosting and unequal B/V completeness can bias faint `B-V` colors blue.

### `[light_curve.*]`, `[extinction_fit]`

LC color index/color term과 extinction fitting 옵션을 정의합니다.
Step 9-11의 일부 interactive 선택은 result directory의 JSON/CSV state로
저장되며 TOML만으로 완전히 표현되지 않습니다.

## Change Impact

| 변경 범주 | 일반적으로 다시 실행할 시작 단계 |
| --- | ---: |
| input frame/filter/crop | Step 1-2 |
| background/FWHM/detection | Step 3-4 |
| WCS solver/refinement/QC | Step 5 |
| master matching | Step 6 |
| aperture/noise/apcorr | Step 7 |
| PSF model | CMD Step 8 |
| calibration/color terms | CMD Step 10 |
| comparison stars | LC Step 8-9 |
| detrend/ensemble/SysRem | LC Step 10 |
| period grid/method/FAP | LC Step 11 |

Cache reuse는 output 존재 여부만이 아니라 입력 signature와 parameter
compatibility를 확인해야 합니다. 자세한 정책은
[Cache Manager Design](cache-manager-design.md)을 참조하십시오.

## Developer Rules

새 설정을 추가할 때:

1. `parameters.example.toml`에 default와 단위를 명확히 추가한다.
2. 사용하는 mode의 parameter map에 runtime attribute를 연결한다.
3. CMD와 LC가 공유하면 두 map을 모두 갱신한다.
4. 숫자 단위는 key 이름에 포함한다: `_px`, `_arcsec`, `_adu`, `_s`.
5. parameter loading과 default behavior를 pytest로 검증한다.
6. 필요하면 다음 명령으로 inventory를 갱신한다.

```powershell
python scripts\inventory_foundation.py --write-docs
```

`docs/parameter-inventory.md`는 생성 문서이므로 수동으로 schema 설명을
추가하지 않습니다.
