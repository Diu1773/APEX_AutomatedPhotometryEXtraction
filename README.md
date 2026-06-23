# APEX - Automated Photometry EXtraction

APEX는 FITS 관측 영상에서 측광 결과까지 이어지는 과정을 단계별로
수행하는 Python/PyQt5 데스크톱 애플리케이션입니다.

APEX is a Python/PyQt5 desktop application for step-by-step astronomical
photometry, from FITS observations to calibrated CMDs, light curves, and
period analysis.

> 현재 배포 대상은 Windows x64이며, 소스 실행에는 Python 3.10 이상이
> 필요합니다. GitHub Actions의 기준 환경은 Python 3.11입니다.

## 분석 모드

- **CMD mode (12 steps)**: source detection, WCS, forced/PSF photometry,
  zeropoint calibration, CMD plotting, and isochrone fitting.
- **LC mode (11 steps)**: forced photometry, target/comparison selection,
  differential light curves, detrending, multi-night merge, and period
  analysis.
- **Shared Steps 1-7**: file selection, crop, sky/QC preview, source
  detection, WCS solving, master catalog build, and forced aperture
  photometry.

## 주요 기능

- FITS 파일 및 필터별 프레임 관리
- sky, FWHM, saturation, elongation 기반 품질 확인
- SEP/DAO 기반 source detection과 master catalog 구성
- 내장 Python WCS solver, ASTAP, local astrometry.net 지원
- Gaia/SIMBAD 조회와 Gaia 기반 WCS refinement/QC
- forced aperture photometry와 선택적 PSF photometry
- CMD 영점·색항 보정 및 PARSEC isochrone fitting
- 비교성 QC, global ensemble, detrending, SysRem, multi-night merge
- Lomb-Scargle, PDM, BLS 및 bootstrap FAP 기반 주기 분석
- artificial-star benchmark, IRAF cross-check, CMD validation 도구

### 측광 검증

APEX 측광은 독립 소프트웨어와 직접 교차검증됩니다 — SExtractor(`sep`)·IRAF/DAOPHOT와
**~3 mmag**(robust MAD) 일치, Gaia DR3 합성측광 대비 **색의존 편향 없음**.
- [측광 교차검증](docs/validation_crosscheck.md) — primary 측정 검증 (`apex validate --suite crosscheck`)
- [실측 성단 재현](docs/validation_realdata.md) — 거리 ~4% (0.85–7.8 kpc)
- [이소크론 피팅 방법론](docs/isochrone_fitting_methodology.md) — §10: gri/BVR 축퇴와 u밴드 해법

## 빠른 시작

### 소스에서 실행

PowerShell:

```powershell
git clone https://github.com/Diu1773/APEX_AutomatedPhotometryEXtraction.git
cd APEX_AutomatedPhotometryEXtraction

python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

python main.py
```

모드를 직접 실행할 수도 있습니다.

```powershell
python apex\cmd\main.py
python apex\lightcurve\main.py
```

루트 launcher 또는 배포된 `APEX.exe`를 처음 실행하면
`parameters.example.toml`을 기반으로 로컬 `parameters.toml`이
생성됩니다. 원본 예제 파일은 기본값과 설정 스키마의 기준이므로,
사용자별 경로와 장비 값은 `parameters.toml`에서 변경하십시오.

설정 항목 설명은 [Configuration Guide](docs/configuration.md)를
참조하십시오.

### Windows 빌드 사용

GitHub Actions의 `Windows Build` 아티팩트에는 다음 파일이 포함됩니다.

- `setup-APEX-<version>.exe`: 사용자별 설치 프로그램
- `APEX-Portable-<version>-x64.zip`: 설치 없이 사용할 수 있는 portable
  bundle

설치본과 portable bundle 모두 Python 런타임과 필수 패키지를
포함합니다. 외부 WCS solver와 해당 catalog/index 데이터는 포함하지
않습니다.

## 명령줄 인터페이스 (CLI)

APEX는 GUI 없이 동작하는 헤드리스 CLI를 제공합니다. `pip` 설치 시 `apex`
명령이 등록되며, 소스 체크아웃에서는 `python -m apex`로 동일하게 실행할 수
있습니다.

```powershell
# 설치 (editable). 환경이 용도별로 분리되어 있습니다:
python -m pip install -e .            # headless core (서버/CI/배치용, GUI 제외)
python -m pip install -e ".[gui]"     # 데스크톱 GUI 포함
python -m pip install -e ".[dev]"     # 개발(GUI+test+build+docs)

# 런타임 진단: Python·의존성·외부 solver·(옵션)네트워크 점검
apex doctor
apex doctor --network

# 설정 파일 관리
apex config init        # parameters.example.toml -> parameters.toml
apex config path
apex config show

# 헤드리스 파이프라인 (공유 Step 1-7)
apex run --mode cmd --dry-run            # 실행 계획만 미리보기
apex run --mode cmd --steps 1            # Step 1(헤더 스캔·타깃 해석) 실행
apex run --mode cmd --steps 4-7 --force  # 일부 단계 강제 재실행
apex run --mode lc  --result-dir D:\out  # 출력 경로 오버라이드

# GUI 실행
apex gui              # 런처
apex gui --mode cmd   # 모드 직접 실행
```

`apex run`은 단계별 산출물 존재 여부로 완료를 판정해 재실행을 건너뛰고
(`--force`로 무시), 실행 결과를 `result_dir/pipeline_run.json` 매니페스트로
남깁니다. 현재 Step 1(scan)이 헤드리스로 동작하며, Step 2-7은 GUI 산출물을
인식·계획에 포함하고 점진적으로 헤드리스 이식 중입니다(권장 순서 4→7→3→6→5→2).

## 워크플로

### Shared Steps 1-7

| Step | Name | 역할 | 출력 폴더 |
| ---: | --- | --- | --- |
| 1 | File Selection | FITS scan, target/filter/frame 선택 | `step1_file_selection/` |
| 2 | Image Crop | 공통 crop 영역 적용 | `step2_crop/` |
| 3 | Sky Preview & QC | sky, FWHM, saturation, frame QC | `step3_sky_preview/` |
| 4 | Source Detection | frame별 source catalog와 품질 정보 생성 | `step4_detection/` |
| 5 | WCS Plate Solving | WCS solve, Gaia refinement, solver QC | `step5_wcs/` |
| 6 | Master Catalog Build | frame 간 source를 master ID로 통합 | `step6_refbuild/` |
| 7 | Forced Aperture Phot | master 위치 기반 frame별 aperture 측광 | `step7_forced_phot/` |

Step 1 can use valid RA/Dec values from the selected FITS headers when a
SIMBAD name lookup is unnecessary or unavailable. Click `Use Header RA/Dec`
to accept those coordinates and enable the next step.

Step 4 builds `detection.sigma_by_filter` from the active FITS `FILTER`
headers. Filter keys follow the canonical APEX convention, including
case-sensitive Johnson-Cousins `R`/`I` and Sloan `r`/`i`.

### CMD Steps 8-12

| Step | Name | 역할 | 출력 폴더 |
| ---: | --- | --- | --- |
| 8 | PSF Photometry | 혼잡장 선택적 PSF 측광 | `cmd_psf/` |
| 9 | Master ID Editor | source/ROI/membership 검토 | `cmd_selection/` |
| 10 | Zeropoint Calibration | frame zeropoint와 color term 보정 | `cmd_zeropoint/` |
| 11 | CMD Plot | calibrated color-magnitude diagram 생성 | `cmd_plot/` |
| 12 | Isochrone Model | PARSEC/BaSTI isochrone 탐색 및 fitting | `cmd_isochrone/` |

Step 12 reads the photometric system and magnitude-column names from the
isochrone header. The selected observation bands must exist in that file; for
example, Johnson `B-V` data require a Johnson-Cousins/Bessell file and cannot
be fitted with an SDSS `ugriz` file. A compact near-solar BaSTI Johnson grid
can be generated with `python tools/download_basti_isochrones.py`.

Step 10 uses the Pancino et al. (2022) Gaia EDR3 dwarf transformation for
Johnson `B`, preserving the main-sequence `B-V` color scale. Johnson `V`,
`R`, and `I` continue to use the Riello et al. (2021) transformations.

Step 12 enables an SNR 20 display and fitting threshold by default. Fainter
low-SNR points remain available by lowering or disabling the filter, but they
should not be used to infer the lower-main-sequence slope.

### LC Steps 8-11

| Step | Name | 역할 | 출력 폴더 |
| ---: | --- | --- | --- |
| 8 | Target/Comparison Selection | target, comparison, check star 선택 | `lc_selection/` |
| 9 | Light Curve Builder | differential light curve와 comparison QC | `lc_lightcurve/` |
| 10 | Detrend & Night Merge | detrending, ensemble correction, night merge | `lc_detrend/` |
| 11 | Period Analysis | LS/PDM/BLS, FAP, phase-folded curve | `lc_period/` |

상세한 데이터 흐름과 모듈 경계는
[Architecture](ARCHITECTURE.md)를 참조하십시오.

## WCS Solver

Step 5에는 세 가지 경로가 있습니다.

1. **Internal (Python)**

   Gaia catalog와 Step 4 source를 quad code로 매칭하고 RANSAC 검증 후
   TAN/SIP WCS를 적합합니다. frame header 좌표, Step 1 target 좌표,
   local-blind 순서로 hint를 재시도합니다. 외부 실행 파일은 필요
   없지만 Gaia catalog cache 또는 네트워크 접근이 필요합니다.
2. **ASTAP (Local)**

   Windows에서 사용할 수 있는 외부 solver입니다. ASTAP 실행 파일과
   영상 시야에 맞는 D80/D50 star database를 별도로 설치해야 합니다.
3. **Astrometry.net (Local)**

   `solve-field`와 시야에 맞는 index files가 필요합니다. Windows에서는
   일반적으로 WSL/Ubuntu를 통해 사용합니다.

외부 설치 자료:

- [ASTAP](https://www.hnsky.org/astap.htm)
- [ASTAP star databases](https://sourceforge.net/projects/astap-program/files/star_databases/)
- [Astrometry.net](https://astrometry.net/doc/readme.html)
- [Astrometry.net index files](https://data.astrometry.net/)

## 개발과 검증

```powershell
python -m compileall apex main.py scripts deploy
python -m pytest tests
python main.py --smoke
```

`--smoke`는 GUI를 열지 않고 배포에 필요한 주요 workflow, Gaia/SIMBAD,
Tools 모듈을 import합니다.

Windows release 전체 빌드:

```powershell
.\build.bat
```

필요 조건과 산출물은 [Deployment Guide](deploy/README.md)에 정리되어
있습니다.

과학 검증:

```powershell
python validation\cmd_step12_synthetic.py
python benchmark\run_benchmark.py --config benchmark\configs\baseline.toml
```

추가 검증 절차는 [Benchmark Guide](benchmark/README.md)와
[Validation Guide](validation/README.md)를 참조하십시오.

## 프로젝트 구조

```text
apex/
  analysis/       astrometry, CMD, light-curve, merge science services
  benchmark/      reusable benchmark and validation services
  cmd/            CMD mode entry point
  config/         TOML-backed CMD/LC parameter models
  core/           project state, file management, instrument model
  gui/            main window, workflow steps, tools, widgets
  lightcurve/     LC mode entry point
  resources/      SVG and packaged GUI assets
  utils/          paths, cache, I/O, astronomy, photometry helpers
benchmark/        benchmark CLI wrappers, configs, and run outputs
deploy/           PyInstaller/Inno Setup release pipeline
docs/             technical design and configuration documents
tests/            pytest suite
validation/       synthetic and real-data validation runners
```

## 문서

- **[사용자 매뉴얼 (User Manual)](docs/manual/README.md)** — 스크린샷과 함께 단계별로
  따라 하는 실사용 가이드(CMD/LC 전 단계 · 도구 · 파라미터 · 문제 해결)
- [Documentation Index](docs/README.md)
- [Architecture](ARCHITECTURE.md)
- [Configuration Guide](docs/configuration.md)
- [Deployment Guide](deploy/README.md)
- [Benchmark Guide](benchmark/README.md)
- [Validation Guide](validation/README.md)
- [Cache Manager Design](docs/cache-manager-design.md)
- [PSF Photometry Analysis](docs/psf-photometry-analysis.md)

## 상태와 제한

- 버전은 `deploy/version.txt`에서 관리합니다.
- Windows installer는 현재 사용자 영역에 설치되며 관리자 권한을
  요구하지 않습니다.
- ASTAP/astrometry.net 실행 파일과 catalog/index 데이터는 배포본에
  포함되지 않습니다.
- 일부 Gaia/SIMBAD 기능은 네트워크 상태와 원격 서비스 가용성에
  영향을 받습니다.
- 라이선스는 MIT입니다([LICENSE](LICENSE) 참조).
