# 실 GUI 구동 계측 — CMD 12스텝 · LC 완주

**2026-07-31 야간** · 요청: 「LC까지 다 돌려서 실 GUI 테스팅, 얼마나 걸리는지까지」

계측 도구: `scripts/gui_step_timing.py` — 실제 창 클래스를 스텝 순서대로 열어
로드 시간·표 모양·창 크기·예외를 기록한다. 사람이 클릭하는 것과 같은 코드경로다
(오프스크린 렌더만 다르다).

## 1. CMD 12스텝 — 전부 열린다

M13 재처리 결과(15프레임, 혼잡장) 기준. **12/12 예외 없음.**

| 스텝 | 창 | 열기 | 로드 | 표 |
|---|---|---:|---:|---|
| 0 File Selection | | 8.7s | — | 15×10 |
| 1 Image Crop | | 8.1s | — | — |
| 2 Sky Preview & QC | | 3.6s | — | — |
| 3 Source Detection | | 8.6s | — | 15×7 |
| 4 WCS Plate Solving | | 1.3s | — | 15×6 |
| 5 Master Catalog Build | | 1.2s | — | 0×0 |
| 6 Forced Aperture Phot | | 0.2s | — | 0×12 |
| **7 PSF Photometry** | | 2.0s | **17.8s** | 15×9 |
| 8 Master ID Editor | | **10.1s** | — | 1500×4 |
| 9 Zeropoint Calibration | | 0.3s | — | — |
| 10 CMD Plot | | 0.6s | — | — |
| 11 Isochrone Model | | 2.9s | 0.3s | — |

메인창 부팅 4.1~5.0s. Step 8 표는 9열이고 `Forced %` 가 채워진다(어제 추가분 확인).

## 2. 찾아서 고친 것

### 2.1 PSF 창을 열 때마다 비교표를 처음부터 다시 만든다 — 10.9초 (`8cf63d4`)

Step 8 창 로드 17.8초를 분해하니 `build_ap_psf_comparison` 이 **10.93초**였다.
`_refresh_qc()` 가 창을 열 때마다 `_cmp_merged_df` 를 버리고 다시 병합한다.

**어제 `mag_ap` 폴백을 넣기 전에는 이 병합이 늘 "All magnitudes are NaN" 으로
즉시 끝나서 비용이 안 보였다.** 버그를 고치자 비용이 드러난 것이다 — 내가 만든
회귀다.

Step 8 산출물이 그대로면 병합 결과도 같으므로 `psf_ap_vs_psf.csv` 가
`photometry_index.csv` 보다 새로울 때만 캐시로 읽게 했다. 헤드리스 export 도
같은 캐시(+ meta json)를 남긴다. 확인: 프로파일 상위에서 사라졌다.

프레임 수에 비례하므로 **프레임이 많을수록 이득이 커진다**(M13 은 15장이다).

### 2.2 `param_file` 을 지정해도 이전 프로젝트가 이긴다 (`6c859e3`)

`MainWindowWorkflow(param_file=M13)` 로 열었는데 창이 NGC6811 의 `result_dir` 을
쓰고 있었다. `_bootstrap_file_selection_state()` 가
`apex/.state/<mode>/project_state.json` — **모드당 하나뿐**이고 마지막에 GUI 로 연
프로젝트가 남아 있는 파일 — 로 `params.P` 를 덮어쓴다.

**평소 GUI 실행에서 겪는 문제는 아니다.** 앱에는 파라미터 파일을 바꾸는 UI 가
없어서 사용자는 늘 같은 `parameters.toml` 을 쓰고 Step 1 에서 폴더를 고른다.
그 복원은 의도된 편의다. 그래서 기본 동작은 그대로 두고 **파일을 명시해서 연
경우에만** 지정한 쪽을 따르게 했다 — 검증 스크립트가 조용히 다른 자료를 재는 것을
막는다.

> 함정: 이 상태 파일은 레포 루트가 아니라 **`apex/.state/`** 안에 있다
> (`project_root = Path(__file__).parent.parent` 가 패키지 디렉터리다).

## 3. LC 모드 완주 — YZ Bootis

| | |
|---|---|
| 자료 | `E:\observe_raw_Analysis\YZbootis\YZbootis_20250430` |
| 기기 | Moravian C3-61000, bin2, 4788×3194 (15.3 MP), 0.393″/px |
| 구성 | Light **283장** (g 97 · r 94 · i 92, 30 s) + dark 30 s ×10 + flat g/r/i ×5 |
| 보정 | bias `E:\bias` · dark 노출 정확히 일치 · flat 자체 폴더 |
| 출력 | `E:\APEX_validation\reprocess\YZBoo_20250430\` |

원본은 `observed_Analysis` 가 아니라 **`observe_raw_Analysis`** 에 있다
(전자에는 결과 폴더만 남고 raw 는 0장).

### 실측 소요

| 단계 | 시간 |
|---|---|
| Step 0 보정 | **장당 37.5 s** (283장이면 약 177분) |
| Step 1–7 | (측정 예정) |
| LC 8–11 | (측정 예정) |

장당 37.5초는 15.3 MP + 우주선 제거(astroscrappy)가 함께 도는 값이다.

283장 전체는 3시간이라 **g 필터 97장으로 범위를 좁혀 완주를 먼저 확보**한다.
광곡선은 필터별로 만들므로 파이프라인 완주 검증에는 한 필터로 충분하다.

### 시도했다가 되돌린 것 — Step 0 보정 스레드 병렬화

`step0_detector_calibration.py` 의 라이트 보정 루프는 **완전 순차**다
(`ThreadPool` 이 없다). 16코어를 놀리고 있길래 마스터 매칭을 먼저 끝내고
보정만 `ThreadPoolExecutor` 4워커로 돌려 봤다.

**효과가 없었다 — 장당 37.5 s → 39.3 s (배속 1.0x).** 작업 프로세스의 CPU 가
**61%**(0.6코어)라 애초에 CPU 병목이 아니었다. 다만 이 측정은 pytest 와
ASTERION 파이썬 8개가 같이 도는 경합 환경에서 잰 값이라 **병렬화가 무용하다는
증명은 아니다**. 검증되지 않은 복잡도는 남기지 않기로 하고 되돌렸다.

> 다시 볼 때의 조건: 다른 작업이 없는 상태에서 순차/병렬을 각각 20장씩 재고,
> 그때도 CPU 가 100% 를 못 넘으면 I/O 쪽을 파야 한다(E 드라이브는 관측 자료용
> 대용량 디스크다).

### LC 스텝 창은 자료 없이도 전부 열린다

Step 0~11 **12/12 예외 없음**(빈 `result_dir` 기준). LC 고유 스텝은
8 Target/Comparison Selection · 9 Light Curve Builder · 10 Detrend & Night Merge ·
11 Period Analysis 이고, 각각 워커가 있다
(`_ComparisonAutomationWorker` · `_LightCurveTaskWorker` · `_DetrendTaskWorker` ·
`PeriodAnalysisWorker`). 워커가 창(owner)을 받으므로 헤드리스로 돌리려면 창을
띄우고 실행 경로를 부르는 방식이 맞다.

**완주 판정 기준: 문헌 주기 P = 0.104092 일을 재현하는가.**

## 4. 오라클

수정 3건(비교표 캐시 · `param_file` 존중 · 병렬화 되돌림) 뒤
**625 passed, 0 failed** (16분 59초). 회귀 없음.

## 5. 아직 안 된 것

- LC Step 8–11 은 **헤드리스 러너가 없다.** `apex.pipeline.registry` 에 등록된
  것은 두 모드 모두 Step 1–7 뿐이고, CMD 8/10/12 만 개별 스크립트가 있다.
  LC 는 GUI 워커로 돌린다.
- 기존 LC 결과(`RESULT_YZbootis_*`)는 **옛 스텝 번호 체계**(`step10_lightcurve`)라
  지금 코드로 그대로 열 수 없다. 그래서 처음부터 돌린다.
