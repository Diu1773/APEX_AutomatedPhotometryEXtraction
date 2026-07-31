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

### 실측 소요 — g 97장 완주

| 단계 | 시간 | 결과 |
|---|---:|---|
| **Step 0 보정** | **50분** (31.0 s/장) | 97장, 마스터 bias/dark 30 s/flat g |
| Step 1 scan | 16.5 s | 97 files, 타깃 좌표 확인 |
| Step 2 crop | 0.2 s | 설정 없어 skip |
| Step 3 sky QC | 0.0 s | 헤드리스 no-op |
| Step 4 detect | 196 s | 10,145 소스, 중앙 FWHM 5.71 px |
| **Step 5 WCS** | **1796 s (30분)** | 97/97 solved (astnet), **WCS-QC pass 42** |
| Step 6 refbuild | 16.8 s | 마스터 **203** 소스 |
| Step 7 forcedphot | 221 s | detected 8,760 / forced 10,931 · 투명도 QC PASS 85 |
| Step 8 선택 | 20 s | 타깃 115 + 비교성 5 |
| Step 9 광곡선 | 35 s | `lightcurve_ID115_raw.csv` 97점 |
| Step 11 주기 | ~10분 | `lc_period/` 6개 산출 |

**Step 1–7 합계 37.4분.** 그중 **WCS 가 80%** 다 — 여기가 유일한 병목이다.

### 완주 판정 — 문헌 주기를 재현한다

| 방법 | 최적 주기(일) | 문헌 대비 | FAP |
|---|---:|---:|---:|
| **Lomb–Scargle** | **0.104792** | **+0.67%** | 1.6e-39 |
| PDM | 0.113360 | +8.90% | — |

문헌 YZ Boo P = 0.104092 일. **LS 는 0.67% 로 재현**했고 FAP 가 1e-39 이라 우연이
아니다. PDM 이 8.9% 벗어난 것은 관측이 **3.05시간 = 1.2주기**뿐이라 위상 빈이
고르게 안 채워졌기 때문이다(PDM 은 빈 기반이라 짧은 관측에 취약하다).

### 주의 — Step 9 의 `mag` 는 기기등급이다

`lightcurve_ID115_raw.csv` 의 `mag` 는 차등이 아니라 **기기등급**이다. 실제로
타깃과 비교성이 **함께** 1 mag 가까이 움직인다 — 별의 변광이 아니라 투명도 변화다.

| | p2p | σ |
|---|---:|---:|
| 타깃 115 | 1.329 | 0.256 |
| 비교성 121 | 0.988 | 0.205 |
| 비교성 36 | 0.978 | 0.205 |
| 비교성 194 | 1.019 | 0.186 |

직접 차등(타깃 − 비교성 3개 평균)을 내면 진폭이 **1.329 → 0.562 mag** 로 줄어
공통 성분이 걷힌다. SX Phe 형 YZ Boo 로 그럴듯한 값이다. 타깃은 포화가 아니다
(`is_saturated` 전부 0, SNR ~900, mag_err 0.0012).

`variable_analysis_bundle.json` 은 릴리스를 **BLOCKED** 로 막는다 — 이유가
「체크성 미선택 · 비교성 안정도 메타 없음」이다. 품질 게이트가 제대로 작동한 것이다.

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

## 4.5 WCS 병목의 진짜 원인 — 엔진을 잘못 골랐다 (2026-07-31 낮)

「WCS 가 30분으로 전체의 80%」는 **엔진 선택 실수**였다. 헤드리스 엔진 우선순위는

```
wcs_engine 명시  →  astap_enable  →  astnet_local_enable  →  internal (기본)
```

인데(`resolve_wcs_engine`), M13 템플릿을 복사할 때 딸려 온
`astnet_local_enable = true` 때문에 astrometry.net(WSL)이 선택됐다. `[wcs] engine`
을 비워 두면 자체 엔진이 기본인데 그게 밀린 것이다.

`engine = "internal"` 로 바꿔 같은 97장을 다시 풀었다.

| Step 5 | 시간 | QC 통과 | match_rate | rms_px |
|---|---:|---:|---:|---:|
| astnet | 1796 s | **42**/97 | 0.971 | — |
| **internal** | **398 s** | **97**/97 | **0.991** | 0.895 |
| internal (캐시 적중) | **53 s** | 97/97 | | |

측광도 따라 좋아진다 — 같은 Step 7 에서

| | detected | forced | 투명도 QC |
|---|---:|---:|---|
| astnet | 8,760 | 10,931 | PASS 85 / REVIEW 11 / **FAIL 1** |
| internal | **10,105** | **8,810** | PASS 90 / REVIEW 7 / **FAIL 0** |

WCS 가 정확해지니 강제 측광이 줄고 실제 검출이 늘었다.

### 그 과정에서 찾은 버그 — 자체 엔진이 `frame_wcs_qc.csv` 를 안 썼다

`frame_wcs_qc.csv` 를 쓰는 곳은 `WcsWorkerBase`(ASTAP)와
`AstrometryNetWorkerBase`(astnet) 둘뿐이고, **내부 엔진 경로
(`_write_internal_summary_csv`)는 요약 CSV 만 쓰고 QC 파일을 남기지 않았다.**
Step 6/7 은 이 파일로 프레임을 거르므로(`qc_utils.filter_files_by_wcs_qc`),
자체 엔진이 97/97 을 통과시킨 뒤에도 **이전 솔버가 남긴 낡은 판정**
(astnet 시절 `{False:55, True:42}`)이 그대로 읽혔다. 판정에 필요한 열은 이미
요약 rows 에 다 있어서, 같은 스키마로 함께 내보내도록 고쳤다.

> **다만 이 자료에서는 하류 결과가 바뀌지 않았다.** 낡은 QC(42통과)로 돌린
> Step 6/7 과 새 QC(97통과)로 돌린 결과가 detected 10,105 / forced 8,810 으로
> 동일하다. 구조적 결함이라 고치는 게 맞지만, **측광이 달라진다고 말할 근거는
> 이 실험에 없다.**

### 틀렸던 가설

「마스터가 203개로 적은 것은 WCS-QC 가 42 라서」— **아니다.** QC 를 97 로 올려도
**195개 그대로**다. 마스터 수는 QC 통과 프레임 수가 아니라 필드의 실제 별 수와
refbuild 조건이 정한다.

## 4.6 엔진 교차검증 — 자료 4개 · 검출기 2종

YZ Boo 한 필드로 「internal 이 훨씬 좋다」고 결론 내면 위험해서, 같은 프레임을
두 엔진으로 각각 풀어 비교했다(`wcs_engine_crosscheck.py`).

| 대상 | 검출기 | 성격 | 엔진 | 시간 | **rms_px** | QC 통과 | n_match |
|---|---|---|---|---:|---:|---:|---:|
| M13 | Moravian C3-61000 | 구상(혼잡) | internal | 33.8 s | **0.427** | **15/15** | — |
| | | | astnet | **13.4 s** | 0.605 | 14/15 | 938 |
| NGC6811 | Moravian C3-61000 | 산개 | internal | **63.5 s** | **0.345** | **21/21** | 881 |
| | | | astnet | 87.2 s | 1.511 | 19/21 | 881 |
| M67 | **LCO QHY600** | 산개 | internal | **25.3 s** | **0.191** | **10/10** | 850 |
| | | | astnet | 53.4 s | 0.253 | 10/10 | 850 |
| YZ Boo | Moravian C3-61000 | 변광성장 | internal | **398 s** | **0.895** | **97/97** | 104 |
| | | | astnet | 1796 s | 2.599 | 42/97 | 102 |

**정확도는 4/4 전부 internal 이 이긴다**(1.3~4.4배). QC 통과율도 4/4 전부 100%.
결정적인 것은 **`n_match` 가 두 엔진에서 같다**는 점이다(NGC6811 881, M67 850) —
같은 별을 같은 수로 매칭하는데 **해의 정밀도만 다르다**. 별을 못 찾아서가 아니다.

속도는 3/4 우세이고 **M13 만 astnet 이 빠르다**(13 s vs 34 s). 다만 시간은 캐시
상태에 민감하다 — M13 internal 이 첫 실행 92 s → 재실행 34 s 였다.

> **검증 도구가 먼저 틀렸던 기록.** 처음 돌렸을 때 astnet 이 3개 자료에서 모두
> "4초 만에 실패"로 나왔다. 엔진 문제가 아니라 **내 스크립트가 두 번째 엔진을
> 설정할 때 `engine` 키를 중복으로 넣어 TOML 이 깨진 것**이었다. 그대로 믿었으면
> 「astnet 은 아예 안 돈다」는 결론을 낼 뻔했다. 지금은 설정 후 `tomllib` 로 다시
> 읽어 값이 맞는지 확인하고, 실패 시 오류 꼬리를 남긴다.

검증에 쓴 자료의 `step5_wcs` 는 전부 원래 상태로 되돌렸다.

## 5. LC Step 8 에서 걸린 것 — 토글과 지연 저장

헤드리스로 몰 때 두 번 헛돌았다. 둘 다 **실제 GUI 사용에서는 문제가 아니지만**
스크립트로 몰 때는 반드시 알아야 한다.

1. **타깃·비교성 지정이 둘 다 토글이다.** `set_target_selected()` 는 같은 별을
   다시 부르면 **해제**하고, `toggle_comparison_selected()` 도 마찬가지다. 이전
   실행이 남긴 `selection_g.json` 을 `load_selections()` 가 복원하므로, 같은
   스크립트를 두 번 돌리면 결과가 반대로 뒤집힌다. 실제로 target=None 이 되거나
   comparison=[] 이 되는 것을 번갈아 봤다. → **선택 파일을 지우고 한 번에** 지정.
2. **선택 저장이 `_queue_selection_save()` 로 지연된다.** 지정 직후 프로세스를
   끝내면 저장이 안 된다. 이벤트 루프를 몇 초 돌려야 타이머가 발화한다.
3. **`build_light_curve()` 는 `runtime_mode=True` 일 때만 동기 실행**된다
   (아니면 워커에 큐잉). 헤드리스에서는 이 플래그를 켜야 결과를 바로 받는다.
4. `auto_select_comparisons()` 는 워커가 12.9초 돌고도 비교성을 하나도 안
   붙였다. 원인 미확인 — **비교성은 수동 지정으로 우회**했다.

## 6. 아직 안 된 것

- LC Step 8–11 은 **헤드리스 러너가 없다.** `apex.pipeline.registry` 에 등록된
  것은 두 모드 모두 Step 1–7 뿐이고, CMD 8/10/12 만 개별 스크립트가 있다.
  LC 는 GUI 워커로 돌린다.
- 기존 LC 결과(`RESULT_YZbootis_*`)는 **옛 스텝 번호 체계**(`step10_lightcurve`)라
  지금 코드로 그대로 열 수 없다. 그래서 처음부터 돌린다.
