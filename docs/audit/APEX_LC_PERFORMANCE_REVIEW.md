# APEX LC·도구 경로 성능 검토

검토일: 2026-08-07
범위: Step 8--11, `apex/analysis/light_curve/`, LC 관련 도구
목적: 사용자가 보고한 "데이터를 불러올 때부터 멈추는" 현상의 원인과 Bottleneck 적용 경계를 구분한다.

## 결론

현재 LC 후반부의 첫 병목은 `numpy.nan*` 통계 함수가 아니라 다음 네 가지가 겹치는
데이터 경로다.

1. 프레임마다 같은 source-id→표시 ID 카탈로그를 다시 찾고 읽는다.
2. 모든 프레임을 전체 DataFrame으로 메모리에 올린 뒤, target·comparison 별로 같은
   행을 반복해서 Boolean 필터링한다.
3. PSF 선택 시 한 프레임의 aperture 표와 PSF 표를 매번 각각 읽고 여러 열을 다시
   매핑한다.
4. Step 10의 플롯은 데이터가 조금만 바뀌어도 세 축을 지우고 전체 점을 다시 그린다.

따라서 우선순위는 **I/O·자료구조·플롯 캐시를 먼저 고치고, 그 다음 실제로 큰 배열인
calibration/ensemble/detrend 통계에만 Bottleneck을 적용**하는 것이다. `Bottleneck을
전체 코드에 적용`하는 것은 이 증상에 대한 정확한 처방이 아니다.

이 문서는 코드 검토 결과이지 속도 개선을 이미 입증한 결과가 아니다. 숫자 성능 주장은
고정 입력, cache warm/cold, worker 수, peak RSS를 기록한 별도 실험 뒤에만 허용한다.

## 관찰된 경로와 병목 후보

| 우선순위 | 위치 | 현재 동작 | 예상 영향 | Bottleneck 적합성 |
|---|---|---|---|---|
| P0 | `apex/utils/photometry_loader.py:92-129,132-230`; Step 9 `step9_lightcurve_builder.py:2365-2391` | `load_frame_photometry(..., sid_map=None)`이면 프레임마다 selection catalog를 glob하고 다시 읽어 source-id→ID map을 만든다. Step 9 preload 경로는 이 map을 전달하지 않는다. | 프레임 수에 비례하는 중복 파일 검색·CSV 파싱. 첫 LC build, 도구의 Step 7 재로드 모두 영향. | 아니오. 한 번 만든 map을 공유하는 캐시/인자 설계 문제. |
| P0 | Step 9 `step9_lightcurve_builder.py:2398-2541,2543-2871,2873-3070` | cache에 전체 DataFrame을 보관하지만 target과 각 comparison을 `df[df["ID"] == ...]` 또는 source-id 선택으로 매번 찾는다. | 대략 `프레임 × 별 수 × 행 수`의 Python/pandas 반복. ensemble·diff·check series에서 중복. | 거의 아니오. `set_index` 1회, compact frame×star matrix 또는 join이 우선. |
| P0 | Step 9 `step9_lightcurve_builder.py:335-401` | 별별 median을 계산할 때 프레임 DataFrame을 읽은 후 모든 star ID에 대해 같은 `df_ids == sid` mask를 다시 만든다. | median table/색지수 계산이 별 수와 프레임 수에 따라 불필요하게 증가. | 아니오. 한 번의 groupby/pivot 또는 indexed lookup이 우선. |
| P0 | `apex/analysis/light_curve/photometry_source_service.py:154-269` | PSF 경로가 프레임마다 aperture CSV/TSV를 읽은 뒤 PSF TSV를 다시 읽고 `frame_uid.map(psf[column])`을 여러 번 수행한다. | PSF를 선택한 첫 Step 8 preview/Step 9 build에서 파일 읽기와 복사가 2배 이상. | 아니오. merged per-frame cache와 필요한 열만 읽는 것이 우선. |
| P1 | Step 9 `step9_lightcurve_builder.py:3745-3838,4714-4785` | build 시작 시 각 comparison을 target으로 삼아 `_build_ensemble_series`를 반복해 QC를 계산한 뒤, 실제 target ensemble과 check-star ensemble을 다시 만든다. | 같은 frame cache를 사용해도 row 선택과 Python loop를 여러 번 반복. | 부분적. 먼저 공통 frame×star 배열을 만들고 QC/build/check가 공유해야 한다. |
| P1 | Step 9 `step9_lightcurve_builder.py:2365-2391` | 모든 파일을 `ThreadPoolExecutor`로 동시에 전체 DataFrame preload하고 캐시는 무제한이다. | 디스크가 빠르지 않으면 동시 CSV 파싱이 오히려 지연을 만들고, 큰 필드에서는 peak RSS가 커진다. | 아니오. bounded cache, column pruning, sequential read benchmark가 필요. |
| P1 | Step 8 `photometry_source_service.py:271-431` | 비교 후보 preview가 필터의 모든 프레임을 순회하여 전체 측광 행을 만들고, 첫 요청 시에만 필터 캐시를 채운다. | 첫 후보를 클릭할 때 긴 UI 대기. 이후 캐시는 유효하므로 cold-start 병목. | 아니오. background prefetch와 compact cache가 우선. |
| P1 | Step 10 `step10_detrend_merge.py:1750-1835,4091-4300` | raw CSV는 한 번 읽지만 copy/concat/보정 후 세 Matplotlib 축을 `clear()`하고 모든 점을 다시 그린다. | 데이터 행 수와 플롯 호출 횟수에 따라 화면 갱신이 느려짐. | 아니오. artist 재사용, 표시용 downsample/binning, `draw_idle` coalescing이 우선. |
| P1 | Step 11 `period_io_service.py:74-195` | 필터를 바꿀 때마다 light-curve CSV 전체를 `pd.read_csv`한 뒤 필터/ID를 선택한다. | 큰 combined CSV에서 필터 콤보 변경이 반복 I/O가 된다. | 아니오. fingerprint 기반 파일 cache와 `usecols`/필터별 materialized cache가 우선. |
| P1 | 도구 `extinction_fit.py:439-650,1461-1620`, `qa_report.py:173-260` | 두 도구 모두 프레임마다 `load_frame_photometry`를 직접 호출하고 전체 표를 복사·정규화한다. `sid_map` 공유가 없다. | LC와 같은 중복 map 읽기 및 전체 표 재구성. | 아니오. 공통 frame loader/cache를 공유해야 한다. |
| P1 | 도구 `variable_star.py:142-230,4056-4160` | 같은 workspace의 primary lightcurve CSV를 경로별로 전체 읽고 series 선택 배열을 여러 개 만든다. | 변수성 분석 창을 여는 순간 파일 수와 CSV 크기에 비례하는 지연. | 아니오. path fingerprint cache와 `usecols`가 우선. |
| P2 | `comparison_stability_service.py:101-213` | leave-one-out residual이 star마다 다른 comparison 열을 다시 구성하고 median을 계산한다. | comparison 수가 커지면 `stars × frames × comparisons`와 DataFrame 복사 증가. | 제한적. leave-one-out median 알고리즘/배열화가 우선이며 Bottleneck은 그 뒤의 큰 median에만 후보. |
| P2 | `global_ensemble.py:237-865` | sparse `lsqr`, pandas groupby/apply, weighted groupby가 중심이다. | 계산량은 있지만 병목은 sparse solver와 groupby일 가능성이 높다. | 제한적. `nan*`를 교체하기 전에 solver·groupby를 측정해야 한다. |
| P2 | `period_analysis_service.py` | Lomb--Scargle/BLS/PDM/Bootstrap은 Astropy·NumPy·SciPy 경로와 별도 worker를 사용한다. | period 계산 자체가 오래 걸릴 수 있으나 CSV 로드 지연과는 별개다. | 대체로 아니오. Bottleneck을 붙여도 LS/BLS kernel 시간은 줄지 않는다. |

## Bottleneck을 실제로 시험할 곳

다음 조건을 모두 만족하는 곳만 `fast_stats` 라우팅 후보로 삼는다.

* 입력이 프레임/별/화소 축을 가진 큰 수치 배열이다.
* 같은 reduction이 한 작업에서 반복되고, pandas groupby·CSV·Matplotlib보다 큰 비중을
  차지한다.
* NumPy와 Bottleneck의 NaN, axis, dtype, `ddof` 결과가 parity test로 고정되어 있다.

우선 후보는 다음 세 곳이다.

1. Step 0 calibration master stack와 활성화된 overscan 영역의 `nanmedian`, `nanmean`,
   `nanstd`, `nansum`. 현재 wrapper 경계와 가장 잘 맞고, 이미지 stack은 충분히 크다.
2. Step 9/10에서 먼저 compact matrix를 만든 뒤의 frame×star ensemble residual 통계.
   현재의 반복 DataFrame 선택을 그대로 둔 채 `np.nanmedian`만 바꾸는 것은 이득이 작다.
3. global ensemble/detrend에서 대규모 residual 행렬을 만든 뒤의 robust summary. sparse
   `lsqr`와 groupby 자체를 Bottleneck 성능으로 포장하면 안 된다.

반대로 다음은 Bottleneck 후보에서 제외한다.

* 10--20개 정도의 per-frame/per-star 배열(호출 overhead가 결과를 지배한다).
* `period_analysis_service`의 Lomb--Scargle/BLS/PDM/Bootstrap kernel.
* pandas `read_csv`, `concat`, `groupby`, `iterrows`, Matplotlib redraw.
* `benchmark/`의 NumPy 기준 구현. 독립 기준선이므로 production wrapper로 바꾸지 않는다.

## 권장 최적화 순서

### 1단계: 측정 계측부터 추가

각 경로에 다음 이벤트를 기록한다(과학 산출물에는 포함하지 않음).

* `load_index_ms`, `load_frame_ms`, `source_map_read_count`, `psf_table_read_count`
* `cache_hit/miss`, `frames_loaded`, `rows_loaded`, `selected_star_count`
* `build_ensemble_ms`, `build_check_ms`, `qc_ms`, `plot_ms`
* process peak RSS, worker 수, input fingerprint, Bottleneck 활성 여부

측정 단위는 Step 8 preview 첫 클릭, Step 9 첫 build, Step 9 comp 변경, Step 10 raw
로드, Step 11 필터 변경으로 고정한다. warm cache와 cold cache를 별도로 측정한다.

### 2단계: 자료 경로를 줄인다

* result directory/filter 단위로 source-id→ID map을 한 번만 만들고 Step 9·extinction·QA에
  전달한다. `load_lightcurve_frame_photometry`에도 optional `sid_map`을 전달할 수 있게
  하되, legacy frame의 source-id/ID fallback semantics는 유지한다.
* frame loader가 필요한 열만 읽거나, 한 번 만든 공통 schema를 cache한다. 원본 CSV/TSV를
  덮어쓰지 말고 cache key에 파일 size+mtime(가능하면 content fingerprint), source mode,
  선택 열을 포함한다.
* 전체 DataFrame cache와 별도로 target/comp/check에 필요한 compact matrix를 만든다.
  중복 source-id가 있을 때 기존의 `first`/`drop_duplicates` 정책을 그대로 보존하고,
  64-bit Gaia source ID를 `int64`로 유지한다.
* PSF 경로는 aperture와 PSF의 공통 `det_uid` join을 한 번만 수행한 뒤 preview/build가
  재사용한다. PSF QC flag와 provenance 열을 버리면 안 된다.

### 3단계: 화면 갱신을 줄인다

* Step 10/9의 기존 Matplotlib artists를 갱신하고, 축을 매번 `clear()`하지 않는다.
* overview는 점 수 상한 또는 시간 bin median/quantile envelope를 사용하고, 확대/선택
  때만 원자료를 그린다. downsample은 저장 CSV와 분석 배열에는 절대 적용하지 않는다.
* filter/date/phase 조작을 debounce하여 짧은 연속 이벤트가 한 번만 render되게 한다.

### 4단계: 마지막으로 Bottleneck을 시험한다

compact matrix와 동일 입력을 사용해 NumPy/Bottleneck을 각각 최소 5회 측정한다. 보고할
값은 median wall time, 변동성, peak RSS, 결과 parity이다. 배열 shape, dtype, axis, `ddof`,
NaN 비율, package versions를 함께 보존한다. 효과가 작은 per-frame 호출은 wrapper에
추가하지 않는다.

## 과학적 불변조건

최적화 후 다음은 byte-for-byte가 아니어도 허용오차와 행 선택 규칙까지 검증해야 한다.

* frame exclusion, filter normalization, night assignment, BJD/JD 선택
* source_id 64-bit 보존과 source-id→display-ID 매핑 우선순위
* 중복 source/frame의 기존 first/quality 선택 규칙
* comparison weight, missing/error 처리, check-star 독립성
* raw/corrected series의 row count와 finite mask
* period 분석 입력 배열과 선택된 correction mode

속도가 빨라졌다는 이유만으로 map cache나 matrix join이 과학적 동일성을 보장하지는
않는다. 위 조건을 고정한 회귀 테스트가 먼저다.

## 현재 판단

현재 코드에서 Bottleneck을 바로 더 심는 것은 우선순위가 아니다. 사용자가 느끼는 LC
초기 지연의 가장 유력한 순서는 **(1) 프레임별 source map 재읽기, (2) PSF의 이중 표
읽기, (3) target/comparison 반복 행 검색, (4) Step 10 전체 플롯 재작성**이다. 이 네
가지를 계측·캐시·인덱스화한 뒤에도 큰 reduction이 남을 때만 `fast_stats`에 추가한다.
