# LC 스텝 9 헤드리스 — 라이트커브가 창 밖으로 나왔다

2026-08-19. 커밋: 이 문서와 같은 작업.

## 무엇이 막고 있었나

LC 는 스텝 7 에서 멈춰 있었다. 이유는 계산이 Qt 에 묶여서가 아니다 —
`apex/analysis/light_curve/` 의 14 개 모듈 약 7,000 줄은 처음부터 Qt 를 안 쓴다.
막고 있던 건 **라이트커브를 만드는 코드 자체가 창의 메서드였다**는 것이다.

`LightCurveBuilderWindow._build_ensemble_series` 271 줄. 프레임마다 대상 별과
비교성 앙상블을 찾아 차등등급을 내는, LC 가 실제로 하는 계산 전부가
`QMainWindow` 자식 클래스 안에 있었다. 창에 의존하는 부분은 하나도 위젯이
아니었다 — 파라미터 모델, 캐시 둘, 밤 배정. 전부 데이터다.

## 무엇을 했나

스텝 8·10 이 받았던 것과 같은 처리를 했다. 계산을
`apex/analysis/light_curve/raw_lightcurve.py` 로 옮기고, **창이 그것을
상속한다.**

```
LightCurveBuilderWindow(StepWindowBase, RawLightCurveBuilder)
```

옮긴 것: 모듈 함수 12 개 + 메서드 17 개, 합쳐 **1,149 줄**.
창은 5,041 → 3,892 줄이 됐다.

동일성은 이제 검사 대상이 아니라 구조다:

```python
>>> LightCurveBuilderWindow._build_ensemble_series is RawLightCurveBuilder._build_ensemble_series
True
```

20 개 메서드 중 17 개가 이렇게 참이고, 나머지 셋은 **일부러** 다르다 —
`log`(창은 로그창에 쓴다)와 진행표시 훅 둘(`_preload_progress`,
`_preload_finished`). 미리읽기가 워커 패널을 갱신하는 게 유일하게 창에 남을
이유가 있던 부분이라, 기반 클래스에서는 아무것도 안 하는 훅으로 두고 창이
자기 Qt 신호로 덮는다. 배치 실행은 진행표시를 건너뛰고 같은 읽기를 한다.

## 대조 — 창이 저장한 커브와 같은 수를 내는가

**YZ Boo 2 밤 364 프레임** (`E:/APEX_validation/reprocess/YZBoo_2n/result`,
Moravian C3-61000, g·r·i, 2025-04-29~30). 창이 2026-08-01 에 저장한
`lightcurve_ID153_raw.csv` 가 기준이다. 대조를 오염시키지 않으려고
`lc_lightcurve/` 는 복사하지 않고 스텝이 직접 만들게 했다.

| 열 | 저장본 대비 최대 차이 |
|---|---|
| `mag` · `mag_err` | **0.000e+00** |
| `comp_avg` · `comp_err` | **0.000e+00** |
| `diff_mag_raw` · `diff_err` · `diff_mag` | **0.000e+00** |
| `airmass` · `JD` · `BJD_TDB` · `rel_time_hr` | **0.000e+00** |
| `color_index` · `color_index_ref` | **0.000e+00** |
| `night_id` | **2** ← 아래 |
| 문자열 열 6 개 (`filter`·`date`·`dataset`·`photometry_source`…) | 전부 일치 |

364 행이 `(파일, 필터)` 로 전부 대응하고 열 22 개가 순서까지 같다. 만들어진
파일도 10 개로 같다 (대상 커브 1 · 체크성 커브 7 · `comp_qc_summary.csv` ·
`comp_selection.json`).

### `night_id` 가 다른 이유 — 새 코드가 맞다

저장본은 두 밤을 **모두 `night_id 0`** 으로 뒀다. 새 실행은 4/29 를 1,
4/30 을 2 로 갈랐다.

이건 GUI 와 헤드리스가 갈린 게 아니라 **대조본이 수정 이전**이라 갈린 것이다.

- 저장본 작성: **2026-08-01 12:02**
- 밤 추론 폴백 커밋 `3d52f7b`: **2026-08-03**, 제목이 정확히
  *"fix(lc): headless 워크스페이스에서 night_id 가 전부 0 이던 것"*

밤 배정은 GUI 스텝 1 의 mixin 이 만든다. 그게 없는 워크스페이스에서는 전 프레임이
0 으로 떨어져 밤별 영점 보정이 두 밤을 한 밤으로 취급했다. 그 수정이 DATE-OBS 에서
밤을 추론하고, 지금 창과 헤드리스가 **둘 다** 그 코드를 쓴다.

## 옮기면서 드러난 결함 셋

이관 자체보다 **실제로 돌려본 것**이 잡았다.

1. **`LcTargetStep` 이 없는 파일을 찾고 있었다.** 마스터 카탈로그를
   `master_sources.csv` 로 적어 뒀는데 스텝 6 이 쓰는 건 `ref_catalog.tsv` 다.
   2026-08-19 에 그 스텝을 만들 때 이름을 확인 안 하고 적었고, 그래서
   **어떤 워크스페이스에서도 항상 `blocked` 였다.** 스텝을 파이프라인으로
   돌려보기 전까지는 안 보였다 — 관문이 막는 게 정상 동작처럼 보이기 때문이다.
   구분자도 쉼표에서 탭으로 고쳤다.

2. **문자열 치환이 메서드 이름을 갈았다.** 모듈을 `ensemble_series` →
   `raw_lightcurve` 로 옮기며 이름을 바꿀 때, `_build_ensemble_series` 안의
   `ensemble_series` 까지 바뀌어 `_build_raw_lightcurve` 가 됐다. 창은 여전히
   옛 이름을 부르니 `AttributeError`. 3 곳.

3. **딸려와야 할 이름 둘이 안 왔다.** `_build_source_to_id_map` 과
   `_load_check_star_meta_by_filter`. 의존성 추적을 `self.X(` 로만 해서
   모듈 수준 호출을 놓쳤다. 체크성 내보내기가 조용히 실패했다
   (`Check star export failed: name ... is not defined`) — 예외를 잡아 로그만
   남기는 자리라 라이트커브는 정상으로 보였다. `ast` 로 미정의 이름을
   훑어서 찾았다.

## 함정 — 레거시 워크스페이스는 빈 커브가 나온다

`_load_selection_ids_by_filter` 는 `step8_selection_dir()`(= `lc_selection/`)만
본다. 옛 이름(`step9_selection/`)을 쓰는 워크스페이스에서는 `{}` 를 돌려주고,
그러면 코드가 `df["ID"] == target_id` 로 떨어진다. **옛 측광 TSV 에는 `ID` 열이
없다** (`det_uid` 로만 식별) — 결과는 0/77 유효점, 조용한 빈 커브.

이건 이번 이관으로 생긴 게 아니다. HEAD 에서도 같은 코드였다. 현행
파이프라인은 스텝 7 이 `ID`·`source_id` 를 측광 TSV 에 직접 쓰므로
(`forced_photometry.py:774-783`) 영향이 없고, GUI 스텝 8 은 레거시를 **읽어서**
새 이름으로 **쓰므로** 창에서 한 번 열면 워크스페이스가 승격된다.

남은 위험은 **레거시 워크스페이스에 헤드리스만 돌리는 경우**다. 두 파일 옆에
레거시 인식 헬퍼(`selection_input_dir`, `step_paths_lc.py:72`)가 이미 있는데
이 읽기만 안 쓴다. 고칠지는 사용자 판단으로 남긴다 — 레거시 지원 범위를
넓히는 결정이라 이관과 분리한다.

## 지금 LC 파이프라인

```
1 scan  2 crop  3 sky  4 detect  5 wcs  6 refbuild  7 forcedphot
8 lctarget   ← 2026-08-19 (설정에서 대상을 읽고, 없으면 막는다)
9 lclightcurve ← 이번
```

남은 것은 10(추세제거)·11(주기분석)과 그림 5 장이다.

---

# LC 스텝 11 헤드리스 — 옮길 게 없었다

2026-08-19, 같은 날. 스텝 9 와 정반대의 작업이었다.

## 왜 쉬웠나

스텝 9 는 계산 1,149 줄이 창 안에 있었다. 스텝 11 은 **창에 계산이 없다.**
`PeriodAnalysisWorker` 를 열어 보면 순수 전달자다 — 배열을 받아
`run_period_analysis()` 를 부르고 돌아온 사전을 신호로 내보낸다. 나머지도 전부
이미 서비스였다:

| 하는 일 | 어디에 있었나 |
|---|---|
| 입력 파일 고르기 | `step_paths_lc.find_best_lightcurve_csv()` |
| 커브 읽기 | `period_io_service.load_period_lightcurve_csv()` |
| 주기 찾기 | `period_analysis_service.run_period_analysis()` |
| 결과 저장 | `period_io_service.save_period_analysis_outputs()` |

없던 건 **스핀박스에 값을 넣어 줄 사람 없이 이것들을 부르는 호출자**뿐이다.
그래서 스텝 파일 170 줄이 전부이고, 창에서는 아무것도 안 뺐다.

## 설정 5 개를 새로 열었다

대상 별과 달리 주기 탐색 값은 **전부 방어 가능한 기본값이 있다** — 창이 쓰던
그 값들이다. 그래서 막지 않고 기본값으로 돈다.

    lightcurve.period_min_days          0.01    탐색 하한 (일)
    lightcurve.period_max_days          10.0    탐색 상한 (일)
    lightcurve.period_samples_per_peak  10      주파수 격자 촘촘함
    lightcurve.period_methods           ls,pdm  쓸 방법
    lightcurve.period_pdm_bins          10      PDM 위상 구간 수

이 중 **실제로 대상마다 정해야 하는 건 탐색 창**이다. 창 밖에 있는 주기는
그냥 안 찾아진다. 0.01–10 일이면 δ Scuti 부터 대부분의 식쌍성까지 든다.

## 대조 — 창이 저장한 주기와 같은가

**YZ Boo 2 밤 364 점** (`reprocess/YZBoo_2n`). 창이 저장한
`period_analysis_all_ID153.json` 기준.

| 방법 | 저장본 (일) | 헤드리스 (일) | 차이 |
|---|---|---|---|
| `raw_ls` | 0.09446790 | 0.09446790 | **0.000e+00** |
| `corr_ls` | 0.09528867 | 0.09528867 | **0.000e+00** |
| `raw_pdm` | 0.10429465 | 0.10429465 | **0.000e+00** |
| `corr_pdm` | 0.10529656 | 0.10529656 | **0.000e+00** |

**첫 대조는 안 맞았고, 그게 코드 탓이 아니었다.** 저장본은 `filter="all"` 로
세 밴드를 합쳐 분석했는데 내 기본값은 밴드별이었다. 조건을 맞추자 전부 0 이
됐다. 입력 파일도 달랐는데(`_current.csv` vs `_offset.csv`) **두 파일은 MD5 가
같았다** — 창이 쓰는 「현재 선택본」이 offset 판의 사본이다.

**문헌 대조**: `raw_pdm` = 0.104295 일, YZ Boo 문헌값 0.104092 일 대비
**0.19 %** (17 초). LS 는 0.0945 일로 별칭을 집었다 — 알려진 거동이고, PDM 이
맞는 값을 준다는 점이 밴드별·방법별로 같이 봐야 하는 이유다.

## 밴드별이 기본인 이유

주기는 별의 성질이지만 **밴드마다 따로 잰다.** 세 밴드가 다른 주기를 주면 그건
결과가 아니라 진단이다 — 합쳐 버리면 그 신호가 사라진다. 그래서 기본은 밴드별
분석이고, 합치려면 `lightcurve.filter = all` 을 명시해야 한다.

## 스텝 10 은 미룬 게 아니라 다른 종류의 일이다

레지스트리에서 10 은 `DeferredStep` 이다. 스텝 9 는 창 의존이 전부 데이터라
통째로 떼어졌지만, 스텝 10 은 **계산이 위젯에서 입력을 읽고 위젯에 결과를 쓴다** —
`fit_and_apply`·`_run_sysrem`·`_run_global_ensemble` 이 모두 `self.target_edit`
를 읽고 `_update_results_table` 을 부른다. 폐포 69 개 메서드 3,314 줄 중 29 개가
GUI 를 만진다. 이건 이관이 아니라 입력읽기·계산·표시를 가르는 재설계다.

**스텝 11 은 10 을 기다리지 않는다.** `find_best_lightcurve_csv` 가 추세제거본이
있으면 그걸, 없으면 스텝 9 의 원시 커브를 쓴다. 추세제거는 답을 좋게 하지
답의 전제가 아니다.

## 지금 LC 파이프라인

```
1 scan  2 crop  3 sky  4 detect  5 wcs  6 refbuild  7 forcedphot
8 lctarget      2026-08-19  설정에서 대상을 읽고, 없으면 막는다
9 lclightcurve  2026-08-19  차등 라이트커브 (창 저장본과 0.0e+00)
10 lcdetrend    deferred — 창에서만 (재설계 필요)
11 lcperiod     2026-08-19  주기분석 (창 저장본과 0.0e+00)
```

## 그림도 같이 나왔다 — 그리고 거기서 결함 셋이 더 나왔다

주기 요약 그림(라이트커브·주기도표·위상접기 3단)은 창의 483 줄짜리 메서드였다.
`apex/analysis/light_curve/period_plot.py` 로 옮기고 창이 상속한다 — 17 개 중
**11 개가 동일 객체**, 나머지 6 개는 창이 위젯에서 읽어야 하는 것들이다
(도형·캔버스폭·재그리기·탐색창·별칭표시·체크성).

창은 2,413 → 1,915 줄.

**레이어 규칙을 지키느라 색을 지연조회로 바꿨다.** 그림이
`apex.gui.theme.Tokens` 에서 색 7 곳을 읽는데, `analysis` 는 `gui` 를 임포트하면
안 된다. 모듈 수준 임포트를 지우고 `_colors()` 함수 안으로 넣었다 — GUI 가 있으면
살아 있는(테마가 바뀐) 토큰을 그대로 쓰고, 없으면 같은 값의 기본 팔레트를 쓴다.
**창의 테마 변경이 그림에 반영되던 동작이 그대로 유지된다.**

### 결함 1 — 별칭 분석을 안 돌려서 엉뚱한 주기로 접었다

첫 그림은 **0.094468 일**로 위상을 접었다. 창 저장본은 **0.104209 일**.

채택 주기는 `alias_analysis["adopted_period"]` 에서 오는데 내 스텝은 별칭 분석을
아예 안 돌렸고, 그러면 코드가 Lomb-Scargle 최고점으로 떨어진다. **두 밤짜리
관측에서 LS 최고점은 흔히 관측 창의 별칭이다.**

    LS 단독        0.094468 일   문헌 대비  9.2 %
    별칭 해소 후   0.104209 일   문헌 대비  0.11 %   ← 창과 동일

설정 `lightcurve.period_resolve_aliases` 를 열고 **기본값을 켜짐**으로 했다.
창은 체크박스로 주지만 배치 실행에는 체크할 사람이 없고, 꺼져 있으면 조용히
별칭을 답으로 낸다.

### 결함 2 — 밤 기준선을 지우는 보정 뒤에 보정본을 썼다

별칭을 붙였더니 이번엔 **0.105359 일** — 창의 0.104209 와 1.15e-3 일 차이.
`input_series` 를 보니 창은 `raw`, 나는 `corrected` 였다.

**야간 영점 보정(nightly offset)은 밤 사이 기준선을 없앤다.** 그건 여러 밤에
걸친 주기 탐색이 쓰는 바로 그 정보다. 창의 워커는
`correction_preserves_nightly_baseline` 을 함께 보는데 나는 「보정본이 있으면
쓴다」로만 짰다. 그 플래그는 이미 `load_period_lightcurve_csv` 가 실어 준다.

    조건 무시   0.105359 일   문헌 대비  1.22 %
    조건 반영   0.104209 일   문헌 대비  0.11 %   ← 창과 0.000e+00

### 결함 3 — 체크성을 내가 넘겨서 오히려 잃었다

`check_star_data` 를 스텝에서 만들어 넘겼는데, 필터가 `"all"` 일 때
`load_check_star_csv(filt="all")` 은 아무것도 못 찾는다. 넘기지 않으면 플로터가
**실제로 그린 필터**로 스스로 찾는다 — 창이 하는 그대로다. 넘기던 걸 뺐더니
체크성 364 점이 창 그림과 똑같이 들어왔다.

### 결함 4 — 이관이 창의 스핀박스를 끊었다

그림을 빼내면서 탐색창·별칭표시·체크성을 접근자로 만들고 기본값을 줬는데,
**창에 그 접근자를 안 덮어 뒀다.** 그 상태로는 창의 스핀박스를 돌려도 그림이
기본값(0.01–10 일)으로 그려진다. 이관 직후 `is` 대조로 잡았고
(`test_the_window_still_reads_its_own_controls`), 그 테스트가 이 결함 때문에
존재한다.

## 최종 대조 — YZ Boo, 창 저장본 기준

| | 창 | 헤드리스 | 차이 |
|---|---|---|---|
| `raw_ls` | 0.09446790 | 0.09446790 | **0.000e+00** |
| `corr_ls` | 0.09528867 | 0.09528867 | **0.000e+00** |
| `raw_pdm` | 0.10429465 | 0.10429465 | **0.000e+00** |
| `corr_pdm` | 0.10529656 | 0.10529656 | **0.000e+00** |
| **채택 주기** | 0.10420857 | 0.10420857 | **0.000e+00** |
| 별칭 판정 | RESOLVED | RESOLVED | 같음 |
| 입력 계열 | raw | raw | 같음 |
| 그림 | 3 단 + 체크성 364 | 3 단 + 체크성 364 | 같음 |

문헌 YZ Boo P = 0.104092 일 대비 **0.11 %**.

라벨만 다르다 — 창은 저장 당시 캔버스가 900 px 미만이라 짧은 라벨(`Diff. mag`)을
썼고 배치는 넓은 캔버스라 전체 라벨을 쓴다. 그게 반응형 레이아웃의 의도다.

## 두 번째 대상이 정반대 결함을 드러냈다 — AE UMa

YZ Boo 하나로 끝냈으면 못 봤을 것이다. **AE UMa**(SX Phe 맥동성, 문헌
P = 0.086017 일, 2 밤 · 514 점)로 같은 스텝을 돌렸다.

    corr_ls    0.08601130 일    문헌 대비  0.01 %
    corr_pdm   0.08604727 일    문헌 대비  0.04 %
    채택       0.08251356 일    문헌 대비  4.07 %   ← 나쁘다

**별칭 해소가 답을 망쳤다.** 그런데 서비스는 잘못이 없었다 —
`status: AMBIGUOUS`, 이유는 *"Leave-one-night-out agreement is 50%"*. 즉
**「나는 이걸 못 가렸다」고 명시했는데 호출하는 쪽이 그 말을 무시하고 1 순위
후보를 채택**했다. 후보 2 순위는 0.086079 일로 정답에 가까웠다.

`analyze_period_aliases` 는 세 상태를 낸다 — `RESOLVED` · `AMBIGUOUS` ·
`INSUFFICIENT`. 처방은 **`RESOLVED` 일 때만 주기도표 최고점을 덮는 것**이다.
`_summary_period` 를 그렇게 고쳤고, **창도 같은 코드를 상속하므로 같이 고쳐졌다.**

    AE UMa   AMBIGUOUS → 최고점 유지    0.086011 일   0.01 %  ✓
    YZ Boo   RESOLVED  → 후보 채택      0.104209 일   0.11 %  ✓

두 대상이 정반대 경우이고, 이제 각각 맞게 처리한다. 그림에도
`alias: AMBIGUOUS` 가 찍혀서 읽는 사람이 그 사정을 안다.

**교훈**: 「기본값이 방어 가능하다」와 「기본 동작이 옳다」는 다른 문제다.
별칭 해소를 기본 켜짐으로 둔 건 맞았지만(YZ Boo 9.2 % → 0.11 %), 그 결과를
**무조건** 믿은 건 틀렸다(AE UMa 0.01 % → 4.07 %). **하나의 대상으로는 이
두 가지가 구분되지 않는다.**

곁가지로 보고 메시지 버그도 나왔다 — 필터 개수를 `written` 길이로 세는데
거기엔 그림과 주기도표 CSV 4 개가 같이 들어 있어서, 한 필터 실행이
「2 filter(s)」로 나왔다.

---

# LC 스텝 10 (추세제거) — 「재설계」가 아니었다

2026-08-19. 마지막 남은 deferred 스텝이고, **미룬 이유가 틀렸다.**

## 내가 잘못 읽었던 것

「계산이 위젯에서 입력을 읽고 위젯에 결과를 쓴다. 폐포 69 개 메서드 3,314 줄
중 29 개가 GUI 를 만진다. 이건 이관이 아니라 재설계다 — 반나절 이상.」

**절반만 맞았다.** 쓰기는 전부 표시가 맞다. 읽기는 아니었다 — 창의
`_sync_state_from_controls` 가 **스핀박스 전부를 한 곳에서 평범한 속성으로
복사**하고 있었고, 계산은 처음부터 위젯이 아니라 속성을 읽고 있었다. 그 통로를
우회하는 직접 읽기는 **네 개**뿐이었다:

    self.target_edit.text()        →  self._target_id_text()
    self.filter_combo.currentText()→  self._filter_selection()
    self.chk_global_k2.isChecked() →  self._use_global_k2()
    self.date_list (in _selected_dates) → 속성 또는 데이터에서 유도

**경계가 이미 있었다.** 내가 위젯 접촉 횟수만 세고 그게 읽기인지 쓰기인지
안 갈랐다.

## 무엇을 옮겼나

`apex/analysis/light_curve/detrend_runner.py` — 메서드 58 개 + 모듈 함수 6 개.
창은 **4,874 → 1,580 줄**. 70 개 중 **59 개가 동일 객체**이고, 나머지 11 개는
접근자 넷과 Qt 배관 넷, 그림 훅 셋이다.

**Qt 를 코드에서 완전히 걷어냈다.** `QMessageBox` 17 곳이 `self._tell_user(...)`
가 됐다 — `PsfPhotometryRunner` 가 `.emit()` 을 버릴 때 한 것과 같은 교환이다.
호출부는 이미 전부 `if not silent:` 나 `if update_ui:` 뒤에 있어서 배치 경로는
닿지도 않았지만, **옮긴 메서드는 정의된 모듈의 전역을 본다** — 창이 상속해도
`detrend_runner` 의 이름을 찾으므로, 모듈에 스텁을 두면 창의 대화상자까지
죽는다. 그래서 채널이 답이다. `refresh()`(테마 스타일 재적용)도 같은 이유로
`self._refresh_style()` 이 됐다.

## 그림 4 장도 같이 나온다

`_update_plots`(210 줄)와 `_plot_global_diagnostics`(91 줄)는 **Qt 를 한 줄도
안 쓴다** — matplotlib 뿐이다. 헤드리스에서 안 나오던 이유는 그리기가 Qt 를
필요로 해서가 아니라 **손에 잡히는 figure 가 `FigureCanvas` 것뿐**이었기
때문이다. `_plot_figure()` 훅 하나로 갈렸다.

### 백지 그림 — 세션에서 두 번째

첫 실행이 PNG 를 정상적으로 저장했고 **내용이 완전히 백지**였다. 모든
`_update_plots()` 호출이 `if update_ui:` 안에 있어서 배치는 그리지 않고,
`savefig` 는 빈 figure 도 불평 없이 쓴다. 저장 직전에 `_ensure_plot_drawn()` 을
부르게 했고, **창은 그걸 no-op 으로 덮는다** — 창의 캔버스는 이미 사용자가 보고
있는 것이라 저장할 때 다시 그리면 보이는 것과 달라진다.

### 패널 수가 다른 건 창의 보기 모드였다

창 저장본은 1 단, 내 것은 3 단이었다. 원인은 `_plot_view_mode` 로,
기본값이 `"corr"`(보정곡선만)다. 캔버스가 조작부와 화면을 나눠 쓰기 때문이다.
**파일은 그런 제약이 없고**, 단일 보기가 감추는 두 패널(원본 곡선, 에어매스 대
Δmag)이 바로 그림을 기록으로 만드는 것들이라 배치 기본값은 `all` 로 뒀다.
설정 `lightcurve.detrend_plot_view` 로 바꿀 수 있다.

## 대조 — 창이 저장한 보정곡선과 같은가

**YZ Boo 2 밤 364 점.** 창이 2026-08-01 에 저장한 것 기준.

| | 결과 |
|---|---|
| 보정곡선 26 열 전부 | **0.000e+00** |
| 적합 파라미터 6 행 | **0.000e+00** |
| 파일 9 개 (그림 2 장 포함) | 이름까지 동일 |

## 세 가지가 달라 보였고, 통제 실험이 전부 입력 탓으로 확정했다

파이프라인으로 8→9→10→11 을 통째로 돌리자 세 가지가 어긋났다.

    night_id            0 vs 1·2          차이 2.0
    적합표의 date        2025-04-29 vs "Night 1"
    별칭 판정            RESOLVED vs AMBIGUOUS

셋 다 **한 뿌리**다. 저장본은 밤 추론 수정(`3d52f7b`, 08-03) 이전이라 night_id 가
전부 0 이고, 지금 스텝 9 는 1·2 로 가른다. 그러면
`_fill_night_id` 가 `date` 를 "Night N" 으로 덮고(원래 설계된 동작),
`n_nights` 가 1→2 가 되어 leave-one-night-out 이 비로소 의미를 가지며 50 % 로
갈려 `AMBIGUOUS` 가 된다.

**통제 실험**: 저장본의 입력(night_id 전부 0)을 그대로 주고 현재 코드를 돌렸다.

    보정곡선  364/364 행   0.000e+00
    적합표      6/6 행     0.000e+00
    별칭       RESOLVED 0.10420857  (n_nights=1) — 저장본과 동일

**차이는 전부 입력이지 이관이 아니다.** 그리고 저장본의 `RESOLVED` 는
night_id 버그의 부산물이었다 — 밤을 제대로 가르면 두 밤으로는 별칭이 안 갈린다는
것이 정직한 답이다.

## 그래서 남은 문제 — 갈리지 않을 때 무엇을 답으로 낼 것인가

밤을 제대로 세면 YZ Boo 는 `AMBIGUOUS` 이고, 그러면 코드가 LS 로 떨어진다.

    LS   0.094468 일   문헌 대비  9.2 %
    PDM  0.104295 일   문헌 대비  0.19 %

**대상 둘로는 규칙을 만들 수 없다.** AE UMa 에서는 LS 가 0.01 % 로 더 좋았다.
그래서 순서를 바꾸지 않고 **불일치를 드러내게** 했다:

    ok: target ID 153, 1 filter(s) — all LS=0.095289 d [alias ambiguous]
        — methods disagree by 11 %

조용히 하나를 고르면 9 % 오차가 숫자로 발표된다. 갈리지 않았다고 말하는 게
이 데이터가 지지하는 전부다. **어느 방법을 우선할지는 사용자 판단으로 남긴다.**

## 지금 LC 파이프라인 — deferred 없음

```
1 scan  2 crop  3 sky  4 detect  5 wcs  6 refbuild  7 forcedphot
8 lctarget      대상 선택 (설정에서 읽고, 없으면 막는다)
9 lclightcurve  차등 라이트커브        창 저장본과 0.000e+00
10 lcdetrend    추세제거 + 그림 2 장    창 저장본과 0.000e+00
11 lcperiod     주기분석 + 그림 1 장    창 저장본과 0.000e+00
```

---

# 사용자 결정 셋을 닫았다 (2026-08-20)

## D-014 · 별칭이 안 갈릴 때 — 「고르지 않는다」로 닫혔다

**내가 질문을 잘못 냈다.** 「이건 관측 판단이니 선생님이 정하세요」로 올렸는데,
사용자는 관측천문·측광 전문가이지 시계열 주기추정 방법론 전공이 아니다.
분야 경계를 잘못 그었고, 하루 전 C-034 와 같은 실수다 (OPERATOR C-036).

사용자 답이 더 낫다:

> 뭐 걍 문헌이랑 제일 비슷한걸 사용자가 고르던가 해야지 뭐, **애초에 정보가
> 부족한건데**, pdm이 맞는지 LS가 맞는지 전공자가 아니라서 몰라

**두 밤으로는 원리상 안 갈린다.** 그러면 코드가 하나를 고르는 것 자체가 틀린
설계다. 규칙을 정하는 게 아니라 **후보를 문헌과 대조할 수 있게 내놓는 것**이 답이다.

### 무엇을 만들었나

`period_candidates_<filter>_ID<n>.csv` — 이 실행이 찾은 주기 전부를, 주기 순으로.

| source | method | series | period_days | period_hours | rank | note |
|---|---|---|---|---|---|---|
| alias candidate | alias resolver | raw | 0.086514 | 2.076 | 8 | |
| periodogram | Lomb-Scargle | raw | 0.094468 | 2.267 | | |
| alias candidate | alias resolver | raw | 0.094605 | 2.271 | 3 | |
| periodogram | Lomb-Scargle | corrected | 0.095289 | 2.287 | | |
| **alias candidate** | alias resolver | raw | **0.104209** | **2.501** | 1 | adopted |
| **periodogram** | **PDM** | raw | **0.104295** | **2.503** | | |
| periodogram | PDM | corrected | 0.105297 | 2.527 | | |
| alias candidate | alias resolver | raw | 0.115791 | 2.779 | 2 | |
| … | | | | | | |

문헌 YZ Boo **0.104092 일 = 2.498 시간**. 표를 보면 **가까운 셋이 2.50~2.53 시간에
모여 있고** 나머지는 2.08·2.27·2.78 시간으로 흩어져 있다 — 사람이 문헌 한 줄만
알면 즉시 고른다. 시간 단위를 같이 넣은 것도 그래서다(문헌은 두 단위로 인용된다).

메시지도 표를 가리킨다:

    all LS=0.095289 d [alias ambiguous] — methods disagree by 11 %;
      candidates in period_candidates_all_ID153.csv

**여전히 자동으로 고르지는 않는다.** 갈렸을 때(`RESOLVED`)만 별칭 후보가 주기도표
최고점을 덮고, 안 갈렸으면 최고점을 두되 표를 남긴다.

## D-015 · 중복 삭제 — 30.4 GB 회수, 이름으로 불리는 것은 하나도 안 건드렸다

사용자: **「중복 지워」**

지운 것: **545 개 / 30.4 GB** (E 여유 105 → **136 GB**).

| 지운 곳 | |
|---|---|
| `reprocess_cr/M3/calibrated/` | 3.42 GB |
| `reprocess/M5/sci` · `M3/sci` · `NGC6811/sci` · `M13/sci` | 5.97 GB |
| `reprocess/*/calibrated/`, `YZBoo_20250430/calibrated/` | 3.6 GB |
| `psf_engines/qfit_ab/qfit_*/sci` | 3.4 GB |
| `masters/` 마스터 프레임 | 5.9 GB |
| IRAF 작업본 `frame.fits` | 2.8 GB |

**남긴 것**: 표(csv·tsv·json) 163 개 — 크기가 0.21 GB 뿐이고 측정값이라 잃으면
비싸다. AstrylStudio 의 NGC6888 중복 176 개 / 12.8 GB — 다른 프로젝트다.

### 안전을 어떻게 확인했나

**정적 도달성 분석이 두 번 틀렸다.** 첫 판은 정규식 `[^"']*` 가 개행을 넘어
docstring 을 통째로 삼켜 「경로 조각 7,406 개」가 대부분 쓰레기였다. 고친 뒤에도
`photometry_*.tsv` 글롭이 안 걸렸다.

그래서 **좁고 결정적인 질문**으로 바꿨다 — 논문 스크립트가 **손으로 이름을 적은
파일 57 개** 중 삭제 목록에 있는 게 몇 개인가. **0 개**였다.

삭제 뒤 확인:

- 손으로 적은 7 개 프레임 전부 존재 (`M67/sci/` 3 개, `sci_nocr/` 4 개)
- `fig_completeness_realvssynth.py` **실제 재생성 성공** — 7 개 프레임 다 읽고
  σ_e·FWHM·m50·S/N50 을 예전과 같은 형태로 출력

`sci_nocr` 는 애초에 삭제 후보에 **0 개**였다. 그 폴더가 중요한 이유가 스크립트
주석에 있다 — M13·NGC 6811 은 우주선 제거로 재환원되어 `sci/` 가 논문이 쓴
프레임이 아니고, **주입 당시 프레임은 `sci_nocr/` 에만 남아 있다.**

**되돌릴 수 있다**: 지운 545 개는 전부 같은 바이트가 다른 경로에 있고,
`restore_map.csv` 에 지운 경로 ↔ 남은 경로가 적혀 있다.

## D-013 · 옛 워크스페이스 — 고치려다 「고칠 수 없다」가 나왔고, 그래서 말하게 했다

사용자: **「13은 알아서 해」**

처음엔 두 줄로 될 줄 알았다. `_load_selection_ids_by_filter` 가 `lc_selection/` 만
보니 레거시 인식 헬퍼(`selection_input_dir`)로 바꾸면 끝이라고 봤다. 바꿨고,
**선택 지도는 실제로 로드됐다** (0 → 3 필터).

그런데 **여전히 0/77 이었다.** 두 번째 원인이 있다:

    옛 측광표 열: det_uid, x_det, y_det, xcenter, ycenter, FILTER, ...
    → source_id 없음, ID 없음

선택 지도는 Gaia `source_id` 로 별을 찾는데 **옛 측광표에는 그 열이 없다.**
`step8_idmatch/frame_sourceid_to_ID.tsv` 에 `source_id ↔ ID` 는 있지만 `det_uid` 와는
**위치로 다시 맞춰야** 이어진다 — 옛 스텝 8 을 재구현하는 일이다.

**비례에 맞지 않는다.** 해당 워크스페이스는 9 개이고 전부 YZ Boo·AE UMa 인데,
둘 다 현행 레이아웃 워크스페이스가 따로 있다. 그래서 **진짜 결함인 「조용함」만
없앴다** — 이제 이렇게 말한다:

    blocked: built 77 rows and none carry a measurement, because the photometry
    tables identify sources by det_uid and have neither an ID nor a source_id
    column. That is a workspace from before the current Step 7 — open it once in
    the GUI, which upgrades it, or re-run Steps 1-7 headless.

`selection_input_dir` 변경은 그대로 뒀다 — 대상 좌표 조회가 레거시 워크스페이스
에서도 되므로 BJD 계산에 도움이 된다.

