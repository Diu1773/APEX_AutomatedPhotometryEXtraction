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
