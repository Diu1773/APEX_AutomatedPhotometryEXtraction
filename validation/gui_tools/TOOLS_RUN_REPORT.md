# Step 12 + 도구 10종 실 GUI 구동 — 캡처와 문제 기록

**2026-08-01 야간** · 요청: 「Step 12 랑 tools 들까지 다 실 GUI 에서 써보고 png 로도
남기고, 문제 생기는 것도 다 기록」

도구: `scripts` 가 아니라 **사람이 메뉴에서 여는 그 진입점**(`main_window` 의
launcher 메서드)을 그대로 호출했다. 창을 띄운 뒤 **주 실행 버튼까지 눌러**
결과가 나온 상태를 `widget.grab()` 으로 저장했다 — `validation/gui_tools/`.

## 1. 결과 한눈에

| 모드 | 대상 | 열림 | 비고 |
|---|---|---|---|
| CMD (M13) | Step 12 + 도구 5 | **6/6** | Step 12 CMD 플롯·이소크론 정상 |
| LC (YZBoo_2n) | 도구 8 | **8/8** | QA Report 364×10 생성 |

## 2. 찾아서 고친 버그 3건

### 2.1 `numpy.bool_` 이 PyQt5 setter 로 샌다 — 창이 통째로 죽는다

`cluster_structure/window.py`

```python
has_mag = bool(...) and np.isfinite(...).sum() > 5   # ← `and` 는 뒤를 그대로 돌려준다
self.chk_bright.setEnabled(has_mag)                  # TypeError: numpy.bool
```

파이썬 `A and B` 는 A 가 참이면 **B 를 그대로** 돌려준다. 뒤쪽 비교가
`numpy.bool_` 이라 PyQt5 가 거부하고, 그 예외로 **`load_inputs` 전체가 실패**했다.
증상은 「입력이 안 잡힌 빈 창」이라 원인이 안 보인다 — 창 로그의

```
ERROR load_inputs: setEnabled(self, a0: bool): argument 1 has unexpected type 'numpy.bool'
```

한 줄로 잡았다. 고친 뒤 같은 창이 1280×1432 · 캔버스 3 으로 정상 표시된다.

**같은 패턴이 하나 더 있었다** — `variable_star.py:3942`
(`np.isfinite(adopted_period) and adopted_period > 0`). 아직 안 터졌을 뿐이라
함께 고쳤다.

### 2.2 airmass 도구를 열면 **뒤이어 여는 도구들이 엉뚱한 폴더를 본다**

`airmass_debug.py` 의 `_restore_file_selection_state()` 가

```python
self.params.P.result_dir = self.params.P.data_dir / "result"   # 항상 파생
```

로 **저장된 `result_dir` 을 무시**한다. `main_window._bootstrap_file_selection_state`
는 저장값을 먼저 쓰는데 이 도구만 규칙이 다르다.

`params` 는 창들이 공유하므로 한 번 열면 그 세션 전체가 오염된다. 실제로

```
Variable Star Analysis
  Status: No lightcurve_*.csv found in E:
  Workspace: …\YZBoo_2n\sci\result      ← 실제는 …\YZBoo_2n\result
```

가 찍혔다. data_dir 과 result_dir 이 다른 워크스페이스에서
**variable_star · transit · eclipsing_binary 가 연쇄로 광곡선을 못 찾는다.**
저장값을 먼저 쓰도록 맞췄다.

**수정 확인** — 같은 창을 다시 캡처하면(`lc_variable_star.png`)

| | 고치기 전 | 고친 뒤 |
|---|---|---|
| Status | `No lightcurve_*.csv found` (빨강) | **`lightcurve_ID153_ra… 124 pts [Raw]`** (초록) |
| Workspace | `…\YZBoo_2n\sci\result` | **`…\YZBoo_2n\result`** |
| Filter | — | g |

하네스에도 감시를 넣었다 — 도구를 하나 열 때마다 `params.P.result_dir` 이
바뀌면 `result_dir_changed` 로 기록한다. 수정 후 재실행에서는 경고가 없다.

> Run 버튼은 여전히 비활성인데 이건 설계된 게이트다 —
> 「Load a validated Step 12 release to run」. LC 도구들은 Step 11 완주 +
> 릴리스가 선행 조건이고, 이번 워크스페이스는 Step 9 까지만 돌렸다.

## 3. 데이터 한계 — 도구 잘못이 아니다

### `extinction_fit` 이 CMD·LC 양쪽에서 실패한다

`Extinction fit failed: No valid per-star extinction fits produced`

airmass 범위를 재보면 이유가 분명하다.

| 자료 | airmass | 폭 |
|---|---|---:|
| M13 (15프레임) | 1.0107 ~ 1.0130 | **0.0023** |
| YZ Boo 2밤 (364프레임) | 0.9997 ~ 1.2165 | 0.2168 |

소광계수는 `m = m0 + k·X` 의 **기울기**다. M13 은 천정 부근에서만 찍어 X 가
사실상 고정이라 기울기를 잴 수가 없다. YZ Boo 의 0.217 도 좁다(보통 X = 1~2 필요).

> **개선 여지**: 메시지가 「per-star fit 이 하나도 없다」라 사용자가 원인을
> 알 수 없다. 「airmass 폭이 0.002 라 소광계수를 결정할 수 없다」로 바꾸면
> 자료를 더 모아야 한다는 판단을 바로 할 수 있다.

### `iraf_photometry` — `IRAF dir not found: …/result/iraf_phot`

IRAF 비교 산출물이 없어서다(그 파이프라인을 안 돌렸다). 도구는 정상이다.

## 4. 기록해 둘 것 — Step 12 그림의 제목이 틀린다

M13 워크스페이스로 Step 12 를 열면 플롯 제목이 **`NGC6811`** 로 나온다
(`cmd_step12_isochrone.png`). `reprocess/M13/parameters.toml` 의
`[target] name` 이 NGC6811 로 남아 있어서다 — 설정 복붙 흔적이고 코드 문제는
아니지만, **그림에 그대로 박히므로** 논문·보고서에 쓰기 전에 고쳐야 한다.

같은 이유로 `cluster_structure` 창 제목도 `NGC6811 - Analyze Cluster Structure` 다.

## 5. 내가 틀렸던 것

처음 캡처에서 네 도구(gaia_3d_viewer · airmass_debug · variable_star ·
eclipsing_binary)가 **100×30 또는 0×0** 으로 나와 「창을 못 연다」고 적었다.
**하네스 문제였다** — launcher 는 창을 반환하지 않고 `self.qa_window` 처럼
도구마다 다른 속성에 넣는데, 내가 최상위 창 목록에서 마지막 것을 집어
숨은 헬퍼 위젯을 캡처했다. `main_window` 의 새 속성을 먼저 보도록 고치니
1906×1111 · 1303×816 · 1304×800 · 1170×1149 로 전부 정상이었다.

또 첫 캡처는 **글자가 하나도 없이 색 블록만** 나왔다. `QT_QPA_PLATFORM=offscreen`
이 시스템 폰트를 못 읽고, `apply_theme(app)` 도 빠뜨렸기 때문이다. 실제 앱은
세 진입점 모두 `apply_theme` 를 부른다 — 그걸 안 부르면 캡처가 실사용과 다르다.

## 6. 후속 (2026-08-01) — 남은 세 가지를 다 처리했다

### 6.1 타깃 이름이 틀린 워크스페이스 3개

좌표로 확인해 보니 **M13 · M3 · M67 셋 다 `[target] name` 이 `NGC6811`** 이었다
(좌표는 각각 250.421 / 205.548 / 132.825 로 정확하다 — NGC6811 설정을 복사해
좌표만 고친 흔적). 그림·창 제목에 그대로 박히므로 각각 제 이름으로 고쳤다.
확인: `cluster_structure` 제목이 `M13 - Analyze Cluster Structure` 로 바뀐다.

### 6.2 `extinction_fit` 에러가 원인을 말해 준다

`No valid per-star extinction fits produced` 뒤에 실제 airmass 범위와 폭을
붙이고, 폭이 0.3 미만이면 「소광계수는 airmass 기울기라 이 폭으로는 결정할 수
없다」고 알려 준다. 자료를 더 모아야 하는지 설정 문제인지 바로 갈린다.

### 6.3 Step 10 → 11 을 완주시키다 **치명적 회귀**를 찾았다

LC 도구들이 요구하는 「validated Step 12 release」를 만들려고 Step 10(디트렌드)
을 돌리고 Step 11 로 넘어갔더니 **창이 아예 안 열렸다.**

```
step11_period_analysis.py:682
    mode_label = _CORR_MODE_LABELS.get(pref, mode_label)
NameError: name '_CORR_MODE_LABELS' is not defined
```

`_CORR_MODE_LABELS` 는 `period_io_service.py` 에 있는데 step11 이 import 하지
않는다(언더스코어라 자동으로 안 온다).

> **발현 조건이 핵심이다.** 그 줄은 `load_detrend_preference()` 가 값을 돌려줄
> 때만 탄다 — 즉 **Step 10 을 거친 뒤에만 죽는다.** 어제 YZBoo_g 는 Step 10 을
> 안 돌려 preference 가 없어 통과했고, 오늘 Step 10 을 돌리자 막혔다.
> **LC 를 정석대로(10 → 11) 진행하면 반드시 걸리는 버그**다.

언더스코어 이름을 남의 모듈에서 끌어쓰는 대신 `corr_mode_label(key, default)`
공개 함수를 만들어 연결했다. 고친 뒤 Step 11 이 열리고
`lc_detrend/lightcurve_ID153_offset.csv` 를 읽어 `lc_period/` 를 낸다.

### 6.4 그런데 파이프라인 주기는 별칭을 집는다 (자료 한계)

Step 10 `offset` 디트렌드 → Step 11 (필터 `all`, 364점) 결과다.

| 방법 | 주기 | 문헌 대비 |
|---|---:|---:|
| raw_ls | 0.094468 | **−9.25%** (1일 별칭) |
| corr_ls | 0.095289 | −8.46% (별칭) |
| raw_pdm | 0.130339 | +25.21% |
| **corr_pdm** | **0.105297** | **+1.16%** |
| (손으로 g+r+i 정규화 후 합침) | 0.104550 | +0.44% |

`GUI_RUN_REPORT.md` §4.8 에 적은 「진짜 피크와 1일 별칭의 power 차이가 0.72%
뿐」이 그대로 발현된 것이다 — **방법이 조금만 달라도 다른 피크를 집는다.**
버그가 아니라 기저선 1.1일이 만드는 한계이고, 밤을 더 넣기 전에는 어느 값도
확정이 아니다. 처리 경로도 다르다(파이프라인은 밤별 offset 디트렌드, 손 계산은
필터별 중앙값 정규화).

## 7. 산출물

```
validation/gui_tools/
  run_tools_capture.py      하네스 (--run 으로 주 실행 버튼까지 누른다)
  cmd_*.png / lc_*.png      창 캡처
  cmd_tools.json            창별 상태·시간·표·다이얼로그
  lc_tools.json
```
