# APEX 산출물 디스크 조사 — 2026-08-28

**아무것도 지우지 않았다.** 어디에 무엇이 얼마나 있는지만 쟀다.

## 한 줄

**86.6 GB 중 83.2 GB(96 %)가 `.fit` 파일 2,901 개다.** 표·그림·기록을 전부 합쳐도
2 GB 가 안 된다. 그런데 **그 83 GB 를 지워도 되는지는 지금 판정할 수 없다** —
아래 「왜 아직 못 지우나」.

---

## 어디에 있나

| 뿌리 | 크기 |
|---|---|
| `E:\APEX_validation\reprocess` | 68.0 GB |
| `E:\APEX_validation_output` (레포의 `validation/` 링크 대상) | 18.6 GB |
| **합계** | **86.6 GB · 24,990 파일** |

`E:\APEX_validation` 과 `E:\APEX_validation_output` 은 **다른 경로**다. 뒤엣것이
레포 안에서 `validation/` 으로 보인다 — 2026-08-20 에 이 링크 때문에 레포 추적
파일이 「삭제 가능」 구역에 들어갔다(F-061).

### 성단·타겟별 (reprocess)

| 폴더 | GB |
|---|---|
| YZBoo_2n | 21.42 |
| NGC6811 | 10.88 |
| M5 | 9.97 |
| M67 | 8.03 |
| M13 | 6.82 |
| YZBoo_g | 5.67 |
| M3 | 3.54 |
| YZBoo_20250430 | 1.65 |

### 산출물 쪽 (validation_output)

| 폴더 | GB |
|---|---|
| **real_gui_run** | **16.29** |
| paper | 0.97 |
| _ndet_m3B | 0.58 |
| _ndet_m67g | 0.57 |
| reports | 0.11 |
| 나머지 열여덟 | 0.05 |

`real_gui_run` 은 파라미터 쓸기 실험이 쌓인 곳이다 — `psf_auto_v3` 4.44 GB,
`m3_fwhm_fixed_ast_fit12/20/24/30_v2` 각 1.15~1.49 GB, `m13_variants` 1.37,
`m13_repeatability` 0.98, `m5_seeing_sweep` 0.59 …

---

## 무엇이 자리를 쓰나

| 확장자 | GB | 파일 수 |
|---|---|---|
| **`.fit`** | **83.19** | 2,901 |
| `.fits` | 1.21 | 289 |
| `.csv` | 1.05 | 6,828 |
| `.tsv` | 0.70 | 1,212 |
| `.npy` | 0.15 | 5,272 |
| `.png` | 0.10 | 394 |
| `.pdf` | 0.08 | 67 |
| `.json` | 0.01 | 5,068 |
| `.ecsv` | 0.01 | 18 |

**논문이 실제로 인용하는 것 — 표·그림·기록 — 을 전부 합치면 약 2.0 GB 다.**
파일 수로는 18,000 개가 넘는데 자리는 2 %도 안 쓴다.

### 캐시는 문제가 아니었다

`cache/` 폴더 전부 합쳐 **0.36 GB**. 「캐시를 지우면 된다」는 처음 가설이었고
**틀렸다.**

### `.fit` 이 어느 폴더에 있나

| 폴더 이름 | GB | 개수 | 무엇 |
|---|---|---|---|
| `cmd_psf` | 38.32 | 2,106 | PSF 단계 산출물 |
| `sci` | 22.85 | 401 | 보정된 science 프레임 |
| `sci_pre20260807` | 7.24 | 127 | **이름이 옛 판본이라고 말한다** |
| `20250430` | 5.53 | 97 | 날짜 이름 — 원본일 가능성 |
| `sci_nocr` | 2.05 | 36 | 우주선 제거 안 한 변형 |
| `20260509` · `20260611` | 2.91 | 51 | 날짜 이름 |
| `_seq_control` | 1.65 | 29 | 통제 실험 |
| `input` · `data` | 2.17 | 38 | **원본일 가능성** |

---

## 왜 아직 못 지우나

`cmd_psf`(38 GB)와 `sci`(23 GB)는 **처리 단계의 산출물이므로 원리상 다시 만들 수
있다.** 합치면 61 GB, 전체의 70 % 다.

**그런데 다시 만들려면 무엇으로 만들었는지를 알아야 하고, 지금 그걸 아는 폴더가
하나도 없다.** 결과 대장이 세어 보면 저널 기록이 있는 폴더는 **0 개**다
(`scripts/result_ledger.py`). 소급도 안 된다.

즉 **디스크 문제와 출처 문제는 같은 문제다.**

- 83 GB 를 안전하게 지우려면 다시 만들 수 있어야 한다
- 다시 만들 수 있으려면 폴더가 무엇으로 만들어졌는지 말해야 한다
- 지금은 말하지 못한다

**그래서 지금 지우면 「용량은 벌었는데 논문 수치를 다시 못 만드는」 상태가 된다.**

---

## 지금 당장 할 수 있는 것

### 1. 2 GB 를 따로 보관한다 — 이건 조건이 없다

표·그림·기록(`.csv`·`.tsv`·`.ecsv`·`.png`·`.pdf`·`.json`·`apex_journal.jsonl`)은
**합쳐서 2 GB 가 안 되는데, 잃으면 다시 만들 수 없는 유일한 부분이다.** 83 GB 의
프레임은 원본이 있으면 다시 만들 수 있지만 표는 그 실행이 사라지면 끝이다.

용량이 작아서 어디에나 들어간다. **이건 판단이 필요 없는 작업이다.**

### 2. 이름이 스스로 옛 판본이라 말하는 것부터 본다

`sci_pre20260807` **7.24 GB**. 이름이 2026-08-07 이전 판본이라고 적고 있고 같은
자리에 현행 `sci` 가 있다. **다만 「이름이 그렇다」와 「실제로 대체됐다」는 다른
사실이므로, 지우기 전에 현행 `sci` 가 같은 프레임을 덮는지 확인해야 한다.**

### 3. `real_gui_run` 의 쓸기 실험들이 무엇을 답했는지 적는다

16.3 GB 가 파라미터 쓸기다 — `fit12/20/24/30`, `seeing_sweep`, `psf_ast_scale_off`,
`group25`. **이 실험들은 질문에 답하려고 만든 것이고, 답이 어딘가에 적혀 있으면
프레임은 안 들고 있어도 된다.** 적혀 있지 않으면 16 GB 를 들고 있으면서도 무엇을
알아냈는지 모르는 상태다.

`apex journal <폴더> --note "이 폴더는 무엇을 물었고 답이 무엇이었나"` 가 그 자리다.

---

## 다시 재는 법

레포 안 도구는 느리다(E 드라이브 순회에 15 분 넘게 걸린다). 위 수치는
PowerShell 로 쟀다.

```powershell
$o = [System.IO.EnumerationOptions]::new(); $o.RecurseSubdirectories=$true
$o.IgnoreInaccessible=$true; $o.AttributesToSkip=[System.IO.FileAttributes]::ReparsePoint
$agg = @{}
foreach ($r in @("E:\APEX_validation\reprocess","E:\APEX_validation_output")) {
  foreach ($f in [System.IO.Directory]::EnumerateFiles($r,'*',$o)) {
    $fi = New-Object System.IO.FileInfo $f
    $e = $fi.Extension.ToLower(); if (-not $agg.ContainsKey($e)) { $agg[$e] = @{n=0;b=0L} }
    $agg[$e].n++; $agg[$e].b += $fi.Length
  }
}
$agg.GetEnumerator() | ForEach-Object {
  [pscustomobject]@{ 확장자=$_.Key; GB=[math]::Round($_.Value.b/1GB,2); 파일수=$_.Value.n }
} | Sort-Object GB -Descending | Format-Table -AutoSize
```

`scripts/disk_survey.py` 는 같은 것을 재면서 **삭제 금지 목록을 함께 만든다** —
git 이 추적하는 파일과 레포가 이름으로 적어 둔 파일. 2026-08-20 의 사고 둘이
정확히 그 두 목록이 없어서 났다(F-061·F-062). 느리지만 **실제로 지우기 전에는
반드시 돌려야 한다.**
