# 7. 문제 해결 / FAQ

[← 파라미터 레퍼런스](06-parameters-reference.md) · [매뉴얼 목차](index.md)

증상별로 빠르게 찾아보세요. 대부분은 "앞 단계를 다시 돌리거나 파라미터 하나를 조정"하면 해결됩니다.

---

## 실행·시작

**Q. 처음 실행했더니 설정이 비어 있습니다.**
처음 실행 시 `parameters.example.toml`을 복사해 `parameters.toml`을 만듭니다. 이 파일에서
`[io] data_dir`/`result_dir`, `[target]` 좌표, `[instrument]` 장비 값을 채우세요. → [시작하기 1.5](01-getting-started.md#15-parameterstoml-기초)

**Q. 모드별 진입점(`apex/cmd/main.py`)을 직접 실행했더니 설정이 없습니다.**
루트 런처(`python main.py`)를 한 번 실행해 `parameters.toml`을 먼저 만든 뒤 모드를 직접 실행하세요.

**Q. 환경이 제대로 깔렸는지 확인하고 싶습니다.**
```powershell
apex doctor            # 파이썬·의존성·외부 솔버 점검
apex doctor --network  # 네트워크(Gaia/SIMBAD)까지 점검
```

---

## 단계 잠김·이동

**Q. 다음 단계가 회색(🔒 Locked)이라 안 열립니다.**
앞 단계의 산출물이 있어야 열립니다. 위에서부터 순서대로 끝내세요. 단계 막대 색 의미:
🟢 완료 · 🔵 진행 가능 · ⚪ 잠김. → [시작하기 1.2](01-getting-started.md#12-런처-화면과-모드-선택)

**Q. `Next Step →`가 비활성(회색)입니다.**
그 단계가 아직 "유효"하지 않습니다. 필요한 동작(예: `Run …`, 타겟 좌표 확정, 별 검출)을
끝내면 자동으로 켜집니다. 각 단계 문서의 "따라하기" 마지막 항목을 확인하세요.

**Q. 무엇을 해야 할지 모르겠습니다.**
오른쪽 위 **`⎘ 가이드`** 버튼을 누르면 그 단계 요약 도움말이 뜹니다.

---

## Step 1~2 (파일·크롭)

**Q. 폴더를 바꿨는데 파일 목록이 그대로입니다.**
**`Rescan Files`** 를 눌러야 반영됩니다.

**Q. SIMBAD 좌표 조회가 안 됩니다.**
인터넷이 필요합니다. 안 되면 **`수동 입력`** 으로 RA/Dec를 직접 넣거나, 헤더 좌표가 있으면
**`Use Header RA/Dec`** 를 쓰세요.

**Q. 크롭했더니 "Re-crop Warning"이 뜹니다.**
Step 4 이후에 다시 크롭하면 픽셀 좌표가 바뀌므로 Step 4부터 다시 돌려야 합니다.

---

## Step 3~4 (스카이·검출)

**Q. 별이 너무 적게/많이 검출됩니다.**
**`Detection Parameters`** → **Detection Sigma** 를 낮추면(예 3.2→2.5) 더 많이, 높이면 더 적게
검출됩니다. 혼잡장은 `Crowded` 프리셋, 희미한 장은 `Faint` 프리셋을 `Apply Preset` 하세요.

**Q. 자동 QC가 제대로 안 됩니다.**
자동 QC(robust z)는 **필터당 프레임 ≥ 10개**가 있어야 합니다. 프레임이 적으면 수동으로
의심 프레임을 `Exclude` 하세요.

**Q. QC에서 elong이 전부 N/A입니다.**
검출 캐시가 오래됐습니다. **`Clear Detection Cache`** 후 다시 `Run Detection`.

---

## Step 5 (WCS)

**Q. WCS가 안 풀립니다.**
- 내장 솔버는 정확한 **픽셀 스케일**(장비 값)·대상 좌표·충분한 검출·Gaia 카탈로그가 필요합니다.
- ASTAP는 실행파일 경로(`wcs.astap_exe`)와 별 DB(D80/D50)가, astrometry.net은 `solve-field`와
  인덱스 파일이 **별도 설치**돼 있어야 합니다.
- 솔버를 바꿔 다시 돌려도 성공한 프레임만 갱신되니, 탭을 바꿔 재시도해 보세요(체인: ASTAP → astrometry.net → Internal).

**Q. "No frames remain after Step 4 QC filtering"**
Step 4에서 프레임을 너무 많이 제외했습니다. Step 4 QC의 제외를 줄이세요.

**Q. Gaia 조회가 자꾸 멈춥니다(타임아웃).**
`gaia.wcs_mag_max`(기본 18)를 너무 어둡게 잡으면 쿼리가 무거워집니다. 큰 시야에서는 등급
상한을 낮추세요. `gaia.hard_deadline_s`로 한 쿼리의 최대 시간을 강제할 수 있습니다.

---

## Step 6~7 (마스터·측광)

**Q. "Run Source Detection first." / "No Master Catalog"**
앞 단계(4 검출, 5 WCS, 6 마스터)를 먼저 끝내세요.

**Q. 측광 결과를 다시 안 만들고 그대로입니다.**
**`Use existing output if complete`** 체크가 켜져 있으면 입력·파라미터가 같을 때 기존 결과를
재사용합니다. 강제로 다시 돌리려면 체크를 끄세요.

**Q. 조리개 보정(apcorr)이 이상합니다.**
밝고 고립된 별(SNR ≥ `photometry.apcorr.min_snr` = 40)이 충분해야 합니다. `Apcorr` 탭에서
현재 조리개(빨강)가 최적 반경 `r_opt`(보라)에 가까운지 보고, `apcorr_summary.csv`의 reject
수를 확인하세요.

**Q. Stats 탭에 REVIEW/CHECK가 많습니다.**
센터링이 불안정합니다. WCS 품질을 점검하거나 **`Match / recenter limit`** 를 조정하세요.

---

## CMD 모드 (8~12)

**Q. PSF가 꼭 필요한가요?**
아니요. **선택 사항**입니다. 넓은 시야/산개성단은 **`Skip PSF →`** 로 넘어가도 됩니다(이후 단계는
Step 7 강제 측광을 사용).

**Q. "Not enough Gaia matches for calibration." (Step 10)**
Gaia 매칭이 `match.min_gaia_matches`(기본 10)보다 적습니다. WCS 품질과 `gaia.wcs_mag_max`를
점검하세요.

**Q. "Isochrone Filter Mismatch" (Step 12)**
관측 밴드와 아이소크론 파일의 측광 시스템이 다릅니다. 예: Johnson B-V는 Johnson/Bessell
파일이라야 합니다(SDSS ugriz 파일로는 안 됨).

**Q. 아이소크론 나이는 맞는데 금속함량이 이상합니다.**
gri/BVR 색만으로는 금속함량이 거의 제약되지 않습니다(사실상 "나이 측정기"). **청색/자외선
밴드**(u-g, U-B)를 쓰거나 MCMC에서 **[M/H] 사전값**을 주세요. → [CMD §Step12](03-cmd-mode.md#step-12--isochrone-model-아이소크론-모델)

**Q. 희미한 별이 CMD에서 너무 파랗게 쏠립니다.**
**SNR-20 게이트**(기본 ON)를 끄지 마세요. 저SNR 점은 색을 왜곡합니다.

---

## LC 모드 (8~11)

**Q. 광도곡선이 들쭉날쭉합니다.**
**변광하는 비교성**이 섞였을 가능성이 큽니다(차등등급 = 타겟 − 비교성평균). Step 9
`Comparison QC`에서 `Run QC` → `Auto Use` 후 **각 비교성 미리보기**를 눈으로 확인하고
흔들리는 비교성의 `Use`를 끄세요. `Gaia Var`·`SIMBAD` 열로 변광 여부를 미리 거를 수 있습니다.

**Q. 한 밤이 두 밤으로 갈립니다(또는 둘이 합쳐집니다).**
LC Step 1에서 **`Night gap:`** 시간을 키우거나 줄이고 **`Rescan Files`** 하세요.

**Q. 여러 밤을 합쳤더니 계단처럼 어긋납니다.**
Step 10 **디트렌드**에서 밤별 영점(ZP₀)을 보정합니다. `Fit && Apply (저장)` 후 **RMS전→RMS후**
가 줄었는지, 체크별 곡선이 평평해지는지 확인하세요.

**Q. 주기가 진짜인지 모르겠습니다.**
LS·PDM·BLS **여러 방법의 결과가 일치**하면 신뢰도가 높습니다. 1일 별칭은 `Alias?` 열·주황
보조선으로 표시됩니다. 짧거나 불균일한 자료는 **Bootstrap FAP**로 유의성을 보세요.

**Q. 더 정밀한 변광/트랜짓 분석을 하고 싶습니다.**
Step 11은 빠른 스캔용입니다. O-C·다중모드·트랜짓·식쌍성 정밀 피팅은 **Tools 메뉴**에
있습니다. → [5. 도구](05-tools.md)

---

## 도구 (Tools)

**Q. Tools 메뉴에 원하는 도구가 안 보입니다.**
도구는 모드별로 다릅니다. 예: Gaia 3D 뷰어·성단 구조는 **CMD 전용**, 변광성·트랜짓·식쌍성·
다중밤 병합은 **LC 전용**입니다. 모드를 맞춰 실행하세요.

**Q. IRAF 도구가 안 돕니다.**
Windows에서는 **WSL + PyRAF/IRAF/DAOPHOT**가 필요합니다. **`Check Environment`** 로 먼저
점검하세요.

**Q. 트랜짓/변광성 도구에서 의존성 오류가 납니다.**
`batman-package`·`emcee`는 선택 의존성입니다. 해당 기능을 쓰려면 추가 설치하세요.

**Q. Gaia 3D 뷰어 MP4 내보내기가 안 됩니다.**
MP4는 `ffmpeg`가 PATH에 있어야 합니다. 없으면 GIF로 내보내세요.

---

## 성능·재실행

**Q. 너무 느립니다.**
- 워커 수는 `parallel.max_workers`(0=자동)로 조절합니다. 0이면 코어 수에 맞춰 자동 결정.
- Gaia 등급 상한을 낮추면 WCS·매칭이 빨라집니다.
- 캐시 재사용 체크박스(`Use detection cache`, `Use existing output if complete`)를 켜 두면
  바뀌지 않은 부분을 건너뜁니다.

**Q. 파라미터를 바꿨는데 어디부터 다시 돌려야 하나요?**
변경 범주별 재실행 시작 단계는 [시작하기 1.5 표](01-getting-started.md#설정을-바꾸면-어디서부터-다시-돌려야-하나)를 보세요.

---

## 검증·재현

```powershell
python -m compileall apex main.py     # 문법 점검
python -m pytest tests                # 테스트
python main.py --smoke                # GUI 없이 주요 모듈 import 점검
```

과학 검증(합성·실데이터)은 `validation/`·`benchmark/` 폴더와 `docs/validation*.md`를 참고하세요.

---

[← 파라미터 레퍼런스](06-parameters-reference.md) · [매뉴얼 목차](index.md)
