<!--
  청람천문대 구글사이트(Google Sites)용 붙여넣기 콘텐츠
  ------------------------------------------------------------
  이 파일은 "웹사이트에 들어갈 내용"을 정리한 것입니다. 구글사이트 편집기에서
  아래 블록을 그대로 복사해 붙여넣으세요. (이 파일 자체는 사이트에 올라가지 않습니다.)

  먼저 알아둘 링크:
    · 다운로드(항상 최신):  https://github.com/Diu1773/APEX_AutomatedPhotometryEXtraction/releases/latest
    · 사용 매뉴얼(웹):      https://diu1773.github.io/APEX_AutomatedPhotometryEXtraction/manual/
    · 프로젝트(소스):       https://github.com/Diu1773/APEX_AutomatedPhotometryEXtraction
-->

# 청람천문대 구글사이트 — APEX 배포 페이지 만들기

구글사이트(사이트) 편집기는 제가 직접 조작할 수 없어서, **붙여넣기만 하면 되는 내용**을
아래에 정리했습니다. 5분이면 페이지 하나가 완성됩니다.

세 가지 핵심 링크부터 기억하세요:

| 용도 | 주소 |
| --- | --- |
| **다운로드** (항상 최신 버전) | `https://github.com/Diu1773/APEX_AutomatedPhotometryEXtraction/releases/latest` |
| **사용 매뉴얼** (웹) | `https://diu1773.github.io/APEX_AutomatedPhotometryEXtraction/manual/` |
| **프로젝트/소스** | `https://github.com/Diu1773/APEX_AutomatedPhotometryEXtraction` |

> 다운로드 링크는 **항상 최신 릴리스**로 자동 연결되므로, 새 버전을 올려도 사이트는 고칠 필요가 없습니다.

---

## 가장 빠른 방법 — 3단계

1. 구글사이트 편집기에서 **삽입 → 버튼** 으로 **다운로드** 버튼을 만들고 위 *다운로드* 주소를 넣습니다.
2. **삽입 → 버튼**(또는 텍스트 링크)으로 **사용 매뉴얼** 버튼을 만들고 *매뉴얼* 주소를 넣습니다.
3. 아래 "소개 문단"과 "스크린샷"을 붙여넣습니다. 끝.

> 더 멋지게 하고 싶으면 아래 **"매뉴얼을 페이지 안에 끼워넣기(임베드)"** 를 보세요.

---

## 붙여넣기용 콘텐츠 (그대로 복사)

### ① 페이지 제목

```
APEX — 천체 측광 프로그램
Automated Photometry EXtraction
```

### ② 소개 문단

```
APEX는 청람천문대에서 만든 천체 측광(photometry) 데스크톱 프로그램입니다.
망원경으로 찍은 FITS 영상에서 별의 밝기를 측정해, 성단의 색-등급도(CMD)와
나이·거리, 또는 변광성·외계행성의 광도곡선과 변광 주기까지 단계별로 구합니다.

천문 측광을 처음 다루는 학생도 화면을 보며 순서대로 따라 하면 결과가 나오도록
설계했습니다. Windows에서 설치형 또는 무설치(포터블)로 바로 실행할 수 있습니다.
```

### ③ 다운로드 버튼 (버튼 2개 권장)

- 버튼 1 — 라벨: **`APEX 내려받기 (Windows)`** · 링크:
  `https://github.com/Diu1773/APEX_AutomatedPhotometryEXtraction/releases/latest`
- 버튼 2 — 라벨: **`사용 매뉴얼 보기`** · 링크:
  `https://diu1773.github.io/APEX_AutomatedPhotometryEXtraction/manual/`

> 다운로드 페이지(Releases)의 **Assets** 에서 두 파일 중 하나를 받습니다:
> `setup-APEX-<버전>.exe`(설치형, 관리자 권한 불필요) 또는
> `APEX-Portable-<버전>-x64.zip`(무설치 — 압축 풀고 `APEX.exe` 실행).

### ④ 주요 기능 (글머리표)

```
· FITS 파일·필터별 프레임 관리, sky·FWHM·포화 기반 품질 확인(QC)
· SEP/DAO 소스 검출과 마스터 카탈로그 구성
· 내장 Python WCS 솔버 + ASTAP/astrometry.net 지원, Gaia/SIMBAD 조회
· 강제 조리개 측광과 선택적 PSF 측광 (독립 소프트웨어와 ~3 mmag 일치 검증)
· CMD 영점·색항 보정과 PARSEC/BaSTI 아이소크론 피팅 (성단 나이·거리)
· 비교성 QC, 디트렌드, SysRem, 다중밤 병합
· Lomb-Scargle·PDM·BLS 주기 분석과 변광성·외계행성·식쌍성 전용 도구
```

### ⑤ 시스템 요구사항

```
· 운영체제: Windows 10/11 (x64)
· 설치 권한: 불필요 (사용자 영역 설치)
· (선택) 외부 WCS 솔버 ASTAP 또는 astrometry.net — 없어도 내장 솔버로 동작
· 소스로 직접 실행 시: Python 3.10 이상
```

### ⑥ 매뉴얼·문의 안내

```
· 사용 매뉴얼(웹): https://diu1773.github.io/APEX_AutomatedPhotometryEXtraction/manual/
· 소스/이슈 신고: https://github.com/Diu1773/APEX_AutomatedPhotometryEXtraction
· 라이선스: MIT (자유롭게 사용·배포 가능)
· 논문에 사용 시 인용: 저장소의 CITATION.cff 참고
```

---

## 스크린샷 — 어떤 걸 올릴까

구글사이트 **삽입 → 이미지** 로 올리세요. 파일은 프로젝트의 `ui_screenshots/` 폴더에 있습니다
(추천 4장):

| 파일 | 설명(캡션으로 사용) |
| --- | --- |
| `ui_screenshots/00_main_cmd.png` | APEX 메인 화면 — 12단계 워크플로 |
| `ui_screenshots/cmd_step11_cmd_plot.png` | 색-등급도(CMD) — 유효온도 색상으로 표시 |
| `ui_screenshots/cmd_step12_isochrone_model.png` | 아이소크론 피팅 — 성단 나이·거리 측정 |
| `ui_screenshots/lc_step11_period_analysis.png` | 광도곡선 주기 분석 (LS·PDM·BLS) |

> 더 많은 화면은 매뉴얼에 모두 들어 있으니, 사이트에는 대표 몇 장만 올리면 충분합니다.

---

## (선택) 매뉴얼을 페이지 안에 끼워넣기 (임베드)

매뉴얼 웹사이트를 구글사이트 안에 통째로 보여줄 수 있습니다.

1. 구글사이트 편집기 → **삽입 → 삽입(Embed) → 전체 페이지 URL** 선택
2. 아래 주소를 붙여넣기:
   ```
   https://diu1773.github.io/APEX_AutomatedPhotometryEXtraction/manual/
   ```
3. 크기를 페이지 폭에 맞게 조정.

> 임베드가 빈 화면으로 보이면(일부 사이트는 iframe을 막습니다), 임베드 대신 위 **버튼 링크**로
> 새 탭에서 열도록 하세요. 가장 안전한 방식은 "버튼 링크"입니다.

---

## 동작 원리 (관리가 거의 필요 없는 이유)

```
  [청람천문대 구글사이트]  ──(버튼/링크)──▶  [GitHub Releases]   ← 실제 다운로드 파일
          │                                  (항상 최신 = /releases/latest)
          └────(버튼/링크 또는 임베드)──▶  [GitHub Pages 매뉴얼]  ← 웹 매뉴얼
```

- 새 버전 배포: 개발 컴퓨터에서 `git tag v0.1.1 && git push --tags` → GitHub가 자동으로
  빌드하고 **새 릴리스를 게시**합니다. 구글사이트는 "최신"을 가리키므로 **고칠 필요 없음.**
- 매뉴얼 수정: `docs/manual/` 의 내용을 고쳐 push 하면 매뉴얼 웹사이트가 **자동 갱신**됩니다.

즉, 한 번 세팅해두면 구글사이트 쪽은 거의 손댈 일이 없습니다.
