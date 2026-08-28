# 의존 구조 대조 — §2.1 「돌리는 데 필요한 것」의 근거

2026-08-28 조사. **시간을 재는 대신 무엇을 따로 깔아야 하는지를 센다** — 사람이
필요 없고, 공개 문서에서 확인되며, 누구나 다시 셀 수 있다.

셀마다 출처를 달았다. **확인 못 한 칸은 「미확인」이라고 적었지 추측으로 채우지
않았다.**

---

## 표 — APEX · IRAF · AstroImageJ

| | APEX | IRAF (community 2.18.1) | AstroImageJ 6.0.10 |
|---|---|---|---|
| **사용자가 따로 설치할 프로그램** | 없음 — 설치본이 Python 런타임과 과학 패키지를 품는다 | **여럿.** 빌드에 gcc·make·flex·bison·zlib dev·readline dev, 실행에 readline 또는 libedit·zlib. 화면 표시는 X 서버 필요(맥은 XQuartz) | 없음 — 설치본이 Java 런타임(25)을 품는다 |
| **이미지 표시 도구** | 내장 | **따로 깐다.** X11IRAF 의 ximtool 이 딸려 오지만 *"most people install SAOImageDS9"* — DS9 은 번들이 아니다 | 내장(ImageJ 계열) |
| **공식 지원 운영체제** | Windows 설치본 배포. `pyproject.toml` 분류자에 Windows·Linux·macOS | Linux · macOS · FreeBSD · GNU HURD. **Windows 네이티브 미지원** — *"Cygwin and big endian architectures like macosx/ppc are not supported anymore"* | Windows · Linux(x64·Arm) · macOS(Intel·Arm) |
| **의존 버전을 못 박았는가** | **둘 다 있다.** `requirements.txt` 는 범위 21 개(설치가 패치 릴리스로 깨지지 않게), `requirements-lock.txt` 는 **91 개 전부 `==`** — 논문 수치를 낸 그 조합 | 미확인 | 미확인 (런타임을 함께 갱신한다고만 밝힘) |
| **인터넷이 필요한 지점** | 내장 측성 솔버의 첫 실행 — 관측 영역의 Gaia 목록을 받는다 (`apex/analysis/wcs_solve.py:2776`) | 미확인 | 미확인 |
| **배포 형태** | PyInstaller 묶음 + Inno Setup 설치본, 그리고 무설치 zip (`deploy/apex_windows.spec` · `installer.iss`) | 소스 빌드 또는 배포판 패키지 | 운영체제별 설치본 |

---

## 이 표를 읽는 법 — 과장하지 않으려면

**첫째, AstroImageJ 도 런타임을 품는다.** 6 판이 Java 25 를 함께 배포하고 갱신도
같이 한다. 그러므로 **「따로 깔 게 없다」는 APEX 만의 성질이 아니다.** 이 축에서
다른 것은 IRAF 하나다. 표를 「우리가 제일 쉽다」로 읽으면 틀리고, 심사에서 바로
걸린다.

**둘째, 런타임을 싼다는 것과 의존이 적다는 것은 다른 말이다.** APEX 안쪽에는
astropy · photutils · SEP · scipy · numpy · matplotlib · PyQt5 를 비롯해 91 개가
있다(`requirements-lock.txt`). **사용자가 깔 것이 없다는 뜻이지 의존이 적다는
뜻이 아니다.** 두 말을 섞어 쓰지 않는다.

**셋째, IRAF 의 Windows 문제는 인상이 아니라 문서에 적힌 사실이다.** 커뮤니티
배포판이 지원 목록에서 Windows 를 빼고 Cygwin 지원 종료를 명시한다. 저자 경험담을
근거로 삼을 필요가 없다.

**넷째, 버전 고정은 우리가 내세울 수 있는 자리다.** 범위와 잠금 파일을 둘 다 두는
구성이고, 잠금 파일이 있는 이유가 구체적이다 — 2026-08-17 에 **동일한 코드의 두
실행 사이에서 PSF 측정 22,305 개 중 다섯이 4.8e-05 mag 움직였고, 원인이
scipy 1.18.0 대 1.17.1 이었다.** 다섯 다 이미 불량으로 걸러진 별이었고 발표 그림에
들어가지 않았지만, 다음번이 그렇게 얌전하리라는 보장이 없어서 잠금 파일을 만들었다
(D-011). **재현 절차도 파일 머리말에 적혀 있다.**

---

## 출처

- IRAF — [IRAF Community Distribution 설치 안내](https://iraf-community.github.io/install.html) (2026-08-28 열람)
- AstroImageJ — [astroimagej.com](https://astroimagej.com/) 및 [6.0.0.00 릴리스 기록](https://astroimagej.com/releases/60000/) (2026-08-28 열람)
- APEX — 이 저장소의 `requirements.txt` · `requirements-lock.txt` · `pyproject.toml` ·
  `deploy/build_release.bat` · `deploy/apex_windows.spec` · `deploy/installer.iss`

## 투고 전에 채울 칸

- IRAF 와 AstroImageJ 의 **의존 버전 고정 여부** — 각 배포 방식을 더 봐야 한다
- 두 도구의 **인터넷이 필요한 지점** — AstroImageJ 는 astrometry.net 웹 API 를
  쓰는 것으로 알려져 있으나 **확인 전에는 적지 않는다**
