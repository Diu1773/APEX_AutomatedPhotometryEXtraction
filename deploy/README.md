# APEX Windows Deployment

APEX의 공식 release pipeline은 Windows x64용 PyInstaller `onedir` bundle,
portable ZIP, Inno Setup installer를 생성합니다.

## Outputs

`build.bat` 실행이 성공하면 `release\Setup\`에 다음 파일이 생성됩니다.

```text
setup-APEX-<version>.exe
APEX-Portable-<version>-x64.zip
```

버전은 `deploy/version.txt`에서 읽습니다. Installer는 기본적으로
`%LOCALAPPDATA%\Programs\APEX`에 설치되며 관리자 권한을 요구하지
않습니다.

## Build Requirements

- Windows x64
- Python 3.10 이상
- Inno Setup 6
- 네트워크 접근: isolated build venv의 Python dependency 설치

Python runtime dependency는 `requirements.txt`, build-only dependency는
`requirements-build.txt`에서 설치됩니다.

## Local Build

저장소 루트에서 실행합니다.

```powershell
.\build.bat
```

pause 없이 자동화하려면:

```powershell
$env:APEX_NO_PAUSE = "1"
.\build.bat
```

`deploy/build_release.bat`를 직접 실행해도 동일한 pipeline을 사용합니다.

## Build Pipeline

1. `.venv-deploy` isolated environment 준비
2. `verify_release.py --source-only` preflight
3. SVG logo에서 Windows ICO 생성
4. 이전 `build/`, `dist/`, `release/Setup/` 정리
5. `compileall`과 전체 pytest 실행
6. `apex_windows.spec`로 PyInstaller bundle 생성
7. `APEX.exe --smoke` 실행 및 stdout/stderr 수집
8. portable ZIP 생성
9. Inno Setup installer 생성
10. 최종 파일과 bundle 구조 검증

`--smoke`는 GUI를 열지 않고 다음 범주를 import합니다.

- shared/CMD workflow의 build-critical modules
- Step 1 target resolution
- `astroquery.gaia`, `astroquery.simbad`, TAP support
- SSL CA bundle support
- Tools-menu registry의 모든 module

## Packaged Configuration

배포본에는 `parameters.example.toml`이 포함됩니다. 첫 실행 시
실행 파일 옆에 `parameters.toml`이 없으면 예제 파일을 복사합니다.

Portable bundle을 read-only 위치에서 실행하면 설정 파일 생성이 실패할
수 있으므로 사용자가 쓸 수 있는 폴더에 압축을 해제해야 합니다.

## External WCS Dependencies

내장 Python solver는 bundle에 포함되지만 Gaia catalog cache 또는
네트워크 접근이 필요합니다. 다음 외부 solver는 포함되지 않습니다.

- ASTAP executable과 D80/D50 star database
- local/WSL astrometry.net `solve-field`
- astrometry.net index files

공식 설치 자료:

- <https://www.hnsky.org/astap.htm>
- <https://sourceforge.net/projects/astap-program/files/star_databases/>
- <https://astrometry.net/doc/readme.html>
- <https://data.astrometry.net/>

## PyInstaller Data Rules

일부 scientific package는 Python import 외에 package data를 런타임에
읽습니다. `deploy/apex_windows.spec`의 collection 목록을 변경할 때는
반드시 frozen smoke test로 확인해야 합니다.

현재 중요한 예:

- `astroquery`와 TAP support
- `pyvo`의 SAMP data
- `certifi` CA bundle
- `astropy`, `photutils`, `matplotlib` data/hooks
- `apex/resources`
- `parameters.example.toml`

소스 환경에서 import가 성공해도 package data가 누락되면 frozen
application에서만 실패할 수 있습니다.

## CI

`.github/workflows/windows-build.yml`은 다음 경우 실행됩니다.

- `main` branch push
- pull request
- manual `workflow_dispatch`

CI는 Python 3.11과 `windows-latest`에서 전체 release build를 수행하고
`release/Setup/*`을 `APEX-Windows` artifact로 업로드합니다.

`.github/workflows/tests.yml`은 Ubuntu/Python 3.11에서 syntax check와
pytest를 실행합니다. Windows packaging 성공 여부는 별도 Windows Build
workflow가 최종 기준입니다.

## Release Checklist

1. `deploy/version.txt`를 확인한다.
2. `python -m compileall apex main.py scripts deploy`를 통과시킨다.
3. `python -m pytest tests`를 통과시킨다.
4. `.\build.bat`를 실행한다.
5. source preflight와 `APEX.exe --smoke` 성공을 확인한다.
6. installer와 portable ZIP이 모두 생성됐는지 확인한다.
7. clean Windows 사용자 계정 또는 VM에서 installer를 실행한다.
8. CMD/LC 모드가 열리고 `parameters.toml`이 생성되는지 확인한다.
9. Internal WCS와 최소 하나의 실제 FITS workflow를 점검한다.
10. 외부 solver 사용 release라면 ASTAP/astrometry.net 경로를 점검한다.

## Troubleshooting

### Wrong Python selected

Build script는 현재 shell의 `python`을 우선 사용하고, 없으면 `py -3`,
`python3` 순으로 탐색합니다. CI와 같은 결과가 필요하면 Python 3.11
environment에서 실행하십시오.

### Smoke test fails only after packaging

`build\smoke-stdout.log`와 `build\smoke-stderr.log`를 확인하십시오.
대부분 hidden import 또는 package data 누락입니다.

### Installer is missing

`release\Setup\setup-APEX-<version>.exe`가 없다면 Inno Setup 단계와
`ISCC.exe` 탐색 로그를 확인하십시오.
