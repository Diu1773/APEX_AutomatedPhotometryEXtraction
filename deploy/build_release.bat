@chcp 65001 >nul
@echo off
setlocal EnableExtensions EnableDelayedExpansion

for %%I in ("%~dp0..") do set "ROOT_DIR=%%~fI"
set "DEPLOY_DIR=%ROOT_DIR%\deploy"
set "BUILD_DIR=%ROOT_DIR%\build"
set "DIST_DIR=%ROOT_DIR%\dist"
set "RELEASE_ROOT=%ROOT_DIR%\release"
set "SETUP_DIR=%RELEASE_ROOT%\Setup"
set "VERSION_FILE=%DEPLOY_DIR%\version.txt"
set "SPEC_FILE=%DEPLOY_DIR%\apex_windows.spec"
set "ISS_FILE=%DEPLOY_DIR%\installer.iss"
set "VENV_DIR=%ROOT_DIR%\.venv-deploy"
set "VENV_PY=%VENV_DIR%\Scripts\python.exe"

if not exist "%VERSION_FILE%" (
    >"%VERSION_FILE%" echo 0.1.0
)
set /p APP_VERSION=<"%VERSION_FILE%"
if "%APP_VERSION%"=="" set "APP_VERSION=0.1.0"

echo.
echo ============================================================
echo   APEX Release Builder
echo   Version: %APP_VERSION%
echo ============================================================
echo.

set "PYTHON_CMD="
where py >nul 2>nul
if not errorlevel 1 (
    set "PYTHON_CMD=py -3"
) else (
    where python >nul 2>nul
    if not errorlevel 1 (
        set "PYTHON_CMD=python"
    ) else (
        where python3 >nul 2>nul
        if not errorlevel 1 set "PYTHON_CMD=python3"
    )
)
if not defined PYTHON_CMD (
    echo [ERROR] Python was not found.
    echo         Install Python 3.10+ for Windows first.
    exit /b 1
)
%PYTHON_CMD% -c "import sys; raise SystemExit(0 if sys.version_info >= (3, 10) else 1)" >nul 2>nul
if errorlevel 1 (
    echo [ERROR] Python 3.10+ runtime was not found.
    echo         The Python launcher may exist, but no usable Python 3.10+ install is registered.
    exit /b 1
)
echo [OK] Python: %PYTHON_CMD%

set "ISCC_EXE="
where ISCC.exe >nul 2>nul
if not errorlevel 1 (
    for /f "delims=" %%I in ('where ISCC.exe') do (
        if not defined ISCC_EXE set "ISCC_EXE=%%I"
    )
)

if not defined ISCC_EXE if exist "%ProgramFiles(x86)%\Inno Setup 6\ISCC.exe" (
    set "ISCC_EXE=%ProgramFiles(x86)%\Inno Setup 6\ISCC.exe"
)
if not defined ISCC_EXE if exist "%ProgramFiles%\Inno Setup 6\ISCC.exe" (
    set "ISCC_EXE=%ProgramFiles%\Inno Setup 6\ISCC.exe"
)
if not defined ISCC_EXE (
    echo [ERROR] Inno Setup 6 compiler was not found.
    echo         Install Inno Setup 6 and try again:
    echo         https://jrsoftware.org/isdl.php
    exit /b 1
)
echo [OK] Inno Setup: %ISCC_EXE%

if not exist "%RELEASE_ROOT%" mkdir "%RELEASE_ROOT%"
if exist "%SETUP_DIR%" rmdir /s /q "%SETUP_DIR%"
mkdir "%SETUP_DIR%"

echo [1/10] Preparing isolated build environment...
if not exist "%VENV_PY%" (
    %PYTHON_CMD% -m venv "%VENV_DIR%"
    if errorlevel 1 exit /b 1
)
"%VENV_PY%" -m pip install --upgrade pip setuptools wheel
if errorlevel 1 exit /b 1
"%VENV_PY%" -m pip install -r "%ROOT_DIR%\requirements.txt" -r "%ROOT_DIR%\requirements-build.txt"
if errorlevel 1 exit /b 1

echo [2/10] Running release source preflight...
"%VENV_PY%" "%DEPLOY_DIR%\verify_release.py" --project-root "%ROOT_DIR%" --source-only
if errorlevel 1 exit /b 1

echo [3/10] Preparing Windows icon...
"%VENV_PY%" "%DEPLOY_DIR%\make_icon.py" --svg "%ROOT_DIR%\apex\resources\logo_base.svg" --output "%ROOT_DIR%\apex\resources\apex.ico"
if errorlevel 1 exit /b 1

echo [4/10] Cleaning old build outputs...
if exist "%BUILD_DIR%" rmdir /s /q "%BUILD_DIR%"
if exist "%DIST_DIR%" rmdir /s /q "%DIST_DIR%"
mkdir "%BUILD_DIR%"
mkdir "%DIST_DIR%"

echo [5/10] Validating Python sources...
pushd "%ROOT_DIR%"
"%VENV_PY%" -m compileall apex main.py scripts deploy
if errorlevel 1 (
    popd
    exit /b 1
)
"%VENV_PY%" -m pytest tests
if errorlevel 1 (
    popd
    exit /b 1
)

echo [6/10] Building PyInstaller app bundle...
"%VENV_PY%" -m PyInstaller --noconfirm --clean "%SPEC_FILE%"
if errorlevel 1 (
    popd
    exit /b 1
)
popd

if not exist "%DIST_DIR%\APEX\APEX.exe" (
    echo [ERROR] Build finished without dist\APEX\APEX.exe
    exit /b 1
)

echo [7/10] Smoke testing PyInstaller app bundle...
start "" /wait "%DIST_DIR%\APEX\APEX.exe" --smoke
if errorlevel 1 exit /b 1

echo [8/10] Creating portable ZIP...
"%VENV_PY%" "%DEPLOY_DIR%\make_portable_zip.py" --source "%DIST_DIR%\APEX" --output "%SETUP_DIR%\APEX-Portable-%APP_VERSION%-x64.zip"
if errorlevel 1 exit /b 1

echo [9/10] Building setup installer...
"%ISCC_EXE%" ^
  "/DMyAppVersion=%APP_VERSION%" ^
  "/DSourceDir=%DIST_DIR%\APEX" ^
  "/DOutputDir=%SETUP_DIR%" ^
  "/DIconFile=%ROOT_DIR%\apex\resources\apex.ico" ^
  "%ISS_FILE%"
if errorlevel 1 exit /b 1

echo [10/10] Verifying release bundle...
"%VENV_PY%" "%DEPLOY_DIR%\verify_release.py" --project-root "%ROOT_DIR%"
if errorlevel 1 exit /b 1

echo.
echo Release output:
echo   %SETUP_DIR%
echo.
dir /b "%SETUP_DIR%"
exit /b 0
