@chcp 65001 >nul
@echo off
setlocal EnableExtensions

REM Development launcher only. Release builds use build.bat.

set "ROOT_DIR=%~dp0"
if "%ROOT_DIR:~-1%"=="\" set "ROOT_DIR=%ROOT_DIR:~0,-1%"
cd /d "%ROOT_DIR%"

set "PYTHON_EXE="
set "PYTHON_ARGS="

if defined APEX_PYTHON (
    set "PYTHON_EXE=%APEX_PYTHON%"
    goto :python_found
)

for /f "usebackq delims=" %%I in (`powershell -NoProfile -ExecutionPolicy Bypass -Command "$ErrorActionPreference='Stop'; $root=(Resolve-Path '.').Path; $settings=Join-Path $root '.vscode/settings.json'; if (Test-Path -LiteralPath $settings) { $json=Get-Content -Raw -LiteralPath $settings | ConvertFrom-Json; $p=[string]$json.'python.defaultInterpreterPath'; if ($p) { $p=$p.Replace('${workspaceFolder}', $root); $p=[Environment]::ExpandEnvironmentVariables($p); if (Test-Path -LiteralPath $p) { (Resolve-Path -LiteralPath $p).Path } } }" 2^>nul`) do (
    if not defined PYTHON_EXE set "PYTHON_EXE=%%I"
)
if defined PYTHON_EXE goto :python_found

if exist "%ROOT_DIR%\.venv\Scripts\python.exe" (
    set "PYTHON_EXE=%ROOT_DIR%\.venv\Scripts\python.exe"
    goto :python_found
)

if exist "%ROOT_DIR%\.venv-deploy\Scripts\python.exe" (
    set "PYTHON_EXE=%ROOT_DIR%\.venv-deploy\Scripts\python.exe"
    goto :python_found
)

where py >nul 2>nul
if not errorlevel 1 (
    set "PYTHON_EXE=py"
    set "PYTHON_ARGS=-3"
    goto :python_found
)

where python >nul 2>nul
if not errorlevel 1 (
    set "PYTHON_EXE=python"
    goto :python_found
)

echo [ERROR] Python was not found.
echo         Set APEX_PYTHON, configure .vscode\settings.json, or create .venv.
exit /b 1

:python_found
if /I not "%PYTHON_EXE%"=="py" if /I not "%PYTHON_EXE%"=="python" if not exist "%PYTHON_EXE%" (
    echo [ERROR] Configured Python does not exist:
    echo         %PYTHON_EXE%
    exit /b 1
)

set "ENTRY=main.py"
set "EXTRA_ARGS="
set "RUN_ONCE=0"

if defined APEX_RUN_ONCE set "RUN_ONCE=1"
if defined APEX_NO_RERUN set "RUN_ONCE=1"

if /I "%~1"=="cmd" (
    set "ENTRY=apex\cmd\main.py"
    shift
) else if /I "%~1"=="lc" (
    set "ENTRY=apex\lightcurve\main.py"
    shift
) else if /I "%~1"=="lightcurve" (
    set "ENTRY=apex\lightcurve\main.py"
    shift
) else if /I "%~1"=="launcher" (
    set "ENTRY=main.py"
    shift
) else if /I "%~1"=="smoke" (
    set "ENTRY=main.py"
    set "EXTRA_ARGS=--smoke"
    set "RUN_ONCE=1"
    shift
)

set "APP_ARGS="
:collect_args
if "%~1"=="" goto :run_app
if /I "%~1"=="--once" (
    set "RUN_ONCE=1"
    shift
    goto :collect_args
)
if /I "%~1"=="once" (
    set "RUN_ONCE=1"
    shift
    goto :collect_args
)
if /I "%~1"=="--no-rerun" (
    set "RUN_ONCE=1"
    shift
    goto :collect_args
)
set APP_ARGS=%APP_ARGS% "%~1"
shift
goto :collect_args

:run_app
if "%RUN_ONCE%"=="1" (
    echo [APEX] Python: %PYTHON_EXE% %PYTHON_ARGS%
    echo [APEX] Entry : %ENTRY% %EXTRA_ARGS% %APP_ARGS%
    "%PYTHON_EXE%" %PYTHON_ARGS% "%ENTRY%" %EXTRA_ARGS% %APP_ARGS%
    exit /b %ERRORLEVEL%
)

set "APEX_DEV_ROOT=%ROOT_DIR%"
set "APEX_DEV_PYTHON_EXE=%PYTHON_EXE%"
set "APEX_DEV_PYTHON_ARGS=%PYTHON_ARGS%"
set "APEX_DEV_ENTRY=%ENTRY%"
set "APEX_DEV_EXTRA_ARGS=%EXTRA_ARGS%"
set "APEX_DEV_APP_ARGS=%APP_ARGS%"

powershell -NoProfile -ExecutionPolicy Bypass -File "%ROOT_DIR%\scripts\dev_run_hotkeys.ps1"
exit /b %ERRORLEVEL%
