@chcp 65001 >nul
@echo off
setlocal EnableExtensions

REM Canonical Windows release entry point.
REM This creates release\Setup\setup.exe plus a portable zip.

call "%~dp0deploy\build_release.bat"
set "EXIT_CODE=%ERRORLEVEL%"

echo.
if not "%EXIT_CODE%"=="0" (
    echo Build failed with exit code %EXIT_CODE%.
) else (
    echo Build finished successfully.
)

if not "%APEX_NO_PAUSE%"=="1" pause
exit /b %EXIT_CODE%
