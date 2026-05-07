$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$Root = Resolve-Path (Join-Path $PSScriptRoot "..")
$env:APEX_NO_PAUSE = "1"

cmd /c "`"$Root\build.bat`""
if ($LASTEXITCODE -ne 0) {
    throw "APEX Windows release build failed with exit code $LASTEXITCODE."
}
