$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$root = $env:APEX_DEV_ROOT
$pythonExe = $env:APEX_DEV_PYTHON_EXE
$pythonArgs = $env:APEX_DEV_PYTHON_ARGS
$entry = $env:APEX_DEV_ENTRY
$extraArgs = $env:APEX_DEV_EXTRA_ARGS
$appArgs = $env:APEX_DEV_APP_ARGS

if (-not $root) {
    $root = Resolve-Path (Join-Path $PSScriptRoot "..")
}
if (-not $pythonExe) {
    throw "APEX_DEV_PYTHON_EXE is not set."
}
if (-not $entry) {
    throw "APEX_DEV_ENTRY is not set."
}

Set-Location -LiteralPath $root

function Quote-CommandArg {
    param([string] $Value)
    if ([string]::IsNullOrWhiteSpace($Value)) {
        return ""
    }
    return '"' + ($Value -replace '"', '\"') + '"'
}

function Get-ArgumentLine {
    $parts = @(
        $pythonArgs,
        (Quote-CommandArg $entry),
        $extraArgs,
        $appArgs
    ) | Where-Object { -not [string]::IsNullOrWhiteSpace([string]$_) }
    return ($parts -join " ")
}

function Stop-ProcessTree {
    param([int] $ProcessId)

    try {
        $children = Get-CimInstance Win32_Process -Filter "ParentProcessId=$ProcessId"
        foreach ($child in $children) {
            Stop-ProcessTree -ProcessId ([int]$child.ProcessId)
        }
    } catch {
        # Best effort; the main process is still stopped below.
    }

    try {
        $proc = Get-Process -Id $ProcessId -ErrorAction SilentlyContinue
        if ($proc -and -not $proc.HasExited) {
            Stop-Process -Id $ProcessId -Force -ErrorAction SilentlyContinue
        }
    } catch {
        # Process may have exited between checks.
    }
}

function Start-ApexProcess {
    $argumentLine = Get-ArgumentLine
    Write-Host "[APEX] Python: $pythonExe $pythonArgs"
    Write-Host "[APEX] Entry : $entry $extraArgs $appArgs"
    Write-Host "[APEX] Hotkeys: Q=quit, R=rerun"
    return Start-Process -FilePath $pythonExe `
        -ArgumentList $argumentLine `
        -WorkingDirectory $root `
        -NoNewWindow `
        -PassThru
}

$lastExitCode = 0
$apexProcess = $null

try {
    $apexProcess = Start-ApexProcess

    while ($true) {
        if ($apexProcess.HasExited) {
            $lastExitCode = $apexProcess.ExitCode
            Write-Host ""
            Write-Host "[APEX] App exited with code $lastExitCode."
            Write-Host "[APEX] Q=quit, R=rerun"

            while ($true) {
                $key = [Console]::ReadKey($true).Key
                if ($key -eq [ConsoleKey]::Q) {
                    exit $lastExitCode
                }
                if ($key -eq [ConsoleKey]::R) {
                    Write-Host ""
                    $apexProcess = Start-ApexProcess
                    break
                }
            }
        }

        if ([Console]::KeyAvailable) {
            $key = [Console]::ReadKey($true).Key
            if ($key -eq [ConsoleKey]::Q) {
                Write-Host ""
                Write-Host "[APEX] Q pressed. Stopping app..."
                Stop-ProcessTree -ProcessId $apexProcess.Id
                exit 0
            }
            if ($key -eq [ConsoleKey]::R) {
                Write-Host ""
                Write-Host "[APEX] R pressed. Restarting app..."
                Stop-ProcessTree -ProcessId $apexProcess.Id
                Start-Sleep -Milliseconds 300
                $apexProcess = Start-ApexProcess
            }
        }

        Start-Sleep -Milliseconds 120
    }
} finally {
    if ($apexProcess -and -not $apexProcess.HasExited) {
        Stop-ProcessTree -ProcessId $apexProcess.Id
    }
}
