# APEX Deployment

APEX is distributed as a Windows installer built from a PyInstaller `onedir`
bundle. The build bundles the Python runtime, PyQt5, scientific Python packages,
`apex/resources`, and `parameters.example.toml`.

On first launch, `APEX.exe` copies `parameters.example.toml` to
`parameters.toml` next to the executable when no local runtime config exists.

External WCS solvers are not bundled:

- ASTAP executable and star database must be installed separately.
- Local astrometry.net / `solve-field` must be installed separately, typically
  through WSL on Windows.
- Users should point `parameters.toml` at those local solver paths after first
  launch.

## Local Windows Build

From the repository root:

```powershell
.\build.bat
```

The distributables are written to:

```text
release\Setup\
```

## Release Checklist

1. Run `python -m compileall apex main.py scripts`.
2. Run `python -m pytest tests`.
3. Build with `build.bat`.
4. Confirm `release\Setup\setup.exe` exists.
5. Install `setup.exe` on a clean Windows user account or VM.
6. Confirm both CMD and LC modes open.
7. Confirm `parameters.toml` is created on first launch.
8. Confirm ASTAP / astrometry.net paths are documented for the release.
