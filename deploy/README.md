# APEX Deployment

APEX is distributed as a Windows installer built from a PyInstaller `onedir`
bundle. The build bundles the Python runtime, PyQt5, scientific Python packages,
`apex/resources`, and `parameters.example.toml`.

On first launch, `APEX.exe` copies `parameters.example.toml` to
`parameters.toml` next to the executable when no local runtime config exists.

External WCS solvers are not bundled. Release users must install at least one
solver and its data before Step 5 WCS solving:

- ASTAP executable and an ASTAP star database must be installed separately.
  APEX currently passes `D80` or `D50` to ASTAP, so the installed database and
  the Step 5 `ASTAP Star DB` setting must match. Use `D80` as the default;
  `D50` is smaller and is normally useful when the image field of view is
  comfortably above about 0.2 deg.
- Local astrometry.net / `solve-field` is an optional fallback and must be
  installed separately, typically through WSL on Windows. Index files matching
  the image field of view are required; missing or mismatched indexes usually
  cause no-solution failures.
- Users should point `parameters.toml` or the Step 5 parameter dialogs at those
  local solver paths after first launch.
- APEX checks ASTAP before Step 5 starts, and checks local/WSL `solve-field`
  when ASTAP is unavailable or fallback is actually needed. Missing D80/D50
  database or astrometry.net index coverage can still appear as solver
  no-solution failures.
- The release smoke test also imports `astroquery.gaia`, `astroquery.simbad`,
  TAP support, and `certifi`. If those imports fail, Gaia attach, SIMBAD
  resolution, WCS refinement, residual medians, and Step 6 Gaia stats can be
  absent in installed builds.
- `[gaia].wcs_mag_max` limits Step 5 Gaia/VizieR WCS queries server-side. The
  release default is `18.0` to avoid large-field TAP timeouts.

Official setup references:

- ASTAP: <https://www.hnsky.org/astap.htm>
- ASTAP star databases: <https://sourceforge.net/projects/astap-program/files/star_databases/>
- Astrometry.net README: <https://astrometry.net/doc/readme.html>
- Astrometry.net index files: <https://data.astrometry.net/>

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

1. Run `python -m compileall apex main.py scripts deploy`.
2. Run `python -m pytest tests`.
3. Build with `build.bat`.
4. Confirm the source preflight and `APEX.exe --smoke` steps pass.
5. Confirm `release\Setup\setup.exe` exists.
6. Install `setup.exe` on a clean Windows user account or VM.
7. Confirm both CMD and LC modes open.
8. Confirm `parameters.toml` is created on first launch.
9. Confirm ASTAP / astrometry.net paths are documented for the release.
