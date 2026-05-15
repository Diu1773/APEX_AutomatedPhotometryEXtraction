# APEX — Automated Photometry EXtraction

APEX is a PyQt5-based GUI toolkit for aperture and PSF photometry of astronomical images. It supports two operational modes:

- **CMD mode** (`apex/cmd/`): Cluster photometry pipeline — detection through CMD diagram and isochrone fitting (12 steps).
- **LC mode** (`apex/lightcurve/`): Light curve analysis pipeline — multi-night photometry, detrending, and period analysis (11 steps).

Both modes share a common core through Step 7: file selection, crop, sky
preview, source detection, WCS plate solving, master catalog build, and forced
aperture photometry. CMD and LC then branch into mode-specific steps.

## Requirements

- Python 3.10+
- PyQt5
- astropy >= 5.0
- photutils >= 1.5
- numpy
- pandas
- scipy
- matplotlib
- tomli / tomllib (Python 3.11+ built-in)
- astroquery (Gaia access)
- certifi (HTTPS CA bundle for packaged Gaia/SIMBAD queries)
- ASTAP with an external star database, or local astrometry.net with index
  files, for Step 5 WCS solving

## Installation

```bash
cd /path/to/Automated_Photometry_EXtraction
pip install -r requirements.txt
```

## WCS Solver Setup

APEX does not bundle external plate solvers or their star/index databases.
Install at least one local solver before running Step 5.

### ASTAP

ASTAP is the recommended first solver. Install ASTAP for Windows and one ASTAP
star database, then set `ASTAP CLI Path` and `ASTAP Star DB` in Step 5 >
ASTAP Parameters.

- ASTAP: <https://www.hnsky.org/astap.htm>
- ASTAP star databases: <https://sourceforge.net/projects/astap-program/files/star_databases/>

APEX currently passes `D80` or `D50` to ASTAP. Use `D80` as the default; `D50`
is smaller and is normally useful when the image field of view is comfortably
above about 0.2 deg. The installed database and the Step 5 `ASTAP Star DB`
setting must match.

### Local astrometry.net

Local astrometry.net is an optional fallback when ASTAP fails or does not leave
a valid WCS header. On Windows, use WSL/Ubuntu and install `solve-field` plus
index files matching the image field of view.

- Astrometry.net README: <https://astrometry.net/doc/readme.html>
- Astrometry.net index files: <https://data.astrometry.net/>

The astrometry.net docs recommend downloading only the index scales needed for
your images. In practice, missing or mismatched index files are the most common
reason for `solve-field` no-solution failures.

Before Step 5 starts solving, APEX checks whether the configured ASTAP
executable is reachable. If ASTAP is unavailable, or if ASTAP fails and local
astrometry.net fallback is enabled, APEX then checks the local/WSL
`solve-field` command. Solver databases and astrometry.net index coverage are
still external data requirements, so an installed executable can still fail if
the selected D80/D50 database or index file scale set is missing. APEX also
checks whether `astroquery.gaia` and the `certifi` CA bundle are importable;
without them, WCS solving can still succeed but Gaia attach, refinement,
residual medians, and Step 6 Gaia stats will be unavailable unless a compatible
`gaia_fov.ecsv` cache already exists.

Step 5 uses `[gaia].wcs_mag_max` as a server-side Gaia/VizieR query cap for WCS
refinement and QC. Keep this near `18`-`20` for large fields; the broader
`[gaia].mag_max` value can remain higher for later catalog products.

## Quickstart

```bash
# CMD mode (cluster photometry)
python apex/cmd/main.py

# LC mode (light curve analysis)
python apex/lightcurve/main.py
```

Or from the project root launcher:
```bash
python main.py
```

## Screenshot

<!-- TODO: add screenshot -->

## License

<!-- TODO: add license -->
