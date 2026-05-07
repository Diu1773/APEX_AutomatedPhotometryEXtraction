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
- astrometry.net (local installation for WCS solving)

## Installation

```bash
cd /path/to/Automated_Photometry_EXtraction
pip install -r requirements.txt
```

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
