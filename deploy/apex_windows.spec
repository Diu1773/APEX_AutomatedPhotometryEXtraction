# -*- mode: python ; coding: utf-8 -*-

from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_data_files, collect_submodules, copy_metadata


ROOT = Path(SPECPATH).resolve().parent


def _dedupe_toc(*tocs):
    merged = []
    seen = set()
    for toc in tocs:
        for item in toc:
            key = str(item[0]).lower() if item else repr(item)
            if key in seen:
                continue
            seen.add(key)
            merged.append(item)
    return merged


def _collect_all(package: str):
    try:
        return collect_all(package)
    except Exception:
        return [], [], []


def _collect_data_files(package: str):
    try:
        return collect_data_files(package)
    except Exception:
        return []


def _copy_metadata(package: str):
    try:
        return copy_metadata(package)
    except Exception:
        return []


def _collect_submodules(package: str):
    try:
        return collect_submodules(package)
    except Exception:
        return []


def _filter_hidden(modules):
    filtered = []
    for module in modules:
        parts = module.split(".")
        if "tests" in parts or any(part.startswith("test_") for part in parts):
            continue
        filtered.append(module)
    return filtered


datas = []
binaries = []
hiddenimports = []

# Scientific packages with compiled extensions and package data. collect_all is
# intentionally broad because the GUI opens many workflow windows lazily.
for package in (
    "numpy",
    "bottleneck",
    "sep",
    "scipy",
    "pandas",
    "matplotlib",
    "astropy",
    "photutils",
    "pyvo",
    "pydantic",
):
    package_datas, package_binaries, package_hidden = _collect_all(package)
    datas += package_datas
    binaries += package_binaries
    hiddenimports += _filter_hidden(package_hidden)

# astroquery imports TAP/VOTable helpers dynamically. Collect the full package
# so installed builds can resolve Gaia/SIMBAD on clean machines.
package_datas, package_binaries, package_hidden = _collect_all("astroquery")
datas += package_datas
binaries += package_binaries
hiddenimports += _filter_hidden(package_hidden)
hiddenimports += [
    "astroquery.gaia",
    "astroquery.simbad",
    "astroquery.nasa_exoplanet_archive",
    "astroquery.utils.tap.core",
    "astroquery.utils.tap.conn.tapconn",
    "astroquery.utils.tap.model.job",
    "astroquery.utils.tap.xmlparser",
    "astroquery.vizier",
]

# Optional GUI tools. These are small compared with the core scientific stack
# and prevent post-install menu actions from failing due to hidden imports.
for package in ("batman", "emcee", "PIL", "requests", "certifi"):
    package_datas, package_binaries, package_hidden = _collect_all(package)
    datas += package_datas
    binaries += package_binaries
    hiddenimports += _filter_hidden(package_hidden)

hiddenimports += [
    "PyQt5.sip",
    "PyQt5.QtSvg",
    "PyQt5.QtPrintSupport",
    "matplotlib.backends.backend_qt5agg",
    "matplotlib.backends.backend_agg",
    "mpl_toolkits.mplot3d",
    "scipy.sparse.csgraph._validation",
    "scipy.ndimage",
    "scipy.ndimage._filters",
    "astropy.visualization",
    "astropy.visualization.interval",
    "astropy.visualization.stretch",
    "certifi",
]

# APEX uses mode-specific windows and tool windows through lazy imports. Collect
# the whole package so sub-widgets keep working in the installed app.
hiddenimports += _collect_submodules("apex")

datas += [(str(ROOT / "parameters.example.toml"), ".")]
datas += [(str(ROOT / "README.md"), ".")]
resources_dir = ROOT / "apex" / "resources"
icon_path = resources_dir / "apex.ico"
if resources_dir.is_dir():
    datas.append((str(resources_dir), "apex/resources"))

excludes = [
    "PyQt6",
    "PySide2",
    "PySide6",
    "pytest",
    "tests",
    "IPython",
    "ipykernel",
    "jupyter",
    "notebook",
    "sphinx",
    "tkinter",
]

a = Analysis(
    [str(ROOT / "main.py")],
    pathex=[str(ROOT)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={"matplotlib": {"backends": ["Qt5Agg", "Agg"]}},
    runtime_hooks=[],
    excludes=excludes,
    noarchive=False,
    optimize=1,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="APEX",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=str(icon_path) if icon_path.is_file() else None,
)

coll = COLLECT(
    exe,
    _dedupe_toc(a.binaries),
    _dedupe_toc(a.datas),
    strip=False,
    upx=False,
    upx_exclude=[],
    name="APEX",
)
