"""Gaia photometric transformation tables shared across the pipeline.

All transforms use ``G - band_mag = poly(G_BP - G_RP)`` with coefficients
stored in ascending power order.

The default table combines Pancino et al. (2022) for Johnson B dwarfs,
Riello et al. (2021) for Johnson-Cousins V/R/I, and Jordi et al. (2010)
for Sloan g/r/i/z. Named source tables can be selected with
``get_gaia_to_band(source=...)``.

Each entry is ``(coefficients, BP-RP_lo, BP-RP_hi, source_label, sigma_mag)``.
"""
from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# Per-source transform tables
# ---------------------------------------------------------------------------

_RIELLO2021: dict[str, tuple] = {
    # Riello+2021 (A&A 649 A3) Table C.2 — Gaia EDR3 → Johnson-Cousins (Vega)
    "V": ([-0.02704,  0.01424, -0.2156,   0.01426],            -0.5, 5.0, "Riello+2021", 0.030),
    "R": ([-0.02275,  0.39610, -0.1243,  -0.01396, 0.003775],   0.0, 4.0, "Riello+2021", 0.040),
    "I": ([ 0.01753,  0.76000, -0.0991],                        0.0, 4.5, "Riello+2021", 0.045),
    "B": ([-0.0295,  -0.7543,  -0.0702,   0.0088],             -0.1, 3.0, "derived",     0.080),
    "U": ([-0.020,   -0.980,   -0.320,    0.050],               0.0, 1.5, "approx",      0.200),
}

# Pancino et al. (2022), A&A 664, A109, Table 9, row 29.
# The published relation is B-G = poly(BP-RP) for dwarfs. APEX stores
# G-band = poly(BP-RP), so the coefficients are sign-reversed here.
_PANCINO2022: dict[str, tuple] = {
    "B": (
        [
             0.00902616221721914,
            -0.703856376726825,
            -0.439617600260539,
             0.201347189673066,
             0.545325598993276,
            -1.3392842024655,
             1.1141015574016,
            -0.435726384509969,
             0.0818759549249881,
            -0.00597884341541752,
        ],
        -0.4,
        3.5,
        "Pancino+2022 dwarf",
        0.0248,
    ),
}

_JORDI2010: dict[str, tuple] = {
    "g": ([ 0.2199,  -0.6365,  -0.1548,   0.0064],           0.3, 3.0, "Jordi+2010", 0.050),
    "r": ([-0.09837,  0.08592,  0.1907,  -0.1701, 0.02263],  0.0, 3.0, "Jordi+2010", 0.050),
    "i": ([-0.293,    0.6404,  -0.09609, -0.002104],          0.5, 2.0, "Jordi+2010", 0.050),
    "z": ([-0.09718,  0.08116,  0.2460,  -0.2294,  0.04764], 0.0, 3.5, "Jordi+2010", 0.060),
}

# Named source registry — extend here when new calibrations are published.
GAIA_TRANSFORM_TABLES: dict[str, dict[str, tuple]] = {
    "pancino2022": _PANCINO2022,
    "riello2021": _RIELLO2021,
    "jordi2010":  _JORDI2010,
    # "evans2018": _EVANS2018,   # placeholder for Gaia DR2 SDSS (Evans+2018 A&A 616 A4)
}

# Default merged table (backward-compatible name kept for existing imports).
GAIA_TO_BAND: dict[str, tuple] = {
    **_RIELLO2021,
    **_PANCINO2022,
    **_JORDI2010,
}


def get_gaia_to_band(source: str | None = None) -> dict[str, tuple]:
    """Return the Gaia→band transform table for the given *source*.

    Parameters
    ----------
    source : str or None
        One of the keys in ``GAIA_TRANSFORM_TABLES`` (e.g. ``"riello2021"``,
        ``"jordi2010"``), or ``"default"`` / ``None`` for the merged table.

    Returns
    -------
    dict mapping canonical filter name → (coeffs, bp_rp_lo, bp_rp_hi, label, sigma)
    """
    if not source or str(source).strip().lower() in ("default", "auto", ""):
        return GAIA_TO_BAND
    key = str(source).strip().lower()
    if key not in GAIA_TRANSFORM_TABLES:
        import warnings
        warnings.warn(
            f"Unknown gaia_transform_source={source!r}; "
            f"available: {sorted(GAIA_TRANSFORM_TABLES)}. Using default.",
            UserWarning, stacklevel=2,
        )
        return GAIA_TO_BAND
    return GAIA_TRANSFORM_TABLES[key]

# Preferred instrumental color index per filter (first available pair wins).
# Cross-system fallbacks are listed last so same-system pairs are always preferred.
FILTER_COLOR_PREF: dict[str, list[tuple[str, str]]] = {
    "U": [("U", "B"), ("U", "V")],
    "B": [("B", "V"), ("B", "R"), ("B", "r")],
    "V": [("B", "V"), ("V", "R"), ("V", "r"), ("V", "I")],
    "R": [("V", "R"), ("B", "R"), ("R", "I"), ("r", "i"), ("B", "V")],
    "I": [("R", "I"), ("V", "I"), ("r", "i"), ("B", "V")],
    "g": [("g", "r"), ("g", "i"), ("B", "V")],
    "r": [("g", "r"), ("r", "i"), ("B", "V"), ("V", "r"), ("B", "r")],
    "i": [("r", "i"), ("g", "i"), ("B", "V"), ("r", "z")],
    "z": [("r", "z"), ("i", "z"), ("r", "i"), ("B", "V")],
}

# Header value → canonical key in GAIA_TO_BAND
BAND_ALIASES: dict[str, str] = {
    "Rc": "R", "RC": "R", "R_c": "R", "Rc_j": "R",
    "Ic": "I", "IC": "I", "I_c": "I", "Ic_j": "I",
    "Bj": "B", "BJ": "B", "B_j": "B",
    "Vj": "V", "VJ": "V", "V_j": "V",
    "Uj": "U", "UJ": "U", "U_j": "U",
}

# ---------------------------------------------------------------------------
# Filter wavelength (nm) — used for sorting bands and building color pairs
# ---------------------------------------------------------------------------
FILTER_WAVELENGTH_NM: dict[str, float] = {
    "U":   365.0,
    "B":   440.0,
    "g":   480.0,
    "Hb":  486.0,
    "V":   550.0,
    "r":   625.0,
    "R":   640.0,
    "Ha":  656.0,
    "i":   760.0,
    "I":   800.0,
    "z":   900.0,
    "Y":   1020.0,
}


def filter_wavelength(filt: str) -> float:
    """Return approximate central wavelength (nm) for a filter name."""
    return FILTER_WAVELENGTH_NM.get(filt, 600.0)


def filter_bands_from_columns(columns, prefix: str) -> list[str]:
    """Return canonical filter names parsed from columns like ``mag_inst_B``.

    Skips derived columns (``mag_inst_err_B``, ``mag_inst_wmean_B`` …) by requiring
    the suffix to contain no underscore.  Result is sorted by central wavelength.
    """
    bare = {
        s for c in columns if c.startswith(prefix)
        for s in [c[len(prefix):]] if "_" not in s and s
    }
    return sorted(bare, key=filter_wavelength)


def build_color_pairs(bands: list[str], adjacent_only: bool = True) -> list[tuple[str, str]]:
    """Return (bluer, redder) color index pairs from *bands*.

    With ``adjacent_only=True`` (default) only consecutive pairs in wavelength
    order are returned — e.g. for BVR: [(B,V), (V,R)].  This matches the
    conventional CMD color indices and keeps the UI uncluttered.

    With ``adjacent_only=False`` all combinations are returned, which is useful
    for ZP calibration fallback logic.
    """
    avail = sorted(bands, key=filter_wavelength)
    if adjacent_only:
        return [(avail[i], avail[i + 1]) for i in range(len(avail) - 1)]
    pairs = [
        (avail[i], avail[j])
        for i in range(len(avail))
        for j in range(i + 1, len(avail))
    ]
    pairs.sort(key=lambda p: abs(filter_wavelength(p[1]) - filter_wavelength(p[0])))
    return pairs


# ---------------------------------------------------------------------------
# Teff vs color index anchors for CMD star-temperature colormap
# Sources: Pecaut & Mamajek 2013 (Johnson), Covey+2007 (SDSS), Gaia DR3
# ---------------------------------------------------------------------------
TEFF_COLOR_ANCHORS: dict[str, dict] = {
    # Johnson-Cousins
    "B-V": {
        "x": np.array([-0.30, -0.16,  0.00,  0.28,  0.58,  0.82,  1.14,  1.60], float),
        "t": np.array([35000, 20000,  9500,  7200,  5900,  5300,  4000,  3200], float),
    },
    "V-R": {
        "x": np.array([-0.15, -0.05,  0.00,  0.22,  0.35,  0.50,  0.75,  1.20], float),
        "t": np.array([35000, 20000, 10000,  7500,  6000,  5000,  4000,  3200], float),
    },
    "V-I": {
        "x": np.array([-0.30, -0.10,  0.00,  0.45,  0.70,  1.05,  1.55,  2.30], float),
        "t": np.array([35000, 20000, 10000,  7500,  6000,  5000,  4000,  3200], float),
    },
    "R-I": {
        "x": np.array([-0.15, -0.05,  0.00,  0.20,  0.33,  0.50,  0.80,  1.10], float),
        "t": np.array([35000, 20000, 10000,  7500,  6000,  5000,  4000,  3200], float),
    },
    "B-R": {
        "x": np.array([-0.45, -0.20,  0.00,  0.55,  0.90,  1.30,  1.90,  2.80], float),
        "t": np.array([35000, 20000, 10000,  7500,  6000,  5000,  4000,  3200], float),
    },
    "U-B": {
        "x": np.array([-1.20, -0.90, -0.45,  0.00,  0.10,  0.60,  1.10,  1.40], float),
        "t": np.array([35000, 20000, 10000,  8000,  7000,  5500,  4500,  3500], float),
    },
    # SDSS
    "g-r": {
        "x": np.array([-0.40, -0.20,  0.00,  0.30,  0.45,  0.80,  1.40,  1.80], float),
        "t": np.array([35000, 20000, 10000,  7500,  6000,  4500,  3200,  2400], float),
    },
    "r-i": {
        "x": np.array([-0.30, -0.20, -0.05,  0.10,  0.20,  0.40,  0.80,  1.10], float),
        "t": np.array([35000, 20000, 10000,  7500,  6000,  4500,  3200,  2400], float),
    },
    "g-i": {
        "x": np.array([-0.70, -0.45,  0.00,  0.50,  0.80,  1.50,  2.60,  3.20], float),
        "t": np.array([35000, 20000, 10000,  7500,  6000,  4500,  3200,  2400], float),
    },
    "i-z": {
        "x": np.array([-0.20, -0.10,  0.00,  0.10,  0.20,  0.40,  0.60,  0.80], float),
        "t": np.array([35000, 20000, 10000,  7500,  6000,  4500,  3200,  2400], float),
    },
    # Gaia
    "BP-RP": {
        "x": np.array([-0.40, -0.20,  0.00,  0.40,  0.80,  1.30,  2.30,  3.50], float),
        "t": np.array([35000, 20000, 10000,  7500,  6000,  4500,  3200,  2400], float),
    },
}


def teff_from_color(color_values: "np.ndarray", band_a: str, band_b: str,
                    teff_min: float = 2400.0, teff_max: float = 40000.0) -> "np.ndarray":
    """Interpolate effective temperature from a color index array.

    Uses ``TEFF_COLOR_ANCHORS`` for the given ``band_a - band_b`` color.
    Returns NaN for unknown color combinations.
    """
    key = f"{band_a}-{band_b}"
    anchor = TEFF_COLOR_ANCHORS.get(key)
    if anchor is None or len(color_values) == 0:
        return np.full_like(color_values, np.nan, dtype=float)
    x = np.asarray(color_values, float)
    teff = np.interp(x, anchor["x"], anchor["t"])
    return np.clip(teff, teff_min, teff_max)
