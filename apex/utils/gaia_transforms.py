"""Gaia photometric transformation tables shared across the pipeline.

All transforms use the form:
    G - band_mag = poly(G_BP - G_RP)
where poly coefficients are in ascending power order (constant, linear, …).

Sources
-------
Riello+2021  A&A 649 A3  — Gaia EDR3, Table 5.7  (V, R_c, I_c)
Jordi+2010   A&A 523 A48 — Gaia pre-launch SDSS   (g, r, i, z)
B            Derived: Riello+2021 V + empirical B-V(BP-RP)     σ ≈ 0.08 mag
U            Approximate only, FGK range (BP-RP < 1.5)         σ ≈ 0.20 mag
"""
from __future__ import annotations

# (coefficients, bpRP_lo, bpRP_hi, source, sigma_mag)
GAIA_TO_BAND: dict[str, tuple] = {
    "V": ([-0.02704,  0.01424, -0.2156,   0.01426], -0.5, 5.0, "Riello+2021",  0.030),
    "R": ([-0.02275,  0.39610, -0.1243,  -0.01396],  0.0, 4.0, "Riello+2021",  0.032),
    "I": ([ 0.01753,  0.76000, -0.0991,   0.03765], -0.5, 4.5, "Riello+2021",  0.045),
    "B": ([-0.0295,  -0.7543,  -0.0702,   0.0088],  -0.1, 3.0, "derived",      0.080),
    "U": ([-0.020,   -0.980,   -0.320,    0.050],    0.0, 1.5, "approx",       0.200),
    "g": ([ 0.2199,  -0.6365,  -0.1548,   0.0064],   0.3, 3.0, "Jordi+2010",   0.050),
    "r": ([-0.09837,  0.08592,  0.1907,  -0.1701, 0.02263], 0.0, 3.0, "Jordi+2010", 0.050),
    "i": ([-0.293,    0.6404,  -0.09609, -0.002104], 0.5, 2.0, "Jordi+2010",   0.050),
    "z": ([-0.09718,  0.08116,  0.2460,  -0.2294,  0.04764], 0.0, 3.5, "Jordi+2010", 0.060),
}

# Preferred instrumental color index per filter (first available pair wins).
FILTER_COLOR_PREF: dict[str, list[tuple[str, str]]] = {
    "U": [("U", "B"), ("U", "V")],
    "B": [("B", "V"), ("B", "R")],
    "V": [("B", "V"), ("V", "R"), ("V", "I")],
    "R": [("V", "R"), ("B", "R"), ("R", "I")],
    "I": [("R", "I"), ("V", "I")],
    "g": [("g", "r"), ("g", "i")],
    "r": [("g", "r"), ("r", "i")],
    "i": [("r", "i"), ("g", "i")],
    "z": [("r", "z"), ("i", "z")],
}

# Header value → canonical key in GAIA_TO_BAND
BAND_ALIASES: dict[str, str] = {
    "Rc": "R", "RC": "R", "R_c": "R", "Rc_j": "R",
    "Ic": "I", "IC": "I", "I_c": "I", "Ic_j": "I",
    "Bj": "B", "BJ": "B", "B_j": "B",
    "Vj": "V", "VJ": "V", "V_j": "V",
    "Uj": "U", "UJ": "U", "U_j": "U",
}
