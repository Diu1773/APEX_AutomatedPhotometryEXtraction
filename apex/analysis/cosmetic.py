"""Cosmetic correction — cosmic rays + hot/cold pixels (Step 0, optional).

Qt-free.  Cosmic-ray rejection is delegated to :mod:`astroscrappy` (the
AstroPy-affiliated C implementation of the L.A.Cosmic algorithm,
van Dokkum 2001, PASP 113, 1420), which is the standard single-frame CR
rejector in CCD photometry pipelines and natively protects saturated stars
from being flagged as cosmic rays.  Persistent hot/cold pixels are handled by a
deterministic bad-pixel mask derived from the master dark/flat and interpolated
over.

This is scientifically standard reduction (astropy CCD reduction guide §6.3-6.4;
IRAF ``cosmicray``/``ccdmask``+``fixpix``), NOT imaging-only cosmetics.  It is
optional and off by default; when on, ``calibrate_light`` calls
:func:`clean_frame` after flat-fielding.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from apex.utils import fast_stats
from apex.utils.constants import MAD_TO_SIGMA

try:
    from astroscrappy import detect_cosmics
    HAS_ASTROSCRAPPY = True
except Exception:                       # pragma: no cover - optional dependency
    HAS_ASTROSCRAPPY = False


def hot_pixel_mask(master_dark: Optional[np.ndarray], sigma: float = 6.0) -> Optional[np.ndarray]:
    """Boolean mask of persistent hot pixels from the master dark (dark >> median)."""
    if master_dark is None:
        return None
    arr = np.asarray(master_dark)
    med, mad = fast_stats.robust_median_mad(arr, 0.0)
    sig = mad * MAD_TO_SIGMA
    if not np.isfinite(sig) or sig <= 0:      # uniform master → fall back to std
        sig = float(np.nanstd(arr))
    if not np.isfinite(sig) or sig <= 0:
        return None
    return arr > (med + sigma * sig)


def star_protect_mask(data: np.ndarray, sigma: float = 4.0) -> Optional[np.ndarray]:
    """Mask of extended bright cores (real stars) to shield from CR flagging.

    L.A.Cosmic can mistake the sharp peak of an under-sampled star for a cosmic
    ray.  A star's core is *extended* (its 3x3 neighbourhood is also bright),
    whereas a thin cosmic ray or hot pixel is not — so requiring a bright local
    mean separates stars from artefacts.  Passed to astroscrappy's ``inmask``.
    """
    try:
        from scipy.ndimage import convolve, binary_dilation
    except Exception:
        return None
    finite = np.nan_to_num(np.asarray(data, dtype=np.float64),
                           nan=float(np.nanmedian(data)))
    med, mad = fast_stats.robust_median_mad(finite, 0.0)
    sig = mad * MAD_TO_SIGMA
    if not np.isfinite(sig) or sig <= 0:      # noiseless synthetic → std fallback
        sig = float(np.nanstd(finite))
    if not np.isfinite(sig) or sig <= 0:
        return None
    bright = finite > med + sigma * sig
    # A star core is surrounded by bright pixels; a hot pixel or thin cosmic-ray
    # streak has few/no bright neighbours. Require >=4 bright neighbours (of 8).
    kernel = np.ones((3, 3), dtype=int)
    kernel[1, 1] = 0
    n_bright_neigh = convolve(bright.astype(int), kernel, mode="constant")
    star_core = bright & (n_bright_neigh >= 4)
    return binary_dilation(star_core, iterations=1)


def _interp_mask(data: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Replace masked pixels with the local 3x3 median (isolated bad pixels)."""
    if mask is None or not np.any(mask):
        return data
    try:
        from scipy.ndimage import median_filter
    except Exception:
        return data
    out = data.copy()
    med = median_filter(np.nan_to_num(data, nan=float(np.nanmedian(data))),
                        size=3, mode="nearest")
    out[mask] = med[mask]
    return out


def clean_frame(data: np.ndarray, *, gain: float = 1.0, readnoise: float = 6.5,
                satlevel: float = 65535.0, sigclip: float = 4.5, objlim: float = 5.0,
                hot_mask: Optional[np.ndarray] = None, protect_stars: bool = True
                ) -> Tuple[np.ndarray, np.ndarray, int]:
    """Remove cosmic rays (L.A.Cosmic / astroscrappy) and interpolate hot pixels.

    ``satlevel``/``objlim`` and an explicit star-protection mask keep saturated
    and real point sources from being clipped.  Dead pixels (NaN, from
    flat-fielding) are preserved.  Returns ``(cleaned, mask, n_pixels)`` where
    ``mask`` marks corrected pixels.
    """
    if not HAS_ASTROSCRAPPY:
        raise RuntimeError("astroscrappy is required for cosmetic correction "
                           "(pip install astroscrappy)")
    arr = np.asarray(data, dtype=np.float32)
    dead = ~np.isfinite(arr)
    fill = float(np.nanmedian(arr)) if np.any(np.isfinite(arr)) else 0.0
    clean_in = np.where(dead, fill, arr)

    inmask = dead.copy()
    if hot_mask is not None:
        inmask = inmask | np.asarray(hot_mask, dtype=bool)
    if protect_stars:
        stars = star_protect_mask(clean_in)
        if stars is not None:
            inmask = inmask | stars

    crmask, cleaned = detect_cosmics(
        np.ascontiguousarray(clean_in), inmask=inmask,
        gain=float(gain), readnoise=float(readnoise), satlevel=float(satlevel),
        sigclip=float(sigclip), objlim=float(objlim), cleantype="meanmask",
    )
    cleaned = np.asarray(cleaned, dtype=np.float64)

    # astroscrappy ignores inmask pixels; interpolate the hot ones ourselves.
    if hot_mask is not None:
        cleaned = _interp_mask(cleaned, np.asarray(hot_mask, dtype=bool))

    corrected = np.asarray(crmask, dtype=bool)
    if hot_mask is not None:
        corrected = corrected | np.asarray(hot_mask, dtype=bool)
    corrected &= ~dead                  # dead pixels were not "corrected"
    cleaned[dead] = np.nan              # preserve dead pixels as NaN
    return cleaned, corrected, int(np.count_nonzero(corrected))
