"""Generate a self-contained synthetic reference FITS frame.

This module exists so the *entire* artificial-star benchmark pipeline
(empirical-PSF extraction + Step 4 detection + injection + matching) can run
with ZERO external data. The frame is a population of Gaussian-PSF stars over a
Poisson+read-noise background with a realistic FITS header (FILTER/EXPTIME/
EGAIN/RDNOISE/...), so:

* it contains plenty of well-detected, isolated, unsaturated stars for
  :func:`apex.benchmark.artificial_stars.extract_empirical_psf`, and
* its noise floor produces a *non-degenerate* completeness curve (bright stars
  recovered ~100%, faint stars ~0%) when artificial stars are injected across
  the magnitude range.

The module is intentionally Qt-free; it only depends on numpy + astropy.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.io import fits

__all__ = ["make_synthetic_reference_frame"]


def _gaussian_sigma_from_fwhm(fwhm_px: float) -> float:
    return float(fwhm_px) / 2.3548200450309493


def _stamp_gaussian(size: int, sigma: float) -> np.ndarray:
    """Unit-sum Gaussian PSF stamp of odd ``size``."""
    half = size // 2
    yy, xx = np.mgrid[-half : half + 1, -half : half + 1]
    stamp = np.exp(-(xx * xx + yy * yy) / (2.0 * sigma * sigma))
    total = float(stamp.sum())
    return stamp / total if total > 0 else stamp


def make_synthetic_reference_frame(
    out_path: str | Path,
    *,
    seed: int,
    n_stars: int = 120,
    size: int = 1024,
    fwhm_px: float = 3.5,
    zeropoint: float = 25.0,
    background: float = 800.0,
    gain: float = 1.5,
    read_noise: float = 8.0,
    mag_min: float = 14.0,
    mag_max: float = 20.5,
    filter_name: str = "r",
    exptime: float = 1.0,
    saturation_adu: float = 65000.0,
    edge_margin_px: int = 24,
    min_separation_px: float | None = None,
) -> Path:
    """Write a realistic synthetic FITS frame and return its path.

    Parameters
    ----------
    out_path:
        Destination ``.fits`` path (parent directories are created).
    seed:
        Deterministic RNG seed.
    n_stars:
        Number of "real" stars baked into the frame (the detector + empirical
        PSF use these). They span ``mag_min..mag_max`` so a good fraction are
        bright and isolated enough for PSF extraction.
    size:
        Square image side in pixels.
    fwhm_px:
        Gaussian PSF FWHM in pixels.
    zeropoint:
        Magnitude zero point: ``flux_e = 10 ** (-0.4 * (mag - zeropoint))``.
    background:
        Sky background level in ADU (Poisson-realized through ``gain``).
    gain:
        Detector gain in e-/ADU (Poisson statistics applied in electrons).
    read_noise:
        Read noise in electrons (Gaussian, added in electrons).
    mag_min, mag_max:
        Magnitude range of the baked-in stars (brightest..faintest).
    filter_name:
        FITS ``FILTER`` keyword. Use a key the detector recognizes (e.g. "r").
    exptime:
        FITS ``EXPTIME`` (the benchmark uses a unit-exposure convention).
    saturation_adu:
        ``SATURATE`` header value and the clip ceiling for the realized image.
    edge_margin_px:
        Keep baked-in stars at least this far from the edges.
    min_separation_px:
        Minimum separation between baked-in stars. Defaults to ``6 * fwhm_px``
        so a healthy isolated subset survives the empirical-PSF isolation cut.

    Returns
    -------
    Path
        The written FITS path.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if size <= 4 * edge_margin_px:
        raise ValueError("size is too small for the configured edge margin")
    if gain <= 0:
        raise ValueError("gain must be positive")

    rng = np.random.default_rng(int(seed))
    sigma = _gaussian_sigma_from_fwhm(fwhm_px)
    stamp_size = int(2 * int(np.ceil(4.0 * fwhm_px)) + 1)
    psf = _stamp_gaussian(stamp_size, sigma)
    half = stamp_size // 2

    if min_separation_px is None:
        min_separation_px = 6.0 * float(fwhm_px)

    # --- place non-overlapping star centers --------------------------------
    lo = edge_margin_px + half
    hi = size - edge_margin_px - half
    if hi <= lo:
        raise ValueError("image too small for the PSF stamp and margins")
    placed: list[tuple[float, float]] = []
    max_attempts = max(20000, n_stars * 400)
    attempts = 0
    while len(placed) < n_stars and attempts < max_attempts:
        attempts += 1
        x = float(rng.uniform(lo, hi))
        y = float(rng.uniform(lo, hi))
        if placed:
            arr = np.asarray(placed, dtype=float)
            if np.min(np.hypot(arr[:, 0] - x, arr[:, 1] - y)) < min_separation_px:
                continue
        placed.append((x, y))
    if len(placed) < max(3, n_stars // 2):
        raise RuntimeError(
            f"Could only place {len(placed)}/{n_stars} stars; reduce n_stars or "
            "min_separation_px, or increase size."
        )
    centers = np.asarray(placed, dtype=float)

    # Brightness distribution: roughly uniform in magnitude with a slight skew
    # toward the bright end so PSF extraction always has isolated bright stars.
    mags = rng.uniform(float(mag_min), float(mag_max), size=len(centers))
    mags = np.sort(mags)  # brightest first is irrelevant; order does not matter
    flux_e = 10.0 ** (-0.4 * (mags - float(zeropoint)))

    # --- build the noiseless electron image --------------------------------
    signal_e = np.zeros((size, size), dtype=np.float64)
    for (x, y), fe in zip(centers, flux_e):
        xi = int(round(x))
        yi = int(round(y))
        dx = x - xi
        dy = y - yi
        # Sub-pixel shift via a fresh sampled Gaussian (cheap, accurate enough).
        yy, xx = np.mgrid[-half : half + 1, -half : half + 1]
        shifted = np.exp(-((xx - dx) ** 2 + (yy - dy) ** 2) / (2.0 * sigma * sigma))
        total = float(shifted.sum())
        if total <= 0:
            continue
        shifted = shifted / total * float(fe)
        signal_e[yi - half : yi + half + 1, xi - half : xi + half + 1] += shifted

    # --- noise model in electrons ------------------------------------------
    background_e = float(background) * float(gain)
    expected_e = signal_e + background_e
    # Poisson photon noise on source + sky, then Gaussian read noise.
    realized_e = rng.poisson(np.clip(expected_e, 0.0, None)).astype(np.float64)
    realized_e += rng.normal(0.0, float(read_noise), size=realized_e.shape)

    # Convert back to ADU (the pipeline works in ADU + gain).
    image_adu = realized_e / float(gain)
    image_adu = np.clip(image_adu, 0.0, float(saturation_adu))
    image = image_adu.astype(np.float32)

    # --- header ------------------------------------------------------------
    header = fits.Header()
    header["BITPIX"] = -32
    header["FILTER"] = (str(filter_name), "Photometric filter (synthetic)")
    header["EXPTIME"] = (float(exptime), "Exposure time [s] (synthetic)")
    header["EGAIN"] = (float(gain), "Gain [e-/ADU] (synthetic)")
    header["GAIN"] = (float(gain), "Gain [e-/ADU] (synthetic)")
    header["RDNOISE"] = (float(read_noise), "Read noise [e-] (synthetic)")
    header["XBINNING"] = (1, "X binning (synthetic)")
    header["YBINNING"] = (1, "Y binning (synthetic)")
    header["SATURATE"] = (float(saturation_adu), "Saturation [ADU] (synthetic)")
    header["BUNIT"] = ("adu", "Pixel units")
    header["OBJECT"] = ("APEX_SYNTH_FIELD", "Synthetic validation field")
    header["APEXSYN"] = (True, "APEX synthetic reference frame")
    header["APEXSEED"] = (int(seed), "Synthetic frame RNG seed")
    header["APEXNST"] = (int(len(centers)), "Number of baked-in stars")
    header["APEXZP"] = (float(zeropoint), "Synthetic frame magnitude zeropoint")
    header["APEXBKG"] = (float(background), "Synthetic sky background [ADU]")
    header["APEXFWHM"] = (float(fwhm_px), "Synthetic PSF FWHM [px]")
    header["DATE"] = (datetime.now(timezone.utc).isoformat(), "Creation time (UTC)")
    header.add_history("Synthetic frame generated by apex.benchmark.synthetic_frame")

    fits.writeto(out_path, image, header, overwrite=True)
    return out_path
