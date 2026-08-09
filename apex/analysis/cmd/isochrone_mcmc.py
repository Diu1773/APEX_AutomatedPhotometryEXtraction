"""Generative-mixture isochrone solver sampled by emcee (Qt-free).

Why this module exists
----------------------
The legacy :mod:`isochrone_fitter_v2` objective is a nearest-neighbour chi2 over
the closest ``fit_fraction`` of stars.  That objective does NOT encode CMD
*density / morphology*: a young isochrone whose main-sequence turn-off sits a bit
bluer can score nearly identical chi2 to an old one because every observed star
is matched to *some* nearest track point regardless of how many model stars
actually live there.  Age is therefore unconstrained and the optimiser rails on
bounds.

This module replaces that with a **generative forward model**.  Each observed
star is assumed to be drawn from a mixture of three populations:

    L_i = (1 - f_field) * [ (1 - f_bin) * L_single_i + f_bin * L_binary_i ]
          + f_field * L_field

where the single-star term is the IMF-weighted Gaussian convolution of the
isochrone track (so a populous main sequence genuinely costs less than a sparse
giant branch), the binary term is the same track shifted along the
equal-/intermediate-mass binary ridge, and the field term is a flat density over
the CMD bounding box that absorbs non-members and outliers (the principled
replacement for the "closest fraction" hack).  Because the likelihood integrates
the IMF density along the track, the turn-off luminosity -- hence age -- is
constrained.

The posterior is sampled with :mod:`emcee` so we get full credible intervals and
a corner-plottable chain instead of a single bound-railed point estimate.

Reuse
-----
Isochrone interpolation, distance/extinction application, IMF weighting, the
grid-scan initialiser and :class:`FitBounds` all come from
:mod:`apex.analysis.cmd.isochrone_fitter_v2`.  Nothing here re-implements grid
interpolation.

Imports are restricted to numpy / scipy / emcee / apex.analysis / apex.utils --
no PyQt5 -- so this is safe to import in headless validation and tests.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from apex.analysis.cmd.isochrone_fitter_v2 import FitBounds

# Type aliases for the two callables we borrow from IsochroneFitterV2.
# interp_fn(log_age, mh, dm, e_color) -> (iso_color, iso_mag, iso_mass)
InterpFn = Callable[[float, float, float, float], Tuple[np.ndarray, np.ndarray, np.ndarray]]
# imf_fn(mass) -> per-point Kroupa density (need not be normalised; we renormalise)
ImfFn = Callable[[np.ndarray], np.ndarray]
# Multi-band interpolator: interp_bands_fn(log_age, mh) -> (bands_dict, mass)
# where bands_dict maps band-name -> ABSOLUTE magnitude array (no dm / extinction
# applied); the caller applies the distance modulus + per-band extinction itself.
InterpBandsFn = Callable[[float, float], Tuple[Dict[str, np.ndarray], np.ndarray]]

# Equal-mass-binary brightening: two equal stars co-add flux ->
# Δm = -2.5*log10(2) = -0.7526 mag, at (essentially) the same colour.
_EQUAL_MASS_DMAG = -2.5 * np.log10(2.0)

# Memory cap for each (N_obs_block, N_iso) Gaussian block.  Kept small (a few
# MB) on purpose: the codebase's other batched numerics cap at ~50 MB, but here
# we additionally want robustness on memory-constrained / fragmented hosts where
# a single large contiguous reservation can fail even with RAM nominally free,
# so blocks are sized to a few MB and buffers are reused in place.
_MAX_BLOCK_BYTES = 512 * 1024  # 512 KB per float32 (B, K) block

# Numerical floors.
_TINY = 1e-300            # additive floor before log() to avoid log(0)
_ERR_FLOOR = 0.02         # systematic/model-error floor (mag), added in quadrature.
# Represents the irreducible CMD width from isochrone-model error, unresolved
# binaries, differential reddening and calibration -- NOT just photon noise. With
# the old 0.005 mag value the combined many-star likelihood was razor-sharp
# (acceptance collapsed to a few percent and walkers stuck at the seed); ~0.02 mag
# is the realistic ground-based floor that makes the posterior samplable and the
# credible intervals honest. General (same floor for all data), not per-target.
_MIN_ISO_POINTS = 10      # below this the track is unusable -> -inf


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------
@dataclass
class McmcResult:
    """Posterior summary from :func:`fit_isochrone_mcmc`.

    Each ``*_p`` field is a 3-tuple ``(p16, p50, p84)``.  Convenience scalars
    (``*_med``, ``*_err_lo``, ``*_err_hi``) expose the median and asymmetric
    1-sigma credible interval.
    """

    # Sampled-parameter percentiles (16/50/84).
    log_age_p: Tuple[float, float, float]
    mh_p: Tuple[float, float, float]
    dm_p: Tuple[float, float, float]
    e_color_p: Tuple[float, float, float]

    # Derived-quantity percentiles.
    age_gyr_p: Tuple[float, float, float]
    distance_pc_p: Tuple[float, float, float]
    e_bv_p: Optional[Tuple[float, float, float]] = None  # requires R coeffs

    # Optional nuisance parameter (binary fraction), if sampled.
    f_bin_p: Optional[Tuple[float, float, float]] = None

    # Chain + diagnostics.
    flat_chain: Optional[np.ndarray] = None       # (n_samples, n_dim) downsampled
    labels: Sequence[str] = field(default_factory=list)
    acceptance_fraction: float = 0.0
    autocorr_time: Optional[np.ndarray] = None
    n_burn: int = 0
    n_steps: int = 0
    n_walkers: int = 0
    n_dim: int = 4
    convergence_ok: bool = False
    # Why convergence was or was not established. Empty means "converged".
    convergence_detail: str = ""
    # Split-R-hat per parameter; near 1.0 means the ensemble halves agree.
    rhat: Optional[np.ndarray] = None
    sample_binary_fraction: bool = False

    # Median parameter vector in the *sampled* space (for overlay plotting).
    median_theta: Optional[np.ndarray] = None

    # -- convenience accessors -------------------------------------------------
    @staticmethod
    def _med(p: Optional[Tuple[float, float, float]]) -> Optional[float]:
        return None if p is None else float(p[1])

    @property
    def log_age_med(self) -> float:
        return self._med(self.log_age_p)

    @property
    def mh_med(self) -> float:
        return self._med(self.mh_p)

    @property
    def dm_med(self) -> float:
        return self._med(self.dm_p)

    @property
    def e_color_med(self) -> float:
        return self._med(self.e_color_p)

    @property
    def age_gyr_med(self) -> float:
        return self._med(self.age_gyr_p)

    @property
    def distance_pc_med(self) -> float:
        return self._med(self.distance_pc_p)

    def contains(self, name: str, true_value: float) -> bool:
        """Return True if ``true_value`` lies within the 16-84 percentile band."""
        p = getattr(self, f"{name}_p", None)
        if p is None:
            return False
        return bool(p[0] <= true_value <= p[2])

    def summary(self) -> str:
        def fmt(p):
            if p is None:
                return "n/a"
            lo, mid, hi = p
            return f"{mid:.4g} (+{hi - mid:.3g} / -{mid - lo:.3g})"

        lines = [
            "=" * 56,
            "Isochrone MCMC posterior",
            "=" * 56,
            f"  log(Age)    = {fmt(self.log_age_p)}",
            f"  Age [Gyr]   = {fmt(self.age_gyr_p)}",
            f"  [M/H]       = {fmt(self.mh_p)}",
            f"  (m-M)       = {fmt(self.dm_p)}",
            f"  distance pc = {fmt(self.distance_pc_p)}",
            f"  E(color)    = {fmt(self.e_color_p)}",
        ]
        if self.e_bv_p is not None:
            lines.append(f"  E(B-V)      = {fmt(self.e_bv_p)}")
        if self.f_bin_p is not None:
            lines.append(f"  f_bin       = {fmt(self.f_bin_p)}")
        lines.extend([
            "-" * 56,
            f"  walkers={self.n_walkers}  burn={self.n_burn}  steps={self.n_steps}",
            f"  acceptance fraction = {self.acceptance_fraction:.3f}",
            f"  autocorr time = "
            + (
                ", ".join(f"{t:.0f}" for t in self.autocorr_time)
                if self.autocorr_time is not None else "unavailable"
            ),
            f"  convergence_ok = {self.convergence_ok}",
            "=" * 56,
        ])
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Hyper-parameters
def split_rhat(chain: np.ndarray) -> Optional[np.ndarray]:
    """Split-R-hat per parameter from an emcee chain ``(steps, walkers, dim)``.

    The autocorrelation criterion is the wrong tool for this posterior. The CMD
    likelihood is a thin correlated valley, so emcee mixes slowly and ``tau``
    does not settle — measured on M67 it grew from 44 at 400 steps to 202 at
    4,000, which only ever says "run longer". Meanwhile the estimates were
    stable to three decimals across 32x2000, 64x4000 and 64x8000.

    Split-R-hat asks the question that actually matters: do independent parts
    of the ensemble agree? Each walker is cut in half so between-chain variance
    can be compared with within-chain variance (Gelman & Rubin 1992; the split
    variant of Vehtari et al. 2021). Values near 1 mean the halves agree;
    > 1.05 means they do not.
    """
    if chain is None or chain.ndim != 3 or chain.shape[0] < 4:
        return None
    n_steps = chain.shape[0] // 2
    if n_steps < 2:
        return None
    # (2*walkers) half-chains of n_steps each
    halves = np.concatenate([chain[:n_steps], chain[n_steps:2 * n_steps]], axis=1)
    m = halves.shape[1]
    means = halves.mean(axis=0)                     # (m, dim)
    variances = halves.var(axis=0, ddof=1)          # (m, dim)
    W = variances.mean(axis=0)
    B = n_steps * means.var(axis=0, ddof=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        var_hat = ((n_steps - 1) / n_steps) * W + B / n_steps
        rhat = np.sqrt(var_hat / W)
    return np.where(np.isfinite(rhat), rhat, np.nan)


# ---------------------------------------------------------------------------
@dataclass
class CmdHyper:
    """Mixture hyper-parameters and CMD geometry passed into the likelihood.

    Attributes
    ----------
    f_bin, f_field :
        Mixture weights for the binary and field populations.
    field_log_density :
        ``log(1 / area_of_CMD_box)`` -- the uniform field density.  Precomputed
        once from the observed bounding box.
    n_binary_ratios :
        Number of intermediate mass-ratio companions used to approximate the
        binary ridge (in addition to the single-star and equal-mass endpoints).
    err_floor :
        Photometric softening (mag) added in quadrature to each star's error to
        keep the Gaussian from collapsing on a perfectly-on-track star.
    mh_prior, e_color_prior :
        Optional Gaussian priors ``(mean, sigma)`` on [M/H] and E(color) -- e.g.
        a spectroscopic [Fe/H] or an independent reddening estimate.
    """

    f_bin: float = 0.3
    f_field: float = 0.1
    field_log_density: float = 0.0
    n_binary_ratios: int = 4
    err_floor: float = _ERR_FLOOR
    mh_prior: Optional[Tuple[float, float]] = None
    e_color_prior: Optional[Tuple[float, float]] = None
    dm_prior: Optional[Tuple[float, float]] = None   # Gaussian (mean, sigma), e.g. Gaia parallax
    obs_weights: Optional[np.ndarray] = None   # per-star likelihood weight (turn-off emphasis)


def _compute_field_log_density(obs_c: np.ndarray, obs_m: np.ndarray) -> float:
    """Uniform field log-density = log(1 / box_area) over the observed CMD."""
    c_lo, c_hi = np.nanmin(obs_c), np.nanmax(obs_c)
    m_lo, m_hi = np.nanmin(obs_m), np.nanmax(obs_m)
    width = max(float(c_hi - c_lo), 1e-3)
    height = max(float(m_hi - m_lo), 1e-3)
    return float(-np.log(width * height))


def _compute_field_log_density_nd(obs: np.ndarray) -> float:
    """Uniform field log-density over the N-D observed bounding box.

    ``obs`` is (N_stars, N_axes).  Returns ``log(1 / Π_axis range)`` so the flat
    field density is on the same magnitude^(-N_axes) scale as the member Gaussian
    (whose per-star prefactor is Π 1/(√(2π) σ_axis)).
    """
    lo = np.nanmin(obs, axis=0)
    hi = np.nanmax(obs, axis=0)
    ranges = np.maximum(hi - lo, 1e-3)
    return float(-np.sum(np.log(ranges)))


# ---------------------------------------------------------------------------
# Multi-band isochrone evaluation (g, r, i simultaneously)
# ---------------------------------------------------------------------------
#
# The single-colour fitter (:class:`IsochroneFitterV2`) only ever carries ONE
# colour (two bands) plus a magnitude band.  To break the age/[M/H]/reddening
# degeneracy we need TWO independent colours (g-r AND r-i), which means we need
# all THREE band magnitudes (g, r, i) of the isochrone at a given (age, [M/H]).
#
# :class:`MultiBandIsochrone` reuses the same compact PARSEC subset and the same
# bilinear-in-mass interpolation strategy as ``_interpolate_bilinear``, but
# returns *absolute* magnitudes per band.  The likelihood then applies the
# distance modulus to every band and the per-band extinction A_band = R_band ·
# E(B-V), with E(B-V) derived from the sampled E(colour) via the colour's
# (R_b1 - R_b2) so the reddening of all three bands stays mutually consistent.
class MultiBandIsochrone:
    """Multi-band (e.g. g, r, i) isochrone interpolator on a PARSEC subset.

    Parameters
    ----------
    iso_data :
        Compact 2-D array.  Columns ``col_age`` and ``col_mh`` index the grid;
        ``col_mass`` is the initial mass; ``band_cols`` maps each band name to
        its absolute-magnitude column.
    col_age, col_mh, col_mass :
        Column indices for log(age), [M/H] and initial mass.
    band_cols :
        ``{"g": col_g, "r": col_r, "i": col_i}`` etc.
    R_band :
        ``{band: R_band}`` total-to-selective extinction ratio A_band/E(B-V)
        (Cardelli+1989), used to redden each band consistently.
    n_mass_pts :
        Number of points on the common mass grid for bilinear interpolation.
    """

    def __init__(
        self,
        iso_data: np.ndarray,
        *,
        col_age: int,
        col_mh: int,
        col_mass: int,
        band_cols: Dict[str, int],
        R_band: Dict[str, float],
        n_mass_pts: int = 200,
    ) -> None:
        arr = np.asarray(iso_data, dtype=float)
        if arr.ndim != 2 or arr.shape[0] == 0:
            raise ValueError("iso_data must be a non-empty 2-D array")
        self.bands: List[str] = list(band_cols.keys())
        self.band_cols = dict(band_cols)
        self.R_band = dict(R_band)
        self.col_age = int(col_age)
        self.col_mh = int(col_mh)
        self.col_mass = int(col_mass)
        self.n_mass_pts = max(50, int(n_mass_pts))

        rounded_age = np.round(arr[:, self.col_age], 2)
        rounded_mh = np.round(arr[:, self.col_mh], 2)
        self.ages = np.unique(rounded_age)
        self.metallicities = np.unique(rounded_mh)

        # Cache per grid point in EEP (file) order, NOT mass-sorted: PARSEC
        # isochrones are EEP-parametrised so the same initial mass recurs at
        # different phases; mass-sorting interleaves phases. Each entry is
        # (u, mass, {band: abs_mag}) where u is the normalised cumulative
        # arc-length along the track (0..1) used to align evolutionary phases
        # across grid corners when interpolating (see :meth:`absolute`).
        self._cache: Dict[Tuple[float, float], Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray]]] = {}
        for age in self.ages:
            for mh in self.metallicities:
                sel = (rounded_age == age) & (rounded_mh == mh)
                if sel.sum() < _MIN_ISO_POINTS:
                    continue
                sub = arr[sel]
                mass = sub[:, self.col_mass]
                cols = {b: sub[:, c] for b, c in self.band_cols.items()}
                finite = np.isfinite(mass)
                for v in cols.values():
                    finite &= np.isfinite(v)
                if finite.sum() < _MIN_ISO_POINTS:
                    continue
                mass = mass[finite]               # file/EEP order preserved
                cols = {b: v[finite] for b, v in cols.items()}
                pts = np.column_stack([cols[b] for b in self.bands])
                step = np.sqrt((np.diff(pts, axis=0) ** 2).sum(axis=1))
                s = np.concatenate([[0.0], np.cumsum(step)])
                if not np.isfinite(s[-1]) or s[-1] <= 0:
                    u = np.linspace(0.0, 1.0, len(mass))
                else:
                    u = np.maximum.accumulate(s / s[-1])
                    u = u + 1e-9 * np.arange(len(u))
                    u = u / u[-1]
                self._cache[(float(age), float(mh))] = (u, mass, cols)

    # -- grid lookup -----------------------------------------------------------
    def _bracket(self, arr: np.ndarray, value: float) -> Tuple[float, float, float]:
        idx = np.searchsorted(arr, value)
        lo = arr[max(0, idx - 1)]
        hi = arr[min(len(arr) - 1, idx)]
        if hi > lo:
            t = float(np.clip((value - lo) / (hi - lo), 0.0, 1.0))
        else:
            t = 0.0
        return float(lo), float(hi), t

    def absolute(
        self, log_age: float, mh: float
    ) -> Tuple[Optional[Dict[str, np.ndarray]], Optional[np.ndarray]]:
        """Return ({band: absolute_mag}, mass) bilinearly interpolated in (age, [M/H]).

        Falls back to the nearest grid point when the four corners are not all
        available or their mass ranges do not overlap.  Returns ``(None, None)``
        if no usable grid point exists.
        """
        age_lo, age_hi, t_age = self._bracket(self.ages, log_age)
        mh_lo, mh_hi, t_mh = self._bracket(self.metallicities, mh)

        corners = {}
        for a in (age_lo, age_hi):
            for m in (mh_lo, mh_hi):
                key = (round(a, 2), round(m, 2))
                if key in self._cache:
                    corners[(a, m)] = self._cache[key]

        if len(corners) < 4:
            return self._nearest(log_age, mh)

        w00 = (1 - t_age) * (1 - t_mh)
        w10 = t_age * (1 - t_mh)
        w01 = (1 - t_age) * t_mh
        w11 = t_age * t_mh

        # Blend the four corners at a COMMON normalised arc-length (EEP) coordinate
        # u∈[0,1], NOT at fixed initial mass. Each corner's track is resampled onto
        # u_grid via its own arc-length parametrisation, so the main sequence aligns
        # with the main sequence and the giant branch with the giant branch — the
        # interpolated track varies smoothly with (age,[M/H]) instead of fanning.
        c00 = corners[(age_lo, mh_lo)]   # (u, mass, cols)
        c10 = corners[(age_hi, mh_lo)]
        c01 = corners[(age_lo, mh_hi)]
        c11 = corners[(age_hi, mh_hi)]
        u_grid = np.linspace(0.0, 1.0, self.n_mass_pts)

        out: Dict[str, np.ndarray] = {}
        for b in self.bands:
            v00 = np.interp(u_grid, c00[0], c00[2][b])
            v10 = np.interp(u_grid, c10[0], c10[2][b])
            v01 = np.interp(u_grid, c01[0], c01[2][b])
            v11 = np.interp(u_grid, c11[0], c11[2][b])
            out[b] = w00 * v00 + w10 * v10 + w01 * v01 + w11 * v11
        m00 = np.interp(u_grid, c00[0], c00[1]); m10 = np.interp(u_grid, c10[0], c10[1])
        m01 = np.interp(u_grid, c01[0], c01[1]); m11 = np.interp(u_grid, c11[0], c11[1])
        mass_blended = w00 * m00 + w10 * m10 + w01 * m01 + w11 * m11
        return out, mass_blended

    def _nearest(
        self, log_age: float, mh: float
    ) -> Tuple[Optional[Dict[str, np.ndarray]], Optional[np.ndarray]]:
        if not self._cache:
            return None, None
        ai = int(np.argmin(np.abs(self.ages - log_age)))
        mi = int(np.argmin(np.abs(self.metallicities - mh)))
        key = (round(float(self.ages[ai]), 2), round(float(self.metallicities[mi]), 2))
        hit = self._cache.get(key)
        if hit is None:
            return None, None
        _u, mass, cols = hit
        return {b: v.copy() for b, v in cols.items()}, mass.copy()

    def apparent(
        self,
        log_age: float,
        mh: float,
        dm: float,
        e_bv: float,
    ) -> Tuple[Optional[Dict[str, np.ndarray]], Optional[np.ndarray]]:
        """Apparent per-band magnitudes: m_band = M_band + dm + R_band·E(B-V)."""
        bands, mass = self.absolute(log_age, mh)
        if bands is None:
            return None, None
        app = {}
        for b, M in bands.items():
            app[b] = M + dm + self.R_band[b] * e_bv
        return app, mass


# ---------------------------------------------------------------------------
# Multi-colour mixture hyper-parameters and likelihood
# ---------------------------------------------------------------------------
@dataclass
class CmdHyperMC:
    """Hyper-parameters + geometry for the multi-colour likelihood.

    The fit axes are ``len(colors)`` colours plus one magnitude, e.g.
    ``colors=[("g","r"), ("r","i")]`` and ``mag="g"`` gives a 3-D space
    (g-r, r-i, g).  ``e_color_to_ebv`` converts the single sampled E(colour)
    (defined on the *first* colour, matching the single-colour fitter's
    ``e_color`` = E(g-r)) into E(B-V); per-band reddening is then R_band·E(B-V).
    """

    colors: Sequence[Tuple[str, str]]
    mag: str
    iso: MultiBandIsochrone
    e_color_to_ebv: float            # E(B-V) = e_color * e_color_to_ebv
    f_bin: float = 0.3
    f_field: float = 0.1
    field_log_density: float = 0.0
    n_binary_ratios: int = 4
    err_floor: float = _ERR_FLOOR
    mh_prior: Optional[Tuple[float, float]] = None
    e_color_prior: Optional[Tuple[float, float]] = None
    dm_prior: Optional[Tuple[float, float]] = None   # Gaussian (mean, sigma), e.g. Gaia parallax
    obs_weights: Optional[np.ndarray] = None   # per-star likelihood weight (turn-off emphasis)


def _gaussian_mixture_density_nd(
    obs: np.ndarray,        # (N_obs, D)
    var: np.ndarray,        # (N_obs, D)  diagonal variances (used when cov is None)
    track: np.ndarray,      # (N_iso, D)
    track_w: np.ndarray,    # (N_iso,)
    cov: Optional[np.ndarray] = None,   # (N_obs, D, D) full covariance (overrides var)
) -> np.ndarray:
    """Σ_k w_k · N_D(obs_i | track_k, Σ_i), per observed star.

    Generalises :func:`_gaussian_mixture_density` to D axes (D = n_colours + 1).
    With ``cov=None`` the per-star covariance is diagonal ``diag(var_i)`` (fast
    path). When ``cov`` is given (e.g. colours that share a band are correlated),
    the full Mahalanobis form ``Δᵀ Σ_i⁻¹ Δ`` with prefactor ``1/√det Σ_i`` is used.
    The per-star prefactor (2π)^(-D/2) keeps the member density on the same
    magnitude^(-D) scale as the flat field term. Built in row blocks bounded to a
    few MB so memory stays flat.
    """
    n_obs, d = obs.shape
    n_iso = track.shape[0]
    if n_iso == 0 or n_obs == 0:
        return np.zeros(n_obs, dtype=float)

    norm = (2.0 * np.pi) ** (-0.5 * d)
    out = np.empty(n_obs, dtype=float)
    tr32 = np.ascontiguousarray(track, dtype=np.float32)          # (K, D)
    tw32 = np.ascontiguousarray(track_w, dtype=np.float32)        # (K,)

    if cov is None:
        # ---- diagonal fast path (one (B,K) buffer live at a time) ----
        bytes_per_row = max(1, n_iso * 4)
        rows_per_block = max(1, int(_MAX_BLOCK_BYTES // bytes_per_row))
        for start in range(0, n_obs, rows_per_block):
            stop = min(start + rows_per_block, n_obs)
            ob = obs[start:stop].astype(np.float32)               # (B, D)
            inv_v = (1.0 / var[start:stop]).astype(np.float32)    # (B, D)
            acc = None
            for dd in range(d):
                diff = ob[:, dd:dd + 1] - tr32[np.newaxis, :, dd]
                diff *= diff
                diff *= inv_v[:, dd:dd + 1]
                if acc is None:
                    acc = diff
                else:
                    acc += diff
                    del diff
            acc *= -0.5
            np.exp(acc, out=acc)
            acc *= tw32
            dens32 = acc.sum(axis=1)
            del acc
            sigma_prod = np.sqrt(np.prod(var[start:stop], axis=1))   # (B,)
            out[start:stop] = dens32.astype(np.float64) * (norm / sigma_prod)
        return out

    # ---- full-covariance path ----
    # Δᵀ Σ⁻¹ Δ = Σ_{d,e} Σ⁻¹[d,e] Δ_d Δ_e (Σ⁻¹ symmetric → diag + 2·upper).
    # We hold D difference buffers (B,K) at once, so shrink the block by ~D.
    bytes_per_row = max(1, n_iso * 4 * d)
    rows_per_block = max(1, int(_MAX_BLOCK_BYTES // bytes_per_row))
    for start in range(0, n_obs, rows_per_block):
        stop = min(start + rows_per_block, n_obs)
        ob = obs[start:stop].astype(np.float32)                   # (B, D)
        cb = np.ascontiguousarray(cov[start:stop], dtype=np.float64)  # (B, D, D)
        cinv = np.linalg.inv(cb)                                  # (B, D, D)
        _sign, logdet = np.linalg.slogdet(cb)                     # (B,)

        diffs = [ob[:, dd:dd + 1] - tr32[np.newaxis, :, dd] for dd in range(d)]  # D×(B,K)
        acc = None
        for dd in range(d):
            term = diffs[dd] * diffs[dd]
            term *= cinv[:, dd, dd][:, None].astype(np.float32)
            acc = term if acc is None else (acc + term)
        for dd in range(d):
            for ee in range(dd + 1, d):
                cross = diffs[dd] * diffs[ee]
                cross *= (2.0 * cinv[:, dd, ee])[:, None].astype(np.float32)
                acc += cross
                del cross
        acc *= -0.5
        np.exp(acc, out=acc)
        acc *= tw32
        dens32 = acc.sum(axis=1)
        del acc, diffs
        prefac = norm * np.exp(-0.5 * logdet)                     # = norm / √det Σ
        out[start:stop] = dens32.astype(np.float64) * prefac

    return out


def _color_covariance(
    err: np.ndarray,
    colors: Sequence[Tuple[str, str]],
    mag: str,
    floor2: float,
) -> Optional[np.ndarray]:
    """Full per-star covariance for fit axes ``[colours…, mag]`` that share bands.

    Colours that share a band (e.g. ``g-r`` and ``r-i`` share ``r``) are
    *correlated*; the diagonal ``var = err²`` misses the cross terms, which
    mis-weights the multi-colour density. The axis→band map ``J`` is square and
    invertible for a colour chain + magnitude, so per-band variances are recovered
    exactly from the per-axis variances and the full covariance is
    ``C = J diag(var_band) Jᵀ + floor²·I``. Its diagonal equals the old ``var``
    (so this is a pure addition of the off-diagonal correlations). Returns
    ``None`` (→ diagonal fallback) if the axis/band system is not square-invertible.
    """
    err = np.asarray(err, dtype=float)
    n_obs, d = err.shape
    axis_bands = [list(c) for c in colors] + [[mag]]
    if len(axis_bands) != d:
        return None
    band_order: List[str] = []
    for ax in axis_bands:
        for b in ax:
            if b not in band_order:
                band_order.append(b)
    if len(band_order) != d:
        return None
    jac = np.zeros((d, d), dtype=float)
    for a, (b1, b2) in enumerate(colors):
        jac[a, band_order.index(b1)] += 1.0
        jac[a, band_order.index(b2)] -= 1.0
    jac[d - 1, band_order.index(mag)] += 1.0
    m = jac ** 2
    try:
        if abs(np.linalg.det(m)) < 1e-9:
            return None
        m_inv = np.linalg.inv(m)
    except np.linalg.LinAlgError:
        return None
    var_band = np.clip((err ** 2) @ m_inv.T, 0.0, None)          # (N, d)
    cov = np.einsum("ak,nk,bk->nab", jac, var_band, jac)         # (N, d, d)
    cov += floor2 * np.eye(d)[None, :, :]
    return cov


def cmd_loglike_multicolor(
    theta: Sequence[float],
    obs: np.ndarray,        # (N_obs, D) = [color_0, color_1, ..., mag]
    err: np.ndarray,        # (N_obs, D) 1-sigma per axis
    *,
    imf_fn: ImfFn,
    hyper: CmdHyperMC,
    f_bin_override: Optional[float] = None,
    obs_weights: Optional[np.ndarray] = None,
) -> float:
    """Multi-colour generative log-likelihood (mirrors :func:`cmd_loglike`).

    The model is identical in structure to the single-colour case -- IMF-weighted
    single-star Gaussians + binary ridge + flat field -- but lives in D = n_colours
    + 1 dimensions so that two colours with different reddening / metallicity
    directions jointly constrain (age, [M/H], dm, E).
    """
    log_age, mh, dm, e_color = (float(theta[0]), float(theta[1]),
                                float(theta[2]), float(theta[3]))
    e_bv = e_color * hyper.e_color_to_ebv

    bands, mass_k = hyper.iso.apparent(log_age, mh, dm, e_bv)
    if bands is None or mass_k is None or len(mass_k) < _MIN_ISO_POINTS:
        return -np.inf

    # Build the track point cloud in fit-axis space: colours then magnitude.
    track_axes = []
    for b1, b2 in hyper.colors:
        track_axes.append(bands[b1] - bands[b2])
    track_axes.append(bands[hyper.mag])
    track = np.column_stack(track_axes)                      # (N_iso, D)
    mass_k = np.asarray(mass_k, dtype=float)

    good = np.isfinite(track).all(axis=1) & np.isfinite(mass_k)
    track = track[good]
    mass_k = mass_k[good]
    if track.shape[0] < _MIN_ISO_POINTS:
        return -np.inf

    w_k = np.asarray(imf_fn(mass_k), dtype=float)
    w_k = np.where(np.isfinite(w_k) & (w_k > 0), w_k, 0.0)
    wsum = w_k.sum()
    if wsum <= 0:
        return -np.inf
    w_k = w_k / wsum

    floor2 = hyper.err_floor ** 2
    var = np.asarray(err, dtype=float) ** 2 + floor2          # (N_obs, D)
    # Full covariance accounting for shared-band colour correlations (g-r & r-i
    # share r). Cached on the hyper object since it depends only on (err, floor),
    # not on theta — so it is built once per fit, not per MCMC step.
    cov = getattr(hyper, "_cov_cache", None)
    if cov is None:
        cov = _color_covariance(err, hyper.colors, hyper.mag, floor2)
        try:
            object.__setattr__(hyper, "_cov_cache", cov if cov is not None else False)
        except Exception:
            pass
    elif cov is False:
        cov = None

    f_bin = hyper.f_bin if f_bin_override is None else float(f_bin_override)
    f_field = hyper.f_field

    track_parts = [track]
    w_parts = [(1.0 - f_bin) * w_k]

    if f_bin > 0:
        ridge_t, ridge_w = _binary_ridge_nd(
            bands, hyper.colors, hyper.mag, mass_k, w_k, imf_fn,
            hyper.n_binary_ratios, good,
        )
        if ridge_t is not None:
            rw_sum = ridge_w.sum()
            if rw_sum > 0:
                ridge_w = (f_bin / rw_sum) * ridge_w
                track_parts.append(ridge_t)
                w_parts.append(ridge_w)
            else:
                w_parts[0] = w_k
        else:
            w_parts[0] = w_k

    track_all = np.vstack(track_parts)
    track_w = np.concatenate(w_parts)

    member = _gaussian_mixture_density_nd(obs, var, track_all, track_w)
    l_field = np.exp(hyper.field_log_density)
    l_i = (1.0 - f_field) * member + f_field * l_field

    l_i = np.where(np.isfinite(l_i) & (l_i > 0), l_i, 0.0)
    logterms = np.log(l_i + _TINY)
    if obs_weights is None:
        obs_weights = getattr(hyper, "obs_weights", None)
    if obs_weights is not None:
        # Per-star weighting. Used to balance the CMD by magnitude density so the
        # densely-populated, age-insensitive faint main sequence does not swamp the
        # sparse but age-defining turn-off / sub-giant / giant region (a Hess-diagram
        # style equalisation; see hyper.obs_weights / fit service).
        loglike = float(np.sum(np.asarray(obs_weights, dtype=float) * logterms))
    else:
        loglike = float(np.sum(logterms))
    if not np.isfinite(loglike):
        return -np.inf
    return loglike


def _binary_ridge_nd(
    bands: Dict[str, np.ndarray],
    colors: Sequence[Tuple[str, str]],
    mag: str,
    mass_k: np.ndarray,
    w_k: np.ndarray,
    imf_fn: ImfFn,
    n_ratios: int,
    good_mask: np.ndarray,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Unresolved-binary ridge in multi-band space.

    Adds a companion of mass ``q·M`` (interpolated on the track) to every primary
    in *each band's flux*, then re-derives the fit-axis colours/magnitude from the
    combined per-band magnitudes -- so the ridge geometry is exact in all bands
    (no per-colour approximation, unlike the single-colour ridge which has to
    assume which band carries the magnitude).
    """
    # All band arrays share the track's mass ordering; restrict to the same
    # finite rows used for the single-star track.
    band_names = list(bands.keys())
    mags = {b: bands[b][good_mask] for b in band_names}
    n = len(mass_k)
    if n < 2:
        return None, None

    order = np.argsort(mass_k)
    ms = mass_k[order]
    mags_s = {b: mags[b][order] for b in band_names}
    w_s = w_k[order]
    flux_s = {b: 10.0 ** (-0.4 * mags_s[b]) for b in band_names}

    qs = np.linspace(0.4, 1.0, max(1, int(n_ratios)))
    color_lists = {idx: [] for idx in range(len(colors))}
    mag_list = []
    w_list = []
    for q in qs:
        comp_mass = q * ms
        comb = {}
        for b in band_names:
            comp_m = np.interp(comp_mass, ms, mags_s[b],
                               left=mags_s[b][-1], right=mags_s[b][-1])
            comp_f = 10.0 ** (-0.4 * comp_m)
            comb[b] = -2.5 * np.log10(flux_s[b] + comp_f)
        for ci, (b1, b2) in enumerate(colors):
            color_lists[ci].append(comb[b1] - comb[b2])
        mag_list.append(comb[mag])
        w_list.append(w_s * imf_fn(comp_mass))

    axes = [np.concatenate(color_lists[ci]) for ci in range(len(colors))]
    axes.append(np.concatenate(mag_list))
    ridge_t = np.column_stack(axes)
    ridge_w = np.concatenate(w_list)
    valid = np.isfinite(ridge_t).all(axis=1) & np.isfinite(ridge_w) & (ridge_w > 0)
    if not valid.any():
        return None, None
    return ridge_t[valid], ridge_w[valid]


# ---------------------------------------------------------------------------
# Binary ridge construction
# ---------------------------------------------------------------------------
def _binary_ridge(
    c_k: np.ndarray,
    m_k: np.ndarray,
    mass_k: np.ndarray,
    w_k: np.ndarray,
    imf_fn: ImfFn,
    n_ratios: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Approximate the unresolved-binary ridge along a single isochrone.

    Approximation
    -------------
    An unresolved binary is a *primary* track star (colour ``c_k``, mag ``m_k``,
    mass ``M_k``) plus a fainter *companion* drawn from the SAME isochrone with
    mass ``q * M_k`` (mass ratio ``q in (0, 1]``).  The pair's combined flux is
    the sum of the two flux contributions, so the binary sits BRIGHTER and
    slightly REDDER than the primary.

    Rather than the full 2-D primary x companion convolution (N_iso^2 points), we
    discretise the mass ratio into ``n_ratios`` values plus the equal-mass
    endpoint and, for each primary point, place a companion by interpolating the
    track at the companion mass.  Companion flux is added in each band; we
    reconstruct per-band magnitudes from ``m_k`` (the mag band) and
    ``c_k = b1 - b2`` (assuming the listed mag band is ``b1``; the small
    colour-term error from this assumption is sub-0.02 mag and documented).

    The weight of each binary point is the IMF weight of the primary times the
    IMF weight of the companion (probability of drawing that pairing),
    integrated over the ratio grid.

    Returns
    -------
    (cb, mb, wb) : flattened binary-ridge colour, mag, and (un-normalised)
    weights.  Returns the equal-mass-only ridge if interpolation is degenerate.
    """
    n = len(mass_k)
    if n < 2:
        return c_k.copy(), m_k.copy(), w_k.copy()

    # Sort by mass for monotone interpolation of the companion.
    order = np.argsort(mass_k)
    ms = mass_k[order]
    cs = c_k[order]
    mgs = m_k[order]

    # Per-band magnitudes of the primary track.  Treat the mag band as b1 so
    # m1 = m_k and m2 = m1 - color (color = b1 - b2).  This pins the binary
    # ridge geometry; the residual error from band choice is small (<~0.02 mag)
    # and is absorbed by the photometric error floor.
    m1 = mgs
    m2 = mgs - cs
    f1 = 10.0 ** (-0.4 * m1)
    f2 = 10.0 ** (-0.4 * m2)

    # Mass-ratio grid: a few intermediate companions plus equal-mass (q=1).
    # Skip very small q (companion negligible -> ridge ~ single star).
    qs = np.linspace(0.4, 1.0, max(1, int(n_ratios)))

    cb_list = []
    mb_list = []
    wb_list = []
    for q in qs:
        comp_mass = q * ms
        # Companion per-band mags by interpolating the track at comp_mass.
        cm1 = np.interp(comp_mass, ms, m1, left=m1[-1], right=m1[-1])
        cm2 = np.interp(comp_mass, ms, m2, left=m2[-1], right=m2[-1])
        cf1 = 10.0 ** (-0.4 * cm1)
        cf2 = 10.0 ** (-0.4 * cm2)
        tot1 = f1 + cf1
        tot2 = f2 + cf2
        bm1 = -2.5 * np.log10(tot1)
        bm2 = -2.5 * np.log10(tot2)
        cb_list.append(bm1 - bm2)
        mb_list.append(bm1)
        # Weight: P(primary) * P(companion at this mass) along the ratio grid.
        comp_w = imf_fn(comp_mass)
        wb_list.append(w_k[order] * comp_w)

    cb = np.concatenate(cb_list)
    mb = np.concatenate(mb_list)
    wb = np.concatenate(wb_list)
    if not np.isfinite(wb).any() or wb.sum() <= 0:
        # Degenerate -> fall back to equal-mass binary ridge (−0.7526 mag).
        return c_k.copy(), m_k + _EQUAL_MASS_DMAG, w_k.copy()
    return cb, mb, wb


# ---------------------------------------------------------------------------
# Vectorised Gaussian mixture density
# ---------------------------------------------------------------------------
def _gaussian_mixture_density(
    obs_c: np.ndarray,
    obs_m: np.ndarray,
    var_c: np.ndarray,
    var_m: np.ndarray,
    track_c: np.ndarray,
    track_m: np.ndarray,
    track_w: np.ndarray,
) -> np.ndarray:
    """Σ_k w_k * N2(obs_i | track_k, diag(var_c_i, var_m_i)), per observed star.

    Builds the (N_obs, N_iso) Gaussian array in row blocks bounded to ~50 MB so
    memory stays flat regardless of catalogue / track size (mirrors the
    codebase's other batched numerics).  ``track_w`` should already sum to 1.

    The 2-D diagonal Gaussian is

        N2 = 1/(2π σc σm) · exp(-0.5[(Δc/σc)² + (Δm/σm)²]).

    Note the per-star normalisation 1/(2π σc σm) is included so that the field
    term (a flat density in the same magnitude-squared units) is on a consistent
    scale -- otherwise the mixture weighting between member and field would be
    arbitrary.
    """
    n_obs = obs_c.shape[0]
    n_iso = track_c.shape[0]
    if n_iso == 0 or n_obs == 0:
        return np.zeros(n_obs, dtype=float)

    # Rows per block so each (B, K) float32 block stays small.  The cap is kept
    # deliberately modest (a few MB) and we reuse buffers in-place so peak
    # working memory per block is ~2 x (B*K*4 bytes), not 5-6 temporaries.  This
    # keeps the sampler robust on memory-constrained / fragmented hosts where a
    # single large contiguous reservation can fail even with RAM nominally free.
    bytes_per_row = max(1, n_iso * 4)            # float32 kernel
    rows_per_block = max(1, int(_MAX_BLOCK_BYTES // bytes_per_row))

    norm = 1.0 / (2.0 * np.pi)  # per-star 1/(σc σm) factored in below
    out = np.empty(n_obs, dtype=float)

    # The (B, K) kernel is evaluated in float32: np.exp on float32 is ~2x faster
    # and a relative density error of ~1e-7 is negligible before the log-sum.
    tc32 = np.ascontiguousarray(track_c, dtype=np.float32)[np.newaxis, :]
    tm32 = np.ascontiguousarray(track_m, dtype=np.float32)[np.newaxis, :]
    tw32 = np.ascontiguousarray(track_w, dtype=np.float32)

    for start in range(0, n_obs, rows_per_block):
        stop = min(start + rows_per_block, n_obs)
        oc = obs_c[start:stop, np.newaxis].astype(np.float32)   # (B, 1)
        om = obs_m[start:stop, np.newaxis].astype(np.float32)
        inv_vc = (1.0 / var_c[start:stop, np.newaxis]).astype(np.float32)
        inv_vm = (1.0 / var_m[start:stop, np.newaxis]).astype(np.float32)

        # exponent = -0.5 * (dc^2/vc + dm^2/vm), built with in-place reuse so
        # only two (B, K) buffers are live at once.
        dc = oc - tc32                            # (B, K) buffer #1
        dc *= dc
        dc *= inv_vc
        dm = om - tm32                            # (B, K) buffer #2
        dm *= dm
        dm *= inv_vm
        dc += dm                                  # reuse #1; free #2 after
        del dm
        dc *= -0.5
        np.exp(dc, out=dc)                        # kernel in-place
        dc *= tw32                                # weight in-place
        dens32 = dc.sum(axis=1)                   # (B,) float32
        del dc
        # Per-star prefactor 1/(2π σc σm); σc σm = sqrt(vc·vm).
        prefac = norm / np.sqrt(var_c[start:stop] * var_m[start:stop])  # (B,) f64
        out[start:stop] = dens32.astype(np.float64) * prefac

    return out


# ---------------------------------------------------------------------------
# Per-cluster log-likelihood
# ---------------------------------------------------------------------------
def cmd_loglike(
    theta: Sequence[float],
    obs_c: np.ndarray,
    obs_m: np.ndarray,
    err_c: np.ndarray,
    err_m: np.ndarray,
    *,
    interp_fn: InterpFn,
    imf_fn: ImfFn,
    hyper: CmdHyper,
    f_bin_override: Optional[float] = None,
) -> float:
    """Marginalised generative log-likelihood of a CMD given ``theta``.

    Parameters
    ----------
    theta : (log_age, mh, dm, e_color)
    obs_c, obs_m, err_c, err_m : observed colour/mag and their 1-sigma errors.
    interp_fn : returns (iso_color, iso_mag, iso_mass) with dm + extinction
        already applied (e.g. ``IsochroneFitterV2._interpolate_isochrone``).
    imf_fn : Kroupa IMF density per mass point.
    hyper : mixture weights, field density and optional priors.
    f_bin_override : if given, use this binary fraction instead of ``hyper.f_bin``
        (used when f_bin is sampled as a nuisance parameter).

    Returns
    -------
    Σ_i log(L_i + tiny).  Returns ``-inf`` for unusable tracks or non-finite
    results so emcee rejects the proposal cleanly.
    """
    log_age, mh, dm, e_color = (float(theta[0]), float(theta[1]),
                                float(theta[2]), float(theta[3]))
    try:
        c_k, m_k, mass_k = interp_fn(log_age, mh, dm, e_color)
    except Exception:
        return -np.inf

    if c_k is None or len(c_k) < _MIN_ISO_POINTS:
        return -np.inf

    c_k = np.asarray(c_k, dtype=float)
    m_k = np.asarray(m_k, dtype=float)
    mass_k = np.asarray(mass_k, dtype=float)
    good = np.isfinite(c_k) & np.isfinite(m_k) & np.isfinite(mass_k)
    c_k, m_k, mass_k = c_k[good], m_k[good], mass_k[good]
    if len(c_k) < _MIN_ISO_POINTS:
        return -np.inf

    # IMF weights = relative number density along the track; normalise to 1.
    w_k = np.asarray(imf_fn(mass_k), dtype=float)
    w_k = np.where(np.isfinite(w_k) & (w_k > 0), w_k, 0.0)
    wsum = w_k.sum()
    if wsum <= 0:
        return -np.inf
    w_k = w_k / wsum

    # Photometric variances with a softening floor in quadrature.
    floor2 = hyper.err_floor ** 2
    var_c = np.asarray(err_c, dtype=float) ** 2 + floor2
    var_m = np.asarray(err_m, dtype=float) ** 2 + floor2

    f_bin = hyper.f_bin if f_bin_override is None else float(f_bin_override)
    f_field = hyper.f_field

    # Build the member track as a single weighted point cloud: the single-star
    # track (weight (1-f_bin)·w_k) plus the binary ridge (weight f_bin·wb_j).
    # Evaluating both in ONE Gaussian pass avoids a second (N_obs × N_iso)
    # exp() over the data and keeps the memory bound intact.
    track_c_parts = [c_k]
    track_m_parts = [m_k]
    track_w_parts = [(1.0 - f_bin) * w_k]

    if f_bin > 0:
        cb, mb, wb = _binary_ridge(c_k, m_k, mass_k, w_k, imf_fn,
                                   hyper.n_binary_ratios)
        wb = np.where(np.isfinite(wb) & np.isfinite(cb) & np.isfinite(mb)
                      & (wb > 0), wb, 0.0)
        wbsum = wb.sum()
        if wbsum > 0:
            wb = (f_bin / wbsum) * wb     # normalise ridge then scale by f_bin
            track_c_parts.append(cb)
            track_m_parts.append(mb)
            track_w_parts.append(wb)
        else:
            # Binary ridge degenerate: hand its weight back to the single track.
            track_w_parts[0] = w_k

    track_c = np.concatenate(track_c_parts)
    track_m = np.concatenate(track_m_parts)
    track_w = np.concatenate(track_w_parts)

    member = _gaussian_mixture_density(obs_c, obs_m, var_c, var_m,
                                       track_c, track_m, track_w)

    # FIELD term: flat density over the CMD bounding box.
    l_field = np.exp(hyper.field_log_density)

    l_i = (1.0 - f_field) * member + f_field * l_field

    # Guard and sum log-densities.
    l_i = np.where(np.isfinite(l_i) & (l_i > 0), l_i, 0.0)
    loglike = float(np.sum(np.log(l_i + _TINY)))
    if not np.isfinite(loglike):
        return -np.inf
    return loglike


# ---------------------------------------------------------------------------
# Prior + posterior
# ---------------------------------------------------------------------------
def _log_prior(
    theta: Sequence[float],
    bounds: FitBounds,
    hyper: CmdHyper,
    sample_f_bin: bool,
) -> float:
    """Uniform-within-bounds prior + optional Gaussian priors on mh / e_color.

    When ``sample_f_bin`` is True, ``theta`` has a 5th element constrained to
    ``[0, 1]``.
    """
    log_age, mh, dm, e_color = theta[0], theta[1], theta[2], theta[3]
    la_lo, la_hi = bounds.log_age
    mh_lo, mh_hi = bounds.metallicity
    dm_lo, dm_hi = bounds.distance_mod
    ec_lo, ec_hi = bounds.extinction_gr

    if not (la_lo <= log_age <= la_hi):
        return -np.inf
    if not (mh_lo <= mh <= mh_hi):
        return -np.inf
    if not (dm_lo <= dm <= dm_hi):
        return -np.inf
    if not (ec_lo <= e_color <= ec_hi):
        return -np.inf

    lp = 0.0
    if sample_f_bin:
        f_bin = theta[4]
        if not (0.0 <= f_bin <= 1.0):
            return -np.inf

    # Optional Gaussian priors (e.g. spectroscopic [Fe/H], reddening map).
    if hyper.mh_prior is not None:
        mu, sig = hyper.mh_prior
        if sig > 0:
            lp += -0.5 * ((mh - mu) / sig) ** 2
    if hyper.e_color_prior is not None:
        mu, sig = hyper.e_color_prior
        if sig > 0:
            lp += -0.5 * ((e_color - mu) / sig) ** 2
    dm_prior = getattr(hyper, "dm_prior", None)
    if dm_prior is not None:
        mu, sig = dm_prior
        if sig > 0:
            lp += -0.5 * ((dm - mu) / sig) ** 2
    return lp


def log_prob(
    theta: Sequence[float],
    obs_c: np.ndarray,
    obs_m: np.ndarray,
    err_c: np.ndarray,
    err_m: np.ndarray,
    *,
    bounds: FitBounds,
    interp_fn: InterpFn,
    imf_fn: ImfFn,
    hyper: CmdHyper,
    sample_f_bin: bool = False,
) -> float:
    """log posterior = log_prior + cmd_loglike.  ``-inf`` outside bounds."""
    lp = _log_prior(theta, bounds, hyper, sample_f_bin)
    if not np.isfinite(lp):
        return -np.inf
    f_bin_override = float(theta[4]) if sample_f_bin else None
    ll = cmd_loglike(
        theta[:4], obs_c, obs_m, err_c, err_m,
        interp_fn=interp_fn, imf_fn=imf_fn, hyper=hyper,
        f_bin_override=f_bin_override,
    )
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll


def log_prob_multicolor(
    theta: Sequence[float],
    obs: np.ndarray,
    err: np.ndarray,
    *,
    bounds: FitBounds,
    imf_fn: ImfFn,
    hyper: CmdHyperMC,
    sample_f_bin: bool = False,
) -> float:
    """log posterior for the multi-colour likelihood.  ``-inf`` outside bounds."""
    lp = _log_prior(theta, bounds, hyper, sample_f_bin)
    if not np.isfinite(lp):
        return -np.inf
    f_bin_override = float(theta[4]) if sample_f_bin else None
    ll = cmd_loglike_multicolor(
        theta[:4], obs, err,
        imf_fn=imf_fn, hyper=hyper, f_bin_override=f_bin_override,
    )
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll


# ---------------------------------------------------------------------------
# Sampler driver
# ---------------------------------------------------------------------------
def _percentiles(samples: np.ndarray) -> Tuple[float, float, float]:
    p16, p50, p84 = np.percentile(samples, [16, 50, 84])
    return (float(p16), float(p50), float(p84))


def fit_isochrone_mcmc(
    obs_c: np.ndarray,
    obs_m: np.ndarray,
    err_c: np.ndarray,
    err_m: np.ndarray,
    *,
    bounds: FitBounds,
    interp_fn: InterpFn,
    imf_fn: ImfFn,
    n_walkers: int = 32,
    n_burn: int = 400,
    n_steps: int = 1500,
    init_from_gridscan: Optional[Sequence[float]] = None,
    priors: Optional[Dict[str, Tuple[float, float]]] = None,
    f_bin: float = 0.3,
    f_field: float = 0.1,
    sample_f_bin: bool = False,
    n_binary_ratios: int = 4,
    err_floor: float = _ERR_FLOOR,
    R_color: Optional[float] = None,
    seed: int = 1234,
    progress: bool = False,
    progress_cb: Optional[Callable[[float, str], None]] = None,
    max_chain_samples: int = 4000,
    # --- multi-colour option ------------------------------------------------
    colors: Optional[Sequence[Tuple[str, str]]] = None,
    mag: Optional[str] = None,
    multiband_iso: Optional[MultiBandIsochrone] = None,
    obs_multicolor: Optional[np.ndarray] = None,
    err_multicolor: Optional[np.ndarray] = None,
    obs_weights: Optional[np.ndarray] = None,
) -> McmcResult:
    """Sample the isochrone posterior with emcee and return credible intervals.

    Parameters
    ----------
    obs_c, obs_m, err_c, err_m : observed CMD and 1-sigma photometric errors
        (single-colour path).  Ignored when ``colors`` is given.
    bounds : :class:`FitBounds` defining the uniform prior box.
    interp_fn, imf_fn : isochrone interpolator and IMF density (borrow from
        :class:`IsochroneFitterV2`).
    n_walkers, n_burn, n_steps : ensemble size and chain lengths.
    init_from_gridscan : ``(log_age, mh, dm, e_color)`` best grid-scan fit used
        to seed a tight walker ball; falls back to the bounds midpoint.
    priors : optional ``{"mh": (mu, sig), "e_color": (mu, sig)}`` Gaussian priors.
    f_bin, f_field : default mixture weights.  ``sample_f_bin=True`` promotes the
        binary fraction to a sampled 5th parameter (uniform on [0, 1]).
    R_color : if given, the colour-extinction coefficient ``R1 - R2`` so E(B-V)
        percentiles can be reported (E(B-V) = E(color) / R_color).
    seed : RNG seed for walker init + sampler reproducibility.
    colors, mag, multiband_iso, obs_multicolor, err_multicolor :
        Multi-colour mode.  When ``colors`` (e.g. ``[("g","r"), ("r","i")]``) is
        provided, the fit uses :func:`cmd_loglike_multicolor` over the D-axis
        space (one axis per colour + ``mag``) instead of the single-colour
        likelihood.  ``multiband_iso`` is a :class:`MultiBandIsochrone` providing
        all band magnitudes; ``obs_multicolor`` / ``err_multicolor`` are
        ``(N_stars, D)`` observed values + errors in the SAME axis order
        (colours then mag).  The sampled E(colour) is interpreted on the FIRST
        colour, matching the single-colour fitter's E(g-r).

    Returns
    -------
    :class:`McmcResult` with 16/50/84 percentiles, a downsampled flat chain,
    acceptance fraction, autocorrelation time and a ``convergence_ok`` flag.
    """
    import emcee  # local import keeps module import cheap and Qt-free

    multicolor = colors is not None
    if multicolor:
        if multiband_iso is None or mag is None:
            raise ValueError("multi-colour mode needs multiband_iso and mag")
        if obs_multicolor is None or err_multicolor is None:
            raise ValueError("multi-colour mode needs obs_multicolor/err_multicolor")
        obs_mc = np.asarray(obs_multicolor, dtype=float)
        err_mc = np.asarray(err_multicolor, dtype=float)
        if obs_mc.ndim != 2 or obs_mc.shape[1] != len(colors) + 1:
            raise ValueError(
                f"obs_multicolor must be (N, {len(colors) + 1}); got {obs_mc.shape}")
        finite = np.isfinite(obs_mc).all(axis=1) & np.isfinite(err_mc).all(axis=1)
        obs_mc, err_mc = obs_mc[finite], err_mc[finite]
        if obs_mc.shape[0] < _MIN_ISO_POINTS:
            raise ValueError("Too few finite stars to run multi-colour MCMC")
        err_mc = np.clip(err_mc, 1e-3, None)
        ow_mc = None
        if obs_weights is not None:
            ow_mc = np.asarray(obs_weights, dtype=float)[finite]

        # E(colour) is on the first colour; convert to E(B-V) via that colour's
        # (R_b1 - R_b2) so all bands redden consistently.
        b1, b2 = colors[0]
        r_diff = multiband_iso.R_band[b1] - multiband_iso.R_band[b2]
        e_color_to_ebv = 1.0 / r_diff if abs(r_diff) > 1e-9 else 0.0

        mh_prior = priors.get("mh") if priors else None
        e_color_prior = priors.get("e_color") if priors else None
        dm_prior = priors.get("dm") if priors else None
        hyper_mc = CmdHyperMC(
            colors=list(colors),
            mag=mag,
            iso=multiband_iso,
            e_color_to_ebv=e_color_to_ebv,
            f_bin=float(f_bin),
            f_field=float(f_field),
            field_log_density=_compute_field_log_density_nd(obs_mc),
            n_binary_ratios=int(n_binary_ratios),
            err_floor=float(err_floor),
            mh_prior=mh_prior,
            e_color_prior=e_color_prior,
            dm_prior=dm_prior,
            obs_weights=ow_mc,
        )
    else:
        obs_c = np.asarray(obs_c, dtype=float)
        obs_m = np.asarray(obs_m, dtype=float)
        err_c = np.asarray(err_c, dtype=float)
        err_m = np.asarray(err_m, dtype=float)

        finite = (np.isfinite(obs_c) & np.isfinite(obs_m)
                  & np.isfinite(err_c) & np.isfinite(err_m))
        obs_c, obs_m = obs_c[finite], obs_m[finite]
        err_c, err_m = err_c[finite], err_m[finite]
        if obs_c.size < _MIN_ISO_POINTS:
            raise ValueError("Too few finite stars to run MCMC")

        err_c = np.clip(err_c, 1e-3, None)
        err_m = np.clip(err_m, 1e-3, None)

        mh_prior = priors.get("mh") if priors else None
        e_color_prior = priors.get("e_color") if priors else None
        dm_prior = priors.get("dm") if priors else None
        hyper = CmdHyper(
            f_bin=float(f_bin),
            f_field=float(f_field),
            field_log_density=_compute_field_log_density(obs_c, obs_m),
            n_binary_ratios=int(n_binary_ratios),
            err_floor=float(err_floor),
            mh_prior=mh_prior,
            e_color_prior=e_color_prior,
            dm_prior=dm_prior,
        )

    n_dim = 5 if sample_f_bin else 4
    e_label = f"E({colors[0][0]}-{colors[0][1]})" if multicolor else "E(color)"
    labels = ["log_age", "[M/H]", "(m-M)", e_label]
    if sample_f_bin:
        labels.append("f_bin")

    # --- walker initialisation ------------------------------------------------
    rng = np.random.default_rng(seed)
    box = bounds.to_list()  # [(lo,hi)] for 4 params
    if init_from_gridscan is not None and len(init_from_gridscan) >= 4:
        center = np.array([float(v) for v in init_from_gridscan[:4]], dtype=float)
    else:
        center = np.array([(lo + hi) * 0.5 for lo, hi in box], dtype=float)
    if sample_f_bin:
        center = np.append(center, float(f_bin))
        box = box + [(0.0, 1.0)]

    # Clamp the (grid-scan) seed into the bounds first: with tightened prior
    # windows (parallax dm, reddening, [M/H]) the seed can land OUTSIDE a window,
    # which would push every walker onto the same bound edge -> a constant column
    # -> linearly-dependent walkers -> emcee "large condition number" error.
    lo_arr = np.array([lo for lo, hi in box], dtype=float)
    hi_arr = np.array([hi for lo, hi in box], dtype=float)
    center = np.clip(center, lo_arr, hi_arr)

    # Tight Gaussian ball: 2% of each parameter's prior width (with an absolute
    # floor so no dimension is degenerate).
    widths = hi_arr - lo_arr
    spread = np.maximum(0.02 * widths, 1e-4)
    p0 = center[np.newaxis, :] + spread[np.newaxis, :] * rng.standard_normal((n_walkers, n_dim))
    # Clip just inside bounds; if a dimension still collapsed onto a single value
    # (seed at the very edge of a tight window), re-seed it uniformly so the
    # walkers remain linearly independent.
    for d, (lo, hi) in enumerate(box):
        w = max(hi - lo, 1e-9)
        eps = 1e-3 * w
        p0[:, d] = np.clip(p0[:, d], lo + eps, hi - eps)
        if float(np.ptp(p0[:, d])) < 1e-9:
            p0[:, d] = rng.uniform(lo + eps, hi - eps, n_walkers)

    # The CMD posterior is narrow and strongly correlated (the age-metallicity-
    # reddening-distance degeneracy is a thin diagonal valley). emcee's default
    # StretchMove mixes poorly there (acceptance collapses to a few percent, so
    # walkers stay stuck at the seed). A DE/DESnooker move mixture is the standard
    # remedy for correlated posteriors (Foreman-Mackey 2013; emcee docs) and is a
    # general, data-independent fix — the same moves are used for every target.
    mcmc_moves = [
        (emcee.moves.DEMove(), 0.8),
        (emcee.moves.DESnookerMove(), 0.2),
    ]

    # emcee's ensemble moves draw from the global NumPy RNG, which is NOT seeded by
    # our local ``rng`` — so the chain (and the test/credible intervals) depended on
    # whatever consumed the global state earlier (test-order-dependent flakiness, and
    # non-reproducible science). Seed the global RNG from ``seed`` so a given fit is
    # reproducible regardless of process history.
    np.random.seed(int(seed) % (2 ** 32))

    if multicolor:
        sampler = emcee.EnsembleSampler(
            n_walkers,
            n_dim,
            log_prob_multicolor,
            args=(obs_mc, err_mc),
            kwargs=dict(
                bounds=bounds,
                imf_fn=imf_fn,
                hyper=hyper_mc,
                sample_f_bin=sample_f_bin,
            ),
            moves=mcmc_moves,
        )
    else:
        sampler = emcee.EnsembleSampler(
            n_walkers,
            n_dim,
            log_prob,
            args=(obs_c, obs_m, err_c, err_m),
            kwargs=dict(
                bounds=bounds,
                interp_fn=interp_fn,
                imf_fn=imf_fn,
                hyper=hyper,
                sample_f_bin=sample_f_bin,
            ),
            moves=mcmc_moves,
        )

    # Run via sampler.sample() so per-iteration progress can be reported to the
    # GUI (the long MCMC otherwise shows no movement). total = burn + production.
    _total = max(1, int(n_burn) + int(n_steps))

    def _run_phase(initial, n_iter, base, label):
        n_iter = int(n_iter)
        every = max(1, n_iter // 25)
        last = initial
        i = 0
        for last in sampler.sample(initial, iterations=n_iter, progress=progress):
            i += 1
            if progress_cb is not None and (i % every == 0 or i == n_iter):
                try:
                    progress_cb((base + i) / _total, label)
                except Exception:
                    pass
        return last

    # Burn-in. If emcee rejects the initial ensemble as linearly dependent
    # (degenerate seed under very tight prior windows), re-seed every walker
    # uniformly across the bounds once and retry — so the user never sees the
    # cryptic "large condition number" error.
    try:
        state = _run_phase(p0, n_burn, 0, "MCMC burn-in")
    except ValueError as exc:
        if "condition number" not in str(exc).lower():
            raise
        p0 = np.empty((n_walkers, n_dim), dtype=float)
        for d, (lo, hi) in enumerate(box):
            eps = 1e-3 * max(hi - lo, 1e-9)
            p0[:, d] = rng.uniform(lo + eps, hi - eps, n_walkers)
        sampler.reset()
        state = _run_phase(p0, n_burn, 0, "MCMC burn-in")
    sampler.reset()
    # Production.
    _run_phase(state, n_steps, n_burn, "MCMC sampling")

    # --- diagnostics ----------------------------------------------------------
    acc = float(np.mean(sampler.acceptance_fraction))
    tau = None
    thin = 1
    convergence_ok = False
    # Why the run is or is not judged converged. A bare ``except: tau = None``
    # used to swallow this, so every fit reported "not converged" with no way
    # to tell an unconverged chain from a diagnostic that never ran.
    convergence_detail = ""
    try:
        tau = sampler.get_autocorr_time(tol=0)
        tau = np.asarray(tau, dtype=float)
        if np.all(np.isfinite(tau)) and np.all(tau > 0):
            max_tau = float(np.max(tau))
            thin = max(1, int(max_tau // 2))
            # Standard emcee convergence heuristic: chain >= 50 autocorr times.
            long_enough = bool(n_steps > 50 * max_tau)
            acc_ok = 0.15 <= acc <= 0.7
            convergence_ok = long_enough and acc_ok
            if not long_enough:
                convergence_detail = (
                    f"chain is {n_steps} steps but the autocorrelation time is "
                    f"{max_tau:.0f}; {int(50 * max_tau)} are needed")
            elif not acc_ok:
                convergence_detail = f"acceptance {acc:.3f} outside 0.15-0.70"
        else:
            convergence_detail = (
                "autocorrelation time is not finite/positive "
                f"({np.array2string(tau, precision=1)})")
            tau = None
    except Exception as exc:  # noqa: BLE001 - reported, not hidden
        convergence_detail = f"{type(exc).__name__}: {exc}"
        tau = None

    # Second opinion. tau only ever says "run longer" on this posterior, so
    # report whether independent halves of the ensemble actually agree.
    rhat = None
    try:
        rhat = split_rhat(sampler.get_chain())
    except Exception:
        rhat = None
    if rhat is not None and np.all(np.isfinite(rhat)):
        worst = float(np.max(rhat))
        agreed = worst < 1.05
        convergence_detail = (
            f"{convergence_detail}; split-R-hat max {worst:.3f} "
            f"({'halves agree' if agreed else 'halves disagree'})"
        ).lstrip("; ")
        # A chain whose halves agree and whose acceptance is healthy is usable
        # even when the autocorrelation criterion cannot be satisfied.
        if not convergence_ok and agreed and 0.15 <= acc <= 0.7:
            convergence_ok = True

    discard = 0  # burn-in already discarded via reset()
    flat = sampler.get_chain(discard=discard, thin=thin, flat=True)

    # Downsample the flat chain for corner plots / storage.
    if flat.shape[0] > max_chain_samples:
        sel = rng.choice(flat.shape[0], size=max_chain_samples, replace=False)
        flat_store = flat[np.sort(sel)]
    else:
        flat_store = flat

    # --- percentiles ----------------------------------------------------------
    log_age_p = _percentiles(flat[:, 0])
    mh_p = _percentiles(flat[:, 1])
    dm_p = _percentiles(flat[:, 2])
    e_color_p = _percentiles(flat[:, 3])
    age_gyr_p = _percentiles(10.0 ** (flat[:, 0] - 9.0))
    distance_pc_p = _percentiles(10.0 ** (1.0 + flat[:, 2] / 5.0))
    e_bv_p = None
    if R_color is not None and abs(R_color) > 1e-9:
        e_bv_p = _percentiles(flat[:, 3] / float(R_color))
    f_bin_p = _percentiles(flat[:, 4]) if sample_f_bin else None

    median_theta = np.array([log_age_p[1], mh_p[1], dm_p[1], e_color_p[1]],
                            dtype=float)

    return McmcResult(
        log_age_p=log_age_p,
        mh_p=mh_p,
        dm_p=dm_p,
        e_color_p=e_color_p,
        age_gyr_p=age_gyr_p,
        distance_pc_p=distance_pc_p,
        e_bv_p=e_bv_p,
        f_bin_p=f_bin_p,
        flat_chain=flat_store,
        labels=labels,
        acceptance_fraction=acc,
        convergence_detail=convergence_detail,
        rhat=rhat,
        autocorr_time=tau,
        n_burn=n_burn,
        n_steps=n_steps,
        n_walkers=n_walkers,
        n_dim=n_dim,
        convergence_ok=convergence_ok,
        sample_binary_fraction=sample_f_bin,
        median_theta=median_theta,
    )
