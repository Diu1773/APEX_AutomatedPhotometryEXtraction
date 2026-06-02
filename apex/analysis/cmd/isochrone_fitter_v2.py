"""
Isochrone Fitting Module v2

Fits (log_age, [M/H], distance_mod, E(color)) to an observed CMD using:
- Nearest grid-point lookup (argmin on precomputed age/metallicity arrays)
- KD-tree nearest-neighbour distance in error-normalised CMD space
- chi²-like objective: sum of squared (Δcol/σ_col)² + (Δmag/σ_mag)² for the
  closest ``fit_fraction`` of stars (robust against field-star contamination)
- Optional EM membership weighting across multiple fit rounds

Note: bilinear interpolation between grid points and perpendicular-distance
computation are not yet implemented.  ``_imf_weight`` is defined but unused
(kept for future IMF-likelihood extension).

Author: KNUEMAO Pipeline
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional, Tuple, Callable, Dict, List

import numpy as np
from scipy.optimize import minimize, differential_evolution
from scipy.spatial import cKDTree


class FitMode(Enum):
    """Fitting mode selection"""
    FAST = "fast"
    HESSIAN = "hessian"
    GRID_SCAN = "grid_scan"


@dataclass
class FitBounds:
    """Parameter bounds for fitting"""
    log_age: Tuple[float, float] = (8.0, 10.0)
    metallicity: Tuple[float, float] = (-0.5, 0.5)
    distance_mod: Tuple[float, float] = (8.0, 13.0)
    extinction_gr: Tuple[float, float] = (0.0, 0.5)

    def to_list(self):
        return [self.log_age, self.metallicity,
                self.distance_mod, self.extinction_gr]


@dataclass
class FitResult:
    """Isochrone fitting result"""
    log_age: float
    metallicity: float
    distance_mod: float
    extinction_gr: float

    log_age_err: Optional[float] = None
    metallicity_err: Optional[float] = None
    distance_mod_err: Optional[float] = None
    extinction_gr_err: Optional[float] = None

    chi2: float = 0.0
    reduced_chi2: float = 0.0
    n_stars: int = 0
    n_params: int = 4

    fit_mode: str = "fast"
    elapsed_sec: float = 0.0
    converged: bool = True
    extinction_label: str = "E(g-r)"
    color_label: str = "g-r"
    mag_label: str = "g"

    @property
    def age_gyr(self) -> float:
        return 10 ** (self.log_age - 9)

    @property
    def distance_pc(self) -> float:
        return 10 ** (1 + self.distance_mod / 5)

    def summary(self) -> str:
        lines = [
            "=" * 50,
            f"Isochrone Fit Result ({self.fit_mode.upper()} mode)",
            "=" * 50,
            "",
            "Best-fit Parameters:",
        ]

        if self.log_age_err:
            lines.append(f"  log(Age) = {self.log_age:.3f} ± {self.log_age_err:.3f}")
            age_gyr_err = self.age_gyr * np.log(10) * self.log_age_err
            lines.append(f"           = {self.age_gyr:.2f} ± {age_gyr_err:.2f} Gyr")
        else:
            lines.append(f"  log(Age) = {self.log_age:.3f}  ({self.age_gyr:.2f} Gyr)")

        if self.metallicity_err:
            lines.append(f"  [M/H]    = {self.metallicity:.3f} ± {self.metallicity_err:.3f}")
        else:
            lines.append(f"  [M/H]    = {self.metallicity:.3f}")

        if self.distance_mod_err:
            lines.append(f"  (m-M)    = {self.distance_mod:.3f} ± {self.distance_mod_err:.3f}")
            dist_pc_err = self.distance_pc * np.log(10) / 5 * self.distance_mod_err
            lines.append(f"           = {self.distance_pc:.0f} ± {dist_pc_err:.0f} pc")
        else:
            lines.append(f"  (m-M)    = {self.distance_mod:.3f}  ({self.distance_pc:.0f} pc)")

        e_label = self.extinction_label or "E(color)"
        if self.extinction_gr_err:
            lines.append(f"  {e_label:<8s} = {self.extinction_gr:.4f} ± {self.extinction_gr_err:.4f}")
        else:
            lines.append(f"  {e_label:<8s} = {self.extinction_gr:.4f}")

        if self.color_label or self.mag_label:
            lines.append(f"  CMD axes : color={self.color_label}, mag={self.mag_label}")

        lines.extend([
            "",
            "Fit Quality:",
            f"  chi²     = {self.chi2:.1f}",
            f"  chi²/dof = {self.reduced_chi2:.3f}",
            f"  N_stars  = {self.n_stars}",
            "",
            f"Elapsed: {self.elapsed_sec:.2f} sec",
            "=" * 50,
        ])

        return "\n".join(lines)


@dataclass
class GridScanResult:
    """Result from grid scan fitting, including chi2 landscape."""
    best_fit: FitResult
    chi2_grid: np.ndarray          # shape (n_ages, n_mhs)
    dm_grid: np.ndarray            # shape (n_ages, n_mhs)
    egr_grid: np.ndarray           # shape (n_ages, n_mhs)
    grid_ages: np.ndarray          # 1D sorted ages within bounds
    grid_mhs: np.ndarray           # 1D sorted metallicities within bounds
    n_evaluated: int = 0
    elapsed_sec: float = 0.0

    def summary(self) -> str:
        lines = [
            "=" * 50,
            "Grid Scan Result",
            "=" * 50,
            "",
            f"Grid: {len(self.grid_ages)} ages x {len(self.grid_mhs)} [M/H]"
            f" = {len(self.grid_ages) * len(self.grid_mhs)} cells",
            f"Evaluated: {self.n_evaluated} cells",
            f"Elapsed: {self.elapsed_sec:.2f} sec",
            "",
        ]
        lines.append(self.best_fit.summary())
        return "\n".join(lines)


class IsochroneFitterV2:
    """KD-tree based isochrone fitter (nearest grid-point, error-normalised CMD space)."""

    # Extinction coefficients R_band = A_band / E(B-V) (Cardelli+1989)
    # Shared with step12 via EXTINCTION_R; kept here for internal use.
    SDSS_R = {"U": 4.902, "B": 4.035, "V": 3.116, "R": 2.634, "I": 1.903,
              "g": 3.303, "r": 2.285, "i": 1.698, "z": 1.263}

    def __init__(
        self,
        isochrone_file: str | Path,
        col_mh: int = 1,
        col_age: int = 2,
        col_g: int = 29,
        col_r: int = 30,
        col_mag: Optional[int] = None,
        col_mass: int = 5,
        fit_fraction: float = 0.7,
        R_band1: float = 3.303,
        R_band2: float = 2.285,
        R_mag: Optional[float] = None,
        iso_data: Optional[np.ndarray] = None,
        use_imf_weight: bool = False,
        use_bilinear: bool = True,
        bilinear_n_points: int = 200,
    ):
        self.iso_file = Path(isochrone_file)
        self.COL_MH = col_mh
        self.COL_AGE = col_age
        # band1, band2 define color = band1 - band2
        self.COL_G = col_g
        self.COL_R = col_r
        # mag band can be independent from color bands.
        self.COL_MAG = int(col_mag) if col_mag is not None else int(col_g)
        self.COL_MASS = col_mass
        self.fit_fraction    = fit_fraction
        self.use_imf_weight  = bool(use_imf_weight)
        self.use_bilinear    = bool(use_bilinear)
        self.bilinear_n_pts  = max(50, int(bilinear_n_points))
        self.R_1 = R_band1
        self.R_2 = R_band2
        self.R_MAG = float(R_mag) if R_mag is not None else float(R_band1)

        self.iso_data = self._load_isochrone() if iso_data is None else self._clean_isochrone_data(iso_data)
        self._validate_isochrone_columns()

        # Build grid structure
        self.ages = np.unique(np.round(self.iso_data[:, self.COL_AGE], 2))
        self.metallicities = np.unique(np.round(self.iso_data[:, self.COL_MH], 2))

        # Pre-compute isochrone curves for each (age, mh) combination
        self._iso_cache: Dict[Tuple[float, float], np.ndarray] = {}
        self._build_iso_cache()

        # Pre-extracted column cache: key -> (raw_color, raw_mag, mass)
        # Avoids repeated 2D array slicing in the inner objective loop.
        self._col_cache: Dict[Tuple[float, float], Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        self._build_col_cache()

        # Fitting state
        self.progress_callback: Optional[Callable[[float, str], None]] = None
        self._fit_obs: Optional[np.ndarray] = None
        self._fit_err: Optional[np.ndarray] = None
        self._fit_weight: Optional[np.ndarray] = None
        self._fit_iteration: int = 0
        self._fit_max_iter: int = 100

    def _load_isochrone(self) -> np.ndarray:
        data = np.genfromtxt(self.iso_file, comments='#')
        return self._clean_isochrone_data(data)

    @staticmethod
    def _clean_isochrone_data(data: np.ndarray) -> np.ndarray:
        arr = np.asarray(data, dtype=float)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        if arr.ndim != 2:
            raise ValueError("Isochrone data must be a 2D numeric array")
        if arr.size == 0:
            return arr
        return arr[~np.isnan(arr).any(axis=1)]

    def _validate_isochrone_columns(self) -> None:
        if self.iso_data.size == 0:
            raise ValueError("Isochrone file is empty")
        max_col = max(self.COL_MH, self.COL_AGE, self.COL_G, self.COL_R, self.COL_MAG, self.COL_MASS)
        if self.iso_data.shape[1] <= max_col:
            raise ValueError(
                f"Isochrone data has {self.iso_data.shape[1]} columns, "
                f"but column index {max_col} is required"
            )

    def _build_iso_cache(self):
        """Pre-cache isochrone data for each grid point."""
        if len(self.iso_data) == 0:
            return

        rounded_ages = np.round(self.iso_data[:, self.COL_AGE], 2)
        rounded_mhs = np.round(self.iso_data[:, self.COL_MH], 2)
        order = np.lexsort((self.iso_data[:, self.COL_MAG], rounded_mhs, rounded_ages))

        self.iso_data = self.iso_data[order]
        rounded_ages = rounded_ages[order]
        rounded_mhs = rounded_mhs[order]

        group_breaks = (
            (rounded_ages[1:] != rounded_ages[:-1])
            | (rounded_mhs[1:] != rounded_mhs[:-1])
        )
        bounds = np.concatenate(([0], np.flatnonzero(group_breaks) + 1, [len(self.iso_data)]))

        for start, stop in zip(bounds[:-1], bounds[1:]):
            if stop - start > 10:
                key = (float(rounded_ages[start]), float(rounded_mhs[start]))
                self._iso_cache[key] = self.iso_data[start:stop]

    def _build_col_cache(self) -> None:
        """Pre-extract band columns for every cached (age, mh) grid point.

        Stores (raw_color, raw_mag, mass) with validity masking already applied,
        so ``_interpolate_isochrone`` only needs to apply dm/E(color) offsets.
        """
        for key, iso_data in self._iso_cache.items():
            b1   = iso_data[:, self.COL_G].copy()
            b2   = iso_data[:, self.COL_R].copy()
            mag  = iso_data[:, self.COL_MAG].copy()
            mass = iso_data[:, self.COL_MASS].copy()
            valid = np.isfinite(b1) & np.isfinite(b2) & np.isfinite(mag) & np.isfinite(mass)
            if valid.sum() >= 10:
                self._col_cache[key] = (b1[valid] - b2[valid], mag[valid], mass[valid])

    def _get_nearby_isochrones(self, log_age: float, mh: float) -> List[Tuple[float, float, np.ndarray]]:
        """Get up to 4 nearest isochrones for interpolation"""
        # Find bracketing ages
        age_idx = np.searchsorted(self.ages, log_age)
        age_lo = self.ages[max(0, age_idx - 1)]
        age_hi = self.ages[min(len(self.ages) - 1, age_idx)]

        # Find bracketing metallicities
        mh_idx = np.searchsorted(self.metallicities, mh)
        mh_lo = self.metallicities[max(0, mh_idx - 1)]
        mh_hi = self.metallicities[min(len(self.metallicities) - 1, mh_idx)]

        nearby = []
        for a in [age_lo, age_hi]:
            for m in [mh_lo, mh_hi]:
                key = (round(a, 2), round(m, 2))
                if key in self._iso_cache:
                    nearby.append((a, m, self._iso_cache[key]))

        return nearby

    def _apply_extinction(
        self, raw_color: np.ndarray, raw_mag: np.ndarray,
        dm: float, e_gr: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Apply distance modulus and Cardelli extinction to raw columns."""
        E_BV  = e_gr / (self.R_1 - self.R_2)
        A_mag = self.R_MAG * E_BV
        return raw_color + e_gr, raw_mag + dm + A_mag

    def _interpolate_nearest(
        self, log_age: float, mh: float, dm: float, e_gr: float
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Nearest grid-point lookup (fast fallback)."""
        age_idx = np.argmin(np.abs(self.ages - log_age))
        mh_idx  = np.argmin(np.abs(self.metallicities - mh))
        key = (round(float(self.ages[age_idx]), 2),
               round(float(self.metallicities[mh_idx]), 2))
        raw = self._col_cache.get(key)
        if raw is None:
            return np.array([]), np.array([]), np.array([])
        raw_color, raw_mag, mass = raw
        color, mag = self._apply_extinction(raw_color, raw_mag, dm, e_gr)
        return color, mag, mass

    def _interpolate_bilinear(
        self, log_age: float, mh: float, dm: float, e_gr: float
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Bilinear interpolation across the 4 bracketing (age, [M/H]) grid points.

        Uses initial stellar mass as the common parameterisation axis so that
        evolutionary stage is matched correctly between isochrones at different
        ages/metallicities.  Falls back to nearest-grid if any corner is absent
        or the mass ranges do not overlap sufficiently.
        """
        # ── Bracketing grid points ──────────────────────────────────────────
        ai = np.searchsorted(self.ages, log_age)
        age_lo = self.ages[max(0, ai - 1)]
        age_hi = self.ages[min(len(self.ages) - 1, ai)]

        mi = np.searchsorted(self.metallicities, mh)
        mh_lo = self.metallicities[max(0, mi - 1)]
        mh_hi = self.metallicities[min(len(self.metallicities) - 1, mi)]

        # Bilinear weights in (age, mh) space
        t_age = float(np.clip(
            (log_age - age_lo) / (age_hi - age_lo) if age_hi > age_lo else 0.0,
            0.0, 1.0,
        ))
        t_mh = float(np.clip(
            (mh - mh_lo) / (mh_hi - mh_lo) if mh_hi > mh_lo else 0.0,
            0.0, 1.0,
        ))

        # ── Load all 4 corners ──────────────────────────────────────────────
        corners: dict[tuple, tuple] = {}
        for a in (age_lo, age_hi):
            for m in (mh_lo, mh_hi):
                key = (round(float(a), 2), round(float(m), 2))
                raw = self._col_cache.get(key)
                if raw is not None:
                    corners[(a, m)] = raw

        if len(corners) < 4:
            return self._interpolate_nearest(log_age, mh, dm, e_gr)

        # ── Common mass grid (intersection of all 4 mass ranges) ───────────
        mass_lo = max(raw[2][0]  for raw in corners.values())
        mass_hi = min(raw[2][-1] for raw in corners.values())
        if mass_hi <= mass_lo + 1e-6:
            return self._interpolate_nearest(log_age, mh, dm, e_gr)

        mass_grid = np.linspace(mass_lo, mass_hi, self.bilinear_n_pts)

        def _resample(raw: tuple) -> Tuple[np.ndarray, np.ndarray]:
            rc, rm, ms = raw
            order = np.argsort(ms)
            return (np.interp(mass_grid, ms[order], rc[order]),
                    np.interp(mass_grid, ms[order], rm[order]))

        c_lolo, g_lolo = _resample(corners[(age_lo, mh_lo)])
        c_hilo, g_hilo = _resample(corners[(age_hi, mh_lo)])
        c_lohi, g_lohi = _resample(corners[(age_lo, mh_hi)])
        c_hihi, g_hihi = _resample(corners[(age_hi, mh_hi)])

        # ── Bilinear combination ────────────────────────────────────────────
        w00 = (1 - t_age) * (1 - t_mh)
        w10 = t_age        * (1 - t_mh)
        w01 = (1 - t_age) * t_mh
        w11 = t_age        * t_mh

        raw_color = w00 * c_lolo + w10 * c_hilo + w01 * c_lohi + w11 * c_hihi
        raw_mag   = w00 * g_lolo + w10 * g_hilo + w01 * g_lohi + w11 * g_hihi

        color, mag = self._apply_extinction(raw_color, raw_mag, dm, e_gr)
        return color, mag, mass_grid

    def _interpolate_isochrone(
        self, log_age: float, mh: float, dm: float, e_gr: float
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (color, mag, mass) with extinction/distance applied.

        Uses bilinear interpolation between grid points when ``use_bilinear=True``
        (default); falls back to nearest-grid when corners are unavailable.
        """
        if self.use_bilinear:
            return self._interpolate_bilinear(log_age, mh, dm, e_gr)
        return self._interpolate_nearest(log_age, mh, dm, e_gr)

    def _fast_distance(
        self,
        obs_color: np.ndarray,
        obs_mag: np.ndarray,
        obs_color_err: np.ndarray,
        obs_mag_err: np.ndarray,
        iso_color: np.ndarray,
        iso_mag: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """KD-tree nearest-neighbour distance in error-normalised CMD space.

        Returns
        -------
        chi_dist : (N,) array of error-normalised distances per observed star
        iso_idx  : (N,) integer array of matched isochrone point indices
        """
        if len(iso_color) < 2:
            inf = np.full(len(obs_color), np.inf)
            return inf, np.zeros(len(obs_color), dtype=int)

        color_scale = max(np.median(obs_color_err), 1e-6)
        mag_scale   = max(np.median(obs_mag_err),   1e-6)

        iso_pts = np.column_stack([iso_color / color_scale, iso_mag / mag_scale])
        obs_pts = np.column_stack([obs_color / color_scale, obs_mag / mag_scale])

        tree = cKDTree(iso_pts)
        _, iso_idx = tree.query(obs_pts)

        closest_iso_c = iso_color[iso_idx]
        closest_iso_m = iso_mag[iso_idx]

        dist_c = (obs_color - closest_iso_c) / obs_color_err
        dist_m = (obs_mag   - closest_iso_m) / obs_mag_err
        chi_dist = np.sqrt(dist_c ** 2 + dist_m ** 2)

        return chi_dist, iso_idx

    def _imf_weight(self, mass: np.ndarray) -> np.ndarray:
        """Kroupa (2001) IMF probability density for each mass value.

        Returns normalised weights (sum = 1) proportional to dN/dM so that
        stars at densely-populated evolutionary stages (main sequence) receive
        higher weight in the chi² objective when ``use_imf_weight=True``.

        Stars below 0.08 M☉ (brown-dwarf regime) are assigned zero weight.
        """
        mass = np.asarray(mass, dtype=float)
        weights = np.zeros_like(mass)
        valid = mass >= 0.08
        m = mass[valid]
        w = np.where(m < 0.5,
                     m ** (-1.3),
                     0.5 ** (-1.3) * (m / 0.5) ** (-2.3))
        weights[valid] = w
        total = weights.sum()
        if total > 0:
            weights /= total
        return weights

    def _effective_dof(self, n_stars: int) -> int:
        """Degrees of freedom matching the chi2 sample size.

        ``_objective`` sums squared distances over only the closest
        ``fit_fraction`` of stars (``n_use``), so reduced_chi2 must divide by
        that same subset size minus the 4 fit parameters — not the full
        ``n_stars`` — otherwise reduced_chi2 is systematically under-estimated.
        """
        n_use = min(int(n_stars), max(20, int(n_stars * self.fit_fraction)))
        return max(1, n_use - 4)

    def _objective(self, params: np.ndarray) -> float:
        """
        Objective function using KD-tree distance with error normalization.
        """
        log_age, mh, dm, e_gr = params

        try:
            iso_c, iso_m, iso_mass = self._interpolate_isochrone(log_age, mh, dm, e_gr)

            if len(iso_c) < 10:
                return 1e10

            obs_c = self._fit_obs[:, 0]
            obs_m = self._fit_obs[:, 1]
            obs_c_err = self._fit_err[:, 0]
            obs_m_err = self._fit_err[:, 1]

            # KD-tree distance in error-normalised CMD space
            dist, iso_idx = self._fast_distance(
                obs_c, obs_m, obs_c_err, obs_m_err, iso_c, iso_m
            )

            # Robust: use only closest fraction
            n_use = max(20, int(len(dist) * self.fit_fraction))
            sorted_idx  = np.argsort(dist)
            closest_idx = sorted_idx[:n_use]

            # Base weights: membership (EM) if available, else uniform
            weights = self._fit_weight
            if weights is not None and len(weights) == len(dist):
                w = np.asarray(weights, dtype=float)[closest_idx]
                w = np.clip(w, 1e-3, None)
                w = w / w.mean()   # normalize so scale is comparable across EM rounds
            else:
                w = np.ones(len(closest_idx), dtype=float)

            # Optional IMF weighting: upweight stars near populous CMD stages.
            # Each observed star is weighted by the Kroupa IMF density at the
            # nearest matched isochrone mass — scaled so the mean weight stays 1.
            if self.use_imf_weight and len(iso_mass) > 0:
                matched_mass = iso_mass[iso_idx[closest_idx]]
                imf_w = self._imf_weight(matched_mass)
                # _imf_weight normalises by sum; rescale to mean=1 for chi² stability
                imf_w = imf_w * len(imf_w) / max(imf_w.sum(), 1e-30)
                w = w * imf_w

            chi2 = np.sum((dist[closest_idx] ** 2) * w) * (len(w) / max(w.sum(), 1e-30))
            return chi2

        except Exception:
            return 1e10

    def _de_callback(self, xk, convergence):
        self._fit_iteration += 1
        progress = min(0.95, self._fit_iteration / self._fit_max_iter)
        if self.progress_callback:
            self.progress_callback(progress, f"Iteration {self._fit_iteration}")

    def fit(
        self,
        obs_color: np.ndarray,
        obs_mag: np.ndarray,
        obs_color_err: np.ndarray,
        obs_mag_err: np.ndarray,
        mode: FitMode = FitMode.FAST,
        bounds: Optional[FitBounds] = None,
        snr_min: float = 5.0,
        **kwargs
    ) -> FitResult:
        """Fit isochrone to observed CMD"""
        t0 = time.time()

        if bounds is None:
            bounds = FitBounds()

        obs_color_err = np.clip(obs_color_err, 0.01, None)
        obs_mag_err = np.clip(obs_mag_err, 0.01, None)

        # Quality mask
        snr_c = 1.0 / obs_color_err
        snr_m = 1.0 / obs_mag_err
        snr = np.minimum(snr_c, snr_m)
        mask = (snr > snr_min) & np.isfinite(obs_color) & np.isfinite(obs_mag)

        color = obs_color[mask]
        mag = obs_mag[mask]
        color_err = obs_color_err[mask]
        mag_err = obs_mag_err[mask]

        n_stars = len(color)

        if n_stars < 20:
            return FitResult(
                log_age=9.0, metallicity=0.0, distance_mod=10.0, extinction_gr=0.0,
                chi2=np.inf, n_stars=n_stars, converged=False,
                fit_mode=mode.value, elapsed_sec=time.time() - t0
            )

        # Store for objective function
        self._fit_obs = np.column_stack([color, mag])
        self._fit_err = np.column_stack([color_err, mag_err])
        self._fit_weight = np.ones(n_stars, dtype=float)
        self._fit_iteration = 0

        em_iters = max(1, int(kwargs.pop("em_iters", 1)))
        membership_sigma_scale = float(kwargs.pop("membership_sigma_scale", 1.3))
        weight_floor = float(kwargs.pop("weight_floor", 0.05))

        best_result = None

        for em_idx in range(em_iters):
            # Scale per-round progress into overall [0, 1]
            original_callback = self.progress_callback
            if original_callback:
                def _round_callback(progress, message, _em=em_idx):
                    p = (float(_em) + float(np.clip(progress, 0.0, 1.0))) / float(em_iters)
                    original_callback(p, f"Round {_em + 1}/{em_iters}: {message}")
                self.progress_callback = _round_callback

            try:
                if mode == FitMode.FAST:
                    result = self._fit_fast(bounds, n_stars, **kwargs)
                elif mode == FitMode.HESSIAN:
                    result = self._fit_hessian(bounds, n_stars, **kwargs)
                else:
                    raise ValueError(f"Unsupported fit mode: {mode}")
            finally:
                self.progress_callback = original_callback

            if (best_result is None) or (float(result.chi2) < float(best_result.chi2)):
                best_result = result

            # Update weights from membership for next EM round.
            if em_idx < em_iters - 1:
                try:
                    prob = self.compute_membership(
                        result,
                        color,
                        mag,
                        sigma_scale=membership_sigma_scale,
                    )
                    if prob is not None and len(prob) == n_stars:
                        self._fit_weight = np.clip(np.asarray(prob, dtype=float), weight_floor, 1.0)
                except Exception:
                    pass

        result = best_result if best_result is not None else result

        result.n_stars = n_stars
        result.elapsed_sec = time.time() - t0
        result.fit_mode = mode.value

        self._fit_obs = None
        self._fit_err = None
        self._fit_weight = None

        return result

    def _fit_fast(self, bounds: FitBounds, n_stars: int, **kwargs) -> FitResult:
        maxiter = kwargs.get('maxiter', 100)
        seed = kwargs.get('seed', 42)
        self._fit_max_iter = maxiter

        if self.progress_callback:
            self.progress_callback(0.05, "Starting optimization...")

        result = differential_evolution(
            self._objective,
            bounds=bounds.to_list(),
            maxiter=maxiter,
            workers=1,
            updating='immediate',
            polish=True,
            seed=seed,
            tol=0.01,
            atol=0.01,
            popsize=10,
            mutation=(0.5, 1.0),
            recombination=0.7,
            callback=self._de_callback
        )

        if self.progress_callback:
            self.progress_callback(1.0, "Complete")

        chi2 = result.fun
        dof = self._effective_dof(n_stars)

        return FitResult(
            log_age=result.x[0],
            metallicity=result.x[1],
            distance_mod=result.x[2],
            extinction_gr=result.x[3],
            chi2=chi2,
            reduced_chi2=chi2 / dof,
            converged=result.success
        )

    def _fit_hessian(self, bounds: FitBounds, n_stars: int, **kwargs) -> FitResult:
        de_maxiter = int(kwargs.get("de_maxiter", kwargs.get("maxiter", 50)))
        local_maxiter = int(kwargs.get("local_maxiter", 100))
        n_starts = max(1, int(kwargs.get("n_starts", 6)))
        de_seed = int(kwargs.get("seed", 42))
        initial_guess = kwargs.get("initial_guess", None)
        bounds_list = bounds.to_list()

        def _clip_guess(guess):
            try:
                g = np.asarray(guess, dtype=float).reshape(-1)
            except Exception:
                return None
            if g.size != 4 or (not np.isfinite(g).all()):
                return None
            clipped = [float(np.clip(val, lo, hi)) for val, (lo, hi) in zip(g, bounds_list)]
            return np.array(clipped, dtype=float)

        original_callback = self.progress_callback

        # ---------------------------------------------------------------------
        # 1) Global multi-start search (Differential Evolution)
        # ---------------------------------------------------------------------
        global_candidates = []
        global_base = 0.05
        global_span = 0.50
        for i in range(n_starts):
            seg_base = global_base + global_span * (i / n_starts)
            seg_span = global_span / n_starts

            if original_callback:
                def _scaled_global_cb(progress, message, _b=seg_base, _s=seg_span, _i=i):
                    p = _b + _s * float(np.clip(progress, 0.0, 1.0))
                    original_callback(p, f"Global {_i + 1}/{n_starts}: {message}")
                self.progress_callback = _scaled_global_cb

            try:
                fast_result = self._fit_fast(
                    bounds,
                    n_stars,
                    maxiter=de_maxiter,
                    seed=de_seed + i * 1009,
                )
            finally:
                self.progress_callback = original_callback

            x_fast = np.array(
                [
                    fast_result.log_age,
                    fast_result.metallicity,
                    fast_result.distance_mod,
                    fast_result.extinction_gr,
                ],
                dtype=float,
            )
            global_candidates.append(
                {
                    "source": f"global#{i + 1}",
                    "x0": x_fast,
                    "chi2": float(fast_result.chi2),
                    "converged": bool(fast_result.converged),
                }
            )

        # ---------------------------------------------------------------------
        # 2) Build local-start candidate set
        # ---------------------------------------------------------------------
        local_candidates = []
        seen = set()

        g0 = _clip_guess(initial_guess)
        if g0 is not None:
            key = tuple(np.round(g0, 6))
            seen.add(key)
            local_candidates.append({"source": "slider", "x0": g0})
            if original_callback:
                original_callback(0.57, "Added slider initial guess candidate")

        for c in sorted(global_candidates, key=lambda x: x["chi2"]):
            key = tuple(np.round(c["x0"], 6))
            if key in seen:
                continue
            seen.add(key)
            local_candidates.append({"source": c["source"], "x0": c["x0"]})

        if not local_candidates and global_candidates:
            best_global = min(global_candidates, key=lambda x: x["chi2"])
            local_candidates.append({"source": best_global["source"], "x0": best_global["x0"]})

        # ---------------------------------------------------------------------
        # 3) Local refinement for each candidate; select minimum chi²
        # ---------------------------------------------------------------------
        best_local_result = None
        best_local_source = None
        local_base = 0.55
        local_span = 0.40
        n_local = max(1, len(local_candidates))

        for j, cand in enumerate(local_candidates):
            seg_base = local_base + local_span * (j / n_local)
            seg_span = local_span / n_local
            iter_box = {"n": 0}

            def _local_callback(_xk, _b=seg_base, _s=seg_span, _j=j, _src=cand["source"]):
                iter_box["n"] += 1
                if original_callback:
                    frac = min(1.0, iter_box["n"] / max(local_maxiter, 1))
                    p = _b + _s * frac
                    original_callback(
                        p,
                        f"Local {_j + 1}/{n_local} [{_src}] {iter_box['n']}/{local_maxiter}",
                    )

            if original_callback:
                original_callback(seg_base, f"Local refine {j + 1}/{n_local} [{cand['source']}]")

            try:
                result = minimize(
                    self._objective,
                    cand["x0"],
                    method='L-BFGS-B',
                    bounds=bounds_list,
                    options={'maxiter': local_maxiter},
                    callback=_local_callback,
                )
            except Exception:
                continue

            if (best_local_result is None) or (float(result.fun) < float(best_local_result.fun)):
                best_local_result = result
                best_local_source = cand["source"]

        # ---------------------------------------------------------------------
        # 4) Compare with best global-only result and finalize
        # ---------------------------------------------------------------------
        best_global = min(global_candidates, key=lambda x: x["chi2"]) if global_candidates else None
        use_local = best_local_result is not None
        use_global = best_global is not None

        errors = [None, None, None, None]
        converged = False

        if use_local and (not use_global or float(best_local_result.fun) <= float(best_global["chi2"])):
            best_x = np.asarray(best_local_result.x, dtype=float)
            chi2 = float(best_local_result.fun)
            converged = bool(best_local_result.success)
            if hasattr(best_local_result, 'hess_inv') and best_local_result.hess_inv is not None:
                try:
                    if hasattr(best_local_result.hess_inv, 'todense'):
                        hess_inv = np.array(best_local_result.hess_inv.todense())
                    else:
                        hess_inv = np.array(best_local_result.hess_inv)
                    # Objective is χ², so Hessian H ≈ 2·F (Fisher info) ⇒ Cov ≈ 2·H⁻¹.
                    # L-BFGS-B's hess_inv is also a BFGS approximation; treat as
                    # an order-of-magnitude estimate only.
                    errors = np.sqrt(np.abs(2.0 * np.diag(hess_inv)))
                except Exception:
                    pass
            selected_source = best_local_source or "local"
        elif use_global:
            best_x = np.asarray(best_global["x0"], dtype=float)
            chi2 = float(best_global["chi2"])
            converged = bool(best_global["converged"])
            selected_source = best_global["source"]
        else:
            # Fallback: center of bounds
            best_x = np.array([(lo + hi) * 0.5 for lo, hi in bounds_list], dtype=float)
            chi2 = float(self._objective(best_x))
            converged = False
            selected_source = "fallback"

        if original_callback:
            original_callback(0.98, f"Selected {selected_source} | chi2={chi2:.2f}")
            original_callback(1.0, "Complete")

        dof = self._effective_dof(n_stars)
        return FitResult(
            log_age=float(best_x[0]),
            metallicity=float(best_x[1]),
            distance_mod=float(best_x[2]),
            extinction_gr=float(best_x[3]),
            log_age_err=errors[0] if errors[0] else None,
            metallicity_err=errors[1] if errors[1] else None,
            distance_mod_err=errors[2] if errors[2] else None,
            extinction_gr_err=errors[3] if errors[3] else None,
            chi2=chi2,
            reduced_chi2=chi2 / dof,
            converged=converged,
        )

    # -----------------------------------------------------------------
    # Grid Scan
    # -----------------------------------------------------------------

    def _grid_objective_2p(self, params_2p: np.ndarray, age: float, mh: float) -> float:
        """2-parameter wrapper: optimize (DM, E(color)) for fixed (age, mh)."""
        dm, e_gr = params_2p
        return self._objective(np.array([age, mh, dm, e_gr]))

    def _fit_grid_scan(self, bounds: FitBounds, n_stars: int, **kwargs) -> GridScanResult:
        local_maxiter = int(kwargs.get("local_maxiter", 50))
        initial_dm = float(kwargs.get("initial_dm", (bounds.distance_mod[0] + bounds.distance_mod[1]) * 0.5))
        initial_egr = float(kwargs.get("initial_egr", (bounds.extinction_gr[0] + bounds.extinction_gr[1]) * 0.5))

        # Filter grid to bounds
        age_mask = (self.ages >= bounds.log_age[0]) & (self.ages <= bounds.log_age[1])
        mh_mask = (self.metallicities >= bounds.metallicity[0]) & (self.metallicities <= bounds.metallicity[1])
        grid_ages = self.ages[age_mask]
        grid_mhs = self.metallicities[mh_mask]

        n_ages = len(grid_ages)
        n_mhs = len(grid_mhs)
        total_cells = n_ages * n_mhs

        if total_cells == 0:
            empty_fit = FitResult(
                log_age=0.0, metallicity=0.0, distance_mod=0.0, extinction_gr=0.0,
                chi2=np.inf, n_stars=n_stars, converged=False, fit_mode="grid_scan",
            )
            return GridScanResult(
                best_fit=empty_fit,
                chi2_grid=np.full((0, 0), np.nan),
                dm_grid=np.full((0, 0), np.nan),
                egr_grid=np.full((0, 0), np.nan),
                grid_ages=grid_ages, grid_mhs=grid_mhs,
            )

        chi2_grid = np.full((n_ages, n_mhs), np.nan)
        dm_grid = np.full((n_ages, n_mhs), np.nan)
        egr_grid = np.full((n_ages, n_mhs), np.nan)

        bounds_2p = [bounds.distance_mod, bounds.extinction_gr]
        x0_2p = np.array([initial_dm, initial_egr], dtype=float)
        n_evaluated = 0

        for i, age in enumerate(grid_ages):
            for j, mh in enumerate(grid_mhs):
                cell_idx = i * n_mhs + j

                # Check cache
                key = (round(float(age), 2), round(float(mh), 2))
                if key not in self._iso_cache:
                    if self.progress_callback:
                        self.progress_callback(
                            (cell_idx + 1) / total_cells,
                            f"Cell {cell_idx + 1}/{total_cells} (skip)"
                        )
                    continue

                try:
                    res = minimize(
                        self._grid_objective_2p,
                        x0_2p,
                        args=(float(age), float(mh)),
                        method='L-BFGS-B',
                        bounds=bounds_2p,
                        options={'maxiter': local_maxiter},
                    )
                    chi2_grid[i, j] = float(res.fun)
                    dm_grid[i, j] = float(res.x[0])
                    egr_grid[i, j] = float(res.x[1])
                    n_evaluated += 1
                except Exception:
                    pass

                if self.progress_callback:
                    self.progress_callback(
                        (cell_idx + 1) / total_cells,
                        f"Cell {cell_idx + 1}/{total_cells} | age={age:.2f} [M/H]={mh:.2f}"
                    )

        # Find best cell
        if n_evaluated == 0:
            empty_fit = FitResult(
                log_age=0.0, metallicity=0.0, distance_mod=0.0, extinction_gr=0.0,
                chi2=np.inf, n_stars=n_stars, converged=False, fit_mode="grid_scan",
            )
            return GridScanResult(
                best_fit=empty_fit,
                chi2_grid=chi2_grid, dm_grid=dm_grid, egr_grid=egr_grid,
                grid_ages=grid_ages, grid_mhs=grid_mhs,
                n_evaluated=n_evaluated,
            )

        flat_idx = np.nanargmin(chi2_grid)
        i_best, j_best = np.unravel_index(flat_idx, chi2_grid.shape)
        best_age = float(grid_ages[i_best])
        best_mh = float(grid_mhs[j_best])
        best_dm = float(dm_grid[i_best, j_best])
        best_egr = float(egr_grid[i_best, j_best])
        best_chi2 = float(chi2_grid[i_best, j_best])

        # Final 4-param refinement for error estimation
        errors = [None, None, None, None]
        try:
            bounds_list = bounds.to_list()
            ref = minimize(
                self._objective,
                np.array([best_age, best_mh, best_dm, best_egr]),
                method='L-BFGS-B',
                bounds=bounds_list,
                options={'maxiter': 200},
            )
            if ref.fun <= best_chi2:
                best_age, best_mh, best_dm, best_egr = (float(v) for v in ref.x)
                best_chi2 = float(ref.fun)
            if hasattr(ref, 'hess_inv') and ref.hess_inv is not None:
                if hasattr(ref.hess_inv, 'todense'):
                    hess_inv = np.array(ref.hess_inv.todense())
                else:
                    hess_inv = np.array(ref.hess_inv)
                errors = list(np.sqrt(np.abs(2.0 * np.diag(hess_inv))))
        except Exception:
            pass

        dof = self._effective_dof(n_stars)
        best_fit = FitResult(
            log_age=best_age,
            metallicity=best_mh,
            distance_mod=best_dm,
            extinction_gr=best_egr,
            log_age_err=errors[0],
            metallicity_err=errors[1],
            distance_mod_err=errors[2],
            extinction_gr_err=errors[3],
            chi2=best_chi2,
            reduced_chi2=best_chi2 / dof,
            n_stars=n_stars,
            converged=True,
            fit_mode="grid_scan",
        )

        return GridScanResult(
            best_fit=best_fit,
            chi2_grid=chi2_grid, dm_grid=dm_grid, egr_grid=egr_grid,
            grid_ages=grid_ages, grid_mhs=grid_mhs,
            n_evaluated=n_evaluated,
        )

    def fit_grid_scan(
        self,
        obs_color: np.ndarray,
        obs_mag: np.ndarray,
        obs_color_err: np.ndarray,
        obs_mag_err: np.ndarray,
        bounds: Optional[FitBounds] = None,
        snr_min: float = 5.0,
        **kwargs,
    ) -> GridScanResult:
        """Exhaustive grid scan over all (age, mh) within bounds."""
        t0 = time.time()
        if bounds is None:
            bounds = FitBounds()

        obs_color_err = np.clip(obs_color_err, 0.01, None)
        obs_mag_err = np.clip(obs_mag_err, 0.01, None)
        snr_c = 1.0 / obs_color_err
        snr_m = 1.0 / obs_mag_err
        snr = np.minimum(snr_c, snr_m)
        mask = (snr > snr_min) & np.isfinite(obs_color) & np.isfinite(obs_mag)

        color = obs_color[mask]
        mag = obs_mag[mask]
        color_err = obs_color_err[mask]
        mag_err = obs_mag_err[mask]
        n_stars = len(color)

        if n_stars < 20:
            empty_fit = FitResult(
                log_age=0.0, metallicity=0.0, distance_mod=0.0, extinction_gr=0.0,
                chi2=np.inf, n_stars=n_stars, converged=False, fit_mode="grid_scan",
            )
            return GridScanResult(
                best_fit=empty_fit,
                chi2_grid=np.full((0, 0), np.nan),
                dm_grid=np.full((0, 0), np.nan),
                egr_grid=np.full((0, 0), np.nan),
                grid_ages=np.array([]), grid_mhs=np.array([]),
            )

        self._fit_obs = np.column_stack([color, mag])
        self._fit_err = np.column_stack([color_err, mag_err])
        self._fit_weight = np.ones(n_stars, dtype=float)

        try:
            grid_result = self._fit_grid_scan(bounds, n_stars, **kwargs)
        finally:
            self._fit_obs = None
            self._fit_err = None
            self._fit_weight = None

        grid_result.elapsed_sec = time.time() - t0
        grid_result.best_fit.n_stars = n_stars
        grid_result.best_fit.elapsed_sec = grid_result.elapsed_sec
        return grid_result

    def get_best_fit_isochrone(self, result: FitResult) -> Tuple[np.ndarray, np.ndarray]:
        """Get best-fit isochrone CMD coordinates"""
        iso_c, iso_m, _ = self._interpolate_isochrone(
            result.log_age, result.metallicity,
            result.distance_mod, result.extinction_gr
        )
        return iso_c, iso_m

    def compute_membership(
        self,
        result: FitResult,
        obs_color: np.ndarray,
        obs_mag: np.ndarray,
        sigma_scale: float = 1.0
    ) -> np.ndarray:
        """
        Compute membership probability for each star.

        Uses perpendicular distance to isochrone for consistency
        with the fitting algorithm.

        Parameters
        ----------
        result : FitResult
            Fitting result with best-fit parameters
        obs_color, obs_mag : arrays
            Observed CMD data
        sigma_scale : float
            Scale factor for membership threshold

        Returns
        -------
        prob : array
            Membership probability [0, 1] for each star
        """
        iso_c, iso_m, _ = self._interpolate_isochrone(
            result.log_age, result.metallicity,
            result.distance_mod, result.extinction_gr
        )

        if len(iso_c) < 10:
            return np.zeros(len(obs_color))

        # Use KD-tree for fast nearest neighbor
        iso_pts = np.column_stack([iso_c, iso_m])
        tree = cKDTree(iso_pts)

        obs_pts = np.column_stack([obs_color, obs_mag])
        dist, _ = tree.query(obs_pts)

        # Adaptive sigma based on median distance
        valid_dist = dist[np.isfinite(dist)]
        if len(valid_dist) == 0:
            return np.zeros(len(obs_color))

        sigma = np.median(valid_dist) * sigma_scale

        # Gaussian membership probability
        prob = np.exp(-0.5 * (dist / sigma) ** 2)
        prob = np.clip(prob, 0, 1)

        return prob
