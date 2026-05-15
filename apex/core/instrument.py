"""
Instrument configuration (Telescope + Camera).
"""

from __future__ import annotations
import math
import re
import numpy as np
import pandas as pd
from pathlib import Path
from astropy.coordinates import SkyCoord
import astropy.units as u
from typing import Optional

from apex.utils.step_paths import step1_dir
from apex.utils.ssl_certificates import configure_ssl_certificates

configure_ssl_certificates()

try:
    from astroquery.simbad import Simbad
    _HAS_SIMBAD = True
except Exception as e:
    print(f"[SIMBAD] astroquery.simbad import failed → SIMBAD unavailable: {e}")
    _HAS_SIMBAD = False


_CATALOG_ALIASES = {
    "M3": ["NGC 5272"],
    "MESSIER3": ["NGC 5272"],
    "M5": ["NGC 5904"],
    "MESSIER5": ["NGC 5904"],
}

_OFFLINE_TARGETS = {
    "M3": {
        "ra_hms": "13:42:11.62",
        "dec_dms": "+28:22:38.2",
        "otype": "GlCl",
        "canonical": "M 3",
    },
    "MESSIER3": {
        "ra_hms": "13:42:11.62",
        "dec_dms": "+28:22:38.2",
        "otype": "GlCl",
        "canonical": "M 3",
    },
    "NGC5272": {
        "ra_hms": "13:42:11.62",
        "dec_dms": "+28:22:38.2",
        "otype": "GlCl",
        "canonical": "NGC 5272",
    },
    "M5": {
        "ra_hms": "15:18:33.22",
        "dec_dms": "+02:04:51.7",
        "otype": "GlCl",
        "canonical": "M 5",
    },
    "MESSIER5": {
        "ra_hms": "15:18:33.22",
        "dec_dms": "+02:04:51.7",
        "otype": "GlCl",
        "canonical": "M 5",
    },
    "NGC5904": {
        "ra_hms": "15:18:33.22",
        "dec_dms": "+02:04:51.7",
        "otype": "GlCl",
        "canonical": "NGC 5904",
    },
}


class InstrumentConfig:
    """
    Telescope and camera configuration.
    Handles pixel scale, FOV calculation, and target resolution.
    """

    def __init__(self, params, binning: Optional[int] = None):
        """
        Args:
            params: Parameters object from config module.
            binning: Binning mode override (default: from params).
        """
        self.params = params

        # Telescope specs — defaults for Planewave CDK500 @ KNUEMAO.
        # Override via params.P if the hardware changes.
        self.telescope_name = getattr(params.P, "telescope_name", "Planewave CDK500")
        self.aperture_mm = float(getattr(params.P, "telescope_aperture_mm", 508.0))
        self.focal_length_mm = float(getattr(params.P, "telescope_focal_mm", 3947.0))
        self.focal_ratio = self.focal_length_mm / self.aperture_mm

        # Camera specs — defaults for Moravian C3-61000.
        self.camera_name = getattr(params.P, "camera_name", "Moravian C3-61000")
        self.pix_size_um = float(getattr(params.P, "camera_pixel_um", 3.76))
        self.sensor_w_mm = float(getattr(params.P, "sensor_w_mm", 36.01))
        self.sensor_h_mm = float(getattr(params.P, "sensor_h_mm", 24.02))
        self.sensor_nx_1x = int(getattr(params.P, "sensor_nx_1x", 9576))
        self.sensor_ny_1x = int(getattr(params.P, "sensor_ny_1x", 6388))

        # Binning
        if binning is not None:
            self.binning = int(binning)
        else:
            try:
                default_bin = getattr(params.P, "binning_default", None)
                if default_bin is None:
                    default_bin = getattr(params.P, "camera_binning", 2)
                self.binning = int(float(default_bin or 2))
            except Exception:
                self.binning = 2

        self.pix_scale_1x = self._pixel_scale_arcsec(1)
        self.pix_scale_bin = self._pixel_scale_arcsec(self.binning)
        self.fov_w_deg = self._fov_deg(self.sensor_w_mm)
        self.fov_h_deg = self._fov_deg(self.sensor_h_mm)

        self.params.P.pixel_scale_arcsec = float(self.pix_scale_bin)
        self.params.P.telescope_focal_mm = float(self.focal_length_mm)
        self.params.P.camera_pixel_um = float(self.pix_size_um)
        self.params.P.binning_default = int(self.binning)

        self._apply_fwhm_conversions()

        self.targets_resolved = []
        self.primary_target = None
        self.primary_coord = None
        self.last_target_attempts = []
        self.last_target_errors = []

    def _pixel_scale_arcsec(self, binning: int = 1) -> float:
        return 206.265 * self.pix_size_um * float(binning) / float(self.focal_length_mm)

    def _fov_deg(self, sensor_mm: float) -> float:
        return 57.2957795 * sensor_mm / float(self.focal_length_mm)

    def _apply_fwhm_conversions(self):
        P = self.params.P
        arc = getattr(P, "fwhm_guess_arcsec", None)
        if arc is not None and np.isfinite(arc) and arc > 0:
            P.fwhm_seed_px = max(2.0, float(arc) / self.pix_scale_bin)
            P._fwhm_seed_from = "arcsec"

        arcmin = getattr(P, "fwhm_arcsec_min", None)
        if arcmin is not None and np.isfinite(arcmin):
            P.fwhm_px_min = max(float(P.fwhm_px_min), float(arcmin) / self.pix_scale_bin)

        arcmax = getattr(P, "fwhm_arcsec_max", None)
        if arcmax is not None and np.isfinite(arcmax):
            P.fwhm_px_max = min(float(P.fwhm_px_max), float(arcmax) / self.pix_scale_bin)

    @staticmethod
    def _target_query_aliases(name: str) -> list[str]:
        raw = " ".join(str(name or "").strip().split())
        aliases: list[str] = []

        def add(value: str) -> None:
            value = " ".join(str(value or "").strip().split())
            if value and value not in aliases:
                aliases.append(value)

        add(raw)
        compact = raw.replace(" ", "")

        messier = re.fullmatch(r"(?i)(?:m|messier)0*([0-9]+[a-z]?)", compact)
        if messier:
            number = messier.group(1).upper()
            add(f"M {number}")
            add(f"Messier {number}")
            for alias in _CATALOG_ALIASES.get(f"M{number}", []):
                add(alias)

        catalog = re.fullmatch(r"(?i)(ngc|ic)0*([0-9]+[a-z]?)", compact)
        if catalog:
            prefix = catalog.group(1).upper()
            number = catalog.group(2).upper()
            add(f"{prefix} {number}")

        return aliases

    @classmethod
    def _offline_target_record(cls, name: str) -> Optional[dict]:
        for alias in cls._target_query_aliases(name):
            key = alias.replace(" ", "").upper()
            rec = _OFFLINE_TARGETS.get(key)
            if not rec:
                continue
            coord = SkyCoord(rec["ra_hms"], rec["dec_dms"], unit=(u.hourangle, u.deg))
            return {
                "name": name,
                "ra_deg": float(coord.ra.deg),
                "dec_deg": float(coord.dec.deg),
                "ra_str": coord.ra.to_string(unit="hour", sep=":", precision=2),
                "dec_str": coord.dec.to_string(unit="deg", sep=":", precision=1, alwayssign=True),
                "vmag": np.nan,
                "otype": rec.get("otype", ""),
                "simbad_query": f"offline:{rec.get('canonical', alias)}",
            }
        return None

    def resolve_targets(self, target_names: Optional[list] = None):
        raw_targets = list(target_names or [])

        if hasattr(self.params.P, "_raw"):
            tn_param = self.params.P._raw.get("target_name", None)
            if tn_param and str(tn_param).strip():
                raw_targets.append(str(tn_param).strip())

        target_list_path = self.params.P.data_dir / "targets.txt"
        if target_list_path.exists():
            for line in target_list_path.read_text(encoding="utf-8").splitlines():
                s = line.strip()
                if s and not s.startswith("#"):
                    raw_targets.append(s)

        if not raw_targets:
            raw_targets = ["M31"]

        seen: set = set()
        targets = []
        for name in raw_targets:
            if name not in seen:
                seen.add(name)
                targets.append(name)

        simbad_ready = bool(_HAS_SIMBAD)
        self.last_target_errors = []
        if simbad_ready:
            try:
                Simbad.reset_votable_fields()
                Simbad.add_votable_fields("ra(d)", "dec(d)", "flux(V)", "otype")
            except Exception as e:
                simbad_ready = False
                self.last_target_errors.append(f"SIMBAD setup: {type(e).__name__}: {e}")
        else:
            self.last_target_errors.append("astroquery.simbad is not importable")

        print("\n=== SIMBAD Target Resolution ===")
        self.last_target_attempts = []
        for name in targets:
            res = None
            query_used = None
            aliases = self._target_query_aliases(name)
            self.last_target_attempts.extend(aliases)
            if simbad_ready:
                for query_name in aliases:
                    try:
                        res = Simbad.query_object(query_name)
                        if res is not None and len(res) > 0:
                            query_used = query_name
                            break
                    except Exception as e:
                        self.last_target_errors.append(
                            f"{query_name}: {type(e).__name__}: {e}"
                        )

            if res is None or len(res) == 0:
                offline = self._offline_target_record(name)
                if offline is not None:
                    self.targets_resolved.append(offline)
                    print(
                        f"[SIMBAD] {name:20s} ({offline['simbad_query']}) -> "
                        f"RA={offline['ra_str']}  DEC={offline['dec_str']}  "
                        f"({offline['ra_deg']:9.5f} deg, {offline['dec_deg']:9.5f} deg)  "
                        f"V~n/a  [{offline['otype']}]"
                    )
                    continue

                err_txt = "; ".join(self.last_target_errors[-3:])
                suffix = f" errors: {err_txt}" if err_txt else ""
                print(f"[SIMBAD] Not found: {name} (tried: {', '.join(aliases)}){suffix}")
                continue

            try:
                ra_deg = float(res["RA_d"][0])
                dec_deg = float(res["DEC_d"][0])
                ra_str = str(res["RA"][0])
                dec_str = str(res["DEC"][0])
                vmag = float(res["FLUX_V"][0]) if "FLUX_V" in res.colnames else np.nan
                otype = str(res["OTYPE"][0]) if "OTYPE" in res.colnames else ""

                rec = dict(name=name, ra_deg=ra_deg, dec_deg=dec_deg,
                           ra_str=ra_str, dec_str=dec_str, vmag=vmag, otype=otype,
                           simbad_query=query_used or name)
                self.targets_resolved.append(rec)
                print(
                    f"[SIMBAD] {name:20s} ({query_used or name}) → RA={ra_str}  DEC={dec_str}  "
                    f"({ra_deg:9.5f}°, {dec_deg:9.5f}°)  V~{vmag if np.isfinite(vmag) else 'n/a'}  [{otype}]"
                )
            except Exception as e:
                self.last_target_errors.append(f"{name}: {type(e).__name__}: {e}")
                print(f"[SIMBAD] ERROR for '{name}': {e}")

        if self.targets_resolved:
            self.primary_target = self.targets_resolved[0]["name"]
            self.primary_coord = SkyCoord(
                self.targets_resolved[0]["ra_deg"],
                self.targets_resolved[0]["dec_deg"],
                unit="deg",
            )
            step1_out = step1_dir(self.params.P.result_dir)
            step1_out.mkdir(parents=True, exist_ok=True)
            out = step1_out / "targets_simbad.tsv"
            pd.DataFrame(self.targets_resolved).to_csv(out, sep="\t", index=False)
            print(f"→ SIMBAD results saved: {out}")
        else:
            print("No targets successfully resolved via SIMBAD")

    def print_summary(self):
        print(f"\n=== Instrument Summary ({self.telescope_name}) ===")
        print(f"Telescope : {self.telescope_name}")
        print(f"  Aperture D = {self.aperture_mm:.1f} mm")
        print(f"  Focal length F = {self.focal_length_mm:.1f} mm  (f/{self.focal_ratio:.2f})")
        print(f"Camera    : {self.camera_name}")
        print(f"  Pixel size (1x1) = {self.pix_size_um:.3f} μm")
        print(
            f"  Sensor = {self.sensor_w_mm:.2f} × {self.sensor_h_mm:.2f} mm  "
            f"({self.sensor_nx_1x} × {self.sensor_ny_1x} px @1x1)"
        )
        print(f"Default binning = {self.binning}×{self.binning}")
        print(f"Pixel scale (1x1)  ≈ {self.pix_scale_1x:.6f} \" / px")
        print(
            f"Pixel scale ({self.binning}x{self.binning}) ≈ {self.pix_scale_bin:.6f} \" / px  "
            f"(pipeline pixel_scale_arcsec={self.params.P.pixel_scale_arcsec:.6f} \" / px)"
        )
        print(f"Full FOV ≈ {self.fov_w_deg:.3f} × {self.fov_h_deg:.3f} deg")

        if self.primary_target is not None and self.primary_coord is not None:
            ra_hms = self.primary_coord.ra.to_string(unit=u.hour, sep=":", precision=2)
            dec_dms = self.primary_coord.dec.to_string(unit=u.deg, sep=":", precision=1, alwayssign=True)
            print(f"\nPRIMARY TARGET: {self.primary_target}  (RA={ra_hms}, DEC={dec_dms})")
