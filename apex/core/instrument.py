"""
Instrument configuration (Telescope + Camera).
"""

from __future__ import annotations
import math
import re
import warnings
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
    from astroquery.simbad import Simbad, conf as _SIMBAD_CONF
    _HAS_SIMBAD = True
except Exception as e:
    print(f"[SIMBAD] astroquery.simbad import failed → SIMBAD unavailable: {e}")
    _SIMBAD_CONF = None
    _HAS_SIMBAD = False


_CATALOG_ALIASES = {
    "M3":       ["NGC 5272"],
    "MESSIER3": ["NGC 5272"],
    "M5":       ["NGC 5904"],
    "MESSIER5": ["NGC 5904"],
    "M13":      ["NGC 6205"],
    "MESSIER13":["NGC 6205"],
    "M92":      ["NGC 6341"],
    "MESSIER92":["NGC 6341"],
    "M31":      ["NGC 224"],
    "MESSIER31":["NGC 224"],
    "M37":      ["NGC 2099"],
    "MESSIER37":["NGC 2099"],
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
        self.last_simbad_servers = []

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

    @staticmethod
    def _table_value(table, *candidates):
        names = {str(name).lower(): name for name in getattr(table, "colnames", [])}
        for candidate in candidates:
            name = names.get(str(candidate).lower())
            if name is not None:
                return table[name][0]
        raise KeyError(f"none of {candidates} in {getattr(table, 'colnames', [])}")

    @classmethod
    def _record_from_simbad_row(cls, name: str, res, query_used: str | None) -> dict:
        try:
            ra_deg = float(cls._table_value(res, "RA_d", "RA(d)", "ra_d"))
            dec_deg = float(cls._table_value(res, "DEC_d", "DEC(d)", "dec_d"))
            coord = SkyCoord(ra_deg, dec_deg, unit="deg")
        except Exception:
            ra_raw = cls._table_value(res, "ra", "RA")
            dec_raw = cls._table_value(res, "dec", "DEC")
            try:
                ra_deg = float(ra_raw)
                dec_deg = float(dec_raw)
                if np.isfinite(ra_deg) and np.isfinite(dec_deg):
                    coord = SkyCoord(ra_deg, dec_deg, unit="deg")
                else:
                    raise ValueError("non-finite numeric coordinates")
            except Exception:
                coord = SkyCoord(str(ra_raw), str(dec_raw), unit=(u.hourangle, u.deg))

        try:
            vmag_raw = cls._table_value(res, "FLUX_V", "flux_V", "V")
            vmag = np.nan if np.ma.is_masked(vmag_raw) else float(vmag_raw)
        except Exception:
            vmag = np.nan
        try:
            otype = str(cls._table_value(res, "OTYPE", "otype"))
        except Exception:
            otype = ""

        return dict(
            name=name,
            ra_deg=float(coord.ra.deg),
            dec_deg=float(coord.dec.deg),
            ra_str=coord.ra.to_string(unit="hour", sep=":", precision=2),
            dec_str=coord.dec.to_string(unit="deg", sep=":", precision=1, alwayssign=True),
            vmag=vmag,
            otype=otype,
            simbad_query=query_used or name,
        )

    def _prepare_simbad(self) -> bool:
        """Prepare optional SIMBAD fields, but keep default query usable on failure."""
        if not _HAS_SIMBAD:
            self.last_target_errors.append("astroquery.simbad is not importable")
            return False
        try:
            Simbad.clear_cache()
        except Exception:
            pass
        try:
            Simbad.reset_votable_fields()
        except Exception as e:
            self.last_target_errors.append(f"SIMBAD reset: {type(e).__name__}: {e}")
        try:
            Simbad.add_votable_fields("otype")
        except Exception as e:
            self.last_target_errors.append(
                f"SIMBAD optional field otype: {type(e).__name__}: {e}"
            )
        try:
            timeout_s = float(getattr(self.params.P, "simbad_timeout_s", 60.0))
            Simbad.TIMEOUT = max(timeout_s, 10.0)
        except Exception:
            pass
        return True

    def _simbad_server_candidates(self) -> list[str]:
        configured = getattr(self.params.P, "simbad_servers", None)
        if configured is None:
            configured = getattr(self.params.P, "simbad_server", None)

        servers: list[str] = []

        def add(value) -> None:
            text = str(value or "").strip().strip("/")
            if text and text not in servers:
                servers.append(text)

        if isinstance(configured, (list, tuple)):
            for item in configured:
                add(item)
        elif configured:
            for item in str(configured).split(","):
                add(item)

        try:
            for item in list(getattr(_SIMBAD_CONF, "servers_list", []) or []):
                add(item)
        except Exception:
            pass

        add("simbad.cds.unistra.fr")
        add("simbad.harvard.edu")
        return servers

    @staticmethod
    def _set_simbad_server(server: str) -> str:
        value = str(server or "").strip().rstrip("/")
        if not value:
            value = "simbad.cds.unistra.fr"
        if value.startswith(("http://", "https://")):
            if value.endswith("/sim-script"):
                url = value
            elif value.endswith("/simbad"):
                url = f"{value}/sim-script"
            else:
                url = f"{value}/simbad/sim-script"
            label = value.split("://", 1)[-1].split("/", 1)[0]
        else:
            label = value
            url = f"https://{value}/simbad/sim-script"
        Simbad.SIMBAD_URL = url
        return label

    def resolve_targets(self, target_names: Optional[list] = None, log_fn=None):
        raw_targets = [
            str(name).strip()
            for name in (target_names or [])
            if str(name).strip()
        ]

        seen: set = set()
        targets = []
        for name in raw_targets:
            if name not in seen:
                seen.add(name)
                targets.append(name)

        _log = log_fn or print

        self.targets_resolved = []
        self.primary_target = None
        self.primary_coord = None
        self.last_target_errors = []
        simbad_ready = self._prepare_simbad()
        simbad_servers = self._simbad_server_candidates() if simbad_ready else []
        self.last_simbad_servers = simbad_servers

        _log("\n=== SIMBAD Target Resolution ===")
        try:
            stale_out = step1_dir(Path(self.params.P.result_dir)) / "targets_simbad.tsv"
            stale_out.unlink(missing_ok=True)
        except Exception:
            pass
        self.last_target_attempts = []
        for name in targets:
            res = None
            query_used = None
            query_server = None
            query_attempted = False
            query_failed = False
            aliases = self._target_query_aliases(name)
            self.last_target_attempts.extend(aliases)
            if simbad_ready:
                for server in simbad_servers:
                    server_label = self._set_simbad_server(server)
                    for query_name in aliases:
                        try:
                            query_attempted = True
                            with warnings.catch_warnings():
                                warnings.simplefilter("ignore")
                                res = Simbad.query_object(query_name)
                            if res is not None and len(res) > 0:
                                query_used = query_name
                                query_server = server_label
                                break
                        except Exception as e:
                            query_failed = True
                            self.last_target_errors.append(
                                f"{query_name}@{server_label}: {type(e).__name__}: {e}"
                            )
                    if res is not None and len(res) > 0:
                        break

            if res is None or len(res) == 0:
                err_txt = "; ".join(self.last_target_errors[-3:])
                suffix = f"  ({err_txt})" if err_txt else ""
                server_txt = f"; servers: {', '.join(simbad_servers)}" if simbad_servers else ""
                if not simbad_ready:
                    _log(f"[SIMBAD] unavailable: {name} (tried: {', '.join(aliases)}{server_txt}){suffix}")
                elif query_failed and not query_attempted:
                    _log(f"[SIMBAD] query failed: {name} (tried: {', '.join(aliases)}{server_txt}){suffix}")
                elif query_failed and err_txt:
                    _log(f"[SIMBAD] query failed: {name} (tried: {', '.join(aliases)}{server_txt}){suffix}")
                else:
                    _log(f"[SIMBAD] Not found: {name} (tried: {', '.join(aliases)}{server_txt}){suffix}")
                continue

            try:
                rec = self._record_from_simbad_row(name, res, query_used)
                rec["simbad_server"] = query_server or ""
                self.targets_resolved.append(rec)
                source_label = query_used or name
                if query_server:
                    source_label = f"{source_label} @ {query_server}"
                _log(
                    f"[SIMBAD] {name:20s} ({source_label}) -> "
                    f"RA={rec['ra_str']}  DEC={rec['dec_str']}  "
                    f"({rec['ra_deg']:9.5f} deg, {rec['dec_deg']:9.5f} deg)  "
                    f"V~{rec['vmag'] if np.isfinite(rec['vmag']) else 'n/a'}  [{rec['otype']}]"
                )
            except Exception as e:
                self.last_target_errors.append(f"{name}: {type(e).__name__}: {e}")
                _log(f"[SIMBAD] parse error for '{name}': {e}")

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
            _log(f"[SIMBAD] results saved: {out}")
        else:
            _log("No targets successfully resolved via SIMBAD")

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
