"""Background target-resolution worker for Step 1 SIMBAD lookups."""
from __future__ import annotations

from PyQt5.QtCore import QThread, pyqtSignal


class TargetResolveWorker(QThread):
    """Resolve a target name without blocking the Qt GUI thread."""

    log_message = pyqtSignal(str)
    result_ready = pyqtSignal(dict)
    error_message = pyqtSignal(str)

    def __init__(self, instrument, target_name: str, parent=None):
        super().__init__(parent)
        self.instrument = instrument
        self.target_name = str(target_name or "").strip()

    def run(self) -> None:
        inst = self.instrument
        try:
            inst.targets_resolved = []
            inst.primary_target = None
            inst.primary_coord = None
            inst.resolve_targets([self.target_name], log_fn=self.log_message.emit)

            tried = list(getattr(inst, "last_target_attempts", []) or [])
            errors = list(getattr(inst, "last_target_errors", []) or [])
            coord = getattr(inst, "primary_coord", None)
            if coord is None:
                self.result_ready.emit({
                    "ok": False,
                    "name": self.target_name,
                    "tried": tried,
                    "errors": errors,
                })
                return

            ra_deg = float(coord.ra.deg)
            dec_deg = float(coord.dec.deg)
            self.result_ready.emit({
                "ok": True,
                "name": self.target_name,
                "ra_deg": ra_deg,
                "dec_deg": dec_deg,
                "ra_hms": coord.ra.to_string(unit="hour", sep=":", precision=2),
                "dec_dms": coord.dec.to_string(
                    unit="deg", sep=":", precision=1, alwayssign=True
                ),
                "tried": tried,
                "errors": errors,
            })
        except Exception as exc:
            self.error_message.emit(f"{type(exc).__name__}: {exc}")


def target_failure_message(result: dict) -> str:
    """Build a compact user-facing message for failed target resolution."""
    name = result.get("name", "")
    tried = result.get("tried") or []
    errors = result.get("errors") or []
    parts = [f"Target not resolved: {name}"]
    if tried:
        parts.append(f"Tried: {', '.join(str(x) for x in tried[:8])}")
    if errors:
        parts.append(f"Last error: {errors[-1]}")
    parts.append("Use Manual RA/Dec, or open SIMBAD in a browser and retry with a catalog alias.")
    return "\n".join(parts)
