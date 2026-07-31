"""Read the ``[calibration]`` TOML table into a flat settings dict.

Shared by the CMD and LC parameter models. Until this existed the section was
documented in ``parameters.example.toml`` but nothing parsed it, so editing it
had no effect — the GUI always started from the dataclass defaults.

The keys produced here are exactly ``CalibrationOptions`` field names, so the
analysis layer can consume them with ``CalibrationOptions.from_mapping`` and
never has to know about the config models.
"""

from __future__ import annotations

from typing import Any, Dict, Mapping, Optional

# ``[calibration.overscan]`` is a sub-table in the TOML but flat fields on
# CalibrationOptions, so its keys are prefixed on the way in.
_OVERSCAN_KEYS = {
    "enable": "overscan_enable",
    "edge": "overscan_edge",
    "width": "overscan_width",
    "trim": "overscan_trim",
}


def read_calibration_section(data: Optional[Mapping[str, Any]]) -> Dict[str, Any]:
    """Flatten ``[calibration]`` (+ ``[calibration.overscan]``) from parsed TOML.

    Unknown keys are passed through unchanged; ``CalibrationOptions`` drops the
    ones it does not recognise. Returns ``{}`` when the section is absent.
    """
    if not isinstance(data, Mapping):
        return {}
    section = data.get("calibration")
    if not isinstance(section, Mapping):
        return {}

    out: Dict[str, Any] = {}
    for key, value in section.items():
        if key == "overscan":
            continue
        out[str(key)] = value

    overscan = section.get("overscan")
    if isinstance(overscan, Mapping):
        for key, value in overscan.items():
            out[_OVERSCAN_KEYS.get(str(key), f"overscan_{key}")] = value
    return out


def calibration_toml_sections(settings: Mapping[str, Any]) -> Dict[tuple, Dict[str, Any]]:
    """Split a flat settings dict back into TOML tables for writing.

    Inverse of :func:`read_calibration_section`; ``enabled`` is preserved by the
    writer because it is not a ``CalibrationOptions`` field.
    """
    inverse = {v: k for k, v in _OVERSCAN_KEYS.items()}
    main: Dict[str, Any] = {}
    overscan: Dict[str, Any] = {}
    for key, value in settings.items():
        if key in inverse:
            overscan[inverse[key]] = value
        else:
            main[str(key)] = value
    return {("calibration",): main, ("calibration", "overscan"): overscan}
