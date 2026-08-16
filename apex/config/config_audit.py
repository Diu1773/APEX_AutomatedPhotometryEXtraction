"""Which settings in a config file actually reach the code.

A config value can die in three places, and until 2026-08-16 all three were
silent — the run succeeded, the log said nothing, and the setting did nothing:

    unmapped   the key is not in the key map, so it never enters `raw`
    dropped    it enters `raw`, but the namespace constructor does not name it
               as a keyword argument, so it is discarded on the way into `P`
    unread     it reaches `P` and no module ever reads that attribute

The third one needs the source tree, so it lives in the test. The first two are
computable from a loaded `Parameters` and are what this module reports.

The case that forced this: `photometry.apcorr.small_scale` was mapped, entered
`raw`, and had no keyword in `parameters_cmd`'s `SimpleNamespace(...)` call, so
`P.apcorr_small_scale` did not exist at all. Step 7 meanwhile read the aperture
scale off `forced_r_ap_scale`, a name nothing sets, and took its literal. Both
halves worked; they were just never connected. Fifty config files on disk
carried an aperture radius that no run had ever used.

    python -X utf8 -m apex.config.config_audit <apex_config.json>
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable

# Keys that carry structure rather than settings, or that the loader consumes
# into a differently-named attribute on purpose.
_STRUCTURAL = ("_meta", "schema_version", "_calibration")


def _flatten(node: dict, known: set[str], prefix: tuple[str, ...] = ()) -> Iterable[str]:
    """Leaf keys, but stop at any subtree the map claims whole.

    `detection.sigma_by_filter` is mapped as one dict, so descending into it
    would report `detection.sigma_by_filter.B` as unmapped when it is not.
    """
    for key, value in node.items():
        here = prefix + (key,)
        dotted = ".".join(here)
        if dotted in known:
            continue
        if isinstance(value, dict):
            yield from _flatten(value, known, here)
        else:
            yield dotted


def dropped_settings(P: Any) -> list[str]:
    """Flat names the file set that never became an attribute of `P`.

    Dotted names are excluded: the `5x.*` HUD settings are deliberately
    raw-only and are read back through `P._raw[...]`, never as attributes.
    """
    raw = getattr(P, "_raw", None) or {}
    return sorted(
        name for name in raw
        if not name.startswith("_") and "." not in name and not hasattr(P, name)
    )


def unmapped_keys(config_path: str | Path, key_map: Iterable[tuple]) -> list[str]:
    """Dotted keys present in the file that the key map does not know."""
    data = json.loads(Path(config_path).read_text(encoding="utf-8"))
    known = {".".join(row[0]) for row in key_map}
    return sorted(
        dotted for dotted in _flatten(data, known)
        if not dotted.startswith(_STRUCTURAL)
    )


def audit(config_path: str | Path, mode: str = "cmd") -> dict[str, list[str]]:
    """Both silent losses for one config file, in one call."""
    if mode == "lc":
        from apex.config.parameters_lc import TOML_KEY_MAP, read_params
    else:
        from apex.config.parameters_cmd import TOML_KEY_MAP, read_params
    P = read_params(config_path).P
    return {
        "unmapped": unmapped_keys(config_path, TOML_KEY_MAP),
        "dropped": dropped_settings(P),
    }


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--mode", choices=("cmd", "lc"), default="cmd")
    args = ap.parse_args(argv)

    result = audit(args.config, args.mode)
    for label, heading in (("unmapped", "매핑이 없다 — 로더가 통째로 무시한다"),
                           ("dropped", "매핑은 됐지만 생성자가 안 받아 사라진다")):
        names = result[label]
        print(f"\n=== {heading} ({len(names)}) ===")
        for name in names:
            print(f"  {name}")
    total = len(result["unmapped"]) + len(result["dropped"])
    print(f"\n적었지만 코드에 닿지 않는 설정 {total} 개")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
