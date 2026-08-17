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
from functools import lru_cache
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


# The third kind, the quiet one. It needs the source tree, so it is separate
# from the two above — but it is not a test-only concern: a window offering a
# knob that changes nothing is a defect a user meets, not a regression a suite
# catches. `specs_from_map` refuses to build a widget for anything this reports.
_CONFIG_LAYER = {
    "apex/config/parameter_map.py",
    "apex/config/parameters_cmd.py",
    "apex/config/parameters_lc.py",
    "apex/config/schema.py",
    "apex/config/config_audit.py",
    # Names a setting only to say that nothing reads it. Counting it as a reader
    # would make the check erase its own findings: writing the list down would
    # shorten the list.
    "tests/test_settings_nobody_reads.py",
}
_SKIP_DIRS = (".venv", ".venv-deploy", "build/", "dist/", "validation/", "benchmark/")


@lru_cache(maxsize=2)
def _python_sources(repo_root: Path) -> dict[str, str]:
    """Every source file that could read a setting.

    `validation/` is skipped on purpose and not because it does not matter: it
    is a junction onto the E: drive holding ~20 GB, and walking it turns this
    from a second into several minutes. Settings used only by a validation
    script are reported here as unread, which is the honest answer anyway —
    they do not reach the pipeline.
    """
    out: dict[str, str] = {}
    for path in repo_root.rglob("*.py"):
        rel = path.relative_to(repo_root).as_posix()
        if rel.startswith(_SKIP_DIRS) or "__pycache__" in rel:
            continue
        try:
            out[rel] = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
    return out


@lru_cache(maxsize=8)
def _unread_cached(root: Path, mode: str) -> tuple[tuple[str, ...], tuple[str, ...]]:
    import re

    from apex.config.parameter_map import toml_key_map_for_mode

    sources = _python_sources(root)
    attrs = sorted({row[1] for row in toml_key_map_for_mode(mode)
                    if len(row) >= 2 and isinstance(row[1], str)})
    dead, gui_only = [], []
    for attr in attrs:
        pattern = re.compile(rf"\b{re.escape(attr)}\b")
        users = [rel for rel, text in sources.items()
                 if rel not in _CONFIG_LAYER and pattern.search(text)]
        if not users:
            dead.append(attr)
        elif all(u.startswith("apex/gui/") for u in users):
            gui_only.append(attr)
    return tuple(dead), tuple(gui_only)


def unread_settings(repo_root: str | Path | None = None,
                    mode: str = "cmd") -> dict[str, list[str]]:
    """Settings that reach `P` and that no module outside the config layer reads.

    Returns {"dead": [...], "gui_only": [...]}. The split matters: a setting the
    desktop reads is doing something even if the pipeline never sees it, but it
    also means a headless run silently ignores it.

    A name counts as read if it appears anywhere in a source file, not only in a
    `getattr`. That is deliberately generous — two call sites build the attribute
    name at run time from a literal list, and a stricter rule would call those
    dead. Being generous here means the list under-reports, never over-reports.
    """
    root = Path(repo_root) if repo_root else Path(__file__).absolute().parents[2]
    dead, gui_only = _unread_cached(root, mode)
    return {"dead": list(dead), "gui_only": list(gui_only)}


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("config", type=Path, nargs="?",
                    help="검사할 apex_config.json (--unread 만 볼 때는 생략 가능)")
    ap.add_argument("--mode", choices=("cmd", "lc"), default="cmd")
    ap.add_argument("--unread", action="store_true",
                    help="세 번째 종류도 검사한다 — P 에 도착하지만 아무도 안 읽는 설정")
    args = ap.parse_args(argv)

    def section(heading: str, names: list[str]) -> None:
        print()
        print(f"=== {heading} ({len(names)}) ===")
        for name in names:
            print(f"  {name}")

    if args.config is not None:
        result = audit(args.config, args.mode)
        section("매핑이 없다 — 로더가 통째로 무시한다", result["unmapped"])
        section("매핑은 됐지만 생성자가 안 받아 사라진다", result["dropped"])
        total = len(result["unmapped"]) + len(result["dropped"])
        print()
        print(f"적었지만 코드에 닿지 않는 설정 {total} 개")

    if args.unread or args.config is None:
        found = unread_settings(mode=args.mode)
        section("도착하지만 아무도 안 읽는다", found["dead"])
        section("GUI 만 읽는다 — 헤드리스는 무시한다", found["gui_only"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
