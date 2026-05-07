"""
Inventory APEX parameter and cache foundations.

This script is intentionally read-only unless --write-docs is passed. It scans
the current TOML files, runtime parameter maps, params.P usages, and obvious
cache/output filename references so refactors can be planned from evidence.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

try:
    import tomllib  # type: ignore
except Exception:
    import tomli as tomllib  # type: ignore


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


PARAM_ATTR_PATTERNS = [
    re.compile(r"\b(?:self\.)?params\.P\.([A-Za-z_][A-Za-z0-9_]*)"),
    re.compile(r"\bP\.([A-Za-z_][A-Za-z0-9_]*)"),
    re.compile(r"getattr\([^,\n]*?\.P,\s*['\"]([A-Za-z_][A-Za-z0-9_]*)['\"]"),
    re.compile(r"getattr\(\s*P,\s*['\"]([A-Za-z_][A-Za-z0-9_]*)['\"]"),
]

CACHE_LITERAL_PATTERN = re.compile(
    r"""["']([^"'\r\n]*(?:cache|detect_|photometry_|apcorr|frame_quality|"""
    r"""ref_catalog|master_catalog|wcs|idmatch|\.json|\.csv|\.tsv|\.npy|\.npz)[^"'\r\n]*)["']""",
    re.IGNORECASE,
)


def load_toml(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("rb") as fh:
        return tomllib.load(fh)


def flatten_toml(data: dict[str, Any], prefix: tuple[str, ...] = ()) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in data.items():
        path = prefix + (str(key),)
        if isinstance(value, dict):
            out.update(flatten_toml(value, path))
        else:
            out[".".join(path)] = value
    return out


def path_key(path: Iterable[str]) -> str:
    return ".".join(str(part) for part in path)


def iter_py_files() -> Iterable[Path]:
    for base in (ROOT / "apex", ROOT / "tests"):
        if not base.exists():
            continue
        yield from base.rglob("*.py")


def scan_params_p_usage() -> tuple[Counter[str], dict[str, list[str]]]:
    counts: Counter[str] = Counter()
    locations: dict[str, list[str]] = defaultdict(list)
    for path in iter_py_files():
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            continue
        rel = path.relative_to(ROOT).as_posix()
        for lineno, line in enumerate(text.splitlines(), start=1):
            for pattern in PARAM_ATTR_PATTERNS:
                for match in pattern.finditer(line):
                    attr = match.group(1)
                    counts[attr] += 1
                    if len(locations[attr]) < 8:
                        locations[attr].append(f"{rel}:{lineno}")
    return counts, locations


def scan_cache_literals() -> dict[str, Counter[str]]:
    by_file: dict[str, Counter[str]] = {}
    for path in iter_py_files():
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            continue
        hits = Counter(match.group(1) for match in CACHE_LITERAL_PATTERN.finditer(text))
        if hits:
            by_file[path.relative_to(ROOT).as_posix()] = hits
    return by_file


def load_runtime_maps() -> dict[str, list[tuple[str, str]]]:
    from apex.config import parameters_cmd, parameters_lc

    return {
        "cmd": [(path_key(path), attr) for path, attr in parameters_cmd.TOML_KEY_MAP],
        "lc": [(path_key(path), attr) for path, attr in parameters_lc.TOML_KEY_MAP],
    }


def markdown_table(headers: list[str], rows: list[list[Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(cell).replace("|", "\\|") for cell in row) + " |")
    return "\n".join(lines)


def build_parameter_inventory() -> str:
    from apex.config.parameter_map import (
        CANONICAL_SCHEMA_VERSION,
        COMMON_FOUNDATION_MAP,
        LEGACY_ALIAS_MAP,
        duplicate_runtime_attrs,
        duplicate_toml_paths,
    )

    current = flatten_toml(load_toml(ROOT / "parameters.toml"))
    example = flatten_toml(load_toml(ROOT / "parameters.example.toml"))
    maps = load_runtime_maps()
    cmd_map = [(tuple(path.split(".")), attr) for path, attr in maps["cmd"]]
    lc_map = [(tuple(path.split(".")), attr) for path, attr in maps["lc"]]
    cmd_paths = {path for path, _ in maps["cmd"]}
    lc_paths = {path for path, _ in maps["lc"]}
    cmd_attrs = {attr for _, attr in maps["cmd"]}
    lc_attrs = {attr for _, attr in maps["lc"]}
    usage_counts, usage_locations = scan_params_p_usage()

    all_toml_keys = sorted(set(current) | set(example))
    toml_rows = []
    for key in all_toml_keys:
        value = current.get(key, example.get(key))
        toml_rows.append([
            key,
            type(value).__name__,
            "yes" if key in cmd_paths else "",
            "yes" if key in lc_paths else "",
            "yes" if key in current else "",
            "yes" if key in example else "",
        ])

    used_attrs = sorted(usage_counts)
    attr_rows = []
    for attr in used_attrs:
        attr_rows.append([
            attr,
            usage_counts[attr],
            "yes" if attr in cmd_attrs else "",
            "yes" if attr in lc_attrs else "",
            ", ".join(usage_locations.get(attr, [])[:3]),
        ])

    foundation_rows = [
        [entry.dotted_path, entry.attr, ",".join(sorted(entry.modes)), "yes" if entry.legacy else "", entry.note]
        for entry in (*COMMON_FOUNDATION_MAP, *LEGACY_ALIAS_MAP)
    ]
    duplicate_rows = []
    for mode, mapping in (("cmd", cmd_map), ("lc", lc_map)):
        for attr, paths in sorted(duplicate_runtime_attrs(mapping).items()):
            duplicate_rows.append([mode, "runtime attr", attr, ", ".join(paths)])
        for path, attrs in sorted(duplicate_toml_paths(mapping).items()):
            duplicate_rows.append([mode, "toml path", path, ", ".join(attrs)])

    lines = [
        "# Parameter Inventory",
        "",
        "Generated by `python scripts/inventory_foundation.py --write-docs`.",
        "",
        f"- canonical schema version: `{CANONICAL_SCHEMA_VERSION}`",
        f"- current TOML keys: `{len(current)}`",
        f"- example TOML keys: `{len(example)}`",
        f"- CMD mapped TOML paths: `{len(cmd_paths)}`",
        f"- LC mapped TOML paths: `{len(lc_paths)}`",
        f"- `params.P` attributes referenced in code: `{len(used_attrs)}`",
        "",
        "## Duplicate Map Diagnostics",
        "",
        markdown_table(["mode", "kind", "key", "mapped values"], duplicate_rows),
        "",
        "## Foundation Map Skeleton",
        "",
        markdown_table(["toml path", "runtime attr", "modes", "legacy", "note"], foundation_rows),
        "",
        "## TOML Key Coverage",
        "",
        markdown_table(["toml key", "type", "cmd map", "lc map", "current", "example"], toml_rows),
        "",
        "## Runtime `params.P` Usage",
        "",
        markdown_table(["attr", "uses", "cmd map", "lc map", "sample locations"], attr_rows),
        "",
    ]
    return "\n".join(lines)


def build_cache_inventory() -> str:
    cache_literals = scan_cache_literals()
    rows = []
    for file_path in sorted(cache_literals):
        for literal, count in sorted(cache_literals[file_path].items()):
            rows.append([file_path, literal, count])

    lines = [
        "# Cache and Output Inventory",
        "",
        "Generated by `python scripts/inventory_foundation.py --write-docs`.",
        "",
        "This is a string-literal inventory, not a final cache schema. It is meant",
        "to expose current cache/output naming before a shared cache manager is",
        "introduced.",
        "",
        f"- files with cache/output literals: `{len(cache_literals)}`",
        f"- distinct file/literal pairs: `{len(rows)}`",
        "",
        markdown_table(["file", "literal", "count"], rows),
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-docs", action="store_true", help="write docs/*.md inventories")
    args = parser.parse_args()

    parameter_doc = build_parameter_inventory()
    cache_doc = build_cache_inventory()

    if args.write_docs:
        docs = ROOT / "docs"
        docs.mkdir(parents=True, exist_ok=True)
        (docs / "parameter-inventory.md").write_text(parameter_doc, encoding="utf-8")
        (docs / "cache-inventory.md").write_text(cache_doc, encoding="utf-8")
    else:
        print(parameter_doc)
        print()
        print(cache_doc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
