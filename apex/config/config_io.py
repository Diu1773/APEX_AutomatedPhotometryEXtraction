"""Workspace configuration I/O — JSON authoritative, legacy TOML migrated.

Why this exists (2026-08-06, user decision "TOML 완전 제거"): the runtime
config used to be a TOML file next to the executable. That design caused two
recurring accidents:

* **one file, many targets** — switching clusters meant editing the same
  ``parameters.toml`` in place, and half-edits left mixed identities
  (io paths → M13 while ``[target]`` → NGC 6811, observed 2026-08-06);
* **lossy, fragile writes** — the TOML round-trip dropped comments and the
  escaped-backslash Windows paths broke ad-hoc edits twice in one session.

The cure is structural, not cosmetic: every workspace keeps its own
``apex_config.json`` and *that* file is the single source of truth for the
GUI, the headless runners and the CLI alike. Legacy TOML files are read once,
converted to JSON alongside, and never consulted again (a newer mtime on the
TOML only earns a warning — silently re-importing it would resurrect the
two-sources-of-truth disease this module exists to kill).

Mapping rule for legacy files — collision-free by construction, because one
directory routinely holds several variant TOMLs (``parameters.toml``,
``parameters_result_psf.toml``, …):

======================================  ====================================
legacy file                             JSON authority
======================================  ====================================
``parameters.toml``                     ``apex_config.json``
``parameters_<suffix>.toml``            ``apex_config_<suffix>.json``
``<other>.toml``                        ``<other>.json``
``<dir>/``                              ``<dir>/apex_config.json``
======================================  ====================================
"""

from __future__ import annotations

import json
import os
import tempfile
import warnings
from pathlib import Path
from typing import Any, Tuple

#: Canonical basename of a workspace configuration file.
CONFIG_BASENAME = "apex_config.json"

#: Legacy basename this module migrates away from.
LEGACY_BASENAME = "parameters.toml"


def _json_name_for_legacy(toml_path: Path) -> str:
    stem = toml_path.stem
    if stem == "parameters":
        return CONFIG_BASENAME
    if stem.startswith("parameters_"):
        return f"apex_config_{stem[len('parameters_'):]}.json"
    return f"{stem}.json"


def resolve_config_path(path: str | Path) -> Path:
    """Return the JSON path that is (or will become) authoritative for ``path``.

    Accepts a directory, a ``.json`` path, or a legacy ``.toml`` path.
    """
    p = Path(path)
    if p.suffix.lower() == ".json":
        return p
    if p.suffix.lower() == ".toml":
        return p.with_name(_json_name_for_legacy(p))
    if p.is_dir() or not p.suffix:
        return p / CONFIG_BASENAME
    return p.with_suffix(".json")


def _legacy_candidate(path: str | Path, json_path: Path) -> Path | None:
    """The TOML file that would seed ``json_path``, if any."""
    p = Path(path)
    if p.suffix.lower() == ".toml":
        return p
    # For a directory / json request, only the canonical legacy name applies.
    cand = json_path.parent / LEGACY_BASENAME
    if json_path.name == CONFIG_BASENAME:
        return cand
    return None


def _read_json(path: Path) -> dict:
    text = path.read_bytes().decode("utf-8-sig")
    data = json.loads(text)
    if not isinstance(data, dict):
        raise ValueError(f"config root must be an object: {path}")
    return data


def save_config_data(json_path: str | Path, data: dict) -> bool:
    """Atomically write ``data`` as the workspace JSON. Returns success."""
    json_path = Path(json_path)
    try:
        json_path.parent.mkdir(parents=True, exist_ok=True)
        payload = json.dumps(data, indent=2, ensure_ascii=False, default=str)
        fd, tmp_name = tempfile.mkstemp(
            dir=str(json_path.parent), prefix=json_path.stem + ".", suffix=".tmp")
        try:
            with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(payload + "\n")
            os.replace(tmp_name, json_path)
        except Exception:
            try:
                os.unlink(tmp_name)
            except OSError:
                pass
            raise
        return True
    except Exception as exc:  # pragma: no cover - disk-level failures
        warnings.warn(f"apex_config write failed ({json_path}): {exc}")
        return False


def load_config_data(path: str | Path) -> Tuple[dict, Path]:
    """Load workspace config, migrating a legacy TOML on first touch.

    Returns ``(data, json_path)`` where ``json_path`` is the authoritative
    file from now on. Missing everything yields ``({}, json_path)`` so the
    caller applies its defaults (same contract the TOML reader had).
    """
    json_path = resolve_config_path(path)
    legacy = _legacy_candidate(path, json_path)

    if json_path.exists():
        if legacy is not None and legacy.exists():
            try:
                if legacy.stat().st_mtime > json_path.stat().st_mtime + 1.0:
                    warnings.warn(
                        f"legacy TOML is newer than its JSON authority and is "
                        f"IGNORED (edit {json_path.name} instead): {legacy}")
            except OSError:
                pass
        return _read_json(json_path), json_path

    if legacy is not None and legacy.exists():
        from apex.utils.io_utils import load_toml  # legacy reader (BOM tolerant)

        data = load_toml(legacy)
        data.setdefault("_meta", {})["migrated_from"] = legacy.name
        if save_config_data(json_path, data):
            return data, json_path
        # Could not persist the migration — still hand back the data so a
        # read-only filesystem does not brick the run.
        return data, json_path

    return {}, json_path


def migrate_config_path(path: str | Path) -> Path:
    """Ensure ``path`` has a JSON authority and return it (data discarded)."""
    _, json_path = load_config_data(path)
    return json_path


def _ident_token(text: str) -> str:
    """Lowercase alnum-only form for fuzzy identity matching (NGC 6811 ≡ ngc6811)."""
    return "".join(ch for ch in str(text).lower() if ch.isalnum())


def check_workspace_identity(data: dict, json_path: Path | str) -> list[str]:
    """Heuristic warnings when a config's identity looks self-inconsistent.

    The historical failure mode: one shared file edited in place until the io
    paths pointed at M13 while ``[target] name`` still said NGC 6811 — the
    window title and the plotted cluster disagreed. This check never blocks a
    run (survey fields legitimately have arbitrary names); it only names the
    mismatch so the user can fix the config instead of chasing ghosts.
    """
    issues: list[str] = []
    io_block = data.get("io") or {}
    target = data.get("target") or {}
    name = str(target.get("name") or "").strip()
    data_dir = str(io_block.get("data_dir") or "")
    result_dir = str(io_block.get("result_dir") or "")

    if not result_dir:
        issues.append(f"io.result_dir is empty ({Path(json_path).name})")
    if name:
        tok = _ident_token(name)
        hay = _ident_token(data_dir) + _ident_token(result_dir)
        if tok and len(tok) >= 3 and tok not in hay:
            issues.append(
                f"target.name={name!r} does not appear in io paths "
                f"(data_dir={data_dir or '?'}, result_dir={result_dir or '?'}) "
                f"— this config may mix two objects")
    return issues
