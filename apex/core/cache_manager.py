"""
Shared cache manifest primitives for workflow steps.

The existing pipeline has many step-specific cache files. This module does not
replace those payloads yet; it provides a small, explicit manifest layer that
new or migrated steps can use to describe cache validity consistently.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Iterable, Mapping


CACHE_MANIFEST_VERSION = 1


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize_cache_path(path: str | Path) -> str:
    """Return a stable absolute path string for cache metadata."""

    return str(Path(path).expanduser().resolve())


def safe_cache_key(key: str) -> str:
    """Convert an arbitrary cache key to a filesystem-safe stem."""

    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(key).strip())
    cleaned = cleaned.strip("._")
    if cleaned:
        return cleaned[:180]
    return hashlib.sha1(str(key).encode("utf-8", errors="ignore")).hexdigest()


def file_sha1(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha1()
    with Path(path).open("rb") as fh:
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


@dataclass(frozen=True)
class CacheInputSignature:
    """File identity used for cache validation."""

    path: str
    size: int | None = None
    mtime_ns: int | None = None
    sha1: str | None = None

    @classmethod
    def from_path(cls, path: str | Path, hash_content: bool = False) -> "CacheInputSignature":
        p = Path(path)
        normalized = normalize_cache_path(p)
        try:
            stat = p.stat()
            size = int(stat.st_size)
            mtime_ns = int(stat.st_mtime_ns)
        except OSError:
            size = None
            mtime_ns = None
        digest = file_sha1(p) if hash_content and p.exists() else None
        return cls(path=normalized, size=size, mtime_ns=mtime_ns, sha1=digest)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "CacheInputSignature":
        return cls(
            path=str(data.get("path", "")),
            size=_optional_int(data.get("size")),
            mtime_ns=_optional_int(data.get("mtime_ns")),
            sha1=str(data["sha1"]) if data.get("sha1") else None,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "size": self.size,
            "mtime_ns": self.mtime_ns,
            "sha1": self.sha1,
        }

    def matches(self, other: "CacheInputSignature", allow_mtime_drift: bool = False) -> bool:
        if self.path != other.path:
            return False
        if self.size != other.size:
            return False
        if self.sha1 or other.sha1:
            return bool(self.sha1 and other.sha1 and self.sha1 == other.sha1)
        if allow_mtime_drift:
            return True
        return self.mtime_ns == other.mtime_ns


@dataclass
class CacheManifest:
    """Standard cache metadata envelope."""

    step_id: str
    cache_schema_version: int
    manifest_version: int = CACHE_MANIFEST_VERSION
    created_at: str = field(default_factory=utc_now_iso)
    parameter_hash: str | None = None
    input_files: list[CacheInputSignature] = field(default_factory=list)
    payload_paths: dict[str, str] = field(default_factory=dict)
    dependency_versions: dict[str, str] = field(default_factory=dict)
    extra: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "CacheManifest":
        input_files = [
            CacheInputSignature.from_dict(item)
            for item in data.get("input_files", [])
            if isinstance(item, Mapping)
        ]
        payload_paths = {
            str(key): str(value)
            for key, value in dict(data.get("payload_paths", {})).items()
        }
        dependency_versions = {
            str(key): str(value)
            for key, value in dict(data.get("dependency_versions", {})).items()
        }
        extra = dict(data.get("extra", {})) if isinstance(data.get("extra", {}), Mapping) else {}
        return cls(
            step_id=str(data.get("step_id", "")),
            cache_schema_version=int(data.get("cache_schema_version", 0) or 0),
            manifest_version=int(data.get("manifest_version", CACHE_MANIFEST_VERSION) or CACHE_MANIFEST_VERSION),
            created_at=str(data.get("created_at") or utc_now_iso()),
            parameter_hash=str(data["parameter_hash"]) if data.get("parameter_hash") else None,
            input_files=input_files,
            payload_paths=payload_paths,
            dependency_versions=dependency_versions,
            extra=extra,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "manifest_version": int(self.manifest_version),
            "cache_schema_version": int(self.cache_schema_version),
            "step_id": self.step_id,
            "created_at": self.created_at,
            "parameter_hash": self.parameter_hash,
            "input_files": [item.to_dict() for item in self.input_files],
            "payload_paths": dict(self.payload_paths),
            "dependency_versions": dict(self.dependency_versions),
            "extra": dict(self.extra),
        }


@dataclass(frozen=True)
class CacheValidationResult:
    valid: bool
    reason: str = "ok"
    manifest: CacheManifest | None = None


class StepCacheManager:
    """Manage cache manifests for one workflow step."""

    def __init__(
        self,
        root_dir: str | Path,
        step_id: str,
        cache_schema_version: int,
        parameter_hash: str | None = None,
        hash_inputs: bool = False,
    ):
        self.root_dir = Path(root_dir)
        self.step_id = str(step_id)
        self.cache_schema_version = int(cache_schema_version)
        self.parameter_hash = parameter_hash
        self.hash_inputs = bool(hash_inputs)

    def manifest_path(self, key: str) -> Path:
        return self.root_dir / "_manifests" / f"{safe_cache_key(key)}.manifest.json"

    def build_manifest(
        self,
        input_paths: Iterable[str | Path] = (),
        payload_paths: Mapping[str, str | Path] | None = None,
        dependency_versions: Mapping[str, Any] | None = None,
        extra: Mapping[str, Any] | None = None,
    ) -> CacheManifest:
        return CacheManifest(
            step_id=self.step_id,
            cache_schema_version=self.cache_schema_version,
            parameter_hash=self.parameter_hash,
            input_files=[
                CacheInputSignature.from_path(path, hash_content=self.hash_inputs)
                for path in input_paths
            ],
            payload_paths={
                str(key): normalize_cache_path(value)
                for key, value in dict(payload_paths or {}).items()
            },
            dependency_versions={
                str(key): str(value)
                for key, value in dict(dependency_versions or {}).items()
            },
            extra=dict(extra or {}),
        )

    def write_manifest(self, key: str, manifest: CacheManifest) -> Path:
        path = self.manifest_path(key)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as fh:
            json.dump(manifest.to_dict(), fh, indent=2, sort_keys=True)
        return path

    def read_manifest(self, key: str) -> CacheManifest | None:
        path = self.manifest_path(key)
        if not path.exists():
            return None
        try:
            with path.open("r", encoding="utf-8") as fh:
                return CacheManifest.from_dict(json.load(fh))
        except Exception:
            return None

    def validate_manifest(
        self,
        manifest: CacheManifest | None,
        input_paths: Iterable[str | Path] = (),
        required_payloads: Iterable[str] = (),
        allow_mtime_drift: bool = False,
    ) -> CacheValidationResult:
        if manifest is None:
            return CacheValidationResult(False, "missing_manifest", manifest)
        if manifest.manifest_version != CACHE_MANIFEST_VERSION:
            return CacheValidationResult(False, "manifest_version_mismatch", manifest)
        if manifest.step_id != self.step_id:
            return CacheValidationResult(False, "step_id_mismatch", manifest)
        if manifest.cache_schema_version != self.cache_schema_version:
            return CacheValidationResult(False, "cache_schema_mismatch", manifest)
        if self.parameter_hash and manifest.parameter_hash != self.parameter_hash:
            return CacheValidationResult(False, "parameter_hash_mismatch", manifest)

        expected = [
            CacheInputSignature.from_path(path, hash_content=self.hash_inputs)
            for path in input_paths
        ]
        if len(expected) != len(manifest.input_files):
            return CacheValidationResult(False, "input_count_mismatch", manifest)
        for saved, current in zip(manifest.input_files, expected):
            if not saved.matches(current, allow_mtime_drift=allow_mtime_drift):
                return CacheValidationResult(False, "input_signature_mismatch", manifest)

        for key in required_payloads:
            payload = manifest.payload_paths.get(str(key))
            if not payload:
                return CacheValidationResult(False, f"missing_payload:{key}", manifest)
            if not Path(payload).exists():
                return CacheValidationResult(False, f"payload_not_found:{key}", manifest)

        return CacheValidationResult(True, "ok", manifest)

    def validate_key(
        self,
        key: str,
        input_paths: Iterable[str | Path] = (),
        required_payloads: Iterable[str] = (),
        allow_mtime_drift: bool = False,
    ) -> CacheValidationResult:
        return self.validate_manifest(
            self.read_manifest(key),
            input_paths=input_paths,
            required_payloads=required_payloads,
            allow_mtime_drift=allow_mtime_drift,
        )


def _optional_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except Exception:
        return None
