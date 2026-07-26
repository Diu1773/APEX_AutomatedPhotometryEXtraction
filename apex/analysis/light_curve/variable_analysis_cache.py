"""Persistent cache for headless advanced variable-analysis results."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from apex.analysis.light_curve.variable_analysis_contract import (
    VariableAnalysisRequest,
    VariableAnalysisResult,
    compute_file_fingerprint,
)
from apex.utils.step_paths_lc import step11_period_dir


CACHE_SCHEMA = "apex.variable-analysis.cache"
CACHE_VERSION = 1
ANALYSIS_ALGORITHM_VERSION = 5


def _resolve_source(request: VariableAnalysisRequest) -> Path | None:
    source = Path(request.bundle.source_file)
    candidates = [source]
    workspace = Path(request.bundle.workspace_dir)
    if not source.is_absolute():
        candidates.append(workspace / source)
    candidates.append(workspace / source.name)
    return next((path for path in candidates if path.is_file()), None)


def request_cache_key(request: VariableAnalysisRequest) -> str:
    source = _resolve_source(request)
    current_source = (
        compute_file_fingerprint(source)
        if source is not None
        else {"source_missing": True}
    )
    payload = {
        "algorithm_version": ANALYSIS_ALGORITHM_VERSION,
        "request": request.to_dict(json_safe=True),
        "current_source": current_source,
    }
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def cache_path(request: VariableAnalysisRequest) -> Path:
    cache_dir = step11_period_dir(Path(request.bundle.workspace_dir)) / "variable_analysis_cache"
    return cache_dir / f"{request_cache_key(request)}.json"


def load_cached_result(request: VariableAnalysisRequest) -> VariableAnalysisResult | None:
    path = cache_path(request)
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None
    if payload.get("schema") != CACHE_SCHEMA or int(payload.get("version", 0)) != CACHE_VERSION:
        return None
    if str(payload.get("cache_key", "")) != request_cache_key(request):
        return None
    result_payload = payload.get("result")
    if not isinstance(result_payload, dict):
        return None
    result = VariableAnalysisResult.from_dict(result_payload)
    return result if result.status == "COMPLETE" else None


def store_cached_result(
    request: VariableAnalysisRequest,
    result: VariableAnalysisResult,
) -> Path | None:
    if result.status != "COMPLETE":
        return None
    path = cache_path(request)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": CACHE_SCHEMA,
        "version": CACHE_VERSION,
        "cache_key": request_cache_key(request),
        "result": result.to_dict(json_safe=True),
    }
    temporary = path.with_suffix(".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False),
        encoding="utf-8",
    )
    temporary.replace(path)
    return path
