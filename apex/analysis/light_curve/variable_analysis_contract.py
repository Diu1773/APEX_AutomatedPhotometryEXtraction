"""Versioned contract between the LC workflow and advanced variable analysis."""

from __future__ import annotations

import copy
import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np


BUNDLE_SCHEMA = "apex.variable-analysis.bundle"
BUNDLE_SCHEMA_VERSION = 1
RELEASE_STATUSES = {"APPROVED", "OVERRIDDEN", "BLOCKED", "UNVERIFIED"}
ANALYSIS_STATUSES = {"COMPLETE", "REVIEW_REQUIRED", "FAILED", "CANCELLED"}
ANALYSIS_BRANCHES = {"auto", "single", "multi"}

_PERIOD_ARRAY_KEYS = {
    "frequency",
    "power",
    "theta",
    "trial_periods",
    "top_periods",
    "top_powers",
    "time",
    "mag",
    "mag_err",
}

_ANALYSIS_ARRAY_KEYS = {
    "time",
    "mag",
    "mag_err",
    "model",
    "residual",
    "coeff",
    "components",
    "component_derivatives",
    "amplitudes",
    "phases",
    "fine_frequency",
    "fine_power",
    "bootstrap_periods",
}


def _json_safe(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return [_json_safe(item) for item in value.tolist()]
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _restore_period_arrays(scan_results: Mapping[str, Any]) -> dict[str, dict]:
    restored: dict[str, dict] = {}
    for method_key, raw_result in scan_results.items():
        if not isinstance(raw_result, Mapping):
            continue
        result = copy.deepcopy(dict(raw_result))
        for key in _PERIOD_ARRAY_KEYS:
            value = result.get(key)
            if isinstance(value, list):
                result[key] = np.asarray(value, dtype=float)
        restored[str(method_key)] = result
    return restored


def _restore_analysis_arrays(value: Any, key: str = "") -> Any:
    if isinstance(value, Mapping):
        return {
            str(child_key): _restore_analysis_arrays(child_value, str(child_key))
            for child_key, child_value in value.items()
        }
    if isinstance(value, list):
        if key in _ANALYSIS_ARRAY_KEYS:
            return np.asarray(value, dtype=float)
        return [_restore_analysis_arrays(item) for item in value]
    return value


def compute_file_fingerprint(path: Path, chunk_size: int = 1024 * 1024) -> dict[str, Any]:
    """Return a stable fingerprint for a released CSV without touching FITS data."""
    source = Path(path)
    stat = source.stat()
    digest = hashlib.sha256()
    with source.open("rb") as handle:
        while True:
            chunk = handle.read(max(int(chunk_size), 4096))
            if not chunk:
                break
            digest.update(chunk)
    return {
        "source_size": int(stat.st_size),
        "source_mtime_ns": int(stat.st_mtime_ns),
        "source_sha256": digest.hexdigest(),
    }


@dataclass
class ValidatedLightCurveBundle:
    """Main-workflow release consumed by the Variable Star Tool.

    Main workflow steps own data-quality decisions. The advanced tool receives
    those decisions and performs candidate-level model analysis only.
    """

    workspace_dir: str
    source_file: str
    target_id: int
    analysis_filter: str
    series_mode: str
    mag_col: str
    correction_mode: str
    correction_preserves_nightly_baseline: bool
    input_signature: dict[str, Any]
    adopted_period: float
    scan_results: dict[str, dict]
    alias_analysis: dict[str, Any]
    multimode_diagnostic: dict[str, Any]
    search: dict[str, Any]
    release_status: str = "UNVERIFIED"
    release_reasons: list[str] = field(default_factory=list)
    main_qc: dict[str, Any] = field(default_factory=dict)
    comparison_provenance: dict[str, Any] = field(default_factory=dict)
    photometry_provenance: dict[str, Any] = field(default_factory=dict)
    source: str = "step12_period_analysis"
    schema_version: int = BUNDLE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        self.release_status = str(self.release_status or "UNVERIFIED").upper()
        if self.release_status not in RELEASE_STATUSES:
            raise ValueError(f"Unknown variable-analysis release status: {self.release_status}")
        self.workspace_dir = str(self.workspace_dir)
        self.source_file = str(self.source_file)
        self.target_id = int(self.target_id)
        self.analysis_filter = str(self.analysis_filter or "__all__")
        self.adopted_period = float(self.adopted_period)
        self.release_reasons = [str(reason) for reason in self.release_reasons if str(reason).strip()]

    @property
    def can_launch(self) -> bool:
        return self.release_status != "BLOCKED"

    @property
    def release_message(self) -> str:
        if self.release_reasons:
            return "; ".join(self.release_reasons)
        return f"Main workflow release status: {self.release_status}"

    def to_dict(self, *, json_safe: bool = False) -> dict[str, Any]:
        payload = {
            "schema": BUNDLE_SCHEMA,
            "schema_version": int(self.schema_version),
            "source": self.source,
            "release": {
                "status": self.release_status,
                "reasons": list(self.release_reasons),
            },
            "target_series": {
                "workspace_dir": self.workspace_dir,
                "source_file": self.source_file,
                "target_id": int(self.target_id),
                "analysis_filter": self.analysis_filter,
                "series_mode": self.series_mode,
                "mag_col": self.mag_col,
                "correction_mode": self.correction_mode,
                "correction_preserves_nightly_baseline": bool(
                    self.correction_preserves_nightly_baseline
                ),
                "input_signature": copy.deepcopy(self.input_signature),
            },
            "main_qc": copy.deepcopy(self.main_qc),
            "provenance": {
                "comparison": copy.deepcopy(self.comparison_provenance),
                "photometry": copy.deepcopy(self.photometry_provenance),
            },
            "period_analysis": {
                "adopted_period": float(self.adopted_period),
                "scan_results": copy.deepcopy(self.scan_results),
                "alias_analysis": copy.deepcopy(self.alias_analysis),
                "multimode_diagnostic": copy.deepcopy(self.multimode_diagnostic),
                "search": copy.deepcopy(self.search),
            },
        }
        return _json_safe(payload) if json_safe else payload

    def to_legacy_handoff(self) -> dict[str, Any]:
        """Return the runtime shape consumed by the existing tool UI."""
        return {
            "version": 2,
            "bundle_schema": BUNDLE_SCHEMA,
            "bundle_schema_version": int(self.schema_version),
            "source": self.source,
            "release_status": self.release_status,
            "release_reasons": list(self.release_reasons),
            "main_qc": copy.deepcopy(self.main_qc),
            "comparison_provenance": copy.deepcopy(self.comparison_provenance),
            "photometry_provenance": copy.deepcopy(self.photometry_provenance),
            "workspace_dir": self.workspace_dir,
            "source_file": self.source_file,
            "target_id": int(self.target_id),
            "analysis_filter": self.analysis_filter,
            "series_mode": self.series_mode,
            "mag_col": self.mag_col,
            "correction_mode": self.correction_mode,
            "correction_preserves_nightly_baseline": bool(
                self.correction_preserves_nightly_baseline
            ),
            "input_signature": copy.deepcopy(self.input_signature),
            "adopted_period": float(self.adopted_period),
            "scan_results": _restore_period_arrays(self.scan_results),
            "alias_analysis": copy.deepcopy(self.alias_analysis),
            "multimode_diagnostic": copy.deepcopy(self.multimode_diagnostic),
            "search": copy.deepcopy(self.search),
        }

    def write_json(self, path: Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(self.to_dict(json_safe=True), indent=2, ensure_ascii=False, allow_nan=False),
            encoding="utf-8",
        )
        return path

    @classmethod
    def read_json(cls, path: Path) -> "ValidatedLightCurveBundle":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "ValidatedLightCurveBundle":
        data = dict(payload)
        if data.get("schema") == BUNDLE_SCHEMA:
            target = dict(data.get("target_series") or {})
            period = dict(data.get("period_analysis") or {})
            release = dict(data.get("release") or {})
            provenance = dict(data.get("provenance") or {})
            return cls(
                workspace_dir=str(target.get("workspace_dir", "")),
                source_file=str(target.get("source_file", "")),
                target_id=int(target.get("target_id", 0) or 0),
                analysis_filter=str(target.get("analysis_filter", "__all__")),
                series_mode=str(target.get("series_mode", "raw")),
                mag_col=str(target.get("mag_col", "")),
                correction_mode=str(target.get("correction_mode", "")),
                correction_preserves_nightly_baseline=bool(
                    target.get("correction_preserves_nightly_baseline", True)
                ),
                input_signature=dict(target.get("input_signature") or {}),
                adopted_period=float(period.get("adopted_period", np.nan)),
                scan_results=_restore_period_arrays(period.get("scan_results") or {}),
                alias_analysis=dict(period.get("alias_analysis") or {}),
                multimode_diagnostic=dict(period.get("multimode_diagnostic") or {}),
                search=dict(period.get("search") or {}),
                release_status=str(release.get("status", "UNVERIFIED")),
                release_reasons=list(release.get("reasons") or []),
                main_qc=dict(data.get("main_qc") or {}),
                comparison_provenance=dict(provenance.get("comparison") or {}),
                photometry_provenance=dict(provenance.get("photometry") or {}),
                source=str(data.get("source", "step12_period_analysis")),
                schema_version=int(data.get("schema_version", BUNDLE_SCHEMA_VERSION)),
            )

        # Compatibility with pre-contract Step 12 handoff dictionaries.
        return cls(
            workspace_dir=str(data.get("workspace_dir", "")),
            source_file=str(data.get("source_file", "")),
            target_id=int(data.get("target_id", 0) or 0),
            analysis_filter=str(data.get("analysis_filter", "__all__")),
            series_mode=str(data.get("series_mode", "raw")),
            mag_col=str(data.get("mag_col", "")),
            correction_mode=str(data.get("correction_mode", "")),
            correction_preserves_nightly_baseline=bool(
                data.get("correction_preserves_nightly_baseline", True)
            ),
            input_signature=dict(data.get("input_signature") or {}),
            adopted_period=float(data.get("adopted_period", np.nan)),
            scan_results=_restore_period_arrays(data.get("scan_results") or {}),
            alias_analysis=dict(data.get("alias_analysis") or {}),
            multimode_diagnostic=dict(data.get("multimode_diagnostic") or {}),
            search=dict(data.get("search") or {}),
            release_status=str(data.get("release_status", "UNVERIFIED")),
            release_reasons=list(data.get("release_reasons") or []),
            main_qc=dict(data.get("main_qc") or {}),
            comparison_provenance=dict(data.get("comparison_provenance") or {}),
            photometry_provenance=dict(data.get("photometry_provenance") or {}),
            source=str(data.get("source", "step12_period_analysis")),
            schema_version=int(data.get("bundle_schema_version", BUNDLE_SCHEMA_VERSION)),
        )


def coerce_validated_bundle(
    value: ValidatedLightCurveBundle | Mapping[str, Any],
) -> ValidatedLightCurveBundle:
    if isinstance(value, ValidatedLightCurveBundle):
        return value
    if isinstance(value, Mapping):
        return ValidatedLightCurveBundle.from_dict(value)
    raise TypeError("Variable-analysis handoff must be a bundle or dictionary.")


@dataclass
class ReviewRequired:
    """A scientific decision that must be made explicitly by the observer."""

    code: str
    message: str
    candidates: list[dict[str, Any]] = field(default_factory=list)
    allowed_actions: list[str] = field(default_factory=list)

    def to_dict(self, *, json_safe: bool = False) -> dict[str, Any]:
        payload = {
            "code": str(self.code),
            "message": str(self.message),
            "candidates": copy.deepcopy(self.candidates),
            "allowed_actions": list(self.allowed_actions),
        }
        return _json_safe(payload) if json_safe else payload

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "ReviewRequired":
        return cls(
            code=str(payload.get("code", "REVIEW_REQUIRED")),
            message=str(payload.get("message", "Observer review is required.")),
            candidates=[dict(row) for row in payload.get("candidates", []) if isinstance(row, Mapping)],
            allowed_actions=[str(action) for action in payload.get("allowed_actions", [])],
        )


@dataclass
class VariableAnalysisRequest:
    """Inputs for one deterministic advanced-analysis run."""

    bundle: ValidatedLightCurveBundle
    adopted_period_override: float | None = None
    secondary_period_override: float | None = None
    analysis_branch: str = "auto"
    bootstrap_resamples: int = 300
    refinement_harmonics: int = 4
    single_harmonics: int = 4
    multimode_harmonics: int = 2
    include_night_offsets: bool = False
    random_seed: int | None = None

    def __post_init__(self) -> None:
        self.bundle = coerce_validated_bundle(self.bundle)
        self.analysis_branch = str(self.analysis_branch or "auto").lower()
        if self.analysis_branch not in ANALYSIS_BRANCHES:
            raise ValueError(f"Unknown variable-analysis branch: {self.analysis_branch}")
        if self.adopted_period_override is not None:
            self.adopted_period_override = float(self.adopted_period_override)
            if not np.isfinite(self.adopted_period_override) or self.adopted_period_override <= 0:
                raise ValueError("The adopted-period override must be positive and finite.")
        if self.secondary_period_override is not None:
            self.secondary_period_override = float(self.secondary_period_override)
            if not np.isfinite(self.secondary_period_override) or self.secondary_period_override <= 0:
                raise ValueError("The secondary-period override must be positive and finite.")
        self.bootstrap_resamples = max(int(self.bootstrap_resamples), 0)
        self.refinement_harmonics = max(int(self.refinement_harmonics), 1)
        self.single_harmonics = max(int(self.single_harmonics), 1)
        self.multimode_harmonics = max(int(self.multimode_harmonics), 1)

    def to_dict(self, *, json_safe: bool = False) -> dict[str, Any]:
        payload = {
            "bundle": self.bundle.to_dict(json_safe=False),
            "adopted_period_override": self.adopted_period_override,
            "secondary_period_override": self.secondary_period_override,
            "analysis_branch": self.analysis_branch,
            "bootstrap_resamples": int(self.bootstrap_resamples),
            "refinement_harmonics": int(self.refinement_harmonics),
            "single_harmonics": int(self.single_harmonics),
            "multimode_harmonics": int(self.multimode_harmonics),
            "include_night_offsets": bool(self.include_night_offsets),
            "random_seed": self.random_seed,
        }
        return _json_safe(payload) if json_safe else payload


@dataclass
class VariableAnalysisResult:
    """Serializable output from the headless advanced-analysis service."""

    status: str
    branch: str = ""
    adopted_period: float | None = None
    refined_period: float | None = None
    local_period_sigma: float | None = None
    data_summary: dict[str, Any] = field(default_factory=dict)
    refinement: dict[str, Any] = field(default_factory=dict)
    model: dict[str, Any] = field(default_factory=dict)
    per_filter_models: dict[str, dict[str, Any]] = field(default_factory=dict)
    diagnostics: dict[str, Any] = field(default_factory=dict)
    review: ReviewRequired | None = None
    error: str = ""

    def __post_init__(self) -> None:
        self.status = str(self.status or "FAILED").upper()
        if self.status not in ANALYSIS_STATUSES:
            raise ValueError(f"Unknown variable-analysis result status: {self.status}")
        self.branch = str(self.branch or "").lower()

    @property
    def is_complete(self) -> bool:
        return self.status == "COMPLETE"

    def to_dict(self, *, json_safe: bool = False) -> dict[str, Any]:
        payload = {
            "status": self.status,
            "branch": self.branch,
            "adopted_period": self.adopted_period,
            "refined_period": self.refined_period,
            "local_period_sigma": self.local_period_sigma,
            "data_summary": copy.deepcopy(self.data_summary),
            "refinement": copy.deepcopy(self.refinement),
            "model": copy.deepcopy(self.model),
            "per_filter_models": copy.deepcopy(self.per_filter_models),
            "diagnostics": copy.deepcopy(self.diagnostics),
            "review": self.review.to_dict(json_safe=False) if self.review else None,
            "error": str(self.error or ""),
        }
        return _json_safe(payload) if json_safe else payload

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "VariableAnalysisResult":
        restored = _restore_analysis_arrays(dict(payload))
        review_payload = restored.get("review")
        return cls(
            status=str(restored.get("status", "FAILED")),
            branch=str(restored.get("branch", "")),
            adopted_period=restored.get("adopted_period"),
            refined_period=restored.get("refined_period"),
            local_period_sigma=restored.get("local_period_sigma"),
            data_summary=dict(restored.get("data_summary") or {}),
            refinement=dict(restored.get("refinement") or {}),
            model=dict(restored.get("model") or {}),
            per_filter_models={
                str(key): dict(value)
                for key, value in dict(restored.get("per_filter_models") or {}).items()
            },
            diagnostics=dict(restored.get("diagnostics") or {}),
            review=(
                ReviewRequired.from_dict(review_payload)
                if isinstance(review_payload, Mapping)
                else None
            ),
            error=str(restored.get("error", "")),
        )
