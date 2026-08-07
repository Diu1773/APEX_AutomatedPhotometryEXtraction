"""Runtime context passed to every pipeline step.

Carries the resolved parameters, paths, logger, and optional ProjectState.
Deliberately Qt-free so it works in CI, on servers, and inside tests.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional


@dataclass
class RunContext:
    mode: str                      # "cmd" | "lc"
    params: Any                    # apex.config.parameters_{cmd,lc}.Parameters
    result_dir: Path
    data_dir: Path
    logger: logging.Logger
    force: bool = False
    dry_run: bool = False
    project_state: Optional[Any] = None
    # Config-driven substitutes for GUI interactions (e.g. crop rectangle).
    extra: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def build(
        cls,
        mode: str,
        config_path: str | Path = "parameters.toml",
        *,
        result_dir: str | Path | None = None,
        data_dir: str | Path | None = None,
        force: bool = False,
        dry_run: bool = False,
        logger: logging.Logger | None = None,
    ) -> "RunContext":
        if mode not in ("cmd", "lc"):
            raise ValueError(f"mode must be 'cmd' or 'lc', got {mode!r}")

        if mode == "cmd":
            from apex.config.parameters_cmd import read_params
        else:
            from apex.config.parameters_lc import read_params

        config_path = Path(config_path)
        if not config_path.exists():
            raise FileNotFoundError(
                f"Config not found: {config_path}. Run 'apex config init' to create one."
            )
        params = read_params(config_path)

        resolved_result = Path(result_dir or getattr(params.P, "result_dir", "") or "result")
        resolved_data = Path(data_dir or getattr(params.P, "data_dir", "") or ".")

        # Make the resolved paths authoritative for the whole run: components
        # such as FileManager read params.P.result_dir / params.P.data_dir
        # directly, so CLI overrides must be reflected there too.
        configured_result = Path(getattr(params.P, "result_dir", "") or "")
        configured_cache = Path(getattr(params.P, "cache_dir", "") or "")
        try:
            params.P.result_dir = resolved_result
            params.P.data_dir = resolved_data
        except Exception:  # noqa: BLE001 - never let path sync break a run
            pass

        # `read_params` has already resolved cache_dir against the *config's*
        # result_dir, so overriding result_dir alone leaves the detection and
        # WCS caches behind in the config's tree — two runs with different
        # --result-dir then share one cache and can reuse each other's
        # results.  Re-root it whenever it lived under the old result dir; an
        # absolute cache_dir configured outside it is a deliberate choice and
        # is left alone.
        if configured_cache and configured_result:
            try:
                relative = configured_cache.relative_to(configured_result)
            except ValueError:
                relative = None
            if relative is not None:
                params.P.cache_dir = resolved_result / relative
                params.P.cache_dir.mkdir(parents=True, exist_ok=True)

        return cls(
            mode=mode,
            params=params,
            result_dir=resolved_result,
            data_dir=resolved_data,
            logger=logger or logging.getLogger("apex.pipeline"),
            force=force,
            dry_run=dry_run,
        )
