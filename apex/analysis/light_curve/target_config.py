"""Which star the light curve is of — read from the config, not clicked.

LC's calculation has been Qt-free all along: fourteen services in
`apex.analysis.light_curve`, some seven thousand lines, none of them importing
Qt. What kept the whole branch in the window is smaller and harder than a port —
nothing in the configuration could say *which star*, and a batch run has nobody
to click it.

That is the same gate the isochrone fit sat behind until 2026-08-17, and it gets
the same answer: config rows, and a refusal rather than a guess. There is no
defensible default here. Picking the brightest star, or the one nearest the
field centre, would produce a confident light curve of the wrong object — and a
light curve does not look wrong, it just belongs to something else.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd


@dataclass
class LcTarget:
    """The target and its comparison ensemble, resolved from the config."""

    target_id: int
    target_name: str = ""
    comparison_ids: list[int] = field(default_factory=list)
    comparison_mode: str = "auto"
    comparison_count: int = 10
    filter_key: str = ""

    @property
    def is_resolved(self) -> bool:
        return self.target_id > 0


def _ids(raw: str) -> list[int]:
    out = []
    for token in str(raw or "").replace(";", ",").split(","):
        token = token.strip()
        if not token:
            continue
        try:
            out.append(int(float(token)))
        except ValueError:
            continue
    return out


def read_target(params) -> LcTarget:
    """The config's answer, unvalidated — `missing_target_settings` judges it."""
    P = getattr(params, "P", params)
    return LcTarget(
        target_id=int(getattr(P, "lc_target_id", -1) or -1),
        target_name=str(getattr(P, "lc_target_name", "") or "").strip(),
        comparison_ids=_ids(getattr(P, "lc_comparison_ids", "")),
        comparison_mode=str(getattr(P, "lc_comparison_mode", "auto") or "auto").strip().lower(),
        comparison_count=int(getattr(P, "lc_comparison_count", 10) or 10),
        filter_key=str(getattr(P, "lc_filter", "") or "").strip(),
    )


def missing_target_settings(params) -> list[str]:
    """What must be written down before a batch light curve means anything."""
    target = read_target(params)
    missing = []
    if not target.is_resolved:
        missing.append("lightcurve.target_id (Step 8 에서 고른 별의 ID)")
    if target.comparison_mode == "manual" and not target.comparison_ids:
        missing.append("lightcurve.comparison_ids "
                       "(comparison_mode=manual 이면 비교성을 적어야 한다)")
    return missing


def resolve_comparisons(target: LcTarget, catalog: pd.DataFrame,
                        stability: pd.DataFrame | None = None) -> list[int]:
    """The comparison ensemble this run should use.

    `manual` takes the config's list verbatim — the user chose those stars and a
    batch run must not quietly substitute others. `auto` ranks by the stability
    metrics Step 8 already computes when they are available, and falls back to
    the catalogue order when they are not, so a run without a stability table
    still produces something reproducible rather than nothing.
    """
    if target.comparison_mode == "manual":
        return list(target.comparison_ids)

    if catalog is None or catalog.empty or "ID" not in catalog.columns:
        return []
    available = [int(v) for v in pd.to_numeric(catalog["ID"], errors="coerce").dropna()
                 if int(v) != target.target_id]

    if stability is not None and not stability.empty and "ID" in stability.columns:
        score_col = next((c for c in ("stability_score", "score", "rms")
                          if c in stability.columns), None)
        if score_col is not None:
            ranked = stability.copy()
            ranked["_score"] = pd.to_numeric(ranked[score_col], errors="coerce")
            ascending = score_col == "rms"        # rms: smaller is better
            ranked = ranked.dropna(subset=["_score"]).sort_values(
                "_score", ascending=ascending)
            ordered = [int(v) for v in ranked["ID"] if int(v) in set(available)]
            if ordered:
                return ordered[:max(target.comparison_count, 1)]

    return available[:max(target.comparison_count, 1)]


def write_selection(result_dir, target: LcTarget, comparisons: list[int],
                    *, check_id: int | None = None,
                    selected_by: str = "config") -> Path:
    """Persist the choice where the downstream LC steps look for it.

    `check_id` and `selected_by` are recorded because a file that says only
    which stars were used cannot tell an ensemble ranked by stability from one
    taken in catalogue order — and those two gave averages 0.68 mag apart on the
    same frames.
    """
    import json

    from apex.utils.step_paths_lc import step8_selection_dir

    out_dir = Path(step8_selection_dir(result_dir))
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "lc_target_selection.json"
    payload = {
        "target_id": target.target_id,
        "target_name": target.target_name,
        "filter": target.filter_key,
        "comparison_mode": target.comparison_mode,
        "comparison_ids": [int(v) for v in comparisons],
        "source": "config",
        "selected_by": selected_by,
    }
    if check_id is not None:
        payload["check_id"] = int(check_id)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False),
                    encoding="utf-8")
    return path


def read_selection(result_dir) -> dict | None:
    """The selection a previous step wrote, whether by config or by window."""
    import json

    from apex.utils.step_paths_lc import step8_selection_dir

    path = Path(step8_selection_dir(result_dir)) / "lc_target_selection.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:                               # noqa: BLE001
        return None
