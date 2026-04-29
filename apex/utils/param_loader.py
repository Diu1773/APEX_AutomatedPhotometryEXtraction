"""Resolve and persist the last-used parameters.toml path for APEX modes."""
from __future__ import annotations
from pathlib import Path
from PyQt5.QtWidgets import QFileDialog, QWidget

_CACHE_FILE = Path(__file__).parent.parent.parent.parent / ".apex_last_param.txt"


def _load_cached() -> Path | None:
    try:
        p = Path(_CACHE_FILE.read_text(encoding="utf-8").strip())
        return p if p.exists() else None
    except Exception:
        return None


def _save_cached(path: Path) -> None:
    try:
        _CACHE_FILE.write_text(str(path), encoding="utf-8")
    except Exception:
        pass


def resolve_param_file(parent: QWidget | None = None, param_file: str | None = None) -> Path:
    """Return a valid parameters.toml path, prompting only when necessary."""
    # 1. explicit arg
    if param_file:
        p = Path(param_file)
        if p.exists():
            _save_cached(p)
            return p

    # 2. parameters.toml next to cwd
    local = Path("parameters.toml")
    if local.exists():
        _save_cached(local.resolve())
        return local.resolve()

    # 3. last used
    cached = _load_cached()
    if cached:
        return cached

    # 4. file dialog
    chosen, _ = QFileDialog.getOpenFileName(
        parent, "APEX — parameters.toml 선택",
        str(Path.home()),
        "TOML Files (*.toml);;All Files (*.*)",
    )
    if not chosen:
        raise RuntimeError("파라미터 파일이 선택되지 않았습니다.")
    p = Path(chosen)
    _save_cached(p)
    return p
