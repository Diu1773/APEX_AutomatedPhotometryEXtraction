"""A proxy that remembers which settings were asked for.

Lives in `utils` rather than beside the run manifest that first needed it,
because two layers record now and they must not import each other. The runner
stands one of these in front of the parameter namespace for the duration of a
step; forced photometry stands another inside each worker **process**, which
is the only way names read on the far side of a process boundary come back —
the parent's proxy is copied into the child by pickle and its additions stay
there. `apex.pipeline.provenance` re-exports this for its existing callers.
"""

from __future__ import annotations

from typing import Any


class RecordingNamespace:
    """Delegates to the real namespace and remembers what was asked for.

    Deliberately not a subclass: `SimpleNamespace` attribute access is what the
    whole codebase uses, so wrapping it is enough, and staying out of the type
    hierarchy means an `isinstance` check somewhere cannot change behaviour.
    """

    def __init__(self, target: Any) -> None:
        object.__setattr__(self, "_target", target)
        object.__setattr__(self, "_seen", set())

    def __getattr__(self, name: str) -> Any:
        target = object.__getattribute__(self, "_target")
        if not name.startswith("_"):
            object.__getattribute__(self, "_seen").add(name)
        return getattr(target, name)

    def __setattr__(self, name: str, value: Any) -> None:
        setattr(object.__getattribute__(self, "_target"), name, value)

    def __dir__(self):
        return dir(object.__getattribute__(self, "_target"))

    @property
    def seen(self) -> set[str]:
        return set(object.__getattribute__(self, "_seen"))
