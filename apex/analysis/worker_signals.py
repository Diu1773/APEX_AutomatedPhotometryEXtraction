"""How long-running calculations report progress without knowing who is asking.

Steps 8 and 10 do their work inside `QThread` subclasses, and the headless
pipeline imported those GUI classes to get at them — so a script could not do
CMD photometry unless PyQt5 was installed. The calculation never needed Qt: the
two workers use no `QThread` method at all, only the base class and a handful of
`pyqtSignal` declarations.

So the calculation lives in `apex.analysis` and announces what it is doing on
named `Channel`s. Whoever drives it subscribes:

    runner.on_progress.subscribe(lambda i, n, name: print(f"{i}/{n} {name}"))

The window subscribes a Qt signal's `emit`, so the line still crosses to the GUI
thread the way it always did. A script subscribes a logger, or nothing at all.
The calculation does not know the difference, which is the point — the same
object runs in both, so "the app and the script agree" is not a claim that needs
testing, it is the same code.

Log lines go through `on_log` rather than the `logging` module. `logging` would
be the obvious home for them, but the caller would then receive each line twice
(once through its own subscription, once through the root logger) unless every
entry point agreed on handler configuration. One channel, one delivery, and a
caller that wants `logging` writes a one-line subscriber (2026-08-16).
"""

from __future__ import annotations

from typing import Any, Callable


class Channel:
    """One named thing a run announces, and the callbacks waiting for it.

    Deliberately the shape of an observer, not of Qt: `subscribe`/`send` rather
    than `connect`/`emit`, so nothing in `apex.analysis` reads as if it were
    talking to a widget toolkit.
    """

    __slots__ = ("name", "_subscribers")

    def __init__(self, name: str = "") -> None:
        self.name = name
        self._subscribers: list[Callable[..., Any]] = []

    def subscribe(self, callback: Callable[..., Any]) -> None:
        self._subscribers.append(callback)

    def unsubscribe(self, callback: Callable[..., Any] | None = None) -> None:
        if callback is None:
            self._subscribers.clear()
        elif callback in self._subscribers:
            self._subscribers.remove(callback)

    def send(self, *args: Any) -> None:
        """Announce. A subscriber that raises must not end the run — losing an
        hour of photometry because a progress bar threw is not a trade anyone
        would make — so failures here are dropped, as a queued Qt slot's would
        be."""
        for callback in list(self._subscribers):
            try:
                callback(*args)
            except Exception:
                pass

    def __bool__(self) -> bool:
        return bool(self._subscribers)


class ReportsProgress:
    """Base for a calculation that announces progress on named channels.

    Subclasses list the channel names in `_CHANNELS`; each becomes a `Channel`
    attribute prefixed with `on_`, so `progress` is reached as `on_progress` and
    cannot collide with a `pyqtSignal` of the same name in a GUI subclass.
    """

    _CHANNELS: tuple[str, ...] = ()

    def __init__(self) -> None:
        super().__init__()
        for name in self._CHANNELS:
            setattr(self, f"on_{name}", Channel(name))

    def stop(self) -> None:
        """Ask the run to end at its next checkpoint."""
        self._stop_requested = True

    @property
    def stop_requested(self) -> bool:
        return bool(getattr(self, "_stop_requested", False))
