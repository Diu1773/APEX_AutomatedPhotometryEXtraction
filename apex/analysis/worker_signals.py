"""Qt-shaped callbacks for worker code that must also run without Qt.

Steps 8 and 10 do their computation inside `QThread` subclasses, and the
headless pipeline reached into `apex.gui` to get at them — so a script could not
do CMD photometry unless PyQt5 was installed, and two steps reported
NOT_IMPLEMENTED on a Qt-free install. The computation itself never needed Qt:
between them the two workers use no `QThread` method at all, only the base class
and twelve `pyqtSignal` declarations they emit on 38 times.

So the workers move down to `apex.analysis` and keep their bodies exactly as
they were. `self.progress.emit(...)` still reads the same, because:

- headless, `SignalHost.__init__` gives each name a `Signal`, which is a plain
  list of callbacks with Qt's `emit`/`connect` shape;
- in the GUI, the subclass is `class Worker(QThread, Runner)` and declares real
  `pyqtSignal`s. Those are class attributes, so `hasattr` already finds them and
  `SignalHost` leaves them alone — the Qt signal is what gets emitted, across
  threads, exactly as before.

One class, two drivers, one code path. That is also what makes the parity claim
in the paper literal rather than approximate: the script and the window are not
running equivalent code, they are running the same object.
"""

from __future__ import annotations

from typing import Any, Callable


class Signal:
    """A callback list with the shape of `pyqtSignal`.

    Only what the workers use: `connect` to register, `emit` to call. A slot
    that raises must not take the run down with it — a progress line failing is
    not a reason to lose an hour of photometry — so exceptions are swallowed
    the way a Qt slot invoked across a thread boundary would swallow them.
    """

    __slots__ = ("_slots",)

    def __init__(self) -> None:
        self._slots: list[Callable[..., Any]] = []

    def connect(self, slot: Callable[..., Any]) -> None:
        self._slots.append(slot)

    def disconnect(self, slot: Callable[..., Any] | None = None) -> None:
        if slot is None:
            self._slots.clear()
        elif slot in self._slots:
            self._slots.remove(slot)

    def emit(self, *args: Any) -> None:
        for slot in list(self._slots):
            try:
                slot(*args)
            except Exception:
                pass

    def __bool__(self) -> bool:
        return bool(self._slots)


class SignalHost:
    """Base for worker code that emits progress but must not require Qt.

    Subclasses list the names they emit on in `_SIGNALS`. Each becomes a
    `Signal` unless the instance already has one — which is how the GUI's
    `pyqtSignal` declarations survive: they are class attributes, so they are
    found first and used unchanged.
    """

    _SIGNALS: tuple[str, ...] = ()

    def __init__(self) -> None:
        super().__init__()
        for name in self._SIGNALS:
            if getattr(self, name, None) is None:
                setattr(self, name, Signal())

    def stop(self) -> None:
        """Ask the run to end at its next checkpoint."""
        self._stop_requested = True

    @property
    def stop_requested(self) -> bool:
        return bool(getattr(self, "_stop_requested", False))
