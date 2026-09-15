"""파라미터 창을 여는 함수가 없는 이름을 참조하면 앱이 통째로 끝난다.

2026-09-15 에 CMD 3·5 번에서 실제로 그랬다.

    5 번   `FittedDialog` 을 세 곳에서 쓰는데 임포트가 없다 → NameError
    3 번   `run_param_dialog(...)` 뒤에 옛 코드 두 줄이 남아 `buttons`·`dialog` 를
           만진다 → NameError

**증상이 예외로 안 보이는 것이 이 결함의 성질이다.** Qt 슬롯 안에서 처리되지 않은
파이썬 예외가 나면 PyQt5 가 프로세스를 끝내므로, 화면에는 앱이 그냥 사라지고 윈도는
`0xC0000409`(스택 버퍼 오버런)로 적는다. `faulthandler` 도 못 잡고 트레이스백도 안
남아서, 처음에는 원인을 측정 도구 쪽으로 잘못 짚었다.

그래서 **개별 두 줄이 아니라 그 종류를 잠근다** — 파라미터 창을 여는 함수가 참조하는
전역 이름이 전부 실제로 있는지 본다. 창을 띄우지 않으므로 빠르고, 창이 열리는 조건
(자료·설정)과 무관하게 돈다.

바이트코드의 `LOAD_GLOBAL` 만 보므로 속성 이름(`dlg.reject` 의 `reject`)은 안 센다.
중첩 함수(`_apply` 같은 것)도 따라 들어간다.
"""
from __future__ import annotations

import builtins
import dis
import importlib
import pkgutil

import pytest

#: 훑을 꾸러미. 창 클래스가 전부 이 아래에 있다.
PACKAGES = ["apex.gui.workflow", "apex.gui.workflow.cmd", "apex.gui.workflow.lc"]
#: 이 이름으로 끝나는 메서드를 본다 — 파라미터 창을 여는 자리들이다.
METHOD_NAMES = ("open_parameters_dialog",)


def _globals_used(code, seen=None) -> set[str]:
    """그 코드와 중첩 함수가 `LOAD_GLOBAL` 로 집는 이름 전부."""
    seen = seen if seen is not None else set()
    for ins in dis.get_instructions(code):
        if ins.opname == "LOAD_GLOBAL" and isinstance(ins.argval, str):
            seen.add(ins.argval)
    for const in code.co_consts:
        if hasattr(const, "co_names"):
            _globals_used(const, seen)
    return seen


def _methods():
    """(모듈 이름, 클래스 이름, 메서드) 목록."""
    out = []
    for pkg_name in PACKAGES:
        try:
            pkg = importlib.import_module(pkg_name)
        except Exception as exc:                      # noqa: BLE001
            pytest.skip(f"{pkg_name} 를 못 불러온다: {exc}")
        for info in pkgutil.iter_modules(pkg.__path__):
            mod_name = f"{pkg_name}.{info.name}"
            try:
                mod = importlib.import_module(mod_name)
            except Exception:                         # noqa: BLE001
                continue                              # 임포트 자체가 안 되는 것은 다른 시험의 몫
            for obj_name in dir(mod):
                obj = getattr(mod, obj_name, None)
                if not isinstance(obj, type):
                    continue
                for meth_name in METHOD_NAMES:
                    meth = obj.__dict__.get(meth_name)
                    if callable(meth):
                        out.append((mod, obj_name, meth_name, meth))
    return out


def test_there_are_parameter_dialog_openers_to_check():
    """훑는 규칙이 망가지면 아래 시험이 조용히 아무것도 안 본다."""
    pytest.importorskip("PyQt5")
    found = _methods()
    assert len(found) >= 5, f"파라미터 창 여는 함수를 {len(found)} 개밖에 못 찾았다"


def test_every_parameter_dialog_opener_resolves_its_names():
    """없는 이름을 집으면 그 단계의 파라미터 창이 앱을 끝낸다."""
    pytest.importorskip("PyQt5")
    missing = []
    for mod, cls_name, meth_name, meth in _methods():
        code = getattr(meth, "__code__", None)
        if code is None:
            continue
        for name in _globals_used(code):
            if name in vars(mod) or hasattr(builtins, name):
                continue
            missing.append(f"{mod.__name__}.{cls_name}.{meth_name} → {name}")
    assert not missing, (
        "파라미터 창을 여는 함수가 없는 이름을 참조한다 — 누르면 앱이 끝난다:\n  "
        + "\n  ".join(sorted(missing)))
