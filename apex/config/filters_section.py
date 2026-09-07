"""설정의 ``[filters]`` 표를 읽어 필터 별칭으로 만든다.

**왜 설정이어야 하나.** 관측소마다 필터를 자기 방식으로 적는다. LCO 는 g' 를
``gp`` 로, 어떤 곳은 ``SDSS-G``, 어떤 곳은 ``Bessell-V``, 어떤 곳은 그냥 번호를
쓴다. APEX 의 기본 별칭표는 흔한 철자를 모아 둔 것일 뿐이고, **세상의 모든
이름을 담으려 드는 것 자체가 틀린 설계다** — 새 관측소를 만날 때마다 우리 코드를
고쳐야 한다는 뜻이 되기 때문이다.

그래서 설정에서 얼마든지 더할 수 있게 한다.

    {
      "filters": {
        "aliases": {
          "Bessell-V": "V",
          "zs": "z",
          "F01": "B"
        }
      }
    }

**사용자 별칭이 기본표를 이긴다.** 어떤 관측소가 우리 표와 다른 뜻으로 이름을
쓴다면 우리가 아니라 그쪽이 옳다.

**오른쪽(정본 키)은 APEX 의 규약을 따른다** — Johnson 은 대문자(B, V, R, I),
SDSS 는 소문자(g, r, i, z), 협대역은 첫 글자만 대문자(Ha, OIII).
"""

from __future__ import annotations

from typing import Any, Dict, Mapping, Optional


def read_filters_section(data: Optional[Mapping[str, Any]]) -> Dict[str, str]:
    """``[filters].aliases`` 를 ``{원래 이름: 정본 키}`` 로 돌려준다.

    ``[filters]`` 아래에 ``aliases`` 표가 있으면 그것을 쓰고, 없으면
    ``[filters]`` 자체를 별칭표로 본다(짧게 쓰고 싶은 사람을 위해).
    빈 값은 버린다 — 설정에서 지우다 만 줄이 필터 이름을 빈칸으로 만들면 안 된다.
    """
    if not isinstance(data, Mapping):
        return {}
    section = data.get("filters")
    if not isinstance(section, Mapping):
        return {}

    table = section.get("aliases")
    if not isinstance(table, Mapping):
        table = {k: v for k, v in section.items() if not isinstance(v, Mapping)}

    out: Dict[str, str] = {}
    for key, value in table.items():
        k, v = str(key).strip(), str(value).strip()
        if k and v:
            out[k] = v
    return out


def filters_toml_section(aliases: Mapping[str, str]) -> Dict[tuple, Dict[str, Any]]:
    """설정을 다시 쓸 때 쓰는 역변환. 비어 있으면 절을 만들지 않는다."""
    clean = {str(k).strip(): str(v).strip()
             for k, v in dict(aliases or {}).items()
             if str(k).strip() and str(v).strip()}
    return {("filters", "aliases"): clean} if clean else {}
