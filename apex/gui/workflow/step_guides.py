"""Per-step quick guides shown by the "❔ 가이드" button in every step window.

Centralised here so :class:`StepWindowBase` stays small and the 23 step
windows don't each need editing.  Look-up order:

    1. mode-specific  (``(mode, step_index)``)
    2. shared          (``step_index``, used for the common Step 1–7)
    3. generic fallback

Each guide is a short, plain-language "what to do here" — not a manual.
"""

from __future__ import annotations

# Shared steps (identical in CMD and LC): index 0–6.
_SHARED: dict[int, str] = {
    0: (
        "<b>파일 선택</b><br>"
        "• 관측 FITS가 들어있는 폴더를 지정합니다.<br>"
        "• 타겟 이름을 입력하고 <b>Resolve SIMBAD</b>로 좌표를 받아옵니다 "
        "(이후 WCS·카탈로그가 이 좌표를 사용).<br>"
        "• 표에서 분석할 프레임만 체크해 둡니다."
    ),
    1: (
        "<b>이미지 크롭</b> (선택 단계)<br>"
        "• 가장자리 노이즈나 비네팅을 잘라낼 때 사용합니다.<br>"
        "• 크롭이 필요 없으면 그대로 다음 단계로 넘어가도 됩니다."
    ),
    2: (
        "<b>스카이 미리보기 & QC</b><br>"
        "• 배경(sky) 추정과 측광 조리개(aperture)의 크기를 미리 확인합니다.<br>"
        "• FWHM·sky 시그마가 비정상이면 파라미터를 조정하세요.<br>"
        "• QC를 통과한 프레임만 이후 단계에서 사용됩니다."
    ),
    3: (
        "<b>소스 검출</b><br>"
        "• 모든 프레임에서 별을 검출하고 결과를 캐시에 저장합니다.<br>"
        "• <b>Run Detection</b> 후 미리보기에서 검출 갯수·FWHM을 확인하세요.<br>"
        "• 검출이 너무 적으면 detection sigma를 낮춰 보세요."
    ),
    4: (
        "<b>WCS 플레이트 솔빙</b><br>"
        "• 각 프레임에 천구 좌표(WCS)를 부여합니다.<br>"
        "• <b>Internal (Python)</b> 탭은 외부 프로그램 없이 동작합니다 "
        "(Gaia 다운로드 필요).<br>"
        "• ASTAP/Astrometry.net이 설치돼 있으면 해당 탭이 더 빠릅니다.<br>"
        "• 솔버를 바꿔 다시 돌려도 성공한 프레임만 갱신됩니다."
    ),
    5: (
        "<b>마스터 카탈로그 생성</b><br>"
        "• WCS가 풀린 프레임들을 합쳐 고정된 별 목록(master)을 만듭니다.<br>"
        "• 이후 모든 측광이 이 마스터 ID를 기준으로 정렬됩니다.<br>"
        "• WCS 매칭이 좋은 프레임이 자동으로 기준 프레임이 됩니다."
    ),
    6: (
        "<b>강제 조리개 측광</b><br>"
        "• 마스터 별 위치를 각 프레임에 투영해 측광합니다.<br>"
        "• 성장곡선(growth curve)으로 조리개 보정을 계산합니다.<br>"
        "• Stats 탭에서 센터링 오차·아웃라이어를 점검하세요."
    ),
}

# CMD-mode steps after the shared ones (index 7–11).
_CMD: dict[int, str] = {
    7: (
        "<b>PSF 측광</b> (선택 단계)<br>"
        "• 밀집한 성단 코어에서 더 정밀한 측광이 필요할 때 사용합니다.<br>"
        "• 필요 없으면 <b>Skip PSF</b>로 넘어가도 됩니다 "
        "(이후 단계는 Step 7 강제 측광 결과를 사용)."
    ),
    8: (
        "<b>마스터 ID 편집기</b><br>"
        "• 별을 클릭해 정보를 확인하고, 잘못 잡힌 별을 정리합니다.<br>"
        "• <b>Members ●</b>를 켜면 Gaia 멤버십 확률이 높은 별에 점이 찍힙니다.<br>"
        "• 필요하면 ROI(관심 영역) 원을 설정해 성단 영역만 남깁니다."
    ),
    9: (
        "<b>영점 보정 (Zeropoint)</b><br>"
        "• 기기 등급을 Gaia/표준 등급으로 변환하는 영점을 맞춥니다.<br>"
        "• 색-등급 관계의 기울기와 잔차를 확인하세요.<br>"
        "• 결과가 다음 CMD/등급 계산의 기준이 됩니다."
    ),
    10: (
        "<b>CMD 플롯</b><br>"
        "• 보정된 등급으로 색-등급도(CMD)를 그립니다.<br>"
        "• 좌측은 기기(Instrumental), 우측은 표준화 등급입니다.<br>"
        "• Membership·SNR·시차(parallax) 필터로 성단 멤버만 강조할 수 있습니다."
    ),
    11: (
        "<b>아이소크론 모델</b><br>"
        "• PARSEC 아이소크론을 불러와 성단의 나이·거리·금속함량을 맞춥니다.<br>"
        "• <b>CMD Viewer</b> 탭에서 슬라이더로 곡선을 데이터에 대략 맞춘 뒤,<br>"
        "&nbsp;&nbsp;<b>Auto Fit</b> 탭에서 그 값을 seed로 가져와 정밀 피팅하세요.<br>"
        "• 첫 실행은 아이소크론 캐시를 만드느라 잠깐 걸립니다(이후 빠름)."
    ),
}

# LC-mode steps after the shared ones (index 7–10).
_LC: dict[int, str] = {
    7: (
        "<b>타겟 / 비교성 선택</b><br>"
        "• 변광을 측정할 타겟별과, 밝기 기준이 될 비교성을 고릅니다.<br>"
        "• 비교성은 변광이 없고 타겟과 밝기가 비슷한 별이 좋습니다.<br>"
        "• 이미지에서 별을 클릭하거나 표에서 역할(Target/Comp)을 지정하세요."
    ),
    8: (
        "<b>라이트커브 생성</b><br>"
        "• 선택한 비교성으로 타겟의 차등 등급(differential)을 만듭니다.<br>"
        "• 비교성 QC로 흔들리는 비교성을 자동 제외할 수 있습니다.<br>"
        "• 여러 폴더(밤)를 합쳐서 사용할 수도 있습니다."
    ),
    9: (
        "<b>디트렌드 & 밤 병합</b><br>"
        "• 밤마다 다른 영점·대기 소광을 보정해 라이트커브를 이어 붙입니다.<br>"
        "• 보정 모드(offset/global)를 고르고 <b>Fit &amp; Apply</b>로 적용하세요.<br>"
        "• RMS 전/후 값으로 보정 효과를 확인합니다."
    ),
    10: (
        "<b>주기 분석</b><br>"
        "• Lomb-Scargle·PDM·BLS로 변광 주기를 찾습니다.<br>"
        "• 여러 방법의 결과가 일치하면 그 주기가 신뢰도 높습니다.<br>"
        "• Phase Plot 탭에서 찾은 주기로 위상 접기(phase-fold)를 확인하세요.<br>"
        "• 정밀 분석(O-C, 트랜짓 피팅)은 Tools 메뉴에 있습니다."
    ),
}

_FALLBACK = "이 단계에 대한 가이드가 아직 준비되지 않았습니다."


def step_guide(mode: str, step_index: int, step_name: str = "") -> str:
    """Return the HTML guide string for a given mode/step."""
    table = _CMD if mode == "cmd" else _LC
    if step_index in table:
        return table[step_index]
    if step_index in _SHARED:
        return _SHARED[step_index]
    return _FALLBACK
