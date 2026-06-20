import pandas as pd

from apex.gui.workflow.step1_header_target import select_header_target


def test_select_header_target_prefers_selected_eligible_frame():
    headers = pd.DataFrame([
        {
            "Filename": "first.fit",
            "OBJECT": "First",
            "RA_DEG": 10.0,
            "DEC_DEG": 20.0,
        },
        {
            "Filename": "selected.fit",
            "OBJECT": "Selected",
            "RA_DEG": 30.0,
            "DEC_DEG": -15.0,
        },
    ])

    result = select_header_target(
        headers,
        eligible_filenames=["first.fit", "selected.fit"],
        preferred_filenames=["selected.fit"],
    )

    assert result == {
        "filename": "selected.fit",
        "object_name": "Selected",
        "ra_deg": 30.0,
        "dec_deg": -15.0,
    }


def test_select_header_target_skips_invalid_or_excluded_frames():
    headers = pd.DataFrame([
        {
            "Filename": "invalid.fit",
            "OBJECT": "Invalid",
            "RA_DEG": float("nan"),
            "DEC_DEG": 10.0,
        },
        {
            "Filename": "excluded.fit",
            "OBJECT": "Excluded",
            "RA_DEG": 50.0,
            "DEC_DEG": 10.0,
        },
        {
            "Filename": "usable.fit",
            "OBJECT": "",
            "RA_DEG": 359.9,
            "DEC_DEG": 89.0,
        },
    ])

    result = select_header_target(
        headers,
        eligible_filenames=["invalid.fit", "usable.fit"],
        preferred_filenames=["excluded.fit", "invalid.fit"],
    )

    assert result == {
        "filename": "usable.fit",
        "object_name": "usable",
        "ra_deg": 359.9,
        "dec_deg": 89.0,
    }


def test_select_header_target_returns_none_without_valid_coordinates():
    headers = pd.DataFrame([
        {
            "Filename": "bad.fit",
            "OBJECT": "Bad",
            "RA_DEG": 400.0,
            "DEC_DEG": -100.0,
        },
    ])

    assert select_header_target(headers, ["bad.fit"]) is None
