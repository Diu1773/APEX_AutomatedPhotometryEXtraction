from pathlib import Path

import pandas as pd

from apex.utils.photometry_loader import load_frame_photometry
from apex.utils.step_paths import step7_forced_phot_dir
from apex.utils.step_paths_lc import step8_selection_dir


def _write_catalog(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_load_frame_photometry_preserves_johnson_filter_case_for_id_map(tmp_path):
    result_dir = tmp_path
    forced_dir = step7_forced_phot_dir(result_dir)
    selection_dir = step8_selection_dir(result_dir)
    forced_dir.mkdir(parents=True)
    selection_dir.mkdir(parents=True)

    _write_catalog(
        selection_dir / "master_catalog_B.tsv",
        [{"ID": 7, "source_id": 123456789012345678}],
    )
    _write_catalog(
        selection_dir / "master_catalog_V.tsv",
        [{"ID": 42, "source_id": 123456789012345678}],
    )
    pd.DataFrame(
        [{
            "file": "frame_V.fits",
            "FILTER": "V",
            "source_id": 123456789012345678,
            "ID": 999,
            "mag": 14.2,
        }]
    ).to_csv(forced_dir / "photometry_frame_V.fits.tsv", sep="\t", index=False)

    df = load_frame_photometry(result_dir, "frame_V.fits")

    assert df is not None
    assert int(df.loc[0, "ID"]) == 42


def test_load_frame_photometry_exposes_mag_inst_as_common_mag(tmp_path):
    forced_dir = step7_forced_phot_dir(tmp_path)
    forced_dir.mkdir(parents=True)
    pd.DataFrame(
        [{"file": "frame.fits", "ID": 3, "source_id": -3, "mag_inst": 15.25}]
    ).to_csv(forced_dir / "photometry_frame.fits.tsv", sep="\t", index=False)

    df = load_frame_photometry(tmp_path, "frame.fits")

    assert df is not None
    assert df.loc[0, "mag"] == 15.25
