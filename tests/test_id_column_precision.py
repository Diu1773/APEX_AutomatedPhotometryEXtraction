"""A 19-digit Gaia identifier must survive the trip to disk and back.

float64 carries about 15-16 significant digits; a Gaia DR3 identifier needs 19.
Put one in a float column and the last three digits are gone, and ``to_csv``
then writes it in scientific notation. That is how ``8.14859719594028e+17``
ended up in real master_catalog files next to the intact
``814859719594028032`` in the neighbouring column (``Main/FAILURES.md`` F-286).

The defence is per column, not per file: every column that can hold such an id
is listed in ``io_utils.ID_LIKE_COLUMNS`` and passed through Int64 on read, on
construction, and before write.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from apex.utils.io_utils import (
    ID_LIKE_COLUMNS,
    coerce_int64_source_id,
    normalize_id_columns,
    read_csv_int64_source_id,
)

# Three real-shaped ids one apart. float64 collapses all three into one value.
# The first happens to sit on a float64 grid point and round-trips unharmed —
# which is why the damage is invisible in some real files and plain in others
# (25 of 180 rows in one merged YZ Boo catalogue). LOSSY is one that does not
# survive; use it wherever the point is the loss itself.
IDS = [814859719594028032, 814859719594028033, 814859719594028034]
LOSSY = IDS[1]


def test_the_hazard_is_real_so_the_test_means_something():
    """Guard the premise: plain float64 really does destroy these ids."""
    as_float = pd.Series(IDS, dtype="float64")
    assert len(set(as_float.astype("int64").tolist())) < 3


def test_every_id_column_is_named_in_one_place():
    assert "source_id" in ID_LIKE_COLUMNS
    assert "gaia_source_id" in ID_LIKE_COLUMNS
    assert "gaia_id" in ID_LIKE_COLUMNS


@pytest.mark.parametrize("column", ["source_id", "gaia_source_id", "gaia_id"])
def test_ids_survive_a_write_and_read(tmp_path, column):
    df = pd.DataFrame({column: IDS, "ra_deg": [1.0, 2.0, 3.0]})
    path = tmp_path / "catalog.tsv"

    normalize_id_columns(df)
    df.to_csv(path, sep="\t", index=False, na_rep="NaN")

    text = path.read_text(encoding="utf-8")
    assert "e+" not in text.lower(), f"scientific notation in the file:\n{text}"

    back = read_csv_int64_source_id(path, sep="\t")
    assert back[column].astype("int64").tolist() == IDS


def test_a_missing_value_does_not_drag_the_column_into_float(tmp_path):
    """The row that has no Gaia match is what made the column float64."""
    df = pd.DataFrame({
        "gaia_source_id": [IDS[0], pd.NA, IDS[2]],
        "source_id": [IDS[0], -1, IDS[2]],
    })
    path = tmp_path / "catalog.tsv"
    normalize_id_columns(df)
    df.to_csv(path, sep="\t", index=False, na_rep="NaN")

    back = read_csv_int64_source_id(path, sep="\t")
    got = back["gaia_source_id"]
    assert int(got.iloc[0]) == IDS[0]
    assert pd.isna(got.iloc[1])
    assert int(got.iloc[2]) == IDS[2]
    assert back["source_id"].astype("int64").tolist() == [IDS[0], -1, IDS[2]]


def test_np_nan_in_row_dicts_is_the_trap_pd_na_is_not():
    """Why the merger builds its rows with pd.NA. Documents the mechanism."""
    lossy = pd.DataFrame([{"gaia_id": LOSSY}, {"gaia_id": np.nan}])
    assert lossy["gaia_id"].dtype == "float64"
    assert int(lossy["gaia_id"].iloc[0]) != LOSSY         # already rounded

    kept = pd.DataFrame([{"gaia_id": LOSSY}, {"gaia_id": pd.NA}])
    assert kept["gaia_id"].dtype == object
    assert int(coerce_int64_source_id(kept["gaia_id"]).iloc[0]) == LOSSY


def test_a_file_already_written_in_scientific_notation_reads_consistently(tmp_path):
    """Old files cannot be repaired, but they must still read the same way
    every time — the rounded value is what they now mean."""
    path = tmp_path / "old.tsv"
    path.write_text("gaia_source_id\n8.148597195940281e+17\n", encoding="utf-8")

    first = read_csv_int64_source_id(path, sep="\t")["gaia_source_id"].iloc[0]
    second = read_csv_int64_source_id(path, sep="\t")["gaia_source_id"].iloc[0]
    assert int(first) == int(second)
    assert int(first) != LOSSY, "the file no longer holds the exact id"
