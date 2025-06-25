import pandas as pd
import pytest
from gelexy.utils import Aligner


@pytest.fixture
def master_row_ids():
    return ["S3", "S1", "S5"]


@pytest.fixture
def master_col_ids():
    return ["yield", "height", "moisture"]


@pytest.fixture
def df_messy():
    data = {
        "fat": [3.9, 4.1, 3.5, 4.2],
        "yield": [112, 105, 100, 108],
        "height": [122, 115, 120, 118],
        "temp": [37, 38, 36, 39],
    }
    df = pd.DataFrame(data, index=["S5", "S1", "S99", "S3"])
    df.index.name = "id"
    return df


def test_align_rows_only(df_messy, master_row_ids):
    aligner = Aligner(index_ids=master_row_ids)
    aligned = aligner.align(df_messy, axis=0)
    assert aligned.index.tolist() == master_row_ids


def test_align_columns_only(df_messy, master_col_ids):
    aligner = Aligner(column_ids=master_col_ids)
    aligned = aligner.align(df_messy, axis=1)
    assert aligned.columns.tolist() == ["yield", "height"]


def test_align_both_axes(df_messy, master_row_ids, master_col_ids):
    aligner = Aligner(index_ids=master_row_ids, column_ids=master_col_ids)
    aligned = aligner.align(df_messy, axis=[0, 1])
    assert aligned.index.tolist() == master_row_ids
    assert aligned.columns.tolist() == ["yield", "height"]
