"""
This module contains unit tests for the phenx model initialization and data cleaning functionality.

The tests cover:
- Model initialization with DataFrame or file input
- Handling of categorical variables
- Formula parsing and model creation
- Data cleaning operations including missing value handling
- Genomic relationship matrix (GRM) validation and alignment

Key test scenarios include:
- Basic initialization with sample data
- File-based initialization using temporary CSV files
- Categorical variable type checking
- Fixed and random effect specification
- Missing value detection and handling
- GRM matching and sample alignment
"""

import numpy as np
import pandas as pd
import pytest
from phenx import make_model


@pytest.fixture
def sample_data():
    return pd.DataFrame(
        {
            "y": [1.0, 2.0, 3.0],
            "x1": [0.5, 1.5, 2.5],
            "x2": ["a", "b", "a"],
            "sample_id": ["s1", "s2", "s3"],
        }
    ).set_index("sample_id")


@pytest.fixture
def grm_a():
    return pd.DataFrame(
        {"s1": [1.0, 0.1, 0.2], "s2": [0.1, 1.0, 0.3], "s3": [0.2, 0.3, 1.0]},
        index=["s1", "s2", "s3"],
    )


@pytest.fixture
def grm_b():
    return pd.DataFrame(
        {"s1": [1.0, 0.1, 0.2], "s2": [0.1, 1.0, 0.3], "s3": [0.2, 0.3, 1.0]},
        index=["s1", "s2", "s3"],
    )


@pytest.fixture
def grm_with_size_mismatch():
    return pd.DataFrame(np.eye(2), index=["s1", "s2"], columns=["s1", "s2"])


@pytest.fixture
def grm_with_index_mismatch():
    return pd.DataFrame(
        np.eye(3),
        index=["s1", "s2", "s4"],  # 's4' causes mismatch
        columns=["s1", "s2", "s4"],
    )


def test_sample_size_mismatch(sample_data, grm_with_size_mismatch):
    model = make_model(sample_data)

    with pytest.raises(ValueError, match="Sample size mismatch for 'grm1' and 'y'."):
        model.make("y ~ x1", {"grm1": grm_with_size_mismatch})


def test_index_mismatch(sample_data, grm_with_index_mismatch):
    model = make_model(sample_data)

    with pytest.raises(ValueError, match="GRM sample indices do not align with 'y'."):
        model.make("y ~ x1", {"grm1": grm_with_index_mismatch})


def test_initialization_nan_in_response(sample_data, grm_a):
    sample_data.iloc[0, 0] = np.nan
    model = make_model(sample_data)
    with pytest.raises(ValueError, match="`y` contains missing value"):
        model.make("y~x1", grm={"a": grm_a})


def test_initialization_with_dataframe(sample_data, grm_a):
    model = make_model(sample_data)
    model.make("y ~ 1", grm={"a": grm_a})
    pd.testing.assert_frame_equal(model.data, sample_data)


def test_initialization_with_series(sample_data, grm_a):
    model = make_model(sample_data["y"])
    model.make("y ~ 1", grm={"a": grm_a})
    pd.testing.assert_frame_equal(model.data, pd.DataFrame(sample_data["y"]))


def test_initialization_with_file(tmp_path, sample_data):
    # Create a temporary file
    file_path = tmp_path / "test_data.csv"
    sample_data.to_csv(file_path, sep="\t", index=True, header=True)
    model = make_model(file_path)
    pd.testing.assert_frame_equal(model.data, sample_data)


def test_initialization_with_categorical(sample_data, grm_a):
    model = make_model(sample_data, categorical="x2")
    lmm = model.make("y ~ x1 + x2", {"grm1": grm_a})
    assert isinstance(model.data["x2"].dtype, pd.CategoricalDtype)
    assert lmm.num_fixed_effects == 3


def test_make_method_missing_values(sample_data, grm_a):
    sample_data.loc["s1", "x1"] = np.nan
    model = make_model(sample_data)

    with pytest.raises(ValueError, match="`x1` contains missing value"):
        model.make("y ~ x1", {"grm1": grm_a})


def test_make_method(sample_data, grm_a):
    model = make_model(sample_data)
    lmm = model.make("y ~ x1", {"grm1": grm_a})

    assert lmm.num_individuals == 3
    assert lmm.num_fixed_effects == 2
    assert lmm.num_random_effects == 2
