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
def mismatch_grm():
    return pd.DataFrame(
        {"s2": [1.0, 0.1, 0.2], "s3": [0.1, 1.0, 0.3], "s4": [0.2, 0.3, 1.0]},
        index=["s2", "s3", "s4"],
    )


def test_initialization_with_dataframe(sample_data):
    model = make_model(sample_data)
    pd.testing.assert_frame_equal(model.data, sample_data)


def test_initialization_with_file(tmp_path, sample_data):
    # Create a temporary file
    file_path = tmp_path / "test_data.csv"
    sample_data.to_csv(file_path, sep="\t", index=True, header=True)
    model = make_model(file_path)
    pd.testing.assert_frame_equal(model.data, sample_data)


def test_initialization_with_categorical(sample_data):
    model = make_model(sample_data, categorical="x2")
    assert isinstance(model.data["x2"].dtype, pd.CategoricalDtype)


def test_make_method(sample_data, grm_a):
    model = make_model(sample_data)
    lmm = model.make("y ~ x1", {"grm1": grm_a})

    assert lmm.n_samples == 3
    assert lmm.n_fixed_effect == 2
    assert lmm.n_random_effect == 2


def test_make_method_missing_values(sample_data, grm_a):
    # Add missing value
    sample_data.loc["s1", "x1"] = np.nan
    model = make_model(sample_data)

    with pytest.raises(
        ValueError, match="fixed effects should not contain missing values"
    ):
        model.make("y ~ x1", {"grm1": grm_a})


def test_clean_data(sample_data):
    model = make_model(sample_data)
    cleaned = model._clean_data("y")
    pd.testing.assert_frame_equal(cleaned, sample_data)  # no NAs in sample data

    # Add NA and test
    sample_data.loc["s1", "y"] = np.nan
    model = make_model(sample_data)
    cleaned = model._clean_data("y")
    assert len(cleaned) == 2  # should drop s1


def test_load_grm(sample_data, grm_a):
    model = make_model(sample_data)
    data, grm_cube, names = model._load_grm({"grm1": grm_a}, sample_data)
    assert grm_cube.shape == (3, 3, 1)
    assert names == ["grm1"]
    pd.testing.assert_frame_equal(data, sample_data, check_names=False)


def test_load_grm_mismatch(sample_data, mismatch_grm):
    # Create GRM with different samples

    model = make_model(sample_data)
    data, _, _ = model._load_grm({"grm1": mismatch_grm}, sample_data)
    sample_data = sample_data.drop("s1")
    pd.testing.assert_frame_equal(data, sample_data, check_names=False)
