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
        np.eye(3),
        index=["s1", "s2", "s3"],
        columns=["s1", "s2", "s3"],
    )


@pytest.fixture
def grm_b():
    return pd.DataFrame(
        np.eye(3),
        index=["s1", "s2", "s3"],
        columns=["s1", "s2", "s3"],
    )


@pytest.fixture
def grm_with_less_individuals():
    return pd.DataFrame(np.eye(2), index=["s1", "s2"], columns=["s1", "s2"])


@pytest.fixture
def grm_with_more_individuals():
    return pd.DataFrame(
        np.eye(4),
        index=["s1", "s2", "s3", "s4"],  # 's4' causes mismatch
        columns=["s1", "s2", "s3", "s4"],
    )


def test_with_less_individuals_grm(sample_data, grm_with_less_individuals):
    model = make_model(sample_data)

    with pytest.raises(
        ValueError,
        match="Some individuals in the `y` are not present in the GRM. Are you sure you are using the correct GRM?",
    ):
        model.make("y ~ x1", {"grm1": grm_with_less_individuals})


def test_with_more_individuals(sample_data, grm_with_more_individuals):
    model_make = make_model(sample_data)
    model = model_make.make("y ~ x1", {"grm1": grm_with_more_individuals})
    assert model._dropped_individuals == ["s4"]


def test_initialization_nan_in_response(sample_data, grm_a):
    sample_data.iloc[0, 0] = np.nan
    model = make_model(sample_data)
    model = model.make("y~x1", {"a": grm_a})
    assert model._dropped_individuals == ["s1"]


def test_initialization_with_dataframe(sample_data, grm_a):
    model = make_model(sample_data)
    model.make("y ~ 1", {"a": grm_a})
    pd.testing.assert_frame_equal(model.data, sample_data)


def test_initialization_with_series(sample_data, grm_a):
    model = make_model(sample_data["y"])
    model.make("y ~ 1", {"a": grm_a})
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
