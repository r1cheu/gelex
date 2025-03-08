from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from phenx import load_grm, make_grm


@pytest.fixture
def test_bed_file():
    """Fixture for test BED file path"""
    return Path(__file__).resolve().parent / "data/test.bed"


def test_make_grm_additive(test_bed_file):
    """Test additive GRM computation"""
    grm = make_grm(test_bed_file, method="add", chunk_size=2, save=False)

    assert isinstance(grm, pd.DataFrame)
    assert grm.shape == (3, 3)  # 3 samples
    assert np.allclose(grm.to_numpy(), grm.to_numpy().T)  # Symmetric
    assert np.all(np.diag(grm.to_numpy()) >= 0)  # Positive diagonal


def test_make_grm_dominance(test_bed_file):
    """Test dominance GRM computation"""
    grm = make_grm(test_bed_file, method="dom", chunk_size=2, save=False)

    assert isinstance(grm, pd.DataFrame)
    assert grm.shape == (3, 3)
    assert np.allclose(grm.to_numpy(), grm.to_numpy().T)
    assert np.all(np.diag(grm.values) >= 0)


def test_make_grm_invalid_method(test_bed_file):
    """Test invalid method raises ValueError"""
    with pytest.raises(ValueError):
        make_grm(test_bed_file, method="invalid")


def test_make_grm_large_chunk_size(test_bed_file):
    """Test with chunk size larger than number of SNPs"""
    grm = make_grm(test_bed_file, method="add", chunk_size=1000, save=False)
    assert isinstance(grm, pd.DataFrame)
    assert grm.shape == (3, 3)


def test_make_grm_save_file(test_bed_file):
    """Test GRM saving to file"""
    output_path: Path = test_bed_file.with_suffix(".add.grm")
    grm = make_grm(test_bed_file, method="add", save=True)
    assert grm.to_numpy().flags["F_CONTIGUOUS"]
    assert output_path.exists()
    output_path.unlink()


def test_load_grm(test_bed_file):
    """Test loading GRM from file, make sure f-contiguous"""
    output_path: Path = test_bed_file.with_suffix(".add.grm")
    _ = make_grm(test_bed_file, method="add", save=True)
    loaded_grm = load_grm(output_path, True)
    assert loaded_grm.flags["F_CONTIGUOUS"]
    output_path.unlink()


@pytest.fixture
def mock_grm_file(tmp_path):
    """Create a mock GRM file for testing."""
    grm_path = tmp_path / "test_grm.h5"

    # Sample data
    mock_grm = np.array([[1.0, 0.5, 0.2], [0.5, 1.0, 0.3], [0.2, 0.3, 1.0]])
    individuals = ["ind1", "ind2", "ind3"]

    # Create the test H5 file
    with h5py.File(grm_path, "w") as f:
        f.create_dataset("grm", data=mock_grm)
        f.create_dataset("individuals", data=individuals)

    return grm_path


def test_load_grm_dataframe(mock_grm_file):
    """Test loading GRM as a pandas DataFrame."""
    result = load_grm(mock_grm_file)

    assert isinstance(result, pd.DataFrame)
    assert result.shape == (3, 3)
    assert list(result.index) == ["ind1", "ind2", "ind3"]
    assert list(result.columns) == ["ind1", "ind2", "ind3"]
    assert result.loc["ind1", "ind2"] == 0.5


def test_load_grm_array(mock_grm_file):
    """Test loading GRM as a NumPy array."""
    result = load_grm(mock_grm_file, return_array=True)

    assert isinstance(result, np.ndarray)
    assert result.shape == (3, 3)
    assert result[0, 1] == 0.5


def test_load_grm_with_dropped_individual(mock_grm_file):
    """Test loading GRM with dropped individuals."""
    result = load_grm(mock_grm_file, dropped_individual=["ind2"])

    assert isinstance(result, pd.DataFrame)
    assert result.shape == (2, 2)
    assert "ind2" not in result.index
    assert "ind2" not in result.columns
    assert list(result.index) == ["ind1", "ind3"]


def test_load_grm_file_not_found():
    """Test that FileNotFoundError is raised for non-existent file."""
    with pytest.raises(FileNotFoundError):
        load_grm("non_existent_file.h5")


def test_load_grm_fcontig(mock_grm_file):
    result = load_grm(mock_grm_file, return_array=True)
    assert result.flags["F_CONTIGUOUS"]

    result = load_grm(mock_grm_file)
    assert result.to_numpy().flags["F_CONTIGUOUS"]
