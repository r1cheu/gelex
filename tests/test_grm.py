from pathlib import Path

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


def test_make_grm_chunking_consistency(test_bed_file):
    """Test that chunked and non-chunked GRM computation gives same result"""
    # Compute GRM with chunking
    grm_chunked = make_grm(
        test_bed_file, method="add", chunk_size=2, save=False
    ).to_numpy()

    # Compute GRM without chunking
    grm_no_chunk = make_grm(
        test_bed_file, method="add", chunk_size=1000, save=False
    ).to_numpy()

    assert np.allclose(grm_chunked, grm_no_chunk)


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
