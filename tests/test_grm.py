import numpy as np
import pandas as pd
import pytest
from bed_reader import to_bed
from phenx import _load_grm, make_grm


@pytest.fixture
def test_bed_file(tmp_path):
    """Create a simple test BED file"""
    bed_path = tmp_path / "test.bed"
    val = np.array(
        [[1.0, 0.0, np.nan, 0.0], [2.0, 0.0, np.nan, 2.0], [0.0, 1.0, 2.0, 0.0]],
        order="F",
    )
    properties = {
        "fid": ["fid1", "fid1", "fid2"],
        "iid": ["iid1", "iid2", "iid3"],
        "father": ["iid23", "iid23", "iid22"],
        "mother": ["iid34", "iid34", "iid33"],
        "sex": [1, 2, 0],
        "pheno": ["red", "red", "blue"],
        "chromosome": ["1", "1", "5", "Y"],
        "sid": ["sid1", "sid2", "sid3", "sid4"],
        "cm_position": [100.4, 2000.5, 4000.7, 7000.9],
        "bp_position": [1, 100, 1000, 1004],
        "allele_1": ["A", "T", "A", "T"],
        "allele_2": ["A", "C", "C", "G"],
    }
    to_bed(bed_path, val, properties=properties)
    return bed_path


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


def test_make_grm_no_chunking(test_bed_file):
    """Test GRM computation without chunking"""
    grm_add = make_grm(test_bed_file, method="add", chunk_size=False, save=False)
    assert isinstance(grm_add, pd.DataFrame)
    assert grm_add.shape == (3, 3)

    grm_dom = make_grm(test_bed_file, method="dom", chunk_size=False, save=False)
    assert isinstance(grm_dom, pd.DataFrame)
    assert grm_dom.shape == (3, 3)


def test_make_grm_save_file(test_bed_file):
    """Test GRM saving to file"""
    output_path = test_bed_file.with_suffix(".add.grm")
    grm = make_grm(test_bed_file, method="add", save=True)
    assert grm.to_numpy().flags["F_CONTIGUOUS"]

    assert output_path.exists()
    loaded_grm = _load_grm(output_path)
    assert loaded_grm.flags["F_CONTIGUOUS"]
    assert np.allclose(grm.to_numpy(), loaded_grm)

    slice_grm = _load_grm(output_path)
    assert slice_grm.flags["F_CONTIGUOUS"]

    slice_grm = _load_grm(
        output_path, samples_col=["iid1", "iid2"], samples_row=["iid1", "iid2"]
    )
    assert slice_grm.flags["F_CONTIGUOUS"]

    slice_grm = _load_grm(output_path, samples_col=["iid1", "iid2"])
    assert slice_grm.flags["F_CONTIGUOUS"]

    slice_grm = _load_grm(output_path, samples_row=["iid1", "iid2"])
    assert slice_grm.flags["F_CONTIGUOUS"]


def test_make_grm_large_chunk_size(test_bed_file):
    """Test with chunk size larger than number of SNPs"""
    grm = make_grm(test_bed_file, method="add", chunk_size=100000, save=False)
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
        test_bed_file, method="add", chunk_size=False, save=False
    ).to_numpy()

    assert np.allclose(grm_chunked, grm_no_chunk)
