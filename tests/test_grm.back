import numpy as np
import pandas as pd
import pytest
from bed_reader import sample_file
from phenx.dataset import grm

rng = np.random.default_rng(42)


@pytest.mark.filterwarnings("ignore::ResourceWarning")
def test_grm_with_bed_file():
    result = grm(str(sample_file("small.bed")))
    assert isinstance(result, pd.DataFrame)


def test_grm_with_h5_file(tmp_path):
    # Create a temporary .h5 file for testing
    h5_path = tmp_path / "test.h5"
    # Write test data to the .h5 file
    # (Note: This is a placeholder. Actual .h5 file creation depends on the format)
    data = pd.DataFrame(rng.normal(size=(10, 10)))
    data.to_hdf(h5_path, key="grm")

    result = grm(str(h5_path))
    assert isinstance(result, pd.DataFrame)


def test_grm_with_dataframe():
    # Test with a pandas DataFrame
    genotype_df = pd.DataFrame(rng.normal(size=(10, 10)))
    result = grm(genotype_df)
    assert isinstance(result, pd.DataFrame)


def test_grm_with_invalid_input():
    # Test with an invalid input type
    with pytest.raises(FileNotFoundError):
        grm("invalid_path.txt")


def test_grm_with_hybird_encoding(tmp_path):
    # Test with hybird encoding method
    genotype_df = pd.DataFrame()
    phenotype_series = pd.Series(rng.normal(size=(10,)))
    save_d_values = tmp_path / "test_d_values.h5"

    result = grm(
        genotype_df,
        phenotype=phenotype_series,
        encode_method="hybird",
        save_d_values=str(save_d_values),
    )
    assert isinstance(result, pd.DataFrame)
    assert save_d_values.exists()


def test_grm_with_load_hybird_encoding(tmp_path):
    # Test with hybird encoding method
    genotype_df = pd.DataFrame(rng.normal(size=(10, 10)))
    phenotype_series = pd.Series(rng.normal(size=(10,)))
    save_d_values = tmp_path / "test_d_values.h5"

    result = grm(
        genotype_df,
        phenotype=phenotype_series,
        encode_method="hybird",
        save_d_values=str(save_d_values),
    )

    expected = grm(
        genotype_df,
        phenotype=save_d_values,
        encode_method="hybird",
    )
    assert isinstance(expected, pd.DataFrame)
    assert np.array_equal(result.to_numpy(), expected.to_numpy())


def test_grm_save_grm(tmp_path):
    # Test saving the GRM to a file
    genotype_df = pd.DataFrame(rng.normal(size=(10, 10)))
    save_grm_path = tmp_path / "test_grm.h5"

    result = grm(genotype_df, save_grm=str(save_grm_path))
    assert isinstance(result, pd.DataFrame)
    assert save_grm_path.exists()
