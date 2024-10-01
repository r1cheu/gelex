import numpy as np
import pandas as pd
import pytest
from phenx.dataset import _hybird_encode


@pytest.fixture
def sample_genotype():
    return pd.DataFrame(
        {"SNP1": [0, 1, 2, 0], "SNP2": [1, 0, 1, 2], "SNP3": [2, 2, 0, 1]},
        index=["sample1", "sample2", "sample3", "sample4"],
    )


@pytest.fixture
def sample_phenotype():
    return pd.Series(
        [0.1, 0.2, 0.3, 0.4], index=["sample1", "sample2", "sample3", "sample4"]
    )


@pytest.fixture
def sample_d_values():
    return pd.DataFrame(
        {"SNP1": [2, 0.2], "SNP2": [0, 0.4], "SNP3": [0, 0.6]},
        index=["value1", "value2"],
    )


def test_hybird_encode_with_phenotype(sample_genotype, sample_phenotype):
    result = _hybird_encode(sample_genotype, sample_phenotype)
    assert isinstance(result, np.ndarray)
    assert result.shape == (4, 3)


def test_hybird_encode_with_d_values_file(sample_genotype, sample_d_values, tmp_path):
    d_values_file = tmp_path / "d_values.h5"
    sample_d_values.to_hdf(d_values_file, key="d_values")
    result = _hybird_encode(sample_genotype, str(d_values_file))
    assert isinstance(result, np.ndarray)
    assert result.shape == (4, 3)


def test_hybird_encode_invalid_input(sample_genotype):
    with pytest.raises(TypeError):
        _hybird_encode(sample_genotype, 123)


def test_hybird_encode_save_d_values(sample_genotype, sample_phenotype, tmp_path):
    save_path = tmp_path / "d_values_saved.h5"
    _hybird_encode(sample_genotype, sample_phenotype, save_path=str(save_path))
    assert save_path.exists()
    saved_d_values = pd.read_hdf(save_path, key="d_values")
    assert saved_d_values.shape == (2, 3)
