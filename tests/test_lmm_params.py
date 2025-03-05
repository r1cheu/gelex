from pathlib import Path

import h5py
import numpy as np
import pytest
from phenx.model import LinearMixedModelParams
from phenx.model.io import load_params


@pytest.fixture
def sample_params():
    # Create sample data
    beta = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    sigma = np.array([0.5, 0.3, 0.2], dtype=np.float64)
    dropped_individuals = ["dropped1", "dropped2"]

    return {
        "beta": beta,
        "sigma": sigma,
        "dropped_individuals": dropped_individuals,
    }


def test_linear_mixed_model_params_init(sample_params):
    # Test initialization
    params = LinearMixedModelParams(
        beta=sample_params["beta"],
        sigma=sample_params["sigma"],
        dropped_individuals=sample_params["dropped_individuals"],
    )

    assert params is not None


def test_linear_mixed_model_params_properties(sample_params):
    # Create instance
    params = LinearMixedModelParams(
        beta=sample_params["beta"],
        sigma=sample_params["sigma"],
        dropped_individuals=sample_params["dropped_individuals"],
    )

    # Test beta property
    np.testing.assert_array_almost_equal(params.beta, sample_params["beta"])

    # Test sigma property
    np.testing.assert_array_almost_equal(params.sigma, sample_params["sigma"])

    # Test dropped_individuals property
    assert params.dropped_individuals == sample_params["dropped_individuals"]


def test_linear_mixed_model_params_invalid_inputs():
    with pytest.raises(TypeError):
        LinearMixedModelParams(
            beta=[1, 2, 3],  # should be numpy array
            sigma=np.array([0.5, 0.3, 0.2]),
            dropped_individuals=["dropped1"],
        )

    with pytest.raises(TypeError):
        LinearMixedModelParams(
            beta=np.array([1.0, 2.0, 3.0]),
            sigma=[0.5, 0.3, 0.2],  # should be numpy array
            dropped_individuals=["dropped1"],
        )


def test_linear_mixed_model_params_empty_arrays():
    # Test with empty arrays
    params = LinearMixedModelParams(
        beta=np.array([]), sigma=np.array([]), dropped_individuals=[]
    )

    assert len(params.beta) == 0
    assert len(params.sigma) == 0
    assert len(params.dropped_individuals) == 0


def test_load_params(tmp_path, sample_params):
    # Create a temporary HDF5 file with test data
    test_file = tmp_path / "test_params.h5"

    with h5py.File(test_file, "w") as f:
        f.create_dataset("beta", data=sample_params["beta"])
        f.create_dataset("sigma", data=sample_params["sigma"])
        f.create_dataset(
            "dropped_individuals", data=sample_params["dropped_individuals"]
        )

    params = load_params(test_file)

    np.testing.assert_array_almost_equal(params.beta, sample_params["beta"])
    np.testing.assert_array_almost_equal(params.sigma, sample_params["sigma"])

    assert params.dropped_individuals == sample_params["dropped_individuals"]

    params = load_params(Path(test_file))
    np.testing.assert_array_almost_equal(params.beta, sample_params["beta"])
