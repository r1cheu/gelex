import h5py
import numpy as np
import pytest
from gelexy.model import LinearMixedModelParams
from gelexy.model.io import load_params


@pytest.fixture
def sample_params():
    # Create sample data
    beta = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    sigma = np.array([0.5, 0.3, 0.2], dtype=np.float64)
    proj_y = np.array([0.5, 0.3, 0.2, 1.0, 2.0], dtype=np.float64)
    dropped_individuals = ["dropped1", "dropped2"]

    return {
        "beta": beta,
        "sigma": sigma,
        "proj_y": proj_y,
        "dropped_individuals": dropped_individuals,
    }


def test_linear_mixed_model_params_init(sample_params):
    # Test initialization
    params = LinearMixedModelParams(
        beta=sample_params["beta"],
        sigma=sample_params["sigma"],
        proj_y=sample_params["proj_y"],
        dropped_individuals=sample_params["dropped_individuals"],
    )
    assert params is not None


def test_linear_mixed_model_params_properties(sample_params):
    # Create instance
    params = LinearMixedModelParams(
        beta=sample_params["beta"],
        sigma=sample_params["sigma"],
        proj_y=sample_params["proj_y"],
        dropped_individuals=sample_params["dropped_individuals"],
    )

    np.testing.assert_array_almost_equal(params.beta, sample_params["beta"])
    np.testing.assert_array_almost_equal(params.sigma, sample_params["sigma"])
    np.testing.assert_array_almost_equal(params.proj_y, sample_params["proj_y"])

    # Test dropped_individuals property
    assert params.dropped_individuals == sample_params["dropped_individuals"]


def test_load_params(tmp_path, sample_params):
    # Create a temporary HDF5 file with test data
    test_file = tmp_path / "test_params.h5"

    with h5py.File(test_file, "w") as f:
        f.create_dataset("beta", data=sample_params["beta"])
        f.create_dataset("sigma", data=sample_params["sigma"])
        f.create_dataset("proj_y", data=sample_params["proj_y"])
        f.create_dataset(
            "dropped_individuals", data=sample_params["dropped_individuals"]
        )
        f.attrs["lhs"] = "phenotype"
        f.attrs["rhs"] = "1"

    params = load_params(test_file)

    np.testing.assert_array_almost_equal(params.beta, sample_params["beta"])
    np.testing.assert_array_almost_equal(params.sigma, sample_params["sigma"])
    np.testing.assert_array_almost_equal(params.proj_y, sample_params["proj_y"])
    assert params.dropped_individuals == sample_params["dropped_individuals"]
