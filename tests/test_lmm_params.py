import numpy as np
import pytest
from phenx._chenx import LinearMixedModelParams


@pytest.fixture
def sample_params():
    # Create sample data
    beta = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    sigma = np.array([0.5, 0.3, 0.2], dtype=np.float64)
    individuals = ["person1", "person2", "person3"]
    dropped_individuals = ["dropped1", "dropped2"]

    return {
        "beta": beta,
        "sigma": sigma,
        "individuals": individuals,
        "dropped_individuals": dropped_individuals,
    }


def test_linear_mixed_model_params_init(sample_params):
    # Test initialization
    params = LinearMixedModelParams(
        beta=sample_params["beta"],
        sigma=sample_params["sigma"],
        individuals=sample_params["individuals"],
        dropped_individuals=sample_params["dropped_individuals"],
    )

    assert params is not None


def test_linear_mixed_model_params_properties(sample_params):
    # Create instance
    params = LinearMixedModelParams(
        beta=sample_params["beta"],
        sigma=sample_params["sigma"],
        individuals=sample_params["individuals"],
        dropped_individuals=sample_params["dropped_individuals"],
    )

    # Test beta property
    np.testing.assert_array_almost_equal(params.beta, sample_params["beta"])

    # Test sigma property
    np.testing.assert_array_almost_equal(params.sigma, sample_params["sigma"])

    # Test individuals property
    assert params.individuals == sample_params["individuals"]

    # Test dropped_individuals property
    assert params.dropped_individuals == sample_params["dropped_individuals"]


def test_linear_mixed_model_params_invalid_inputs():
    with pytest.raises(TypeError):
        LinearMixedModelParams(
            beta=[1, 2, 3],  # should be numpy array
            sigma=np.array([0.5, 0.3, 0.2]),
            individuals=["person1"],
            dropped_individuals=["dropped1"],
        )

    with pytest.raises(TypeError):
        LinearMixedModelParams(
            beta=np.array([1.0, 2.0, 3.0]),
            sigma=[0.5, 0.3, 0.2],  # should be numpy array
            individuals=["person1"],
            dropped_individuals=["dropped1"],
        )


def test_linear_mixed_model_params_empty_arrays():
    # Test with empty arrays
    params = LinearMixedModelParams(
        beta=np.array([]), sigma=np.array([]), individuals=[], dropped_individuals=[]
    )

    assert len(params.beta) == 0
    assert len(params.sigma) == 0
    assert len(params.individuals) == 0
    assert len(params.dropped_individuals) == 0
