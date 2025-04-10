import h5py
import numpy as np
import pytest
from gelexy.model import GBLUP
from gelexy.model.io import load_params


@pytest.fixture
def test_data():
    """Create test data for LMM"""
    rng = np.random.default_rng(42)
    n_samples = 100
    n_random = 2

    # Generate test data
    y = rng.normal(0, 1, (n_samples, 1))
    X = np.stack(
        [
            np.ones(n_samples),
            rng.normal(0, 1, n_samples),
            rng.normal(0, 1, n_samples),
        ],
        axis=-1,
    )
    covar_cube = np.stack([np.eye(n_samples) for _ in range(n_random)], axis=-1)
    names = ["random1", "random2"]
    return (
        np.asfortranarray(y),
        np.asfortranarray(X),
        np.asfortranarray(covar_cube),
        names,
    )


@pytest.fixture
def model():
    response = np.asfortranarray(np.array([1.0, 2.0, 3.0]).reshape(3, 1))
    design_matrix = np.asfortranarray(np.ones((3, 1)))
    grm_cube = np.asfortranarray(np.stack([np.eye(3)], axis=-1))

    random_effect_names = ["effect1"]
    model = GBLUP(
        response, design_matrix, grm_cube, random_effect_names
    )
    model._dropped_individuals = ["ind1"]
    model._lhs = "phenotype"
    model._rhs = "1"
    return model


def test_lmm_initialization(model):
    """Test GBLUP initialization"""

    assert model is not None
    assert model.num_individuals == 3
    assert model.num_fixed_effects == 1
    assert model.num_random_effects == 2


def test_lmm_properties(model):
    """Test GBLUP properties"""

    beta = model.beta
    assert isinstance(beta, np.ndarray)
    assert beta.shape == (1,)

    sigma = model.sigma
    assert isinstance(sigma, np.ndarray)
    assert sigma.shape == (2,)


def test_lmm_reset(test_data):
    """Test reset functionality"""
    y, X, covar_cube, names = test_data
    model = GBLUP(y, X, covar_cube, names)

    # Store initial values
    initial_beta = model.beta.copy()
    initial_sigma = model.sigma.copy()

    # Modify model state
    model.reset()

    # Verify reset
    np.testing.assert_array_equal(model.beta, initial_beta)
    np.testing.assert_array_equal(model.sigma, initial_sigma)


def test_lmm_invalid_input(test_data):
    """Test invalid input handling"""
    y, X, covar_cube, names = test_data

    # Test mismatched dimensions
    with pytest.raises(RuntimeError):
        GBLUP(y[:-1], X, covar_cube, names)


def test_lmm_noconvert(test_data):
    """Test noconvert functionality"""
    y, X, covar_cube, names = test_data

    y_noncontig = y[::2]
    X_noncontig = X[::2, ::2]
    cube_noncontig = covar_cube[::2, ::2, ::2]

    with pytest.raises(TypeError):
        GBLUP(y_noncontig, X_noncontig, cube_noncontig, names)


def test_lmm_save(tmp_path, model):
    save_path = tmp_path / "test.h5"
    model.save(save_path)

    with h5py.File(save_path, "r") as f:
        assert "beta" in f
        assert "sigma" in f
        assert "proj_y" in f
        assert "dropped_individuals" in f


def test_lmm_params_from_model(model):
    model_params = load_params(model)

    np.testing.assert_array_almost_equal(model_params.beta, model.beta)
    np.testing.assert_array_almost_equal(model_params.sigma, model.sigma)
    np.testing.assert_array_equal(
        model_params.dropped_individuals, model._dropped_individuals
    )
