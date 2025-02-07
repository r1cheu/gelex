import numpy as np
import pytest
from phenx._chenx import LinearMixedModel


@pytest.fixture
def test_data():
    """Create test data for LMM"""
    rng = np.random.default_rng(42)
    n_samples = 100
    n_random = 2

    # Generate test data
    y = rng.normal(0, 1, (n_samples, 1))
    X = np.stack(
        [np.ones(n_samples), rng.normal(0, 1, n_samples), rng.normal(0, 1, n_samples)],
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


def test_lmm_initialization(test_data):
    """Test LinearMixedModel initialization"""
    y, X, covar_cube, names = test_data
    model = LinearMixedModel(y, X, covar_cube, names)

    assert model is not None
    assert model.n_samples == y.shape[0]
    assert model.n_fixed_effect == X.shape[1]
    assert model.n_random_effect == len(names) + 1


def test_lmm_properties(test_data):
    """Test LinearMixedModel properties"""
    y, X, covar_cube, names = test_data
    model = LinearMixedModel(y, X, covar_cube, names)

    # Test beta property
    beta = model.beta
    assert isinstance(beta, np.ndarray)
    assert beta.shape == (X.shape[1],)

    # Test sigma property
    sigma = model.sigma
    assert isinstance(sigma, np.ndarray)
    assert sigma.shape == (len(names) + 1,)


def test_lmm_reset(test_data):
    """Test reset functionality"""
    y, X, covar_cube, names = test_data
    model = LinearMixedModel(y, X, covar_cube, names)

    # Store initial values
    initial_beta = model.beta.copy()
    initial_sigma = model.sigma.copy()

    # Modify model state
    model.reset()

    # Verify reset
    assert np.array_equal(model.beta, initial_beta)
    assert np.array_equal(model.sigma, initial_sigma)


def test_lmm_repr(test_data):
    """Test string representation"""
    y, X, covar_cube, names = test_data
    model = LinearMixedModel(y, X, covar_cube, names)

    rep = repr(model)
    assert "Linear Mixed Model" in rep
    assert str(y.shape[0]) in rep
    assert str(X.shape[1]) in rep
    assert all(name in rep for name in names)


def test_lmm_invalid_input(test_data):
    """Test invalid input handling"""
    y, X, covar_cube, names = test_data

    # Test mismatched dimensions
    with pytest.raises(RuntimeError):
        LinearMixedModel(y[:-1], X, covar_cube, names)


def test_lmm_noconvert(test_data):
    """Test noconvert functionality"""
    y, X, covar_cube, names = test_data

    # Test with non-contiguous arrays
    y_noncontig = y[::2]
    X_noncontig = X[::2, ::2]
    cube_noncontig = covar_cube[::2, ::2, ::2]

    with pytest.raises(TypeError):
        LinearMixedModel(y_noncontig, X_noncontig, cube_noncontig, names)
