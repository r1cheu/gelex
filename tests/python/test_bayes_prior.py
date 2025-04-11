import numpy as np
from gelexy._core import Priors, sigma_prior


def test_sigma_prior():
    """Test the sigma_prior class bindings."""
    sp = sigma_prior(nu=5.0, s2=2.0)

    assert sp.nu == 5.0
    assert sp.s2 == 2.0

    sp.nu = 10.0
    sp.s2 = 3.5
    assert sp.nu == 10.0
    assert sp.s2 == 3.5


def test_priors_default_constructor():
    """Test the Priors class default constructor."""
    priors = Priors()

    sigma_a = priors.sigma_a()
    assert isinstance(sigma_a, sigma_prior)

    sigma_r = priors.sigma_r()
    assert isinstance(sigma_r, sigma_prior)

    sigma_e = priors.sigma_e()
    assert isinstance(sigma_e, sigma_prior)


def test_priors_pi_constructor():
    """Test the Priors class constructor with pi parameter."""
    probabilities = np.array([0.2, 0.3, 0.5], dtype=np.float64)

    priors = Priors(probabilities)

    sigma_a = priors.sigma_a()
    assert isinstance(sigma_a, sigma_prior)

    sigma_r = priors.sigma_r()
    assert isinstance(sigma_r, sigma_prior)

    sigma_e = priors.sigma_e()
    assert isinstance(sigma_e, sigma_prior)

    np.testing.assert_allclose(priors.pi(), probabilities)


def test_priors_pi_modification():
    probabilities = np.array([0.2, 0.3, 0.5], dtype=np.float64)
    priors = Priors(probabilities)

    priors.pi()[0] = 0.1
    np.testing.assert_allclose(priors.pi(), [0.1, 0.3, 0.5])


def test_sigma_prior_modification():
    """Test modifying sigma_prior objects returned by Priors methods."""
    priors = Priors()

    sigma_a = priors.sigma_a()

    initial_nu = sigma_a.nu
    sigma_a.nu = 20.0

    updated_sigma_a = priors.sigma_a()
    assert updated_sigma_a.nu == 20.0
    assert updated_sigma_a.nu != initial_nu
