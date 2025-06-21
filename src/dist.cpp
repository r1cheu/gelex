#include "gelex/dist.h"

namespace gelex
{
arma::dvec dirichlet(const arma::uvec& alphas, std::mt19937_64& rng)
{
    arma::dvec pi(alphas.n_elem, arma::fill::zeros);
    double sum = 0.0;
    for (size_t i = 0; i < alphas.size(); ++i)
    {
        std::gamma_distribution<double> gamma_dist(alphas[i], 1.0);
        pi.at(i) = gamma_dist(rng);
    }
    return pi / arma::sum(pi);
}
ScaledInvChiSq::ScaledInvChiSq(const ScaledInvChiSqParams& prior_params)
    : params_(prior_params)
{
}
ScaledInvChiSq::ScaledInvChiSq(double initial_nu, double initial_s2)
    : params_{initial_nu, initial_s2}
{
}

void ScaledInvChiSq::update(
    double sum_of_squared_errors,
    size_t num_observations)
{
    if (num_observations <= 0)
    {
        return;
    }

    const double old_nu = params_.nu;
    const double old_s2 = params_.s2;

    const double posterior_nu = old_nu + static_cast<double>(num_observations);
    const double posterior_s2
        = ((old_nu * old_s2) + sum_of_squared_errors) / posterior_nu;
    params_ = {posterior_nu, posterior_s2};
}

void ScaledInvChiSq::update(double single_observation_squared_error)
{
    update(single_observation_squared_error, 1);
}

double ScaledInvChiSq::sample(std::mt19937_64& rng) const
{
    return sample_scale_inv_chi_squared(rng, params_.nu, params_.s2);
}
const ScaledInvChiSqParams& ScaledInvChiSq::get_params() const
{
    return params_;
}

}  // namespace gelex
