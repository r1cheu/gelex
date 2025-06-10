#pragma once
#include <armadillo>
#include <random>

namespace gelex
{

arma::dvec dirichlet(const arma::uvec& alphas, std::mt19937_64& rng);

inline double
sample_scale_inv_chi_squared(std::mt19937_64& rng, double nu, double s2 = 1.0)
{
    std::chi_squared_distribution<double> chisq{nu};
    return s2 / chisq(rng);
}
}  // namespace gelex
