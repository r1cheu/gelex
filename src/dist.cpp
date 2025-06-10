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
}  // namespace gelex
