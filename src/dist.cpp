#include "gelex/dist.h"

namespace gelex
{
dvec dirichlet(const uvec& alphas, std::mt19937_64& rng)
{
    dvec pi(alphas.n_elem, arma::fill::zeros);
    double sum = 0.0;
    for (size_t i = 0; i < alphas.size(); ++i)
    {
        std::gamma_distribution<double> gamma_dist(alphas[i], 1.0);
        pi.at(i) = gamma_dist(rng);
    }
    return pi / arma::sum(pi);
}
}  // namespace gelex
