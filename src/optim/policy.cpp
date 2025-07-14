#include "gelex/optim/policy.h"

namespace gelex
{
arma::dvec ExpectationMaximizationPolicy::apply(Optimizer& optim, GBLUP& model)
{
    optim.compute_proj(model);
    arma::dvec sigma = model.var_comp();
    arma::dvec sigma_2 = sigma % sigma;
    auto n = static_cast<double>(model.n_individuals());

    sigma.at(0)
        = arma::as_scalar(
              (sigma_2.at(0) * optim.proj_y_.t() * optim.proj_y_)
              + (arma::trace(-sigma_2.at(0) * optim.proj_)) + (sigma.at(0) * n))
          / n;

    auto compute_sigma = [&](const auto& effect, size_t idx)
    {
        sigma.at(idx++) = arma::as_scalar(
                              (sigma_2.at(idx) * optim.proj_y_.t()
                               * effect.covariance_matrix * optim.proj_y_)
                              + ((arma::trace(
                                     -sigma_2.at(idx) * optim.proj_
                                     * effect.covariance_matrix))
                                 + (sigma.at(idx) * n)))
                          / n;
    };
    optim.visitor_effects(model, compute_sigma, 1);
    return optim.constrain(sigma, optim.phenotype_var_);
}

arma::dvec AverageInformationPolicy::apply(Optimizer& optim, GBLUP& model)
{
    optim.compute_proj(model);
    optim.compute_first_grad(model);

    arma::dvec sigma = model.var_comp();
    size_t n = sigma.n_elem;
    arma::dmat hess(n, n, arma::fill::zeros);
    for (size_t i{0}; i < n; ++i)
    {
        for (size_t j{i}; j < n; ++j)
        {
            hess.at(i, j) = -0.5
                            * arma::as_scalar(
                                optim.dvpy_.unsafe_col(i).t() * optim.proj_
                                * optim.dvpy_.unsafe_col(j));
            if (i != j)
            {
                hess.at(j, i) = hess.at(i, j);
            }
        }
    }

    arma::dmat hess_inv;
    if (!arma::pinv(hess_inv, hess))
    {
        throw std::runtime_error("Hessian matrix is not invertible!");
    }
    optim.hess_inv_ = hess_inv;
    arma::dvec delta = -hess_inv * optim.first_grad_;
    sigma += delta;
    return optim.constrain(sigma, optim.phenotype_var_);
}
}  // namespace gelex
