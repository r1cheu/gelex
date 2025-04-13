#include "gelex/optim/optimizers.h"

#include <armadillo>
#include <cstdint>

namespace gelex
{

void ExpectationMaximizationOptimizer::step_inner(GBLUP& model)
{
    prepare_proj(model);
    dvec sigma = model.sigma();
    dvec sigma_2 = {sigma % sigma};
    uint64_t idx{};
    auto n = static_cast<double>(n_individuals());
    for (const sp_dmat& mat : model.env_cov_mats())
    {
        sigma.at(idx) = as_scalar(
                            sigma_2.at(idx) * proj_y().t() * mat * proj_y()
                            + ((arma::trace(-sigma_2.at(idx) * proj() * mat))
                               + (sigma.at(idx) * n)))
                        / n;
        ++idx;
    }
    for (const dmat& mat : model.genetic_cov_mats())
    {
        sigma.at(idx) = as_scalar(
                            sigma_2.at(idx) * proj_y().t() * mat * proj_y()
                            + ((arma::trace(-sigma_2.at(idx) * proj() * mat))
                               + (sigma.at(idx) * n)))
                        / n;
        ++idx;
    }
    sigma.back()
        = as_scalar(
              sigma_2.back() * proj_y().t() * proj_y()
              + ((arma::trace(-sigma_2.back() * proj())) + (sigma.back() * n)))
          / n;
    sigma = constrain(sigma, phenotype_var());
    model.set_sigma(sigma);
};

void AverageInformationOptimizer::step_inner(GBLUP& model)
{
    prepare_proj(model);
    compute_first_grad(model);
    dvec sigma = model.sigma();
    dmat hess = compute_hess();
    dmat hess_inv;
    if (!pinv(hess_inv, hess))
    {
        throw std::runtime_error("Hessian matrix is not invertible!");
    }
    dvec delta = -hess_inv * first_grad();
    if (loglike_diff() > 1)
    {
        sigma += 0.316 * delta;
    }
    else
    {
        sigma += delta;
    }
    sigma = constrain(sigma, phenotype_var());

    model.set_sigma(OptimizerBase::constrain(sigma, phenotype_var()));
};

dmat AverageInformationOptimizer::compute_hess()
{
    uint64_t n{n_sigma()};
    dmat hess(n, n, arma::fill::zeros);
    for (uint64_t i{0}; i < n; ++i)
    {
        for (uint64_t j{i}; j < n; ++j)
        {
            hess.at(i, j)
                = -0.5
                  * arma::as_scalar(
                      dvpy().unsafe_col(i).t() * proj() * dvpy().unsafe_col(j));
            if (i != j)
            {
                hess.at(j, i) = hess.at(i, j);
            }
        }
    }
    return hess;
}
};  // namespace gelex
