#pragma once

#include "gelex/dist.h"
#include "gelex/estimator/gibbs/base.h"

namespace gelex
{
using arma::dmat;
using arma::dvec;

template <typename ModelType, bool EstimatePi>
void BayesCKernel(
    ModelType& model,
    dvec& y_adj,
    const dmat& genotype_mat,
    const dvec& cols_norm,
    const dvec& cols_var,
    Normal& normal,
    ScaleInvChiSq& chisq,
    Uniform& uniform,
    uvec& snp_tracker,
    uvec& fold,
    std::mt19937_64& gen)
{
    dvec& a = model.a();
    double sigma_a = model.sigma_a();
    const double sigma_e = model.sigma_e();
    const double log_pi0 = std::log(model.pi().at(0));
    const double log_pi1 = std::log(model.pi().at(1));

    double var_a{};
    const dvec inv_scaler = 1.0 / (cols_norm + sigma_e / sigma_a);

    for (size_t i = 0; i < a.n_elem; ++i)
    {
        if (cols_var.at(i) == 0.0)
        {
            continue;
        }
        const double old_i = a.at(i);
        const dvec& col_i = genotype_mat.unsafe_col(i);
        const double col_norm = cols_norm.at(i);
        const double inv_scaler_i = inv_scaler.at(i);

        double rhs
            = arma::dot(col_i, y_adj) + (old_i != 0 ? col_norm * old_i : 0);
        double logdetV = log((sigma_a * col_norm / sigma_e) + 1);
        double uhat = rhs * inv_scaler_i;

        const double s1_minus_s0
            = (-0.5 * (logdetV - uhat * rhs / sigma_e)) + log_pi1 - log_pi0;
        const double accept_prob = 1 / (1 + exp(s1_minus_s0));

        const bool has_effect = (uniform() >= accept_prob);
        snp_tracker.at(i) = has_effect ? 1 : 0;

        double new_i{};
        if (has_effect)
        {
            new_i = normal(uhat, sqrt(sigma_e * inv_scaler_i));
            daxpy_auto(y_adj, col_i, old_i - new_i);
            var_a += new_i * new_i;
        }
        else if (old_i != 0.0)
        {
            daxpy_auto(y_adj, col_i, old_i);
        }
        a.at(i) = new_i;
    }
    fold.at(1) = arma::sum(snp_tracker);
    fold.at(0) = a.n_elem - model.n_var_0() - fold.at(1);
    model.set_sigma_a(chisq.update(static_cast<double>(fold.at(1)), var_a));
    if constexpr (EstimatePi)
    {
        model.set_pi(dirichlet(fold, gen));
    }
}
}  // namespace gelex
