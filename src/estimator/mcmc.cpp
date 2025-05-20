#include <chrono>
#include <cstddef>

#include "gelex/estimator/gibbs/base.h"
#include "gelex/estimator/mcmc.h"
#include "gelex/estimator/mcmc_logger.h"
#include "gelex/estimator/mcmc_result.h"
#include "gelex/model/bayes_policy.h"
#include "gelex/utils.h"

namespace gelex
{

MCMC::MCMC(MCMCParams params) : params_(params), rng_(params.seed) {}

const MCMCResult& MCMC::run(Bayes& model, size_t log_freq)
{
    y_adj_ = model.phenotype();
    y_adj_ -= model.mu().value;

    storage_ = std::make_unique<MCMCStorage>(params_);
    storage_->initialize(model);
    size_t record_idx = 0;

    logger_.log_model_information(model, params_);

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < params_.iter; ++i)
    {
        double sigma_e = model.residual().value;
        sample_mu(model);
        if (model.fixed().exist)
        {
            sample_fixed_effect(model.fixed(), sigma_e);
        }
        for (auto& eff : model.random())
        {
            sample_random_effect(eff, sigma_e);
        }

        for (auto& eff : model.genetic())
        {
            sample_genetic_effect(eff, sigma_e);
        }

        const double ssq = arma::dot(y_adj_, y_adj_);
        auto& residual = model.residual();
        residual.value = sample_scale_inv_chi_squared(
            rng_,
            residual.prior.nu + static_cast<double>(model.n_individuals()),
            (ssq + (residual.prior.s2 * residual.prior.nu)));

        if (i == params_.n_burnin)
        {
            logger_.log_burnin_finished();
        }

        if (i >= params_.n_burnin && i % params_.n_thin == 0)
        {
            storage_->store(model, record_idx++);
        }

        if (log_freq > 0 && i % log_freq == 0)
        {
            logger_.log_iteration(
                i, model, compute_time_left(start, i, params_.iter));
        }
    }
    logger_.log_iteration(
        params_.iter,
        model,
        compute_time_left(start, params_.iter, params_.iter));
    result_ = compute_mcmc_result(*storage_, model);
    return result_;
}

void MCMC::sample_mu(Bayes& model)
{
    const auto n = static_cast<double>(model.phenotype().n_elem);
    std::normal_distribution<double> adjustment_sampler{
        arma::sum(y_adj_) / n, sqrt(model.residual().value / n)};
    double adj = -adjustment_sampler(rng_);
    model.mu().value -= adj;
    y_adj_ += adj;
}

void MCMC::sample_fixed_effect(bayes::FixedEffect& effect, double sigma_e)
{
    dvec& coeff = effect.coeff;
    const dvec& cols_norm = effect.cols_norm;
    const dmat& design_mat = effect.design_mat;
    const dvec& post_sigma = arma::sqrt(sigma_e / cols_norm);
    std::normal_distribution<double> normal{0, 1};
    double* y_adj = y_adj_.memptr();

    const int n = static_cast<int>(design_mat.n_rows);
    for (size_t i = 0; i < coeff.n_elem; ++i)
    {
        const double old_i = coeff.at(i);
        const double* col_i = design_mat.colptr(i);
        const double norm = cols_norm.at(i);

        const double rhs = ddot_ptr(n, col_i, y_adj) + (norm * old_i);
        const double new_i = (normal(rng_) * post_sigma.at(i)) + (rhs / norm);
        coeff.at(i) = new_i;
        daxpy_ptr(n, old_i - new_i, col_i, y_adj);
    }
}

void MCMC::sample_random_effect(bayes::RandomEffect& effect, double sigma_e)
{
    dvec& coeff = effect.coeff;
    const dvec& cols_norm = effect.cols_norm;
    const dmat& design_mat = effect.design_mat;
    const double sigma = arma::as_scalar(effect.sigma);
    const dvec inv_scaler = 1 / (cols_norm + sigma_e / sigma);
    const dvec post_sigma = arma::sqrt(sigma_e * inv_scaler);

    std::normal_distribution<double> normal{0, 1};

    const int n = static_cast<int>(design_mat.n_rows);
    double* y_adj = y_adj_.memptr();

    for (size_t i = 0; i < coeff.n_elem; ++i)
    {
        const double old_i = coeff.at(i);
        const double* col_i = design_mat.memptr();
        const double norm = cols_norm.at(i);

        const double rhs = ddot_ptr(n, col_i, y_adj) + (norm * old_i);
        const double new_i
            = (normal(rng_) * post_sigma.at(i)) + (rhs * inv_scaler.at(i));

        coeff.at(i) = new_i;
        daxpy_ptr(n, old_i - new_i, col_i, y_adj);
    }
    const double ssq = arma::dot(effect.coeff, effect.coeff);
    effect.sigma = sample_scale_inv_chi_squared(
        rng_,
        effect.prior.nu + static_cast<double>(coeff.n_elem),
        (ssq + (effect.prior.s2 * effect.prior.nu)));
}

void MCMC::sample_genetic_effect(bayes::GeneticEffect& effect, double sigma_e)
{
    bayes_trait_sample.at(to_index(effect.bayes.type))(
        effect, y_adj_.memptr(), snp_tracker_, sigma_e, rng_);
}

}  // namespace gelex
