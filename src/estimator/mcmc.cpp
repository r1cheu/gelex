#include <chrono>
#include <cstddef>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#include <barkeep.h>
#include <fmt/format.h>

#include <omp.h>

#include "gelex/estimator/gibbs/base.h"
#include "gelex/estimator/mcmc.h"
#include "gelex/estimator/mcmc_logger.h"
#include "gelex/estimator/mcmc_result.h"
#include "gelex/estimator/mcmc_samples.h"
#include "gelex/model/bayes.h"
#include "gelex/model/bayes_effects.h"
#include "gelex/model/bayes_policy.h"
#include "gelex/utils.h"

namespace gelex
{
using arma::dmat;
using arma::dvec;
using arma::uvec;

namespace bk = barkeep;
MCMC::MCMC(MCMCParams params) : params_(params) {}

void MCMC::run(const BayesModel& model, size_t seed)
{
    omp_set_num_threads(1);
    samples_ = std::make_unique<MCMCSamples>(params_, model);

    std::vector<std::thread> threads;
    const size_t n_chains = params_.n_chains;
    std::vector<size_t> idxs(n_chains);
    auto bars = MCMCLogger::progress_bar(idxs, params_.iter);
    logger_.log_model_information(model, params_);

    threads.reserve(n_chains);

    bars->show();
    for (size_t i = 0; i < n_chains; ++i)
    {
        threads.emplace_back([this, &model, seed, i, &idxs]
                             { run_one_chain(model, i, seed + i, idxs[i]); });
    }
    for (auto& t : threads)
    {
        t.join();
    }
    bars->done();
}

void MCMC::run_one_chain(
    const BayesModel& model,
    size_t chain,
    size_t seed,
    size_t& iter)
{
    std::mt19937_64 rng(seed);

    arma::dvec y_adj = model.phenotype();
    BayesStatus status(model);
    y_adj -= model.mu().value;

    size_t record_idx = 0;
    uvec snp_tracker(model.genetic()[0].design_mat.n_cols, arma::fill::zeros);

    for (; iter < params_.iter; ++iter)
    {
        double sigma_e = status.residual.value;
        sample_mu(status.mu, y_adj, sigma_e, rng);
        if (model.fixed())
        {
            sample_fixed_effect(
                *model.fixed(), status.fixed, y_adj.memptr(), sigma_e, rng);
        }
        for (size_t i = 0; i < model.random().size(); ++i)
        {
            sample_random_effect(
                model.random()[i],
                status.random[i],
                y_adj.memptr(),
                sigma_e,
                rng);
        }
        for (size_t i = 0; i < model.genetic().size(); ++i)
        {
            sample_genetic_effect(
                model.genetic()[i],
                status.genetic[i],
                y_adj.memptr(),
                snp_tracker,
                sigma_e,
                rng);
        }

        const double ssq = arma::dot(y_adj, y_adj);
        auto& residual = status.residual;
        residual.value = sample_scale_inv_chi_squared(
            rng,
            residual.prior.nu + static_cast<double>(model.n_individuals()),
            (ssq + (residual.prior.s2 * residual.prior.nu)));

        if (iter >= params_.n_burnin
            && (iter + 1 - params_.n_burnin) % params_.n_thin == 0)
        {
            samples_->store(status, record_idx++, chain);
        }
    }
}

void MCMC::sample_mu(Mu& mu, dvec& y_adj, double sigma_e, std::mt19937_64& rng)
{
    const auto n = static_cast<double>(y_adj.n_elem);
    std::normal_distribution<double> adjustment_sampler{
        arma::mean(y_adj), sqrt(sigma_e / n)};
    double adj = adjustment_sampler(rng);
    mu.value += adj;
    y_adj -= adj;
}

void MCMC::sample_fixed_effect(
    const FixedEffectDesign& design,
    FixedEffectState& state,
    double* y_adj,

    double sigma_e,
    std::mt19937_64& rng)
{
    dvec& coeff = state.coeff;
    const dvec& cols_norm = design.cols_norm;
    const dmat& design_mat = design.design_mat;
    const dvec& post_sigma = arma::sqrt(sigma_e / cols_norm);
    std::normal_distribution<double> normal{0, 1};

    const int n = static_cast<int>(design_mat.n_rows);
    for (size_t i = 0; i < coeff.n_elem; ++i)
    {
        const double old_i = coeff.at(i);
        const double* col_i = design_mat.colptr(i);
        const double norm = cols_norm.at(i);

        const double rhs = ddot_ptr(n, col_i, y_adj) + (norm * old_i);
        const double new_i = (normal(rng) * post_sigma.at(i)) + (rhs / norm);
        coeff.at(i) = new_i;
        daxpy_ptr(n, old_i - new_i, col_i, y_adj);
    }
}

void MCMC::sample_random_effect(
    const RandomEffectDesign& design,
    RandomEffectState& state,
    double* y_adj,
    double sigma_e,
    std::mt19937_64& rng)
{
    dvec& coeff = state.coeff;
    const dvec& cols_norm = design.cols_norm;
    const dmat& design_mat = design.design_mat;
    const double sigma = arma::as_scalar(state.sigma);
    const dvec inv_scaler = 1 / (cols_norm + sigma_e / sigma);
    const dvec post_sigma = arma::sqrt(sigma_e * inv_scaler);

    std::normal_distribution<double> normal{0, 1};

    const int n = static_cast<int>(design_mat.n_rows);
    for (size_t i = 0; i < coeff.n_elem; ++i)
    {
        const double old_i = coeff.at(i);
        const double* col_i = design_mat.colptr(i);
        const double norm = cols_norm.at(i);

        const double rhs = ddot_ptr(n, col_i, y_adj) + (norm * old_i);
        const double new_i
            = (normal(rng) * post_sigma.at(i)) + (rhs * inv_scaler.at(i));

        coeff.at(i) = new_i;
        daxpy_ptr(n, old_i - new_i, col_i, y_adj);
    }
    const double ssq = arma::dot(state.coeff, state.coeff);
    state.sigma = sample_scale_inv_chi_squared(
        rng,
        design.prior.nu + static_cast<double>(coeff.n_elem),
        (ssq + (design.prior.s2 * design.prior.nu)));
}

void MCMC::sample_genetic_effect(
    const GeneticEffectDesign& design,
    GeneticEffectState& state,
    double* y_adj,
    uvec& snp_tracker,
    double sigma_e,
    std::mt19937_64& rng)
{
    bayes_trait_sample.at(to_index(design.type))(
        design, state, y_adj, snp_tracker, sigma_e, rng);
}

}  // namespace gelex
