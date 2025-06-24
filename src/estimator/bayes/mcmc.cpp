#include <atomic>
#include <cstddef>
#include <memory>
#include <thread>
#include <vector>

#include <fmt/base.h>
#include <fmt/format.h>
#include <gelex/barkeep.h>

#include <omp.h>

#include "armadillo"
#include "gelex/dist.h"
#include "gelex/estimator/bayes/base.h"
#include "gelex/estimator/bayes/indicator.h"
#include "gelex/estimator/bayes/logger.h"
#include "gelex/estimator/bayes/mcmc.h"
#include "gelex/estimator/bayes/result.h"
#include "gelex/estimator/bayes/samples.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/policy.h"
#include "gelex/utils/formatter.h"

namespace gelex
{
using arma::dmat;
using arma::dvec;
using arma::uvec;

namespace bk = barkeep;
MCMC::MCMC(MCMCParams params) : params_(params) {}

const MCMCResult& MCMC::run(const BayesModel& model, size_t seed)
{
    samples_ = std::make_unique<MCMCSamples>(params_, model);
    result_ = std::make_unique<MCMCResult>(*samples_);

    const size_t n_chains = params_.n_chains;
    std::vector<std::atomic<size_t>> idxs(n_chains);
    for (auto& i : idxs)
    {
        i = 0;
    }

    auto status_names = Indicator::create_status_names(model);
    Indicator indicator(n_chains, params_.n_iters, idxs, status_names);

    logger_.log_model_information(model, params_);

    indicator.show();
    std::vector<std::thread> threads;
    threads.reserve(n_chains);

    for (size_t i = 0; i < n_chains; ++i)
    {
        threads.emplace_back(
            [this, &model, seed, i, &idxs, &indicator]
            {
                omp_set_num_threads(1);
                run_one_chain(model, i, seed + i, idxs[i], indicator);
            });
    }
    for (auto& t : threads)
    {
        t.join();
    }

    indicator.done();
    result_->compute_summary_statistics(*samples_, 0.9);
    logger_.log_result(*result_, model);

    return *result_;
}

void MCMC::run_one_chain(
    const BayesModel& model,
    size_t chain,
    size_t seed,
    std::atomic_size_t& iter,
    Indicator& indicator)
{
    std::mt19937_64 rng(seed);

    arma::dvec y_adj = model.phenotype();
    BayesStatus status(model);

    size_t record_idx = 0;
    uvec snp_tracker(model.genetic()[0].design_mat.n_cols, arma::fill::zeros);

    ScaledInvChiSq residuel_sample{model.residual().prior};

    for (; iter < params_.n_iters; ++iter)
    {
        double sigma_e = status.residual.value;

        sample_fixed_effect(
            model.fixed(), status.fixed, y_adj.memptr(), sigma_e, rng);
        for (size_t i = 0; i < model.random().size(); ++i)
        {
            sample_random_effect(
                model.random()[i],
                status.random[i],
                y_adj.memptr(),
                sigma_e,
                rng);
            indicator.update(
                chain,
                sigma_squared("_" + model.random()[i].name),
                status.random[i].sigma.at(0));
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
            indicator.update(
                chain,
                sigma_squared("_" + model.genetic()[i].name),
                status.genetic[i].genetic_var);

            if (bayes_trait_estimate_pi[to_index(model.genetic()[i].type)])
            {
                for (size_t j = 0; j < model.genetic()[i].pi.size(); ++j)
                {
                    indicator.update(
                        chain,
                        fmt::format("π{}_{}", j, model.genetic()[i].name),
                        status.genetic[i].pi.prop.at(j));
                }
            }
        }

        auto& residual = status.residual;
        residuel_sample.update(arma::dot(y_adj, y_adj), model.n_individuals());
        residual.value = residuel_sample.sample(rng);
        status.compute_heritability();

        for (size_t i = 0; i < model.genetic().size(); ++i)
        {
            indicator.update(
                chain,
                h2("_" + model.genetic()[i].name),
                status.genetic[i].heritability);
        }

        indicator.update(chain, sigma_squared("_e"), status.residual.value);

        if (iter >= params_.n_burnin
            && (iter + 1 - params_.n_burnin) % params_.n_thin == 0)
        {
            std::lock_guard<std::mutex> lock(samples_mutex_);
            samples_->store(status, record_idx++, chain);
        }
    }
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
    ScaledInvChiSq chi_squared{design.prior};
    chi_squared.update(arma::dot(state.coeff, state.coeff), coeff.n_elem);
    state.sigma.at(0) = chi_squared.sample(rng);
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
