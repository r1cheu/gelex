#include <atomic>
#include <chrono>
#include <cstddef>
#include <memory>
#include <thread>
#include <vector>

#include <fmt/base.h>
#include <fmt/format.h>
#include <gelex/barkeep.h>
#include <omp.h>
#include <Eigen/Core>

#include "../src/estimator/bayes/indicator.h"
#include "../src/logger/bayes_logger.h"
#include "../src/model/bayes/distribution.h"
#include "gelex/estimator/bayes/mcmc.h"
#include "gelex/estimator/bayes/result.h"
#include "gelex/estimator/bayes/samples.h"
#include "gelex/logger.h"
#include "gelex/model/bayes/model.h"
#include "gelex/utils/formatter.h"

namespace gelex
{
using Eigen::Index;
using Eigen::Ref;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

namespace bk = barkeep;
MCMC::MCMC(MCMCParams params) : params_(params) {}

const MCMCResult& MCMC::run(const BayesModel& model, Index seed)
{
    MCMCSamples samples(params_, model);

    const size_t n_chains = params_.n_chains;
    std::vector<std::atomic<size_t>> idxs(n_chains);
    detail::Indicator indicator(model, params_.n_iters, idxs);

    logger_.log_model_information(model, params_);

    indicator.show();
    std::vector<std::thread> threads;
    threads.reserve(n_chains);
    int nthreads = Eigen::nbThreads();

    auto start = std::chrono::high_resolution_clock::now();
    for (Index i = 0; i < n_chains; ++i)
    {
        threads.emplace_back(
            [this, &model, &samples, seed, i, &idxs, &indicator]
            {
                omp_set_num_threads(1);
                Eigen::setNbThreads(1);
                run_one_chain(model, samples, i, seed + i, idxs[i], indicator);
            });
    }
    for (auto& t : threads)
    {
        t.join();
    }
    Eigen::setNbThreads(nthreads);
    indicator.done();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration
        = std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
              .count();

    logger_.info("");
    logger_.info(
        "MCMC sampling completed, Total time: {:.2f} s.",
        static_cast<double>(duration) / 1000.0);

    result_ = std::make_unique<MCMCResult>(std::move(samples), model, 0.9);
    result_->compute();
    logger_.log_result(*result_, model);

    return *result_;
}

void MCMC::run_one_chain(
    const BayesModel& model,
    MCMCSamples& samples,
    Index chain,
    Index seed,
    std::atomic_size_t& iter,
    detail::Indicator& indicator)
{
    auto logger = logging::get();
    std::mt19937_64 rng(seed);

    VectorXd y_adj = model.phenotype().array() - model.phenotype().mean();
    BayesStatus status(model);

    Index record_idx = 0;

    detail::ScaledInvChiSq residual_sample{model.residual().prior};

    for (; iter < params_.n_iters; ++iter)
    {
        double sigma_e = status.residual.value;

        if (model.fixed())
        {
            sample_fixed_effect(
                *model.fixed(), *status.fixed, y_adj, sigma_e, rng);
        }

        for (Index i = 0; i < model.random().size(); ++i)
        {
            sample_random_effect(
                model.random()[i], status.random[i], y_adj, sigma_e, rng);
            indicator.update(
                chain,
                sigma_squared("_" + model.random()[i].name),
                status.random[i].effect_variance(0));
        }

        if (model.additive())
        {
            sample_additive_effect(
                *model.additive(),
                *status.additive,
                y_adj,
                sigma_e,
                rng,
                model.trait());

            indicator.update(
                chain,
                sigma_squared("_" + model.additive()->name),
                status.additive->effect_variance);

            if (model.trait()->estimate_pi())
            {
                for (int j = 0; j < status.additive->pi.prop.size(); ++j)
                {
                    indicator.update(
                        chain,
                        fmt::format("π_{}", j),
                        status.additive->pi.prop(j));
                }
            }
        }
        auto& residual = status.residual;
        residual_sample.compute(y_adj.squaredNorm(), model.n_individuals());
        residual.value = residual_sample(rng);
        status.compute_heritability();

        if (model.additive())
        {
            indicator.update(
                chain,
                h2("_" + model.additive()->name),
                status.additive->heritability);
        }

        indicator.update(chain, sigma_squared("_e"), status.residual.value);

        if (iter >= params_.n_burnin
            && (iter + 1 - params_.n_burnin) % params_.n_thin == 0)
        {
            std::lock_guard<std::mutex> lock(samples_mutex_);
            samples.store(status, record_idx++, chain);
        }
    }
}

void MCMC::sample_fixed_effect(
    const bayes::FixedEffect& effect,
    bayes::FixedStatus& status,
    Ref<VectorXd> y_adj,
    double sigma_e,
    std::mt19937_64& rng)
{
    // for convenience
    VectorXd& coeff = status.coeff;
    const VectorXd& cols_norm = effect.cols_norm;
    const MatrixXd& design_matrix = effect.design_matrix;

    const int n = static_cast<int>(design_matrix.rows());

    std::normal_distribution<double> normal{0, 1};

    for (int i = 0; i < coeff.size(); ++i)
    {
        // for convenience
        const double old_i = coeff(i);
        const auto& col = design_matrix.col(i);
        const double norm = cols_norm(i);

        // calculate the posterior mean
        const double rhs = col.dot(y_adj) + (norm * old_i);
        const double post_mean = rhs / norm;
        const double post_stddev = std::sqrt(sigma_e / norm);

        // sample a new coefficient
        const double new_i = (normal(rng) * post_stddev) + post_mean;
        coeff(i) = new_i;

        // update the y_adj vector
        const double diff = old_i - new_i;
        y_adj.array() += diff * col.array();
    }
}

void MCMC::sample_random_effect(
    const bayes::RandomEffect& effect,
    bayes::RandomStatus& status,
    Ref<VectorXd> y_adj,
    double sigma_e,
    std::mt19937_64& rng)
{
    // for convenience
    VectorXd& coeff = status.coeff;
    const VectorXd& cols_norm = effect.cols_norm;
    const MatrixXd& design_matrix = effect.design_matrix;
    const double sigma = status.effect_variance(0);

    const int n = static_cast<int>(design_matrix.rows());

    // calculate precision kernel and posterior standard deviation
    const VectorXd inv_scaler = 1.0 / (cols_norm.array() + sigma_e / sigma);
    const VectorXd post_stddev = (sigma_e * inv_scaler.array()).sqrt();

    // Setup distributions for sampling
    std::normal_distribution<double> normal{0, 1};

    for (int i = 0; i < coeff.size(); ++i)
    {
        // for convenience
        const double old_i = coeff(i);
        const auto& col = design_matrix.col(i);
        const double norm = cols_norm(i);

        // calculate the posterior mean
        const double rhs = col.dot(y_adj) + (norm * old_i);
        const double post_mean = rhs * inv_scaler(i);

        // sample a new coefficient
        const double new_i = (normal(rng) * post_stddev(i)) + post_mean;
        coeff(i) = new_i;

        // update the y_adj vector
        const double diff = old_i - new_i;
        y_adj.array() += col.array() * diff;
    }

    // sample a new variance
    detail::ScaledInvChiSq chi_squared{effect.prior};
    chi_squared.compute(coeff.squaredNorm(), coeff.size());
    status.effect_variance(0) = chi_squared(rng);
}

void MCMC::sample_additive_effect(
    const bayes::AdditiveEffect& effect,
    bayes::AdditiveStatus& status,
    Ref<VectorXd> y_adj,
    double sigma_e,
    std::mt19937_64& rng,
    const std::unique_ptr<GeneticTrait>& trait)
{
    (*trait)(effect, status, y_adj, sigma_e, rng);
}

}  // namespace gelex
