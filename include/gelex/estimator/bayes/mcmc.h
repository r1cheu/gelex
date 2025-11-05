#pragma once
#include <atomic>
#include <chrono>
#include <memory>
#include <ranges>
#include <thread>
#include <vector>

#include <fmt/format.h>
#include <gelex/barkeep.h>
#include <omp.h>
#include <Eigen/Core>

#include "../src/estimator/bayes/indicator.h"
#include "../src/logger/bayes_logger.h"
#include "../src/utils/formatter.h"
#include "estimator/bayes/posterior_calculator.h"
#include "gelex/estimator/bayes/params.h"
#include "gelex/logger.h"
#include "gelex/model/bayes/model.h"
#include "gelex/types/mcmc_results.h"
#include "gelex/types/mcmc_samples.h"

namespace gelex
{
template <typename TraitSampler>
class MCMC
{
   public:
    explicit MCMC(MCMCParams params, TraitSampler trait_sampler);
    MCMCResult run(const BayesModel& model, Eigen::Index seed = 42);

   private:
    void run_one_chain(
        const BayesModel& model,
        MCMCSamples& samples,
        Eigen::Index chain,
        Eigen::Index seed,
        std::atomic_ptrdiff_t& iter,
        detail::Indicator& indicator);

    void update_indicators(
        const BayesModel& model,
        const BayesState& states,
        Eigen::Index chain,
        detail::Indicator& indicator);

    detail::MCMCLogger logger_;

    MCMCParams params_;
    TraitSampler trait_sampler_;
};

template <typename TraitSampler>
MCMC<TraitSampler>::MCMC(MCMCParams params, TraitSampler trait_sampler)
    : params_(params), trait_sampler_(std::move(trait_sampler))
{
}

template <typename TraitSampler>
MCMCResult MCMC<TraitSampler>::run(const BayesModel& model, Eigen::Index seed)
{
    MCMCSamples samples(params_, model);

    const Eigen::Index n_chains = params_.n_chains;
    std::vector<std::atomic_ptrdiff_t> idxs(n_chains);
    detail::Indicator indicator(params_.n_iters, idxs);

    logger_.log_model_information(model);

    indicator.show();
    std::vector<std::thread> threads;
    threads.reserve(n_chains);

    const detail::EigenThreadGuard guard;

    auto start = std::chrono::high_resolution_clock::now();
    for (Eigen::Index i = 0; i < n_chains; ++i)
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
    indicator.done();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration
        = std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
              .count();

    logger_.info("");
    logger_.info("MCMC sampling completed.");
    logger_.info(
        "  - Samples collected per parameter: {}.",
        params_.n_records * params_.n_chains);
    logger_.info(
        "  - Total time elapsed: {:.2f} seconds.",
        static_cast<double>(duration) / 1000.0);

    MCMCResult result(std::move(samples), model, 0.9);
    result.compute();
    logger_.log_result(result, model);

    return result;
}

template <typename TraitSampler>
void MCMC<TraitSampler>::run_one_chain(
    const BayesModel& model,
    MCMCSamples& samples,
    Eigen::Index chain,
    Eigen::Index seed,
    std::atomic_ptrdiff_t& iter,
    detail::Indicator& indicator)
{
    auto logger = logging::get();
    std::mt19937_64 rng(seed);

    BayesState status{model};

    Eigen::Index record_idx = 0;

    for (; iter < params_.n_iters; ++iter)
    {
        // Sample all effects using trait sampler
        trait_sampler_(model, status, rng);

        // Compute heritability
        status.compute_heritability();

        // Update indicators for monitoring
        update_indicators(model, status, chain, indicator);

        // Record samples after burnin with thinning
        if (iter >= params_.n_burnin
            && (iter + 1 - params_.n_burnin) % params_.n_thin == 0)
        {
            samples.store(status, record_idx++, chain);
        }
    }
}

template <typename TraitSampler>
void MCMC<TraitSampler>::update_indicators(
    const BayesModel& model,
    const BayesState& states,
    Eigen::Index chain,
    detail::Indicator& indicator)
{
    // Update additive effect indicators
    auto update_effect_indicators
        = [&](const auto* state, const auto* effect, std::string_view prefix)
    {
        if (state != nullptr)
        {
            indicator.update(
                chain, fmt::format("{}_variance", prefix), state->variance);
            indicator.update(
                chain,
                fmt::format("{}_heritability", prefix),
                state->heritability);
            if (effect->estimate_pi)
            {
                for (Eigen::Index j = 0; j < state->pi.prop.size(); ++j)
                {
                    indicator.update(
                        chain,
                        fmt::format("mixture_proportion_{}_{}", prefix, j),
                        state->pi.prop(j));
                }
            }
        }
    };

    update_effect_indicators(states.additive(), model.additive(), "additive");
    update_effect_indicators(states.dominant(), model.dominant(), "dominant");

    // Update random effects indicators
    for (auto&& [effect, state] :
         std::views::zip(model.random(), states.random()))
    {
        // TODO(rlchen): add correct name according to effect name
        indicator.update(chain, sigma_squared("_"), state.variance);
    }

    const auto& residual = states.residual();
    indicator.update(chain, "residual_variance", residual.variance);

    indicator.flush_status(chain);
}

}  // namespace gelex
