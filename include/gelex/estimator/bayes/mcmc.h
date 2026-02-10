/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef GELEX_ESTIMATOR_BAYES_MCMC_H_
#define GELEX_ESTIMATOR_BAYES_MCMC_H_
#include <atomic>
#include <chrono>
#include <string_view>

#include <barkeep.h>
#include <omp.h>
#include <Eigen/Core>

#include "../src/estimator/bayes/posterior_calculator.h"
#include "../src/logger/bayes_logger.h"
#include "gelex/detail/indicator.h"
#include "gelex/estimator/bayes/params.h"
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
    MCMCResult run(
        const BayesModel& model,
        Eigen::Index seed = 42,
        std::string_view sample_prefix = "");

   private:
    void run_impl(
        const BayesModel& model,
        MCMCSamples& samples,
        Eigen::Index seed,
        std::atomic_ptrdiff_t& iter,
        detail::Indicator& indicator);

    void update_indicators(
        const BayesState& states,
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
MCMCResult MCMC<TraitSampler>::run(
    const BayesModel& model,
    Eigen::Index seed,
    std::string_view sample_prefix)
{
    MCMCSamples samples(params_, model, sample_prefix);

    std::atomic_ptrdiff_t iter{0};
    detail::Indicator indicator(params_.n_iters, iter);

    logger_.log_model_information(model);

    indicator.show();

    const detail::EigenThreadGuard guard;
    auto start = std::chrono::high_resolution_clock::now();
    omp_set_num_threads(1);
    run_impl(model, samples, seed, iter, indicator);

    indicator.done();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration
        = std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
              .count();

    MCMCResult result(std::move(samples), model, 0.9);
    result.compute();
    logger_.log_result(
        result,
        model,
        static_cast<double>(duration) / 1000.0,
        params_.n_records);

    return result;
}

template <typename TraitSampler>
void MCMC<TraitSampler>::run_impl(
    const BayesModel& model,
    MCMCSamples& samples,
    Eigen::Index seed,
    std::atomic_ptrdiff_t& iter,
    detail::Indicator& indicator)
{
    BayesState status{model};

    std::mt19937_64 rng(seed);
    Eigen::Index record_idx = 0;

    for (; iter < params_.n_iters; ++iter)
    {
        trait_sampler_(model, status, rng);

        status.compute_heritability();

        update_indicators(status, indicator);

        if (iter >= params_.n_burnin
            && (iter + 1 - params_.n_burnin) % params_.n_thin == 0)
        {
            samples.store(status, record_idx++);
        }
    }
}

template <typename TraitSampler>
void MCMC<TraitSampler>::update_indicators(
    const BayesState& states,
    detail::Indicator& indicator)
{
    using sm = detail::Indicator::StatusMetric;

    indicator.update(
        sm::additive_heritability, states.additive()->heritability);
    indicator.update(
        sm::dominant_heritability, states.dominant()->heritability);
    indicator.update(sm::residual_variance, states.residual().variance);

    indicator.flush_status();
}

}  // namespace gelex

#endif  // GELEX_ESTIMATOR_BAYES_MCMC_H_
