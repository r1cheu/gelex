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
#include <string_view>

#include <omp.h>
#include <Eigen/Core>

#include "gelex/algo/infer/params.h"
#include "gelex/algo/infer/posterior_calculator.h"
#include "gelex/algo/report/bayes_logger.h"
#include "gelex/infra/logging/fit_event.h"
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
        std::string_view sample_prefix = "",
        const FitObserver& observer = {});

   private:
    void run_impl(
        const BayesModel& model,
        MCMCSamples& samples,
        Eigen::Index seed,
        const FitObserver& observer);

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
    std::string_view sample_prefix,
    const FitObserver& observer)
{
    MCMCSamples samples(params_, model, sample_prefix);

    logger_.log_model_information(model);

    const detail::EigenThreadGuard guard;
    omp_set_num_threads(1);
    run_impl(model, samples, seed, observer);

    if (observer)
    {
        observer(
            FitMcmcProgressEvent{
                .current = static_cast<size_t>(params_.n_iters),
                .total = static_cast<size_t>(params_.n_iters),
                .done = true,
                .h2 = std::nullopt,
                .h2_dom = std::nullopt,
                .sigma2_e = std::nullopt,
            });
    }

    MCMCResult result(std::move(samples), model, 0.9);
    result.compute();
    logger_.log_result(result, model, params_.n_records);

    return result;
}

template <typename TraitSampler>
void MCMC<TraitSampler>::run_impl(
    const BayesModel& model,
    MCMCSamples& samples,
    Eigen::Index seed,
    const FitObserver& observer)
{
    BayesState status{model};

    std::mt19937_64 rng(seed);
    Eigen::Index record_idx = 0;

    for (Eigen::Index iter = 0; iter < params_.n_iters; ++iter)
    {
        trait_sampler_(model, status, rng);

        status.compute_heritability();

        if (observer)
        {
            observer(
                FitMcmcProgressEvent{
                    .current = static_cast<size_t>(iter + 1),
                    .total = static_cast<size_t>(params_.n_iters),
                    .done = false,
                    .h2 = status.additive()
                              ? std::optional{status.additive()->heritability}
                              : std::nullopt,
                    .h2_dom
                    = status.dominant()
                          ? std::optional{status.dominant()->heritability}
                          : std::nullopt,
                    .sigma2_e = status.residual().variance,
                });
        }

        if (iter >= params_.n_burnin
            && (iter + 1 - params_.n_burnin) % params_.n_thin == 0)
        {
            samples.store(status, record_idx++);
        }
    }
}

}  // namespace gelex

#endif  // GELEX_ESTIMATOR_BAYES_MCMC_H_
