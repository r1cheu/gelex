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
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/notify.h"
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
    MCMCSamples samples(model, sample_prefix, params_.n_records);

    notify(observer, FitModelReadyEvent{&model});

    const detail::EigenThreadGuard guard;
    omp_set_num_threads(1);
    run_impl(model, samples, seed, observer);
    samples.finalize();

    notify(
        observer,
        FitMcmcProgressEvent{
            .current = static_cast<size_t>(params_.n_iters),
            .total = static_cast<size_t>(params_.n_iters),
            .done = true,
            .genetic_heritabilities = {},
            .dom_positive_prob = std::nullopt,
            .sigma2_e = std::nullopt,
        });

    MCMCResult result(std::move(samples), model);
    result.compute();

    notify(observer, FitMcmcCompleteEvent{&result, &model, params_.n_records});

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

    std::vector<std::pair<GeneticEffectType, double>> gen_h2;
    gen_h2.reserve(status.genetics().size());

    for (Eigen::Index iter = 0; iter < params_.n_iters; ++iter)
    {
        trait_sampler_(model, status, rng);

        status.compute_heritability();

        gen_h2.clear();
        std::optional<std::vector<double>> dom_pi;
        std::optional<double> dom_positive_prob;
        for (const auto& gen : status.genetics())
        {
            gen_h2.emplace_back(gen.type, gen.heritability);
            if (gen.type == GeneticEffectType::Dom && gen.mixture)
            {
                const auto& prop = gen.mixture->pi.proportion;
                dom_pi.emplace(std::vector<double>(prop.begin(), prop.end()));
            }
            if (gen.type == GeneticEffectType::Dom && gen.sign)
            {
                dom_positive_prob = gen.sign->positive_prob;
            }
        }
        notify(
            observer,
            FitMcmcProgressEvent{
                .current = static_cast<size_t>(iter + 1),
                .total = static_cast<size_t>(params_.n_iters),
                .done = false,
                .genetic_heritabilities = gen_h2,
                .dom_positive_prob = dom_positive_prob,
                .sigma2_e = status.residual().variance,
            });

        if (iter >= params_.n_burnin
            && (iter + 1 - params_.n_burnin) % params_.n_thin == 0)
        {
            samples.store(status);
        }
    }
}

}  // namespace gelex

#endif  // GELEX_ESTIMATOR_BAYES_MCMC_H_
