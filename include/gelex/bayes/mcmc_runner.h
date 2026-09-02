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

#ifndef GELEX_BAYES_MCMC_RUNNER_H_
#define GELEX_BAYES_MCMC_RUNNER_H_

#include <cstddef>
#include <cstdint>
#include <random>

#include "gelex/bayes/draws.h"
#include "gelex/bayes/kernel.h"
#include "gelex/bayes/mcmc_progress.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/infra/notify.h"

namespace gelex
{

class MCMCRunner
{
   public:
    MCMCRunner(int iterations, int burn_in, int thin);

    [[nodiscard]] auto draw_count() const noexcept -> std::uint64_t
    {
        return static_cast<std::uint64_t>((iterations_ - burn_in_) / thin_);
    }

    template <typename GeneticPrior>
    auto run(
        const BayesModel& model,
        const BayesPrior<GeneticPrior>& prior,
        BayesDraws<GeneticPrior>& draws,
        int seed = 42,
        const MCMCObserver& observer = {}) -> void
    {
        auto state = make_state(prior, model);
        BayesKernel kernel(prior);
        auto rng
            = std::mt19937_64{static_cast<std::mt19937_64::result_type>(seed)};

        for (int iteration = 0; iteration < iterations_; ++iteration)
        {
            kernel.step(model, state, rng);
            if (iteration >= burn_in_
                && (iteration + 1 - burn_in_) % thin_ == 0)
            {
                draws.append(state);
            }
            notify(
                observer,
                MCMCProgressEvent{
                    .current = static_cast<std::size_t>(iteration + 1),
                    .total = static_cast<std::size_t>(iterations_)});
        }

        notify(
            observer,
            MCMCProgressEvent{
                .current = static_cast<std::size_t>(iterations_),
                .total = static_cast<std::size_t>(iterations_),
                .done = true});
    }

   private:
    int iterations_;
    int burn_in_;
    int thin_;
};

}  // namespace gelex

#endif  // GELEX_BAYES_MCMC_RUNNER_H_
