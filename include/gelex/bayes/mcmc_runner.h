// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_MCMC_RUNNER_H_
#define GELEX_BAYES_MCMC_RUNNER_H_

#include <cstddef>
#include <cstdint>
#include <functional>
#include <random>

#include "gelex/bayes/draws.h"
#include "gelex/bayes/kernel.h"
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
        const std::function<void(std::size_t)>& observer = {}) -> void
    {
        auto state = make_state(prior, model);
        auto kernel = make_kernel(prior);
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
            notify(observer, static_cast<std::size_t>(iteration + 1));
        }
    }

   private:
    int iterations_;
    int burn_in_;
    int thin_;
};

}  // namespace gelex

#endif  // GELEX_BAYES_MCMC_RUNNER_H_
