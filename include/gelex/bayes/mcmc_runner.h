// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_MCMC_RUNNER_H_
#define GELEX_BAYES_MCMC_RUNNER_H_

#include <cstddef>
#include <functional>
#include <random>
#include <string_view>

#include "gelex/bayes/draws.h"
#include "gelex/bayes/kernel.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/sampling_plan.h"
#include "gelex/bayes/state.h"
#include "gelex/infra/notify.h"

namespace gelex
{

class MCMCRunner
{
   public:
    explicit MCMCRunner(SamplingPlan plan) : plan_{plan} {}

    [[nodiscard]] auto plan() const noexcept -> const SamplingPlan&
    {
        return plan_;
    }

    template <typename GeneticPrior>
    auto run(
        const BayesModel& model,
        const BayesPrior<GeneticPrior>& prior,
        std::string_view output_path,
        const std::function<void(std::size_t)>& observer = {}) -> void
    {
        auto state = make_state(prior, model);
        auto draws = BayesDraws{state, model, output_path, plan_.draw_count()};
        auto kernel = make_kernel(prior);
        auto rng = std::mt19937_64{
            static_cast<std::mt19937_64::result_type>(plan_.seed())};

        for (int iteration = 0; iteration < plan_.iterations(); ++iteration)
        {
            kernel.step(model, state, rng);
            if (plan_.retains(iteration))
            {
                draws << state;
            }
            notify(observer, static_cast<std::size_t>(iteration + 1));
        }
        draws.close();
    }

   private:
    SamplingPlan plan_;
};

}  // namespace gelex

#endif  // GELEX_BAYES_MCMC_RUNNER_H_
