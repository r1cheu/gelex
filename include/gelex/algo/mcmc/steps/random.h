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

#ifndef GELEX_ALGO_MCMC_STEPS_RANDOM_H_
#define GELEX_ALGO_MCMC_STEPS_RANDOM_H_

#include <random>
#include <span>

#include "gelex/bayes/design.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/infra/stats/normal_sampler.h"
#include "gelex/infra/stats/scaled_inv_chi2_sampler.h"

namespace gelex::mcmc
{

class RandomStep final
{
   public:
    RandomStep(
        const bayes::RandomPrior& prior,
        std::span<const bayes::RandomDesign> designs,
        std::span<bayes::RandomState> states,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    auto step() -> void;

   private:
    std::span<const bayes::RandomDesign> designs_;
    std::span<bayes::RandomState> states_;
    bayes::ResidualState& residual_;
    std::mt19937_64& rng_;
    stats::NormalSampler<double> normal_{0.0};
    stats::ScaledInvChi2Sampler<double> variance_sampler_;
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_MCMC_STEPS_RANDOM_H_
