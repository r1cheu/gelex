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

#ifndef GELEX_ALGO_INFER_MCMC_STEPS_RANDOM_VARIANCE_H_
#define GELEX_ALGO_INFER_MCMC_STEPS_RANDOM_VARIANCE_H_

#include <random>
#include <span>

#include "gelex/algo/infer/mcmc/step.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/state.h"

namespace gelex::mcmc
{

class RandomVarianceStep final : public Step
{
   public:
    RandomVarianceStep(
        const bayes::RandomPrior& prior,
        std::span<bayes::RandomState> states,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    std::span<bayes::RandomState> states_;
    std::mt19937_64& rng_;
    stats::ScaledInvChi2Sampler<double> sampler_;
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_STEPS_RANDOM_VARIANCE_H_
