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

#ifndef GELEX_ALGO_MCMC_STEPS_FIXED_H_
#define GELEX_ALGO_MCMC_STEPS_FIXED_H_

#include <random>

#include "gelex/bayes/legacy_state.h"
#include "gelex/infra/stats/normal_sampler.h"
#include "gelex/types/fixed_designs.h"

namespace gelex
{

class FixedStep final
{
   public:
    FixedStep(
        const FixedDesign& design,
        bayes::FixedState& state,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    auto step() -> void;

   private:
    const FixedDesign& design_;
    bayes::FixedState& state_;
    bayes::ResidualState& residual_;
    std::mt19937_64& rng_;
    NormalSampler<double> normal_{0.0};
};

}  // namespace gelex

#endif  // GELEX_ALGO_MCMC_STEPS_FIXED_H_
