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

#ifndef GELEX_ALGO_INFER_MCMC_STEPS_RANDOM_COEFFICIENT_H_
#define GELEX_ALGO_INFER_MCMC_STEPS_RANDOM_COEFFICIENT_H_

#include <random>
#include <span>

#include "gelex/algo/infer/mcmc/kernels/random_coefficient.h"
#include "gelex/algo/infer/mcmc/step.h"
#include "gelex/algo/infer/mcmc/sweeps/random_coefficient.h"
#include "gelex/model/bayes/designs.h"
#include "gelex/model/bayes/state.h"

namespace gelex::mcmc
{

class RandomCoefficientStep final : public Step
{
   public:
    RandomCoefficientStep(
        std::span<const bayes::RandomDesign> designs,
        std::span<bayes::RandomState> states,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    RandomCoefficientSweep sweep_;
    RandomCoefficientKernel kernel_;
    std::mt19937_64& rng_;
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_STEPS_RANDOM_COEFFICIENT_H_
