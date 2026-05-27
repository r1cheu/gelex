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

#ifndef GELEX_ALGO_INFER_MCMC_STEPS_FIXED_COEFFICIENT_H_
#define GELEX_ALGO_INFER_MCMC_STEPS_FIXED_COEFFICIENT_H_

#include <random>

#include "gelex/algo/infer/mcmc/kernels/fixed_coefficient.h"
#include "gelex/algo/infer/mcmc/step.h"
#include "gelex/algo/infer/mcmc/sweeps/fixed_coefficient.h"
#include "gelex/model/bayes/state.h"
#include "gelex/types/fixed_designs.h"

namespace gelex::mcmc
{

class FixedCoefficientStep final : public Step
{
   public:
    FixedCoefficientStep(
        const FixedDesign& design,
        bayes::FixedState& state,
        bayes::ResidualState& residual,
        std::mt19937_64& rng)
        : sweep_(design, state, residual), rng_(rng)
    {
    }

    auto step() -> void override;

   private:
    FixedCoefficientSweep sweep_;
    FixedCoefficientKernel kernel_;
    std::mt19937_64& rng_;
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_STEPS_FIXED_COEFFICIENT_H_
