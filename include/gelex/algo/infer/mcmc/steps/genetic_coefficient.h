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

#ifndef GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_COEFFICIENT_H_
#define GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_COEFFICIENT_H_

#include <random>

#include "gelex/algo/infer/mcmc/step.h"
#include "gelex/algo/infer/mcmc/sweeps/genetic_coefficient.h"
#include "gelex/model/bayes/designs.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/state.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
template <typename Transition>
class SingleGeneticCoefficientStep final : public Step
{
   public:
    SingleGeneticCoefficientStep(
        const bayes::GeneticDesign& design,
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng)
        : sweep_(design), transition_(prior, block, residual), rng_(rng)
    {
    }

    auto step() -> void override { sweep_.run(transition_, rng_); }

   private:
    SingleGeneticMarkerSweep sweep_;
    Transition transition_;
    std::mt19937_64& rng_;
};

template <typename Transition>
class JointGeneticCoefficientStep final : public Step
{
   public:
    JointGeneticCoefficientStep(
        const bayes::GeneticDesign& additive,
        const bayes::GeneticDesign& dominance,
        const bayes::JointGeneticPrior& prior,
        bayes::JointGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng)
        : sweep_(additive, dominance),
          transition_(prior, block, residual),
          rng_(rng)
    {
    }

    auto step() -> void override { sweep_.run(transition_, rng_); }

   private:
    JointGeneticMarkerSweep sweep_;
    Transition transition_;
    std::mt19937_64& rng_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_COEFFICIENT_H_
