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

#ifndef GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_PROPORTION_H_
#define GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_PROPORTION_H_

#include <random>

#include "gelex/algo/infer/mcmc/step.h"
#include "gelex/bayes/state.h"
#include "gelex/infra/stats/conjugate_prior.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class SingleProportionStep final : public Step
{
   public:
    SingleProportionStep(
        bayes::SingleGeneticBlockState& block,
        const bayes::SingleGeneticPrior& prior,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::SampledProportionState& proportion_;
    stats::DirichletSampler<double> dirichlet_;
    std::mt19937_64& rng_;
};

class JointProportionStep final : public Step
{
   public:
    JointProportionStep(
        bayes::JointGeneticBlockState& block,
        const bayes::JointGeneticPrior& prior,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::SampledProportionState& proportion_;
    stats::DirichletSampler<double> dirichlet_;
    std::mt19937_64& rng_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_PROPORTION_H_
