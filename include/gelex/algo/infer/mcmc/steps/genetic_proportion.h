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
#include "gelex/exception.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/capabilities.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/state.h"
#include "gelex/model/bayes/state_capabilities.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class SingleGeneticProportionStep final : public Step
{
   public:
    SingleGeneticProportionStep(
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        std::mt19937_64& rng)
        : proportion_(block.prior_state()
                          .require<bayes::SingleProportionStateCap>()
                          .proportion()),
          dirichlet_(
              prior.query<bayes::SingleMixtureProportionCap>() != nullptr
                  ? prior.query<bayes::SingleMixtureProportionCap>()
                        ->proportion()
                        .parameter()
                        .prior()
                        .concentration()
                  : throw GelexException(
                        "SingleGeneticProportionStep: prior lacks mixture "
                        "proportion capability")),
          rng_(rng)
    {
    }

    auto step() -> void override
    {
        dirichlet_.reset();
        if (proportion_.update == bayes::UpdatePolicy::sampled)
        {
            proportion_.value = dirichlet_(proportion_.count, rng_);
        }
    }

   private:
    bayes::ProportionState& proportion_;
    stats::DirichletSampler<double> dirichlet_;
    std::mt19937_64& rng_;
};

class JointGeneticProportionStep final : public Step
{
   public:
    JointGeneticProportionStep(
        const bayes::JointGeneticPrior& prior,
        bayes::JointGeneticBlockState& block,
        std::mt19937_64& rng)
        : proportion_(block.prior_state()
                          .require<bayes::JointProportionStateCap>()
                          .proportion()),
          dirichlet_(
              prior.query<bayes::JointMixtureProportionCap>() != nullptr
                  ? prior.query<bayes::JointMixtureProportionCap>()
                        ->proportion()
                        .parameter()
                        .prior()
                        .concentration()
                  : throw GelexException(
                        "JointGeneticProportionStep: prior lacks mixture "
                        "proportion capability")),
          rng_(rng)
    {
    }

    auto step() -> void override
    {
        dirichlet_.reset();
        if (proportion_.update == bayes::UpdatePolicy::sampled)
        {
            proportion_.value = dirichlet_(proportion_.count, rng_);
        }
    }

   private:
    bayes::ProportionState& proportion_;
    stats::DirichletSampler<double> dirichlet_;
    std::mt19937_64& rng_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_PROPORTION_H_
