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

#include "gelex/algo/infer/mcmc/steps/genetic_proportion.h"
#include <random>

#include "gelex/model/bayes/capabilities.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/state.h"
#include "gelex/model/bayes/state_capabilities.h"

namespace gelex::mcmc
{

SingleProportionStep::SingleProportionStep(
    const bayes::SingleGeneticPrior& prior,
    bayes::SingleGeneticBlockState& block,
    std::mt19937_64& rng)
    : proportion_(block.prior_state()
                      .get<bayes::SingleSampledMixtureProportionStateCap>()
                      .proportion()),
      dirichlet_(prior.get<bayes::SingleSampledMixtureProportionCap>()
                     .proportion()
                     .parameter()
                     .prior()
                     .concentration()),
      rng_(rng)
{
}

auto SingleProportionStep::step() -> void
{
    dirichlet_.reset();
    proportion_.value = dirichlet_(proportion_.assignment.count, rng_);
}

JointProportionStep::JointProportionStep(
    const bayes::JointGeneticPrior& prior,
    bayes::JointGeneticBlockState& block,
    std::mt19937_64& rng)
    : proportion_(block.prior_state()
                      .get<bayes::JointSampledMixtureProportionStateCap>()
                      .proportion()),
      dirichlet_(prior.get<bayes::JointSampledMixtureProportionCap>()
                     .proportion()
                     .parameter()
                     .prior()
                     .concentration()),
      rng_(rng)
{
}

auto JointProportionStep::step() -> void
{
    dirichlet_.reset();
    proportion_.value = dirichlet_(proportion_.assignment.count, rng_);
}

}  // namespace gelex::mcmc
