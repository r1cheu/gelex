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

#ifndef GELEX_ALGO_INFER_MCMC_CHAIN_H_
#define GELEX_ALGO_INFER_MCMC_CHAIN_H_

#include <random>
#include <variant>
#include <vector>

#include "gelex/algo/infer/mcmc/steps/fixed.h"
#include "gelex/algo/infer/mcmc/steps/joint_genetic_step.h"
#include "gelex/algo/infer/mcmc/steps/random.h"
#include "gelex/algo/infer/mcmc/steps/residual.h"
#include "gelex/algo/infer/mcmc/steps/single_genetic_step.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"

namespace gelex::mcmc
{

using SingleGeneticStep = std::variant<
    SingleSharedGaussianStep,
    SinglePerMarkerGaussianStep,
    SingleSharedSpikeSlabStep,
    SinglePerMarkerSpikeSlabStep,
    SingleScaledMixtureStep>;

using JointGeneticStep = std::variant<JointGaussianMixtureStep>;

class Chain
{
   public:
    Chain(
        FixedStep fixed,
        RandomStep random,
        std::vector<SingleGeneticStep> single_genetics,
        std::vector<JointGeneticStep> joint_genetics,
        ResidualStep residual,
        BayesState& state);

    auto step() -> void;

    static auto make(
        const BayesModel& model,
        const bayes::BayesPrior& prior,
        BayesState& state,
        std::mt19937_64& rng) -> Chain;

   private:
    FixedStep fixed_;
    RandomStep random_;
    std::vector<SingleGeneticStep> single_genetics_;
    std::vector<JointGeneticStep> joint_genetics_;
    ResidualStep residual_;
    BayesState& state_;
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_CHAIN_H_
