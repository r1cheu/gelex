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

#ifndef GELEX_ALGO_MCMC_STEPS_RESIDUAL_H_
#define GELEX_ALGO_MCMC_STEPS_RESIDUAL_H_

#include <random>

#include <Eigen/Core>

#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/infra/stats/scaled_inv_chi2_sampler.h"

namespace gelex
{

class ResidualStep final
{
   public:
    ResidualStep(
        Eigen::Index num_individuals,
        const bayes::ResidualPrior& prior,
        bayes::ResidualState& state,
        std::mt19937_64& rng);

    auto step() -> void;

   private:
    Eigen::Index num_individuals_{};
    bayes::ResidualState& state_;
    std::mt19937_64& rng_;
    ScaledInvChi2Sampler<double> sampler_;
};

}  // namespace gelex

#endif  // GELEX_ALGO_MCMC_STEPS_RESIDUAL_H_
