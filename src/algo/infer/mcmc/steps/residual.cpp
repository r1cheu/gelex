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

#include "gelex/algo/infer/mcmc/steps/residual.h"

#include <random>

#include <Eigen/Core>

#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"

namespace gelex::mcmc
{

ResidualStep::ResidualStep(
    Eigen::Index num_individuals,
    const bayes::ResidualPrior& prior,
    bayes::ResidualState& state,
    std::mt19937_64& rng)
    : num_individuals_(num_individuals),
      state_(state),
      rng_(rng),
      sampler_(prior.prior())
{
}

auto ResidualStep::step() -> void
{
    sampler_.reset();
    state_.variance = sampler_(
        {.n = num_individuals_, .sum_squares = state_.y_adj.squaredNorm()},
        rng_);
}

}  // namespace gelex::mcmc
