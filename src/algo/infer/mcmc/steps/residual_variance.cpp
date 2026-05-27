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

#include "gelex/algo/infer/mcmc/steps/residual_variance.h"

namespace gelex::mcmc
{

ResidualVarianceStep::ResidualVarianceStep(
    Eigen::Index num_individuals,
    const bayes::ResidualPrior& prior,
    bayes::ResidualState& state,
    std::mt19937_64& rng)
    : num_individuals_(num_individuals),
      state_(state),
      rng_(rng),
      kernel_(prior)
{
}

auto ResidualVarianceStep::step() -> void
{
    kernel_.prepare();
    state_.variance = kernel_.sample(num_individuals_, state_.y_adj, rng_);
}

}  // namespace gelex::mcmc
