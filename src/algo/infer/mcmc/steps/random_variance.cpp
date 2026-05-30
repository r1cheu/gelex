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

#include "gelex/algo/infer/mcmc/steps/random_variance.h"

#include <random>
#include <span>

#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"

namespace gelex::mcmc
{

RandomVarianceStep::RandomVarianceStep(
    const bayes::RandomPrior& prior,
    std::span<bayes::RandomState> states,
    std::mt19937_64& rng)
    : states_(states), rng_(rng), sampler_(prior.prior())
{
}

auto RandomVarianceStep::step() -> void
{
    sampler_.reset();
    for (auto& state : states_)
    {
        state.variance = sampler_(
            {.n = state.coeffs.size(),
             .sum_squares = state.coeffs.squaredNorm()},
            rng_);
    }
}

}  // namespace gelex::mcmc
