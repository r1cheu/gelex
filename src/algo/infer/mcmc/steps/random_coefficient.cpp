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

#include "gelex/algo/infer/mcmc/steps/random_coefficient.h"

#include <fmt/format.h>

#include "gelex/exception.h"

namespace gelex::mcmc
{

RandomCoefficientStep::RandomCoefficientStep(
    std::span<const bayes::RandomDesign> designs,
    std::span<bayes::RandomState> states,
    bayes::ResidualState& residual,
    std::mt19937_64& rng)
    : sweep_(designs, states, residual), rng_(rng)
{
    if (designs.size() != states.size())
    {
        throw GelexException(
            fmt::format(
                "RandomCoefficientStep: design/state size mismatch: {} != {}",
                designs.size(),
                states.size()));
    }
}

auto RandomCoefficientStep::step() -> void
{
    sweep_.run(kernel_, rng_);
}

}  // namespace gelex::mcmc
