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

#include "gelex/model/bayes/samplers/detail/additive.h"

#include <random>

#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/samplers/detail/gibbs/a.h"
#include "gelex/model/bayes/samplers/detail/gibbs/b.h"
#include "gelex/model/bayes/samplers/detail/gibbs/c.h"
#include "gelex/model/bayes/samplers/detail/gibbs/r.h"
#include "gelex/model/bayes/samplers/detail/gibbs/rr.h"

namespace gelex::detail::AdditiveSampler
{

auto A::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* effect = model.additive();
    auto* state = states.additive();
    auto& residual = states.residual();
    Gibbs::A(*effect, *state, residual, rng);
}

auto B::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* effect = model.additive();
    auto* state = states.additive();
    auto& residual = states.residual();
    Gibbs::B(*effect, *state, residual, rng);
}

auto C::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* effect = model.additive();
    auto* state = states.additive();
    auto& residual = states.residual();
    Gibbs::C(*effect, *state, residual, rng);
}

auto R::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* effect = model.additive();
    auto* state = states.additive();
    auto& residual = states.residual();
    Gibbs::R(*effect, *state, residual, rng);
}

auto RR::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* effect = model.additive();
    auto* state = states.additive();
    auto& residual = states.residual();
    Gibbs::RR(*effect, *state, residual, rng);
}

}  // namespace gelex::detail::AdditiveSampler
