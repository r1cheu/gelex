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

#include "../src/model/bayes/samplers/additive.h"

#include <random>

#include "../src/model/bayes/samplers/gibbs/a.h"
#include "../src/model/bayes/samplers/gibbs/b.h"
#include "../src/model/bayes/samplers/gibbs/c.h"
#include "../src/model/bayes/samplers/gibbs/r.h"
#include "../src/model/bayes/samplers/gibbs/rr.h"
#include "../src/model/bayes/samplers/mh/rrd.h"
#include "gelex/model/bayes/model.h"

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

auto RRD::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* add_effect = model.additive();
    auto* add_state = states.additive();
    const auto* dom_effect = model.dominant();
    auto* dom_state = states.dominant();
    auto& residual = states.residual();

    MH::RRD(*add_effect, *add_state, *dom_effect, *dom_state, residual, rng);
}

}  // namespace gelex::detail::AdditiveSampler
