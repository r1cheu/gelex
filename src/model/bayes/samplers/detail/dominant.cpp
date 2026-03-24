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

#include "gelex/model/bayes/samplers/detail/dominant.h"

#include <random>

#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/samplers/detail/gibbs/a.h"
#include "gelex/model/bayes/samplers/detail/gibbs/at.h"
#include "gelex/model/bayes/samplers/detail/gibbs/b.h"
#include "gelex/model/bayes/samplers/detail/gibbs/c.h"
#include "gelex/model/bayes/samplers/detail/gibbs/r.h"
#include "gelex/model/bayes/samplers/detail/gibbs/rr.h"

namespace gelex::detail::DominantSampler
{

auto A::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    Gibbs::A(
        *model.genetic(GeneticKind::Dom),
        *model.priors().genetic(GeneticKind::Dom),
        *states.genetic(GeneticKind::Dom),
        states.residual(),
        rng);
}

auto B::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    Gibbs::B(
        *model.genetic(GeneticKind::Dom),
        *model.priors().genetic(GeneticKind::Dom),
        *states.genetic(GeneticKind::Dom),
        states.residual(),
        rng);
}

auto C::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    Gibbs::C(
        *model.genetic(GeneticKind::Dom),
        *model.priors().genetic(GeneticKind::Dom),
        *states.genetic(GeneticKind::Dom),
        states.residual(),
        rng);
}

auto R::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    Gibbs::R(
        *model.genetic(GeneticKind::Dom),
        *model.priors().genetic(GeneticKind::Dom),
        *states.genetic(GeneticKind::Dom),
        states.residual(),
        rng);
}

auto RR::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    Gibbs::RR(
        *model.genetic(GeneticKind::Dom),
        *model.priors().genetic(GeneticKind::Dom),
        *states.genetic(GeneticKind::Dom),
        states.residual(),
        rng);
}

auto AT::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    Gibbs::R<Gibbs::policy::AT>(
        *model.genetic(GeneticKind::Dom),
        *model.priors().genetic(GeneticKind::Dom),
        *states.genetic(GeneticKind::Dom),
        states.residual(),
        rng);
}

}  // namespace gelex::detail::DominantSampler
