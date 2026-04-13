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

#include "gelex/algo/infer/mcmc/samplers/dominant.h"

#include <random>

#include "gelex/algo/infer/mcmc/samplers/gibbs/a.h"
#include "gelex/algo/infer/mcmc/samplers/gibbs/at.h"
#include "gelex/algo/infer/mcmc/samplers/gibbs/b.h"
#include "gelex/algo/infer/mcmc/samplers/gibbs/c.h"
#include "gelex/algo/infer/mcmc/samplers/gibbs/r.h"
#include "gelex/algo/infer/mcmc/samplers/gibbs/rr.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"

namespace gelex::detail::DominantSampler
{

auto A::operator()(
    const BayesModel& model,
    const bayes::Priors& priors,
    mcmc::State& states,
    std::mt19937_64& rng) const -> void
{
    Gibbs::A(
        *model.genetic(GeneticMode::D),
        *priors.genetic(GeneticMode::D),
        *states.genetic(GeneticMode::D),
        states.residual(),
        rng);
}

auto B::operator()(
    const BayesModel& model,
    const bayes::Priors& priors,
    mcmc::State& states,
    std::mt19937_64& rng) const -> void
{
    Gibbs::B(
        *model.genetic(GeneticMode::D),
        *priors.genetic(GeneticMode::D),
        *states.genetic(GeneticMode::D),
        states.residual(),
        rng);
}

auto C::operator()(
    const BayesModel& model,
    const bayes::Priors& priors,
    mcmc::State& states,
    std::mt19937_64& rng) const -> void
{
    Gibbs::C(
        *model.genetic(GeneticMode::D),
        *priors.genetic(GeneticMode::D),
        *states.genetic(GeneticMode::D),
        states.residual(),
        rng);
}

auto R::operator()(
    const BayesModel& model,
    const bayes::Priors& priors,
    mcmc::State& states,
    std::mt19937_64& rng) const -> void
{
    Gibbs::R(
        *model.genetic(GeneticMode::D),
        *priors.genetic(GeneticMode::D),
        *states.genetic(GeneticMode::D),
        states.residual(),
        rng);
}

auto RR::operator()(
    const BayesModel& model,
    const bayes::Priors& priors,
    mcmc::State& states,
    std::mt19937_64& rng) const -> void
{
    Gibbs::RR(
        *model.genetic(GeneticMode::D),
        *priors.genetic(GeneticMode::D),
        *states.genetic(GeneticMode::D),
        states.residual(),
        rng);
}

auto AT::operator()(
    const BayesModel& model,
    const bayes::Priors& priors,
    mcmc::State& states,
    std::mt19937_64& rng) const -> void
{
    Gibbs::R<Gibbs::policy::AT>(
        *model.genetic(GeneticMode::D),
        *priors.genetic(GeneticMode::D),
        *states.genetic(GeneticMode::D),
        states.residual(),
        rng);
}

}  // namespace gelex::detail::DominantSampler
