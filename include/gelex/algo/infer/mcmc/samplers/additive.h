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

#ifndef GELEX_ALGO_INFER_MCMC_SAMPLERS_ADDITIVE_H_
#define GELEX_ALGO_INFER_MCMC_SAMPLERS_ADDITIVE_H_

#include <random>

#include "gelex/model/bayes/model.h"

namespace gelex::bayes
{
class Priors;
}

namespace gelex::detail::AdditiveSampler
{

struct A
{
    auto operator()(
        const BayesModel& model,
        const bayes::Priors& priors,
        mcmc::State& states,
        std::mt19937_64& rng) const -> void;
};

struct B
{
    auto operator()(
        const BayesModel& model,
        const bayes::Priors& priors,
        mcmc::State& states,
        std::mt19937_64& rng) const -> void;
};

struct C
{
    auto operator()(
        const BayesModel& model,
        const bayes::Priors& priors,
        mcmc::State& states,
        std::mt19937_64& rng) const -> void;
};

struct R
{
    auto operator()(
        const BayesModel& model,
        const bayes::Priors& priors,
        mcmc::State& states,
        std::mt19937_64& rng) const -> void;
};

struct RR
{
    auto operator()(
        const BayesModel& model,
        const bayes::Priors& priors,
        mcmc::State& states,
        std::mt19937_64& rng) const -> void;
};

}  // namespace gelex::detail::AdditiveSampler

#endif  // GELEX_ALGO_INFER_MCMC_SAMPLERS_ADDITIVE_H_
