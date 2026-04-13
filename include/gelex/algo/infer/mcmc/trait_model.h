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

#ifndef GELEX_ALGO_INFER_MCMC_TRAIT_MODEL_H_
#define GELEX_ALGO_INFER_MCMC_TRAIT_MODEL_H_

#include <random>

#include "gelex/algo/infer/mcmc/samplers/additive.h"
#include "gelex/algo/infer/mcmc/samplers/common.h"
#include "gelex/algo/infer/mcmc/samplers/dominant.h"
#include "gelex/algo/infer/mcmc/samplers/pi.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"

namespace gelex::mcmc
{

template <typename T>
concept Sampler = requires(
    const T& op,
    const BayesModel& model,
    const bayes::Priors& priors,
    mcmc::State& status,
    std::mt19937_64& rng) {
    { op(model, priors, status, rng) } -> std::same_as<void>;
};

template <Sampler... Samplers>
class TraitModel
{
   public:
    auto operator()(
        const BayesModel& model,
        const bayes::Priors& priors,
        mcmc::State& status,
        std::mt19937_64& rng) const -> void
    {
        std::apply(
            [&](const auto&... sampler)
            { (sampler(model, priors, status, rng), ...); },
            samplers_);
    }

   private:
    std::tuple<Samplers...> samplers_{};
};

template <Sampler... Samplers>
using TraitBasicDefault = TraitModel<
    detail::CommonSampler::Fixed,
    detail::CommonSampler::Random,
    Samplers...,
    detail::CommonSampler::Residual>;
using RR = TraitBasicDefault<detail::AdditiveSampler::RR>;
using A = TraitBasicDefault<detail::AdditiveSampler::A>;
using B = TraitBasicDefault<detail::AdditiveSampler::B>;
using C = TraitBasicDefault<detail::AdditiveSampler::C>;

using Bpi = TraitBasicDefault<
    detail::AdditiveSampler::B,
    detail::AdditiveSampler::Pi>;
using Cpi = TraitBasicDefault<
    detail::AdditiveSampler::C,
    detail::AdditiveSampler::Pi>;
using R = TraitBasicDefault<
    detail::AdditiveSampler::R,
    detail::AdditiveSampler::Pi>;

using RRd = TraitBasicDefault<
    detail::AdditiveSampler::RR,
    detail::DominantSampler::RR>;

using Ad
    = TraitBasicDefault<detail::AdditiveSampler::A, detail::DominantSampler::A>;

using Bd
    = TraitBasicDefault<detail::AdditiveSampler::B, detail::DominantSampler::B>;

using Bdpi = TraitBasicDefault<
    detail::AdditiveSampler::B,
    detail::AdditiveSampler::Pi,
    detail::DominantSampler::B,
    detail::DominantSampler::Pi>;

using Cd
    = TraitBasicDefault<detail::AdditiveSampler::C, detail::DominantSampler::C>;

using Cdpi = TraitBasicDefault<
    detail::AdditiveSampler::C,
    detail::AdditiveSampler::Pi,
    detail::DominantSampler::C,
    detail::DominantSampler::Pi>;

using Rd = TraitBasicDefault<
    detail::AdditiveSampler::R,
    detail::AdditiveSampler::Pi,
    detail::DominantSampler::R,
    detail::DominantSampler::Pi>;

using RdAt = TraitBasicDefault<
    detail::AdditiveSampler::R,
    detail::AdditiveSampler::Pi,
    detail::DominantSampler::AT,
    detail::DominantSampler::Pi>;
}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_TRAIT_MODEL_H_
