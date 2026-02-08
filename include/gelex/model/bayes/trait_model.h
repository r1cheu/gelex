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

#ifndef GELEX_MODEL_BAYES_TRAIT_MODEL_H_
#define GELEX_MODEL_BAYES_TRAIT_MODEL_H_

#include <random>

#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/samplers/detail/additive.h"
#include "gelex/model/bayes/samplers/detail/common.h"
#include "gelex/model/bayes/samplers/detail/dominant.h"
#include "gelex/model/bayes/samplers/detail/pi.h"

namespace gelex
{

template <typename T>
concept Sampler = requires(
    const T& op,
    const BayesModel& model,
    BayesState& status,
    std::mt19937_64& rng) {
    { op(model, status, rng) } -> std::same_as<void>;
};

template <Sampler... Samplers>
class TraitModel
{
   public:
    auto operator()(
        const BayesModel& model,
        BayesState& status,
        std::mt19937_64& rng) const -> void
    {
        std::apply(
            [&](const auto&... sampler) { (sampler(model, status, rng), ...); },
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

using BayesRR = TraitBasicDefault<detail::AdditiveSampler::RR>;
using BayesA = TraitBasicDefault<detail::AdditiveSampler::A>;
using BayesB = TraitBasicDefault<detail::AdditiveSampler::B>;
using BayesC = TraitBasicDefault<detail::AdditiveSampler::C>;

using BayesBpi = TraitBasicDefault<
    detail::AdditiveSampler::B,
    detail::AdditiveSampler::Pi>;
using BayesCpi = TraitBasicDefault<
    detail::AdditiveSampler::C,
    detail::AdditiveSampler::Pi>;
using BayesR = TraitBasicDefault<
    detail::AdditiveSampler::R,
    detail::AdditiveSampler::Pi>;

using BayesRRd = TraitBasicDefault<
    detail::AdditiveSampler::RR,
    detail::DominantSampler::RR>;

using BayesAd
    = TraitBasicDefault<detail::AdditiveSampler::A, detail::DominantSampler::A>;

using BayesBd
    = TraitBasicDefault<detail::AdditiveSampler::B, detail::DominantSampler::B>;

using BayesBdpi = TraitBasicDefault<
    detail::AdditiveSampler::B,
    detail::AdditiveSampler::Pi,
    detail::DominantSampler::B,
    detail::DominantSampler::Pi>;

using BayesCd
    = TraitBasicDefault<detail::AdditiveSampler::C, detail::DominantSampler::C>;

using BayesCdpi = TraitBasicDefault<
    detail::AdditiveSampler::C,
    detail::AdditiveSampler::Pi,
    detail::DominantSampler::C,
    detail::DominantSampler::Pi>;

using BayesRd = TraitBasicDefault<
    detail::AdditiveSampler::R,
    detail::AdditiveSampler::Pi,
    detail::DominantSampler::R,
    detail::DominantSampler::Pi>;

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_TRAIT_MODEL_H_
