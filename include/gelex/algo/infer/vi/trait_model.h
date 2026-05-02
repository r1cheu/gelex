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

#ifndef GELEX_ALGO_INFER_VI_TRAIT_MODEL_H_
#define GELEX_ALGO_INFER_VI_TRAIT_MODEL_H_

#include "gelex/algo/infer/vi/updaters/additive.h"
#include "gelex/algo/infer/vi/updaters/common.h"
#include "gelex/algo/infer/vi/updaters/dominant.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"

namespace gelex::vi
{

template <typename T>
concept Updater = requires(
    const T& op,
    const BayesModel& model,
    const bayes::Priors& priors,
    vi::State& state) {
    { op(model, priors, state) } -> std::same_as<void>;
};

template <Updater... Updaters>
class TraitModel
{
   public:
    auto operator()(
        const BayesModel& model,
        const bayes::Priors& priors,
        vi::State& state) const -> void
    {
        std::apply(
            [&](const auto&... updater)
            { (updater(model, priors, state), ...); },
            updaters_);
    }

   private:
    std::tuple<Updaters...> updaters_{};
};

template <Updater... Updaters>
using TraitBasicDefault = TraitModel<
    ::gelex::vi::detail::Fixed,
    ::gelex::vi::detail::Random,
    Updaters...,
    ::gelex::vi::detail::Residual>;

using RR = TraitBasicDefault<::gelex::vi::detail::AdditiveRR>;
using RRd = TraitBasicDefault<
    ::gelex::vi::detail::AdditiveRR,
    ::gelex::vi::detail::DominantRR>;

}  // namespace gelex::vi

#endif  // GELEX_ALGO_INFER_VI_TRAIT_MODEL_H_
