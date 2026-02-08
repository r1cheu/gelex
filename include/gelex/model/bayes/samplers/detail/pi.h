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

#ifndef GELEX_MODEL_BAYES_SAMPLERS_DETAIL_PI_H_
#define GELEX_MODEL_BAYES_SAMPLERS_DETAIL_PI_H_

#include <random>

#include "gelex/model/bayes/model.h"

namespace gelex::detail
{

/**
 * @brief Sampler for pi parameter in Bayesian models
 *
 * This sampler handles the Dirichlet update for pi parameters.
 * Currently supports 2-component pi (like BayesCpi, BayesBpi),
 * but designed to be extensible for multi-component pi (like BayesR).
 */
namespace AdditiveSampler
{
struct Pi
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};
}  // namespace AdditiveSampler

namespace DominantSampler
{
struct Pi
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};
}  // namespace DominantSampler

}  // namespace gelex::detail

#endif  // GELEX_MODEL_BAYES_SAMPLERS_DETAIL_PI_H_
