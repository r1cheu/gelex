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

#ifndef GELEX_BAYES_GENETIC_DETAIL_PRIOR_SUPPORT_H_
#define GELEX_BAYES_GENETIC_DETAIL_PRIOR_SUPPORT_H_

#include <utility>

#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/parameter.h"

namespace gelex::detail
{

template <MixtureWeightUpdate Update, typename T, typename Prior>
auto make_parameter(T initial, Prior prior)
{
    if constexpr (Update == MixtureWeightUpdate::Disabled)
    {
        return FixedParameter<T>{.initial = std::move(initial)};
    }
    else
    {
        return Parameter<T, Prior>{
            .initial = std::move(initial), .prior = std::move(prior)};
    }
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_DETAIL_PRIOR_SUPPORT_H_
