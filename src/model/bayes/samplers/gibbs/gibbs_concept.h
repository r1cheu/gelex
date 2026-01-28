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

#ifndef GELEX_MODEL_BAYES_SAMPLERS_GIBBS_GIBBS_CONCEPT_H_
#define GELEX_MODEL_BAYES_SAMPLERS_GIBBS_GIBBS_CONCEPT_H_

#include <concepts>

namespace gelex::bayes
{
struct AdditiveEffect;
struct AdditiveState;
struct DominantEffect;
struct DominantState;
}  // namespace gelex::bayes

namespace gelex::detail::Gibbs
{

template <typename E, typename S>
concept IsValidEffectStatePair
    = (std::same_as<std::remove_cvref_t<E>, bayes::AdditiveEffect>
       && std::same_as<std::remove_cvref_t<S>, bayes::AdditiveState>)
      || (std::same_as<std::remove_cvref_t<E>, bayes::DominantEffect>
          && std::same_as<std::remove_cvref_t<S>, bayes::DominantState>);

}  // namespace gelex::detail::Gibbs

#endif  // GELEX_MODEL_BAYES_SAMPLERS_GIBBS_GIBBS_CONCEPT_H_
