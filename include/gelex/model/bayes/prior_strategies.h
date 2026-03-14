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

#ifndef GELEX_MODEL_BAYES_PRIOR_STRATEGIES_H_
#define GELEX_MODEL_BAYES_PRIOR_STRATEGIES_H_

#include <optional>

#include "gelex/model/bayes/prior_strategy.h"
#include "gelex/types/bayes_method.h"

namespace gelex
{

inline auto create_prior_strategy(const BayesMethodConfig& method)
    -> std::optional<PriorSetter>
{
    using enum PriorType;
    using enum VarianceScope;
    using enum GeneticKind;

    auto non_mixture = [](VarianceScope scope, bool has_dominant) -> PriorSpec
    {
        PriorSpec spec{{{Add, {NonMixture, scope, false}}}};
        if (has_dominant)
        {
            spec.genetics.push_back({Dom, {NonMixture, scope, false}});
        }
        return spec;
    };

    auto pi_mixture
        = [](VarianceScope scope, bool estimate, bool has_dominant) -> PriorSpec
    {
        PriorSpec spec{{{Add, {PiMixture, scope, estimate}}}};
        if (has_dominant)
        {
            spec.genetics.push_back({Dom, {PiMixture, scope, estimate}});
        }
        return spec;
    };

    auto scale_mixture = [](bool has_dominant, bool asymmetric) -> PriorSpec
    {
        PriorSpec spec{{{Add, {ScaleMixture, Shared, true}}}};
        if (has_dominant)
        {
            spec.genetics.push_back(
                {Dom, {ScaleMixture, Shared, true, asymmetric}});
        }
        return spec;
    };

    bool dom = method.dominance;

    switch (method.base)
    {
        case BayesBase::A:
        case BayesBase::RR:
            return PriorSetter(non_mixture(PerMarker, dom));
        case BayesBase::B:
            return PriorSetter(pi_mixture(PerMarker, method.estimate_pi, dom));
        case BayesBase::C:
            return PriorSetter(pi_mixture(Shared, method.estimate_pi, dom));
        case BayesBase::R:
            return PriorSetter(scale_mixture(dom, method.asymmetric));
    }
    return std::nullopt;
}

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_PRIOR_STRATEGIES_H_
