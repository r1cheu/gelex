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
#include "gelex/model/effects.h"

namespace gelex
{

inline auto create_prior_strategy(BayesAlphabet type)
    -> std::optional<PriorSetter>
{
    using enum PriorType;
    using enum VarianceScope;
    using bt = BayesAlphabet;

    auto non_mixture = [](VarianceScope scope, bool has_dominant) -> PriorSpec
    {
        PriorSpec spec{{NonMixture, scope, false}, std::nullopt};
        if (has_dominant)
        {
            spec.dominant = {NonMixture, scope, false};
        }
        return spec;
    };

    auto pi_mixture
        = [](VarianceScope scope, bool estimate, bool has_dominant) -> PriorSpec
    {
        PriorSpec spec{{PiMixture, scope, estimate}, std::nullopt};
        if (has_dominant)
        {
            spec.dominant = {PiMixture, scope, estimate};
        }
        return spec;
    };

    auto scale_mixture = [](bool has_dominant) -> PriorSpec
    {
        PriorSpec spec{{ScaleMixture, Shared, true}, std::nullopt};
        if (has_dominant)
        {
            spec.dominant = {ScaleMixture, Shared, true};
        }
        return spec;
    };

    switch (type)
    {
        case bt::A:
            return PriorSetter(non_mixture(PerMarker, false));
        case bt::Ad:
            return PriorSetter(non_mixture(PerMarker, true));
        case bt::RR:
            return PriorSetter(non_mixture(PerMarker, false));
        case bt::RRd:
            return PriorSetter(non_mixture(PerMarker, true));
        case bt::B:
            return PriorSetter(pi_mixture(PerMarker, false, false));
        case bt::Bpi:
            return PriorSetter(pi_mixture(PerMarker, true, false));
        case bt::Bd:
            return PriorSetter(pi_mixture(PerMarker, false, true));
        case bt::Bdpi:
            return PriorSetter(pi_mixture(PerMarker, true, true));
        case bt::C:
            return PriorSetter(pi_mixture(Shared, false, false));
        case bt::Cpi:
            return PriorSetter(pi_mixture(Shared, true, false));
        case bt::Cd:
            return PriorSetter(pi_mixture(Shared, false, true));
        case bt::Cdpi:
            return PriorSetter(pi_mixture(Shared, true, true));
        case bt::R:
            return PriorSetter(scale_mixture(false));
        case bt::Rd:
            return PriorSetter(scale_mixture(true));
        default:
            return std::nullopt;
    }
}

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_PRIOR_STRATEGIES_H_
