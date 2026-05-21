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

#include "gelex/model/bayes/legacy_prior.h"

#include <optional>

#include <Eigen/Core>

#include "gelex/model/bayes/legacy_bayes_policy.h"
#include "gelex/model/bayes/legacy_prior_constants.h"

namespace gelex::bayes
{

auto OldVarianceSpec::make(double phenotype_variance) -> OldVarianceSpec
{
    constexpr double kResidualVarianceProportion = 0.3;
    return OldVarianceSpec{
        .scope = MarkerVarianceScope::per_effect,
        .init = kResidualVarianceProportion * phenotype_variance,
        .prior = ScaledInvChiSqPrior{
            prior_constants::RESIDUAL_VARIANCE_SHAPE,
            prior_constants::RESIDUAL_VARIANCE_SCALE,
        },
    };
}

auto OldVarianceSpec::make(double init, MarkerVarianceScope scope)
    -> OldVarianceSpec
{
    return OldVarianceSpec{
        .scope = scope,
        .init = init,
        .prior = ScaledInvChiSqPrior{
            prior_constants::MARKER_VARIANCE_SHAPE,
            prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER * init,
        },
    };
}

auto Mixture::make(const BayesPolicy& policy, bool estimate_pi)
    -> std::optional<Mixture>
{
    if (!policy.mixture)
    {
        return std::nullopt;
    }
    const auto& m = *policy.mixture;
    return Mixture{
        m.strategy,
        CategoricalSpec{
            m.default_proportion,
            OldDirichletPrior{
                Eigen::VectorXi::Ones(m.default_proportion.size())},
            estimate_pi,
        },
    };
}

}  // namespace gelex::bayes
