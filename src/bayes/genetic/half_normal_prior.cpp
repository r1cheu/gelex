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

#include "gelex/bayes/genetic/half_normal_prior.h"

#include <array>
#include <utility>

#include <Eigen/Core>

#include "gelex/bayes/genetic/half_normal_prior_state.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/exception.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

JointHalfNormalMixturePrior::JointHalfNormalMixturePrior(
    JointSharedMarkerVariance variance,
    MixtureProportion proportion,
    ProbabilityParameter dominance_positive_probability)
    : variances_(std::move(variance)),
      proportion_(std::move(proportion)),
      dominance_positive_probability_(std::move(dominance_positive_probability))
{
    if (this->proportion().size() != 4)
    {
        throw GelexException(
            "JointHalfNormalMixturePrior: proportion must have size 4");
    }
}

auto JointHalfNormalMixturePrior::variance(GeneticMode mode)
    -> SharedMarkerVariance&
{
    return variances_.variance(mode);
}

auto JointHalfNormalMixturePrior::variance(GeneticMode mode) const
    -> const SharedMarkerVariance&
{
    return variances_.variance(mode);
}

auto JointHalfNormalMixturePrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    variances().visit(visitor);
    proportion().visit(visitor);
    dominance_positive_probability().visit(visitor);
}

auto JointHalfNormalMixturePrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const -> JointHalfNormalMixtureState
{
    return JointHalfNormalMixtureState{
        std::array{
            variance(GeneticMode::A).parameter().initial_value(),
            variance(GeneticMode::D).parameter().initial_value()},
        proportion(),
        dominance_positive_probability(),
        num_markers};
}

}  // namespace gelex::bayes
