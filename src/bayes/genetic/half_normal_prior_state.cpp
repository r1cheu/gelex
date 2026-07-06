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

#include "gelex/bayes/genetic/half_normal_prior_state.h"

#include <array>
#include <utility>

#include <fmt/format.h>

#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/genetic/prior_state_values.h"
#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::bayes
{

JointHalfNormalMixtureState::JointHalfNormalMixtureState(
    std::array<double, 2> variances,
    const MixtureProportion& proportion,
    const ProbabilityParameter& dominance_positive_probability,
    Eigen::Index num_markers)
    : variances_(std::move(variances)),
      mixture_(proportion, num_markers),
      dominance_sign_(dominance_positive_probability, num_markers)
{
}

auto JointHalfNormalMixtureState::variance(GeneticMode mode) -> double&
{
    return variances_[std::to_underlying(mode)];
}

auto JointHalfNormalMixtureState::variance(GeneticMode mode) const -> const
    double&
{
    return variances_[std::to_underlying(mode)];
}

auto JointHalfNormalMixtureState::visit(FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    for (const auto mode : modes)
    {
        auto mode_scope = visitor.scope(fmt::format("{}", mode));
        visitor.on(
            "variance",
            variance(mode),
            FieldFlag::checkpoint | FieldFlag::trace);
        visitor.on("variance_name", "σ²_marker", FieldFlag::report);
    }
    mixture_.visit(visitor);
    dominance_sign_.visit(visitor);
}

}  // namespace gelex::bayes
