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

#include "gelex/model/bayes/prior_state.h"

#include <cstddef>
#include <ranges>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/model/bayes/prior_parameters.h"

namespace gelex::bayes
{

MixtureAssignmentState::MixtureAssignmentState(
    Eigen::Index num_components,
    Eigen::Index num_markers)
{
    count = Eigen::VectorXi::Zero(num_components);
    count(0) = static_cast<int>(num_markers);
    assignment = Eigen::VectorXi::Zero(num_markers);
}

auto MixtureAssignmentState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on(
        "assignment", assignment, FieldFlag::checkpoint | FieldFlag::trace);
    visitor.on("count", count, FieldFlag::checkpoint);
}

SampledMixtureProportionState::SampledMixtureProportionState(
    const SampledMixtureProportion& proportion,
    Eigen::Index num_markers)
    : assignment(proportion.size(), num_markers),
      value(proportion.parameter().initial_value())
{
}

auto SampledMixtureProportionState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    assignment.visit(visitor);
    visitor.on("value", value, FieldFlag::checkpoint | FieldFlag::trace);
}

ComponentState::ComponentState(
    Eigen::Index num_components,
    Eigen::Index num_individuals)
{
    gebv.assign(
        static_cast<std::size_t>(num_components),
        Eigen::VectorXd::Zero(num_individuals));
    gebv_var = Eigen::VectorXd::Zero(num_components);
}

auto ComponentState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("gebv_var", gebv_var, FieldFlag::checkpoint | FieldFlag::trace);
    for (auto [i, value] : std::views::enumerate(gebv))
    {
        visitor.on(fmt::format("gebv_{}", i), value, FieldFlag::checkpoint);
    }
}

}  // namespace gelex::bayes
