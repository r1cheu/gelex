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

#include "gelex/bayes/genetic/prior_state_values.h"

#include <cstddef>
#include <ranges>
#include <span>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/bayes/genetic/parameters.h"
#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"

namespace gelex::bayes
{

MixtureState::MixtureState(
    const MixtureProportion& proportion_parameter,
    Eigen::Index num_markers)
    : proportion(proportion_parameter.initial_value()),
      assignment(num_markers, proportion_parameter.initial_value().size())
{
}

auto MixtureState::visit(FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);

    visitor.on(
        "proportion", proportion, FieldFlag::checkpoint | FieldFlag::trace);
    std::vector<std::string> proportion_names;
    proportion_names.reserve(static_cast<std::size_t>(proportion.size()));
    for (Eigen::Index i = 0; i < proportion.size(); ++i)
    {
        proportion_names.push_back(fmt::format("π[{}]", i));
    }
    visitor.on(
        "proportion_names",
        std::span<const std::string>{proportion_names},
        FieldFlag::report);
    visitor.on(
        "assignment", assignment, FieldFlag::checkpoint | FieldFlag::trace);
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

auto ComponentState::visit(FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("gebv_var", gebv_var, FieldFlag::checkpoint | FieldFlag::trace);
    std::vector<std::string> gebv_var_names;
    gebv_var_names.reserve(static_cast<std::size_t>(gebv_var.size()));
    for (Eigen::Index i = 0; i < gebv_var.size(); ++i)
    {
        gebv_var_names.push_back(fmt::format("σ²_component[{}]", i));
    }
    visitor.on(
        "gebv_var_names",
        std::span<const std::string>{gebv_var_names},
        FieldFlag::report);
    for (auto [i, value] : std::views::enumerate(gebv))
    {
        visitor.on(fmt::format("gebv_{}", i), value, FieldFlag::checkpoint);
    }
}

DominanceSignState::DominanceSignState(
    const ProbabilityParameter& probability,
    Eigen::Index num_markers)
    : positive_probability(probability.initial_value()), sign(num_markers, 2)
{
}

auto DominanceSignState::visit(FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on(
        "positive_probability",
        positive_probability,
        FieldFlag::checkpoint | FieldFlag::trace);
    visitor.on("positive_probability_name", "p_s", FieldFlag::report);
    visitor.on("sign", sign, FieldFlag::checkpoint | FieldFlag::trace);
}

}  // namespace gelex::bayes
