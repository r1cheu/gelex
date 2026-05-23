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

#include "gelex/infra/record_visitor.h"
#include "gelex/model/bayes/prior_specs.h"
#include "gelex/model/bayes/state_capabilities.h"

namespace gelex::bayes
{

ProportionState::ProportionState(
    const ProportionSpec& spec,
    Eigen::Index num_markers)
{
    count = Eigen::VectorXi::Zero(spec.size());
    count(0) = static_cast<int>(num_markers);
    assignment = Eigen::VectorXi::Zero(num_markers);
    value = spec.initial_value();
    update = spec.update();
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

auto VarianceStateCap::visit_sample_records(infra::RecordVisitor& visitor) const
    -> void
{
    for (auto [i, value] : std::views::enumerate(variance()))
    {
        visitor.emit("variance", static_cast<std::size_t>(i), "value", value);
    }
}

auto VarianceStateCap::visit_checkpoint_records(
    infra::RecordVisitor& visitor) const -> void
{
    for (auto [i, value] : std::views::enumerate(variance()))
    {
        visitor.emit("variance", static_cast<std::size_t>(i), "value", value);
    }
}

auto VarianceStateCap::visit_checkpoint_records(
    infra::MutableRecordVisitor& visitor) -> void
{
    for (auto [i, value] : std::views::enumerate(variance()))
    {
        visitor.emit("variance", static_cast<std::size_t>(i), "value", value);
    }
}

auto ProportionStateCap::visit_sample_records(
    infra::RecordVisitor& visitor) const -> void
{
    for (auto [i, state] : std::views::enumerate(proportion()))
    {
        const auto slot = static_cast<std::size_t>(i);
        visitor.emit("proportion", slot, "assignment", state.assignment);
        if (state.update == ProportionUpdate::sampled)
        {
            visitor.emit("proportion", slot, "value", state.value);
        }
    }
}

auto ProportionStateCap::visit_checkpoint_records(
    infra::RecordVisitor& visitor) const -> void
{
    for (auto [i, state] : std::views::enumerate(proportion()))
    {
        const auto slot = static_cast<std::size_t>(i);
        visitor.emit("proportion", slot, "assignment", state.assignment);
        visitor.emit("proportion", slot, "count", state.count);
        visitor.emit("proportion", slot, "value", state.value);
    }
}

auto ProportionStateCap::visit_checkpoint_records(
    infra::MutableRecordVisitor& visitor) -> void
{
    for (auto [i, state] : std::views::enumerate(proportion()))
    {
        const auto slot = static_cast<std::size_t>(i);
        visitor.emit("proportion", slot, "assignment", state.assignment);
        visitor.emit("proportion", slot, "count", state.count);
        visitor.emit("proportion", slot, "value", state.value);
    }
}

auto ComponentStateCap::visit_sample_records(
    infra::RecordVisitor& visitor) const -> void
{
    for (auto [i, state] : std::views::enumerate(component()))
    {
        visitor.emit(
            "component",
            static_cast<std::size_t>(i),
            "gebv_var",
            state.gebv_var);
    }
}

auto ComponentStateCap::visit_checkpoint_records(
    infra::RecordVisitor& visitor) const -> void
{
    for (auto [i, state] : std::views::enumerate(component()))
    {
        const auto slot = static_cast<std::size_t>(i);
        visitor.emit("component", slot, "gebv_var", state.gebv_var);
        for (auto [k, value] : std::views::enumerate(state.gebv))
        {
            visitor.emit("component", slot, fmt::format("gebv_{}", k), value);
        }
    }
}

auto ComponentStateCap::visit_checkpoint_records(
    infra::MutableRecordVisitor& visitor) -> void
{
    for (auto [i, state] : std::views::enumerate(component()))
    {
        const auto slot = static_cast<std::size_t>(i);
        visitor.emit("component", slot, "gebv_var", state.gebv_var);
        for (auto [k, value] : std::views::enumerate(state.gebv))
        {
            visitor.emit("component", slot, fmt::format("gebv_{}", k), value);
        }
    }
}

}  // namespace gelex::bayes
