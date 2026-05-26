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
#include <string>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/infra/record_visitor.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/state_capabilities.h"

namespace gelex::bayes
{

ProportionState::ProportionState(
    const MixtureProportion& proportion,
    Eigen::Index num_markers)
    : value(proportion.parameter().initial_value()), update(proportion.update())
{
    count = Eigen::VectorXi::Zero(proportion.size());
    count(0) = static_cast<int>(num_markers);
    assignment = Eigen::VectorXi::Zero(num_markers);
}

auto ProportionState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on(
        "assignment", assignment, FieldFlag::checkpoint | FieldFlag::trace);
    visitor.on("count", count, FieldFlag::checkpoint);
    visitor.on("value", value, FieldFlag::checkpoint | FieldFlag::trace);
    visitor.on("update", update, FieldFlag::report);
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

auto VarianceStateCap::visit_records(StateRecordSet, infra::RecordSink& sink)
    const -> void
{
    for (auto [i, value] : std::views::enumerate(variance()))
    {
        sink.emit("variance", static_cast<std::size_t>(i), "value", value);
    }
}

auto VarianceStateCap::visit_records(
    StateRecordSet set,
    infra::MutableRecordSink& sink) -> void
{
    if (set != StateRecordSet::checkpoint)
    {
        throw GelexException(
            "VarianceStateCap: mutable visit_records requires checkpoint set");
    }
    for (auto [i, value] : std::views::enumerate(variance()))
    {
        sink.emit("variance", static_cast<std::size_t>(i), "value", value);
    }
}

auto ProportionStateCap::visit_records(
    StateRecordSet set,
    infra::RecordSink& sink) const -> void
{
    for (auto [i, state] : std::views::enumerate(proportion()))
    {
        const auto slot = static_cast<std::size_t>(i);
        switch (set)
        {
            case StateRecordSet::sample:
                sink.emit("proportion", slot, "assignment", state.assignment);
                if (state.update == UpdatePolicy::sampled)
                {
                    sink.emit("proportion", slot, "value", state.value);
                }
                break;
            case StateRecordSet::checkpoint:
                sink.emit("proportion", slot, "assignment", state.assignment);
                sink.emit("proportion", slot, "count", state.count);
                sink.emit("proportion", slot, "value", state.value);
                break;
        }
    }
}

auto ProportionStateCap::visit_records(
    StateRecordSet set,
    infra::MutableRecordSink& sink) -> void
{
    if (set != StateRecordSet::checkpoint)
    {
        throw GelexException(
            "ProportionStateCap: mutable visit_records requires "
            "checkpoint set");
    }
    for (auto [i, state] : std::views::enumerate(proportion()))
    {
        const auto slot = static_cast<std::size_t>(i);
        sink.emit("proportion", slot, "assignment", state.assignment);
        sink.emit("proportion", slot, "count", state.count);
        sink.emit("proportion", slot, "value", state.value);
    }
}

auto ComponentStateCap::visit_records(
    StateRecordSet set,
    infra::RecordSink& sink) const -> void
{
    for (auto [i, state] : std::views::enumerate(component()))
    {
        const auto slot = static_cast<std::size_t>(i);
        switch (set)
        {
            case StateRecordSet::sample:
                sink.emit("component", slot, "gebv_var", state.gebv_var);
                break;
            case StateRecordSet::checkpoint:
                sink.emit("component", slot, "gebv_var", state.gebv_var);
                for (auto [k, value] : std::views::enumerate(state.gebv))
                {
                    sink.emit(
                        "component", slot, fmt::format("gebv_{}", k), value);
                }
                break;
        }
    }
}

auto ComponentStateCap::visit_records(
    StateRecordSet set,
    infra::MutableRecordSink& sink) -> void
{
    if (set != StateRecordSet::checkpoint)
    {
        throw GelexException(
            "ComponentStateCap: mutable visit_records requires "
            "checkpoint set");
    }
    for (auto [i, state] : std::views::enumerate(component()))
    {
        const auto slot = static_cast<std::size_t>(i);
        sink.emit("component", slot, "gebv_var", state.gebv_var);
        for (auto [k, value] : std::views::enumerate(state.gebv))
        {
            sink.emit("component", slot, fmt::format("gebv_{}", k), value);
        }
    }
}

}  // namespace gelex::bayes
