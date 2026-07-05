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

#include "gelex/bayes/genetic/prior.h"

#include <variant>

#include <Eigen/Core>

#include "gelex/bayes/genetic/prior_state.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::bayes
{

auto mode(const SingleGeneticPrior& prior) -> GeneticMode
{
    return std::visit([](const auto& value) { return value.mode(); }, prior);
}

auto visit(SingleGeneticPrior& prior, infra::FieldVisitor& visitor) -> void
{
    std::visit([&visitor](auto& value) { value.visit(visitor); }, prior);
}

auto visit(JointGeneticPrior& prior, infra::FieldVisitor& visitor) -> void
{
    std::visit([&visitor](auto& value) { value.visit(visitor); }, prior);
}

auto make_state(
    const SingleGeneticPrior& prior,
    Eigen::Index num_markers,
    Eigen::Index num_individuals) -> SingleGeneticPriorState
{
    return std::visit(
        [num_markers,
         num_individuals](const auto& value) -> SingleGeneticPriorState
        { return value.make_state(num_markers, num_individuals); },
        prior);
}

auto make_state(
    const JointGeneticPrior& prior,
    Eigen::Index num_markers,
    Eigen::Index num_individuals) -> JointGeneticPriorState
{
    return std::visit(
        [num_markers,
         num_individuals](const auto& value) -> JointGeneticPriorState
        { return value.make_state(num_markers, num_individuals); },
        prior);
}

}  // namespace gelex::bayes
