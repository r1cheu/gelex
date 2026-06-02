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

#include "gelex/bayes/state.h"

#include <ranges>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/prior_state.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/exception.h"
#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

FixedState::FixedState(const FixedDesign& design)
    : coeffs(Eigen::VectorXd::Zero(design.X.cols()))
{
}

FixedState::FixedState(Eigen::VectorXd coeffs) : coeffs(std::move(coeffs)) {}

auto FixedState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("coeffs", coeffs, FieldFlag::checkpoint | FieldFlag::trace);
}

RandomState::RandomState(const RandomDesign& design, const RandomPrior& prior)
    : coeffs(Eigen::VectorXd::Zero(design.X.cols())),
      variance{prior.initial_value()}
{
}

auto RandomState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("coeffs", coeffs, FieldFlag::checkpoint | FieldFlag::trace);
    visitor.on("variance", variance, FieldFlag::checkpoint | FieldFlag::trace);
}

GeneticState::GeneticState(
    Eigen::Index num_markers,
    Eigen::Index num_individuals)
    : coeffs(Eigen::VectorXd::Zero(num_markers)),
      u(Eigen::VectorXd::Zero(num_individuals))
{
}

auto GeneticState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("coeffs", coeffs, FieldFlag::checkpoint | FieldFlag::trace);
    visitor.on("u", u, FieldFlag::checkpoint);
    visitor.on("variance", variance, FieldFlag::checkpoint | FieldFlag::trace);
    visitor.on(
        "heritability", heritability, FieldFlag::checkpoint | FieldFlag::trace);
}

SingleGeneticBlockState::SingleGeneticBlockState(
    const GeneticDesign& design,
    const SingleGeneticPrior& prior)
    : state_{design.X.cols(), design.X.rows()},
      prior_state_(make_state(prior, design.X.cols(), design.X.rows()))
{
    if (design.type != gelex::bayes::mode(prior))
    {
        throw GelexException(
            fmt::format(
                "SingleGeneticBlockState: design mode {} != prior mode {}",
                design.type,
                gelex::bayes::mode(prior)));
    }
}

auto SingleGeneticBlockState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    state_.visit(visitor);
    auto prior_scope = visitor.scope("prior_state");
    std::visit([&visitor](auto& state) { state.visit(visitor); }, prior_state_);
}

JointGeneticBlockState::JointGeneticBlockState(
    const GeneticDesign& additive,
    const GeneticDesign& dominance,
    const JointGeneticPrior& prior)
    : additive_{additive.X.cols(), additive.X.rows()},
      dominance_{dominance.X.cols(), dominance.X.rows()},
      prior_state_(make_state(prior, additive.X.cols(), additive.X.rows()))
{
    if (additive.type != GeneticMode::A)
    {
        throw GelexException(
            fmt::format(
                "JointGeneticBlockState: additive design has mode {}",
                additive.type));
    }
    if (dominance.type != GeneticMode::D)
    {
        throw GelexException(
            fmt::format(
                "JointGeneticBlockState: dominance design has mode {}",
                dominance.type));
    }
    if (additive.X.rows() != dominance.X.rows()
        || additive.X.cols() != dominance.X.cols())
    {
        throw GelexException(
            fmt::format(
                "JointGeneticBlockState: genetic designs must share shape; "
                "A is {}x{}, D is {}x{}",
                additive.X.rows(),
                additive.X.cols(),
                dominance.X.rows(),
                dominance.X.cols()));
    }
}

auto JointGeneticBlockState::state(GeneticMode mode) -> GeneticState&
{
    switch (mode)
    {
        case GeneticMode::A:
            return additive_;
        case GeneticMode::D:
            return dominance_;
    }
    throw GelexException(
        fmt::format("JointGeneticBlockState: missing genetic mode {}", mode));
}

auto JointGeneticBlockState::state(GeneticMode mode) const
    -> const GeneticState&
{
    switch (mode)
    {
        case GeneticMode::A:
            return additive_;
        case GeneticMode::D:
            return dominance_;
    }
    throw GelexException(
        fmt::format("JointGeneticBlockState: missing genetic mode {}", mode));
}

auto JointGeneticBlockState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    {
        auto additive_scope = visitor.scope("A");
        additive_.visit(visitor);
    }
    {
        auto dominance_scope = visitor.scope("D");
        dominance_.visit(visitor);
    }
    auto prior_scope = visitor.scope("prior_state");
    std::visit([&visitor](auto& state) { state.visit(visitor); }, prior_state_);
}

auto ResidualState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("y_adj", y_adj, FieldFlag::checkpoint);
    visitor.on("variance", variance, FieldFlag::checkpoint | FieldFlag::trace);
}

}  // namespace gelex::bayes

namespace gelex
{

BayesState::BayesState(const BayesModel& model, const bayes::BayesPrior& prior)
    : fixed_(model.fixed()),
      residual_{
          .y_adj = model.phenotype(),
          .variance = prior.residual().initial_value(),
          .old_y_adj = Eigen::VectorXd::Zero(model.phenotype().size())}
{
    const auto& random_designs = model.random();
    random_.reserve(random_designs.size());
    for (const auto& design : random_designs)
    {
        random_.emplace_back(design, prior.random());
    }

    genetics_.reserve(prior.genetics().size());
    for (const auto& block : prior.genetics())
    {
        std::visit(
            [this, &model](const auto& genetic_prior)
            {
                if constexpr (
                    std::is_same_v<
                        std::decay_t<decltype(genetic_prior)>,
                        bayes::SingleGeneticPrior>)
                {
                    const auto prior_mode = bayes::mode(genetic_prior);
                    const auto* design = model.genetic(prior_mode);
                    if (design == nullptr)
                    {
                        throw GelexException(
                            fmt::format(
                                "BayesState: missing genetic design for "
                                "mode {}",
                                prior_mode));
                    }
                    genetics_.emplace_back(
                        bayes::SingleGeneticBlockState{*design, genetic_prior});
                }
                else
                {
                    const auto* additive = model.genetic(GeneticMode::A);
                    if (additive == nullptr)
                    {
                        throw GelexException(
                            "BayesState: missing genetic design for mode A");
                    }

                    const auto* dominance = model.genetic(GeneticMode::D);
                    if (dominance == nullptr)
                    {
                        throw GelexException(
                            "BayesState: missing genetic design for mode D");
                    }

                    genetics_.emplace_back(
                        bayes::JointGeneticBlockState{
                            *additive, *dominance, genetic_prior});
                }
            },
            block);
    }
}

auto BayesState::compute_heritability() -> void
{
    double total_variance = residual_.variance;
    for (const auto& state : random_)
    {
        total_variance += state.variance;
    }

    for (const auto& block : genetics_)
    {
        std::visit(
            [&total_variance](const auto& value)
            {
                if constexpr (
                    std::is_same_v<
                        std::decay_t<decltype(value)>,
                        bayes::SingleGeneticBlockState>)
                {
                    total_variance += value.state().variance;
                }
                else
                {
                    total_variance += value.state(GeneticMode::A).variance;
                    total_variance += value.state(GeneticMode::D).variance;
                }
            },
            block);
    }

    for (auto& block : genetics_)
    {
        std::visit(
            [total_variance](auto& value)
            {
                if constexpr (
                    std::is_same_v<
                        std::decay_t<decltype(value)>,
                        bayes::SingleGeneticBlockState>)
                {
                    value.state().heritability
                        = value.state().variance / total_variance;
                }
                else
                {
                    value.state(GeneticMode::A).heritability
                        = value.state(GeneticMode::A).variance / total_variance;
                    value.state(GeneticMode::D).heritability
                        = value.state(GeneticMode::D).variance / total_variance;
                }
            },
            block);
    }
}

auto BayesState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    fixed_.visit(visitor);
    for (auto [i, state] : std::views::enumerate(random_))
    {
        auto state_scope = visitor.scope(fmt::format("random_{}", i));
        visitor.on(
            "coeffs", state.coeffs, FieldFlag::checkpoint | FieldFlag::trace);
        visitor.on(
            "variance",
            state.variance,
            FieldFlag::checkpoint | FieldFlag::trace);
    }
    for (auto [i, block] : std::views::enumerate(genetics_))
    {
        auto block_scope = visitor.scope(fmt::format("genetic_{}", i));
        std::visit([&visitor](auto& value) { value.visit(visitor); }, block);
    }
    residual_.visit(visitor);
}

}  // namespace gelex
