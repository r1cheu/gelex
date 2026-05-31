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

#ifndef GELEX_BAYES_STATE_H_
#define GELEX_BAYES_STATE_H_

#include <algorithm>
#include <cstddef>
#include <span>
#include <string_view>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/prior_state.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class BayesModel;

namespace infra
{
class FieldVisitor;
}

namespace bayes
{

class BayesPrior;
class RandomPrior;

struct FixedState
{
    static constexpr std::string_view name = "fixed";

    explicit FixedState(const FixedDesign& design);
    explicit FixedState(Eigen::VectorXd coeffs);

    auto visit(infra::FieldVisitor& visitor) -> void;

    Eigen::VectorXd coeffs;
};

struct RandomState
{
    static constexpr std::string_view name = "random";

    RandomState(const RandomDesign& design, const RandomPrior& prior);

    auto visit(infra::FieldVisitor& visitor) -> void;

    Eigen::VectorXd coeffs;
    double variance{0.0};
};

struct ResidualState
{
    static constexpr std::string_view name = "residual";

    auto visit(infra::FieldVisitor& visitor) -> void;

    Eigen::VectorXd y_adj;
    double variance{0.0};
    Eigen::VectorXd old_y_adj;
};

struct GeneticState
{
    static constexpr std::string_view name = "genetic";

    GeneticState(
        GeneticMode type,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto visit(infra::FieldVisitor& visitor) -> void;

    GeneticMode type;
    Eigen::VectorXd coeffs;
    Eigen::VectorXd u;
    double variance{};
    double heritability{};
};

class SingleGeneticBlockState
{
   public:
    static constexpr std::string_view name = "single";

    SingleGeneticBlockState(
        const GeneticDesign& design,
        const SingleGeneticPrior& prior);

    auto mode() const -> GeneticMode { return state_.type; }
    auto contains(GeneticMode target_mode) const -> bool
    {
        return target_mode == mode();
    }

    auto state() -> GeneticState& { return state_; }
    auto state() const -> const GeneticState& { return state_; }

    auto prior_state() -> SingleGeneticPriorState& { return prior_state_; }
    auto prior_state() const -> const SingleGeneticPriorState&
    {
        return prior_state_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void;

   private:
    GeneticState state_;
    SingleGeneticPriorState prior_state_;
};

class JointGeneticBlockState
{
   public:
    static constexpr std::string_view name = "joint";

    JointGeneticBlockState(
        const GeneticDesign& additive,
        const GeneticDesign& dominance,
        const JointGeneticPrior& prior);

    auto contains(GeneticMode mode) const -> bool;
    auto state(GeneticMode mode) -> GeneticState&;
    auto state(GeneticMode mode) const -> const GeneticState&;

    auto prior_state() -> JointGeneticPriorState& { return prior_state_; }
    auto prior_state() const -> const JointGeneticPriorState&
    {
        return prior_state_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void;

   private:
    GeneticState additive_;
    GeneticState dominance_;
    JointGeneticPriorState prior_state_;
};

using GeneticPriorBlockState
    = std::variant<SingleGeneticBlockState, JointGeneticBlockState>;

}  // namespace bayes

class BayesState
{
   public:
    static constexpr std::string_view name = "state";

    BayesState(const BayesModel& model, const bayes::BayesPrior& prior);

    auto fixed() -> bayes::FixedState& { return fixed_; }
    auto fixed() const -> const bayes::FixedState& { return fixed_; }

    auto random() -> std::vector<bayes::RandomState>& { return random_; }
    auto random() const -> const std::vector<bayes::RandomState>&
    {
        return random_;
    }

    auto genetics() -> std::vector<bayes::GeneticPriorBlockState>&
    {
        return genetics_;
    }
    auto genetics() const -> const std::vector<bayes::GeneticPriorBlockState>&
    {
        return genetics_;
    }

    auto genetic(GeneticMode type) -> bayes::GeneticState*;
    auto genetic(GeneticMode type) const -> const bayes::GeneticState*;

    auto genetic_block_for(GeneticMode type) -> bayes::GeneticPriorBlockState*;
    auto genetic_block_for(GeneticMode type) const
        -> const bayes::GeneticPriorBlockState*;

    auto residual() -> bayes::ResidualState& { return residual_; }
    auto residual() const -> const bayes::ResidualState& { return residual_; }

    auto compute_heritability() -> void;
    auto visit(infra::FieldVisitor& visitor) -> void;

   private:
    bayes::FixedState fixed_;
    std::vector<bayes::RandomState> random_;
    std::vector<bayes::GeneticPriorBlockState> genetics_;
    bayes::ResidualState residual_;
};

}  // namespace gelex

#endif  // GELEX_BAYES_STATE_H_
