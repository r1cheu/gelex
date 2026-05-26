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

#ifndef GELEX_MODEL_BAYES_STATE_H_
#define GELEX_MODEL_BAYES_STATE_H_

#include <algorithm>
#include <cstddef>
#include <memory>
#include <span>
#include <string_view>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/infra/record_visitor.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/prior_state.h"
#include "gelex/model/bayes/state_record_set.h"
#include "gelex/types/fixed_effects.h"
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
class BayesPriorV2;
class JointGeneticPrior;
class RandomPrior;
class SingleGeneticPrior;

struct FixedState
{
    static constexpr std::string_view name = "fixed";

    explicit FixedState(const FixedEffect& effect);
    explicit FixedState(Eigen::VectorXd coeffs);

    auto visit(infra::FieldVisitor& visitor) -> void;

    auto visit_records(StateRecordSet set, infra::RecordSink& sink) const
        -> void;
    auto visit_records(StateRecordSet set, infra::MutableRecordSink& sink)
        -> void;

    Eigen::VectorXd coeffs;
};

struct RandomState
{
    static constexpr std::string_view name = "random";

    RandomState(const RandomEffect& effect, double variance);
    RandomState(const RandomEffect& effect, const RandomPrior& prior);
    RandomState(Eigen::VectorXd coeffs, double variance);

    auto visit(infra::FieldVisitor& visitor) -> void;

    auto visit_records(StateRecordSet set, infra::RecordSink& sink) const
        -> void;
    auto visit_records(StateRecordSet set, infra::MutableRecordSink& sink)
        -> void;

    Eigen::VectorXd coeffs;
    double variance{0.0};
};

struct ResidualState
{
    static constexpr std::string_view name = "residual";

    auto visit(infra::FieldVisitor& visitor) -> void;

    auto visit_records(StateRecordSet set, infra::RecordSink& sink) const
        -> void;
    auto visit_records(StateRecordSet set, infra::MutableRecordSink& sink)
        -> void;

    Eigen::VectorXd y_adj;
    double variance{0.0};
};

struct GeneticState
{
    static constexpr std::string_view name = "genetic";

    GeneticState(
        GeneticMode type,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto visit(infra::FieldVisitor& visitor) -> void;

    auto visit_records(StateRecordSet set, infra::RecordSink& sink) const
        -> void;
    auto visit_records(StateRecordSet set, infra::MutableRecordSink& sink)
        -> void;

    GeneticMode type;
    Eigen::VectorXd coeffs;
    Eigen::VectorXd u;
    double variance{};
    double heritability{};
};

class GeneticBlockState
{
   public:
    GeneticBlockState(
        std::vector<GeneticMode> modes,
        std::vector<std::size_t> genetic_indices,
        std::unique_ptr<GeneticPriorState> prior_state);

    auto modes() const -> std::span<const GeneticMode> { return modes_; }
    auto genetic_indices() const -> std::span<const std::size_t>
    {
        return genetic_indices_;
    }

    auto prior_state() -> GeneticPriorState& { return *prior_state_; }
    auto prior_state() const -> const GeneticPriorState&
    {
        return *prior_state_;
    }

    auto contains(GeneticMode mode) const -> bool
    {
        return std::ranges::contains(modes_, mode);
    }
    auto slot(GeneticMode mode) const -> std::size_t;

    auto visit_records(StateRecordSet set, infra::RecordSink& sink) const
        -> void;
    auto visit_records(StateRecordSet set, infra::MutableRecordSink& sink)
        -> void;

   private:
    std::vector<GeneticMode> modes_;
    std::vector<std::size_t> genetic_indices_;
    std::unique_ptr<GeneticPriorState> prior_state_;
};

class SingleGeneticBlockState
{
   public:
    static constexpr std::string_view name = "single";

    SingleGeneticBlockState(
        const GeneticEffect& effect,
        const SingleGeneticPrior& prior);
    SingleGeneticBlockState(
        GeneticState genetic,
        std::unique_ptr<SingleGeneticPriorState> prior_state);

    auto mode() const -> GeneticMode { return state_.type; }
    auto contains(GeneticMode target_mode) const -> bool
    {
        return target_mode == mode();
    }

    auto state() -> GeneticState& { return state_; }
    auto state() const -> const GeneticState& { return state_; }

    auto prior_state() -> SingleGeneticPriorState& { return *prior_state_; }
    auto prior_state() const -> const SingleGeneticPriorState&
    {
        return *prior_state_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void;

   private:
    GeneticState state_;
    std::unique_ptr<SingleGeneticPriorState> prior_state_;
};

class JointGeneticBlockState
{
   public:
    static constexpr std::string_view name = "joint";

    JointGeneticBlockState(
        const GeneticEffect& additive,
        const GeneticEffect& dominance,
        const JointGeneticPrior& prior);
    JointGeneticBlockState(
        GeneticState additive,
        GeneticState dominance,
        std::unique_ptr<JointGeneticPriorState> prior_state);

    auto contains(GeneticMode mode) const -> bool;
    auto state(GeneticMode mode) -> GeneticState&;
    auto state(GeneticMode mode) const -> const GeneticState&;

    auto prior_state() -> JointGeneticPriorState& { return *prior_state_; }
    auto prior_state() const -> const JointGeneticPriorState&
    {
        return *prior_state_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void;

   private:
    GeneticState additive_;
    GeneticState dominance_;
    std::unique_ptr<JointGeneticPriorState> prior_state_;
};

using GeneticPriorBlockState
    = std::variant<SingleGeneticBlockState, JointGeneticBlockState>;

}  // namespace bayes

class BayesState
{
   public:
    BayesState(const BayesModel& model, const bayes::BayesPrior& prior);

    auto fixed() -> bayes::FixedState& { return fixed_; }
    auto fixed() const -> const bayes::FixedState& { return fixed_; }

    auto random() -> std::vector<bayes::RandomState>& { return random_; }
    auto random() const -> const std::vector<bayes::RandomState>&
    {
        return random_;
    }

    auto genetics() -> std::vector<bayes::GeneticState>& { return genetics_; }
    auto genetics() const -> const std::vector<bayes::GeneticState>&
    {
        return genetics_;
    }

    auto genetic(GeneticMode type) -> bayes::GeneticState*
    {
        auto it
            = std::ranges::find(genetics_, type, &bayes::GeneticState::type);
        return it != genetics_.end() ? &*it : nullptr;
    }
    auto genetic(GeneticMode type) const -> const bayes::GeneticState*
    {
        auto it
            = std::ranges::find(genetics_, type, &bayes::GeneticState::type);
        return it != genetics_.end() ? &*it : nullptr;
    }

    auto genetic_blocks() -> std::vector<bayes::GeneticBlockState>&
    {
        return genetic_blocks_;
    }
    auto genetic_blocks() const -> const std::vector<bayes::GeneticBlockState>&
    {
        return genetic_blocks_;
    }

    auto genetic_block(std::size_t index) -> bayes::GeneticBlockState&
    {
        return genetic_blocks_.at(index);
    }
    auto genetic_block(std::size_t index) const
        -> const bayes::GeneticBlockState&
    {
        return genetic_blocks_.at(index);
    }

    auto genetic_block_for(GeneticMode type) -> bayes::GeneticBlockState*
    {
        auto it = std::ranges::find_if(
            genetic_blocks_,
            [type](const bayes::GeneticBlockState& block)
            { return block.contains(type); });
        return it != genetic_blocks_.end() ? &*it : nullptr;
    }
    auto genetic_block_for(GeneticMode type) const
        -> const bayes::GeneticBlockState*
    {
        auto it = std::ranges::find_if(
            genetic_blocks_,
            [type](const bayes::GeneticBlockState& block)
            { return block.contains(type); });
        return it != genetic_blocks_.end() ? &*it : nullptr;
    }

    auto residual() -> bayes::ResidualState& { return residual_; }
    auto residual() const -> const bayes::ResidualState& { return residual_; }

    auto compute_heritability() -> void;
    auto visit_records(bayes::StateRecordSet set, infra::RecordSink& sink) const
        -> void;
    auto visit_records(
        bayes::StateRecordSet set,
        infra::MutableRecordSink& sink) -> void;

   private:
    bayes::FixedState fixed_;
    std::vector<bayes::RandomState> random_;
    std::vector<bayes::GeneticState> genetics_;
    std::vector<bayes::GeneticBlockState> genetic_blocks_;
    bayes::ResidualState residual_;
};

class BayesStateV2
{
   public:
    static constexpr std::string_view name = "state";

    BayesStateV2(const BayesModel& model, const bayes::BayesPriorV2& prior);

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

#endif  // GELEX_MODEL_BAYES_STATE_H_
