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
#include <vector>

#include <Eigen/Core>

#include "gelex/infra/record_visitor.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/prior_state.h"
#include "gelex/model/bayes/state_record_set.h"
#include "gelex/types/fixed_effects.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class BayesModel;

namespace bayes
{

class BayesPrior;

struct FixedState
{
    explicit FixedState(const FixedEffect& effect);
    explicit FixedState(Eigen::VectorXd coeffs);

    auto visit_records(StateRecordSet set, infra::RecordSink& sink) const
        -> void;
    auto visit_records(StateRecordSet set, infra::MutableRecordSink& sink)
        -> void;

    Eigen::VectorXd coeffs;
};

struct RandomState
{
    RandomState(const RandomEffect& effect, double variance);
    RandomState(const RandomEffect& effect, const VarianceParameter& parameter);
    RandomState(Eigen::VectorXd coeffs, double variance);

    auto visit_records(StateRecordSet set, infra::RecordSink& sink) const
        -> void;
    auto visit_records(StateRecordSet set, infra::MutableRecordSink& sink)
        -> void;

    Eigen::VectorXd coeffs;
    double variance{0.0};
};

struct GeneticState
{
    GeneticState(
        GeneticMode type,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

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

struct ResidualState
{
    auto visit_records(StateRecordSet set, infra::RecordSink& sink) const
        -> void;
    auto visit_records(StateRecordSet set, infra::MutableRecordSink& sink)
        -> void;

    Eigen::VectorXd y_adj;
    double variance{0.0};
};

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

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_STATE_H_
