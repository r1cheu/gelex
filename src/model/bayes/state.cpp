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

#include "gelex/model/bayes/state.h"

#include <cstddef>
#include <memory>
#include <ranges>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/infra/record_visitor.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

namespace
{

template <bool Mutable>
class PrefixedRecordSink final : public infra::BasicRecordSink<Mutable>
{
   public:
    using Base = infra::BasicRecordSink<Mutable>;
    using VectorXfRecord = typename Base::template RecordType<Eigen::VectorXf>;
    using VectorXdRecord = typename Base::template RecordType<Eigen::VectorXd>;
    using VectorXiRecord = typename Base::template RecordType<Eigen::VectorXi>;
    using DoubleRecord = std::conditional_t<Mutable, double&, const double&>;

    PrefixedRecordSink(std::string prefix, Base& inner)
        : prefix_(std::move(prefix)), inner_(inner)
    {
    }

    auto visit(std::string_view path, VectorXfRecord value) -> void override
    {
        inner_.visit(make_path(path), value);
    }
    auto visit(std::string_view path, VectorXdRecord value) -> void override
    {
        inner_.visit(make_path(path), value);
    }
    auto visit(std::string_view path, VectorXiRecord value) -> void override
    {
        inner_.visit(make_path(path), value);
    }
    auto visit(std::string_view path, DoubleRecord value) -> void override
    {
        inner_.visit(make_path(path), value);
    }

   private:
    auto make_path(std::string_view path) const -> std::string
    {
        if (path.empty())
        {
            return prefix_;
        }
        return fmt::format("{}/{}", prefix_, path);
    }

    std::string prefix_;
    Base& inner_;
};

template <bool Mutable, typename Fn>
auto visit_with_prefix(
    infra::BasicRecordSink<Mutable>& sink,
    std::string prefix,
    Fn&& fn) -> void
{
    PrefixedRecordSink<Mutable> prefixed(std::move(prefix), sink);
    std::forward<Fn>(fn)(prefixed);
}

}  // namespace

FixedState::FixedState(const FixedEffect& effect)
    : coeffs(Eigen::VectorXd::Zero(effect.X.cols()))
{
}

FixedState::FixedState(Eigen::VectorXd coeffs) : coeffs(std::move(coeffs)) {}

auto FixedState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("coeffs", coeffs, FieldFlag::checkpoint | FieldFlag::trace);
}

auto FixedState::visit_records(StateRecordSet, infra::RecordSink& sink) const
    -> void
{
    sink.visit("coeffs", coeffs);
}

auto FixedState::visit_records(
    StateRecordSet set,
    infra::MutableRecordSink& sink) -> void
{
    if (set != StateRecordSet::checkpoint)
    {
        throw GelexException(
            "FixedState: mutable visit_records requires checkpoint set");
    }
    sink.visit("coeffs", coeffs);
}

RandomState::RandomState(const RandomEffect& effect, double variance)
    : coeffs(Eigen::VectorXd::Zero(effect.X.cols())), variance{variance}
{
}

RandomState::RandomState(const RandomEffect& effect, const RandomPrior& prior)
    : RandomState(effect, prior.initial_value())
{
}

RandomState::RandomState(Eigen::VectorXd coeffs, double variance)
    : coeffs(std::move(coeffs)), variance{variance}
{
}

auto RandomState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("coeffs", coeffs, FieldFlag::checkpoint | FieldFlag::trace);
    visitor.on("variance", variance, FieldFlag::checkpoint | FieldFlag::trace);
}

auto RandomState::visit_records(StateRecordSet, infra::RecordSink& sink) const
    -> void
{
    sink.visit("coeffs", coeffs);
    sink.visit("variance", variance);
}

auto RandomState::visit_records(
    StateRecordSet set,
    infra::MutableRecordSink& sink) -> void
{
    if (set != StateRecordSet::checkpoint)
    {
        throw GelexException(
            "RandomState: mutable visit_records requires checkpoint set");
    }
    sink.visit("coeffs", coeffs);
    sink.visit("variance", variance);
}

GeneticState::GeneticState(
    GeneticMode type,
    Eigen::Index num_markers,
    Eigen::Index num_individuals)
    : type(type),
      coeffs(Eigen::VectorXd::Zero(num_markers)),
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

auto GeneticState::visit_records(StateRecordSet set, infra::RecordSink& sink)
    const -> void
{
    switch (set)
    {
        case StateRecordSet::sample:
            sink.visit("coeffs", coeffs);
            sink.visit("variance", variance);
            sink.visit("heritability", heritability);
            break;
        case StateRecordSet::checkpoint:
            sink.visit("coeffs", coeffs);
            sink.visit("u", u);
            sink.visit("variance", variance);
            sink.visit("heritability", heritability);
            break;
    }
}

auto GeneticState::visit_records(
    StateRecordSet set,
    infra::MutableRecordSink& sink) -> void
{
    if (set != StateRecordSet::checkpoint)
    {
        throw GelexException(
            "GeneticState: mutable visit_records requires checkpoint set");
    }
    sink.visit("coeffs", coeffs);
    sink.visit("u", u);
    sink.visit("variance", variance);
    sink.visit("heritability", heritability);
}

GeneticBlockState::GeneticBlockState(
    std::vector<GeneticMode> modes,
    std::vector<std::size_t> genetic_indices,
    std::unique_ptr<GeneticPriorState> prior_state)
    : modes_(std::move(modes)),
      genetic_indices_(std::move(genetic_indices)),
      prior_state_(std::move(prior_state))
{
    if (modes_.size() != genetic_indices_.size())
    {
        throw GelexException(
            "GeneticBlockState: modes and genetic indices differ in size");
    }
    if (prior_state_ == nullptr)
    {
        throw GelexException("GeneticBlockState: null prior state");
    }
}

auto GeneticBlockState::slot(GeneticMode mode) const -> std::size_t
{
    for (std::size_t i = 0; i < modes_.size(); ++i)
    {
        if (modes_[i] == mode)
        {
            return i;
        }
    }
    throw GelexException(
        fmt::format("GeneticBlockState: missing genetic mode {}", mode));
}

auto GeneticBlockState::visit_records(
    StateRecordSet set,
    infra::RecordSink& sink) const -> void
{
    prior_state_->visit_records(set, sink);
}

auto GeneticBlockState::visit_records(
    StateRecordSet set,
    infra::MutableRecordSink& sink) -> void
{
    if (set != StateRecordSet::checkpoint)
    {
        throw GelexException(
            "GeneticBlockState: mutable visit_records requires checkpoint set");
    }
    prior_state_->visit_records(set, sink);
}

SingleGeneticBlockState::SingleGeneticBlockState(
    const GeneticEffect& effect,
    const SingleGeneticPrior& prior)
    : SingleGeneticBlockState(
          GeneticState{effect.type, effect.X.cols(), effect.X.rows()},
          prior.make_state(effect.X.cols(), effect.X.rows()))
{
    if (effect.type != prior.mode())
    {
        throw GelexException(
            fmt::format(
                "SingleGeneticBlockState: effect mode {} != prior mode {}",
                effect.type,
                prior.mode()));
    }
}

SingleGeneticBlockState::SingleGeneticBlockState(
    GeneticState genetic,
    std::unique_ptr<SingleGeneticPriorState> prior_state)
    : state_(std::move(genetic)), prior_state_(std::move(prior_state))
{
    if (prior_state_ == nullptr)
    {
        throw GelexException("SingleGeneticBlockState: null prior state");
    }
}

auto SingleGeneticBlockState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    state_.visit(visitor);
    auto prior_scope = visitor.scope("prior_state");
    prior_state_->visit(visitor);
}

JointGeneticBlockState::JointGeneticBlockState(
    const GeneticEffect& additive,
    const GeneticEffect& dominance,
    const JointGeneticPrior& prior)
    : JointGeneticBlockState(
          GeneticState{GeneticMode::A, additive.X.cols(), additive.X.rows()},
          GeneticState{GeneticMode::D, dominance.X.cols(), dominance.X.rows()},
          prior.make_state(additive.X.cols(), additive.X.rows()))
{
    if (additive.type != GeneticMode::A)
    {
        throw GelexException(
            fmt::format(
                "JointGeneticBlockState: additive effect has mode {}",
                additive.type));
    }
    if (dominance.type != GeneticMode::D)
    {
        throw GelexException(
            fmt::format(
                "JointGeneticBlockState: dominance effect has mode {}",
                dominance.type));
    }
    if (additive.X.rows() != dominance.X.rows()
        || additive.X.cols() != dominance.X.cols())
    {
        throw GelexException(
            fmt::format(
                "JointGeneticBlockState: genetic effects must share shape; "
                "A is {}x{}, D is {}x{}",
                additive.X.rows(),
                additive.X.cols(),
                dominance.X.rows(),
                dominance.X.cols()));
    }
}

JointGeneticBlockState::JointGeneticBlockState(
    GeneticState additive,
    GeneticState dominance,
    std::unique_ptr<JointGeneticPriorState> prior_state)
    : additive_(std::move(additive)),
      dominance_(std::move(dominance)),
      prior_state_(std::move(prior_state))
{
    if (additive_.type != GeneticMode::A)
    {
        throw GelexException(
            fmt::format(
                "JointGeneticBlockState: additive state has mode {}",
                additive_.type));
    }
    if (dominance_.type != GeneticMode::D)
    {
        throw GelexException(
            fmt::format(
                "JointGeneticBlockState: dominance state has mode {}",
                dominance_.type));
    }
    if (prior_state_ == nullptr)
    {
        throw GelexException("JointGeneticBlockState: null prior state");
    }
}

auto JointGeneticBlockState::contains(GeneticMode mode) const -> bool
{
    return mode == GeneticMode::A || mode == GeneticMode::D;
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
    prior_state_->visit(visitor);
}

auto ResidualState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("y_adj", y_adj, FieldFlag::checkpoint);
    visitor.on("variance", variance, FieldFlag::checkpoint | FieldFlag::trace);
}

auto ResidualState::visit_records(StateRecordSet set, infra::RecordSink& sink)
    const -> void
{
    switch (set)
    {
        case StateRecordSet::sample:
            sink.visit("variance", variance);
            break;
        case StateRecordSet::checkpoint:
            sink.visit("y_adj", y_adj);
            sink.visit("variance", variance);
            break;
    }
}

auto ResidualState::visit_records(
    StateRecordSet set,
    infra::MutableRecordSink& sink) -> void
{
    if (set != StateRecordSet::checkpoint)
    {
        throw GelexException(
            "ResidualState: mutable visit_records requires checkpoint set");
    }
    sink.visit("y_adj", y_adj);
    sink.visit("variance", variance);
}

}  // namespace gelex::bayes

namespace gelex
{

namespace
{

auto require_model_effect(const BayesModel& model, GeneticMode mode)
    -> const bayes::GeneticEffect&
{
    const auto* effect = model.genetic(mode);
    if (effect == nullptr)
    {
        throw GelexException(
            fmt::format(
                "BayesState: missing genetic effect for mode {}", mode));
    }
    return *effect;
}

auto require_shared_shape(
    const bayes::GeneticEffect& first,
    const bayes::GeneticEffect& current) -> void
{
    if (first.X.rows() != current.X.rows()
        || first.X.cols() != current.X.cols())
    {
        throw GelexException(
            fmt::format(
                "BayesState: genetic prior block modes must share shape; "
                "{} is {}x{}, {} is {}x{}",
                first.type,
                first.X.rows(),
                first.X.cols(),
                current.type,
                current.X.rows(),
                current.X.cols()));
    }
}

}  // namespace

BayesState::BayesState(const BayesModel& model, const bayes::BayesPrior& prior)
    : fixed_(model.fixed()),
      residual_{
          .y_adj = model.phenotype(),
          .variance = prior.residual().initial_value()}
{
    const auto& random_effects = model.random();
    random_.reserve(random_effects.size());
    for (const auto& effect : random_effects)
    {
        random_.emplace_back(effect, prior.random());
    }

    for (const auto& block : prior.genetics())
    {
        const auto modes = block.modes();
        const auto& first_effect = require_model_effect(model, modes.front());
        auto prior_state
            = block.make_state(first_effect.X.cols(), model.num_individuals());

        std::vector<GeneticMode> block_modes;
        std::vector<std::size_t> genetic_indices;
        block_modes.reserve(modes.size());
        genetic_indices.reserve(modes.size());

        for (const auto mode : modes)
        {
            const auto& effect = require_model_effect(model, mode);
            require_shared_shape(first_effect, effect);
            block_modes.push_back(mode);
            genetic_indices.push_back(genetics_.size());
            genetics_.emplace_back(mode, effect.X.cols(), effect.X.rows());
        }

        genetic_blocks_.emplace_back(
            std::move(block_modes),
            std::move(genetic_indices),
            std::move(prior_state));
    }
}

auto BayesState::compute_heritability() -> void
{
    double total_variance = residual_.variance;
    for (const auto& state : random_)
    {
        total_variance += state.variance;
    }
    for (const auto& state : genetics_)
    {
        total_variance += state.variance;
    }

    for (auto& state : genetics_)
    {
        state.heritability = state.variance / total_variance;
    }
}

auto BayesState::visit_records(
    bayes::StateRecordSet set,
    infra::RecordSink& sink) const -> void
{
    bayes::visit_with_prefix<false>(
        sink,
        "fixed/0",
        [set, this](infra::RecordSink& prefixed)
        { fixed_.visit_records(set, prefixed); });

    for (std::size_t i = 0; i < random_.size(); ++i)
    {
        bayes::visit_with_prefix<false>(
            sink,
            fmt::format("random/{}", i),
            [set, this, i](infra::RecordSink& prefixed)
            { random_[i].visit_records(set, prefixed); });
    }

    for (std::size_t i = 0; i < genetics_.size(); ++i)
    {
        bayes::visit_with_prefix<false>(
            sink,
            fmt::format("genetic/{}", i),
            [set, this, i](infra::RecordSink& prefixed)
            { genetics_[i].visit_records(set, prefixed); });
    }

    for (std::size_t i = 0; i < genetic_blocks_.size(); ++i)
    {
        bayes::visit_with_prefix<false>(
            sink,
            fmt::format("genetic_block/{}/prior_state", i),
            [set, this, i](infra::RecordSink& prefixed)
            { genetic_blocks_[i].visit_records(set, prefixed); });
    }

    bayes::visit_with_prefix<false>(
        sink,
        "residual/0",
        [set, this](infra::RecordSink& prefixed)
        { residual_.visit_records(set, prefixed); });
}

auto BayesState::visit_records(
    bayes::StateRecordSet set,
    infra::MutableRecordSink& sink) -> void
{
    if (set != bayes::StateRecordSet::checkpoint)
    {
        throw GelexException(
            "BayesState: mutable visit_records requires checkpoint set");
    }

    bayes::visit_with_prefix<true>(
        sink,
        "fixed/0",
        [this](infra::MutableRecordSink& prefixed)
        { fixed_.visit_records(bayes::StateRecordSet::checkpoint, prefixed); });

    for (std::size_t i = 0; i < random_.size(); ++i)
    {
        bayes::visit_with_prefix<true>(
            sink,
            fmt::format("random/{}", i),
            [this, i](infra::MutableRecordSink& prefixed)
            {
                random_[i].visit_records(
                    bayes::StateRecordSet::checkpoint, prefixed);
            });
    }

    for (std::size_t i = 0; i < genetics_.size(); ++i)
    {
        bayes::visit_with_prefix<true>(
            sink,
            fmt::format("genetic/{}", i),
            [this, i](infra::MutableRecordSink& prefixed)
            {
                genetics_[i].visit_records(
                    bayes::StateRecordSet::checkpoint, prefixed);
            });
    }

    for (std::size_t i = 0; i < genetic_blocks_.size(); ++i)
    {
        bayes::visit_with_prefix<true>(
            sink,
            fmt::format("genetic_block/{}/prior_state", i),
            [this, i](infra::MutableRecordSink& prefixed)
            {
                genetic_blocks_[i].visit_records(
                    bayes::StateRecordSet::checkpoint, prefixed);
            });
    }

    bayes::visit_with_prefix<true>(
        sink,
        "residual/0",
        [this](infra::MutableRecordSink& prefixed)
        {
            residual_.visit_records(
                bayes::StateRecordSet::checkpoint, prefixed);
        });
}

BayesStateV2::BayesStateV2(
    const BayesModel& model,
    const bayes::BayesPriorV2& prior)
    : fixed_(model.fixed()),
      residual_{
          .y_adj = model.phenotype(),
          .variance = prior.residual().initial_value()}
{
    const auto& random_effects = model.random();
    random_.reserve(random_effects.size());
    for (const auto& effect : random_effects)
    {
        random_.emplace_back(effect, prior.random());
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
                        std::unique_ptr<bayes::SingleGeneticPrior>>)
                {
                    const auto mode = genetic_prior->mode();
                    const auto* effect = model.genetic(mode);
                    if (effect == nullptr)
                    {
                        throw GelexException(
                            fmt::format(
                                "BayesStateV2: missing genetic effect for "
                                "mode {}",
                                mode));
                    }
                    genetics_.emplace_back(
                        bayes::SingleGeneticBlockState{
                            *effect, *genetic_prior});
                }
                else
                {
                    const auto* additive = model.genetic(GeneticMode::A);
                    if (additive == nullptr)
                    {
                        throw GelexException(
                            "BayesStateV2: missing genetic effect for mode A");
                    }

                    const auto* dominance = model.genetic(GeneticMode::D);
                    if (dominance == nullptr)
                    {
                        throw GelexException(
                            "BayesStateV2: missing genetic effect for mode D");
                    }

                    genetics_.emplace_back(
                        bayes::JointGeneticBlockState{
                            *additive, *dominance, *genetic_prior});
                }
            },
            block);
    }
}

auto BayesStateV2::genetic(GeneticMode type) -> bayes::GeneticState*
{
    for (auto& block : genetics_)
    {
        auto* state = std::visit(
            [type](auto& value) -> bayes::GeneticState*
            {
                if constexpr (
                    std::is_same_v<
                        std::decay_t<decltype(value)>,
                        bayes::SingleGeneticBlockState>)
                {
                    return value.contains(type) ? &value.state() : nullptr;
                }
                else
                {
                    return value.contains(type) ? &value.state(type) : nullptr;
                }
            },
            block);
        if (state != nullptr)
        {
            return state;
        }
    }
    return nullptr;
}

auto BayesStateV2::genetic(GeneticMode type) const -> const bayes::GeneticState*
{
    return const_cast<BayesStateV2*>(this)->genetic(type);
}

auto BayesStateV2::genetic_block_for(GeneticMode type)
    -> bayes::GeneticPriorBlockState*
{
    for (auto& block : genetics_)
    {
        if (std::visit(
                [type](const auto& value) { return value.contains(type); },
                block))
        {
            return &block;
        }
    }
    return nullptr;
}

auto BayesStateV2::genetic_block_for(GeneticMode type) const
    -> const bayes::GeneticPriorBlockState*
{
    return const_cast<BayesStateV2*>(this)->genetic_block_for(type);
}

auto BayesStateV2::compute_heritability() -> void
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

auto BayesStateV2::visit(infra::FieldVisitor& visitor) -> void
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
