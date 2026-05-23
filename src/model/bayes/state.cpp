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
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

FixedState::FixedState(const FixedEffect& effect)
    : coeffs(Eigen::VectorXd::Zero(effect.X.cols()))
{
}

FixedState::FixedState(Eigen::VectorXd coeffs) : coeffs(std::move(coeffs)) {}

RandomState::RandomState(const RandomEffect& effect, double variance)
    : coeffs(Eigen::VectorXd::Zero(effect.X.cols())), variance{variance}
{
}

RandomState::RandomState(const RandomEffect& effect, const VarianceSpec& spec)
    : RandomState(effect, spec.initial_value())
{
}

RandomState::RandomState(Eigen::VectorXd coeffs, double variance)
    : coeffs(std::move(coeffs)), variance{variance}
{
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

}  // namespace gelex
