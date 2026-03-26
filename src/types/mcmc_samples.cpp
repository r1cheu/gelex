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

#include "gelex/types/mcmc_samples.h"

#include <memory>
#include <ranges>
#include <string_view>
#include <variant>

#include <Eigen/Core>

#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/states.h"
#include "gelex/model/bayes/writer/mcmc_writer.h"
#include "gelex/types/fixed_effects.h"

namespace gelex
{
using Eigen::Index;

MCMCSamples::~MCMCSamples() = default;
MCMCSamples::MCMCSamples(MCMCSamples&&) noexcept = default;
auto MCMCSamples::operator=(MCMCSamples&&) noexcept -> MCMCSamples& = default;

FixedSamples::FixedSamples(const FixedEffect& effect)
    : names(effect.names), levels(effect.levels), n_coeffs_(effect.X.cols())
{
}

RandomSamples::RandomSamples(const bayes::RandomEffect& effect)
    : levels(effect.levels), n_coeffs_(effect.X.cols())
{
}

void FixedSamples::store(const bayes::FixedState& state)
{
    coeffs_stats_.update(state.coeffs);
}

void RandomSamples::store(const bayes::RandomState& state)
{
    coeffs_stats_.update(state.coeffs);
    variance_stats_.update(state.variance);
}

void ResidualSamples::store(const bayes::ResidualState& state)
{
    variance_stats_.update(state.variance);
}

AssignmentSamples::AssignmentSamples(
    Eigen::Index n_snps,
    Eigen::Index n_proportions,
    bool estimate_pi)
    : estimate_pi_(estimate_pi),
      comp_counts_(Eigen::MatrixXd::Zero(n_snps, n_proportions))
{
}

void AssignmentSamples::store(const bayes::Assignment& alloc)
{
    if (estimate_pi_)
    {
        proportion_stats_.update(alloc.proportion);
    }
    ++n_samples_;
    for (Index i = 0; i < alloc.tracker.size(); ++i)
    {
        comp_counts_(i, alloc.tracker(i)) += 1.0;
    }
}

ComponentSamples::ComponentSamples(
    Eigen::Index n_snps,
    Eigen::Index n_proportions,
    bool estimate_pi)
    : assignment(n_snps, n_proportions, estimate_pi)
{
}

void ComponentSamples::store(const bayes::ComponentAllocation& alloc)
{
    assignment.store(alloc.assignment);
    comp_var_stats_.update(alloc.component_variance);
}

auto GeneticSamples::make_group_samples(
    const bayes::GeneticEffect& effect,
    const bayes::GeneticPrior& prior) -> std::optional<MixtureSamples>
{
    const auto n_snps = bayes::get_cols(effect.X);
    if (const auto* sp = std::get_if<bayes::SpikePrior>(&prior.marker))
    {
        return AssignmentSamples{
            n_snps, sp->proportion.init.size(), sp->proportion.estimate};
    }
    if (const auto* mp = std::get_if<bayes::MixturePrior>(&prior.marker))
    {
        return ComponentSamples{
            n_snps, mp->proportion.init.size(), mp->proportion.estimate};
    }
    return std::nullopt;
}

GeneticSamples::GeneticSamples(
    const bayes::GeneticEffect& effect,
    const bayes::GeneticPrior& prior)
    : type(effect.type),
      group(make_group_samples(effect, prior)),
      n_coeffs_(bayes::get_cols(effect.X))
{
    if (prior.sign)
    {
        const auto n_snps = bayes::get_cols(effect.X);
        sign.emplace(n_snps, 3, true);
    }
}

void GeneticSamples::store(const bayes::GeneticState& state)
{
    coeffs_stats_.update(state.coeffs);
    variance_stats_.update(state.variance);
    heritability_stats_.update(state.heritability);

    if (group && state.group)
    {
        std::visit(
            [&](auto& g, const auto& s)
            {
                using G = std::decay_t<decltype(g)>;
                using S = std::decay_t<decltype(s)>;
                if constexpr (
                    (std::is_same_v<G, AssignmentSamples>
                     && std::is_same_v<S, bayes::Assignment>)
                    || (std::is_same_v<G, ComponentSamples>
                        && std::is_same_v<S, bayes::ComponentAllocation>))
                {
                    g.store(s);
                }
            },
            *group,
            *state.group);
    }
    if (sign && state.sign)
    {
        sign->store(*state.sign);
    }
}

MCMCSamples::MCMCSamples(
    const BayesModel& model,
    const bayes::Priors& priors,
    std::string_view sample_prefix,
    Eigen::Index n_records)
    : fixed_(model.fixed())
{
    if (const auto& effects = model.random(); !effects.empty())
    {
        random_.reserve(effects.size());
        for (const auto& effect : effects)
        {
            random_.emplace_back(effect);
        }
    }

    for (const auto& effect : model.genetics())
    {
        const auto* prior = priors.genetic(effect.type);
        genetics_.emplace_back(effect, *prior);
    }

    if (!sample_prefix.empty())
    {
        writer_ = std::make_unique<MCMCWriter>(
            model, priors, sample_prefix, n_records);
    }
}

void MCMCSamples::store(const BayesState& states)
{
    fixed_.store(states.fixed());

    for (auto&& [sample, state] : std::views::zip(random_, states.random()))
    {
        sample.store(state);
    }

    for (auto&& [sample, state] : std::views::zip(genetics_, states.genetics()))
    {
        sample.store(state);
    }

    residual_.store(states.residual());

    if (writer_)
    {
        writer_->write(states);
    }
}

void MCMCSamples::finalize()
{
    if (writer_)
    {
        writer_->finalize();
    }
}

}  // namespace gelex
