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

#include "gelex/algo/infer/mcmc/samples.h"

#include <memory>
#include <ranges>
#include <string_view>
#include <variant>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/exception.h"
#include "gelex/io/mcmc/sample_writer.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/state_capabilities.h"

namespace gelex::mcmc
{
using Eigen::Index;

Samples::~Samples() = default;
Samples::Samples(Samples&&) noexcept = default;
auto Samples::operator=(Samples&&) noexcept -> Samples& = default;

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

void AssignmentSamples::store(const bayes::ProportionState& state)
{
    if (estimate_pi_)
    {
        proportion_stats_.update(state.value);
    }
    ++n_samples_;
    for (Index i = 0; i < state.assignment.size(); ++i)
    {
        comp_counts_(i, state.assignment(i)) += 1.0;
    }
}

ComponentSamples::ComponentSamples(
    Eigen::Index n_snps,
    Eigen::Index n_proportions,
    bool estimate_pi)
    : assignment(n_snps, n_proportions, estimate_pi)
{
}

void ComponentSamples::store(
    const bayes::ComponentState& component,
    const bayes::ProportionState& proportion)
{
    assignment.store(proportion);
    comp_var_stats_.update(component.gebv_var);
}

auto GeneticSamples::make_group_samples(
    const bayes::GeneticEffect& effect,
    const bayes::GeneticPrior& prior,
    const bayes::GeneticBlockState& block,
    std::size_t slot) -> std::optional<MixtureSamples>
{
    const auto* proportion_cap = prior.query<bayes::ProportionSpecCap>();
    if (proportion_cap == nullptr)
    {
        return std::nullopt;
    }
    const auto n_snps = effect.X.cols();
    const auto specs = proportion_cap->proportion();
    const auto n_pi = specs[slot].size();
    const auto estimate_pi = specs[slot].sampled();
    const bool has_components
        = block.prior_state().query<bayes::ComponentStateCap>() != nullptr;
    if (has_components)
    {
        return ComponentSamples{n_snps, n_pi, estimate_pi};
    }
    return AssignmentSamples{n_snps, n_pi, estimate_pi};
}

GeneticSamples::GeneticSamples(
    const bayes::GeneticEffect& effect,
    const bayes::GeneticPrior& prior,
    const bayes::GeneticBlockState& block,
    std::size_t block_index,
    GeneticMode mode)
    : type(mode),
      group(make_group_samples(effect, prior, block, block.slot(mode))),
      block_index_(block_index),
      slot_(block.slot(mode)),
      n_coeffs_(effect.X.cols())
{
}

void GeneticSamples::store(const BayesState& state)
{
    const auto* genetic = state.genetic(type);
    if (genetic == nullptr)
    {
        throw GelexException("GeneticSamples: missing genetic state");
    }
    coeffs_stats_.update(genetic->coeffs);
    variance_stats_.update(genetic->variance);
    heritability_stats_.update(genetic->heritability);

    if (group)
    {
        auto& prior_state = state.genetic_block(block_index_).prior_state();
        auto& proportion = prior_state.require<bayes::ProportionStateCap>()
                               .proportion()[slot_];
        std::visit(
            [&](auto& g)
            {
                using G = std::decay_t<decltype(g)>;
                if constexpr (std::is_same_v<G, AssignmentSamples>)
                {
                    g.store(proportion);
                }
                else
                {
                    const auto& component
                        = prior_state.require<bayes::ComponentStateCap>()
                              .component()[slot_];
                    g.store(component, proportion);
                }
            },
            *group);
    }
}

Samples::Samples(
    const BayesModel& model,
    const bayes::BayesPrior& prior,
    const BayesState& state,
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

    std::vector<const bayes::GeneticPrior*> prior_blocks;
    for (const auto& block : prior.genetics())
    {
        prior_blocks.push_back(&block);
    }

    for (auto&& [block_index, block] :
         std::views::enumerate(state.genetic_blocks()))
    {
        const auto idx = static_cast<std::size_t>(block_index);
        for (const auto mode : block.modes())
        {
            const auto* effect = model.genetic(mode);
            if (effect == nullptr)
            {
                throw GelexException("Samples: missing genetic effect");
            }
            genetics_.emplace_back(
                *effect, *prior_blocks.at(idx), block, idx, mode);
        }
    }

    if (!sample_prefix.empty())
    {
        writer_
            = std::make_unique<mcmc::Writer>(state, sample_prefix, n_records);
    }
}

void Samples::store(const mcmc::State& states)
{
    fixed_.store(states.fixed());

    for (auto&& [sample, state] : std::views::zip(random_, states.random()))
    {
        sample.store(state);
    }

    for (auto& sample : genetics_)
    {
        sample.store(states);
    }

    residual_.store(states.residual());

    if (writer_)
    {
        writer_->write(states);
    }
}

void Samples::finalize()
{
    writer_.reset();
}

}  // namespace gelex::mcmc
