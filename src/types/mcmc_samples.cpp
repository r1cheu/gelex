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

#include <Eigen/Core>

#include "../src/types/bayes_effects.h"
#include "gelex/estimator/bayes/params.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{
using Eigen::Index;

FixedSamples::FixedSamples(const MCMCParams& params, const FixedEffect& effect)
{
    coeffs.emplace_back(effect.X.cols(), params.n_records);
}
RandomSamples::RandomSamples(
    const MCMCParams& params,
    const bayes::RandomEffect& effect)
    : RandomSamples(params, effect.X.cols()) {};

RandomSamples::RandomSamples(const MCMCParams& params, Eigen::Index n_coeffs)
{
    coeffs.emplace_back(n_coeffs, params.n_records);
    variance.emplace_back(1, params.n_records);
}

BaseMarkerSamples::BaseMarkerSamples(
    const MCMCParams& params,
    const bayes::GeneticEffect& effect)
    : RandomSamples(params, bayes::get_cols(effect.X))
{
    heritability.emplace_back(1, params.n_records);

    if (effect.init_pi)  // mixture model
    {
        const Eigen::Index num_snp = bayes::get_cols(effect.X);
        n_proportions = effect.init_pi->size();
        tracker.emplace_back(num_snp, params.n_records);

        if (n_proportions > 2)
        {
            component_variance.emplace_back(
                n_proportions - 1, params.n_records);
        }
    }
    if (effect.estimate_pi)
    {
        const Eigen::Index n_props = effect.init_pi->size();
        mixture_proportion.emplace_back(n_props, params.n_records);
    }
}

AdditiveSamples::AdditiveSamples(
    const MCMCParams& params,
    const bayes::AdditiveEffect& effect)
    : BaseMarkerSamples(params, effect)
{
}

DominantSamples::DominantSamples(
    const MCMCParams& params,
    const bayes::DominantEffect& effect)
    : BaseMarkerSamples(params, effect)
{
}

ResidualSamples::ResidualSamples(const MCMCParams& params)
{
    variance.emplace_back(1, params.n_records);
}

MCMCSamples::MCMCSamples(
    const MCMCParams& params,
    const BayesModel& model,
    std::string_view sample_prefix)
    : residual_(params)
{
    if (const auto* effect = model.fixed(); effect)
    {
        fixed_.emplace(params, *effect);
    }

    if (const auto& effects = model.random(); !effects.empty())
    {
        random_.reserve(effects.size());
        for (const auto& effect : effects)
        {
            random_.emplace_back(params, effect);
        }
    }

    if (const auto* effect = model.additive(); effect)
    {
        additive_.emplace(params, *effect);
        add_writer_
            = sample_prefix.empty()
                  ? nullptr
                  : std::make_unique<detail::BinaryMatrixWriter>(
                        std::filesystem::path(
                            std::format("{}.add.sample", sample_prefix)));
    }

    if (const auto* effect = model.dominant(); effect)
    {
        dominant_.emplace(params, *effect);
        dom_writer_
            = sample_prefix.empty()
                  ? nullptr
                  : std::make_unique<detail::BinaryMatrixWriter>(
                        std::filesystem::path(
                            std::format("{}.dom.sample", sample_prefix)));
    }
}

void MCMCSamples::store(const BayesState& states, Eigen::Index record_idx)
{
    constexpr Eigen::Index chain_idx = 0;
    if (const auto* state = states.fixed(); fixed_ && state != nullptr)
    {
        fixed_->coeffs[chain_idx].col(record_idx) = state->coeffs;
    }

    for (auto&& [sample, state] : std::views::zip(random_, states.random()))
    {
        sample.coeffs[chain_idx].col(record_idx) = state.coeffs;
        sample.variance[chain_idx](0, record_idx) = state.variance;
    }

    if (const auto* state = states.additive(); additive_ && state != nullptr)
    {
        additive_->coeffs[chain_idx].col(record_idx) = state->coeffs;
        if (add_writer_)
        {
            add_writer_->write(state->coeffs);
        }
        additive_->variance[chain_idx](0, record_idx) = state->variance;
        additive_->heritability[chain_idx](0, record_idx) = state->heritability;

        if (!additive_->mixture_proportion.empty())
        {
            additive_->mixture_proportion[chain_idx].col(record_idx)
                = state->pi.prop;
        }
        if (!additive_->tracker.empty())
        {
            additive_->tracker[chain_idx].col(record_idx) = state->tracker;
        }

        if (!additive_->component_variance.empty())
        {
            additive_->component_variance[chain_idx].col(record_idx)
                = state->component_variance;
        }
    }

    if (const auto* state = states.dominant(); dominant_ && state != nullptr)
    {
        dominant_->coeffs[chain_idx].col(record_idx) = state->coeffs;
        if (dom_writer_)
        {
            dom_writer_->write(state->coeffs);
        }
        dominant_->variance[chain_idx](0, record_idx) = state->variance;
        dominant_->heritability[chain_idx](0, record_idx) = state->heritability;

        if (!dominant_->mixture_proportion.empty()
            && state->pi.prop.size() != 0)
        {
            dominant_->mixture_proportion[chain_idx].col(record_idx)
                = state->pi.prop;
        }
        if (!dominant_->tracker.empty() && state->tracker.size() != 0)
        {
            dominant_->tracker[chain_idx].col(record_idx) = state->tracker;
        }

        if (!dominant_->component_variance.empty())
        {
            dominant_->component_variance[chain_idx].col(record_idx)
                = state->component_variance;
        }
    }

    residual_.variance[chain_idx](0, record_idx) = states.residual().variance;
}

}  // namespace gelex
