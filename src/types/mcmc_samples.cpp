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

#include <Eigen/Core>

#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/model.h"
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

MixtureSamples::MixtureSamples(const bayes::GeneticEffect& effect)
    : n_snps_(bayes::get_cols(effect.X)),
      n_proportions_(effect.mixture->init_proportion.size()),
      estimate_pi_(effect.mixture->estimate_pi)
{
    comp_counts_ = Eigen::MatrixXd::Zero(n_snps_, n_proportions_);
}

void MixtureSamples::store(const bayes::MixtureState& state)
{
    if (estimate_pi_)
    {
        proportion_stats_.update(state.pi.proportion);
    }
    if (n_proportions_ > 2)
    {
        comp_var_stats_.update(state.component_variance);
    }
    ++n_samples_;
    for (Index i = 0; i < state.tracker.size(); ++i)
    {
        comp_counts_(i, state.tracker(i)) += 1.0;
    }
}

void SignSamples::store(const bayes::SignState& state)
{
    positive_prob_stats_.update(state.positive_prob);
}

GeneticSamples::GeneticSamples(const bayes::GeneticEffect& effect)
    : type(effect.type), n_coeffs_(bayes::get_cols(effect.X))
{
    if (effect.mixture)
    {
        mixture.emplace(effect);
    }
    if (effect.sign)
    {
        sign.emplace();
    }
}

void GeneticSamples::store(const bayes::GeneticState& state)
{
    coeffs_stats_.update(state.coeffs);
    variance_stats_.update(state.variance);
    heritability_stats_.update(state.heritability);

    if (mixture && state.mixture)
    {
        mixture->store(*state.mixture);
    }
    if (sign && state.sign)
    {
        sign->store(*state.sign);
    }
}

MCMCSamples::MCMCSamples(
    const BayesModel& model,
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
        genetics_.emplace_back(effect);
    }

    if (!sample_prefix.empty())
    {
        writer_ = std::make_unique<MCMCWriter>(model, sample_prefix, n_records);
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
