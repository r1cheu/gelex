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

#include "gelex/types/mcmc_results.h"

#include <optional>
#include <ranges>

#include <Eigen/Core>

#include "../src/estimator/bayes/posterior_calculator.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{

using Eigen::Index;
using Eigen::VectorXd;

MCMCResult::MCMCResult(
    MCMCSamples&& samples,
    const BayesModel& model,
    double prob)
    : samples_(std::move(samples)),
      residual_(1),
      prob_(prob),
      phenotype_var_(model.phenotype_variance())
{
    if (const auto* effect = model.additive(); effect)
    {
        p_freq = bayes::get_means(effect->X).array() / 2;
    }

    if (const auto* sample = samples_.fixed(); sample)
    {
        fixed_.emplace(*sample);
    }

    for (const auto& sample : samples_.random())
    {
        random_.emplace_back(sample);
    }

    if (const auto* sample = samples_.additive(); sample)
    {
        additive_.emplace(*sample);
    }

    if (const auto* sample = samples_.dominant(); sample)
    {
        dominant_.emplace(*sample);
    }
}

void MCMCResult::compute(std::optional<double> prob)
{
    if (prob)
    {
        prob_ = prob.value();
    }

    // Use PosteriorCalculator for all computations
    if (const auto* sample = samples_.fixed(); fixed_ && sample != nullptr)
    {
        fixed_->coeffs = detail::PosteriorCalculator::compute_param_summary(
            samples_.fixed()->coeffs, prob_);
    }

    for (auto&& [result, sample] : std::views::zip(random_, samples_.random()))
    {
        result.coeffs = detail::PosteriorCalculator::compute_param_summary(
            sample.coeffs, prob_);
        result.variance = detail::PosteriorCalculator::compute_param_summary(
            sample.variance, prob_);
    }

    auto compute_summary = [&](auto& effect, const auto* sample)
    {
        effect->coeffs = detail::PosteriorCalculator::compute_param_summary(
            sample->coeffs, prob_);
        effect->variance = detail::PosteriorCalculator::compute_param_summary(
            sample->variance, prob_);
        effect->heritability
            = detail::PosteriorCalculator::compute_param_summary(
                sample->heritability, prob_);

        if (effect->mixture_proportion.size() > 0)
        {
            effect->mixture_proportion
                = detail::PosteriorCalculator::compute_param_summary(
                    sample->mixture_proportion, prob_);
        }
        if (effect->component_variance.size() > 0)
        {
            effect->component_variance
                = detail::PosteriorCalculator::compute_param_summary(
                    sample->component_variance, prob_);
        }
        if (effect->pip.size() > 0)
        {
            const auto n_comp = effect->comp_probs.cols();
            effect->comp_probs
                = detail::PosteriorCalculator::compute_component_probs(
                    sample->tracker, n_comp);
            effect->pip
                = effect->comp_probs.rightCols(n_comp - 1).rowwise().sum();
        }

        detail::PosteriorCalculator::compute_pve(
            effect->pve, sample->coeffs, phenotype_var_);
    };

    if (const auto* sample = samples_.additive();
        additive_ && sample != nullptr)
    {
        compute_summary(additive_, sample);
    }

    if (const auto* sample = samples_.dominant();
        dominant_ && sample != nullptr)
    {
        compute_summary(dominant_, sample);
    }

    // pi and snp_tracker functionality is now handled within AdditiveSummary
    // and DominantSummary

    residual_ = detail::PosteriorCalculator::compute_param_summary(
        samples_.residual().variance, prob_);
}
}  // namespace gelex
