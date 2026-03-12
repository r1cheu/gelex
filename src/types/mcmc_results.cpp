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

#include <ranges>

#include <Eigen/Core>

#include "gelex/algo/infer/posterior_calculator.h"
#include "gelex/infra/stats/running_stats.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{

PosteriorSummary::PosteriorSummary(RunningStatsResult result)
    : mean(std::move(result.mean)), stddev(std::move(result.stddev))
{
}

MCMCResult::MCMCResult(MCMCSamples&& samples, const BayesModel& model)
    : samples_(std::move(samples)),
      fixed_(samples_.fixed()),
      residual_(1),
      phenotype_var_(model.phenotype_variance())
{
    if (const auto* effect = model.genetic(GeneticEffectType::Add); effect)
    {
        p_freq_ = bayes::get_means(effect->X).array() / 2;
    }

    for (const auto& sample : samples_.random())
    {
        random_.emplace_back(sample);
    }

    for (const auto& sample : samples_.genetics())
    {
        genetics_.emplace_back(sample);
    }
}

void FixedSummary::compute(const FixedSamples& sample)
{
    coeffs = PosteriorSummary(sample.coeffs());
}

void RandomSummary::compute(const RandomSamples& sample)
{
    coeffs = PosteriorSummary(sample.coeffs());
    variance = PosteriorSummary(sample.variance());
}

void MixtureSummary::compute(const MixtureSamples& sample)
{
    if (mixture_proportion.size() > 0)
    {
        mixture_proportion = PosteriorSummary(sample.proportion());
    }
    if (component_variance.size() > 0)
    {
        component_variance = PosteriorSummary(sample.comp_var());
    }
    comp_probs = sample.component_probs();
    pip = Eigen::VectorXd::Ones(comp_probs.rows()) - comp_probs.col(0);
}

SignSummary::SignSummary(const SignSamples& /*samples*/) : positive_prob(1) {}

void SignSummary::compute(const SignSamples& samples)
{
    positive_prob = PosteriorSummary(samples.positive_prob());
}

void GeneticSummary::compute(const GeneticSamples& sample, double phenotype_var)
{
    coeffs = PosteriorSummary(sample.coeffs());
    variance = PosteriorSummary(sample.variance());
    heritability = PosteriorSummary(sample.heritability());

    if (mixture && sample.mixture)
    {
        mixture->compute(*sample.mixture);
    }
    if (sign && sample.sign)
    {
        sign->compute(*sample.sign);
    }

    detail::PosteriorCalculator::compute_pve(pve, coeffs.mean, phenotype_var);
}

void MCMCResult::compute()
{
    fixed_.compute(samples_.fixed());

    for (auto&& [result, sample] : std::views::zip(random_, samples_.random()))
    {
        result.compute(sample);
    }

    for (auto&& [summary, sample] :
         std::views::zip(genetics_, samples_.genetics()))
    {
        summary.compute(sample, phenotype_var_);
    }

    residual_ = PosteriorSummary(samples_.residual().variance());
}
}  // namespace gelex
