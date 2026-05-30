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

#include "gelex/algo/infer/mcmc/result.h"

#include <ranges>

#include <Eigen/Core>
#include <type_traits>
#include <utility>
#include <variant>

#include "gelex/algo/detail/posterior_calculator.h"
#include "gelex/algo/infer/mcmc/samples.h"
#include "gelex/bayes/model.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

mcmc::Result::Result(mcmc::Samples&& samples, const BayesModel& model)
    : samples_(std::move(samples)),
      fixed_(samples_.fixed()),
      residual_(1),
      phenotype_var_(model.phenotype_variance())
{
    if (const auto* design = model.genetic(GeneticMode::A); design)
    {
        p_freq_ = design->X.mean().array() / 2;
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

AssignmentSummary::AssignmentSummary(const mcmc::AssignmentSamples& samples)
    : pip(Eigen::VectorXd::Zero(samples.n_snps())),
      comp_probs(
          Eigen::MatrixXd::Zero(samples.n_snps(), samples.n_proportions()))
{
    if (samples.estimate_pi())
    {
        mixture_proportion = PosteriorSummary(samples.n_proportions());
    }
}

void AssignmentSummary::compute(const mcmc::AssignmentSamples& sample)
{
    if (mixture_proportion.size() > 0)
    {
        mixture_proportion = PosteriorSummary(sample.proportion());
    }
    comp_probs = sample.component_probs();
    pip = Eigen::VectorXd::Ones(comp_probs.rows()) - comp_probs.col(0);
}

ComponentSummary::ComponentSummary(const mcmc::ComponentSamples& samples)
    : assignment(samples.assignment),
      component_variance(samples.assignment.n_proportions() - 1)
{
}

void ComponentSummary::compute(const mcmc::ComponentSamples& sample)
{
    assignment.compute(sample.assignment);
    component_variance = PosteriorSummary(sample.comp_var());
}

GeneticSummary::GeneticSummary(const mcmc::GeneticSamples& samples)
    : type(samples.type),
      coeffs(samples.n_coeffs()),
      variance(1),
      heritability(1),
      pve(samples.n_coeffs())
{
    if (samples.group)
    {
        std::visit(
            [&](const auto& s)
            {
                using T = std::decay_t<decltype(s)>;
                if constexpr (std::is_same_v<T, mcmc::AssignmentSamples>)
                {
                    group.emplace(std::in_place_type<AssignmentSummary>, s);
                }
                else
                {
                    group.emplace(std::in_place_type<ComponentSummary>, s);
                }
            },
            *samples.group);
    }
    if (samples.sign)
    {
        sign.emplace(*samples.sign);
    }
}

void GeneticSummary::compute(
    const mcmc::GeneticSamples& sample,
    double phenotype_var)
{
    coeffs = PosteriorSummary(sample.coeffs());
    variance = PosteriorSummary(sample.variance());
    heritability = PosteriorSummary(sample.heritability());

    if (group && sample.group)
    {
        std::visit(
            [](auto& g, const auto& s)
            {
                using G = std::decay_t<decltype(g)>;
                using S = std::decay_t<decltype(s)>;
                if constexpr (
                    (std::is_same_v<G, AssignmentSummary>
                     && std::is_same_v<S, mcmc::AssignmentSamples>)
                    || (std::is_same_v<G, ComponentSummary>
                        && std::is_same_v<S, mcmc::ComponentSamples>))
                {
                    g.compute(s);
                }
            },
            *group,
            *sample.group);
    }
    if (sign && sample.sign)
    {
        sign->compute(*sample.sign);
    }

    posterior::detail::compute_pve(pve, coeffs.mean, phenotype_var);
}

void mcmc::Result::compute()
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
