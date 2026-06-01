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

#include <ranges>
#include <string_view>
#include <type_traits>
#include <variant>

#include <Eigen/Core>

#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/prior_state.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/exception.h"
#include "gelex/io/mcmc/sample_writer.h"
#include "gelex/types/genetic_effect_type.h"

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

void AssignmentSamples::store(const Eigen::VectorXi& assignment)
{
    ++n_samples_;
    for (Index i = 0; i < assignment.size(); ++i)
    {
        comp_counts_(i, assignment(i)) += 1.0;
    }
}

void AssignmentSamples::store(const bayes::MixtureState& state)
{
    if (estimate_pi_)
    {
        proportion_stats_.update(state.proportion);
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
    const Eigen::VectorXi& mixture_assignment)
{
    assignment.store(mixture_assignment);
    comp_var_stats_.update(component.gebv_var);
}

void ComponentSamples::store(
    const bayes::ComponentState& component,
    const bayes::MixtureState& mixture)
{
    assignment.store(mixture);
    comp_var_stats_.update(component.gebv_var);
}

GeneticSamples::GeneticSamples(
    const bayes::GeneticDesign& design,
    GeneticMode mode)
    : type(mode), n_coeffs_(design.X.cols())
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
}

Samples::Samples(
    const BayesModel& model,
    const bayes::BayesPrior& prior,
    const BayesState& state,
    std::string_view sample_prefix,
    Eigen::Index n_records)
    : fixed_(model.fixed())
{
    if (const auto& designs = model.random(); !designs.empty())
    {
        random_.reserve(designs.size());
        for (const auto& design : designs)
        {
            random_.emplace_back(design);
        }
    }

    genetics_.reserve(model.genetics().size());
    for (const auto& block : prior.genetics())
    {
        std::visit(
            [&](const auto& genetic_prior)
            {
                if constexpr (
                    std::is_same_v<
                        std::decay_t<decltype(genetic_prior)>,
                        bayes::SingleGeneticPrior>)
                {
                    const auto mode = bayes::mode(genetic_prior);
                    const auto* design = model.genetic(mode);
                    if (design == nullptr)
                    {
                        throw GelexException(
                            "Samples: missing single genetic design");
                    }
                    genetics_.emplace_back(*design, mode);
                }
                else
                {
                    const auto* additive = model.genetic(GeneticMode::A);
                    if (additive == nullptr)
                    {
                        throw GelexException(
                            "Samples: missing additive genetic design");
                    }
                    genetics_.emplace_back(*additive, GeneticMode::A);

                    const auto* dominance = model.genetic(GeneticMode::D);
                    if (dominance == nullptr)
                    {
                        throw GelexException(
                            "Samples: missing dominance genetic design");
                    }
                    genetics_.emplace_back(*dominance, GeneticMode::D);
                }
            },
            block);
    }

    if (!sample_prefix.empty())
    {
        throw GelexException(
            "MCMC sample writer is not implemented after Bayes prior/state "
            "cleanup");
    }
    static_cast<void>(state);
    static_cast<void>(n_records);
}

void Samples::store(const BayesState& states)
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
