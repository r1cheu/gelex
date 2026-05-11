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

#include "gelex/algo/infer/detail/genetic_binding.h"
#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/io/mcmc/sample_writer.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/method.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"

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
    if (!prior.mixture)
    {
        return std::nullopt;
    }
    const auto n_snps = effect.X.cols();
    const auto n_pi = prior.mixture->proportions.init.size();
    const auto estimate_pi = prior.mixture->proportions.estimate;
    if (std::holds_alternative<bayes::SpikeSlab>(prior.mixture->strategy))
    {
        return AssignmentSamples{n_snps, n_pi, estimate_pi};
    }
    return ComponentSamples{n_snps, n_pi, estimate_pi};
}

GeneticSamples::GeneticSamples(
    const bayes::GeneticEffect& effect,
    const bayes::GeneticPrior& prior,
    GeneticMode mode)
    : type(mode),
      group(make_group_samples(effect, prior)),
      n_coeffs_(effect.X.cols())
{
    const auto* spec = std::get_if<bayes::GeneticSpec>(&prior.spec);
    const bool has_sign = (spec != nullptr && spec->sign.has_value()) || [&]
    {
        if (const auto* js = std::get_if<bayes::JointSpec>(&prior.spec))
        {
            return (js->additive.mode == mode && js->additive.sign.has_value())
                   || (js->dominance.mode == mode
                       && js->dominance.sign.has_value());
        }
        return false;
    }();

    if (has_sign)
    {
        const auto n_snps = effect.X.cols();
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

Samples::Samples(
    const BayesModel& model,
    const bayes::BayesMethod& method,
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
        const auto* prior
            = infer::detail::find_prior_for_mode(method, effect.type);
        genetics_.emplace_back(effect, *prior, effect.type);
    }

    if (!sample_prefix.empty())
    {
        writer_ = std::make_unique<mcmc::Writer>(
            model, method, sample_prefix, n_records);
    }
}

void Samples::store(const mcmc::State& states)
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

void Samples::finalize()
{
    writer_.reset();
}

}  // namespace gelex::mcmc
