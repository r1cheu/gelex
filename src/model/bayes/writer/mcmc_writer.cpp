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

#include "gelex/model/bayes/writer/mcmc_writer.h"

#include <fmt/format.h>
#include <cstdint>
#include <ranges>
#include <string_view>
#include <variant>

#include <Eigen/Core>

#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/states.h"

namespace gelex
{

namespace
{

auto has_group_prior(const bayes::MarkerPrior& marker) -> bool
{
    return !std::holds_alternative<bayes::ContinuousPrior>(marker);
}

auto group_prior_estimate_pi(const bayes::MarkerPrior& marker) -> bool
{
    if (const auto* sp = std::get_if<bayes::SpikePrior>(&marker))
    {
        return sp->proportion.estimate;
    }
    return std::get<bayes::MixturePrior>(marker).proportion.estimate;
}

auto group_prior_n_proportions(const bayes::MarkerPrior& marker) -> Eigen::Index
{
    if (const auto* sp = std::get_if<bayes::SpikePrior>(&marker))
    {
        return sp->proportion.init.size();
    }
    return std::get<bayes::MixturePrior>(marker).proportion.init.size();
}

}  // namespace

mcmc::Writer::Writer(
    const BayesModel& model,
    const bayes::Priors& priors,
    std::string_view prefix,
    Eigen::Index n_records)
    : writer_(fmt::format("{}.samples", prefix))
{
    const auto cols = n_records;

    // Fixed
    fixed_coeffs_ = writer_.reserve<double>(
        fmt::format("{}/coeff", EffectType::fixed()),
        model.fixed().X.cols(),
        cols);
    writer_.write_strings(
        fmt::format("{}/names", EffectType::fixed()), model.fixed().names);

    for (const auto& [name, levels] :
         std::views::zip(model.fixed().names, model.fixed().levels))
    {
        if (!levels)
        {
            writer_.write_strings(
                fmt::format("{}/levels/{}", EffectType::fixed(), name),
                std::views::single(std::string_view{name}));
        }
        else
        {
            writer_.write_strings(
                fmt::format("{}/levels/{}", EffectType::fixed(), name),
                *levels);
        }
    }

    // Random
    for (const auto& random : model.random())
    {
        std::string_view name = random.name;
        writer_.write_strings(
            fmt::format("{}/levels/{}", EffectType::random(), name),
            random.levels);
        auto coeffs = writer_.reserve<double>(
            fmt::format("{}/coeff/{}", EffectType::random(), name),
            random.X.cols(),
            cols);
        auto variance = writer_.reserve<double>(
            fmt::format("{}/variance/{}", EffectType::random(), name), 1, cols);
        random_.push_back({.coeffs = coeffs, .variance = variance});
    }

    // Genetic
    for (const auto& effect : model.genetics())
    {
        auto sect = EffectType::from_genetic(effect.type);
        const auto n_snps = bayes::get_cols(effect.X);
        const auto* prior = priors.genetic(effect.type);

        GeneticHandles gh;
        gh.section_effect = sect;
        gh.coeffs = writer_.reserve<double>(
            fmt::format("{}/coeff", sect), n_snps, cols);
        gh.variance = writer_.reserve<double>(
            fmt::format("{}/variance", sect), 1, cols);

        if ((prior != nullptr) && has_group_prior(prior->marker))
        {
            gh.mixture_tracker = writer_.reserve<int8_t>(
                fmt::format("{}/group/assignment", sect), n_snps, cols);
            if (group_prior_estimate_pi(prior->marker))
            {
                const auto n_pi = group_prior_n_proportions(prior->marker);
                gh.pi = writer_.reserve<double>(
                    fmt::format("{}/group/proportion", sect), n_pi, cols);
            }
        }
        if ((prior != nullptr) && prior->sign)
        {
            gh.sign_tracker = writer_.reserve<int8_t>(
                fmt::format("{}/sign/assignment", sect), n_snps, cols);
            gh.positive_prob = writer_.reserve<double>(
                fmt::format("{}/sign/proportion", sect), 1, cols);
        }

        genetic_.push_back(gh);
    }

    // Residual
    residual_variance_ = writer_.reserve<double>(
        fmt::format("{}/variance", EffectType::residual()), 1, cols);
}

void mcmc::Writer::write(const mcmc::State& state)
{
    // Fixed
    writer_.write(fixed_coeffs_, state.fixed().coeffs);

    // Random
    for (auto&& [handles, rs] : std::views::zip(random_, state.random()))
    {
        writer_.write(handles.coeffs, rs.coeffs);
        writer_.write(handles.variance, rs.variance);
    }

    // Genetic
    for (auto&& [gh, gs] : std::views::zip(genetic_, state.genetics()))
    {
        writer_.write(gh.coeffs, gs.coeffs);
        writer_.write(gh.variance, gs.variance);

        if (gh.mixture_tracker && gs.group)
        {
            std::visit(
                [&](const auto& alloc)
                {
                    using T = std::decay_t<decltype(alloc)>;
                    const auto& assignment = [&]() -> const bayes::Assignment&
                    {
                        if constexpr (std::is_same_v<T, bayes::Assignment>)
                        {
                            return alloc;
                        }
                        else
                        {
                            return alloc.assignment;
                        }
                    }();
                    writer_.write(*gh.mixture_tracker, assignment.tracker);
                    if (gh.pi)
                    {
                        writer_.write(*gh.pi, assignment.proportion);
                    }
                },
                *gs.group);
        }
        if (gh.sign_tracker && gs.sign)
        {
            writer_.write(*gh.sign_tracker, gs.sign->tracker);
            const double pos_prob = gs.sign->proportion(1);
            writer_.write(*gh.positive_prob, pos_prob);
        }
    }

    // Residual
    writer_.write(residual_variance_, state.residual().variance);
}

void mcmc::Writer::finalize()
{
    writer_.finalize();
}

}  // namespace gelex
