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

#include "reporter.h"

#include <Eigen/Core>
#include <cstddef>
#include <fmt/base.h>
#include <fmt/format.h>
#include <iterator>
#include <ranges>
#include <type_traits>
#include <variant>

#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/labels.h"
#include "gelex/bayes/state.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/types/genetic_mode.h"

#include "cli/formatter.h"
#include "cli/progress_bar.h"
#include "cli/report_printer.h"

namespace cli
{

auto McmcReporter::show_prior(const gelex::bayes::BayesPrior& prior) -> void
{
    prior_ = &prior;
    FitReporter::show_prior(prior);
}

auto McmcReporter::on_event(const gelex::MCMCProgressEvent& event) -> void
{
    if (!init_progress_)
    {
        init_progress_ = true;
        cli::printer().block(gelex::section("MCMC Sampling:"));
        bar_ = cli::create_progress_bar(
            iter_, event.total, "{bar} {value}/{total} [{speed:.1f}/s]");
        bar_.display->show();
    }

    if (event.done)
    {
        bar_.display->done();
        cli::printer().on_progress_finished();
        return;
    }

    iter_ = event.current;
    bar_.before_bar->message("");

    const auto* state = event.state;
    stats_.clear();
    if (prior_ != nullptr)
    {
        const auto priors = prior_->genetics();
        for (auto [i, gen] : std::views::enumerate(state->genetics()))
        {
            const auto& prior = priors[static_cast<std::size_t>(i)];
            std::visit(
                [&](const auto& prior_block, const auto& block)
                {
                    using Prior = std::decay_t<decltype(prior_block)>;
                    using Block = std::decay_t<decltype(block)>;
                    if constexpr (
                        std::is_same_v<Prior, gelex::bayes::SingleGeneticPrior>
                        && std::is_same_v<
                            Block,
                            gelex::bayes::SingleGeneticBlockState>)
                    {
                        const auto mode = gelex::bayes::mode(prior_block);
                        const auto& genetic = block.state();
                        fmt::format_to(
                            std::back_inserter(stats_),
                            "{}{}:{:.3f}",
                            stats_.empty() ? "" : " | ",
                            gelex::bayes::to_heritability_label(mode),
                            genetic.heritability);
                        std::visit(
                            [&](const auto& prior_state)
                            {
                                using PriorState
                                    = std::decay_t<decltype(prior_state)>;
                                if constexpr (
                                    std::is_same_v<
                                        PriorState,
                                        gelex::bayes::
                                            SingleSharedSpikeSlabGaussianState>
                                    || std::is_same_v<
                                        PriorState,
                                        gelex::bayes::
                                            SinglePerMarkerSpikeSlabGaussianState>
                                    || std::is_same_v<
                                        PriorState,
                                        gelex::bayes::
                                            SingleScaledMixtureGaussianState>)
                                {
                                    fmt::format_to(
                                        std::back_inserter(stats_), " | π:[");
                                    const auto& proportion
                                        = prior_state.proportion();
                                    for (Eigen::Index j = 0;
                                         j < proportion.size();
                                         ++j)
                                    {
                                        fmt::format_to(
                                            std::back_inserter(stats_),
                                            "{}{:.3f}",
                                            j == 0 ? "" : ",",
                                            proportion(j));
                                    }
                                    fmt::format_to(
                                        std::back_inserter(stats_), "]");
                                }
                            },
                            block.prior_state());
                    }
                    else if constexpr (
                        std::is_same_v<Prior, gelex::bayes::JointGeneticPrior>
                        && std::is_same_v<
                            Block,
                            gelex::bayes::JointGeneticBlockState>)
                    {
                        for (const auto mode :
                             {gelex::GeneticMode::A, gelex::GeneticMode::D})
                        {
                            const auto& genetic = block.state(mode);
                            fmt::format_to(
                                std::back_inserter(stats_),
                                "{}{}:{:.3f}",
                                stats_.empty() ? "" : " | ",
                                gelex::bayes::to_heritability_label(mode),
                                genetic.heritability);
                        }
                        std::visit(
                            [&](const auto& prior_state)
                            {
                                using PriorState
                                    = std::decay_t<decltype(prior_state)>;
                                if constexpr (
                                    std::is_same_v<
                                        PriorState,
                                        gelex::bayes::JointGaussianMixtureState>
                                    || std::is_same_v<
                                        PriorState,
                                        gelex::bayes::
                                            JointHalfNormalMixtureState>)
                                {
                                    fmt::format_to(
                                        std::back_inserter(stats_), " | π:[");
                                    const auto& proportion
                                        = prior_state.proportion();
                                    for (Eigen::Index j = 0;
                                         j < proportion.size();
                                         ++j)
                                    {
                                        fmt::format_to(
                                            std::back_inserter(stats_),
                                            "{}{:.3f}",
                                            j == 0 ? "" : ",",
                                            proportion(j));
                                    }
                                    fmt::format_to(
                                        std::back_inserter(stats_), "]");
                                }
                                if constexpr (
                                    std::is_same_v<
                                        PriorState,
                                        gelex::bayes::
                                            JointHalfNormalMixtureState>)
                                {
                                    fmt::format_to(
                                        std::back_inserter(stats_),
                                        " | p_s:{:.3f}",
                                        prior_state.dominance_sign()
                                            .positive_probability);
                                }
                            },
                            block.prior_state());
                    }
                },
                prior,
                gen);
        }
    }
    fmt::format_to(
        std::back_inserter(stats_),
        "{}σ²_e: {:.3f}",
        stats_.empty() ? "" : " | ",
        state->residual().variance);

    if (bar_.after_bar)
    {
        bar_.after_bar->message(stats_);
    }
}

auto McmcReporter::show_complete(std::ptrdiff_t samples_collected) -> void
{
    cli::printer().block(gelex::section("MCMC Complete:"));
    cli::printer().line("  {:<12}: {}", "Samples", samples_collected);
}

auto McmcReporter::on_event(const gelex::MCMCCheckpointSavedEvent& /*event*/)
    -> void
{
    bar_.before_bar->message(fmt::format(" ckpt saved"));
}

}  // namespace cli
