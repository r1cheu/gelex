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

#include "mcmc_reporter.h"

#include <cstddef>
#include <iterator>
#include <ranges>
#include <type_traits>
#include <variant>

#include <fmt/base.h>
#include <fmt/format.h>

#include <Eigen/Core>

#include "cli/report_printer.h"
#include "config.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/labels.h"
#include "gelex/bayes/state.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/formatter.h"
#include "gelex/infra/logging/progress_bar.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::cli
{

auto McmcReporter::on_event(const MCMCBannerEvent& /*event*/) -> void
{
    cli::printer().block(
        gelex::command_banner(PROJECT_VERSION, "Model Fitting (MCMC)"));
}

auto McmcReporter::on_event(const MCMCConfigEvent& event) -> void
{
    cli::printer().block(gelex::section("[Config]"));
    cli::printer().line(
        "  {:<12}: {}", "Method", fmt::format("{}", event.recipe_scheme));
    cli::printer().line(
        "  {:<12}: {} iters ({} burn-in, {} sampling)",
        "Chain",
        event.n_iters,
        event.n_burn_in,
        event.n_iters - event.n_burn_in);
    cli::printer().line("  {:<12}: {}", "Seed", event.seed);
}

auto McmcReporter::on_event(const FitPriorSetEvent& event) -> void
{
    prior_ = event.prior;
    FitReporter::on_event(event);
}

auto McmcReporter::on_event(const MCMCProgressEvent& event) -> void
{
    if (!init_progress_)
    {
        init_progress_ = true;
        cli::printer().block(gelex::section("[MCMC Sampling]"));
        bar_ = create_progress_bar(
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
                        std::is_same_v<Prior, bayes::SingleGeneticPrior>
                        && std::
                            is_same_v<Block, bayes::SingleGeneticBlockState>)
                    {
                        const auto mode = bayes::mode(prior_block);
                        const auto& genetic = block.state();
                        fmt::format_to(
                            std::back_inserter(stats_),
                            "{}{}:{:.3f}",
                            stats_.empty() ? "" : " | ",
                            bayes::to_heritability_label(mode),
                            genetic.heritability);
                        std::visit(
                            [&](const auto& prior_state)
                            {
                                using PriorState
                                    = std::decay_t<decltype(prior_state)>;
                                if constexpr (
                                    std::is_same_v<
                                        PriorState,
                                        bayes::
                                            SingleSharedSpikeSlabGaussianState>
                                    || std::is_same_v<
                                        PriorState,
                                        bayes::
                                            SinglePerMarkerSpikeSlabGaussianState>
                                    || std::is_same_v<
                                        PriorState,
                                        bayes::
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
                        std::is_same_v<Prior, bayes::JointGeneticPrior>
                        && std::is_same_v<Block, bayes::JointGeneticBlockState>)
                    {
                        for (const auto mode : {GeneticMode::A, GeneticMode::D})
                        {
                            const auto& genetic = block.state(mode);
                            fmt::format_to(
                                std::back_inserter(stats_),
                                "{}{}:{:.3f}",
                                stats_.empty() ? "" : " | ",
                                bayes::to_heritability_label(mode),
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
                                        bayes::JointGaussianMixtureState>
                                    || std::is_same_v<
                                        PriorState,
                                        bayes::JointHalfNormalMixtureState>)
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
                                        bayes::JointHalfNormalMixtureState>)
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

auto McmcReporter::on_event(const MCMCCompleteEvent& event) -> void
{
    cli::printer().block(gelex::section("[MCMC Complete]"));
    cli::printer().line("  {:<12}: {}", "Samples", event.samples_collected);
}

auto McmcReporter::on_event(const FitResultsSavedEvent& event) -> void
{
    cli::printer().block(
        gelex::success(
            "Results saved to '{}' (.param, .summary, .snpeff, .log)",
            event.out_prefix));
}

auto McmcReporter::on_event(const FitCheckpointSavedEvent& /*event*/) -> void
{
    bar_.before_bar->message(fmt::format(" ckpt saved"));
}

}  // namespace gelex::cli
