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

#include <Eigen/Core>
#include <cstddef>
#include <iterator>
#include <type_traits>
#include <variant>

#include <fmt/base.h>
#include <fmt/format.h>

#include "cli/report_printer.h"
#include "config.h"
#include "gelex/algo/infer/mcmc/result.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/formatter.h"
#include "gelex/infra/logging/progress_bar.h"
#include "gelex/model/bayes/labels.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/state.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::cli
{

namespace
{
const int kTableWidth = 40;
}  // namespace

auto McmcReporter::on_event(const MCMCBannerEvent& /*event*/) -> void
{
    cli::printer().block(
        gelex::command_banner(PROJECT_VERSION, "Model Fitting (MCMC)"));
}

auto McmcReporter::on_event(const MCMCConfigEvent& event) -> void
{
    cli::printer().block(gelex::section("[Config]"));
    cli::printer().line(
        "  {:<12}: {}", "Method", fmt::format("{}", event.preset));
    cli::printer().line(
        "  {:<12}: {} iters ({} burn-in, {} sampling)",
        "Chain",
        event.n_iters,
        event.n_burn_in,
        event.n_iters - event.n_burn_in);
    cli::printer().line("  {:<12}: {}", "Seed", event.seed);
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
    for (const auto& gen : state->genetics())
    {
        std::visit(
            [&](const auto& block)
            {
                using Block = std::decay_t<decltype(block)>;
                if constexpr (
                    std::is_same_v<Block, bayes::SingleGeneticBlockState>)
                {
                    const auto& genetic = block.state();
                    fmt::format_to(
                        std::back_inserter(stats_),
                        "{}{}:{:.3f}",
                        stats_.empty() ? "" : " | ",
                        bayes::to_heritability_label(genetic.type),
                        genetic.heritability);
                }
                else
                {
                    for (const auto mode : {GeneticMode::A, GeneticMode::D})
                    {
                        const auto& genetic = block.state(mode);
                        fmt::format_to(
                            std::back_inserter(stats_),
                            "{}{}:{:.3f}",
                            stats_.empty() ? "" : " | ",
                            bayes::to_heritability_label(genetic.type),
                            genetic.heritability);
                    }
                }
            },
            gen);
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
    print_fixed_summary(*event.result, event.samples_collected);
    for (const auto& summary : event.result->genetics())
    {
        print_genetic_summary(
            &summary, event.model->genetic(summary.type), summary.type);
    }
    print_residual_summary(*event.result);
}

auto McmcReporter::on_event(const FitCheckpointSavedEvent& /*event*/) -> void
{
    bar_.before_bar->message(fmt::format(" ckpt saved"));
}

auto McmcReporter::print_fixed_summary(
    const mcmc::Result& result,
    std::ptrdiff_t samples_collected) -> void
{
    const auto& fixed = result.fixed();
    auto& p = cli::printer();

    p.block(gelex::section("[Posterior Summary]"));
    p.line("  Samples collected per parameter: {}", samples_collected);
    p.line("  {:<8} {:>8} {:>8}", "Parameter", "Mean", "SD");
    p.line(gelex::table_separator(kTableWidth));

    fixed.for_each_term([&](const std::string& term, Eigen::Index i)
                        { print_summary_row(term, fixed.coeffs, i); });
}

auto McmcReporter::print_genetic_summary(
    const GeneticSummary* summary,
    const bayes::GeneticDesign* /*design*/,
    GeneticMode type) -> void
{
    if (summary == nullptr)
    {
        return;
    }

    std::string h_name{bayes::to_heritability_label(type)};
    auto& p = cli::printer();

    p.line(gelex::named_section(fmt::format("{}", type), kTableWidth, 2));
    print_summary_row("σ²", summary->variance);
    print_summary_row(h_name, summary->heritability);

    if (summary->group)
    {
        const auto& base = assignment(*summary->group);
        for (Eigen::Index i = 0; i < base.mixture_proportion.size(); ++i)
        {
            print_summary_row(
                fmt::format("π[{}]", i), base.mixture_proportion, i);
        }
        if (const auto* comp = std::get_if<ComponentSummary>(&*summary->group))
        {
            for (Eigen::Index i = 0; i < comp->component_variance.size(); ++i)
            {
                print_summary_row(
                    fmt::format("σ²[{}]", i + 1), comp->component_variance, i);
            }
        }
    }

    if (summary->sign)
    {
        print_summary_row("p(+)", summary->sign->mixture_proportion, 1);
    }
}

auto McmcReporter::print_residual_summary(const mcmc::Result& result) -> void
{
    auto& p = cli::printer();
    p.line(gelex::named_section("Residual", kTableWidth, 2));
    print_summary_row("σ²", result.residual());
    p.line(gelex::table_separator(kTableWidth));
}

}  // namespace gelex::cli
