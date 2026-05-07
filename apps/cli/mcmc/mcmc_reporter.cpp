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

#include <iterator>
#include <ranges>
#include <span>

#include <fmt/format.h>

#include "config.h"
#include "gelex/algo/infer/mcmc/result.h"
#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/formatter.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/method.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::cli
{

namespace
{
const int kTableWidth = 40;
}  // namespace

McmcReporter::McmcReporter() : FitReporter() {}

auto McmcReporter::on_event(const MCMCBannerEvent& /*event*/) const -> void
{
    logger_->info(
        gelex::command_banner(PROJECT_VERSION, "Model Fitting (MCMC)"));
    logger_->info("");
}

auto McmcReporter::on_event(const MCMCConfigEvent& event) const -> void
{
    logger_->info(gelex::section("[Config]"));
    logger_->info("  {:<12}: {}", "Method", fmt::format("{}", event.method));
    logger_->info(
        "  {:<12}: {} iters ({} burn-in, {} sampling)",
        "Chain",
        event.n_iters,
        event.n_burn_in,
        event.n_iters - event.n_burn_in);
    logger_->info("  {:<12}: {}", "Seed", event.seed);
    logger_->info("");
}

auto McmcReporter::on_event(const MCMCProgressEvent& event) -> void
{
    if (!init_progress_)
    {
        init_progress_ = true;
        logger_->info("");
        logger_->info(gelex::section("[MCMC Sampling]"));
        bar_ = create_progress_bar(
            iter_, event.total, "{bar} {value}/{total} [{speed:.1f}/s]");
        bar_.display->show();
    }

    if (event.done)
    {
        bar_.display->done();
        logger_->info("");
        return;
    }

    iter_ = event.current;
    bar_.before_bar->message("");

    const auto* state = event.state;
    stats_.clear();
    for (const auto& gen : state->genetics())
    {
        fmt::format_to(
            std::back_inserter(stats_),
            "{}{}:{:.3f}",
            stats_.empty() ? "" : " | ",
            genetic_mode::to_heritability_label(gen.type),
            gen.heritability);
    }
    fmt::format_to(
        std::back_inserter(stats_),
        "{}σ²_e: {:.3f}",
        stats_.empty() ? "" : " | ",
        state->residual().variance);

    const auto* dom = state->genetic(GeneticMode::D);
    if (dom != nullptr && dom->sign)
    {
        fmt::format_to(
            std::back_inserter(stats_),
            " | p: {:.3f}",
            dom->sign->proportion(1));
    }

    if (bar_.after_bar)
    {
        bar_.after_bar->message(stats_);
    }
}

auto McmcReporter::on_event(const MCMCCompleteEvent& event) const -> void
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
    std::ptrdiff_t samples_collected) const -> void
{
    const auto& fixed = result.fixed();

    logger_->info("");
    logger_->info(gelex::section("[Posterior Summary]"));
    logger_->info("  Samples collected per parameter: {}", samples_collected);
    logger_->info("");

    logger_->info("  {:<8} {:>8} {:>8}", "Parameter", "Mean", "SD");
    logger_->info(gelex::table_separator(kTableWidth));

    fixed.for_each_term([&](const std::string& term, Eigen::Index i)
                        { print_summary_row(term, fixed.coeffs, i); });
}

auto McmcReporter::print_genetic_summary(
    const GeneticSummary* summary,
    const bayes::GeneticEffect* /*effect*/,
    GeneticMode type) const -> void
{
    if (summary == nullptr)
    {
        return;
    }

    std::string h_name{genetic_mode::to_heritability_label(type)};

    logger_->info(
        gelex::named_section(fmt::format("{}", type), kTableWidth, 2));
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

auto McmcReporter::print_residual_summary(const mcmc::Result& result) const
    -> void
{
    logger_->info(gelex::named_section("Residual", kTableWidth, 2));
    print_summary_row("σ²", result.residual());
    logger_->info(gelex::table_separator(kTableWidth));
    logger_->info("");
}

}  // namespace gelex::cli
