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

#include "fit_reporter.h"

#include <iterator>

#include <fmt/format.h>

#include "config.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/utils/formatter.h"

namespace gelex::cli
{

FitReporter::FitReporter() : logger_(gelex::logging::get()) {}

auto FitReporter::on_event(const FitConfigLoadedEvent& event) const -> void
{
    logger_->info(
        gelex::command_banner(PROJECT_VERSION, "Model Fitting (MCMC)"));
    logger_->info("");
    logger_->info(gelex::section("[Config]"));
    logger_->info(
        "  {:<12}: {} ({})",
        "Method",
        fmt::format("{}", event.method),
        event.use_dominance ? "Additive + Dominance" : "Additive");
    logger_->info(
        "  {:<12}: {} iters ({} burn-in, {} sampling)",
        "Chain",
        event.n_iters,
        event.n_burnin,
        event.n_iters - event.n_burnin);
    logger_->info("  {:<12}: {}", "Seed", event.seed);
    logger_->info("  {:<12}: {}", "Threads", event.threads);
    logger_->info("");
}

auto FitReporter::on_event(const FitMcmcProgressEvent& event) -> void
{
    if (!init_progress_)
    {
        init_progress_ = true;
        bar_ = detail::create_progress_bar(
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

    stats_.clear();
    if (event.h2)
    {
        fmt::format_to(std::back_inserter(stats_), "h²: {:.3f}", *event.h2);
    }
    if (event.h2_dom)
    {
        fmt::format_to(
            std::back_inserter(stats_),
            "{}δ²: {:.3f}",
            stats_.empty() ? "" : " | ",
            *event.h2_dom);
    }
    if (event.sigma2_e)
    {
        fmt::format_to(
            std::back_inserter(stats_),
            "{}σ²_e: {:.3f}",
            stats_.empty() ? "" : " | ",
            *event.sigma2_e);
    }

    if (bar_.after_bar)
    {
        bar_.after_bar->message(stats_);
    }
}

auto FitReporter::on_event(const FitResultsSavedEvent& event) const -> void
{
    logger_->info(
        gelex::success(
            "Results saved to '{}' (.param, .snp.eff, .log)",
            event.out_prefix));
}

}  // namespace gelex::cli
