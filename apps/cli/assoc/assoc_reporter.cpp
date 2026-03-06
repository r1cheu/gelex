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

#include "assoc_reporter.h"

#include <fmt/format.h>

#include "cli/reml_reporter.h"
#include "config.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/assoc_event.h"
#include "gelex/infra/utils/formatter.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::cli
{

AssocReporter::AssocReporter() : logger_(gelex::logging::get()), eta_(1) {}

auto AssocReporter::on_event(const AssocConfigLoadedEvent& event) const -> void
{
    logger_->info(gelex::command_banner(PROJECT_VERSION, "GWAS Analysis"));
    logger_->info("");
    logger_->info(gelex::section("[Config]"));
    logger_->info(
        "  {:<12}: {}",
        "Model",
        event.model_type == gelex::ModelType::A ? "Additive" : "Dominance");
    logger_->info("  {:<12}: {}", "LOCO", event.loco ? "Yes" : "No");

    logger_->info("  {:<12}: {}", "Geno Method", event.geno_method);

    logger_->info("  {:<12}: {}", "Max Iter", event.max_iter);
    logger_->info("  {:<12}: {}", "Tolerance", event.tol);
    logger_->info("");
}

auto AssocReporter::on_event(const AssocRemlStartedEvent& event) const -> void
{
    logger_->info(gelex::section(""));
    if (event.chr_name.empty())
    {
        logger_->info(gelex::section("[Variance Component Estimation]"));
    }
    else
    {
        logger_->info(
            gelex::section(
                "[Variance Component Estimation — Chr {}]", event.chr_name));
    }
}

auto AssocReporter::on_event(const AssocScanSummaryEvent& event) -> void
{
    eta_.reset(event.total_snps);

    logger_->info("");
    logger_->info(gelex::section("[Association Scan]"));
    logger_->info(gelex::task("SNPs to test : {}", event.total_snps));
    logger_->info(gelex::task("Chunk size   : {}", event.chunk_size));
    if (event.loco)
    {
        logger_->info(gelex::task("Mode         : LOCO"));
    }
    logger_->info("");

    bar_ = detail::create_progress_bar(progress_, event.total_snps);
    bar_.display->show();
    bar_active_ = true;
}

auto AssocReporter::on_event(const AssocScanProgressEvent& event) -> void
{
    progress_ = event.current;
    if (bar_.after_bar)
    {
        bar_.after_bar->message(
            fmt::format(
                "{:.1f}% ({}/{}) | ETA: {}",
                static_cast<double>(event.current)
                    / static_cast<double>(event.total) * 100.0,
                AbbrNumber(event.current),
                AbbrNumber(event.total),
                eta_.get_eta(event.current)));
    }
}

auto AssocReporter::on_event(const AssocLocoPhaseEvent& event) -> void
{
    if (bar_.before_bar)
    {
        auto color = event.phase == "REML" ? fmt::color::yellow
                                           : fmt::color::light_green;
        bar_.before_bar->message(
            fmt::format(
                " {} [Chr {}]",
                fmt::format(fmt::fg(color), "{}", event.phase),
                event.chr_name));
    }
}

auto AssocReporter::on_event(const AssocLocoRemlSummaryEvent& event) -> void
{
    if (bar_active_)
    {
        bar_.display->done();
        bar_active_ = false;
    }
    cli::print_loco_reml_summary(event.results);
}

auto AssocReporter::on_event(const AssocCompleteEvent& event) -> void
{
    if (bar_active_)
    {
        bar_.display->done();
        bar_active_ = false;
    }

    logger_->info("");
    logger_->info(
        gelex::success("Results saved to : {}.gwas.tsv", event.out_prefix));
}

}  // namespace gelex::cli
