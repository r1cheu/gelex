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

#include <fmt/color.h>
#include <fmt/format.h>

#include "cli/reml_reporter.h"
#include "cli/report_printer.h"
#include "gelex/algo/gwas/assoc_type.h"
#include "gelex/infra/logging/assoc_event.h"
#include "gelex/infra/logging/formatter.h"
#include "gelex/infra/logging/progress_bar.h"
#include "version.h"

namespace cli
{

AssocReporter::AssocReporter() : eta_(1) {}

auto AssocReporter::on_event(const gelex::AssocBannerEvent& /*event*/) const
    -> void
{
    cli::printer().block(
        gelex::command_banner(PROJECT_VERSION, "GWAS Analysis"));
}

auto AssocReporter::on_event(const gelex::AssocConfigLoadedEvent& event) const
    -> void
{
    cli::printer().block(gelex::section("[Config]"));
    cli::printer().line("  {:<12}: {}", "Mode", event.mode);
    cli::printer().line(
        "  {:<12}: {}",
        "Test",
        event.test_type == gelex::AssocType::Single ? "Single" : "Joint");
    cli::printer().line("  {:<12}: {}", "LOCO", event.loco ? "Yes" : "No");
    cli::printer().line("  {:<12}: {}", "Geno Method", event.geno_method);
    cli::printer().line("  {:<12}: {}", "Max Iter", event.max_iter);
    cli::printer().line("  {:<12}: {}", "Tolerance", event.tol);
}

auto AssocReporter::on_event(const gelex::AssocRemlStartedEvent& event) const
    -> void
{
    if (event.chr_name.empty())
    {
        cli::printer().block(gelex::section("[Variance Component Estimation]"));
    }
    else
    {
        cli::printer().block(
            gelex::section(
                "[Variance Component Estimation — Chr {}]", event.chr_name));
    }
}

auto AssocReporter::on_event(const gelex::AssocScanSummaryEvent& event) -> void
{
    eta_.reset(event.total_snps);

    cli::printer().block(gelex::section("[Association Scan]"));
    cli::printer().line("   SNPs to test : {}", event.total_snps);
    cli::printer().line("   Chunk size   : {}", event.chunk_size);
    if (event.loco)
    {
        cli::printer().line("   Mode         : LOCO");
    }

    bar_ = gelex::create_progress_bar(progress_, event.total_snps);
    bar_.display->show();
    bar_active_ = true;
}

auto AssocReporter::on_event(const gelex::AssocScanProgressEvent& event) -> void
{
    progress_ = event.current;
    if (bar_.after_bar)
    {
        bar_.after_bar->message(
            fmt::format(
                "{:.1f}% ({}/{}) | ETA: {}",
                static_cast<double>(event.current)
                    / static_cast<double>(event.total) * 100.0,
                gelex::AbbrNumber(event.current),
                gelex::AbbrNumber(event.total),
                eta_.get_eta(event.current)));
    }
}

auto AssocReporter::on_event(const gelex::AssocLocoPhaseEvent& event) -> void
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

auto AssocReporter::on_event(const gelex::AssocLocoRemlSummaryEvent& event)
    -> void
{
    if (bar_active_)
    {
        bar_.display->done();
        bar_active_ = false;
        cli::printer().on_progress_finished();
    }
    cli::print_loco_reml_summary(event.results);
}

auto AssocReporter::on_event(const gelex::AssocCompleteEvent& event) -> void
{
    if (bar_active_)
    {
        bar_.display->done();
        bar_active_ = false;
        cli::printer().on_progress_finished();
    }
    cli::printer().block(
        gelex::success("Results saved to : {}.gwas.tsv", event.out_prefix));
}

}  // namespace cli
