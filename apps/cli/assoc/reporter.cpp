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

#include <cstddef>
#include <string_view>
#include <vector>

#include <fmt/color.h>
#include <fmt/format.h>

#include "cli/formatter.h"
#include "cli/reml_reporter.h"
#include "cli/report_printer.h"
#include "gelex/algo/gwas/assoc_type.h"
#include "gelex/infra/logging/progress_bar.h"
#include "version.h"

namespace cli
{

AssocReporter::AssocReporter() : eta_(1) {}

auto AssocReporter::show_banner() const -> void
{
    cli::printer().block(
        gelex::command_banner(PROJECT_VERSION, "GWAS Analysis"));
}

auto AssocReporter::show_config(
    gelex::GeneticMode mode,
    gelex::AssocType test_type,
    bool loco,
    gelex::GenotypeMethod geno_method,
    int max_iter,
    double tol) const -> void
{
    cli::printer().block(gelex::section("[Config]"));
    cli::printer().line("  {:<12}: {}", "Mode", mode);
    cli::printer().line(
        "  {:<12}: {}",
        "Test",
        test_type == gelex::AssocType::Single ? "Single" : "Joint");
    cli::printer().line("  {:<12}: {}", "LOCO", loco ? "Yes" : "No");
    cli::printer().line("  {:<12}: {}", "Geno Method", geno_method);
    cli::printer().line("  {:<12}: {}", "Max Iter", max_iter);
    cli::printer().line("  {:<12}: {}", "Tolerance", tol);
}

auto AssocReporter::show_reml_started(std::string_view chr_name) const -> void
{
    if (chr_name.empty())
    {
        cli::printer().block(gelex::section("[Variance Component Estimation]"));
    }
    else
    {
        cli::printer().block(
            gelex::section(
                "[Variance Component Estimation — Chr {}]", chr_name));
    }
}

auto AssocReporter::start_scan(size_t total_snps, int chunk_size, bool loco)
    -> void
{
    eta_.reset(total_snps);

    cli::printer().block(gelex::section("[Association Scan]"));
    cli::printer().line("   SNPs to test : {}", total_snps);
    cli::printer().line("   Chunk size   : {}", chunk_size);
    if (loco)
    {
        cli::printer().line("   Mode         : LOCO");
    }

    bar_ = gelex::create_progress_bar(progress_, total_snps);
    bar_.display->show();
    bar_active_ = true;
}

auto AssocReporter::update_scan_progress(size_t current, size_t total) -> void
{
    progress_ = current;
    if (bar_.after_bar)
    {
        bar_.after_bar->message(
            fmt::format(
                "{:.1f}% ({}/{}) | ETA: {}",
                static_cast<double>(current) / static_cast<double>(total)
                    * 100.0,
                gelex::AbbrNumber(current),
                gelex::AbbrNumber(total),
                eta_.get_eta(current)));
    }
}

auto AssocReporter::show_loco_phase(
    std::string_view chr_name,
    std::string_view phase) -> void
{
    if (bar_.before_bar)
    {
        auto color
            = phase == "REML" ? fmt::color::yellow : fmt::color::light_green;
        bar_.before_bar->message(
            fmt::format(
                " {} [Chr {}]",
                fmt::format(fmt::fg(color), "{}", phase),
                chr_name));
    }
}

auto AssocReporter::show_loco_reml_summary(
    const std::vector<gelex::LocoRemlResult>& results) -> void
{
    if (bar_active_)
    {
        bar_.display->done();
        bar_active_ = false;
        cli::printer().on_progress_finished();
    }
    cli::print_loco_reml_summary(results);
}

auto AssocReporter::show_complete(std::string_view out_prefix) -> void
{
    if (bar_active_)
    {
        bar_.display->done();
        bar_active_ = false;
        cli::printer().on_progress_finished();
    }
    cli::printer().block(
        gelex::success("Results saved to : {}.gwas.tsv", out_prefix));
}

}  // namespace cli
