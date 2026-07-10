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
#include <fmt/color.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <ranges>
#include <string_view>
#include <vector>

#include "gelex/freq/model.h"

#include "cli/formatter.h"
#include "cli/progress_bar.h"
#include "cli/reml_reporter.h"
#include "cli/report_printer.h"

namespace cli
{

AssocReporter::AssocReporter() : eta_(1) {}

auto AssocReporter::show_dataset_summary(
    const gelex::FreqModel& model,
    Eigen::Index n_snps) -> void
{
    auto& p = cli::printer();
    p.block(cli::section("Dataset Summary:"));
    p.line(cli::field("Analyzed Samples", "{}", model.num_individuals()));
    p.line(cli::field("Covariates", "{}", model.fixed().X.cols()));
    auto names = model.random()
                 | std::views::transform([](const auto& d)
                                         { return std::string_view(d.name); });
    p.line(cli::field("Random effects", "{}", fmt::join(names, ", ")));
    p.line(cli::field("SNPs", "{} markers", n_snps));
}

auto AssocReporter::start_scan(size_t total_snps, int chunk_size, bool loco)
    -> void
{
    eta_.reset(total_snps);

    cli::printer().block(cli::section("Association Scan:"));
    cli::printer().line("   SNPs to test : {}", total_snps);
    cli::printer().line("   Chunk size   : {}", chunk_size);
    if (loco)
    {
        cli::printer().line("   Mode         : LOCO");
    }

    bar_ = cli::create_progress_bar(progress_, total_snps);
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
                cli::AbbrNumber(current),
                cli::AbbrNumber(total),
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

auto AssocReporter::finish_scan() -> void
{
    if (bar_active_)
    {
        bar_.display->done();
        bar_active_ = false;
        cli::printer().on_progress_finished();
    }
}

}  // namespace cli
