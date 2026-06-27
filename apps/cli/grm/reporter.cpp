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

#include <fmt/format.h>

#include "cli/formatter.h"
#include "cli/progress_bar.h"
#include "cli/report_printer.h"
#include "gelex/infra/logging/grm_event.h"

namespace cli
{

GrmReporter::GrmReporter() : eta_(1) {}

auto GrmReporter::show_data_loaded(size_t num_samples, size_t num_snps) -> void
{
    cli::printer().block(gelex::section("[Dataset Summary]"));
    cli::printer().line("   Samples    : {} samples", num_samples);
    cli::printer().line("   SNPs       : {} markers", num_snps);
}

auto GrmReporter::start_compute(size_t total_snps) -> void
{
    global_total_ = total_snps;
    accumulated_base_ = 0;
    progress_ = 0;
    eta_.reset(global_total_);

    bar_ = cli::create_progress_bar(progress_, global_total_);
    bar_.display->show();
    bar_active_ = true;
}

auto GrmReporter::on_event(const gelex::GrmProgressEvent& event) -> void
{
    if (event.done)
    {
        bar_.display->done();
        bar_active_ = false;
        cli::printer().on_progress_finished();
        return;
    }
    progress_ = accumulated_base_ + event.current;
    if (bar_.after_bar)
    {
        bar_.after_bar->message(
            fmt::format(
                "{:.1f}% ({}/{}) | ETA: {}",
                static_cast<double>(progress_)
                    / static_cast<double>(global_total_) * 100.0,
                gelex::AbbrNumber(progress_),
                gelex::AbbrNumber(global_total_),
                eta_.get_eta(progress_)));
    }
    if (event.current == event.total)
    {
        accumulated_base_ += event.total;
    }
}

auto GrmReporter::finish_progress() -> void
{
    if (!bar_active_)
    {
        return;
    }
    bar_.display->done();
    bar_active_ = false;
    cli::printer().on_progress_finished();
}

auto GrmReporter::show_files_written(
    size_t num_files,
    std::string_view output_dir,
    std::string_view file_pattern) -> void
{
    cli::printer().block(gelex::section("[File Summary]"));
    cli::printer().line("  Num Files : {}", num_files);
    cli::printer().line("  Output Dir  : {}", output_dir);
    cli::printer().line("  Pattern     : {}", file_pattern);
}

}  // namespace cli
