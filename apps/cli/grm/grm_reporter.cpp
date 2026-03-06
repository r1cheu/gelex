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

#include "grm_reporter.h"

#include "config.h"
#include "gelex/infra/detail/indicator.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/grm_event.h"
#include "gelex/infra/utils/formatter.h"

namespace gelex::cli
{

GrmReporter::GrmReporter() : logger_(gelex::logging::get()), eta_(1) {}

auto GrmReporter::on_event(const GrmConfigLoadedEvent& event) const -> void
{
    logger_->info(gelex::command_banner(PROJECT_VERSION, "GRM Computation"));
    logger_->info("");
    logger_->info(gelex::section("[Config]"));
    logger_->info("  {:<12}: {}", "Method", event.method);
    logger_->info("  {:<12}: {}", "Mode", fmt::format("{}", event.mode));
    logger_->info("  {:<12}: {}", "LOCO", event.do_loco ? "yes" : "no");
    logger_->info("");
}

auto GrmReporter::on_event(const GrmDataLoadedEvent& event) const -> void
{
    logger_->info(gelex::section("[Dataset Summary]"));
    logger_->info(gelex::success("Samples    : {} samples", event.num_samples));
    logger_->info(gelex::success("SNPs       : {} markers", event.num_snps));
    logger_->info("");
}

auto GrmReporter::on_event(const GrmComputeStartedEvent& event) -> void
{
    global_total_ = event.total_snps;
    accumulated_base_ = 0;
    progress_ = 0;
    eta_.reset(global_total_);

    bar_ = detail::create_progress_bar(progress_, global_total_);
    bar_.display->show();
    bar_active_ = true;
}

auto GrmReporter::on_event(const GrmProgressEvent& event) -> void
{
    if (event.done)
    {
        bar_.display->done();
        bar_active_ = false;
        logger_->info("");
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
                AbbrNumber(progress_),
                AbbrNumber(global_total_),
                eta_.get_eta(progress_)));
    }
    if (event.current == event.total)
    {
        accumulated_base_ += event.total;
    }
}

auto GrmReporter::on_event(const GrmFilesWrittenEvent& event) const -> void
{
    logger_->info(gelex::section("[File Summary]"));
    logger_->info("  Num Files : {}", event.num_files);
    logger_->info("  Output Dir  : {}", event.output_dir);
    logger_->info("  Pattern     : {}", event.file_pattern);
}

}  // namespace gelex::cli
