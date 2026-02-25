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

#include "gelex/infra/detail/indicator.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/grm_event.h"
#include "gelex/infra/utils/formatter.h"

namespace gelex::cli
{

GrmReporter::GrmReporter()
    : logger_(gelex::logging::get()), info_(detail::create_progress_info())
{
}

auto GrmReporter::on_event(const GrmConfigLoadedEvent& event) const -> void
{
    logger_->info(gelex::section("[GRM Parameters]"));
    logger_->info(gelex::task("Method     : {}", event.method));
    logger_->info(gelex::task("Mode       : {}", event.mode));
    logger_->info(gelex::task("LOCO       : {}", event.do_loco ? "yes" : "no"));
    logger_->info(gelex::task("Chunk Size : {}", event.chunk_size));
    logger_->info(gelex::task("Threads    : {}", event.threads));
    logger_->info("");
}

auto GrmReporter::on_event(const GrmDataLoadedEvent& event) const -> void
{
    logger_->info(gelex::section("[Loading Data]"));
    logger_->info(gelex::success("Samples    : {} samples", event.num_samples));
    logger_->info(gelex::success("SNPs       : {} markers", event.num_snps));
    logger_->info("");
}

auto GrmReporter::on_event(const GrmProgressEvent& event) -> void
{
    if (!init_progress_)
    {
        init_progress_ = true;
        info_.display->show();
    }

    info_.progress_info->message(
        fmt::format(
            " Computing GRM {}/{} SNPs...",
            gelex::AbbrNumber(event.current),
            gelex::AbbrNumber(event.total)));
    if (event.done)
    {
        info_.display->done();
        logger_->info("");
    }
}

auto GrmReporter::on_event(const GrmFilesWrittenEvent& event) const -> void
{
    logger_->info(gelex::named_section("Computation Summary", 70));
    logger_->info(
        gelex::success("{:<14}: {}", "Time elapsed", event.time_elapsed));
    logger_->info("  Total Files : {}", event.file_paths.size());
    logger_->info("  Output Dir  : {}", event.output_dir);
    logger_->info("  Pattern     : {}", event.file_pattern);
    logger_->info(gelex::separator());
    for (const auto& path : event.file_paths)
    {
        logger_->info(gelex::success("{}", path));
    }
}

}  // namespace gelex::cli
