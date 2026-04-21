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

#include "simulator_reporter.h"

#include "config.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/formatter.h"
#include "gelex/infra/logging/progress_bar.h"
#include "gelex/infra/logging/simulate_event.h"

namespace gelex::cli
{

SimulatorReporter::SimulatorReporter()
    : logger_(gelex::logging::get()), info_(create_progress_info())
{
}

auto SimulatorReporter::on_event(const SimulateBannerEvent& /*event*/) const
    -> void
{
    logger_->info(
        gelex::command_banner(PROJECT_VERSION, "Phenotype Simulation"));
    logger_->info("");
}

auto SimulatorReporter::on_event(const SimulateConfigLoadedEvent& event) const
    -> void
{
    std::string mode_str;
    if (event.add_heritability && event.dom_heritability)
    {
        mode_str = "AD";
    }
    else if (event.add_heritability)
    {
        mode_str = "A";
    }
    else
    {
        mode_str = "D";
    }

    logger_->info(gelex::section("[Config]"));
    logger_->info("  {:<12}: {}", "Mode", mode_str);
    if (event.add_heritability)
    {
        logger_->info("  {:<12}: {:.4f}", "h\u00b2", *event.add_heritability);
    }
    if (event.dom_heritability)
    {
        logger_->info("  {:<12}: {:.4f}", "d\u00b2", *event.dom_heritability);
    }
    logger_->info("  {:<12}: {}", "Seed", event.seed);
    logger_->info("");
}

auto SimulatorReporter::on_event(
    const SimulateVarianceSummaryEvent& event) const -> void
{
    logger_->info(gelex::section("[Realized]"));
    if (event.realized_h2)
    {
        logger_->info("  {:<12}: {:.4f}", "h\u00b2", *event.realized_h2);
    }
    if (event.realized_d2)
    {
        logger_->info("  {:<12}: {:.4f}", "d\u00b2", *event.realized_d2);
    }
    logger_->info("");
}

auto SimulatorReporter::on_event(const SimulateProgressEvent& event) -> void
{
    if (!init_progress_)
    {
        init_progress_ = true;
        info_.display->show();
    }

    info_.progress_info->message(
        fmt::format(
            " Simulating {}/{} SNPs...",
            gelex::AbbrNumber(event.current),
            gelex::AbbrNumber(event.total)));
    if (event.done)
    {
        info_.display->done();
        logger_->info("");
    }
}

}  // namespace gelex::cli
