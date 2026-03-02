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

#include "simulater_reporter.h"

#include "config.h"
#include "gelex/infra/detail/indicator.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/simulate_event.h"
#include "gelex/infra/utils/formatter.h"

namespace gelex::cli
{

SimulaterReporter::SimulaterReporter()
    : logger_(gelex::logging::get()), info_(detail::create_progress_info())
{
}

auto SimulaterReporter::on_event(const SimulateConfigLoadedEvent& event) const
    -> void
{
    std::string mode_str = event.dom_heritability ? "AD" : "A";

    logger_->info(
        gelex::command_banner(PROJECT_VERSION, "Phenotype Simulation"));
    logger_->info("");
    logger_->info(gelex::section("[Config]"));
    logger_->info("  {:<12}: {}", "Mode", mode_str);
    logger_->info("  {:<12}: {:.4f}", "h\u00b2", event.add_heritability);
    if (event.dom_heritability)
    {
        logger_->info("  {:<12}: {:.4f}", "d\u00b2", *event.dom_heritability);
    }
    if (event.intercept != 0.0)
    {
        logger_->info("  {:<12}: {:.4f}", "Intercept", event.intercept);
    }
    logger_->info("  {:<12}: {}", "Seed", event.seed);
    logger_->info("");
}

auto SimulaterReporter::on_event(const SimulateProgressEvent& event) -> void
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

auto SimulaterReporter::on_event(const HeritabilityGeneratedEvent& event) const
    -> void
{
    logger_->info(gelex::section("[True Heritability]"));
    logger_->info(gelex::success("h²                : {:.4f}", event.additive));
    if (event.dominance)
    {
        logger_->info(
            gelex::success("δ²                : {:.4f}", *event.dominance));
    }
    logger_->info("");
}

auto SimulaterReporter::on_event(const OutputsWrittenEvent& event) const -> void
{
    logger_->info(gelex::section("[File Summary]"));
    logger_->info(
        gelex::success("{:<24}: {}", "Phenotype file", event.phenotype_path));
    logger_->info(
        gelex::success(
            "{:<24}: {}", "Snp effects file", event.snp_effect_path));
}

}  // namespace gelex::cli
