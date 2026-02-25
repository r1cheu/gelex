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

auto SimulaterReporter::on_event(const ParameterLoadedEvent& event) const
    -> void
{
    logger_->info(gelex::section("[Simulation Parameters]"));
    logger_->info(
        gelex::task("Heritability (h2)      : {:.2f}", event.add_heritability));
    if (event.dom_heritability)
    {
        logger_->info(
            gelex::task(
                "Dom-Heritability (d2)  : {:.2f}", *event.dom_heritability));
    }
    if (event.intercept != 0.0)
    {
        logger_->info(
            gelex::task("Intercept             : {:.2f}", event.intercept));
    }
    logger_->info(gelex::task("Seed                  : {}", event.seed));
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
    logger_->info(gelex::section("[Generating Phenotypes]"));
    logger_->info(
        gelex::task("True h²                : {:.4f}", event.additive));
    if (event.dominance)
    {
        logger_->info(
            gelex::task("True d²                : {:.4f}", *event.dominance));
    }
    logger_->info("");
}

auto SimulaterReporter::on_event(const OutputsWrittenEvent& event) const -> void
{
    logger_->info(
        gelex::success(
            "{:<24}: {}", "Phenotypes saved to", event.phenotype_path));
    logger_->info(
        gelex::success(
            "{:<24}: {}", "Snp effects saved to", event.snp_effect_path));
}

}  // namespace gelex::cli
