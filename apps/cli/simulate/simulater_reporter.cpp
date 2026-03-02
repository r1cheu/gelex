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
    std::string title
        = fmt::format("gelex v{} :: Phenotype Simulation", PROJECT_VERSION);
    std::string mode_str
        = event.dom_heritability ? "Additive + Dominance" : "Additive";
    std::vector<std::pair<std::string, std::string>> items
        = {{"Mode", mode_str},
           {"h\u00b2", fmt::format("{:.4f}", event.add_heritability)}};
    if (event.dom_heritability)
    {
        items.emplace_back(
            "d\u00b2", fmt::format("{:.4f}", *event.dom_heritability));
    }
    if (event.intercept != 0.0)
    {
        items.emplace_back("Intercept", fmt::format("{:.4f}", event.intercept));
    }
    items.emplace_back("Seed", fmt::format("{}", event.seed));
    for (const auto& line : gelex::header_box(title, items, 70))
    {
        logger_->info(line);
    }
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
