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

#include "post_reporter.h"

#include <fmt/format.h>

#include "config.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/post_event.h"
#include "gelex/infra/utils/formatter.h"

namespace gelex::cli
{

PostReporter::PostReporter() : logger_(gelex::logging::get()) {}

auto PostReporter::on_event(const PostStartEvent& event) const -> void
{
    const auto n_chains = static_cast<Eigen::Index>(event.in_prefixes.size());
    std::string input_str
        = event.in_prefixes[0]
          + (n_chains > 1 ? fmt::format(" (+{} more)", n_chains - 1) : "");

    std::string title
        = fmt::format("gelex v{} :: MCMC Posterior Analysis", PROJECT_VERSION);
    std::vector<std::pair<std::string, std::string>> items
        = {{"Chains", fmt::format("{}", n_chains)}, {"Input", input_str}};
    logger_->info(gelex::header_box(title, items, 70));
    logger_->info("");
}

auto PostReporter::on_event(const DiagnosticsReadyEvent& event) const -> void
{
    logger_->info(gelex::named_section("MCMC Summary", 80));
    logger_->info("");

    logger_->info(
        "   {:<8}  {:>8}  {:>8}  {:>21}  {:>8}  {:>8}",
        "Parameter",
        "Mean",
        "SD",
        "90% HPDI",
        "ESS",
        "R-hat");
    logger_->info(gelex::table_separator(78));

    for (const auto& d : event.diags)
    {
        logger_->info(
            "   {:<8}  {:>8.4f}  {:>8.4f}  {:>21}  {:>8.1f}  "
            "{:>8.3f}",
            d.name,
            d.mean,
            d.sd,
            fmt::format("[{:8.4f}, {:8.4f}]", d.hpdi_lo, d.hpdi_hi),
            d.ess,
            d.rhat);
    }

    logger_->info(gelex::separator(80));
}

}  // namespace gelex::cli
