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

#include <Eigen/Core>
#include <algorithm>
#include <string>

#include <fmt/format.h>

#include "cli/report_printer.h"
#include "gelex/infra/logging/formatter.h"
#include "gelex/infra/logging/post_event.h"
#include "version.h"

namespace cli
{

auto PostReporter::on_event(const gelex::PostBannerEvent& /*event*/) const
    -> void
{
    cli::printer().block(
        gelex::command_banner(PROJECT_VERSION, "MCMC Posterior Analysis"));
}

auto PostReporter::on_event(const gelex::PostStartEvent& event) const -> void
{
    const auto n_chains = static_cast<Eigen::Index>(event.in_prefixes.size());
    std::string input_str
        = event.in_prefixes[0]
          + (n_chains > 1 ? fmt::format(" (+{} more)", n_chains - 1) : "");

    cli::printer().block(gelex::section("[Config]"));
    cli::printer().line("  {:<12}: {}", "Chains", n_chains);
    cli::printer().line("  {:<12}: {}", "Input", input_str);
}

auto PostReporter::on_event(const gelex::DiagnosticsReadyEvent& event) const
    -> void
{
    auto section_order = [](const std::string& section) -> int
    {
        if (section == "Parameter")
        {
            return 0;
        }
        if (section == "Additive")
        {
            return 1;
        }
        if (section == "Dominance")
        {
            return 2;
        }
        if (section == "Genetic")
        {
            return 3;
        }
        if (section == "Residual")
        {
            return 4;
        }
        return 5;
    };

    auto sorted_diags = event.diags;
    std::ranges::stable_sort(
        sorted_diags,
        [&](const auto& lhs, const auto& rhs)
        { return section_order(lhs.section) < section_order(rhs.section); });

    constexpr int TABLE_WIDTH = 92;

    double lo_pct = (1.0 - event.hdpi_prob) / 2.0 * 100.0;
    double hi_pct = (1.0 + event.hdpi_prob) / 2.0 * 100.0;
    auto hpdi_label = fmt::format("[{:g}%, {:g}%]", lo_pct, hi_pct);

    auto& p = cli::printer();

    p.block(gelex::section("[MCMC Summary]"));
    p.line(
        "   {:<12}  {:>8}  {:>8}  {:>8}  {:>21}  {:>8}  {:>8}",
        "Parameter",
        "Mean",
        "Median",
        "SD",
        hpdi_label,
        "ESS",
        "R-hat");
    p.line(gelex::table_separator(TABLE_WIDTH));

    std::string current_section;
    for (const auto& d : sorted_diags)
    {
        if (d.section != current_section)
        {
            current_section = d.section;
            if (!d.section.empty() && d.section != "Parameter")
            {
                p.line(gelex::named_section(current_section, TABLE_WIDTH, 3));
            }
        }
        p.line(
            "   {:<12}  {:>8.4f}  {:>8.4f}  {:>8.4f}  {:>21}  {:>8.1f}  "
            "{:>8.3f}",
            d.name,
            d.mean,
            d.median,
            d.sd,
            fmt::format("[{:8.4f}, {:8.4f}]", d.hpdi_lo, d.hpdi_hi),
            d.ess,
            d.rhat);
    }
    p.line(gelex::table_separator(TABLE_WIDTH));
}

}  // namespace cli
