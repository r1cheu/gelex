// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "diagnostics_reporter.h"

#include <algorithm>
#include <cmath>
#include <fmt/format.h>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/bayes/genetic/diagnostics_traits.h"
#include "gelex/bayes/serialization_ids.h"

#include "cli/formatter.h"
#include "cli/report_printer.h"
#include "cli/table.h"
#include "cli/theme.h"

namespace
{

constexpr double rhat_threshold = 1.05;

auto is_model_level(const gelex::DiagnosticEntry& entry) -> bool
{
    const auto suffix = fmt::format("/{}", gelex::coefficients_id);
    return !entry.id.ends_with(suffix);
}

auto multi_row_ids(std::span<const gelex::DiagnosticEntry> entries)
    -> std::vector<std::string_view>
{
    std::vector<std::string_view> ids;
    for (const auto& entry : entries)
    {
        if (entry.index > 0 && std::ranges::find(ids, entry.id) == ids.end())
        {
            ids.push_back(entry.id);
        }
    }
    return ids;
}

auto parameter_name(
    const gelex::DiagnosticEntry& entry,
    std::span<const std::string_view> multi_row) -> std::string
{
    if (std::ranges::find(multi_row, entry.id) == multi_row.end())
    {
        return entry.id;
    }
    return fmt::format("{}[{}]", entry.id, entry.index);
}

auto number(double value, std::string_view spec) -> std::string
{
    if (std::isnan(value))
    {
        return "-";
    }
    return fmt::vformat(spec, fmt::make_format_args(value));
}

auto flagged(std::string text, bool warn) -> std::string
{
    return warn ? cli::colorize(cli::ColorRole::warning, std::move(text))
                : std::move(text);
}

}  // namespace

namespace cli
{

auto show_diagnostics(std::span<const gelex::DiagnosticEntry> entries) -> void
{
    const auto multi_row = multi_row_ids(entries);

    Table table;
    table.column("Parameter", Align::left);
    table.column("Mean", Align::right);
    table.column("SD", Align::right);
    table.column("ESS", Align::right);
    table.column("R-hat", Align::right);

    for (const auto& entry : entries)
    {
        if (!is_model_level(entry))
        {
            continue;
        }
        const auto name = parameter_name(entry, multi_row);
        const auto& s = entry.stats;
        table.row(
            {name,
             number(s.mean, "{:.4f}"),
             number(s.sd, "{:.4f}"),
             number(s.ess, "{:.1f}"),
             flagged(
                 number(s.split_rhat, "{:.3f}"),
                 s.split_rhat > rhat_threshold)});
    }

    auto& p = cli::printer();
    p.block(cli::section("MCMC Summary:"));
    p.line(table.render());
}

}  // namespace cli
