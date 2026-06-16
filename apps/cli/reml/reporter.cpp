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

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <CLI/CLI.hpp>
#include <Eigen/Core>

#include "cli/report_printer.h"
#include "gelex/algo/reml/statistics.h"
#include "gelex/freq/model.h"
#include "gelex/infra/logging/formatter.h"
#include "gelex/infra/logging/reml_event.h"
#include "version.h"

namespace cli
{

auto RemlCommandReporter::show_banner() const -> void
{
    cli::printer().block(
        gelex::command_banner(
            PROJECT_VERSION, "REML Variance Component Estimation"));
}

auto RemlCommandReporter::show_config(const CLI::App& cmd) const -> void
{
    auto& p = cli::printer();
    p.block(gelex::section("Options in effect:"));

    std::vector<const CLI::Option*> printed;
    for (const auto* option : cmd.parse_order())
    {
        if (option->count() == 0)
        {
            continue;
        }
        if (std::ranges::find(printed, option) != printed.end())
        {
            continue;
        }
        printed.push_back(option);

        std::string name;
        const auto& long_names = option->get_lnames();
        const auto& short_names = option->get_snames();
        if (!long_names.empty())
        {
            name = "--" + long_names.front();
        }
        else if (!short_names.empty())
        {
            name = "-" + short_names.front();
        }
        if (name.empty())
        {
            continue;
        }

        std::vector<std::string> values;
        for (const auto& value : option->results())
        {
            if (!value.empty() && value != "true")
            {
                values.push_back(value);
            }
        }
        if (values.empty())
        {
            p.line("  {}", gelex::rebecca_purple(name));
        }
        else
        {
            p.line(
                "  {} {}", gelex::rebecca_purple(name), fmt::join(values, " "));
        }
    }
}

auto RemlCommandReporter::on_event(const gelex::RemlEmInitEvent& e) -> void
{
    header_printed_ = false;
    cli::printer().line("   Initializing (EM)...");
    cli::printer().line(
        "    LogL: {:.2f} | Init variance: [{}]",
        e.loglike,
        gelex::rebecca_purple(
            fmt::format("{:.2f}", fmt::join(e.init_variances, ", "))));
}

auto RemlCommandReporter::on_event(const gelex::RemlIterationEvent& e) -> void
{
    if (!header_printed_)
    {
        std::string var_header;
        for (const auto& label : e.labels)
        {
            var_header += fmt::format("{:>12}", label);
        }
        cli::printer().block("  {:<4} {:>12} {}", "Iter", "LogL", var_header);
        cli::printer().line(gelex::table_separator(55));
        header_printed_ = true;
    }

    std::string var_str;
    for (const auto& v : e.variances)
    {
        var_str += fmt::format("{:>12.2f}", v);
    }
    cli::printer().line("  {:<4} {:>12.2f}{}", e.iter, e.loglike, var_str);
}

auto RemlCommandReporter::show_result(
    const gelex::FreqModel& model,
    const gelex::FreqState& state,
    bool converged,
    size_t iter_count,
    size_t max_iter,
    double loglike) const -> void
{
    auto& p = cli::printer();

    p.line(gelex::table_separator(55));
    p.block(gelex::named_section("REML Results", 70));

    if (converged)
    {
        p.line(
            gelex::success(
                "Converged successfully in {} iterations", iter_count));
    }
    else
    {
        p.warn("  ! REML did not converge ({} iterations)", max_iter);
        p.warn(
            "    Try to increase max_iter or check the model specification.");
    }

    p.block("  Model Fit:");
    p.line("  - AIC : {:.2f}", gelex::reml::compute_aic(model, loglike));
    p.line("  - BIC : {:.2f}", gelex::reml::compute_bic(model, loglike));

    p.block("  Variance Components:");
    p.line(
        "  {:12} {:>12} {:>12} {:>15} {:>12}",
        "Component",
        "Estimate",
        "SE",
        "Ratio",
        "SE");
    p.line(gelex::table_separator(69));

    for (size_t i = 0; i < state.random().size(); ++i)
    {
        const auto& r = state.random()[i];
        p.line(
            "  {:12} {:>12.3f} {:>12.3f} {:>15.3f} {:>12.3f}",
            model.random()[i].name,
            r.variance,
            r.variance_se,
            r.variance_ratio,
            r.variance_ratio_se);
    }

    p.line(
        "  {:12} {:>12.3f} {:>12.3f} {:>15} {:>12}",
        "Residual",
        state.residual().variance,
        state.residual().variance_se,
        "-",
        "-");

    p.line(gelex::separator(70));
}

}  // namespace cli
