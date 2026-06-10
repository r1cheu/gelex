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

#include "reml_reporter.h"

#include <cstddef>
#include <string>
#include <vector>

#include <fmt/color.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <Eigen/Core>

#include "cli/report_printer.h"
#include "gelex/algo/reml/statistics.h"
#include "gelex/freq/model.h"
#include "gelex/infra/logging/formatter.h"
#include "gelex/infra/logging/reml_event.h"

namespace cli
{

auto RemlReporter::on_event(const gelex::RemlEmInitEvent& e) -> void
{
    header_printed_ = false;
    cli::printer().line("   Initializing (EM)...");
    cli::printer().line(
        "    LogL: {:.2f} | Init variance: [{}]",
        e.loglike,
        gelex::rebecca_purple(
            fmt::format("{:.2f}", fmt::join(e.init_variances, ", "))));
}

auto RemlReporter::on_event(const gelex::RemlIterationEvent& e) -> void
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

auto RemlReporter::on_event(const gelex::RemlCompleteEvent& e) const -> void
{
    const auto& model = *e.model;
    const auto& state = *e.state;
    auto& p = cli::printer();

    p.line(gelex::table_separator(55));
    p.block(gelex::named_section("REML Results", 70));

    if (e.converged)
    {
        p.line(
            gelex::success(
                "Converged successfully in {} iterations", e.iter_count));
    }
    else
    {
        p.warn("  ! REML did not converge ({} iterations)", e.max_iter);
        p.warn(
            "    Try to increase max_iter or check the model specification.");
    }

    // model fit
    p.block("  Model Fit:");
    p.line("  - AIC : {:.2f}", gelex::reml::compute_aic(model, e.loglike));
    p.line("  - BIC : {:.2f}", gelex::reml::compute_bic(model, e.loglike));

    // fixed effects
    p.block("  Fixed Effects:");
    p.line("  {:12} {:>12} {:>12}", "Effect", "Estimate", "SE");
    p.line(gelex::table_separator(40));
    for (Eigen::Index i = 0; i < state.fixed().coeffs.size(); ++i)
    {
        std::string name = fmt::format("X{}", i);
        if (static_cast<size_t>(i) < model.fixed().names.size())
        {
            name = model.fixed().names[i];
        }
        p.line(
            "  {:12} {:>12.3f} {:>12.3f}",
            name,
            state.fixed().coeffs(i),
            state.fixed().se(i));
    }

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

void print_loco_reml_summary(const std::vector<gelex::LocoRemlResult>& results)
{
    if (results.empty())
    {
        return;
    }

    auto& p = cli::printer();

    size_t num_random = results[0].random.size();
    auto format_variances = [](const auto& values) -> std::string
    {
        std::string row;
        for (const auto& v : values)
        {
            row += fmt::format("  {:>10.4f}", v);
        }
        return row;
    };

    std::string header = fmt::format("  {:>5}  {:>10}", "Chr", "LogL");
    for (const auto& r : results[0].random)
    {
        header += fmt::format("  {:>10}", fmt::format("V({})", r.name));
    }
    header += fmt::format("  {:>10}  {:>4}", "V(e)", "Conv");

    p.block(gelex::named_section("LOCO REML Summary", 70));
    p.line("{}", header);
    p.line("{}", gelex::table_separator());

    std::vector<double> sum_random_variance(num_random, 0.0);
    std::vector<double> sum_random_ratio(num_random, 0.0);
    double sum_ve = 0.0;

    for (const auto& r : results)
    {
        std::string conv_mark
            = r.converged ? fmt::format(fmt::fg(fmt::color::light_green), "✓")
                          : fmt::format(fmt::fg(fmt::color::orange_red), "✗");

        std::vector<double> variances(num_random);
        for (size_t i = 0; i < num_random; ++i)
        {
            variances[i] = (i < r.random.size()) ? r.random[i].variance : 0.0;
            sum_random_variance[i] += variances[i];
            sum_random_ratio[i]
                += (i < r.random.size()) ? r.random[i].variance_ratio : 0.0;
        }

        p.line(
            "  {:>5}  {:>10.2f}{}  {:>10.4f}    {}",
            r.chr_name,
            r.loglike,
            format_variances(variances),
            r.residual_variance,
            conv_mark);

        sum_ve += r.residual_variance;
    }

    p.line("{}", gelex::table_separator());

    auto n = static_cast<double>(results.size());
    std::vector<double> mean_random_variance(num_random);
    std::vector<double> mean_random_ratio(num_random);
    for (size_t i = 0; i < num_random; ++i)
    {
        mean_random_variance[i] = sum_random_variance[i] / n;
        mean_random_ratio[i] = sum_random_ratio[i] / n;
    }

    p.line(
        "  {:>5}  {:>10}{}  {:>10.4f}",
        "Mean",
        "",
        format_variances(mean_random_variance),
        sum_ve / n);
    p.line(
        "  {:>5}  {:>10}{}", "Ratio", "", format_variances(mean_random_ratio));
    p.line(gelex::separator());
}

}  // namespace cli
