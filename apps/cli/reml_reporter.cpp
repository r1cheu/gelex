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

#include <Eigen/Core>
#include <cstddef>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/algo/reml/statistics.h"
#include "gelex/algo/reml/summary.h"
#include "gelex/freq/model.h"
#include "gelex/infra/logging/reml_event.h"

#include "cli/formatter.h"
#include "cli/report_printer.h"
#include "theme.h"

namespace cli
{

auto RemlReporter::show_dataset_summary(
    const gelex::FreqModel& model,
    std::string_view pheno_name) const -> void
{
    auto& p = cli::printer();
    p.block(cli::section("Dataset Summary:"));
    p.line(cli::field("Trait", "{}", pheno_name));
    p.line(cli::field("Analyzed Samples", "{}", model.num_individuals()));
    p.line(cli::field("Covariates", "{}", model.fixed().X.cols()));

    auto names = model.random()
                 | std::views::transform([](const auto& d)
                                         { return std::string_view(d.name); });
    p.line(cli::field("Random effects", "{}", fmt::join(names, ", ")));
}

auto RemlReporter::on_event(const gelex::RemlIterationEvent& e) -> void
{
    if (!header_printed_)
    {
        cli::printer().block(cli::section("REML Iterations:"));
        iter_table_.column("Iter", Align::right, 4);
        iter_table_.column("LogL", Align::right, 10);
        for (const auto& label : e.labels)
        {
            iter_table_.column(label, Align::right, 10);
        }
        cli::printer().line(iter_table_.stream_header());
        header_printed_ = true;
    }

    std::vector<std::string> cells;
    cells.reserve(2 + e.variances.size());
    cells.push_back(fmt::format("{}", e.iter));
    cells.push_back(fmt::format("{:.2f}", e.loglike));
    for (const auto& v : e.variances)
    {
        cells.push_back(fmt::format("{:.2f}", v));
    }
    cli::printer().line(iter_table_.stream_row(cells));
}

auto RemlReporter::on_event(const gelex::RemlConstrainedEvent& e) -> void
{
    auto& p = cli::printer();
    if (e.num_constrained == e.num_total)
    {
        p.warn(
            "  ! All {} variance components were constrained; the estimate is "
            "not reliable.",
            e.num_total);
    }
    else
    {
        p.warn(
            "  ! {} of {} variance components were constrained; the estimate "
            "may be unreliable.",
            e.num_constrained,
            e.num_total);
    }
}

auto RemlReporter::show_result(
    const gelex::FreqModel& model,
    const gelex::RemlSummary& summary,
    size_t max_iter) const -> void
{
    auto& p = cli::printer();

    p.block(cli::section("REML Summary:"));

    if (summary.converged)
    {
        p.line(
            cli::success(
                "Converged successfully in {} iterations", summary.iter_count));
    }
    else
    {
        p.warn("  ! REML did not converge ({} iterations)", max_iter);
        p.warn(
            "    Try to increase max_iter or check the model specification.");
    }

    p.ensure_blank();

    Table t;
    t.column("Component", Align::left);
    t.column("Estimate", Align::right);
    t.column("Est.SE", Align::right);
    t.column("Ratio", Align::right);
    t.column("Ratio.SE", Align::right);
    t.column("Pval", Align::right);
    for (const auto& r : summary.random)
    {
        t.row(
            {r.name,
             fmt::format("{:.3f}", r.variance),
             fmt::format("{:.3f}", r.variance_se),
             fmt::format("{:.3f}", r.variance_ratio),
             fmt::format("{:.3f}", r.variance_ratio_se),
             fmt::format(
                 "{:.3g}",
                 gelex::wald_p_onesided(r.variance / r.variance_se))});
    }
    t.row(
        {"Residual",
         fmt::format("{:.3f}", summary.residual_variance),
         fmt::format("{:.3f}", summary.residual_variance_se),
         "-",
         "-",
         "-"});
    p.line(t.render());

    p.ensure_blank();
    p.line(
        "  logL {:.2f}      AIC {:.2f}      BIC {:.2f}",
        summary.loglike,
        gelex::compute_aic(model, summary.loglike),
        gelex::compute_bic(model, summary.loglike));
}

void print_loco_reml_summary(const std::vector<gelex::LocoRemlResult>& results)
{
    if (results.empty())
    {
        return;
    }

    auto& p = cli::printer();

    size_t num_random = results[0].summary.random.size();

    Table t;
    t.column("Chr", Align::right);
    t.column("LogL", Align::right);
    for (const auto& r : results[0].summary.random)
    {
        t.column(fmt::format("V({})", r.name), Align::right);
    }
    t.column("V(e)", Align::right);
    t.column("Conv", Align::right);

    std::vector<double> sum_random_variance(num_random, 0.0);
    std::vector<double> sum_random_ratio(num_random, 0.0);
    double sum_ve = 0.0;

    for (const auto& result : results)
    {
        const auto& s = result.summary;
        std::string conv_mark = s.converged ? colorize(ColorRole::success, "✓")
                                            : colorize(ColorRole::error, "✗");

        std::vector<std::string> cells;
        cells.reserve(num_random + 4);
        cells.push_back(result.chr_name);
        cells.push_back(fmt::format("{:.2f}", s.loglike));
        for (size_t i = 0; i < num_random; ++i)
        {
            double v = (i < s.random.size()) ? s.random[i].variance : 0.0;
            cells.push_back(fmt::format("{:.4f}", v));
            sum_random_variance[i] += v;
            sum_random_ratio[i]
                += (i < s.random.size()) ? s.random[i].variance_ratio : 0.0;
        }
        cells.push_back(fmt::format("{:.4f}", s.residual_variance));
        cells.push_back(std::move(conv_mark));
        t.row(std::move(cells));

        sum_ve += s.residual_variance;
    }

    t.rule();

    auto n = static_cast<double>(results.size());
    std::vector<std::string> mean_cells{"Mean", ""};
    std::vector<std::string> ratio_cells{"Ratio", ""};
    mean_cells.reserve(num_random + 4);
    ratio_cells.reserve(num_random + 2);
    for (size_t i = 0; i < num_random; ++i)
    {
        mean_cells.push_back(fmt::format("{:.4f}", sum_random_variance[i] / n));
        ratio_cells.push_back(fmt::format("{:.4f}", sum_random_ratio[i] / n));
    }
    mean_cells.push_back(fmt::format("{:.4f}", sum_ve / n));
    t.row(std::move(mean_cells));
    t.row(std::move(ratio_cells));

    p.block(cli::section("LOCO REML Summary:"));
    p.line(t.render());
}

}  // namespace cli
