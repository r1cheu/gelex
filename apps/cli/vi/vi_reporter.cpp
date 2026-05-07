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

#include "vi_reporter.h"

#include <string>

#include <fmt/format.h>

#include "config.h"
#include "gelex/algo/infer/vi/result.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/formatter.h"
#include "gelex/model/bayes/method.h"

namespace gelex::cli
{

namespace
{
const int kTableWidth = 40;
}  // namespace

ViReporter::ViReporter() : FitReporter() {}

auto ViReporter::on_event(const VIBannerEvent& /*event*/) const -> void
{
    logger_->info(
        gelex::command_banner(PROJECT_VERSION, "Model Fitting (CAVI)"));
    logger_->info("");
}

auto ViReporter::on_event(const VIConfigEvent& event) const -> void
{
    logger_->info(gelex::section("[Config]"));
    logger_->info("  {:<12}: {}", "Method", fmt::format("{}", event.method));
    logger_->info("  {:<12}: {}", "Max iters", event.max_iters);
    logger_->info("  {:<12}: {:.1e}", "Tolerance", event.tol);
    logger_->info("");
}

auto ViReporter::on_event(const VIProgressEvent& event) -> void
{
    if (!init_progress_)
    {
        init_progress_ = true;
        logger_->info("");
        logger_->info(gelex::section("[CAVI Optimization]"));
        cavi_info_ = create_progress_info();
        cavi_info_.display->show();
    }

    if (event.done)
    {
        cavi_info_.display->done();
        logger_->info("");
        return;
    }

    cavi_info_.progress_info->message(
        fmt::format(
            "iter {:>4} | ELBO: {:.2f} | Δ: {:.2e}",
            event.current,
            event.elbo,
            event.delta));
}

auto ViReporter::on_event(const VICompleteEvent& event) const -> void
{
    if (event.result == nullptr)
    {
        return;
    }

    const auto& result = *event.result;
    logger_->info(gelex::section("[Variational Posterior Summary]"));
    logger_->info("");
    logger_->info("  {:<8} {:>8} {:>8}", "Parameter", "Mean", "SD");
    logger_->info(gelex::table_separator(kTableWidth));

    const auto& fixed = result.fixed();
    fixed.for_each_term([&](const std::string& term, Eigen::Index i)
                        { print_summary_row(term, fixed.coeffs, i); });

    logger_->info(gelex::table_separator(kTableWidth));
    logger_->info("");
}

}  // namespace gelex::cli
