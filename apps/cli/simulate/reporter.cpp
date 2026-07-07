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

#include <optional>

#include <fmt/format.h>

#include "cli/formatter.h"
#include "cli/progress_bar.h"
#include "cli/report_printer.h"

namespace cli
{

SimulatorReporter::SimulatorReporter() : info_(cli::create_progress_info()) {}

auto SimulatorReporter::show_variance_summary(
    std::optional<double> realized_h2,
    std::optional<double> realized_d2) const -> void
{
    cli::printer().block(gelex::section("Realized:"));
    if (realized_h2)
    {
        cli::printer().line("  {:<12}: {:.4f}", "h²", *realized_h2);
    }
    if (realized_d2)
    {
        cli::printer().line("  {:<12}: {:.4f}", "d²", *realized_d2);
    }
}

auto SimulatorReporter::on_event(const gelex::SimulateProgressEvent& event)
    -> void
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
        cli::printer().on_progress_finished();
    }
}

}  // namespace cli
