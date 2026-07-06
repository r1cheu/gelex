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

#include "geno_reporter.h"

#include <stdio.h>
#include <unistd.h>
#include <cstdint>

#include <fmt/format.h>
#include <string>

#include "cli/formatter.h"
#include "cli/progress_bar.h"
#include "cli/report_printer.h"
#include "gelex/types/genetic_mode.h"

namespace cli
{

GenoReporter::GenoReporter() : progress_info_(cli::create_progress_info()) {}

auto GenoReporter::show_loaded(
    gelex::GeneticMode mode,
    int64_t num_snps,
    int64_t invalid_snps) const -> void
{
    const auto effective_snps = num_snps - invalid_snps;
    const std::string label
        = (mode == gelex::GeneticMode::D) ? "Dominance" : "Additive";
    const std::string msg = fmt::format(
        "   {:<13}: {} SNPs ({} invalid excluded)",
        label,
        gelex::AbbrNumber(effective_snps),
        gelex::AbbrNumber(invalid_snps));

    if (isatty(fileno(stdout)) != 0)
    {
        cli::printer().line("{}", "\033[A\r" + msg + "\033[K");
    }
    else
    {
        cli::printer().line("{}", msg);
    }
}

auto GenoReporter::on_event(const gelex::GenotypeProgressEvent& event) -> void
{
    if (!init_progress_)
    {
        init_progress_ = true;
        progress_info_ = cli::create_progress_info();
        progress_info_.display->show();
    }

    progress_info_.progress_info->message(
        fmt::format(
            "  {}/{} SNPs",
            gelex::AbbrNumber(event.current),
            gelex::AbbrNumber(event.total)));

    if (event.done)
    {
        progress_info_.display->done();
        cli::printer().on_progress_finished();
        init_progress_ = false;
    }
}

}  // namespace cli
