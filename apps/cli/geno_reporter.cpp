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

#include <unistd.h>

#include <fmt/format.h>
#include <string>

#include "cli/report_printer.h"
#include "gelex/infra/logging/formatter.h"
#include "gelex/infra/logging/progress_bar.h"

namespace gelex::cli
{

GenoReporter::GenoReporter() : progress_info_(create_progress_info()) {}

auto GenoReporter::on_event(const GenotypeLoadedEvent& event) const -> void
{
    const auto effective_snps = event.num_snps - event.monomorphic_snps;
    const std::string label
        = (event.mode == GeneticMode::D) ? "Dominance" : "Additive";
    const std::string msg = fmt::format(
        "   {:<13}: {} SNPs ({} monomorphic excluded)",
        label,
        gelex::AbbrNumber(effective_snps),
        gelex::AbbrNumber(event.monomorphic_snps));

    if (isatty(fileno(stdout)) != 0)
    {
        cli::printer().line("{}", "\033[A\r" + msg + "\033[K");
    }
    else
    {
        cli::printer().line("{}", msg);
    }
}

auto GenoReporter::on_event(const GenotypeProgressEvent& event) -> void
{
    if (!init_progress_)
    {
        init_progress_ = true;
        progress_info_ = create_progress_info();
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

}  // namespace gelex::cli
