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

#include <cstddef>

#include "gelex/data/grm/progress.h"

#include "cli/formatter.h"
#include "cli/report_printer.h"

namespace cli
{

GrmReporter::GrmReporter(std::size_t total)
    : progress_{"", total, "SNP"},
      estimate_rate_{cli::make_rate()},
      estimate_eta_{cli::make_eta(total)}
{
}

auto GrmReporter::show_data_loaded(size_t num_samples, size_t num_snps) -> void
{
    cli::printer().block(cli::section("Dataset Summary:"));
    cli::printer().line(
        cli::field("Analyzed Samples", "{} samples", num_samples));
    cli::printer().line(cli::field("SNPs", "{} markers", num_snps));
}

auto GrmReporter::on_event(const gelex::GrmProgressEvent& event) -> void
{
    if (event.done)
    {
        progress_.finish();
        cli::printer().on_progress_finished();
        return;
    }

    progress_.update(
        {.current = event.current,
         .rate = estimate_rate_(event.current),
         .eta = estimate_eta_(event.current)});
}

}  // namespace cli
