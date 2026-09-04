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

#include "progress.h"

#include <cstddef>
#include <fmt/color.h>
#include <fmt/format.h>
#include <string>

#include "cli/report_printer.h"
#include "cli/theme.h"

namespace cli
{

AssocProgress::AssocProgress(std::size_t total_snps)
    : progress_{"", total_snps, "SNP"},
      estimate_rate_{cli::make_rate()},
      estimate_eta_{cli::make_eta(total_snps)}
{
}

auto AssocProgress::operator()(std::size_t current) -> void
{
    progress_.update(
        {.current = current,
         .rate = estimate_rate_(current),
         .eta = estimate_eta_(current)});
}

auto AssocProgress::start_reml() -> void
{
    progress_.set_prefix(
        fmt::format(
            "[{}]",
            colorize(
                style_for(ColorRole::warning) | fmt::emphasis::bold,
                std::string{"REML"})));
}

auto AssocProgress::start_scan() -> void
{
    progress_.set_prefix(
        fmt::format(
            "[{}]",
            colorize(
                style_for(ColorRole::success) | fmt::emphasis::bold,
                std::string{"SCAN"})));
}

auto AssocProgress::finish() -> void
{
    progress_.finish();
    cli::printer().on_progress_finished();
}

}  // namespace cli
