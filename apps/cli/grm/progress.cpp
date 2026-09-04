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

#include "cli/report_printer.h"

namespace cli
{

GrmProgress::GrmProgress(std::size_t total)
    : progress_{"", total, "SNP"},
      estimate_rate_{cli::make_rate()},
      estimate_eta_{cli::make_eta(total)}
{
}

auto GrmProgress::operator()(std::size_t current) -> void
{
    progress_.update(
        {.current = current,
         .rate = estimate_rate_(current),
         .eta = estimate_eta_(current)});
}

auto GrmProgress::finish() -> void
{
    progress_.finish();
    cli::printer().on_progress_finished();
}

}  // namespace cli
