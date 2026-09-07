// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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
