// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "progress.h"

#include <cstddef>
#include <fmt/color.h>
#include <fmt/format.h>
#include <string>
#include <string_view>

#include "cli/report_printer.h"
#include "cli/theme.h"

namespace cli
{
namespace
{

auto phase_prefix(std::string_view phase, ColorRole color) -> std::string
{
    return fmt::format(
        "[{}]",
        colorize(style_for(color) | fmt::emphasis::bold, std::string{phase}));
}

}  // namespace

GenotypeProgress::GenotypeProgress(std::size_t total)
    : progress_{"", total, "SNP"},
      estimate_rate_{cli::make_rate()},
      estimate_eta_{cli::make_eta(total)}
{
}

auto GenotypeProgress::operator()(std::size_t current) -> void
{
    progress_.update(
        {.current = current,
         .rate = estimate_rate_(current),
         .eta = estimate_eta_(current)});
}

auto GenotypeProgress::finish() -> void
{
    progress_.finish();
    cli::printer().on_progress_finished();
}

McmcProgress::McmcProgress(std::size_t total, int burn_in)
    : burn_in_{static_cast<std::size_t>(burn_in)},
      progress_{
          burn_in == 0 ? phase_prefix("SAMPLE", ColorRole::success)
                       : phase_prefix("BURN-IN", ColorRole::warning),
          total,
          "iter"},
      estimate_rate_{cli::make_rate()},
      estimate_eta_{cli::make_eta(total)}
{
}

auto McmcProgress::operator()(std::size_t current) -> void
{
    progress_.update(
        {.current = current,
         .rate = estimate_rate_(current),
         .eta = estimate_eta_(current)});
    if (current == burn_in_)
    {
        progress_.set_prefix(phase_prefix("SAMPLE", ColorRole::success));
    }
}

auto McmcProgress::finish() -> void
{
    progress_.finish();
    cli::printer().on_progress_finished();
}

}  // namespace cli
