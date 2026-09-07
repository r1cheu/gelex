// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "summary.h"

#include <string_view>

#include "cli/formatter.h"
#include "cli/logging.h"

namespace cli
{

Summary::Summary(std::string_view title)
{
    lines_.push_back(cli::section("{}:", title));
}

auto Summary::show() const -> void
{
    auto& logger = cli::logging::get();
    logger->info("");
    for (const auto& line : lines_)
    {
        logger->info("{}", line);
    }
}

}  // namespace cli
