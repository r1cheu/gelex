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
