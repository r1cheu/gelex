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

#include <catch2/catch_test_macros.hpp>
#include <memory>
#include <spdlog/logger.h>
#include <spdlog/sinks/ostream_sink.h>
#include <sstream>
#include <string>
#include <utility>

#include "cli/logging.h"
#include "cli/summary.h"

namespace
{

auto capture_output(const cli::Summary& summary) -> std::string
{
    std::ostringstream output;
    auto sink = std::make_shared<spdlog::sinks::ostream_sink_mt>(output, true);
    sink->set_pattern("%v");
    auto logger = std::make_shared<spdlog::logger>("summary_test", sink);
    auto previous_logger = std::exchange(cli::logging::get(), logger);

    summary.show();
    logger->flush();

    cli::logging::get() = std::move(previous_logger);
    return output.str();
}

}  // namespace

TEST_CASE("Summary shows ordered fields", "[cli][summary]")
{
    cli::Summary summary{"Dataset Summary"};
    summary.field("Trait", "{}", "Yield")
        .field("Samples", "{}", 1000)
        .field("SNPs", "{} markers", 50000);

    REQUIRE(
        capture_output(summary)
        == "\n"
           " Dataset Summary:\n"
           "   Trait           : Yield\n"
           "   Samples         : 1000\n"
           "   SNPs            : 50000 markers\n");
}

TEST_CASE("Summary can show only a heading", "[cli][summary]")
{
    const cli::Summary summary{"Dataset Summary"};

    REQUIRE(capture_output(summary) == "\n Dataset Summary:\n");
}
