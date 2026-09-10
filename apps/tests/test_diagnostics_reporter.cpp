// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <catch2/catch_test_macros.hpp>
#include <limits>
#include <memory>
#include <spdlog/logger.h>
#include <spdlog/sinks/ostream_sink.h>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "gelex/bayes/genetic/diagnostics_traits.h"
#include "gelex/bayes/stats/diagnostics.h"

#include "cli/logging.h"
#include "cli/mcmc/diagnostics_reporter.h"

namespace
{

auto capture_output(const std::vector<gelex::DiagnosticEntry>& entries)
    -> std::string
{
    std::ostringstream output;
    auto sink = std::make_shared<spdlog::sinks::ostream_sink_mt>(output, true);
    sink->set_pattern("%v");
    auto logger = std::make_shared<spdlog::logger>("diagnostics_test", sink);
    auto previous_logger = std::exchange(cli::logging::get(), logger);

    cli::show_diagnostics(entries);
    logger->flush();

    cli::logging::get() = std::move(previous_logger);
    return output.str();
}

auto stats(double mean, double ess, double rhat) -> gelex::ChainDiagnostics
{
    return {
        .mean = mean,
        .sd = 0.1,
        .median = mean,
        .hpdi_lower = mean - 0.2,
        .hpdi_upper = mean + 0.2,
        .ess = ess,
        .mcse = 0.01,
        .split_rhat = rhat};
}

}  // namespace

TEST_CASE(
    "Diagnostics reporter prints model-level entries in one table",
    "[cli][mcmc][diagnostics]")
{
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<gelex::DiagnosticEntry> entries{
        {.id = "fixed/coefficients",
         .index = 0,
         .stats = stats(1.0, 50.0, 1.0)},
        {.id = "random/Group/coefficients",
         .index = 1,
         .stats = stats(1.0, 50.0, 1.0)},
        {.id = "genetic/A/variance",
         .index = 0,
         .stats = stats(0.5, 812.3, 1.002)},
        {.id = "genetic/A/probabilities",
         .index = 0,
         .stats = stats(0.25, 40.0, 1.2)},
        {.id = "genetic/A/probabilities",
         .index = 1,
         .stats = stats(0.75, nan, nan)},
        {.id = "residual/variance",
         .index = 0,
         .stats = stats(2.0, 900.0, 1.0)}};

    const auto output = capture_output(entries);

    REQUIRE(output.find("MCMC Summary:") != std::string::npos);
    REQUIRE(output.find("HPDI") == std::string::npos);
    REQUIRE(output.find("MCSE") == std::string::npos);
    REQUIRE(output.find("coefficients") == std::string::npos);
    REQUIRE(output.find("genetic/A/variance ") != std::string::npos);
    REQUIRE(output.find("genetic/A/probabilities[0]") != std::string::npos);
    REQUIRE(output.find("genetic/A/probabilities[1]") != std::string::npos);
    REQUIRE(output.find("0.5000   0.1000") != std::string::npos);
    REQUIRE(output.find("812.3") != std::string::npos);
    REQUIRE(output.find("1.002") != std::string::npos);
}
