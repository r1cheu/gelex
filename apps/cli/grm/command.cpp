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

#include "command.h"

#include <string>
#include <utility>

#include <fmt/format.h>
#include <CLI/CLI.hpp>

#include "cli/cli_helper.h"
#include "config.h"
#include "gelex/engine/grm.h"
#include "gelex/infra/logging/grm_event.h"
#include "reporter.h"

auto grm_execute(CLI::App& cmd) -> int
{
    auto config = cli::make_grm_config(cmd);
    cli::GrmReporter reporter;

    const auto method_name = fmt::format("{}", config.method);

    auto threads = cmd.get_option("--threads")->as<int>();
    cli::setup_parallelization(threads);

    cli::GrmReporter::on_event(gelex::GrmBannerEvent{});
    cli::GrmReporter::on_event(
        gelex::GrmConfigLoadedEvent{
            .method = std::string(method_name),
            .requested_effects = config.requested_effects,
            .do_loco = config.do_loco,
        });

    gelex::GrmEngine engine(std::move(config));

    engine.compute(reporter.as_observer());

    return 0;
}
