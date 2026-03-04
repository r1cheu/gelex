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

#include "grm_command.h"

#include <argparse.h>
#include <fmt/format.h>
#include <variant>

#include "cli/cli_helper.h"
#include "gelex/infra/logging/grm_event.h"
#include "gelex/pipeline/grm_engine.h"
#include "grm_config.h"
#include "grm_reporter.h"

auto grm_execute(argparse::ArgumentParser& cmd) -> int
{
    auto config = gelex::cli::make_grm_config(cmd);
    gelex::cli::GrmReporter reporter;

    const auto method_name = fmt::format("{}", config.method);

    auto threads = cmd.get<int>("--threads");
    gelex::cli::setup_parallelization(threads);

    reporter.on_event(
        gelex::GrmConfigLoadedEvent{
            .method = std::string(method_name),
            .mode = config.mode,
            .do_loco = config.do_loco,
            .chunk_size = config.chunk_size,
            .threads = threads,
        });

    gelex::GrmEngine engine(std::move(config));

    engine.compute(
        [&reporter](const gelex::GrmEvent& event)
        {
            std::visit(
                [&reporter](const auto& e) { reporter.on_event(e); }, event);
        });

    return 0;
}
