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

#include "post_command.h"

#include <argparse.h>

#include <variant>

#include "gelex/infra/logger.h"
#include "gelex/infra/logging/post_event.h"
#include "gelex/pipeline/posterior_analysis_engine.h"
#include "post_config.h"
#include "post_reporter.h"

auto post_execute(argparse::ArgumentParser& post) -> int
{
    auto config = gelex::cli::make_post_config(post);
    auto logger = gelex::logging::get();

    gelex::cli::PostReporter reporter;
    reporter.on_event(gelex::PostStartEvent{.in_prefixes = config.in_prefixes});

    gelex::PosteriorAnalysisEngine engine(
        gelex::PosteriorAnalysisEngine::Config{
            .in_prefixes = config.in_prefixes});

    engine.run(
        [&reporter](const gelex::PostEvent& event)
        {
            std::visit(
                [&reporter](const auto& ev) { reporter.on_event(ev); }, event);
        });
    return 0;
}
