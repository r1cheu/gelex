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

#include <optional>
#include <utility>

#include <CLI/CLI.hpp>

#include "config.h"
#include "gelex/engine/simulation.h"
#include "gelex/infra/logging/simulate_event.h"
#include "reporter.h"

auto simulate_execute(CLI::App& sim) -> int
{
    auto config = cli::make_simulate_config(sim);
    cli::SimulatorReporter reporter;

    reporter.on_event(gelex::SimulateBannerEvent{});
    reporter.on_event(
        gelex::SimulateConfigLoadedEvent{
            .add_heritability
            = config.additive
                  ? std::make_optional(config.additive->heritability)
                  : std::nullopt,
            .dom_heritability
            = config.dominance
                  ? std::make_optional(config.dominance->heritability)
                  : std::nullopt,
            .seed = config.seed,
        });

    gelex::SimulationEngine simulator(std::move(config));
    simulator.run(reporter.as_observer());
    return 0;
}
