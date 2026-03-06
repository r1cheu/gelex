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

#include "simulate_command.h"

#include <argparse.h>

#include "gelex/pipeline/phenotype_simulation_engine.h"
#include "simulate_config.h"
#include "simulator_reporter.h"

auto simulate_execute(argparse::ArgumentParser& sim) -> int
{
    auto config = gelex::cli::make_simulate_config(sim);
    gelex::cli::SimulatorReporter reporter;

    reporter.on_event(
        gelex::SimulateConfigLoadedEvent{
            .intercept = config.intercept,
            .add_heritability = config.add_heritability,
            .dom_heritability = config.dom_heritability,
            .seed = config.seed,
        });

    gelex::PhenotypeSimulationEngine simulator(std::move(config));
    simulator.run(reporter.as_observer());
    return 0;
}
