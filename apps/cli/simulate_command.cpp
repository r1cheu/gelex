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
#include <optional>
#include <variant>
#include <vector>

#include "cli_helper.h"
#include "config_factory.h"
#include "gelex/algo/sim/effect_sampler.h"
#include "gelex/pipeline/sim/phenotype_simulation_engine.h"
#include "simulate_config.h"
#include "simulater_reporter.h"

namespace
{

auto create_effectsize_vec(
    std::span<const double> variances,
    std::span<const double> proportions) -> std::vector<gelex::EffectSizeClass>
{
    std::vector<gelex::EffectSizeClass> classes(variances.size());
    for (size_t i = 0; i < variances.size(); ++i)
    {
        classes[i] = {proportions[i], variances[i]};
    }
    return classes;
}

}  // namespace

int simulate_execute(argparse::ArgumentParser& sim)
{
    auto config = gelex::cli::make_config<SimulateConfig>(sim);
    gelex::cli::SimulaterReporter logger;

    gelex::PhenotypeSimulationEngine::Config simulator_config{
        .bed_path = config.bed_path,
        .output_path = config.output_path,

        .intercept = config.intercept,

        .add_heritability = config.add_heritability,
        .add_effect_classes = create_effectsize_vec(
            config.additive_variances, config.additive_proportions),

        .dom_heritability = config.dom_heritability,
        .dom_effect_classes
        = config.dom_heritability
              ? create_effectsize_vec(
                    config.dominance_variances, config.dominance_proportions)
              : std::vector<gelex::EffectSizeClass>{},
        .seed = config.seed,
    };

    gelex::cli::print_simulate_header(config.dom_heritability.has_value());

    logger.on_event(
        gelex::ParameterLoadedEvent{
            .intercept = config.intercept,
            .add_heritability = config.add_heritability,
            .dom_heritability = config.dom_heritability,
            .seed = config.seed,
        });

    gelex::PhenotypeSimulationEngine simulator(simulator_config);
    simulator.simulate(
        [&logger](const gelex::SimulateEvent& event)
        {
            std::visit(
                [&logger](const auto& concrete_event)
                { logger.on_event(concrete_event); },
                event);
        });
    return 0;
}
