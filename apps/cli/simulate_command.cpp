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

#include <format>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "cli_helper.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/simulate.h"
#include "gelex/exception.h"
#include "gelex/logger.h"
#include "utils/formatter.h"

namespace
{
auto parse_effect_classes(
    argparse::ArgumentParser& sim,
    std::string_view var_flag,
    std::string_view prop_flag) -> std::vector<gelex::EffectSizeClass>
{
    auto variances = sim.get<std::vector<double>>(std::string(var_flag));
    if (!sim.is_used(std::string(prop_flag)))
    {
        throw gelex::ArgumentValidationException(
            std::format(
                "{} is required when {} is specified", prop_flag, var_flag));
    }
    auto proportions = sim.get<std::vector<double>>(std::string(prop_flag));
    if (variances.size() != proportions.size())
    {
        throw gelex::ArgumentValidationException(
            std::format(
                "{} and {} must have the same number of values",
                var_flag,
                prop_flag));
    }

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
    auto logger = gelex::logging::get();
    std::filesystem::path bed
        = gelex::BedPipe::format_bed_path(sim.get("--bfile"));

    gelex::PhenotypeSimulator::Config config{
        .bed_path = bed,
        .add_heritability = sim.get<double>("--h2"),
        .dom_heritability = sim.get<double>("--d2"),
    };

    if (sim.is_used("--add-var"))
    {
        config.add_effect_classes
            = parse_effect_classes(sim, "--add-var", "--add-prop");
    }

    if (sim.is_used("--dom-var"))
    {
        config.dom_effect_classes
            = parse_effect_classes(sim, "--dom-var", "--dom-prop");
    }

    config.intercept = sim.get<double>("--intercept");
    config.seed = sim.get<int>("--seed");
    config.output_path = sim.get("--out");

    gelex::cli::print_simulate_header(config.dom_heritability > 0.0);

    gelex::PhenotypeSimulator simulator(config);
    simulator.simulate();
    logger->info(gelex::separator());
    return 0;
}
