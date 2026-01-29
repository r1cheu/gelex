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
#include "gelex/data/bed_pipe.h"
#include "gelex/data/simulate.h"

#include "gelex/logger.h"

int simulate_execute(argparse::ArgumentParser& sim)
{
    auto logger = gelex::logging::get();
    std::filesystem::path bed
        = gelex::BedPipe::format_bed_path(sim.get("--bfile"));

    gelex::PhenotypeSimulator::Config config{
        .bed_path = bed,
        .causal_variants_path = sim.get("--causal"),
        .heritability = sim.get<double>("--h2"),
        .seed = sim.get<int>("--seed"),
        .output_path = sim.get("--out")};

    gelex::PhenotypeSimulator simulator(config);
    simulator.simulate();
    logger->info("Phenotype simulation completed successfully");
    return 0;
}
