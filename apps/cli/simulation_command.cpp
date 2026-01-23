#include "simulation_command.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/simulation.h"

#include "gelex/logger.h"

int simulate_execute(argparse::ArgumentParser& sim)
{
    std::string out_prefix = sim.get("--out");

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
