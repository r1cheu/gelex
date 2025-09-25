#include "gelex/cli/simulation_command.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/simulation.h"
#include "gelex/error.h"
#include "gelex/logger.h"

void simulate_command(argparse::ArgumentParser& cmd)
{
    cmd.add_description(
        "Simulate phenotypes based on genetic data and specified parameters");

    cmd.add_argument("--bfile").help("PLINK BED file path").required();

    cmd.add_argument("--causal")
        .help(
            "File containing causal variants (one per line, optional effects)")
        .required();

    cmd.add_argument("--h2")
        .help("Heritability (0 < h2 < 1)")
        .default_value(0.5)
        .scan<'g', double>();

    cmd.add_argument("--seed")
        .help("Random seed (-1 for random)")
        .default_value(-1)
        .scan<'i', int>();

    cmd.add_argument("--out")
        .help("Output phenotype file path")
        .default_value("sim.phen");
}

int simulate_execute(argparse::ArgumentParser& sim)
{
    std::string out_prefix = sim.get("--out");

    gelex::logging::initialize(out_prefix);
    auto logger = gelex::logging::get();
    auto bed = gelex::valid_bed(sim.get("--bfile"));
    if (!bed)
    {
        logger->error(bed.error().message);
        return 1;
    }

    gelex::PhenotypeSimulator::Config config{
        .bed_path = bed.value(),
        .causal_variants_path = sim.get("--causal"),
        .heritability = sim.get<double>("--h2"),
        .seed = sim.get<int>("--seed"),
        .output_path = sim.get("--out")};

    auto simulator_result = gelex::PhenotypeSimulator::create(config);

    if (!simulator_result)
    {
        logger->error(simulator_result.error().message);
        return 1;
    }

    auto& simulator = *simulator_result;

    if (auto result = simulator.simulate(); !result)
    {
        logger->error(result.error().message);
        return 1;
    }

    logger->info("Phenotype simulation completed successfully");
    return 0;
}
