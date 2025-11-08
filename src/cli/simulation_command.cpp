#include "gelex/cli/simulation_command.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/simulation.h"
#include "gelex/error.h"
#include "gelex/logger.h"

void simulate_command(argparse::ArgumentParser& cmd)
{
    cmd.add_description(
        "Simulate phenotypes based on genetic data and specified parameters");

    // ================================================================
    // Data Files
    // ================================================================
    cmd.add_group("Data Files");
    cmd.add_argument("-b", "--bfile")
        .help("PLINK binary file prefix (.bed/.bim/.fam)")
        .metavar("<BFILE>")
        .required();
    cmd.add_argument("-c", "--causal")
        .help(
            "Causal variants file (TSV: SNP ID per line, optional effect size)")
        .metavar("<CAUSAL>")
        .required();
    cmd.add_argument("-o", "--out")
        .help("Output file prefix for simulated phenotypes")
        .metavar("<OUT>")
        .default_value("sim.phen");

    // ================================================================
    // Simulation Parameters
    // ================================================================
    cmd.add_group("Simulation Parameters");
    cmd.add_argument("--h2")
        .help("Narrow-sense heritability (range: 0-1)")
        .default_value(0.5)
        .scan<'g', double>();
    cmd.add_argument("--seed")
        .help("Random seed for reproducibility (-1 for time-based)")
        .default_value(-1)
        .scan<'i', int>();
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
