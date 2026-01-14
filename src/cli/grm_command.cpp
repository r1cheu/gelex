#include "gelex/cli/grm_command.h"

#include <string>

#include <omp.h>

#include "gelex/data/bed_pipe.h"
#include "gelex/data/grm.h"
#include "gelex/data/grm_code_policy.h"
#include "gelex/logger.h"

#include "data/grm_bin_writer.h"
#include "data/grm_id_writer.h"

void grm_command(argparse::ArgumentParser& cmd)
{
    cmd.add_description(
        "Compute genomic relationship matrix (GRM) from PLINK "
        "binary files and output in GCTA format");

    // ================================================================
    // Data Files
    // ================================================================
    cmd.add_group("Data Files");
    cmd.add_argument("-b", "--bfile")
        .help("PLINK binary file prefix (.bed/.bim/.fam)")
        .metavar("<BFILE>")
        .required();
    cmd.add_argument("-o", "--out")
        .help("Output file prefix")
        .metavar("<OUT>")
        .default_value(std::string("grm"));

    // ================================================================
    // GRM Options
    // ================================================================
    cmd.add_group("GRM Options");
    cmd.add_argument("--method")
        .help("GRM computation method: su, yang, zeng, vitezica")
        .metavar("<METHOD>")
        .default_value(std::string("su"));
    cmd.add_argument("--chunk")
        .help("Chunk size for memory-efficient computation")
        .metavar("<SIZE>")
        .default_value(10000)
        .scan<'i', int>();
    cmd.add_argument("--threads")
        .help("Number of threads (-1 for all cores)")
        .metavar("<N>")
        .default_value(-1)
        .scan<'i', int>();
    cmd.add_argument("--additive").help("Compute additive GRM").flag();
    cmd.add_argument("--dominant").help("Compute dominance GRM").flag();
}

namespace
{

template <typename CodePolicy>
auto compute_grm(gelex::GRM& grm, int chunk_size, bool additive)
    -> Eigen::MatrixXd
{
    return grm.compute<CodePolicy>(chunk_size, additive);
}

auto compute_grm_with_method(
    gelex::GRM& grm,
    std::string_view method,
    int chunk_size,
    bool additive) -> Eigen::MatrixXd
{
    if (method == "su")
    {
        return compute_grm<gelex::grm::Su>(grm, chunk_size, additive);
    }
    if (method == "yang")
    {
        return compute_grm<gelex::grm::Yang>(grm, chunk_size, additive);
    }
    if (method == "zeng")
    {
        return compute_grm<gelex::grm::Zeng>(grm, chunk_size, additive);
    }
    if (method == "vitezica")
    {
        return compute_grm<gelex::grm::Vitezica>(grm, chunk_size, additive);
    }
    throw std::invalid_argument(
        "Unknown GRM method: " + std::string(method)
        + ". Valid options: su, yang, zeng, vitezica");
}

auto write_grm_files(
    const Eigen::MatrixXd& G,
    const std::vector<std::string>& sample_ids,
    const std::string& out_prefix,
    std::shared_ptr<spdlog::logger> logger) -> void
{
    std::string bin_path = out_prefix + ".grm.bin";
    std::string id_path = out_prefix + ".grm.id";

    gelex::detail::GrmBinWriter(bin_path).write(G);
    logger->info("GRM binary written to: {}", bin_path);

    gelex::detail::GrmIdWriter(id_path).write(sample_ids);
    logger->info("Sample IDs written to: {}", id_path);
}

}  // namespace

auto grm_execute(argparse::ArgumentParser& cmd) -> int
{
    auto logger = gelex::logging::get();

    // Set thread count
    int threads = cmd.get<int>("--threads");
    if (threads > 0)
    {
        omp_set_num_threads(threads);
    }
    logger->info("Using {} threads", omp_get_max_threads());

    // Parse arguments
    std::filesystem::path bed_path
        = gelex::BedPipe::format_bed_path(cmd.get("--bfile"));
    std::string out_prefix = cmd.get("--out");
    std::string method = cmd.get("--method");
    int chunk_size = cmd.get<int>("--chunk");
    bool do_additive = cmd.get<bool>("--additive");
    bool do_dominant = cmd.get<bool>("--dominant");

    // Default to additive if neither specified
    if (!do_additive && !do_dominant)
    {
        do_additive = true;
    }

    logger->info("Input: {}", bed_path.string());
    logger->info("Output prefix: {}", out_prefix);
    logger->info("Method: {}", method);
    logger->info("Chunk size: {}", chunk_size);

    // Create GRM object
    gelex::GRM grm(bed_path);
    const auto& sample_ids = grm.sample_ids();
    logger->info("Number of samples: {}", sample_ids.size());

    // Compute and write GRM
    if (do_additive && do_dominant)
    {
        logger->info("Computing additive GRM...");
        Eigen::MatrixXd G_add
            = compute_grm_with_method(grm, method, chunk_size, true);
        write_grm_files(G_add, sample_ids, out_prefix + ".add", logger);

        logger->info("Computing dominance GRM...");
        Eigen::MatrixXd G_dom
            = compute_grm_with_method(grm, method, chunk_size, false);
        write_grm_files(G_dom, sample_ids, out_prefix + ".dom", logger);
    }
    else if (do_additive)
    {
        logger->info("Computing additive GRM...");
        Eigen::MatrixXd G
            = compute_grm_with_method(grm, method, chunk_size, true);
        write_grm_files(G, sample_ids, out_prefix, logger);
    }
    else
    {
        logger->info("Computing dominance GRM...");
        Eigen::MatrixXd G
            = compute_grm_with_method(grm, method, chunk_size, false);
        write_grm_files(G, sample_ids, out_prefix, logger);
    }

    logger->info("GRM computation completed successfully");
    return 0;
}
