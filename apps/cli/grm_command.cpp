#include "grm_command.h"

#include <atomic>
#include <chrono>
#include <string>

#include <barkeep.h>
#include <fmt/format.h>
#include <omp.h>

#include "cli_helper.h"
#include "config.h"
#include "data/grm_bin_writer.h"
#include "data/grm_id_writer.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/grm.h"
#include "gelex/data/grm_code_policy.h"
#include "gelex/logger.h"
#include "utils/formatter.h"

namespace bk = barkeep;

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

const barkeep::BarParts BAR_STYLE{
    .left = "[",
    .right = "]",
    .fill = {"\033[1;36m━\033[0m"},
    .empty = {"-"}};

template <typename CodePolicy>
auto compute_grm(
    gelex::GRM& grm,
    int chunk_size,
    bool additive,
    std::function<void(Eigen::Index, Eigen::Index)> progress_callback = nullptr)
    -> Eigen::MatrixXd
{
    return grm.compute<CodePolicy>(chunk_size, additive, progress_callback);
}

auto compute_grm_with_method(
    gelex::GRM& grm,
    std::string_view method,
    int chunk_size,
    bool additive,
    std::function<void(Eigen::Index, Eigen::Index)> progress_callback = nullptr)
    -> Eigen::MatrixXd
{
    if (method == "su")
    {
        return compute_grm<gelex::grm::Su>(
            grm, chunk_size, additive, progress_callback);
    }
    if (method == "yang")
    {
        return compute_grm<gelex::grm::Yang>(
            grm, chunk_size, additive, progress_callback);
    }
    if (method == "zeng")
    {
        return compute_grm<gelex::grm::Zeng>(
            grm, chunk_size, additive, progress_callback);
    }
    if (method == "vitezica")
    {
        return compute_grm<gelex::grm::Vitezica>(
            grm, chunk_size, additive, progress_callback);
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
    logger->info(gelex::success("GRM binary written to : {}", bin_path));

    gelex::detail::GrmIdWriter(id_path).write(sample_ids);
    logger->info(gelex::success("Sample IDs written to : {}", id_path));
}

auto compute_grm_with_progress(
    gelex::GRM& grm,
    std::string_view method,
    int chunk_size,
    bool additive,
    Eigen::Index num_snps) -> Eigen::MatrixXd
{
    std::atomic<std::ptrdiff_t> progress{0};

    auto pbar = bk::ProgressBar(
        &progress,
        {.total = num_snps,
         .format = "  {bar} {value}/{total} SNPs [{speed:.1f} snp/s]",
         .speed = 0.1,
         .style = BAR_STYLE,
         .show = false});

    pbar->show();

    auto progress_callback = [&progress](Eigen::Index current, Eigen::Index)
    { progress.store(current, std::memory_order_relaxed); };

    Eigen::MatrixXd result = compute_grm_with_method(
        grm, method, chunk_size, additive, progress_callback);

    pbar->done();

    return result;
}

}  // namespace

auto grm_execute(argparse::ArgumentParser& cmd) -> int
{
    auto logger = gelex::logging::get();

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

    // Set thread count
    int threads = cmd.get<int>("--threads");
    if (threads > 0)
    {
        omp_set_num_threads(threads);
    }
    int actual_threads = omp_get_max_threads();

    // Print header banner
    gelex::cli::print_grm_header(
        method, do_additive, do_dominant, chunk_size, actual_threads);

    // Loading Data section
    logger->info("");
    logger->info(gelex::section("Loading Data..."));
    logger->info(gelex::success("Input      : {}", bed_path.string()));

    gelex::GRM grm(bed_path);
    const auto& sample_ids = grm.sample_ids();
    Eigen::Index num_snps = grm.num_snps();
    logger->info(gelex::success("Samples    : {} samples", sample_ids.size()));
    logger->info(gelex::success("SNPs       : {} markers", num_snps));

    // Computing GRM section
    logger->info("");
    logger->info(gelex::section("Computing GRM..."));

    auto start_time = std::chrono::high_resolution_clock::now();

    if (do_additive && do_dominant)
    {
        logger->info(gelex::task("Additive:"));
        Eigen::MatrixXd G_add = compute_grm_with_progress(
            grm, method, chunk_size, true, num_snps);
        write_grm_files(G_add, sample_ids, out_prefix + ".add", logger);

        logger->info(gelex::task("Dominance:"));
        Eigen::MatrixXd G_dom = compute_grm_with_progress(
            grm, method, chunk_size, false, num_snps);
        write_grm_files(G_dom, sample_ids, out_prefix + ".dom", logger);
    }
    else if (do_additive)
    {
        logger->info(gelex::task("Additive:"));
        Eigen::MatrixXd G = compute_grm_with_progress(
            grm, method, chunk_size, true, num_snps);
        write_grm_files(G, sample_ids, out_prefix, logger);
    }
    else
    {
        logger->info(gelex::task("Dominance:"));
        Eigen::MatrixXd G = compute_grm_with_progress(
            grm, method, chunk_size, false, num_snps);
        write_grm_files(G, sample_ids, out_prefix, logger);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double>(end_time - start_time);

    logger->info("");
    logger->info(gelex::success("Time elapsed: {:.2f}s", duration.count()));
    logger->info(
        fmt::format(
            fmt::emphasis::bold | fmt::fg(fmt::color::light_cyan),
            "───────────────────────────────────"
            "───────────────────────────────────"));

    return 0;
}
