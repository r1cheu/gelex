#include "grm_command.h"

#include <atomic>
#include <chrono>
#include <string>

#include <barkeep.h>
#include <fmt/format.h>
#include <omp.h>

#include "cli_helper.h"
#include "data/grm_bin_writer.h"
#include "data/grm_id_writer.h"
#include "data/loader/bim_loader.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/grm.h"
#include "gelex/data/grm_code_policy.h"
#include "gelex/logger.h"
#include "utils/formatter.h"

namespace bk = barkeep;

auto grm_command(argparse::ArgumentParser& cmd) -> void
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
    cmd.add_argument("--loco").help("Compute GRM for each chromosome").flag();
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
    const std::vector<std::pair<Eigen::Index, Eigen::Index>>& ranges,
    int chunk_size,
    bool additive,
    std::function<void(Eigen::Index, Eigen::Index)> progress_callback = nullptr)
    -> gelex::GrmResult
{
    return grm.compute<CodePolicy>(
        ranges, chunk_size, additive, progress_callback);
}

auto compute_grm_with_method(
    gelex::GRM& grm,
    const std::vector<std::pair<Eigen::Index, Eigen::Index>>& ranges,
    std::string_view method,
    int chunk_size,
    bool additive,
    std::function<void(Eigen::Index, Eigen::Index)> progress_callback = nullptr)
    -> gelex::GrmResult
{
    if (method == "su")
    {
        return compute_grm<gelex::grm::Su>(
            grm, ranges, chunk_size, additive, progress_callback);
    }
    if (method == "yang")
    {
        return compute_grm<gelex::grm::Yang>(
            grm, ranges, chunk_size, additive, progress_callback);
    }
    if (method == "zeng")
    {
        return compute_grm<gelex::grm::Zeng>(
            grm, ranges, chunk_size, additive, progress_callback);
    }
    if (method == "vitezica")
    {
        return compute_grm<gelex::grm::Vitezica>(
            grm, ranges, chunk_size, additive, progress_callback);
    }
    throw std::invalid_argument(
        "Unknown GRM method: " + std::string(method)
        + ". Valid options: su, yang, zeng, vitezica");
}

auto write_grm_files(
    const gelex::GrmResult& result,
    const std::vector<std::string>& sample_ids,
    const std::string& out_prefix,
    std::shared_ptr<spdlog::logger> logger) -> void
{
    std::string bin_path = out_prefix + ".bin";
    std::string id_path = out_prefix + ".id";

    gelex::detail::GrmBinWriter(bin_path).write(result.grm, result.denominator);
    logger->info(gelex::success("GRM binary written to : {}", bin_path));

    gelex::detail::GrmIdWriter(id_path).write(sample_ids);
    logger->info(gelex::success("Sample IDs written to : {}", id_path));
}

auto compute_grm_with_progress(
    gelex::GRM& grm,
    const std::vector<std::pair<Eigen::Index, Eigen::Index>>& ranges,
    std::string_view method,
    int chunk_size,
    bool additive) -> gelex::GrmResult
{
    std::atomic<std::ptrdiff_t> progress{0};
    Eigen::Index num_snps = 0;
    for (const auto& [start, end] : ranges)
    {
        num_snps += (end - start);
    }

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

    gelex::GrmResult result = compute_grm_with_method(
        grm, ranges, method, chunk_size, additive, progress_callback);

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
    bool do_loco = cmd.get<bool>("--loco");

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

    struct ChrGroup
    {
        std::string name;
        std::vector<std::pair<Eigen::Index, Eigen::Index>> ranges;
    };

    std::vector<ChrGroup> groups;

    if (do_loco)
    {
        auto bim_path = bed_path;
        bim_path.replace_extension(".bim");
        gelex::detail::BimLoader bim_loader(bim_path);
        const auto& snp_info = bim_loader.info();

        std::string current_chr;
        Eigen::Index range_start = 0;

        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(snp_info.size());
             ++i)
        {
            if (snp_info[i].chrom != current_chr)
            {
                if (!current_chr.empty())
                {
                    groups.push_back({current_chr, {{range_start, i}}});
                }
                current_chr = snp_info[i].chrom;
                range_start = i;
            }
        }
        if (!current_chr.empty())
        {
            groups.push_back(
                {current_chr,
                 {{range_start, static_cast<Eigen::Index>(snp_info.size())}}});
        }
    }
    else
    {
        groups.push_back({"all", {{0, num_snps}}});
    }

    struct GrmTask
    {
        std::string name;
        std::string label;
        bool is_additive;
    };

    std::vector<GrmTask> tasks;
    if (do_additive)
    {
        tasks.push_back({"add", "Additive", true});
    }
    if (do_dominant)
    {
        tasks.push_back({"dom", "Dominance", false});
    }

    for (const auto& group : groups)
    {
        for (const auto& task : tasks)
        {
            logger->info(gelex::task("{} (chr{}):", task.label, group.name));

            auto result = compute_grm_with_progress(
                grm, group.ranges, method, chunk_size, task.is_additive);

            std::string path;
            if (do_loco || tasks.size() > 1)
            {
                std::string suffix
                    = do_loco ? fmt::format(".chr{}", group.name) : ".grm";
                path = fmt::format("{}.{}{}", out_prefix, task.name, suffix);
            }
            else
            {
                path = out_prefix + ".grm";
            }

            write_grm_files(result, sample_ids, path, logger);
        }
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
