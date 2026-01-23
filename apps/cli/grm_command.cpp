#include "grm_command.h"

#include <atomic>
#include <string>
#include <string_view>
#include <vector>

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
#include "utils/utils.h"

namespace bk = barkeep;

namespace
{

template <typename CodePolicy>
auto compute_grm_impl(
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
        return compute_grm_impl<gelex::grm::Su>(
            grm, ranges, chunk_size, additive, progress_callback);
    }
    if (method == "yang")
    {
        return compute_grm_impl<gelex::grm::Yang>(
            grm, ranges, chunk_size, additive, progress_callback);
    }
    if (method == "zeng")
    {
        return compute_grm_impl<gelex::grm::Zeng>(
            grm, ranges, chunk_size, additive, progress_callback);
    }
    if (method == "vitezica")
    {
        return compute_grm_impl<gelex::grm::Vitezica>(
            grm, ranges, chunk_size, additive, progress_callback);
    }
    throw std::invalid_argument(
        fmt::format(
            "Unknown GRM method: {}. Valid: su, yang, zeng, vitezica", method));
}

auto write_grm_files(
    const gelex::GrmResult& result,
    const std::vector<std::string>& sample_ids,
    const std::string& out_prefix) -> std::vector<std::string>
{
    auto bin_path = out_prefix + ".bin";
    auto id_path = out_prefix + ".id";

    gelex::detail::GrmBinWriter(bin_path).write(result.grm);
    gelex::detail::GrmIdWriter(id_path).write(sample_ids);

    return {bin_path, id_path};
}

}  // namespace

auto grm_execute(argparse::ArgumentParser& cmd) -> int
{
    auto logger = gelex::logging::get();

    auto bed_path = gelex::BedPipe::format_bed_path(cmd.get("--bfile"));
    auto out_prefix = cmd.get("--out");
    auto method = cmd.get("--method");
    auto chunk_size = cmd.get<int>("--chunk");
    auto do_additive = cmd.get<bool>("--additive");
    auto do_dominant = cmd.get<bool>("--dominant");
    auto do_loco = cmd.get<bool>("--loco");
    auto threads = cmd.get<int>("--threads");

    if (!do_additive && !do_dominant)
    {
        do_additive = true;
    }

    gelex::cli::setup_parallelization(threads);
    int actual_threads = omp_get_max_threads();

    gelex::cli::print_grm_header(
        method, do_additive, do_dominant, chunk_size, actual_threads);

    logger->info("");
    logger->info(gelex::section("Loading Data..."));
    logger->info(gelex::success("Input      : {}", bed_path.string()));

    gelex::GRM grm(bed_path);
    const auto& sample_ids = grm.sample_ids();
    auto num_snps = grm.num_snps();
    logger->info(gelex::success("Samples    : {} samples", sample_ids.size()));
    logger->info(gelex::success("SNPs       : {} markers", num_snps));

    logger->info("");
    logger->info(gelex::section("Computing GRM..."));

    auto bim_path = bed_path;
    bim_path.replace_extension(".bim");
    gelex::detail::BimLoader bim_loader(bim_path);

    auto groups = gelex::cli::build_chr_groups(do_loco, bim_loader.info());

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

    size_t total_work_snps = 0;
    for (const auto& group : groups)
    {
        total_work_snps += static_cast<size_t>(group.total_snps) * tasks.size();
    }

    std::atomic<size_t> global_progress{0};
    auto pbar = gelex::cli::create_progress_bar(
        global_progress, total_work_snps, "{bar}");
    pbar.display->show();

    gelex::SmoothEtaCalculator eta_calculator(total_work_snps);
    std::vector<std::string> generated_files;
    size_t completed_snps_base = 0;

    for (const auto& group : groups)
    {
        for (const auto& task : tasks)
        {
            auto progress_callback = [&](Eigen::Index current, Eigen::Index)
            {
                auto current_total
                    = completed_snps_base + static_cast<size_t>(current);
                global_progress.store(current_total, std::memory_order_relaxed);

                pbar.status->message(
                    fmt::format(
                        "{:.1f}% ({}/{}) | {}",
                        static_cast<double>(current_total) / total_work_snps
                            * 100,
                        gelex::HumanReadable(current_total),
                        gelex::HumanReadable(total_work_snps),
                        eta_calculator.get_eta(current_total)));
            };

            auto result = compute_grm_with_method(
                grm,
                group.ranges,
                method,
                chunk_size,
                task.is_additive,
                progress_callback);

            auto path = out_prefix;
            if (do_loco || tasks.size() > 1)
            {
                auto suffix = do_loco ? fmt::format(".chr{}", group.name) : "";
                path = fmt::format("{}.{}{}", out_prefix, task.name, suffix);
            }

            auto files = write_grm_files(result, sample_ids, path);
            generated_files.insert(
                generated_files.end(), files.begin(), files.end());

            completed_snps_base += static_cast<size_t>(group.total_snps);
        }
    }

    pbar.display->done();

    logger->info("");
    logger->info(
        gelex::success(
            "Time elapsed: {}", eta_calculator.total_time_consumed()));
    logger->info("--------------------------------------------------");
    logger->info("Computation Summary:");
    logger->info("  Total Files : {}", generated_files.size());
    logger->info(
        "  Output Dir  : {}",
        std::filesystem::absolute(std::filesystem::path(out_prefix))
            .parent_path()
            .string());
    if (!generated_files.empty())
    {
        logger->info("  Pattern     : {}...", generated_files[0]);
    }
    logger->info(
        fmt::format(
            fmt::emphasis::bold | fmt::fg(fmt::color::light_cyan),
            "──────────────────────────────────────────────────────────────────"
            "──"));

    return 0;
}
