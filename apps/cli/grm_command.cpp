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
#include "gelex/data/genotype_processor.h"
#include "gelex/data/grm.h"
#include "gelex/logger.h"
#include "grm_args.h"
#include "utils/formatter.h"
#include "utils/utils.h"

namespace bk = barkeep;

namespace
{

template <typename Method>
auto dispatch_grm(
    gelex::GRM& grm,
    const std::vector<std::pair<Eigen::Index, Eigen::Index>>& ranges,
    int chunk_size,
    bool additive,
    const std::function<void(Eigen::Index, Eigen::Index)>& progress_callback
    = nullptr) -> gelex::GrmResult
{
    if (additive)
    {
        return grm.compute<typename Method::Additive>(
            ranges, chunk_size, progress_callback);
    }
    return grm.compute<typename Method::Dominant>(
        ranges, chunk_size, progress_callback);
}

auto compute_grm_with_method(
    gelex::GRM& grm,
    const std::vector<std::pair<Eigen::Index, Eigen::Index>>& ranges,
    std::string_view method,
    int chunk_size,
    bool additive,
    const std::function<void(Eigen::Index, Eigen::Index)>& progress_callback
    = nullptr) -> gelex::GrmResult
{
    if (method == "1")
    {
        return dispatch_grm<gelex::grm::OrthStandardized>(
            grm, ranges, chunk_size, additive, progress_callback);
    }
    if (method == "2")
    {
        return dispatch_grm<gelex::grm::Centered>(
            grm, ranges, chunk_size, additive, progress_callback);
    }
    if (method == "3")
    {
        return dispatch_grm<gelex::grm::OrthCentered>(
            grm, ranges, chunk_size, additive, progress_callback);
    }
    throw std::invalid_argument(
        fmt::format("Unknown GRM method: {}. Valid: 1, 2, 3", method));
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
    GrmConfig config{
        .bed_path = gelex::BedPipe::format_bed_path(cmd.get("--bfile")),
        .out_prefix = cmd.get("--out"),
        .method = cmd.get("--method"),
        .chunk_size = cmd.get<int>("--chunk-size"),
        .do_additive = cmd.get<bool>("--add"),
        .do_dominant = cmd.get<bool>("--dom"),
        .do_loco = cmd.get<bool>("--loco"),
        .threads = cmd.get<int>("--threads")};

    if (!config.do_additive && !config.do_dominant)
    {
        config.do_additive = true;
    }

    gelex::cli::setup_parallelization(config.threads);

    gelex::cli::print_grm_header(
        config.method,
        config.do_additive,
        config.do_dominant,
        config.chunk_size,
        config.threads);

    logger->info(gelex::section("Loading Data..."));
    logger->info(gelex::success("Input      : {}", config.bed_path.string()));

    gelex::GRM grm(config.bed_path);
    const auto& sample_ids = grm.sample_ids();
    auto num_snps = grm.num_snps();
    logger->info(gelex::success("Samples    : {} samples", sample_ids.size()));
    logger->info(gelex::success("SNPs       : {} markers", num_snps));

    logger->info("");
    logger->info(gelex::section("Computing GRM..."));

    auto bim_path = config.bed_path;
    bim_path.replace_extension(".bim");
    gelex::detail::BimLoader bim_loader(bim_path);

    auto groups
        = gelex::cli::build_chr_groups(config.do_loco, bim_loader.info());

    struct GrmTask
    {
        std::string name;
        std::string label;
        bool is_additive;
    };

    std::vector<GrmTask> tasks;
    if (config.do_additive)
    {
        tasks.push_back({"add", "Additive", true});
    }
    if (config.do_dominant)
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

                pbar.after->message(
                    fmt::format(
                        "{:.1f}% ({}/{}) | {}",
                        static_cast<double>(current_total)
                            / static_cast<double>(total_work_snps) * 100,
                        gelex::HumanReadable(current_total),
                        gelex::HumanReadable(total_work_snps),
                        eta_calculator.get_eta(current_total)));
            };

            auto result = compute_grm_with_method(
                grm,
                group.ranges,
                config.method,
                config.chunk_size,
                task.is_additive,
                progress_callback);

            auto path = config.out_prefix;
            if (config.do_loco || tasks.size() > 1)
            {
                auto suffix
                    = config.do_loco ? fmt::format(".chr{}", group.name) : "";
                path = fmt::format(
                    "{}.{}{}", config.out_prefix, task.name, suffix);
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
        fmt::format(
            fmt::fg(fmt::color::light_cyan),
            "── Computation Summary {}",
            gelex::separator(70 - 23)));
    logger->info(
        gelex::success(
            "Time elapsed: {}", eta_calculator.total_time_consumed()));

    logger->info("  Total Files : {}", generated_files.size());
    logger->info(
        "  Output Dir  : {}",
        std::filesystem::absolute(std::filesystem::path(config.out_prefix))
            .parent_path()
            .string());

    if (config.do_loco)
    {
        logger->info(
            "  Pattern     : {}.{{add|dom}}.chr{{1..{}}}.{{bin|id}}",
            config.out_prefix,
            groups.size());
    }
    else
    {
        logger->info(
            "  Pattern     : {}.{{add|dom}}.{{bin|id}}", config.out_prefix);
    }
    logger->info(
        fmt::format(
            fmt::emphasis::bold | fmt::fg(fmt::color::light_cyan),
            "──────────────────────────────────────────────────────────────"
            "────────"));
    return 0;
}
