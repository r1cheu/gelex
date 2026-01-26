#include "assoc_command.h"

#include <filesystem>
#include <memory>
#include <ranges>

#include <fmt/format.h>
#include <Eigen/Core>

#include "cli/cli_helper.h"
#include "cli/gwas_runner.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/data_pipe.h"
#include "gelex/logger.h"

#include "data/loader/bim_loader.h"

auto assoc_execute(argparse::ArgumentParser& cmd) -> int
{
    auto logger = gelex::logging::get();
    std::string out_prefix = cmd.get("--out");

    gelex::cli::setup_parallelization(cmd.get<int>("--threads"));

    auto grm_paths = std::ranges::to<std::vector<std::filesystem::path>>(
        cmd.get<std::vector<std::string>>("--grm"));

    auto bed_path = gelex::BedPipe::format_bed_path(cmd.get("bfile"));

    gelex::DataPipe::Config config{
        .phenotype_path = cmd.get("--pheno"),
        .phenotype_column = cmd.get<int>("--pheno-col"),
        .bed_path = bed_path,
        .use_dominance_effect = false,
        .use_mmap = false,
        .chunk_size = cmd.get<int>("--chunk-size"),
        .qcovar_path = cmd.get("--qcovar"),
        .dcovar_path = cmd.get("--dcovar"),
        .iid_only = cmd.get<bool>("--iid-only"),
        .output_prefix = cmd.get("--out"),
        .grm_paths = grm_paths};

    gelex::cli::print_assoc_header(cmd.get<int>("--threads"));

    auto data_pipe = gelex::load_data_for_reml(config);
    auto sample_manager = data_pipe.sample_manager();

    gelex::BedPipe bed_pipe(bed_path, sample_manager);
    auto bim_path = bed_path;
    bim_path.replace_extension(".bim");
    auto snp_effects
        = std::move(gelex::detail::BimLoader(bim_path)).take_info();

    gelex::cli::GwasRunner::Config runner_config{
        .max_iter = cmd.get<int>("--max-iter"),
        .tol = cmd.get<double>("--tol"),
        .chunk_size = cmd.get<int>("--chunk-size"),
        .loco = cmd.get<bool>("--loco"),
        .additive = cmd.get("--model") == "a",
        .grm_paths = grm_paths,
        .out_prefix = out_prefix};

    gelex::cli::GwasRunner runner(
        runner_config,
        std::move(data_pipe),
        std::move(bed_pipe),
        std::move(snp_effects));

    runner.run();

    return 0;
}
