#ifndef GELEX_CLI_GWAS_RUNNER_H
#define GELEX_CLI_GWAS_RUNNER_H

#include <atomic>
#include <filesystem>
#include <functional>
#include <vector>

#include <Eigen/Core>

#include "barkeep.h"
#include "cli/cli_helper.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/data_pipe.h"
#include "gelex/data/grm_code_policy.h"
#include "gelex/data/loco_grm_loader.h"
#include "gelex/estimator/freq/reml.h"
#include "gelex/gwas/gwas_writer.h"
#include "gelex/types/assoc_input.h"
#include "gelex/types/snp_info.h"
#include "utils/utils.h"

namespace gelex::cli
{

class GwasRunner
{
   public:
    struct Config
    {
        int max_iter;
        double tol;
        int chunk_size;
        bool loco;
        bool additive;
        std::vector<std::filesystem::path> grm_paths;
        std::string out_prefix;
    };

    GwasRunner(
        Config config,
        DataPipe data_pipe,
        BedPipe bed_pipe,
        SnpEffects snp_effects);

    auto run() -> void;

   private:
    auto print_summary() const -> void;

    auto run_normal() -> void;

    auto run_loco() -> void;

    auto scan_chromosome(
        AssocInput& input,
        const ChrGroup& group,
        std::atomic<size_t>& progress_counter,
        size_t total_snps,
        const std::function<void(size_t, size_t)>& progress_callback) -> void;

    Config config_;
    DataPipe data_pipe_;
    BedPipe bed_pipe_;
    gwas::GwasWriter writer_;
    SnpEffects snp_effects_;

    grm::Vitezica encoder_;
    SmoothEtaCalculator eta_calculator_;
    std::vector<ChrGroup> chr_groups_;
};

}  // namespace gelex::cli

#endif  // GELEX_CLI_GWAS_RUNNER_H
