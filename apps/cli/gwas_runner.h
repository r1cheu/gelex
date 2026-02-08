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

#ifndef GELEX_CLI_GWAS_RUNNER_H
#define GELEX_CLI_GWAS_RUNNER_H

#include <atomic>
#include <filesystem>
#include <functional>
#include <vector>

#include <Eigen/Core>

#include "cli/cli_helper.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/data_pipe.h"
#include "gelex/gwas/gwas_writer.h"
#include "gelex/types/assoc_input.h"
#include "gelex/types/snp_info.h"
#include "logger/loco_reml_logger.h"
#include "utils/utils.h"

namespace gelex
{
class FreqModel;
class FreqState;
}  // namespace gelex

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
        GenotypeProcessMethod method;
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
    auto print_assoc_summary() const -> void;

    auto run_normal() -> void;

    auto run_loco() -> void;

    auto update_assoc_input(
        const FreqModel& model,
        const FreqState& state,
        Eigen::MatrixXd&& v_inv) -> void;

    auto scan_chromosome(
        const ChrGroup& group,
        std::atomic<size_t>& progress_counter,
        size_t total_snps_to_report,
        size_t total_processed_before,
        const std::function<void(size_t, size_t, size_t)>& progress_callback)
        -> void;

    auto print_scan_summary() const -> void;

    Config config_;
    DataPipe data_pipe_;
    BedPipe bed_pipe_;
    gwas::GwasWriter writer_;
    SnpEffects snp_effects_;

    SmoothEtaCalculator eta_calculator_;
    std::vector<ChrGroup> chr_groups_;

    AssocInput assoc_input_;
    AssocOutput assoc_output_;
    Eigen::VectorXd freqs_;

    std::vector<detail::LocoRemlResult> loco_results_;
};

}  // namespace gelex::cli

#endif  // GELEX_CLI_GWAS_RUNNER_H
