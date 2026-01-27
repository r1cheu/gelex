#ifndef GELEX_CLI_CLI_HELPER_H_
#define GELEX_CLI_CLI_HELPER_H_

#include <argparse.h>
#include <atomic>
#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include <barkeep.h>
#include <Eigen/Core>

#include "gelex/types/snp_info.h"

namespace gelex::cli
{
struct ChrGroup
{
    std::string name;
    std::vector<std::pair<Eigen::Index, Eigen::Index>> ranges;
    Eigen::Index total_snps;
};

struct ProgressBarDisplay
{
    std::shared_ptr<barkeep::CompositeDisplay> display;
    std::shared_ptr<barkeep::StatusDisplay> status;
};

struct RemlSpinnerDisplay
{
    std::shared_ptr<barkeep::CompositeDisplay> display;
    std::shared_ptr<barkeep::StatusDisplay> status;
};

auto is_tty() -> bool;

auto setup_parallelization(int num_threads) -> void;

auto build_chr_groups(bool do_loco, const gelex::SnpEffects& snp_effects)
    -> std::vector<ChrGroup>;

auto create_progress_bar(
    std::atomic<size_t>& counter,
    size_t total,
    std::string_view format = "{bar}") -> ProgressBarDisplay;

auto create_reml_spinner() -> RemlSpinnerDisplay;

auto print_gelex_banner_message(std::string_view version) -> void;

auto print_fit_header(
    std::string_view model_name,
    bool has_dominance,
    int iters,
    int burn_in,
    int threads) -> void;

auto print_grm_header(
    std::string_view method,
    bool do_additive,
    bool do_dominant,
    int chunk_size,
    int threads) -> void;

auto print_assoc_header(int threads) -> void;

}  // namespace gelex::cli

#endif  // GELEX_CLI_CLI_HELPER_H_
