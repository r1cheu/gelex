#ifndef GELEX_CLI_CLI_HELPER_H_
#define GELEX_CLI_CLI_HELPER_H_

#include <argparse.h>
#include <string_view>

namespace gelex::cli
{
/**
 * @brief Check if the output is a terminal (TTY)
 * @return true if stdout is a terminal, false otherwise
 */
bool is_tty();

void setup_parallelization(int num_threads);

void print_gelex_banner_message(std::string_view version);

void print_fit_header(
    std::string_view model_name,
    bool has_dominance,
    int iters,
    int burn_in,
    int threads);

void print_grm_header(
    std::string_view method,
    bool do_additive,
    bool do_dominant,
    int chunk_size,
    int threads);

void print_assoc_header(int threads);

}  // namespace gelex::cli

#endif  // GELEX_CLI_CLI_HELPER_H_
