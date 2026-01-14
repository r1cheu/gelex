#pragma once

#include <string>
#include <string_view>

#include <argparse.h>

namespace gelex::cli
{
/**
 * @brief Check if the output is a terminal (TTY)
 * @return true if stdout is a terminal, false otherwise
 */
bool is_tty();

std::string repeat(size_t n, std::string_view str);

void print_banner_message(std::string_view version);

void print_fit_header(
    std::string_view version,
    std::string_view model_name,
    bool has_dominance,
    int iters,
    int burnin,
    int threads);

void print_grm_header(
    std::string_view version,
    std::string_view method,
    bool do_additive,
    bool do_dominant,
    int chunk_size,
    int threads);

}  // namespace gelex::cli
