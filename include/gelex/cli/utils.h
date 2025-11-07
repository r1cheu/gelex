#pragma once

#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <gelex/argparse.h>

namespace gelex::cli
{
/**
 * @brief Check if the output is a terminal (TTY)
 * @return true if stdout is a terminal, false otherwise
 */
bool is_tty();

std::string repeat(int n, std::string_view str);

std::vector<std::string> parse_command(int argc, char* args[]);

void log_command(
    const argparse::ArgumentParser& subcommand,
    std::span<const std::string> cmd);

}  // namespace gelex::cli
