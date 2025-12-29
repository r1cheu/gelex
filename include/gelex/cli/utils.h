#pragma once

#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <argparse.h>

namespace gelex::cli
{
/**
 * @brief Check if the output is a terminal (TTY)
 * @return true if stdout is a terminal, false otherwise
 */
bool is_tty();

std::string repeat(int n, std::string_view str);

void print_banner_message(std::string_view version);

}  // namespace gelex::cli
