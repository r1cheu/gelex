#include "gelex/cli/utils.h"

#include <unistd.h>
#include <string_view>

#include <gelex/argparse.h>

#include "gelex/logger.h"

namespace gelex::cli
{

bool is_tty()
{
    return isatty(fileno(stdout)) != 0;
}

void print_banner(std::string_view version)
{
    auto logger = gelex::logging::get();
    logger->info("╔═══════════════════════════════════════════╗");
    logger->info(
        "║  G E L E X  v{}{}║",
        version,
        std::string(29 - version.length(), ' '));
    logger->info("║  ▔▔▔▔▔▔▔▔▔▔                               ║");
    logger->info("║  Genomic Prediction with Bayes Models     ║");
    logger->info("╚═══════════════════════════════════════════╝");
}

std::string repeat(size_t n, std::string_view str)
{
    if (n <= 0 || str.empty())
    {
        return "";
    }
    std::string result;
    result.reserve(n * str.length());
    for (int i = 0; i < n; ++i)
    {
        result.append(str);
    }
    return result;
}

std::vector<std::string> parse_command(int argc, char* args[])
{
    std::vector<std::string> command;
    command.emplace_back("gelex");
    for (size_t i = 1; i < argc; ++i)
    {
        std::string arg = args[i];
        if (arg.starts_with("--") || arg.starts_with("-"))
        {
            std::string option = "  " + arg;
            if (i + 1 < argc)
            {
                std::string next_arg = args[i + 1];
                if (!next_arg.starts_with("--") && !next_arg.starts_with("-"))
                {
                    option += " " + next_arg;
                    i++;
                }
            }
            command.push_back(option);
        }
        else
        {
            command.push_back(arg);
        }
    }
    return command;
}

void log_command(
    const argparse::ArgumentParser& subcommand,
    std::string_view version,
    std::span<const std::string> cmd)
{
    gelex::logging::initialize(subcommand.get("--out"));
    print_banner(version);

    auto logger = gelex::logging::get();

    logger->info("");
    size_t max_length = 0;
    for (std::string_view part : cmd)
    {
        max_length = std::max(max_length, part.length());
    }
    max_length += 2;
    max_length = std::max(max_length, static_cast<size_t>(20));

    std::string top_border = " ── Command " + repeat(max_length - 10, "─");
    std::string bottom_border = " " + repeat(max_length, "─");

    logger->info("{}", top_border);
    for (std::string_view part : cmd)
    {
        logger->info("   {:<{}} ", part, max_length);
    }
    logger->info("{}", bottom_border);
    logger->info("");
}

}  // namespace gelex::cli
