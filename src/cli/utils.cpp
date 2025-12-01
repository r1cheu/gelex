#include "gelex/cli/utils.h"

#include <unistd.h>
#include <string_view>

#include <argparse.h>

#include "gelex/logger.h"

namespace gelex::cli
{

bool is_tty()
{
    return isatty(fileno(stdout)) != 0;
}

std::string repeat(int n, std::string_view str)
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
    for (int i = 1; i < argc; ++i)
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
    std::span<const std::string> cmd)
{
    gelex::logging::initialize(subcommand.get("--out"));
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

void print_banner_message(std::string_view version)
{
    std::cout << "Gelex [version " << version
              << "] - High-Performance Genomic Prediction with Bayesian and "
                 "Frequentist Models\n\n";
    std::cout
        << R"(Gelex is a specialized CLI tool designed for genomic selection and prediction in breeding
programs and quantitative genetics research. Built on modern C++23 with memory-mapped I/O
and BLAS/LAPACK acceleration, Gelex offers seamless integration with PLINK binary formats
and efficient processing of large-scale genomic data.

Basic Usage:
    Fit a Bayesian model:
    $ gelex fit --bfile genotypes --pheno phenotypes.tsv --method RR --out results

    Run simulations:
    $ gelex simulate [options]

Found a Bug or Have a Feature Request?
    Open an issue at: https://github.com/r1cheu/gelex/issues

For more information, see the documentation at: https://github.com/r1cheu/gelex
)";
}

}  // namespace gelex::cli
