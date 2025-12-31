#include "gelex/cli/utils.h"

#include <unistd.h>
#include <algorithm>

#include <fmt/format.h>

#include "gelex/logger.h"

namespace gelex::cli
{

bool is_tty()
{
    return isatty(fileno(stdout)) != 0;
}

std::string repeat(size_t n, std::string_view str)
{
    if (n <= 0 || str.empty())
    {
        return "";
    }
    std::string result;
    result.reserve(n * str.length());
    for (size_t i = 0; i < n; ++i)
    {
        result.append(str);
    }
    return result;
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

void print_fit_header(
    std::string_view version,
    std::string_view model_name,
    bool has_dominance,
    int iters,
    int burnin,
    int threads)
{
    auto logger = gelex::logging::get();

    std::string title
        = fmt::format("gelex v{} :: Model Fitting (MCMC)", version);
    std::string model_str = fmt::format(
        "Bayes{} ({})",
        model_name,
        has_dominance ? "Additive + Dominance" : "Additive");
    std::string chain_str = fmt::format(
        "{} iters ({} burn-in, {} sampling)", iters, burnin, iters - burnin);
    std::string comp_str = fmt::format("{} threads (OpenMP)", threads);

    size_t min_width = 69;
    size_t content_prefix_len = 17;
    size_t max_content_len
        = std::max({model_str.length(), chain_str.length(), comp_str.length()});

    size_t required_width = content_prefix_len + max_content_len + 1 + 1;
    size_t width = std::max(min_width, required_width);

    width = std::max(width, title.length() + 8);

    size_t top_dashes = width - 7 - title.length();
    logger->info("");  // empty line
    logger->info(" ┌── {} {}┐", title, repeat(top_dashes, "─"));

    // logger->info(" │{:^{}}│", "", width - 3); empty line

    auto log_row = [&](std::string_view label, std::string_view val)
    {
        size_t pad = width - 18 - val.length();
        logger->info(" │  {:<11}: {}{:<{}}│", label, val, "", pad);
    };

    log_row("Model", model_str);
    log_row("Chain", chain_str);
    log_row("Computing", comp_str);
    logger->info(" └{}┘", repeat(width - 3, "─"));
}

}  // namespace gelex::cli
