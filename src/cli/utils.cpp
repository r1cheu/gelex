#include "gelex/cli/utils.h"

#include <unistd.h>

#include <fmt/format.h>
#include <omp.h>
#include <Eigen/Core>
#include "config.h"

#include "../src/utils/formatter.h"
#include "gelex/logger.h"

namespace gelex::cli
{

bool is_tty()
{
    return isatty(fileno(stdout)) != 0;
}

void setup_parallelization(int num_threads)
{
    if (num_threads > 0)
    {
        omp_set_num_threads(num_threads);
        Eigen::setNbThreads(num_threads);
    }
}

void print_gelex_banner_message(std::string_view version)
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
    std::string_view model_name,
    bool has_dominance,
    int iters,
    int burn_in,
    int threads)
{
    auto logger = gelex::logging::get();

    std::string title
        = fmt::format("gelex v{} :: Model Fitting (MCMC)", PROJECT_VERSION);
    std::string model_str = fmt::format(
        "Bayes{} ({})",
        model_name,
        has_dominance ? "Additive + Dominance" : "Additive");
    std::string chain_str = fmt::format(
        "{} iters ({} burn-in, {} sampling)", iters, burn_in, iters - burn_in);
    std::string comp_str = fmt::format("{}", threads);

    std::vector<std::pair<std::string, std::string>> items
        = {{"Model", model_str}, {"Chain", chain_str}, {"Threads", comp_str}};

    logger->info(gelex::header_box(title, items, 70));
    logger->info("");
}

void print_grm_header(
    std::string_view method,
    bool do_additive,
    bool do_dominant,
    int chunk_size,
    int threads)
{
    auto logger = gelex::logging::get();

    std::string title
        = fmt::format("gelex v{} :: GRM Computation", PROJECT_VERSION);

    std::string mode_str;
    if (do_additive && do_dominant)
    {
        mode_str = "Additive + Dominance";
    }
    else if (do_additive)
    {
        mode_str = "Additive";
    }
    else
    {
        mode_str = "Dominance";
    }

    std::string method_str = std::string(method);
    std::string chunk_str = fmt::format("{}", chunk_size);
    std::string comp_str = fmt::format("{}", threads);

    std::vector<std::pair<std::string, std::string>> items
        = {{"Method", method_str},
           {"Mode", mode_str},
           {"Chunk Size", chunk_str},
           {"Threads", comp_str}};

    logger->info(gelex::header_box(title, items, 70));
    logger->info("");
}

void print_assoc_header(int threads)
{
    auto logger = gelex::logging::get();

    std::string title
        = fmt::format("gelex v{} :: GWAS Analysis", PROJECT_VERSION);
    std::vector<std::pair<std::string, std::string>> header_items
        = {{"Method", "AI-REML (Average Information)"},
           {"Threads", fmt::format("{}", threads)}};
    logger->info(gelex::header_box(title, header_items, 70));
    logger->info("");
}

}  // namespace gelex::cli
