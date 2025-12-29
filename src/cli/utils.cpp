#include "gelex/cli/utils.h"

#include <unistd.h>

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
