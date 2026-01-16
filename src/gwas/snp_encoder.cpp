#include "gelex/gwas/snp_encoder.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace gelex::gwas
{

namespace
{

// Encode additive effect: standardize {0, 1, 2}
auto encode_additive(Eigen::Ref<Eigen::VectorXd> geno) -> double
{
    const auto n = static_cast<double>(geno.size());
    const double mean = geno.mean();
    const double p = mean / 2.0;  // allele frequency

    const double var = (geno.array() - mean).square().sum() / (n - 1.0);
    const double sd = std::sqrt(var);

    if (sd < std::numeric_limits<double>::epsilon())
    {
        return p;  // monomorphic, but return MAF anyway
    }

    geno.array() = (geno.array() - mean) / sd;
    return p;
}

// Encode dominance effect using orthogonal HWE encoding (Vitezica et al.)
// This ensures Cov(a, d) ≈ 0 under HWE
// Encoding: 0 → 0, 1 → 2p, 2 → 4p - 2
// Then standardize
auto encode_dominance_orthogonal(Eigen::Ref<Eigen::VectorXd> geno, double p)
    -> void
{
    const auto n = static_cast<double>(geno.size());

    // Expected values under HWE
    const double mu = 2.0 * p * p;           // E[d]
    const double var = 2.0 * p * (1.0 - p);  // Var[d] under HWE

    if (var < std::numeric_limits<double>::epsilon())
    {
        geno.setZero();
        return;
    }

    // Transform genotype to dominance coding
    const double one_encode = 2.0 * p;
    const double two_encode = 4.0 * p - 2.0;

    geno = geno.unaryExpr(
        [one_encode, two_encode](double x) -> double
        {
            if (x == 1.0)
            {
                return one_encode;
            }
            if (x == 2.0)
            {
                return two_encode;
            }
            return 0.0;  // x == 0
        });

    // Standardize
    const double actual_var = (geno.array() - mu).square().sum() / (n - 1.0);
    const double sd = std::sqrt(actual_var);

    if (sd > std::numeric_limits<double>::epsilon())
    {
        geno.array() = (geno.array() - mu) / sd;
    }
    else
    {
        geno.setZero();
    }
}

}  // namespace

auto encode_snp(Eigen::Ref<Eigen::VectorXd> raw, GwasModel model) -> EncodedSNP
{
    EncodedSNP result;
    const auto n = raw.size();

    // Calculate MAF first (before any transformation)
    const double p = raw.mean() / 2.0;
    result.maf = std::min(p, 1.0 - p);

    // Check for monomorphic
    const double var = (raw.array() - raw.mean()).square().sum()
                       / (static_cast<double>(n) - 1.0);
    if (var < std::numeric_limits<double>::epsilon())
    {
        result.is_valid = false;
        result.Z.resize(n, model == GwasModel::AdditiveDominance ? 2 : 1);
        result.Z.setZero();
        return result;
    }

    switch (model)
    {
        case GwasModel::Additive:
        {
            result.Z.resize(n, 1);
            result.Z.col(0) = raw;
            encode_additive(result.Z.col(0));
            break;
        }
        case GwasModel::Dominance:
        {
            result.Z.resize(n, 1);
            result.Z.col(0) = raw;
            encode_dominance_orthogonal(result.Z.col(0), p);
            break;
        }
        case GwasModel::AdditiveDominance:
        {
            result.Z.resize(n, 2);
            // Column 0: additive
            result.Z.col(0) = raw;
            encode_additive(result.Z.col(0));
            // Column 1: dominance (use original raw for frequency)
            result.Z.col(1) = raw;
            encode_dominance_orthogonal(result.Z.col(1), p);
            break;
        }
    }

    return result;
}

auto parse_gwas_model(std::string_view model_str) -> GwasModel
{
    if (model_str == "a" || model_str == "add" || model_str == "additive")
    {
        return GwasModel::Additive;
    }
    if (model_str == "d" || model_str == "dom" || model_str == "dominance")
    {
        return GwasModel::Dominance;
    }
    if (model_str == "ad" || model_str == "a+d" || model_str == "full")
    {
        return GwasModel::AdditiveDominance;
    }
    throw std::invalid_argument(
        "Invalid GWAS model: " + std::string(model_str)
        + ". Valid options: a, d, ad");
}

}  // namespace gelex::gwas
