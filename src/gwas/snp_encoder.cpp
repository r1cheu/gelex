#include "gelex/gwas/snp_encoder.h"

#include <stdexcept>

#include "gelex/data/grm_code_policy.h"

namespace gelex::gwas
{

auto encode_snp(Eigen::Ref<Eigen::VectorXd> raw, AssocMode model) -> EncodedSNP
{
    EncodedSNP result;
    const auto n = raw.size();
    const double p = raw.mean() / 2.0;
    result.maf = std::min(p, 1.0 - p);

    switch (model)
    {
        case AssocMode::A:
        {
            result.Z.resize(n, 1);
            result.Z.col(0) = raw;
            gelex::grm::yang_add_impl(result.Z.col(0));
            break;
        }
        case AssocMode::D:
        {
            result.Z.resize(n, 1);
            result.Z.col(0) = raw;
            gelex::grm::yang_dom_impl(result.Z.col(0));
            break;
        }
        case AssocMode::AD:
        {
            result.Z.resize(n, 2);
            // Column 0: additive
            result.Z.col(0) = raw;
            gelex::grm::yang_add_impl(result.Z.col(0));
            // Column 1: dominance (use original raw for frequency)
            result.Z.col(1) = raw;
            gelex::grm::yang_dom_impl(result.Z.col(1));
            break;
        }
    }

    return result;
}

auto parse_assoc_mode(std::string_view model_str) -> AssocMode
{
    if (model_str == "a" || model_str == "add" || model_str == "additive")
    {
        return AssocMode::A;
    }
    if (model_str == "d" || model_str == "dom" || model_str == "dominance")
    {
        return AssocMode::D;
    }
    if (model_str == "ad" || model_str == "a+d" || model_str == "full")
    {
        return AssocMode::AD;
    }
    throw std::invalid_argument(
        "Invalid GWAS model: " + std::string(model_str)
        + ". Valid options: a, d, ad");
}

}  // namespace gelex::gwas
