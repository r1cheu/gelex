#ifndef GELEX_GWAS_ASSOCIATION_TEST_H
#define GELEX_GWAS_ASSOCIATION_TEST_H

#include <cmath>

#include <Eigen/Core>

#include "gelex/gwas/snp_encoder.h"

namespace gelex::gwas
{

enum class TestType
{
    Joint,    // Joint test H₀: a = d = 0 (for a+d model)
    Separate  // Separate tests for a and d
};

struct AssociationResult
{
    double beta_a{std::nan("")};
    double se_a{std::nan("")};
    double beta_d{std::nan("")};
    double se_d{std::nan("")};
    double stat{};    // Wald statistic
    double pvalue{};  // Joint p-value or single effect p-value
    double pvalue_a{
        std::nan("")};  // Separate p-value for a (if TestType::Separate)
    double pvalue_d{
        std::nan("")};  // Separate p-value for d (if TestType::Separate)
    int df{};           // Degrees of freedom (1 or 2)
};

// Perform Wald test for a single SNP
// residual: y - Xβ̂ (phenotype minus fixed effects)
// v_inv: pre-computed V⁻¹ from null model
auto wald_test(
    const EncodedSNP& snp,
    const Eigen::Ref<const Eigen::VectorXd>& residual,
    const Eigen::Ref<const Eigen::MatrixXd>& v_inv,
    GwasModel model,
    TestType test_type) -> AssociationResult;

// Parse test type string to enum
auto parse_test_type(std::string_view test_str) -> TestType;

}  // namespace gelex::gwas

#endif  // GELEX_GWAS_ASSOCIATION_TEST_H
