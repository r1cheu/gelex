#ifndef GELEX_GWAS_SNP_ENCODER_H
#define GELEX_GWAS_SNP_ENCODER_H

#include <Eigen/Core>

namespace gelex::gwas
{

enum class GwasModel
{
    Additive,          // a only
    Dominance,         // d only
    AdditiveDominance  // a + d
};

struct EncodedSNP
{
    Eigen::MatrixXd Z;  // n × 1 (a or d) or n × 2 (a+d)
    double maf{};
    bool is_valid{true};
};

// Encode raw genotype {0, 1, 2} based on model type
// Uses orthogonal HWE encoding for dominance
// raw genotype will be modified in place for efficiency
auto encode_snp(Eigen::Ref<Eigen::VectorXd> raw, GwasModel model) -> EncodedSNP;

// Parse model string to enum
auto parse_gwas_model(std::string_view model_str) -> GwasModel;

}  // namespace gelex::gwas

#endif  // GELEX_GWAS_SNP_ENCODER_H
